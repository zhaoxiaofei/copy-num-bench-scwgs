#!/usr/bin/env python
# A new script of https://github.com/zhaoxiaofei/copy-num-bench-scwgs

"""
Benchmark ploidy estimation across CNV calling methods and datasets, in the style of
bench_results/scWGS-performances-eval.py, but on the ploidy-evaluation outputs of
ploidy_eval.py / ploidy_tools.py instead of the CNV-caller performance table.

What is plotted
---------------
The main figure is a grid:
    rows    = CNV calling methods (one row per caller / ploidy-inference tool)
    columns = datasets (one column per donor / sampleType / avgSpotLen group)
    entry   = per-cell ploidy-estimation ERRORS of one method on one dataset:
              one point per cell (observed - expected ploidy), drawn as one jittered
              error cluster per cell sample, exactly the column-wise entries of
              ploidy_eval.plot_ploidy (ploidy_eval.py line 515) translated into error
              space: the grey band is the +/- ploidy window around the expected ploidy
              (error = 0), the crimson line marks the expected ploidy itself, and the
              blue stars mark the 2x / 0.5x scaling-error references.
The row / column labelling follows plot_grid_main of scWGS-performances-eval.py (its
line 153): the method name is the y-label of the first column, the dataset name is the
title of the first row, inner tick labels are hidden, and every entry shares the y
scale so methods and datasets stay comparable.
Two auxiliary figures summarise the same errors as box plots (pooled over datasets,
and one row per dataset), mirroring plot_main / plot_multirow_main of
scWGS-performances-eval.py.

Input
-----
The main input files are the ploidy-evaluation summaries written by ploidy_eval.py and
ploidy_tools.py:
    *_ploidy_eval_summary.json      (CNV callers, e.g. via data_tumor.py --ploidy-file)
    *_ploidy_tool_eval_summary.json (ploidy-inference tools, e.g. scAbsolute)
each carrying the run identity (tool, donor, sampleType, avgSpotLen, ploidy_window,
max_cn, ...).  The per-cell errors themselves live in the sibling file that shares the
summary's prefix, `<prefix>_percell.tsv` (columns `sample`, `ploidy_error`,
`expected_ploidy`, ...), which this script resolves automatically; both writers
(ploidy_eval.py, ploidy_tools.py) always emit the pair.  When a caller was evaluated
both by ploidy_eval.py and by `ploidy_tools.py eval` (identical per-cell numbers), the
`*_ploidy_eval_*` run is kept and the duplicate is dropped with a warning.

Usage
-----
    python bench_results/scWGS-ploidy-performances-eval.py \
        -i '../data/*/4from2_*_ploidy_eval_summary.json' \
           '../data/*/4from2_*_ploidy_tool_eval_summary.json' \
        -o bench_results/ploidy-performances
    # or feed the file list (or glob patterns) on stdin, like scWGS-performances-eval.py:
    find ../data -name '*_ploidy*_eval*_summary.json' | \
        python bench_results/scWGS-ploidy-performances-eval.py -o bench_results/ploidy-performances

Outputs `<output>_ploidy_error_grid.pdf/.png` (the main grid),
`<output>_ploidy_error_main.pdf/.png` and `<output>_ploidy_error_multirow.pdf/.png`
(the box-plot summaries), plus a per-run statistics table on stderr.
"""

import argparse
import collections
import glob
import json
import logging
import os
import re
import sys

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')  # Set non-interactive backend before importing pyplot

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import seaborn as sns

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(filename)s %(levelname)s %(message)s')

# ---------------------------------------------------------------------------------------------
# Constants (kept in step with ploidy_eval.py so the two scripts cannot drift apart)
# ---------------------------------------------------------------------------------------------

SUMMARY_SUFFIX = '_summary.json'
PERCELL_SUFFIX = '_percell.tsv'
# scAbsolute's window around the experimental point estimate (ploidy_eval.py).
DEFAULT_PLOIDY_WINDOW = 0.5
# Cap applied to the per-segment copy number before it is averaged into a ploidy
# (ploidy_eval.py).  The pipeline evaluates every run at both the default cap and `inf`
# and writes them to separate prefixes, hence the --max-cn selection below.
DEFAULT_MAX_CN = 10.0

# The row order of scWGS-performances-eval.py / cnv_gather_results.py, with the
# ploidy-inference tool of ploidy_tools.py first.  Unknown tools are appended after
# these, alphabetically.
CALLER_ORDER = ['scabsolute', 'hmmcopy', 'ginkgo', 'copynumber', 'secnv',
                'sccnv', 'scyn', 'chisel', 'aneufinder', 'flcna']

caller2desc = {
    'aneufinder': 'AneuFinder  Genome Biology               2016',
    'flcna'     : 'FLCNA       Genome Research              2024',
    'chisel'    : 'Chisel      Nature Biotechnology         2021',
    'copynumber': 'CopyNumber  BMC Genomics                 2012',
    'ginkgo'    : 'Ginkgo      Nature Methods               2015',
    'hmmcopy'   : 'HMMcopy     Bioinformatics               2006',
    'secnv'     : 'SeCNV       Briefings in Bioinformatics  2022',
    'sccnv'     : 'SCCNV       Frontiers in Genetics        2020',
    'scyn'      : 'SCYN/SCOPE  Cell Systems                 2020',
    'scabsolute': 'scAbsolute  Genome Biology               2024',
}
# Annotated row labels, off by default exactly as in scWGS-performances-eval.py:
# empty this dict to keep it off, or leave it filled to switch the annotations on.
caller2desc = {}

# Fallback parser for the pipeline file names, used only when the summary JSON carries
# no donor / sampleType / avgSpotLen / tool (i.e. standalone runs without --tool etc.):
#   4from2_2_<donor>_3_<sampleType>_<avgSpotLen>_4_step<n>_<tool>_ploidy_eval[_maxcn_<v>]
#   4from3_2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>_4_step<n>_<tool>_ploidy_eval[...]
# The simulated <cellLine> is deliberately NOT part of the dataset identity: on
# simulated data it is the per-cell sample label, so it shows up as one error cluster
# per cell line inside the dataset's column, like every other cell sample.
STEM_RE = re.compile(
    r'^(?:4from[23]_2_)?'
    r'(?P<donor>.+?)_3_(?P<sampleType>.+?)_(?P<avgSpotLen>\d{2,})'
    r'(?:_(?P<cellLine>.+?))?_4_step\d+_(?P<tool>.+?)_ploidy(?:_tool)?_eval'
    r'(?:_maxcn_(?P<maxcn>[^_]+))?$')

parser1 = argparse.ArgumentParser()
parser1.add_argument('-i', '--input', nargs='+', default=None,
                     help='Ploidy-evaluation summary files or globs (typically '
                          '*_ploidy*_eval*_summary.json); each must sit next to its '
                          '<prefix>_percell.tsv sibling.  Read from stdin instead when '
                          'omitted or "-". ')
parser1.add_argument('-t', '--type', type=int, default=0,
                     help='Output type. 0: all figures. 1: testing (first method and '
                          'first dataset only). 2: only plot the main grid. ')
parser1.add_argument('-o', '--output', default='scWGS-ploidy-performances')
parser1.add_argument('--max-cn', default=str(int(DEFAULT_MAX_CN)),
                     help='Which copy-number cap to plot, given that the pipeline '
                          'evaluates every run at several caps: a number (e.g. 10), '
                          'inf, or all (one extra dataset column per additional cap). '
                          'Runs at other caps are skipped. ')
parser1.add_argument('--methods', nargs='+', default=None,
                     help='Restrict (and order) the rows to these tools, by exact name. ')
parser1.add_argument('--datasets', nargs='+', default=None,
                     help='Restrict (and order) the columns to the datasets whose label '
                          'contains any of these substrings. ')
parser1.add_argument('--sharey', default='all', choices=['all', 'row', 'none'],
                     help='Share the y scale of the grid entries across the whole figure '
                          '(all, the default), within each method row, or not at all. ')
parser1.add_argument('--y-quantile', type=float, default=None,
                     help='Clip the y axis to this central fraction of the per-cell '
                          'errors (e.g. 0.99 keeps the central 99%%); points beyond the '
                          'limits are drawn but not visible. ')
parser1.add_argument('--ylim', type=float, nargs=2, default=None, metavar=('LOW', 'HIGH'),
                     help='Fixed y limits of the main grid, overriding --sharey/--y-quantile. ')
parser1.add_argument('--jitter', type=float, default=0.28,
                     help='Half-width of the horizontal jitter of each per-cell error cluster. ')
parser1.add_argument('--dpi', type=int, default=300, help='Raster resolution of the PNG outputs. ')

args = parser1.parse_args()

# ---------------------------------------------------------------------------------------------
# Input: the *_ploidy*_eval*_summary.json files and their <prefix>_percell.tsv siblings
# ---------------------------------------------------------------------------------------------

def norm_max_cn(value):
    """A copy-number cap as a float, mapping the uncapped spellings of ploidy_eval.py
    (inf / none / nan / '') onto float('inf') so that caps compare by ==."""
    if value is None:
        return float('inf')
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return float(value)
    s = str(value).strip().lower()
    if s in ('inf', 'infinity', 'none', 'no', 'nan', ''):
        return float('inf')
    return float(s)

def fmt_max_cn(value):
    v = norm_max_cn(value)
    return 'inf' if np.isinf(v) else F'{v:g}'

def is_caller_eval_stem(stem):
    """True for the *_ploidy_eval summaries of ploidy_eval.py (CNV callers), False for
    the *_ploidy_tool_eval summaries of ploidy_tools.py.  A caller named in --ploidy-tools
    gets both on real tumors, with identical per-cell numbers; the caller-side file is
    the canonical one, so duplicates of it are dropped."""
    return ('_ploidy_eval' in stem) and ('_ploidy_tool_eval' not in stem)

def expand_inputs(patterns):
    """File paths from a list of paths / glob patterns, de-duplicated, in a stable order."""
    paths, seen = [], set()
    for pat in patterns:
        pat = pat.strip()
        if not pat or pat.startswith('#'):
            continue
        matched = sorted(glob.glob(pat)) if any(c in pat for c in '*?[') else [pat]
        if not matched:
            logging.warning('input pattern matched no file: %s', pat)
        for path in matched:
            apath = os.path.abspath(path)
            if apath in seen:
                continue
            seen.add(apath)
            paths.append(apath)
    return paths

def load_run(path):
    """One ploidy-evaluation run: its identity, read from the summary JSON, plus the
    per-cell errors of the sibling <prefix>_percell.tsv.  Returns None when the run
    cannot be used, with a warning saying why."""
    base = os.path.basename(path)
    if not base.endswith(SUMMARY_SUFFIX):
        logging.warning('%s: not a *%s file; skipping', path, SUMMARY_SUFFIX)
        return None
    stem = base[:-len(SUMMARY_SUFFIX)]
    if not ('ploidy' in stem and 'eval' in stem):
        logging.warning('%s: does not follow the *_ploidy*_eval*%s naming; skipping', path, SUMMARY_SUFFIX)
        return None
    try:
        with open(path) as fh:
            js = json.load(fh)
    except (OSError, ValueError) as exc:
        logging.warning('%s: unreadable JSON (%s); skipping', path, exc)
        return None
    if not isinstance(js, dict):
        logging.warning('%s: the summary is not a JSON object; skipping', path)
        return None
    percell_path = path[:-len(SUMMARY_SUFFIX)] + PERCELL_SUFFIX
    if not (os.path.isfile(percell_path) and os.path.getsize(percell_path) > 0):
        logging.warning('%s: no usable sibling %s, so no per-cell errors; skipping', path, PERCELL_SUFFIX)
        return None
    try:
        df = pd.read_csv(percell_path, sep='\t')
    except (OSError, ValueError) as exc:
        logging.warning('%s: unreadable per-cell table (%s); skipping', percell_path, exc)
        return None
    missing = [c for c in ('sample', 'ploidy_error') if c not in df.columns]
    if missing:
        logging.warning('%s: per-cell table lacks the column(s) %s; skipping', percell_path, missing)
        return None
    m = STEM_RE.match(stem) or {}
    # Method (row) identity: the summary's --tool, else the per-cell table's tool column
    # (both writers carry it through), else the file name.
    tool = str(js.get('tool') or '').strip()
    if not tool and 'tool' in df.columns:
        nonempty = df['tool'].dropna().astype(str)
        nonempty = nonempty[nonempty.str.strip() != '']
        if len(nonempty):
            tool = str(nonempty.iloc[0]).strip()
    if not tool:
        tool = str(m.get('tool') or '').strip()
    if not tool:
        logging.warning('%s: no tool identity (summary key "tool", per-cell column or '
                        'file name); skipping', path)
        return None
    # Dataset (column) identity: donor / sampleType / avgSpotLen, which is exactly what
    # one summary file covers.  The cell samples inside the dataset (the error clusters)
    # are the `sample` values of the per-cell table: the ploidy-file samples on real
    # tumors, the simulated cell line on simulated data.
    donor = str(js.get('donor') or '').strip()
    sample_type = str(js.get('sampleType') or '').strip()
    avg_spot_len = str(js.get('avgSpotLen') or '').strip()
    if not (donor or sample_type or avg_spot_len):
        # standalone runs: fall back to the provenance columns both writers carry through
        for col, field in (('donor', 'donor'), ('sampleType', 'sample_type'),
                           ('avgSpotLen', 'avg_spot_len')):
            if col in df.columns:
                vals = [str(v).strip() for v in df[col].dropna().astype(str)
                        if str(v).strip() != '']
                if vals:
                    if field == 'donor': donor = vals[0]
                    elif field == 'sample_type': sample_type = vals[0]
                    elif field == 'avg_spot_len': avg_spot_len = vals[0]
    if not (donor or sample_type or avg_spot_len) and m:
        donor, sample_type, avg_spot_len = m['donor'], m['sampleType'], m['avgSpotLen']
    ploidy_error = pd.to_numeric(df['ploidy_error'], errors='coerce').to_numpy(dtype=float)
    expected = (pd.to_numeric(df['expected_ploidy'], errors='coerce').to_numpy(dtype=float)
                if 'expected_ploidy' in df.columns
                else np.full(len(df), np.nan))
    n_finite = int(np.isfinite(ploidy_error).sum())
    if n_finite == 0:
        logging.warning('%s: no cell has a finite ploidy error; skipping', path)
        return None
    window = js.get('ploidy_window', DEFAULT_PLOIDY_WINDOW)
    try:
        window = float(window)
    except (TypeError, ValueError):
        window = DEFAULT_PLOIDY_WINDOW
    return {
        'path': path, 'percell_path': percell_path, 'stem': stem,
        'tool': tool, 'donor': donor, 'sampleType': sample_type, 'avgSpotLen': avg_spot_len,
        'max_cn': norm_max_cn(js.get('max_cn', DEFAULT_MAX_CN)), 'window': window,
        'sample': df['sample'].astype(str).to_numpy(),
        'ploidy_error': ploidy_error, 'expected_ploidy': expected,
        'n_cells': int(len(df)), 'n_cells_finite': n_finite,
        'caller_eval': is_caller_eval_stem(stem),
        # fallback dataset label for standalone runs without donor/sampleType/avgSpotLen
        'stem_label': re.sub(r'_4_step\d+_.+?_ploidy(_tool)?_eval(_maxcn_\S+)?$', '', stem),
    }

def dataset_key(run):
    if run['donor'] or run['sampleType'] or run['avgSpotLen']:
        return (run['donor'], run['sampleType'], run['avgSpotLen'])
    return (run['stem_label'], '', '')

def dataset_label(run, max_cn_tag=''):
    parts = []
    if run['donor']:       parts.append(F"donor={run['donor']}")
    if run['sampleType']:  parts.append(F"sampleType={run['sampleType']}")
    if run['avgSpotLen']:  parts.append(F"avgSpotLen={run['avgSpotLen']}")
    label = ' | '.join(parts) if parts else (run['stem_label'] or 'dataset')
    return label + max_cn_tag

def _avg_spot_len_sort(run_or_key):
    key = run_or_key if isinstance(run_or_key, tuple) else dataset_key(run_or_key)
    try:
        return (0, float(key[2])) if key[2] else (1, 0.0)
    except ValueError:
        return (1, 0.0)

def order_methods(tools):
    known = [t for t in CALLER_ORDER if t in tools]
    rest = sorted(set(tools) - set(CALLER_ORDER))
    return known + rest

def dataset_key_str(key):
    return ' | '.join(p for p in key if p)

# ---------------------------------------------------------------------------------------------
# Assemble the runs into one tidy per-cell table
# ---------------------------------------------------------------------------------------------

def gather_runs():
    patterns = list(args.input or [])
    if not patterns or patterns == ['-']:
        patterns = [line.rstrip('\n') for line in sys.stdin]
    if not patterns:
        parser1.error('no input: pass *_ploidy*_eval*_summary.json files or globs with -i, '
                      'or feed them on stdin')
    paths = expand_inputs(patterns)
    if not paths:
        parser1.error('no existing input file among the given patterns')
    runs = []
    for path in paths:
        run = load_run(path)
        if run is not None:
            runs.append(run)
            logging.info('loaded %s: tool=%s donor=%s sampleType=%s avgSpotLen=%s '
                         'max-cn=%s window=%g cells=%d (finite: %d)',
                         os.path.basename(path), run['tool'], run['donor'],
                         run['sampleType'], run['avgSpotLen'], fmt_max_cn(run['max_cn']),
                         run['window'], run['n_cells'], run['n_cells_finite'])
    if not runs:
        sys.exit('ploidy-performances-eval: no usable *_ploidy*_eval*_summary.json input. ')
    # --max-cn selection: the pipeline evaluates every run at several caps.
    wanted = str(args.max_cn).strip().lower()
    if wanted in ('all', 'any', 'both'):
        wanted = 'all'
    else:
        try:
            wanted = float(wanted)
            if not (wanted > 0):
                raise ValueError
        except ValueError:
            parser1.error(F'--max-cn must be a positive number, inf, or all (got {args.max_cn})')
        before = len(runs)
        runs = [r for r in runs if r['max_cn'] == wanted]
        logging.info('--max-cn %s: kept %d of %d runs (the rest were evaluated at other caps)',
                     fmt_max_cn(wanted), len(runs), before)
        if not runs:
            sys.exit(F'ploidy-performances-eval: every input run is at a copy-number cap other '
                     F'than {fmt_max_cn(wanted)}; retry with --max-cn all (or inf). ')
    # De-duplicate (tool, dataset, max_cn): on real tumors a caller named in --ploidy-tools
    # is evaluated both by ploidy_eval.py and by ploidy_tools.py, with identical per-cell
    # numbers; keep the caller-side *_ploidy_eval run and drop the duplicate.
    runs.sort(key=lambda r: (not r['caller_eval'], r['path']))
    seen, kept, dropped = {}, [], []
    for run in runs:
        key = (run['tool'], dataset_key(run), run['max_cn'])
        if key in seen:
            dropped.append(run['path'])
            continue
        seen[key] = run
        kept.append(run)
    for path in dropped:
        logging.warning('%s: duplicate evaluation of the same (method, dataset, max-cn); dropped '
                        'in favour of the caller-side *_ploidy_eval run', os.path.basename(path))
    runs = kept
    # When several caps are shown at once, tell them apart in the column labels.
    caps_by_dataset = {}
    for run in runs:
        caps_by_dataset.setdefault(dataset_key(run), set()).add(run['max_cn'])
    for run in runs:
        key = dataset_key(run)
        tag = (F' [max-cn={fmt_max_cn(run["max_cn"])}]'
               if len(caps_by_dataset[key]) > 1 else '')
        run['dataset_key'] = key
        run['dataset_label'] = dataset_label(run, tag)
    return runs

runs = gather_runs()

frames = []
for run in runs:
    frames.append(pd.DataFrame({
        'method': run['tool'],
        'dataset': run['dataset_label'],
        'sample': run['sample'],
        'ploidy_error': run['ploidy_error'],
        'expected_ploidy': run['expected_ploidy'],
    }))
cells = pd.concat(frames, ignore_index=True)

# Row / column selection and ordering.
the_methods = order_methods({r['tool'] for r in runs})
the_datasets = []
for r in sorted(runs, key=lambda r: (r['dataset_key'][0], r['dataset_key'][1])
                                  + _avg_spot_len_sort(r) + (r['dataset_label'],)):
    if r['dataset_label'] not in the_datasets:
        the_datasets.append(r['dataset_label'])
if args.methods:
    unknown = [m for m in args.methods if m not in the_methods]
    if unknown:
        parser1.error(F'--methods: not among the evaluated tools {the_methods}: {unknown}')
    the_methods = args.methods
if args.datasets:
    the_datasets = [d for d in the_datasets if any(s in d for s in args.datasets)]
    if not the_datasets:
        parser1.error(F'--datasets: no dataset label contains any of {args.datasets}')
if (args.type & 0x1):  # testing: one method, one dataset
    the_methods = the_methods[:1]
    the_datasets = the_datasets[:1]
    logging.info('testing mode: only %s on %s', the_methods, the_datasets)
cells = cells[cells['method'].isin(the_methods) & cells['dataset'].isin(the_datasets)]
if cells.empty:
    sys.exit('ploidy-performances-eval: no per-cell ploidy error left after the selection. ')

window_by_entry = {(r['tool'], r['dataset_label']): r['window'] for r in runs}
n_finite_cells = int(np.isfinite(cells['ploidy_error'].to_numpy(dtype=float)).sum())
logging.info('%d method(s) x %d dataset(s); %d cells, %d with a finite ploidy error '
             '(%d failed to yield one)',
             len(the_methods), len(the_datasets), len(cells), n_finite_cells,
             len(cells) - n_finite_cells)

# Per-run statistics, in the spirit of the per-sample table ploidy_eval.py prints.
stat_rows = []
for run in runs:
    if run['tool'] not in the_methods or run['dataset_label'] not in the_datasets:
        continue
    errs = run['ploidy_error']
    finite = errs[np.isfinite(errs)]
    stat_rows.append({
        'method': run['tool'], 'dataset': run['dataset_label'],
        'max_cn': fmt_max_cn(run['max_cn']), 'window': run['window'],
        'n_samples': len(set(run['sample'])), 'n_cells': run['n_cells'],
        'n_cells_finite': len(finite),
        'mean_ploidy_error': float(np.mean(finite)),
        'mean_abs_ploidy_error': float(np.mean(np.abs(finite))),
        'pct_outliers': float(100.0 * np.mean(np.abs(finite) > run['window'])),
    })
if stat_rows:
    sys.stderr.write(pd.DataFrame(stat_rows).to_string(index=False) + '\n')

# ---------------------------------------------------------------------------------------------
# Y limits of the grid entries
# ---------------------------------------------------------------------------------------------

def entry_ylim(sub, windows):
    """Y limits for one grid entry (or for any set of entries): the range of the per-cell
    errors and of the 2x / 0.5x reference stars, padded, and always covering the
    expected-ploidy line (error = 0) and the ploidy window around it.  With --y-quantile,
    only the central fraction of the errors bounds the axis (the reference stars may then
    fall outside and be invisible)."""
    errs = sub['ploidy_error'].to_numpy(dtype=float)
    errs = errs[np.isfinite(errs)]
    stars = []
    for _, s in sub.groupby('sample', sort=False):
        exp = s['expected_ploidy'].to_numpy(dtype=float)
        exp = exp[np.isfinite(exp)]
        if len(exp):
            e = float(np.mean(exp))
            if e > 0:
                stars += [e, -0.5 * e]      # 2x / 0.5x the expected ploidy, in error space
    w = max(windows) if windows else DEFAULT_PLOIDY_WINDOW
    if args.y_quantile is not None and len(errs):
        lo = float(np.quantile(errs, (1.0 - args.y_quantile) / 2.0))
        hi = float(np.quantile(errs, 1.0 - (1.0 - args.y_quantile) / 2.0))
    elif len(errs):
        lo, hi = float(errs.min()), float(errs.max())
        if stars:
            lo, hi = min(lo, min(stars)), max(hi, max(stars))
    else:
        lo = min(stars + [-w]) if stars else -w
        hi = max(stars + [w]) if stars else w
    lo, hi = min(lo, -w, 0.0), max(hi, w, 0.0)
    pad = max(0.05 * (hi - lo), 0.1)
    return (lo - pad, hi + pad)

def modal_window(sub_runs):
    """The most common ploidy window of a set of runs (they normally all agree)."""
    ws = [r['window'] for r in sub_runs]
    return collections.Counter(ws).most_common(1)[0][0] if ws else DEFAULT_PLOIDY_WINDOW

def wrap_label(label, width=34):
    """Split a long dataset label at its ' | ' separators onto at most two lines."""
    if len(label) <= width:
        return label
    parts = label.split(' | ')
    if len(parts) < 2:
        return label
    half = (len(parts) + 1) // 2
    return ' | '.join(parts[:half]) + '\n' + ' | '.join(parts[half:])

def save_fig(fig, stem):
    """Write one figure as PDF and PNG, creating the output directory first."""
    out_dir = os.path.dirname(os.path.abspath(stem))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    fig.savefig(stem + '.pdf', dpi=args.dpi)
    fig.savefig(stem + '.png', dpi=args.dpi)

# ---------------------------------------------------------------------------------------------
# Main figure: the method (rows) x dataset (columns) grid of per-cell ploidy errors
# ---------------------------------------------------------------------------------------------

def draw_error_panel(ax, sub, samples, window, rng, show_x_labels):
    """One column-wise entry of the main grid: the per-cell ploidy-estimation errors of
    one method on one dataset -- one point per cell, one jittered error cluster per cell
    sample.  This is the column-wise entry of ploidy_eval.plot_ploidy (ploidy_eval.py
    line 515) translated into error space: the grey rectangle is the experimental ploidy
    window (expected +/- window, i.e. error in [-window, +window]), the crimson line is
    the expected ploidy itself (error = 0), and the blue stars mark 2x / 0.5x the
    expected ploidy, the characteristic whole-genome doubling / halving failures."""
    if sub.empty:
        ax.text(0.5, 0.5, 'no data', transform=ax.transAxes, ha='center', va='center',
                fontsize=8, color='0.45')
    for k, sample in enumerate(samples):
        ssub = sub[sub['sample'] == sample]
        errs = ssub['ploidy_error'].to_numpy(dtype=float)
        errs = errs[np.isfinite(errs)]
        exp = ssub['expected_ploidy'].to_numpy(dtype=float)
        exp = exp[np.isfinite(exp)]
        exp = float(np.mean(exp)) if len(exp) else float('nan')
        if not sub.empty:
            ax.add_patch(plt.Rectangle((k - 0.42, -window), 0.84, 2 * window,
                                       facecolor='0.82', edgecolor='none', zorder=1))
        if len(errs):  # one error per cell
            ax.scatter(k + rng.uniform(-args.jitter, args.jitter, size=len(errs)), errs,
                       s=5, alpha=0.45, color='#31688e', linewidths=0, zorder=2,
                       rasterized=True)
        if np.isfinite(exp) and exp > 0:
            for ref in (exp, -0.5 * exp):  # 2x / 0.5x the expected ploidy, in error space
                ax.scatter([k], [ref], marker='*', s=55, color='#1f77b4', zorder=3)
    if not sub.empty:
        ax.axhline(0.0, color='crimson', linewidth=1.5, zorder=4)  # expected: error = 0
    ax.set_xlim(-0.6, max(len(samples) - 0.4, 0.6))
    ax.set_xticks(range(len(samples)))
    if show_x_labels and samples:
        ax.set_xticklabels(samples, rotation=30, ha='right', fontsize=8)
    else:
        ax.set_xticklabels([])
    ax.grid(axis='y', linestyle=':', linewidth=0.6, alpha=0.6)
    ax.tick_params(axis='y', labelsize=8)

def plot_error_grid():
    """The main figure: rows = CNV calling methods, columns = datasets, and in every
    entry one error cluster per cell sample.  The row / column labelling follows
    plot_grid_main of scWGS-performances-eval.py (its line 153): the method name is the
    y-label of the first column, the dataset name the title of the first row, and the
    inner tick labels are hidden."""
    n_rows, n_cols = len(the_methods), len(the_datasets)
    samples_by_dataset = {ds: sorted(cells[cells['dataset'] == ds]['sample'].unique())
                          for ds in the_datasets}
    plotted_runs = [r for r in runs
                    if r['tool'] in the_methods and r['dataset_label'] in the_datasets]
    rng = np.random.default_rng(0)
    col_widths = [max(1.8, 1.15 * max(1, len(samples_by_dataset[ds])) + 1.0)
                  for ds in the_datasets]
    fig_w = max(6.0, 1.2 + sum(col_widths))
    fig_h = max(4.0, 1.6 + 2.5 * n_rows)
    fig1 = plt.figure(figsize=(fig_w, fig_h), constrained_layout=True)
    gs = gridspec.GridSpec(1 + n_rows, n_cols, height_ratios=[3] + [10] * n_rows,
                           figure=fig1, wspace=0, hspace=0.1)
    legend_ax = fig1.add_subplot(gs[0, :])
    legend_ax.set_axis_off()
    # Y limits: fixed by --ylim, else shared as --sharey asks, else per entry.
    ylim_all = tuple(args.ylim) if args.ylim else entry_ylim(cells, [r['window'] for r in plotted_runs])
    for rowidx, method in enumerate(the_methods):
        row_cells = cells[cells['method'] == method]
        row_runs = [r for r in plotted_runs if r['tool'] == method]
        ylim_row = (tuple(args.ylim) if args.ylim
                    else entry_ylim(row_cells, [r['window'] for r in row_runs]))
        for colidx, dataset in enumerate(the_datasets):
            ax2 = fig1.add_subplot(gs[rowidx + 1, colidx])
            sub = row_cells[row_cells['dataset'] == dataset]
            window = window_by_entry.get((method, dataset), DEFAULT_PLOIDY_WINDOW)
            draw_error_panel(ax2, sub, samples_by_dataset[dataset], window, rng,
                             show_x_labels=(rowidx == n_rows - 1))
            if args.ylim:
                ylim = tuple(args.ylim)
            elif args.sharey == 'all':
                ylim = ylim_all
            elif args.sharey == 'row':
                ylim = ylim_row
            else:
                ylim = entry_ylim(sub, [window])
            ax2.set_ylim(*ylim)
            # First column only: the method name (the row label).
            if colidx == 0:
                ax2.set_ylabel(caller2desc.get(method, method), fontsize=10, labelpad=8)
            else:
                ax2.set_ylabel('')
            same_ylim_within_row = bool(args.ylim) or args.sharey in ('all', 'row')
            ax2.tick_params(labelleft=(colidx == 0 or not same_ylim_within_row))
            # First row only: the dataset name (the column label).
            if rowidx == 0:
                ax2.set_title(wrap_label(dataset), fontsize=10)
    windows = sorted({r['window'] for r in plotted_runs})
    w_txt = F'+/- {windows[0]:g}' if len(windows) == 1 else 'varies by run'
    handles = [
        Patch(facecolor='0.82', edgecolor='none',
              label=F'experimental ploidy window ({w_txt})'),
        Line2D([], [], marker='o', linestyle='none', color='#31688e', alpha=0.45,
               markersize=5, label='per-cell ploidy error (one point per cell)'),
        Line2D([], [], color='crimson', linewidth=1.5, label='expected ploidy (error = 0)'),
        Line2D([], [], marker='*', linestyle='none', color='#1f77b4', markersize=9,
               label='2x / 0.5x expected ploidy (scaling-error references)'),
    ]
    legend_ax.legend(handles=handles, loc='center', fontsize=10, ncol=2, frameon=False)
    caps = ' & '.join(sorted({fmt_max_cn(r['max_cn']) for r in plotted_runs})) or '?'
    fig1.suptitle(F'Per-cell ploidy-estimation errors: {n_rows} method(s) (rows) x '
                  F'{n_cols} dataset(s) (columns), {n_finite_cells} cells'
                  F' | copy-number cap: {caps}', fontsize=13)
    fig1.supylabel('per-cell ploidy error (observed - expected)', fontsize=16)
    save_fig(fig1, args.output + '_ploidy_error_grid')
    plt.close(fig1)
    logging.info('wrote %s_ploidy_error_grid.pdf/.png', args.output)

# ---------------------------------------------------------------------------------------------
# Box-plot summaries of the same errors, in the style of scWGS-performances-eval.py
# ---------------------------------------------------------------------------------------------

def plot_main():
    """Pooled per-cell ploidy errors, one box per method (plot_main of
    scWGS-performances-eval.py, on the ploidy error instead of a CNV metric)."""
    fig1 = plt.figure(figsize=(1 * 7, 1 * 5), constrained_layout=True)
    ax = sns.boxplot(data=cells, x='method', y='ploidy_error', order=the_methods,
                     color='#457B9D', linewidth=1.0,
                     flierprops=dict(markersize=1.5, alpha=0.4))
    w = modal_window([r for r in runs
                      if r['tool'] in the_methods and r['dataset_label'] in the_datasets])
    ax.axhspan(-w, w, color='0.82', zorder=0)
    ax.axhline(0.0, color='crimson', linewidth=1.5, zorder=1)
    plt.tick_params(axis='both', which='major', labelsize=10)
    for label in ax.get_xticklabels():
        label.set_rotation(20)
        label.set_ha('right')
    ax.set_xlabel('Copy-number calling method')
    ax.set_ylabel('per-cell ploidy error (observed - expected)')
    ax.set_title(F'Per-cell ploidy-estimation errors pooled over {len(the_datasets)} '
                 F'dataset(s) (window +/- {w:g})', fontsize=10)
    save_fig(fig1, args.output + '_ploidy_error_main')
    plt.close(fig1)
    logging.info('wrote %s_ploidy_error_main.pdf/.png', args.output)

def plot_multirow_main():
    """One row per dataset, each a box plot of the per-cell ploidy errors of every method
    (plot_multirow_main of scWGS-performances-eval.py, with datasets for metrics)."""
    n = len(the_datasets)
    fig, axes = plt.subplots(nrows=n, ncols=1,
                             figsize=(max(7.0, 1.1 * len(the_methods) + 2.0), 1.9 * n + 1.5),
                             constrained_layout=True, squeeze=False)
    for rowidx, (ax, dataset) in enumerate(zip(axes[:, 0], the_datasets)):
        sub = cells[cells['dataset'] == dataset]
        sns.boxplot(data=sub, x='method', y='ploidy_error', order=the_methods, ax=ax,
                    color='steelblue', linewidth=1.2,
                    flierprops=dict(markersize=2, alpha=0.5))
        w = modal_window([r for r in runs if r['dataset_label'] == dataset
                          and r['tool'] in the_methods])
        ax.axhspan(-w, w, color='0.82', zorder=0)
        ax.axhline(0.0, color='crimson', linewidth=1.5, zorder=1)
        ax.set_title(wrap_label(dataset), fontsize=12, weight='bold')
        ax.set_ylabel('ploidy error', fontsize=10)
        ax.grid(axis='y', alpha=0.3)
        ax.tick_params(axis='y', labelsize=10)
        if rowidx == n - 1:
            ax.tick_params(axis='x', labelsize=10)
            for label in ax.get_xticklabels():
                label.set_rotation(15)
                label.set_ha('right')
            ax.set_xlabel('Copy-number calling method')
        else:
            ax.set_xticklabels([])
            ax.set_xlabel('')
    save_fig(fig, args.output + '_ploidy_error_multirow')
    plt.close(fig)
    logging.info('wrote %s_ploidy_error_multirow.pdf/.png', args.output)

# ---------------------------------------------------------------------------------------------
# Go
# ---------------------------------------------------------------------------------------------

plot_error_grid()
if (args.type & 0x2):
    logging.info('--type %d: only the main grid was requested', args.type)
    sys.exit(0)
plot_main()
plot_multirow_main()
sys.stderr.write(F'Wrote {args.output}_ploidy_error_grid.pdf/.png, '
                 F'{args.output}_ploidy_error_main.pdf/.png and '
                 F'{args.output}_ploidy_error_multirow.pdf/.png\n')
