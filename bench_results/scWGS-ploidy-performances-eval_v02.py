#!/usr/bin/env python
# A new script of https://github.com/zhaoxiaofei/copy-num-bench-scwgs
# Revised to emit main-text balloon plots in the style of Song et al.,
# Advanced Science 2025 (doi:10.1002/advs.202507839) Figure 4.

"""
Benchmark ploidy estimation across CNV calling methods and datasets.

What is plotted (main text)
---------------------------
Four balloon / dot-grid figures, one per biological group:

    1. COLO-829   germline-derived (emulated cell-line) data
    2. HCC1395    germline-derived (emulated cell-line) data
    3. HeLa       germline-derived (emulated cell-line) data
    4. ACT        real cancer-derived samples (TN1, TN2, ...)

In every figure:

    rows     = datasets
    columns  = methods  (each CNV caller contributes TWO columns:
                         copy-number cap at 10, and no cap)
    entry    = a filled circle whose SIZE and COLOUR encode the percentage
               of cells whose |observed − expected ploidy| is within the
               0.5 tolerance window.  A red cross marks a missing result
               (runtime error, empty per-cell table, or no finite ploidy).

Dataset identity
----------------
Germline-derived (S01, S02, 234HS, ...): a dataset is the combination of
average-spot-length, emulated cell-line, and original germline sample name.
The three cell-lines are split across three figures, so the row label inside
each figure is ``<donor> · <sampleType> · <avgSpotLen> bp``.

Real cancer-derived (ACT): a dataset is the original sample name (TN1, TN2,
... and the ACT cell-line samples).  There is a single ACT figure; one
ploidy-eval summary that covers many samples is split on the per-cell
``sample`` column.

Input
-----
The main input files are the ploidy-evaluation summaries written by
ploidy_eval.py and ploid_tools.py:

    *_ploidy_eval_summary.json
    *_ploidy_tool_eval_summary.json

each sitting next to its ``<prefix>_percell.tsv`` sibling.  Runs that have a
summary but no usable per-cell table are kept and drawn as red crosses.

Usage
-----
    python bench_results/scWGS-ploidy-performances-eval_v01.py \\
        -i '../data/*/4from2_*_ploidy_eval_summary.json' \\
           '../data/*/4from3_*_ploidy_eval_summary.json' \\
           '../data/*/4from2_*_ploidy_tool_eval_summary.json' \\
        -o bench_results/ploidy-performances

    # style preview with in-script synthetic data (no input files needed):
    python bench_results/scWGS-ploidy-performances-eval_v01.py --demo \\
        -o bench_results/ploidy-performances-demo

Outputs (PDF + PNG at --dpi):

    <output>_main_COLO-829.{pdf,png}
    <output>_main_HCC1395.{pdf,png}
    <output>_main_HeLa.{pdf,png}
    <output>_main_ACT.{pdf,png}
    <output>_main_four.{pdf,png}          (the four panels combined)
    <output>_pct_within_long.tsv          (the numbers behind the dots)

Pass ``--legacy`` to also write the older per-cell error-grid / box-plot
figures (``_ploidy_error_grid``, ``_ploidy_error_main``,
``_ploidy_error_multirow``).
"""

from __future__ import annotations

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
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.cm import ScalarMappable
import seaborn as sns

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(filename)s %(levelname)s %(message)s')

# ---------------------------------------------------------------------------------------------
# Constants (kept in step with ploidy_eval.py so the two scripts cannot drift apart)
# ---------------------------------------------------------------------------------------------

SUMMARY_SUFFIX = '_summary.json'
PERCELL_SUFFIX = '_percell.tsv'
DEFAULT_PLOIDY_WINDOW = 0.5
DEFAULT_MAX_CN = 10.0

CALLER_ORDER = ['scabsolute', 'hmmcopy', 'ginkgo', 'copynumber', 'secnv',
                'sccnv', 'scyn', 'chisel', 'aneufinder', 'flcna']

TOOL_PRETTY = {
    'aneufinder': 'AneuFinder',
    'flcna': 'FLCNA',
    'chisel': 'CHISEL',
    'copynumber': 'CopyNumber',
    'ginkgo': 'Ginkgo',
    'hmmcopy': 'HMMcopy',
    'secnv': 'SeCNV',
    'sccnv': 'SCCNV',
    'scyn': 'SCYN',
    'scabsolute': 'scAbsolute',
}

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
# Annotated row labels, off by default exactly as in scWGS-performances-eval.py.
caller2desc = {}

# The three emulated cell-lines of the germline-derived (4from3) arm, and the
# spellings that show up in file names / JSON / the per-cell `sample` column.
GERMLINE_CELL_LINE_ORDER = ['COLO-829', 'HCC1395', 'HeLa']
_CELL_LINE_ALIASES = {
    'colo829': 'COLO-829', 'colo-829': 'COLO-829', 'colo_829': 'COLO-829',
    'colo 829': 'COLO-829', 'col0829': 'COLO-829',
    'hcc1395': 'HCC1395', 'hcc-1395': 'HCC1395', 'hcc_1395': 'HCC1395',
    'hela': 'HeLa', 'hela-s3': 'HeLa', 'helas3': 'HeLa', 'hela_s3': 'HeLa',
}

# ACT (Minussi 2021, PRJNA629885) sample ids, used to recognise and order rows
# of the cancer-derived figure.  Unknown ACT-like names are appended after these.
ACT_SAMPLE_ORDER = [
    'TN1', 'TN2', 'TN3', 'TN4', 'TN5', 'TN6', 'TN7', 'TN8',
    'MDAMB231c28', 'MDAMB231c8', 'MDAMB231_popp31', 'mb157', 'BT20', 'mb453',
]
_ACT_SAMPLE_ALIASES = {
    'tn1': 'TN1', 'tn-1': 'TN1', 'tn01': 'TN1',
    'tn2': 'TN2', 'tn-2': 'TN2', 'tn02': 'TN2',
    'tn3': 'TN3', 'tn-3': 'TN3', 'tn03': 'TN3',
    'tn4': 'TN4', 'tn-4': 'TN4', 'tn04': 'TN4',
    'tn5': 'TN5', 'tn-5': 'TN5', 'tn05': 'TN5',
    'tn6': 'TN6', 'tn-6': 'TN6', 'tn06': 'TN6',
    'tn7': 'TN7', 'tn-7': 'TN7', 'tn07': 'TN7',
    'tn8': 'TN8', 'tn-8': 'TN8', 'tn08': 'TN8',
    'mdamb231c28': 'MDAMB231c28', 'mdamb231ex1': 'MDAMB231c28',
    'mdamb231c8': 'MDAMB231c8', 'mdamb231ex2': 'MDAMB231c8',
    'mdamb231popp31': 'MDAMB231_popp31', 'mdamb231': 'MDAMB231_popp31',
    'mdamb231parental': 'MDAMB231_popp31', 'mb231': 'MDAMB231_popp31',
    'mb-231': 'MDAMB231_popp31',
    'mb157': 'mb157', 'mdamb157': 'mb157', 'mb-157': 'mb157',
    'bt20': 'BT20', 'bt-20': 'BT20', 'bt_20': 'BT20',
    'mb453': 'mb453', 'mdamb453': 'mb453', 'mb-453': 'mb453',
}

# Song et al. Adv Sci 2025 Fig. 4 sequential scale: cream -> gold -> brown -> near-black.
BALLOON_CMAP = LinearSegmentedColormap.from_list('song_ylorbr', [
    '#fff7bc', '#fee391', '#fec44f', '#fe9929',
    '#ec7014', '#cc4c02', '#993404', '#662506', '#3d1500',
])

# Fallback parser for pipeline file names, used only when the summary JSON
# carries no donor / sampleType / avgSpotLen / tool / cellLine:
#   4from2_2_<donor>_3_<sampleType>_<avgSpotLen>_4_step<n>_<tool>_ploidy_eval[_maxcn_<v>]
#   4from3_2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>_4_step<n>_<tool>_ploidy_eval[...]
STEM_RE = re.compile(
    r'^(?:(?P<source>4from[23])_2_)?'
    r'(?P<donor>.+?)_3_(?P<sampleType>.+?)_(?P<avgSpotLen>\d{2,})'
    r'(?:_(?P<cellLine>.+?))?_4_step\d+_(?P<tool>.+?)_ploidy(?:_tool)?_eval'
    r'(?:_maxcn_(?P<maxcn>[^_]+))?$')

PLOT_SPECS = [
    # (plot_id, title, kind)  kind is 'germline' or 'act'
    ('COLO-829', 'COLO-829  (germline-derived, emulated cell-line)', 'germline'),
    ('HCC1395',  'HCC1395  (germline-derived, emulated cell-line)',  'germline'),
    ('HeLa',     'HeLa  (germline-derived, emulated cell-line)',     'germline'),
    ('ACT',      'ACT  (real cancer-derived samples)',               'act'),
]


# ---------------------------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------------------------

def build_parser():
    p = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__)
    p.add_argument('-i', '--input', nargs='+', default=None,
                   help='Ploidy-evaluation summary files or globs (typically '
                        '*_ploidy*_eval*_summary.json); each must sit next to its '
                        '<prefix>_percell.tsv sibling.  Read from stdin instead when '
                        'omitted or "-". ')
    p.add_argument('-t', '--type', type=int, default=0,
                   help='Output type. 0: the four main-text balloon plots. '
                        '1: testing (first few rows / columns). '
                        '2: only the combined 2x2 figure. ')
    p.add_argument('-o', '--output', default='scWGS-ploidy-performances')
    p.add_argument('--max-cn', default='all',
                   help='Which copy-number cap(s) to plot.  Default "all" draws both '
                        'the cap-at-10 and the uncapped method of every caller.  A '
                        'number (e.g. 10) or "inf" restricts the figure to that cap.')
    p.add_argument('--methods', nargs='+', default=None,
                   help='Restrict (and order) the callers to these tools, by exact name. ')
    p.add_argument('--datasets', nargs='+', default=None,
                   help='Restrict (and order) the rows to the datasets whose label '
                        'contains any of these substrings. ')
    p.add_argument('--window', type=float, default=None,
                   help='Ploidy-error tolerance used for the balloon '
                        '(default: each run\'s own ploidy_window, typically 0.5).')
    p.add_argument('--dpi', type=int, default=300,
                   help='Raster resolution of the PNG outputs. ')
    p.add_argument('--legacy', action='store_true',
                   help='Also write the older per-cell error-grid and box-plot figures.')
    p.add_argument('--sharey', default='all', choices=['all', 'row', 'none'],
                   help='(legacy) Share the y scale of the error-grid entries.')
    p.add_argument('--y-quantile', type=float, default=None,
                   help='(legacy) Clip the error-grid y axis to this central fraction.')
    p.add_argument('--ylim', type=float, nargs=2, default=None, metavar=('LOW', 'HIGH'),
                   help='(legacy) Fixed y limits of the error grid.')
    p.add_argument('--jitter', type=float, default=0.28,
                   help='(legacy) Half-width of the per-cell error-cluster jitter.')
    p.add_argument('--demo', action='store_true',
                   help='Ignore -i and synthesise a realistic-looking table so the '
                        'four main-text figures can be previewed without benchmark files.')
    p.add_argument('--dot-scale', type=float, default=1.0,
                   help='Multiply every balloon area by this factor (default 1).')
    return p


# ---------------------------------------------------------------------------------------------
# Small helpers
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


def cap_tick_label(value):
    v = norm_max_cn(value)
    return 'no cap' if np.isinf(v) else '≤10' if abs(v - 10.0) < 1e-9 else F'≤{v:g}'


def pretty_tool(tool):
    return TOOL_PRETTY.get(tool, tool)


def _fold(name):
    return re.sub(r'[^a-z0-9]', '', str(name).strip().lower())


def canon_cell_line(name):
    if name is None:
        return None
    s = str(name).strip()
    if not s or s.lower() in ('tumor', 'na', 'nan', 'none', 'unknown'):
        return None
    folded = _fold(s)
    if folded in _CELL_LINE_ALIASES:
        return _CELL_LINE_ALIASES[folded]
    # already a canonical name?
    for canon in GERMLINE_CELL_LINE_ORDER:
        if _fold(canon) == folded:
            return canon
    return None


def canon_act_sample(name):
    if name is None:
        return None
    s = str(name).strip()
    if not s:
        return None
    folded = _fold(s)
    if folded in _ACT_SAMPLE_ALIASES:
        return _ACT_SAMPLE_ALIASES[folded]
    m = re.match(r'^tn0*(\d+)$', folded)
    if m:
        return F'TN{int(m.group(1))}'
    return s


def is_act_sample_name(name):
    if name is None:
        return False
    folded = _fold(name)
    if folded in _ACT_SAMPLE_ALIASES:
        return True
    if re.match(r'^tn\d+$', folded):
        return True
    return False


def is_germline_donor(donor):
    d = str(donor or '').strip()
    if re.match(r'^S\d+$', d, re.I):
        return True
    if re.search(r'\d+HS', d, re.I):
        return True
    return False


def is_caller_eval_stem(stem):
    """True for the *_ploidy_eval summaries of ploidy_eval.py (CNV callers), False for
    the *_ploidy_tool_eval summaries of ploidy_tools.py."""
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


def order_methods(tools):
    known = [t for t in CALLER_ORDER if t in tools]
    rest = sorted(set(tools) - set(CALLER_ORDER))
    return known + rest


def wrap_label(label, width=34):
    """Split a long dataset label at its ' | ' / ' · ' separators onto at most two lines."""
    if len(label) <= width:
        return label
    for sep in (' · ', ' | '):
        parts = label.split(sep)
        if len(parts) >= 2:
            half = (len(parts) + 1) // 2
            return sep.join(parts[:half]) + '\n' + sep.join(parts[half:])
    return label


def save_fig(fig, stem, dpi):
    """Write one figure as PDF and PNG, creating the output directory first."""
    out_dir = os.path.dirname(os.path.abspath(stem))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    fig.savefig(stem + '.pdf', dpi=dpi, bbox_inches='tight')
    fig.savefig(stem + '.png', dpi=dpi, bbox_inches='tight')


def _first_nonempty(*values):
    for v in values:
        if v is None:
            continue
        s = str(v).strip()
        if s and s.lower() not in ('nan', 'none'):
            return s
    return ''


def _column_mode(df, col):
    if df is None or col not in df.columns:
        return ''
    vals = [str(v).strip() for v in df[col].dropna().astype(str) if str(v).strip() != '']
    return vals[0] if vals else ''


# ---------------------------------------------------------------------------------------------
# Input: the *_ploidy*_eval*_summary.json files and their <prefix>_percell.tsv siblings
# ---------------------------------------------------------------------------------------------

def load_run(path):
    """One ploidy-evaluation run.  Returns None only when the file cannot even be
    identified (wrong name / unreadable JSON).  A run with no usable per-cell
    table is returned with failed=True so it can be drawn as a red cross."""
    base = os.path.basename(path)
    if not base.endswith(SUMMARY_SUFFIX):
        logging.warning('%s: not a *%s file; skipping', path, SUMMARY_SUFFIX)
        return None
    stem = base[:-len(SUMMARY_SUFFIX)]
    if not ('ploidy' in stem and 'eval' in stem):
        logging.warning('%s: does not follow the *_ploidy*_eval*%s naming; skipping',
                        path, SUMMARY_SUFFIX)
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

    m = STEM_RE.match(stem)
    mg = m.groupdict() if m else {}
    percell_path = path[:-len(SUMMARY_SUFFIX)] + PERCELL_SUFFIX
    df = None
    failed_reason = None
    if not (os.path.isfile(percell_path) and os.path.getsize(percell_path) > 0):
        failed_reason = 'no usable sibling per-cell table'
    else:
        try:
            df = pd.read_csv(percell_path, sep='\t')
        except (OSError, ValueError) as exc:
            failed_reason = F'unreadable per-cell table ({exc})'
            df = None
        if df is not None:
            missing = [c for c in ('sample', 'ploidy_error') if c not in df.columns]
            if missing:
                failed_reason = F'per-cell table lacks the column(s) {missing}'
                df = None

    tool = _first_nonempty(js.get('tool'),
                           _column_mode(df, 'tool'),
                           mg.get('tool'))
    if not tool:
        logging.warning('%s: no tool identity; skipping', path)
        return None

    donor = _first_nonempty(js.get('donor'),
                            _column_mode(df, 'donor'),
                            mg.get('donor'))
    sample_type = _first_nonempty(js.get('sampleType') or js.get('sample_type'),
                                  _column_mode(df, 'sampleType'),
                                  _column_mode(df, 'sample_type'),
                                  mg.get('sampleType'))
    avg_spot_len = _first_nonempty(js.get('avgSpotLen') or js.get('avg_spot_len'),
                                   _column_mode(df, 'avgSpotLen'),
                                   _column_mode(df, 'avg_spot_len'),
                                   mg.get('avgSpotLen'))
    cell_line = _first_nonempty(js.get('cellLine') or js.get('cell_line'),
                                _column_mode(df, 'cellLine'),
                                _column_mode(df, 'cell_line'),
                                mg.get('cellLine'))
    source = _first_nonempty(mg.get('source'),
                             '4from3' if cell_line and canon_cell_line(cell_line) else '',
                             '4from2')

    if df is not None:
        sample = df['sample'].astype(str).to_numpy()
        ploidy_error = pd.to_numeric(df['ploidy_error'], errors='coerce').to_numpy(dtype=float)
        expected = (pd.to_numeric(df['expected_ploidy'], errors='coerce').to_numpy(dtype=float)
                    if 'expected_ploidy' in df.columns
                    else np.full(len(df), np.nan))
        if 'is_outlier' in df.columns:
            is_outlier = df['is_outlier'].to_numpy()
        else:
            is_outlier = None
    else:
        sample = np.array([], dtype=object)
        ploidy_error = np.array([], dtype=float)
        expected = np.array([], dtype=float)
        is_outlier = None

    n_finite = int(np.isfinite(ploidy_error).sum()) if len(ploidy_error) else 0
    if failed_reason is None and n_finite == 0:
        failed_reason = 'no cell has a finite ploidy error'

    window = js.get('ploidy_window', DEFAULT_PLOIDY_WINDOW)
    try:
        window = float(window)
    except (TypeError, ValueError):
        window = DEFAULT_PLOIDY_WINDOW

    if failed_reason:
        logging.warning('%s: %s; kept as a missing (red-cross) result',
                        os.path.basename(path), failed_reason)

    return {
        'path': path, 'percell_path': percell_path, 'stem': stem,
        'tool': tool, 'donor': donor, 'sampleType': sample_type,
        'avgSpotLen': avg_spot_len, 'cellLine': cell_line, 'source': source,
        'max_cn': norm_max_cn(js.get('max_cn', DEFAULT_MAX_CN)), 'window': window,
        'sample': sample, 'ploidy_error': ploidy_error, 'expected_ploidy': expected,
        'is_outlier': is_outlier,
        'n_cells': int(len(sample)), 'n_cells_finite': n_finite,
        'failed': bool(failed_reason),
        'caller_eval': is_caller_eval_stem(stem),
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


def gather_runs(args):
    patterns = list(args.input or [])
    if not patterns or patterns == ['-']:
        if sys.stdin.isatty():
            parser_error = getattr(gather_runs, 'parser_error', None)
            msg = ('no input: pass *_ploidy*_eval*_summary.json files or globs with -i, '
                   'or feed them on stdin (or pass --demo)')
            if parser_error:
                parser_error(msg)
            sys.exit('ploidy-performances-eval: ' + msg)
        patterns = [line.rstrip('\n') for line in sys.stdin]
    if not patterns:
        sys.exit('ploidy-performances-eval: no input: pass files with -i, or --demo')
    paths = expand_inputs(patterns)
    if not paths:
        sys.exit('ploidy-performances-eval: no existing input file among the given patterns')
    runs = []
    for path in paths:
        run = load_run(path)
        if run is not None:
            runs.append(run)
            logging.info('loaded %s: tool=%s donor=%s sampleType=%s avgSpotLen=%s '
                         'cellLine=%s max-cn=%s window=%g cells=%d (finite: %d)%s',
                         os.path.basename(path), run['tool'], run['donor'],
                         run['sampleType'], run['avgSpotLen'], run['cellLine'] or '-',
                         fmt_max_cn(run['max_cn']), run['window'], run['n_cells'],
                         run['n_cells_finite'],
                         ' [FAILED]' if run['failed'] else '')
    if not runs:
        sys.exit('ploidy-performances-eval: no usable *_ploidy*_eval*_summary.json input. ')

    wanted = str(args.max_cn).strip().lower()
    if wanted in ('all', 'any', 'both'):
        wanted = 'all'
    else:
        try:
            wanted = float(wanted) if wanted not in ('inf', 'infinity') else float('inf')
            if not (wanted > 0):
                raise ValueError
        except ValueError:
            sys.exit(F'ploidy-performances-eval: --max-cn must be a positive number, inf, '
                     F'or all (got {args.max_cn})')
        before = len(runs)
        runs = [r for r in runs if r['max_cn'] == wanted]
        logging.info('--max-cn %s: kept %d of %d runs (the rest were evaluated at other caps)',
                     fmt_max_cn(wanted), len(runs), before)
        if not runs:
            sys.exit(F'ploidy-performances-eval: every input run is at a copy-number cap other '
                     F'than {fmt_max_cn(wanted)}; retry with --max-cn all (or inf). ')

    # De-duplicate (tool, dataset, max_cn, cellLine): on real tumors a caller named in
    # --ploidy-tools is evaluated both by ploidy_eval.py and by ploidy_tools.py; keep
    # the caller-side *_ploidy_eval run and drop the duplicate.
    runs.sort(key=lambda r: (not r['caller_eval'], r['path']))
    seen, kept, dropped = {}, [], []
    for run in runs:
        key = (run['tool'], dataset_key(run), run['max_cn'],
               canon_cell_line(run['cellLine']) or run['cellLine'] or '')
        if key in seen:
            dropped.append(run['path'])
            continue
        seen[key] = run
        kept.append(run)
    for path in dropped:
        logging.warning('%s: duplicate evaluation of the same (method, dataset, max-cn); dropped '
                        'in favour of the caller-side *_ploidy_eval run', os.path.basename(path))
    runs = kept
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


# ---------------------------------------------------------------------------------------------
# Balloon-plot table: one row per (plot, dataset, tool, max_cn)
# ---------------------------------------------------------------------------------------------

def germline_row_label(run):
    """Row label inside a cell-line figure: donor · sampleType · avgSpotLen.
    sampleType is omitted when it is empty or the generic 'tumor' tag of data_tumor.py."""
    parts = []
    if run['donor']:
        parts.append(str(run['donor']))
    st = str(run['sampleType'] or '').strip()
    if st and st.lower() not in ('tumor', 'na', 'nan'):
        # drop a trailing _ILLUMINA-style platform suffix that data_tumor.py appends
        st = re.sub(r'_(ILLUMINA|PACBIO|ONT|BGISEQ|MGI)$', '', st, flags=re.I)
        if st and st.lower() not in ('tumor',):
            parts.append(st)
    if run['avgSpotLen']:
        parts.append(F"{run['avgSpotLen']} bp")
    return ' · '.join(parts) if parts else (run['stem_label'] or 'dataset')


def _pct_within(err, window, is_outlier=None):
    """Percentage of finite cells inside the ploidy window.  NaN if none are finite."""
    err = np.asarray(err, dtype=float)
    finite = np.isfinite(err)
    if not finite.any():
        return float('nan'), 0
    if is_outlier is not None and len(is_outlier) == len(err):
        # honour ploidy_eval.py's own flag when the window was not overridden
        try:
            out = np.asarray(is_outlier, dtype=bool)[finite]
            pct = float(100.0 * np.mean(~out))
            return pct, int(finite.sum())
        except (TypeError, ValueError):
            pass
    pct = float(100.0 * np.mean(np.abs(err[finite]) <= window))
    return pct, int(finite.sum())


def _entry(plot_id, dataset, run, err, n_cells, is_outlier, window, failed=False):
    if failed or n_cells == 0:
        pct, n_fin = float('nan'), 0
        mean_abs = float('nan')
    else:
        pct, n_fin = _pct_within(err, window, is_outlier=is_outlier)
        finite = np.asarray(err, dtype=float)
        finite = finite[np.isfinite(finite)]
        mean_abs = float(np.mean(np.abs(finite))) if len(finite) else float('nan')
        if not np.isfinite(pct):
            failed = True
    return {
        'plot': plot_id,
        'dataset': dataset,
        'tool': run['tool'],
        'max_cn': run['max_cn'],
        'method': F"{run['tool']}|{fmt_max_cn(run['max_cn'])}",
        'window': window,
        'n_cells': int(n_cells),
        'n_cells_finite': int(n_fin),
        'pct_within': pct,
        'mean_abs_ploidy_error': mean_abs,
        'failed': bool(failed or not np.isfinite(pct)),
        'donor': run.get('donor', ''),
        'sampleType': run.get('sampleType', ''),
        'avgSpotLen': run.get('avgSpotLen', ''),
        'cellLine': run.get('cellLine', ''),
    }


def expand_runs_to_entries(runs, window_override=None):
    """Turn raw runs into one balloon-plot entry per (plot, dataset, tool, max_cn).

    Two germline layouts are accepted:

      A. cellLine is in the file name / JSON (4from3_..._<cellLine>_...).  The
         whole run is one dataset of that cell-line plot.
      B. cellLine is the per-cell ``sample`` label (the original v01 convention).
         The run is split across the three cell-line plots.

    ACT / 4from2 runs are always split on the per-cell ``sample`` column.
    """
    entries = []
    for run in runs:
        window = (float(window_override) if window_override is not None else run['window'])
        cl_from_name = canon_cell_line(run.get('cellLine'))
        samples = [s for s in pd.unique(run['sample']) if str(s).strip() not in ('', 'nan')]
        sample_cls = {s: canon_cell_line(s) for s in samples}
        n_cell_line_samples = sum(1 for s in samples if sample_cls[s])
        n_act_samples = sum(1 for s in samples if is_act_sample_name(s))

        def mask_of(sample):
            return np.asarray(run['sample'], dtype=str) == str(sample)

        def outlier_of(mask):
            if run['is_outlier'] is None:
                return None
            return np.asarray(run['is_outlier'])[mask]

        # Failed run with no cells: still register the method (and, when we know
        # the cell-line / dataset, a red-cross entry at that row).
        if run['failed'] and len(samples) == 0:
            if cl_from_name:
                entries.append(_entry(cl_from_name, germline_row_label(run), run,
                                      [], 0, None, window, failed=True))
            elif is_germline_donor(run['donor']):
                for cl in GERMLINE_CELL_LINE_ORDER:
                    entries.append(_entry(cl, germline_row_label(run), run,
                                          [], 0, None, window, failed=True))
            else:
                entries.append(_entry('ACT', None, run, [], 0, None, window, failed=True))
            continue

        # Layout A: filename/JSON names a germline cell-line.
        if cl_from_name:
            entries.append(_entry(
                cl_from_name, germline_row_label(run), run,
                run['ploidy_error'], run['n_cells'], run['is_outlier'], window,
                failed=run['failed']))
            continue

        # Layout B: per-cell sample labels ARE the emulated cell-lines.
        if n_cell_line_samples and n_cell_line_samples >= n_act_samples:
            for sample in samples:
                cl = sample_cls.get(sample)
                if not cl:
                    continue
                mask = mask_of(sample)
                entries.append(_entry(
                    cl, germline_row_label(run), run,
                    run['ploidy_error'][mask], int(mask.sum()), outlier_of(mask),
                    window, failed=run['failed']))
            continue

        # ACT / real cancer: one row per original sample name.
        if samples:
            for sample in samples:
                mask = mask_of(sample)
                entries.append(_entry(
                    'ACT', canon_act_sample(sample), run,
                    run['ploidy_error'][mask], int(mask.sum()), outlier_of(mask),
                    window, failed=run['failed']))
        else:
            # no sample column values: one ACT row labelled by donor
            label = run['donor'] or run['stem_label'] or 'ACT'
            entries.append(_entry(
                'ACT', canon_act_sample(label), run,
                run['ploidy_error'], run['n_cells'], run['is_outlier'],
                window, failed=run['failed']))
    return pd.DataFrame(entries)


def _donor_sort_key(donor):
    d = str(donor or '')
    m = re.match(r'^S(\d+)$', d, re.I)
    if m:
        return (0, int(m.group(1)), d)
    m = re.match(r'^(\d+)HS', d, re.I)
    if m:
        return (1, int(m.group(1)), d)
    return (2, 0, d.lower())


def _spot_sort_key(spot):
    try:
        return (0, int(float(spot)))
    except (TypeError, ValueError):
        return (1, 0)


def order_germline_rows(sub):
    """Stable order: donor (S01, S02, 234HS, ...), then avgSpotLen, then sampleType."""
    rows = []
    seen = set()
    recs = (sub[['dataset', 'donor', 'sampleType', 'avgSpotLen']]
            .drop_duplicates().to_dict('records'))
    recs.sort(key=lambda r: (_donor_sort_key(r['donor']),
                             _spot_sort_key(r['avgSpotLen']),
                             str(r['sampleType'] or ''),
                             str(r['dataset'])))
    for r in recs:
        if r['dataset'] and r['dataset'] not in seen:
            seen.add(r['dataset'])
            rows.append(r['dataset'])
    return rows


def order_act_rows(sub):
    names = [d for d in sub['dataset'].dropna().unique() if d]
    rank = {n: i for i, n in enumerate(ACT_SAMPLE_ORDER)}

    def key(n):
        if n in rank:
            return (0, rank[n], n)
        m = re.match(r'^TN(\d+)$', str(n), re.I)
        if m:
            return (0, 1000 + int(m.group(1)), n)
        return (1, 0, str(n).lower())

    return sorted(names, key=key)


def method_columns(tools, max_cn_mode):
    """List of (tool, max_cn) pairs that become the x-axis.

    Default (max_cn_mode == 'all'): every caller contributes the cap-at-10
    column and the uncapped column, in that order, even if one of the two is
    entirely missing (those cells become red crosses).
    """
    tools = order_methods(tools)
    if max_cn_mode == 'all':
        cols = []
        for t in tools:
            cols.append((t, 10.0))
            cols.append((t, float('inf')))
        return cols
    cap = norm_max_cn(max_cn_mode)
    return [(t, cap) for t in tools]


# ---------------------------------------------------------------------------------------------
# Balloon / dot-grid drawing (Song et al. Adv Sci 2025 Fig. 4)
# ---------------------------------------------------------------------------------------------

def _apply_style():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'DejaVu Sans', 'Helvetica', 'Liberation Sans'],
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
        'axes.linewidth': 0.8,
        'xtick.major.size': 0,
        'ytick.major.size': 0,
        'axes.spines.top': True,
        'axes.spines.right': True,
    })


def _dot_size(pct, s_min, s_max):
    """Circle area.  0% is a small pale dot; 100% is the full marker."""
    if not np.isfinite(pct):
        return 0.0
    p = min(max(float(pct), 0.0), 100.0) / 100.0
    return s_min + (s_max - s_min) * p


def draw_balloon_ax(ax, matrix, row_labels, col_pairs, cmap, norm,
                    s_min=28.0, s_max=420.0, show_xlabel=True,
                    tick_fontsize=7.5, ytick_fontsize=8.0):
    """One Song-style balloon panel.

    matrix[i, j] is pct_within in [0, 100], or NaN for a red cross.
    col_pairs is a list of (tool, max_cn).  X-axis is grouped: each caller
    occupies two ticks (≤10, no cap) with the caller name centred underneath.
    """
    n_rows, n_cols = matrix.shape
    for i in range(n_rows):
        for j in range(n_cols):
            val = matrix[i, j]
            if not np.isfinite(val):
                ax.plot(j, i, marker='x', color='#d62728', markersize=9.5,
                        markeredgewidth=1.8, linestyle='none', zorder=4,
                        clip_on=False)
            else:
                ax.scatter([j], [i],
                           s=_dot_size(val, s_min, s_max),
                           c=[cmap(norm(val / 100.0))],
                           edgecolors='#3d1500', linewidths=0.35,
                           zorder=3, clip_on=False)

    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(n_rows - 0.5, -0.5)          # first dataset at the top
    ax.set_xticks(range(n_cols))
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels([wrap_label(r, width=42) for r in row_labels],
                       fontsize=ytick_fontsize)
    if show_xlabel:
        ax.set_xticklabels(
            [F'{pretty_tool(t)}  {cap_tick_label(c)}' for t, c in col_pairs],
            fontsize=tick_fontsize, rotation=90, ha='right', va='center',
            rotation_mode='anchor')
    else:
        ax.set_xticklabels([])

    # Light cell grid, plus a slightly stronger vertical rule between callers.
    ax.set_xticks(np.arange(-0.5, n_cols, 1.0), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1.0), minor=True)
    ax.grid(which='minor', color='#e6e6e6', linestyle='-', linewidth=0.6, zorder=0)
    ax.tick_params(which='minor', bottom=False, left=False)
    ax.set_axisbelow(True)
    for j in range(1, n_cols):
        if col_pairs[j][0] != col_pairs[j - 1][0]:
            ax.axvline(j - 0.5, color='#8a8a8a', linewidth=0.7, zorder=1)
    for spine in ax.spines.values():
        spine.set_color('#444444')
        spine.set_linewidth(0.8)
    ax.set_facecolor('white')


def _add_colorbar(fig, cax, cmap, norm, label):
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label(label, fontsize=8.5)
    cb.ax.tick_params(labelsize=7.5, length=2.5, width=0.6)
    cb.outline.set_linewidth(0.6)
    return cb


def _add_size_legend(ax, s_min, s_max, values=(25, 50, 75, 100)):
    """Vertical size legend, drawn in its own axes to the right of the colourbar."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    ax.text(0.5, 0.96, '% cells', ha='center', va='top', fontsize=8)
    ys = np.linspace(0.82, 0.12, len(values))
    for y, v in zip(ys, values):
        ax.scatter([0.28], [y],
                   s=_dot_size(v, s_min, s_max),
                   facecolors='#b35806', edgecolors='#3d1500',
                   linewidths=0.35, clip_on=False, zorder=3)
        ax.text(0.52, y, F'{v:g}', ha='left', va='center', fontsize=8)


def plot_one_balloon(entries, plot_id, title, row_labels, col_pairs, args,
                     panel_letter=None, s_min=28.0, s_max=420.0):
    """Build, save and return one main-text balloon figure."""
    sub = entries[entries['plot'] == plot_id]
    n_rows, n_cols = len(row_labels), len(col_pairs)
    if n_rows == 0 or n_cols == 0:
        logging.warning('plot %s: nothing to draw (rows=%d cols=%d)', plot_id, n_rows, n_cols)
        return None

    lookup = {}
    for rec in sub.itertuples(index=False):
        if rec.dataset is None or (isinstance(rec.dataset, float) and np.isnan(rec.dataset)):
            continue
        lookup[(rec.dataset, rec.tool, rec.max_cn)] = rec

    matrix = np.full((n_rows, n_cols), np.nan, dtype=float)
    for i, ds in enumerate(row_labels):
        for j, (tool, cap) in enumerate(col_pairs):
            rec = lookup.get((ds, tool, cap))
            if rec is None or rec.failed or not np.isfinite(rec.pct_within):
                matrix[i, j] = np.nan
            else:
                matrix[i, j] = rec.pct_within

    s_min *= args.dot_scale
    s_max *= args.dot_scale
    cmap = BALLOON_CMAP
    norm = Normalize(vmin=0.0, vmax=1.0)

    fig_w = max(9.0, 0.70 * n_cols + 4.2)
    fig_h = max(4.6, 0.42 * n_rows + 3.2)
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = gridspec.GridSpec(
        1, 3, figure=fig,
        width_ratios=[1.0, 0.028, 0.10],
        wspace=0.10,
        left=0.16, right=0.98, top=0.88, bottom=0.28)
    ax = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[0, 1])
    lax = fig.add_subplot(gs[0, 2])

    draw_balloon_ax(ax, matrix, row_labels, col_pairs, cmap, norm,
                    s_min=s_min, s_max=s_max, show_xlabel=True)
    _add_colorbar(fig, cax, cmap, norm, '% cells within ±0.5 of truth')
    _add_size_legend(lax, s_min, s_max)

    letter = F'{panel_letter}  ' if panel_letter else ''
    ax.set_title(letter + title, fontsize=12, loc='left', pad=10, fontweight='bold')
    ax.set_ylabel('Dataset', fontsize=10, labelpad=8)
    ax.tick_params(axis='x', pad=2)

    fig.text(0.16, 0.02,
             'Each caller is shown twice: copy-number cap at 10 (left) and no cap (right).  '
             'Dot size and colour: % of cells with |ploidy error| ≤ 0.5.  '
             'Red ×: no result (runtime error or empty output).',
             ha='left', va='bottom', fontsize=7.5, color='#333333')
    return fig


def plot_combined_four(figs_spec, args):
    """A 2x2 page of the four main-text panels, sharing one colour scale."""
    # figs_spec: list of (plot_id, title, row_labels, matrix, col_pairs)
    nonempty = [s for s in figs_spec if s[2] and s[3] is not None]
    if not nonempty:
        return None
    n_cols = max(s[3].shape[1] for s in nonempty)
    n_rows_max = max(s[3].shape[0] for s in nonempty)
    s_min, s_max = 14.0 * args.dot_scale, 180.0 * args.dot_scale
    cmap = BALLOON_CMAP
    norm = Normalize(vmin=0.0, vmax=1.0)

    fig_w = max(16.0, 0.62 * n_cols * 2 + 6.0)
    fig_h = max(10.0, 0.36 * n_rows_max * 2 + 5.0)
    fig = plt.figure(figsize=(fig_w, fig_h))
    outer = gridspec.GridSpec(2, 2, figure=fig, wspace=0.16, hspace=0.32,
                              left=0.08, right=0.90, top=0.93, bottom=0.12)
    letters = 'ABCD'
    for k, (plot_id, title, row_labels, matrix, pairs) in enumerate(PLOT_SPECS_TO_SPEC(figs_spec)):
        ax = fig.add_subplot(outer[k // 2, k % 2])
        if matrix is None or not row_labels:
            ax.set_axis_off()
            ax.set_title(F'{letters[k]}  {title}', fontsize=10, loc='left',
                         fontweight='bold', color='0.5')
            ax.text(0.5, 0.5, 'no data', transform=ax.transAxes,
                    ha='center', va='center', color='0.5')
            continue
        draw_balloon_ax(ax, matrix, row_labels, pairs, cmap, norm,
                        s_min=s_min, s_max=s_max,
                        show_xlabel=(k >= 2),
                        tick_fontsize=6.0, ytick_fontsize=6.5)
        ax.set_title(F'{letters[k]}  {title}', fontsize=10, loc='left',
                     fontweight='bold', pad=6)

    cax = fig.add_axes([0.92, 0.25, 0.012, 0.55])
    _add_colorbar(fig, cax, cmap, norm, '% cells within ±0.5')
    fig.suptitle('Ploidy-estimation accuracy  ·  % cells within ±0.5 of truth',
                 fontsize=13, fontweight='bold', y=0.98)
    fig.text(0.10, 0.03,
             'Each caller is shown twice: copy-number cap at 10 (left) and no cap (right).  '
             'Red × marks a missing result.',
             fontsize=8, color='#333333')
    return fig


def PLOT_SPECS_TO_SPEC(figs_spec):
    """Yield 4 slots aligned with PLOT_SPECS, filling missing plots with empties."""
    by_id = {s[0]: s for s in figs_spec}
    for plot_id, title, _kind in PLOT_SPECS:
        if plot_id in by_id:
            yield by_id[plot_id]
        else:
            yield (plot_id, title, [], None, [])


def build_matrix(entries, plot_id, row_labels, col_pairs):
    sub = entries[entries['plot'] == plot_id]
    lookup = {}
    for rec in sub.itertuples(index=False):
        if rec.dataset is None or (isinstance(rec.dataset, float) and np.isnan(rec.dataset)):
            continue
        lookup[(rec.dataset, rec.tool, rec.max_cn)] = rec
    matrix = np.full((len(row_labels), len(col_pairs)), np.nan, dtype=float)
    for i, ds in enumerate(row_labels):
        for j, (tool, cap) in enumerate(col_pairs):
            rec = lookup.get((ds, tool, cap))
            if rec is not None and (not rec.failed) and np.isfinite(rec.pct_within):
                matrix[i, j] = rec.pct_within
    return matrix


# ---------------------------------------------------------------------------------------------
# Demo / synthetic table, so the four figures can be previewed without benchmark files
# ---------------------------------------------------------------------------------------------

def make_demo_runs(rng=None):
    """A compact, visually plausible table covering all four main-text figures.

    Numbers are NOT real benchmark results; they only exist so the layout,
    colour scale, grouped x-axis and red-cross encoding can be inspected.
    """
    rng = np.random.default_rng(1 if rng is None else rng)
    tools = list(CALLER_ORDER)
    # Rough per-caller accuracy used to draw a readable figure (ginkgo high, etc.).
    tool_acc = {
        'scabsolute': 0.88, 'hmmcopy': 0.70, 'ginkgo': 0.92, 'copynumber': 0.55,
        'secnv': 0.62, 'sccnv': 0.48, 'scyn': 0.58, 'chisel': 0.66,
        'aneufinder': 0.40, 'flcna': 0.52,
    }
    donors = ['S01', 'S02', '234HS']
    spots = ['50', '75', '100']
    sample_type = 'PB1'
    cell_lines = GERMLINE_CELL_LINE_ORDER
    act_samples = ['TN1', 'TN2', 'TN3', 'TN4', 'TN5', 'TN6', 'TN7', 'TN8',
                   'mb157', 'BT20', 'mb453']

    # Combinations that should appear as red crosses (runtime failures).
    missing = {
        ('flcna', 'HeLa', 'inf'),
        ('aneufinder', 'COLO-829', 'inf'),
        ('chisel', 'ACT', 'inf'),
        ('sccnv', 'HCC1395', '10'),
        ('copynumber', 'ACT', '10'),
        ('scyn', 'HeLa', '10'),
    }

    runs = []

    def _one(tool, donor, sample_type, spot, cell_line, max_cn, samples, acc, tag):
        n = 40
        failed = (tool, tag, fmt_max_cn(max_cn)) in missing
        if failed:
            err = np.array([], dtype=float)
            samp = np.array([], dtype=object)
            exp = np.array([], dtype=float)
        else:
            # Mix in-window cells with a few outliers / 2x failures.
            sigma = 0.18 + 0.55 * (1.0 - acc)
            err = rng.normal(0.0, sigma, size=n * len(samples))
            # push a fraction outside the window
            n_out = int(round((1.0 - acc) * len(err)))
            if n_out:
                err[:n_out] = rng.choice([-1.0, 1.0, 2.2, -1.4]) * rng.uniform(0.6, 2.5, size=n_out)
            samp = np.repeat(np.asarray(samples, dtype=object), n)
            exp = np.full(len(err), 3.2)
        return {
            'path': F'demo/{tool}/{donor}/{spot}/{cell_line}/{fmt_max_cn(max_cn)}',
            'percell_path': '', 'stem': 'demo',
            'tool': tool, 'donor': donor, 'sampleType': sample_type,
            'avgSpotLen': spot, 'cellLine': cell_line,
            'source': '4from3' if cell_line else '4from2',
            'max_cn': float(max_cn), 'window': DEFAULT_PLOIDY_WINDOW,
            'sample': samp, 'ploidy_error': err, 'expected_ploidy': exp,
            'is_outlier': None,
            'n_cells': int(len(samp)), 'n_cells_finite': int(np.isfinite(err).sum()) if len(err) else 0,
            'failed': bool(failed),
            'caller_eval': True,
            'stem_label': 'demo',
            'dataset_key': (donor, sample_type, spot),
            'dataset_label': F'donor={donor} | sampleType={sample_type} | avgSpotLen={spot}',
        }

    for tool in tools:
        acc0 = tool_acc[tool]
        for donor in donors:
            for spot in spots:
                for cl in cell_lines:
                    for cap, acc_delta in ((10.0, 0.04), (float('inf'), -0.06)):
                        acc = min(0.98, max(0.08, acc0 + acc_delta + 0.03 * (spot == '100')))
                        runs.append(_one(tool, donor, sample_type, spot, cl, cap,
                                         samples=[cl], acc=acc, tag=cl))
        for cap, acc_delta in ((10.0, 0.02), (float('inf'), -0.08)):
            acc = min(0.97, max(0.10, acc0 + acc_delta))
            runs.append(_one(tool, 'ACT', 'tumor', '50', '', cap,
                             samples=act_samples, acc=acc, tag='ACT'))
    return runs


# ---------------------------------------------------------------------------------------------
# Legacy error-grid (the original v01 figures), kept behind --legacy
# ---------------------------------------------------------------------------------------------

def entry_ylim(sub, windows, args):
    errs = sub['ploidy_error'].to_numpy(dtype=float)
    errs = errs[np.isfinite(errs)]
    stars = []
    for _, s in sub.groupby('sample', sort=False):
        exp = s['expected_ploidy'].to_numpy(dtype=float)
        exp = exp[np.isfinite(exp)]
        if len(exp):
            e = float(np.mean(exp))
            if e > 0:
                stars += [e, -0.5 * e]
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
    ws = [r['window'] for r in sub_runs]
    return collections.Counter(ws).most_common(1)[0][0] if ws else DEFAULT_PLOIDY_WINDOW


def draw_error_panel(ax, sub, samples, window, rng, args, show_x_labels):
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
        if len(errs):
            ax.scatter(k + rng.uniform(-args.jitter, args.jitter, size=len(errs)), errs,
                       s=5, alpha=0.45, color='#31688e', linewidths=0, zorder=2,
                       rasterized=True)
        if np.isfinite(exp) and exp > 0:
            for ref in (exp, -0.5 * exp):
                ax.scatter([k], [ref], marker='*', s=55, color='#1f77b4', zorder=3)
    if not sub.empty:
        ax.axhline(0.0, color='crimson', linewidth=1.5, zorder=4)
    ax.set_xlim(-0.6, max(len(samples) - 0.4, 0.6))
    ax.set_xticks(range(len(samples)))
    if show_x_labels and samples:
        ax.set_xticklabels(samples, rotation=30, ha='right', fontsize=8)
    else:
        ax.set_xticklabels([])
    ax.grid(axis='y', linestyle=':', linewidth=0.6, alpha=0.6)
    ax.tick_params(axis='y', labelsize=8)


def plot_legacy(runs, cells, the_methods, the_datasets, args):
    window_by_entry = {(r['tool'], r['dataset_label']): r['window'] for r in runs}
    n_finite_cells = int(np.isfinite(cells['ploidy_error'].to_numpy(dtype=float)).sum())
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
    ylim_all = tuple(args.ylim) if args.ylim else entry_ylim(cells, [r['window'] for r in plotted_runs], args)
    for rowidx, method in enumerate(the_methods):
        row_cells = cells[cells['method'] == method]
        row_runs = [r for r in plotted_runs if r['tool'] == method]
        ylim_row = (tuple(args.ylim) if args.ylim
                    else entry_ylim(row_cells, [r['window'] for r in row_runs], args))
        for colidx, dataset in enumerate(the_datasets):
            ax2 = fig1.add_subplot(gs[rowidx + 1, colidx])
            sub = row_cells[row_cells['dataset'] == dataset]
            window = window_by_entry.get((method, dataset), DEFAULT_PLOIDY_WINDOW)
            draw_error_panel(ax2, sub, samples_by_dataset[dataset], window, rng, args,
                             show_x_labels=(rowidx == n_rows - 1))
            if args.ylim:
                ylim = tuple(args.ylim)
            elif args.sharey == 'all':
                ylim = ylim_all
            elif args.sharey == 'row':
                ylim = ylim_row
            else:
                ylim = entry_ylim(sub, [window], args)
            ax2.set_ylim(*ylim)
            if colidx == 0:
                ax2.set_ylabel(caller2desc.get(method, method), fontsize=10, labelpad=8)
            else:
                ax2.set_ylabel('')
            same_ylim_within_row = bool(args.ylim) or args.sharey in ('all', 'row')
            ax2.tick_params(labelleft=(colidx == 0 or not same_ylim_within_row))
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
    save_fig(fig1, args.output + '_ploidy_error_grid', args.dpi)
    plt.close(fig1)
    logging.info('wrote %s_ploidy_error_grid.pdf/.png', args.output)

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
    save_fig(fig1, args.output + '_ploidy_error_main', args.dpi)
    plt.close(fig1)

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
    save_fig(fig, args.output + '_ploidy_error_multirow', args.dpi)
    plt.close(fig)
    logging.info('wrote legacy error-grid / box-plot figures next to %s', args.output)


# ---------------------------------------------------------------------------------------------
# Go
# ---------------------------------------------------------------------------------------------

def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    gather_runs.parser_error = parser.error

    if args.demo:
        logging.info('--demo: synthesising a preview table (not real benchmark numbers)')
        runs = make_demo_runs()
    else:
        runs = gather_runs(args)

    # --- balloon-plot table ---
    entries = expand_runs_to_entries(runs, window_override=args.window)
    if entries.empty:
        sys.exit('ploidy-performances-eval: no balloon-plot entries could be built from the input.')

    tools = order_methods(set(entries['tool'].dropna().astype(str)))
    if args.methods:
        unknown = [m for m in args.methods if m not in tools]
        if unknown:
            parser.error(F'--methods: not among the evaluated tools {tools}: {unknown}')
        tools = args.methods
    max_cn_mode = str(args.max_cn).strip().lower()
    if max_cn_mode not in ('all', 'any', 'both'):
        max_cn_mode = args.max_cn
    else:
        max_cn_mode = 'all'
    col_pairs = method_columns(tools, max_cn_mode)

    if args.datasets:
        keep = entries['dataset'].fillna('').astype(str).apply(
            lambda s: any(k in s for k in args.datasets))
        entries = entries[keep | entries['dataset'].isna()]
        if entries['dataset'].dropna().empty:
            parser.error(F'--datasets: no dataset label contains any of {args.datasets}')

    # Drop placeholder rows (failed ACT run with no sample names) from the row list
    # but keep them so the tool still appears on the x-axis.
    entries_for_rows = entries[entries['dataset'].notna()
                               & (entries['dataset'].astype(str) != 'None')]

    if (args.type & 0x1):
        # testing: first 4 rows of each plot, first 2 callers
        col_pairs = col_pairs[:4]
        logging.info('testing mode: first 2 callers, first 4 datasets per plot')

    _apply_style()
    letters = 'ABCD'
    figs_spec = []
    written = []
    for k, (plot_id, title, kind) in enumerate(PLOT_SPECS):
        sub = entries_for_rows[entries_for_rows['plot'] == plot_id]
        if kind == 'germline':
            row_labels = order_germline_rows(sub)
        else:
            row_labels = order_act_rows(sub)
        if (args.type & 0x1):
            row_labels = row_labels[:4]
        if not row_labels:
            logging.warning('no datasets for plot %s; skipping the individual figure', plot_id)
            figs_spec.append((plot_id, title, [], None, col_pairs))
            continue
        matrix = build_matrix(entries, plot_id, row_labels, col_pairs)
        figs_spec.append((plot_id, title, row_labels, matrix, col_pairs))
        if (args.type & 0x2) == 0:
            fig = plot_one_balloon(entries, plot_id, title, row_labels, col_pairs, args,
                                   panel_letter=letters[k])
            if fig is not None:
                stem = args.output + '_main_' + plot_id
                save_fig(fig, stem, args.dpi)
                plt.close(fig)
                written.append(stem + '.pdf/.png')
                logging.info('wrote %s.pdf/.png', stem)

    fig4 = plot_combined_four(figs_spec, args)
    if fig4 is not None:
        stem = args.output + '_main_four'
        save_fig(fig4, stem, args.dpi)
        plt.close(fig4)
        written.append(stem + '.pdf/.png')
        logging.info('wrote %s.pdf/.png', stem)

    # Numbers behind the dots, for the paper / SI table.
    out_dir = os.path.dirname(os.path.abspath(args.output))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    tsv_path = args.output + '_pct_within_long.tsv'
    cols = ['plot', 'dataset', 'tool', 'max_cn', 'method', 'window',
            'n_cells', 'n_cells_finite', 'pct_within', 'mean_abs_ploidy_error',
            'failed', 'donor', 'sampleType', 'avgSpotLen', 'cellLine']
    tab = entries.copy()
    tab['max_cn'] = tab['max_cn'].map(fmt_max_cn)
    tab = tab[[c for c in cols if c in tab.columns]]
    tab.to_csv(tsv_path, sep='\t', index=False, float_format='%.4f')
    logging.info('wrote %s', tsv_path)

    # stderr summary: one line per (plot, dataset, method)
    show = entries_for_rows.copy()
    if not show.empty:
        show['max_cn'] = show['max_cn'].map(fmt_max_cn)
        show['pct_within'] = show['pct_within'].map(
            lambda x: '' if not np.isfinite(x) else F'{x:.1f}')
        sys.stderr.write(show[['plot', 'dataset', 'tool', 'max_cn',
                               'n_cells_finite', 'pct_within', 'failed']]
                         .to_string(index=False) + '\n')

    if args.legacy:
        frames = []
        for run in runs:
            if run['failed'] or len(run['sample']) == 0:
                continue
            frames.append(pd.DataFrame({
                'method': run['tool'],
                'dataset': run['dataset_label'],
                'sample': run['sample'],
                'ploidy_error': run['ploidy_error'],
                'expected_ploidy': run['expected_ploidy'],
            }))
        if frames:
            cells = pd.concat(frames, ignore_index=True)
            the_methods = order_methods({r['tool'] for r in runs if not r['failed']})
            the_datasets = []
            for r in sorted((x for x in runs if not x['failed']),
                            key=lambda r: (r['dataset_key'][0], r['dataset_key'][1])
                                          + _avg_spot_len_sort(r) + (r['dataset_label'],)):
                if r['dataset_label'] not in the_datasets:
                    the_datasets.append(r['dataset_label'])
            if args.methods:
                the_methods = [m for m in args.methods if m in the_methods]
            cells = cells[cells['method'].isin(the_methods) & cells['dataset'].isin(the_datasets)]
            if not cells.empty:
                plot_legacy(runs, cells, the_methods, the_datasets, args)

    sys.stderr.write('Wrote ' + ', '.join(written) + F' and {tsv_path}\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())

