#!/usr/bin/env python3
# A script to plot the benchmarking results of scWGS-based ploidy estimation
# Self-contained revision: the former plotting front-end and the
# scWGS-ploidy-performances-eval_v02.py base engine are merged into this ONE
# file, so no sibling v02 script is required any more.  stat_tests.py
# (optional; provides the pooled statistical tests) is still loaded from next
# to this file when present -- without it the figures are written and the
# tests are skipped with a warning.

"""
Benchmark ploidy estimation across CNV calling methods and datasets: the four
main-text balloon figures, the pooled statistical tests and the default LaTeX
table.

This file is self-contained: data loading, filtering, row ordering, benchmark
calculations, TSV export, the plotting layer and the statistics all live here
(the former scWGS-ploidy-performances-eval_v02.py base engine is merged in).
Place it anywhere and invoke it with the same CLI arguments as before.

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
               of cells whose |observed - expected ploidy| is within the
               0.5 tolerance window.  A red cross marks a missing result
               (runtime error, empty per-cell table, or no finite ploidy).

Visual encoding of this revision
--------------------------------
* Balloon DIAMETER = percentage of cells whose ploidy estimate is within
  +/-0.5 of the ground truth.
* Balloon COLOUR INTENSITY = log10 of the number of finite/evaluable cells
  contributing to that percentage (``n_cells_finite`` by default).
* Red x = a result was expected but is missing / failed.
* Grey hatched block = method is not applicable to that panel.

Method-specific rules
---------------------
* scAbsolute: shown once, because copy-number cap vs no-cap does not apply; it
  produces a final ploidy estimate rather than a CNV profile.
* CHISEL in ACT: shown as not applicable because phased genotypes are not
  readily available for the ACT samples.

The combined 2x2 figure follows common manuscript conventions: lowercase bold
panel letters (a-d) at the upper-left outside each axes, no embedded figure
headline/caption, fully horizontal, vertically interleaved method labels,
explicit CapAt10 labels, one shared legend below the panels, and
collision-free two-tier x-axis labels.

Statistical tests (ONE pooled family over all donors)
-----------------------------------------------------
After the figures, the ploidy benchmark is tested ONCE over the pooled set of
ALL donors from ALL four panels (COLO-829, HCC1395, HeLa, ACT).  Every
donor's per-dataset evaluations -- across the panels, sample types and read
lengths in which that donor occurs -- are aggregated into ONE observation per
method (the median, non-failed rows preferred), so each donor is exactly one
independent sample and the pooled table keeps one row per (donor, method).
The Friedman omnibus and the two-sided Wilcoxon signed-rank post-hoc tests
(paired by donor, Holm-corrected over the whole donor pool, effect sizes, BCa
bootstrap CIs) are run on this pooled table and written to a single set of
``<output>.pooled.stats.*.tsv`` files.  The LaTeX table is written BY DEFAULT
to ``<output>.pooled.stats.pairwise.tex`` after every successful pooled-stats
run (``--no-latex-table`` disables this); ``--latex-table`` prints it to
stdout from an existing ``<output>.pooled.stats.pairwise.tsv``.
``--stats-only`` re-runs the pooled tests on the existing
``<output>_pct_within_long.tsv`` without redrawing the figures.

The tests need ``stat_tests.py`` next to this script; without it the figures
are still produced and the tests are skipped with a warning.

Dataset identity
----------------
Germline-derived (S01, S02, 234HS, ...): a dataset is the combination of
average-spot-length, emulated cell-line, and original germline sample name.
The three cell-lines are split across three figures, so the row label inside
each figure is ``<donor> - <sampleType> - <avgSpotLen> bp``.

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
    python bench_results/scWGS-ploidy-performances-eval.py \
        -i '../data/*/4from2_*_ploidy_eval_summary.json' \
           '../data/*/4from3_*_ploidy_eval_summary.json' \
           '../data/*/4from2_*_ploidy_tool_eval_summary.json' \
        -o bench_results/ploidy-performances

    # style preview with in-script synthetic data (no input files needed):
    python bench_results/scWGS-ploidy-performances-eval.py --demo \
        -o bench_results/ploidy-performances-demo

Outputs (PDF + PNG at --dpi):

    <output>_main_COLO-829.{pdf,png}
    <output>_main_HCC1395.{pdf,png}
    <output>_main_HeLa.{pdf,png}
    <output>_main_ACT.{pdf,png}
    <output>_main_four.{pdf,png}          (the four panels combined)
    <output>_pct_within_long.tsv          (the numbers behind the dots)
    <output>.pooled.long.tsv              (pooled donor-level table)
    <output>.pooled.stats.*.tsv           (pooled statistical tests)
    <output>.pooled.stats.pairwise.tex    (LaTeX table, written by default)

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
from collections import OrderedDict

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.cm import ScalarMappable
import seaborn as sns

# ---------------------------------------------------------------------------------------------
# Constants (kept in step with ploidy_eval.py so the two scripts cannot drift apart)
# ---------------------------------------------------------------------------------------------

SUMMARY_SUFFIX = '_summary.json'
PERCELL_SUFFIX = '_percell.tsv'
DEFAULT_PLOIDY_WINDOW = 0.5
DEFAULT_MAX_CN = 10.0

# Methods are ordered chronologically by their publication date (exact dates are
# documented in bench_results/scWGS-performances-eval.py, CALLER_PUBLICATION):
# HMMcopy 2006, Copynumber 2012, Ginkgo 2015, AneuFinder 2016, SCCNV 2020,
# CHISEL 2021, SCYN 2021, SeCNV 2022, FLCNA 2024, scAbsolute 2024.
CALLER_ORDER = ['hmmcopy', 'copynumber', 'ginkgo', 'aneufinder', 'sccnv',
                'chisel', 'scyn', 'secnv', 'flcna', 'scabsolute']

TOOL_PRETTY = {
    'aneufinder': 'AneuFinder',
    'flcna': 'FLCNA',
    'chisel': 'CHISEL',
    'copynumber': 'Copynumber',
    'ginkgo': 'Ginkgo',
    'hmmcopy': 'HMMcopy',
    'secnv': 'SeCNV',
    'sccnv': 'SCCNV',
    'scyn': 'SCYN',
    'scabsolute': 'scAbsolute',
}

# Publication year shown after each method name, matching the main caller figures
# (see CALLER_PUBLICATION in bench_results/scWGS-performances-eval.py).
TOOL_PUBLICATION_YEAR = {
    'hmmcopy'   : 2006,
    'copynumber': 2012,
    'ginkgo'    : 2015,
    'aneufinder': 2016,
    'sccnv'     : 2020,
    'chisel'    : 2021,
    'scyn'      : 2021,
    'secnv'     : 2022,
    'flcna'     : 2024,
    'scabsolute': 2024,
}

# Method labels are built by pretty_tool() from TOOL_PRETTY + TOOL_PUBLICATION_YEAR
# (manuscript spelling + publication year); the old raw caller2desc mapping is kept
# empty for backwards compatibility with callers that expect the symbol to exist.
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


def pretty_tool(tool):
    """Manuscript spelling plus publication year, e.g. 'Ginkgo 2015'.  Kept on one
    line because the method labels are vertically staggered under the panels; the
    main caller figures use the two-line 'Ginkgo\\n2015' form instead."""
    name = TOOL_PRETTY.get(tool, tool)
    year = TOOL_PUBLICATION_YEAR.get(tool)
    return F'{name} {year}' if year else name


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


def PLOT_SPECS_TO_SPEC(figs_spec):
    """Yield 4 slots aligned with PLOT_SPECS, filling missing plots with empties."""
    by_id = {s[0]: s for s in figs_spec}
    for plot_id, title, _kind in PLOT_SPECS:
        if plot_id in by_id:
            yield by_id[plot_id]
        else:
            yield (plot_id, title, [], None, [])


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
                ax2.set_ylabel(pretty_tool(method), fontsize=10, labelpad=8)
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

def _base_main(argv=None):
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


# [REV v4] Statistical tests for the ploidy benchmark (Friedman omnibus +
# two-sided Wilcoxon signed-rank post-hoc, paired by DONOR with Holm
# correction, effect sizes, BCa bootstrap CIs), run as ONE POOLED family over
# all donors from all four balloon subplots (COLO-829, HCC1395, HeLa, ACT):
# every donor's per-dataset evaluations are aggregated (median, non-failed
# rows preferred) across all panels in which that donor occurs, so each donor
# is exactly one independent sample and the pooled long table keeps one row
# per (donor, method) -- this both restores the per-sample pivot inside
# stat_tests (which a naive panel concatenation would break with duplicate
# (dataset, method) keys) and enlarges the single test family.  Loaded lazily
# in _run_ploidy_stats() because the sibling module only needs scipy/pandas.


# ======================================================================================
# User-facing style/configuration knobs
# ======================================================================================
COUNT_FIELD = "n_cells_finite"        # use "n_cells" for total loaded cells instead
COUNT_LEGEND_LABEL = r"$\log_{10}$(evaluable cells, n)"
PLOIDY_TOLERANCE_TEXT = "±0.5"

# Academic figures usually place the explanatory title/caption outside the figure.
SHOW_COMBINED_SUPTITLE = False
SHOW_COMBINED_FOOTNOTE = False
SHOW_INDIVIDUAL_FOOTNOTE = False

# Perceptually ordered, colour-blind-safe-ish sequential blue scale.
COUNT_CMAP = LinearSegmentedColormap.from_list(
    "sample_size_blues",
    ["#eff6ff", "#c6dbef", "#9ecae1", "#6baed6", "#3182bd", "#08519c", "#08306b"],
)

TEXT = "#20262B"
MUTED = "#66717A"
GRID = "#E2E7EB"
GROUP_GRID = "#B8C1C8"
ALT_BAND = "#F7F9FA"
EDGE = "#21384E"
MISSING = "#D84A4A"
NA_FACE = "#F0F2F4"
NA_EDGE = "#A5AEB5"

# Methods for which cap/no-cap is conceptually not applicable.
SINGLE_COLUMN_TOOLS = {"scabsolute"}

# Entire caller is not applicable in these panels.
NOT_APPLICABLE_BY_PLOT = {
    "ACT": {"chisel"},
}

# Marker diameter (points).  We square these before passing to scatter because
# matplotlib's ``s`` is marker area in pt^2.  Thus the user-visible *diameter*
# is linear in the percentage, exactly as requested.
D_MIN_INDIV = 2.2
D_MAX_INDIV = 8.0
D_MIN_COMBINED = 1.5
D_MAX_COMBINED = 6.0

# The actual maximum diameter is also limited by the densest panel row pitch so
# 100% circles can never touch the rows above/below after manuscript reduction.
MAX_DIAMETER_ROW_FRACTION = 0.72

# Designed around a full-width/two-column manuscript figure.  Vector PDF output
# remains crisp if the publisher rescales it modestly.
COMBINED_FIGSIZE = (8.25, 9.25)
INDIVIDUAL_FIGSIZE = (8.1, 5.9)

# Cache display payloads because the base engine passes only the percentage
# matrix to its combined-figure function.
_PANEL_CACHE = {}

# The four balloon subplots of PLOT_SPECS.  The ploidy
# statistical tests are run ONCE over the pooled donors of all these panels
# (the constant is kept for panel-level bookkeeping and documentation).
PANEL_ORDER = ('COLO-829', 'HCC1395', 'HeLa', 'ACT')


# ======================================================================================
# Style
# ======================================================================================
def _apply_style():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "Liberation Sans"],
        "font.size": 7.0,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.55,
        "axes.edgecolor": GROUP_GRID,
        "axes.labelcolor": TEXT,
        "xtick.color": TEXT,
        "ytick.color": TEXT,
        "text.color": TEXT,
        "xtick.major.size": 0,
        "ytick.major.size": 0,
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
    })


# ======================================================================================
# Data/display helpers
# ======================================================================================
def _pct_to_area(pct, d_min, d_max):
    if not np.isfinite(pct):
        return 0.0
    p = np.clip(float(pct), 0.0, 100.0) / 100.0
    diameter = d_min + (d_max - d_min) * p
    return float(diameter ** 2)


def _count_to_log10(n):
    """Transform a positive cell count to log10(n); non-positive values are NaN."""
    try:
        x = float(n)
    except Exception:
        return float("nan")
    if not np.isfinite(x) or x <= 0:
        return float("nan")
    return float(np.log10(x))


def _count_norm(entries):
    """Shared linear normalization *after* log10-transforming cell counts."""
    if entries is None or COUNT_FIELD not in entries.columns:
        return Normalize(vmin=0.0, vmax=1.0)
    raw = np.asarray(entries[COUNT_FIELD], dtype=float)
    raw = raw[np.isfinite(raw) & (raw > 0)]
    if raw.size == 0:
        return Normalize(vmin=0.0, vmax=1.0)
    logx = np.log10(raw)
    vmax = float(np.max(logx))
    # Start at log10(1)=0 to keep the colour scale interpretable and stable.
    if not np.isfinite(vmax) or vmax <= 0:
        vmax = 1.0
    return Normalize(vmin=0.0, vmax=vmax)


def _nice_log_ticks(norm):
    """Integer log10 ticks (0,1,2,...) within the shared colour scale."""
    hi = float(norm.vmax)
    ticks = np.arange(0, int(np.floor(hi)) + 1, dtype=float)
    if ticks.size == 0:
        ticks = np.array([0.0])
    return ticks


def _limit_diameter_to_rows(fig, panel_height_frac, n_rows, requested_dmax):
    """Prevent circles from touching neighbouring rows while keeping one scale."""
    if n_rows <= 0:
        return requested_dmax
    row_pitch_pt = panel_height_frac * fig.get_figheight() * 72.0 / float(n_rows)
    return min(float(requested_dmax), MAX_DIAMETER_ROW_FRACTION * row_pitch_pt)


def _choose_single_cap(pairs):
    """Prefer cap=10 for scAbsolute if both duplicate variants exist."""
    if not pairs:
        return None
    for _tool, cap in pairs:
        try:
            if float(cap) == 10.0:
                return cap
        except Exception:
            pass
    return pairs[0][1]


def _cap_state_label(cap):
    """Publication label for a CN-calling variant.

    The repository defines the default copy-number cap as 10, so the capped
    variant is labelled ``CapAt10`` rather than the symbolic ``<=10``.
    """
    try:
        x = float(cap)
        if np.isfinite(x) and np.isclose(x, 10.0):
            return "CapAt10"
        if not np.isfinite(x):
            return "NC"
    except Exception:
        pass
    text = str(cap).strip().lower()
    if text in {"inf", "+inf", "infinity", "none", "no cap", "nocap", "nc"}:
        return "NC"
    # In normal benchmark output only 10 and infinity occur.  Retain an informative
    # fallback for user-filtered/legacy inputs rather than silently mislabelling.
    try:
        return f"CapAt{float(cap):g}"
    except Exception:
        return str(cap)


def _display_colspecs(plot_id, raw_col_pairs):
    """Build the *displayed* columns from the raw tool×cap column pairs."""
    grouped = OrderedDict()
    for tool, cap in raw_col_pairs:
        grouped.setdefault(tool, []).append((tool, cap))

    displayed = []
    na_tools = set(NOT_APPLICABLE_BY_PLOT.get(plot_id, set()))

    for tool, pairs in grouped.items():
        if tool in SINGLE_COLUMN_TOOLS:
            displayed.append({
                "tool": tool,
                "cap": _choose_single_cap(pairs),
                "cap_label": "",
                "applicable": True,
            })
            continue

        if tool in na_tools:
            # Keep the pair width so ACT remains directly comparable with the other
            # panels.  The whole pair is rendered as one hatched NA block later.
            if pairs:
                for _t, cap in pairs:
                    displayed.append({
                        "tool": tool,
                        "cap": cap,
                        "cap_label": _cap_state_label(cap),
                        "applicable": False,
                    })
            else:
                displayed.append({
                    "tool": tool,
                    "cap": None,
                    "cap_label": "",
                    "applicable": False,
                })
            continue

        for _t, cap in pairs:
            displayed.append({
                "tool": tool,
                "cap": cap,
                "cap_label": _cap_state_label(cap),
                "applicable": True,
            })

    return displayed


def _column_groups(colspecs):
    """Return [(tool, left_col, right_col, applicable), ...]."""
    if not colspecs:
        return []
    groups = []
    start = 0
    tool = colspecs[0]["tool"]
    applicable = colspecs[0]["applicable"]
    for j in range(1, len(colspecs) + 1):
        if j == len(colspecs) or colspecs[j]["tool"] != tool:
            groups.append((tool, start, j - 1, applicable))
            if j < len(colspecs):
                start = j
                tool = colspecs[j]["tool"]
                applicable = colspecs[j]["applicable"]
    return groups


def _build_payload(entries, plot_id, row_labels, raw_col_pairs):
    colspecs = _display_colspecs(plot_id, raw_col_pairs)
    sub = entries[entries["plot"] == plot_id]

    lookup = {}
    for rec in sub.itertuples(index=False):
        if rec.dataset is None or (isinstance(rec.dataset, float) and np.isnan(rec.dataset)):
            continue
        lookup[(rec.dataset, rec.tool, rec.max_cn)] = rec

    nr, nc = len(row_labels), len(colspecs)
    pct = np.full((nr, nc), np.nan, dtype=float)
    count = np.full((nr, nc), np.nan, dtype=float)
    status = np.full((nr, nc), "missing", dtype=object)  # valid / missing / na

    for i, dataset in enumerate(row_labels):
        for j, spec in enumerate(colspecs):
            if not spec["applicable"]:
                status[i, j] = "na"
                continue

            rec = lookup.get((dataset, spec["tool"], spec["cap"]))
            if rec is None:
                status[i, j] = "missing"
                continue

            n = getattr(rec, COUNT_FIELD, np.nan)
            if np.isfinite(n):
                count[i, j] = float(n)

            if (not rec.failed) and np.isfinite(rec.pct_within):
                pct[i, j] = float(rec.pct_within)
                status[i, j] = "valid"
            else:
                status[i, j] = "missing"

    payload = {
        "plot_id": plot_id,
        "row_labels": list(row_labels),
        "raw_col_pairs": list(raw_col_pairs),
        "colspecs": colspecs,
        "pct": pct,
        "count": count,
        "status": status,
    }
    _PANEL_CACHE[plot_id] = payload
    return payload


# ======================================================================================
# Figure annotations / legends
# ======================================================================================
def _panel_letter(fig, ax, letter, fontsize=13.5):
    """Large bold lowercase panel label, just outside the upper-left axes corner."""
    if not letter:
        return
    bbox = ax.get_position()
    fig.text(
        bbox.x0 - 0.030,
        bbox.y1 + 0.008,
        str(letter).lower(),
        ha="left",
        va="bottom",
        fontsize=fontsize,
        fontweight="bold",
        color=TEXT,
    )


def _size_legend(ax, d_min, d_max, values=(0, 25, 50, 75, 100), fontsize=6.3):
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.0, 0.98, f"Cells within {PLOIDY_TOLERANCE_TEXT} of truth",
            ha="left", va="top", fontsize=fontsize + 0.3, fontweight="bold", color=TEXT)
    ys = np.linspace(0.80, 0.24, len(values))
    ref_colour = COUNT_CMAP(0.55)
    for y, v in zip(ys, values):
        ax.scatter([0.20], [y], s=_pct_to_area(v, d_min, d_max),
                   facecolor=ref_colour, edgecolor=EDGE, linewidth=0.32)
        ax.text(0.43, y, f"{v}%", ha="left", va="center", fontsize=fontsize, color=TEXT)


def _status_legend(ax, fontsize=6.3):
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.scatter([0.18], [0.68], marker="x", s=28, color=MISSING, linewidth=1.0)
    ax.text(0.38, 0.68, "Missing result", va="center", fontsize=fontsize, color=TEXT)

    rect = Rectangle((0.10, 0.20), 0.16, 0.22, facecolor=NA_FACE,
                     edgecolor=NA_EDGE, hatch="///", linewidth=0.6)
    ax.add_patch(rect)
    ax.text(0.38, 0.31, "Not applicable", va="center", fontsize=fontsize, color=TEXT)
    ax.text(0.10, -0.03, "CapAt10 = CN cap at 10;  NC = no cap",
            va="top", ha="left", fontsize=max(fontsize - 0.35, 4.8), color=MUTED)


def _legend_card(ax, title, fontsize=6.2):
    """Style a small boxed legend region with a clear section title."""
    ax.set_facecolor("#FCFDFE")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#CDD5DB")
        spine.set_linewidth(0.65)
    ax.text(0.045, 0.91, title, ha="left", va="top",
            fontsize=fontsize, fontweight="bold", color=TEXT)


def _colourbar_in_card(card_ax, norm, fontsize=5.8):
    cax = card_ax.inset_axes([0.08, 0.34, 0.84, 0.20])
    sm = ScalarMappable(norm=norm, cmap=COUNT_CMAP)
    sm.set_array([])
    cb = card_ax.figure.colorbar(sm, cax=cax, orientation="horizontal")
    cb.outline.set_linewidth(0.45)
    cb.outline.set_edgecolor(GROUP_GRID)
    cb.ax.tick_params(labelsize=fontsize, length=1.8, width=0.45, color=TEXT, pad=1.2)
    cb.set_ticks(_nice_log_ticks(norm))
    cb.set_label(COUNT_LEGEND_LABEL, fontsize=fontsize + 0.1, color=TEXT, labelpad=2.0)
    return cb


def _size_legend_horizontal(ax, d_min, d_max, values=(0, 25, 50, 75, 100), fontsize=6.0, embedded=False):
    """Compact horizontal size legend for the shared combined-figure legend."""
    if not embedded:
        ax.set_axis_off()
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.text(0.0, 0.95, f"Cells within {PLOIDY_TOLERANCE_TEXT} of truth",
                ha="left", va="top", fontsize=fontsize + 0.25, fontweight="bold", color=TEXT)
        y_circle, y_text = 0.53, 0.16
    else:
        ax.text(0.045, 0.72, f"Cells within {PLOIDY_TOLERANCE_TEXT} of truth",
                ha="left", va="top", fontsize=fontsize, color=MUTED)
        y_circle, y_text = 0.39, 0.13
    xs = np.linspace(0.09, 0.91, len(values))
    ref_colour = COUNT_CMAP(0.55)
    for x, v in zip(xs, values):
        ax.scatter([x], [y_circle], s=_pct_to_area(v, d_min, d_max),
                   facecolor=ref_colour, edgecolor=EDGE, linewidth=0.30)
        ax.text(x, y_text, f"{v}%", ha="center", va="center", fontsize=fontsize, color=TEXT)


def _status_legend_horizontal(ax, fontsize=6.0, embedded=False):
    if not embedded:
        ax.set_axis_off()
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        y = 0.68
    else:
        y = 0.56
    ax.scatter([0.09], [y], marker="x", s=25, color=MISSING, linewidth=0.95)
    ax.text(0.17, y, "Missing", va="center", fontsize=fontsize, color=TEXT)
    rect = Rectangle((0.51, y - 0.10), 0.08, 0.20, facecolor=NA_FACE,
                     edgecolor=NA_EDGE, hatch="///", linewidth=0.55)
    ax.add_patch(rect)
    ax.text(0.63, y, "Not applicable", va="center", fontsize=fontsize, color=TEXT)
    ax.text(0.05, 0.16, "CapAt10 = CN cap at 10    |    NC = no cap",
            va="center", ha="left", fontsize=max(fontsize - 0.25, 4.8), color=MUTED)


def _colourbar_horizontal(fig, rect, norm, fontsize=6.0):
    cax = fig.add_axes(rect)
    sm = ScalarMappable(norm=norm, cmap=COUNT_CMAP)
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.outline.set_linewidth(0.5)
    cb.outline.set_edgecolor(GROUP_GRID)
    cb.ax.tick_params(labelsize=fontsize, length=2.0, width=0.5, color=TEXT, pad=1.5)
    cb.set_ticks(_nice_log_ticks(norm))
    cb.set_label(COUNT_LEGEND_LABEL, fontsize=fontsize + 0.2, color=TEXT, labelpad=2.5)
    return cb


def _colourbar(fig, rect, norm, fontsize=6.2):
    cax = fig.add_axes(rect)
    sm = ScalarMappable(norm=norm, cmap=COUNT_CMAP)
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cax)
    cb.outline.set_linewidth(0.5)
    cb.outline.set_edgecolor(GROUP_GRID)
    cb.ax.tick_params(labelsize=fontsize, length=2.0, width=0.5, color=TEXT, pad=2)
    cb.set_ticks(_nice_log_ticks(norm))
    cb.set_label(COUNT_LEGEND_LABEL, fontsize=fontsize + 0.2, color=TEXT, labelpad=4)
    return cb


# ======================================================================================
# Plotting layer: these are the LIVE implementations of the main-text figures
# (the former base-engine plotting functions are not carried over into this
# merged file; _base_main() reaches them through module-global lookups)
# ======================================================================================
def build_matrix(entries, plot_id, row_labels, col_pairs):
    # _base_main() only needs a matrix-like object to decide panel availability;
    # the actual displayed payload is cached and used by our custom plotters.
    payload = _build_payload(entries, plot_id, row_labels, col_pairs)
    return payload["pct"]


def draw_panel(ax, payload, norm, d_min, d_max, *,
               method_font=5.2, row_font=5.1, cap_font=4.35,
               method_rotation=0.0, show_ylabels=True,
               stagger_method_labels=False):
    pct = payload["pct"]
    count = payload["count"]
    status = payload["status"]
    rows = payload["row_labels"]
    cols = payload["colspecs"]
    groups = _column_groups(cols)
    nr, nc = pct.shape

    # Alternating caller bands improve pair grouping without a boxed-table look.
    for gi, (_tool, left, right, applicable) in enumerate(groups):
        if not applicable:
            ax.axvspan(left - 0.5, right + 0.5, facecolor=NA_FACE,
                       edgecolor=NA_EDGE, linewidth=0.0, hatch="///", zorder=0)
        elif gi % 2 == 1:
            ax.axvspan(left - 0.5, right + 0.5, facecolor=ALT_BAND,
                       edgecolor="none", zorder=0)

    # Quiet grid.
    for y in np.arange(-0.5, nr, 1.0):
        ax.axhline(y, color=GRID, linewidth=0.36, zorder=0.5)
    for x in np.arange(-0.5, nc, 1.0):
        ax.axvline(x, color=GRID, linewidth=0.28, zorder=0.5)
    for _tool, left, _right, _app in groups[1:]:
        ax.axvline(left - 0.5, color=GROUP_GRID, linewidth=0.65, zorder=1)

    xs, ys, sizes, colours = [], [], [], []
    mx, my = [], []
    for i in range(nr):
        for j in range(nc):
            if status[i, j] == "valid":
                xs.append(j)
                ys.append(i)
                sizes.append(_pct_to_area(pct[i, j], d_min, d_max))
                n = count[i, j]
                logn = _count_to_log10(n)
                colours.append(COUNT_CMAP(norm(0.0 if not np.isfinite(logn) else logn)))
            elif status[i, j] == "missing":
                mx.append(j)
                my.append(i)
            # status == 'na' is intentionally represented by the group-wide hatch.

    if xs:
        ax.scatter(xs, ys, s=sizes, c=colours, edgecolor=EDGE,
                   linewidth=0.28, zorder=3)
    if mx:
        ax.scatter(mx, my, marker="x", s=19, color=MISSING,
                   linewidth=0.85, zorder=4)

    # Not-applicable callers are indicated by the hatched band only; the
    # shared legend explains the hatch.  Avoiding text inside the data grid
    # prevents collisions and keeps the ACT panel visually quiet.

    ax.set_xlim(-0.5, nc - 0.5)
    ax.set_ylim(nr - 0.5, -0.5)

    ax.set_yticks(np.arange(nr))
    if show_ylabels:
        ax.set_yticklabels([wrap_label(str(r), width=34) for r in rows], fontsize=row_font)
        ax.tick_params(axis="y", pad=2.5)
    else:
        # Panels a-c share the same germline-derived row ordering, so the
        # upper-right panel can omit repeated labels.  This is standard
        # multi-panel practice and prevents labels from spilling into panel a.
        ax.set_yticklabels([])
        ax.tick_params(axis="y", pad=0.0)

    # Two-tier x axis.  Cap-state labels are vertical because ``CapAt10``
    # is necessarily much wider than a subcolumn.  Method names are the
    # semantic x-axis labels; they are fully horizontal and alternate
    # between two vertical levels so adjacent long names remain readable.
    ax.set_xticks(np.arange(nc))
    ax.set_xticklabels([])
    ax.tick_params(axis="x", pad=0.0)
    for j, spec in enumerate(cols):
        if not spec["cap_label"]:
            continue
        ax.text(j, -0.014, spec["cap_label"],
                transform=ax.get_xaxis_transform(),
                ha="center", va="top", rotation=90,
                fontsize=cap_font, clip_on=False, color=MUTED)

    for gi, (tool, left, right, _applicable) in enumerate(groups):
        center = 0.5 * (left + right)
        # A shallow two-level stagger is used only in the combined figure.
        # It preserves the requested near-horizontal labels while ensuring
        # adjacent long names never collide after manuscript-scale reduction.
        if stagger_method_labels:
            y_method = -0.158 if gi % 2 == 0 else -0.205
        else:
            y_method = -0.165
        # scAbsolute occupies one column rather than a two-column pair; a
        # tiny left offset gives it visual breathing room from HMMcopy.
        if tool == "scabsolute":
            center -= 0.10
        ax.text(center, y_method, pretty_tool(tool),
                transform=ax.get_xaxis_transform(),
                ha="center", va="top", rotation=0,
                fontsize=method_font,
                clip_on=False, color=TEXT)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(GROUP_GRID)
    ax.spines["bottom"].set_color(GROUP_GRID)
    ax.spines["left"].set_linewidth(0.55)
    ax.spines["bottom"].set_linewidth(0.55)
    ax.set_facecolor("white")
    return ax


def plot_one_balloon(entries, plot_id, title, row_labels, col_pairs, args,
                     panel_letter=None, s_min=None, s_max=None):
    payload = _build_payload(entries, plot_id, row_labels, col_pairs)
    norm = _count_norm(entries)
    scale = np.sqrt(max(float(args.dot_scale), 1e-6))
    d_min, d_max = D_MIN_INDIV * scale, D_MAX_INDIV * scale

    fig = plt.figure(figsize=INDIVIDUAL_FIGSIZE)
    ax_rect = [0.155, 0.275, 0.675, 0.605]
    d_max = _limit_diameter_to_rows(fig, ax_rect[3], len(row_labels), d_max)
    d_min = min(d_min, 0.34 * d_max)
    ax = fig.add_axes(ax_rect)
    draw_panel(ax, payload, norm, d_min, d_max,
               method_font=6.0, row_font=6.2, cap_font=5.0, method_rotation=0,
               show_ylabels=True, stagger_method_labels=False)
    _panel_letter(fig, ax, panel_letter, fontsize=14.5)
    ax.set_title(title, loc="left", fontsize=8.8, fontweight="bold", pad=6, color=TEXT)
    ax.set_ylabel("Dataset", fontsize=7.0, labelpad=6)

    # Compact legend block at right.
    size_ax = fig.add_axes([0.848, 0.55, 0.14, 0.27])
    _size_legend(size_ax, d_min, d_max, fontsize=6.2)
    status_ax = fig.add_axes([0.848, 0.405, 0.14, 0.105])
    _status_legend(status_ax, fontsize=6.2)
    _colourbar(fig, [0.905, 0.225, 0.014, 0.15], norm, fontsize=6.1)

    if SHOW_INDIVIDUAL_FOOTNOTE:
        fig.text(0.155, 0.035,
                 "Balloon diameter: accuracy; colour: sample size; red ×: missing; hatched: not applicable.",
                 fontsize=6.0, color=MUTED)
    return fig


def plot_combined_four(figs_spec, args):
    nonempty = [s for s in figs_spec if s[2] and s[3] is not None]
    if not nonempty:
        return None

    entries = _LAST_ENTRIES_FOR_REVISED
    norm = _count_norm(entries)
    scale = np.sqrt(max(float(args.dot_scale), 1e-6))
    d_min, d_max = D_MIN_COMBINED * scale, D_MAX_COMBINED * scale

    fig = plt.figure(figsize=COMBINED_FIGSIZE)

    # Use nearly the full page width for the panels.  The shared legend is
    # moved below the 2x2 grid instead of consuming a right-hand column.
    # This gives long method names (Copynumber, AneuFinder, scAbsolute, ...)
    # enough horizontal room to remain only slightly tilted without overlap.
    left, right = 0.082, 0.988
    bottom, top = 0.245, 0.965
    wspace, hspace = 0.095, 0.155
    panel_w = (right - left - wspace) / 2.0
    panel_h = (top - bottom - hspace) / 2.0

    max_rows = max((len(s[2]) for s in nonempty), default=1)
    d_max = _limit_diameter_to_rows(fig, panel_h, max_rows, d_max)
    d_min = min(d_min, 0.34 * d_max)

    positions = [
        [left,                  bottom + panel_h + hspace, panel_w, panel_h],
        [left + panel_w+wspace, bottom + panel_h + hspace, panel_w, panel_h],
        [left,                  bottom,                    panel_w, panel_h],
        [left + panel_w+wspace, bottom,                    panel_w, panel_h],
    ]

    letters = "abcd"
    for k, (plot_id, title, row_labels, matrix, col_pairs) in enumerate(
            PLOT_SPECS_TO_SPEC(figs_spec)):
        ax = fig.add_axes(positions[k])

        payload = _PANEL_CACHE.get(plot_id)
        if payload is None and row_labels:
            payload = _build_payload(entries, plot_id, row_labels, col_pairs)

        if payload is None or not row_labels:
            ax.set_axis_off()
            _panel_letter(fig, ax, letters[k], fontsize=15.5)
            ax.text(0.0, 1.01, title, transform=ax.transAxes,
                    fontsize=7.3, fontweight="bold", ha="left", va="bottom", color=MUTED)
            ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                    ha="center", va="center", fontsize=6.5, color=MUTED)
            continue

        # b repeats exactly the germline-derived row labels shown in a;
        # suppressing them prevents center-gutter overlap.  d keeps its
        # distinct ACT sample labels, at a slightly smaller size.
        show_y = (k != 1)
        row_fs = 4.55 if k != 3 else 4.15
        draw_panel(ax, payload, norm, d_min, d_max,
                   method_font=5.05, row_font=row_fs, cap_font=3.8,
                   method_rotation=0, show_ylabels=show_y,
                   stagger_method_labels=True)
        _panel_letter(fig, ax, letters[k], fontsize=15.5)
        ax.set_title(title, loc="left", fontsize=7.35, fontweight="bold", pad=5, color=TEXT)

    # Shared legend strip: three visually separated cards.  This avoids the
    # previous "clustered" appearance and makes each visual encoding obvious.
    acc_ax = fig.add_axes([0.082, 0.040, 0.285, 0.125])
    _legend_card(acc_ax, "Balloon diameter")
    _size_legend_horizontal(acc_ax, d_min, d_max, fontsize=5.9, embedded=True)

    stat_ax = fig.add_axes([0.390, 0.040, 0.240, 0.125])
    _legend_card(stat_ax, "Result status")
    _status_legend_horizontal(stat_ax, fontsize=5.9, embedded=True)

    count_ax = fig.add_axes([0.653, 0.040, 0.335, 0.125])
    _legend_card(count_ax, "Balloon colour")
    _colourbar_in_card(count_ax, norm, fontsize=5.8)

    if SHOW_COMBINED_SUPTITLE:
        fig.suptitle("Ploidy-estimation accuracy", fontsize=9.2, fontweight="bold", y=0.993)
    if SHOW_COMBINED_FOOTNOTE:
        fig.text(left, 0.018,
                 "Diameter = % cells within ±0.5 of truth; colour = evaluable-cell count; "
                 "red × = missing; hatch = not applicable.",
                 fontsize=5.8, color=MUTED)
    return fig


# The former monkey-patching of the separate v02 module is gone: the four
# functions above are the live implementations reached by _base_main().
#
# Store the entries table of the current run so that plot_combined_four
# (whose signature is fixed by the base engine) can normalise colours and
# rebuild panel payloads: wrap the base loader transparently.
_LAST_ENTRIES_FOR_REVISED = None
_expand_runs_to_entries_base = expand_runs_to_entries


def expand_runs_to_entries(*args, **kwargs):
    global _LAST_ENTRIES_FOR_REVISED
    out = _expand_runs_to_entries_base(*args, **kwargs)
    _LAST_ENTRIES_FOR_REVISED = out
    return out


# ======================================================================================
# LaTeX export of the pairwise statistical table
# ======================================================================================
# [FIX] The pooled ploidy statistical tests write <output>.pooled.stats.pairwise.tsv
# with the PLOIDY schema -- identity columns (plot, method_a, method_b; plot is the
# constant 'POOLED') plus pvalue_holm and the effect-size columns -- and NOT the
# (scenario, metric, caller_b, p_value_holm) schema of the CNV-caller benchmark that
# the previous version of this block assumed (that exporter could never run: it
# always aborted with "column 'scenario' missing").  This block now turns the pooled
# pairwise rows into a copy-paste-ready booktabs LaTeX table whose caption describes
# the pooled design.  By default the same table is ALSO written to
# <output>.pooled.stats.pairwise.tex after every successful pooled-stats run
# (--no-latex-table disables this); --latex-table additionally prints it to stdout.

_CALLER_DISPLAY = {
    'aneufinder': 'AneuFinder',
    'flcna'     : 'FLCNA',
    'chisel'    : 'CHISEL',
    'copynumber': 'Copynumber',
    'ginkgo'    : 'Ginkgo',
    'hmmcopy'   : 'HMMcopy',
    'secnv'     : 'SeCNV',
    'sccnv'     : 'SCCNV',
    'scyn'      : 'SCYN',
    'scabsolute': 'scAbsolute',
}


def _format_reference(reference):
    """Human-readable name of the reference caller/method (e.g. 'ginkgo|10')."""
    ref = str(reference)
    if '|' in ref:
        tool, cap = ref.split('|', 1)
        return F'{_CALLER_DISPLAY.get(tool, tool.capitalize())} capped at CN {cap}'
    return _CALLER_DISPLAY.get(ref, ref.capitalize())


def _ploidy_legend(ref_desc, show_panel_col=False, show_ref_col=False):
    """Table legend for the pooled ploidy-estimation pairwise table.

    [FIX] The previous caption described the per-scenario stratification of the
    CNV-caller benchmark (Hap_0 vs Hap_1), which the pooled ploidy table does
    not have; it now describes the actual pooled design.
    [REV v4] Describes exactly the four reported statistics (n, p, r, 95% CI
    of r), matching the table columns.
    """
    ref_clause = (F'({ref_desc})'
                  if not show_ref_col else
                  F'($a$; see the Method $a$ column, {ref_desc} by default)')
    panel_clause = ('per panel, ' if show_panel_col else '')
    return (
        'Pairwise comparison of ploidy-estimation accuracy between the reference method '
        F'{ref_clause} and each other method ($b$), {panel_clause}pooled over all four '
        'benchmark panels (COLO-829, HCC1395, HeLa, ACT). Every donor contributes ONE '
        'observation per method -- the median of that donor\'s per-dataset evaluations '
        'across all panels in which it occurs, preferring non-failed rows -- so each '
        'donor is one independent sample (the same germline donors underlie the three '
        'emulated cell-line panels). The metric is the percentage of cells whose ploidy '
        'estimate is within $\\pm$0.5 of the ground truth. Each row reports, in this '
        'order: $n$, the effective sample size (number of donors paired for that '
        'comparison); $p$, the two-sided Wilcoxon signed-rank test on the paired '
        'per-donor differences, Holm--Bonferroni-adjusted within the single pooled '
        'family, bold at the 0.05 family-wise level; $r$, the matched-pairs '
        'rank-biserial correlation on the per-donor differences (positive means the '
        'reference outperforms method $b$); and the 95\\% percentile-bootstrap CI of '
        '$r$ obtained by resampling the donors.'
    )


def _tex_escape(s):
    """Escape LaTeX special characters in a table text cell."""
    return (str(s)
            .replace('\\', r'\textbackslash{}')
            .replace('&', r'\&')
            .replace('%', r'\%')
            .replace('_', r'\_')
            .replace('#', r'\#')
            .replace('$', r'\$'))


def _fmt_pvalue_holm(x, alpha=0.05):
    """Format one Holm-adjusted P value for a LaTeX table cell.

    Missing values become '--'; underflow (P = 0) becomes '<2.2e-16'; values
    >= 1e-3 use three decimals, smaller ones scientific notation.  Values
    significant at level `alpha` are wrapped in \\textbf.
    """
    try:
        x = float(x)
    except (TypeError, ValueError):
        return '--'
    if x != x or x in (float('inf'), float('-inf')):  # NaN / inf
        return '--'
    if x == 0.0:
        s = '$<2.2\\times10^{-16}$'
    elif x < 0.001:
        mant, exp = F'{x:.2e}'.split('e')
        s = F'${mant}\\times10^{{{int(exp)}}}$'
    else:
        s = F'{x:.3f}'
    return F'\\textbf{{{s}}}' if x < alpha else s


def _fmt_signed_float(x, decimals=1):
    """Format an effect-size cell; missing values become '--'."""
    try:
        x = float(x)
    except (TypeError, ValueError):
        return '--'
    if x != x or x in (float('inf'), float('-inf')):
        return '--'
    return F'{x:.{decimals}f}'


def _fmt_effect(x, decimals=2):
    """Format a rank-biserial effect-size cell; missing values become '--'."""
    return _fmt_signed_float(x, decimals=decimals)


def _fmt_ci(low, high, decimals=2):
    """Format a 95% confidence-interval cell as '[low, high]'; missing -> '--'."""
    try:
        lo, hi = float(low), float(high)
    except (TypeError, ValueError):
        return '--'
    if lo != lo or hi != hi or lo in (float('inf'), float('-inf')) \
            or hi in (float('inf'), float('-inf')):
        return '--'
    return F'[{lo:.{decimals}f}, {hi:.{decimals}f}]'


def _ploidy_latex_lines(tsv_path, reference='ginkgo|10',
                        table_label='tab:pairwise', alpha=0.05,
                        caption_note=None):
    """Build the booktabs LaTeX table from a pooled *.stats.pairwise.tsv file.

    [FIX] Reads the PLOIDY schema written by stat_tests.run_ploidy_benchmark_stats:
    required columns are (method_b, pvalue_holm); the optional columns
    (plot, method_a, n_pairs) are shown whenever they add information -- the
    Panel column only when more than one panel occurs (the pooled file has the
    constant 'POOLED'), and the Method $a$ column only when several reference
    methods occur (e.g. a hand-made all-pairs table).
    [REV v4] After the identity columns, every row carries exactly these four
    statistics, in this order and nothing else: the effective sample size n
    (number of donors), the Holm-adjusted two-sided p, the effect size r
    (matched-pairs rank-biserial correlation) and the 95% bootstrap CI of r.
    The previous 'Median diff. (pp)' column (a statistic outside this set) is
    gone: it remains available in the pairwise TSV for exploration.
    Returns the list of table lines, or None on any problem (a reason is logged).
    """
    import pandas as pd
    if not os.path.isfile(tsv_path):
        logging.error('pooled pairwise stats file not found: %s', tsv_path)
        logging.error('run the script with the statistical tests enabled (default) to generate it first')
        return None
    tab = pd.read_csv(tsv_path, sep='\t')
    n_col = next((c for c in ('n_pairs', 'n_clusters') if c in tab.columns), None)
    missing = [c for c in ('method_b', 'pvalue_holm', 'rank_biserial_r',
                           'ci95_r_low', 'ci95_r_high') if c not in tab.columns]
    if missing:
        logging.error('column(s) %s missing from %s (available: %s)',
                      missing, tsv_path, ', '.join(map(str, tab.columns)))
        logging.error('rerun the pooled statistical tests with the current '
                      'stat_tests.py to obtain n, p, r and the 95%% CI of r')
        return None
    sub = tab.dropna(subset=['method_b']).copy()
    if sub.empty:
        logging.error('no pairwise rows in %s', tsv_path)
        return None
    # Compared methods keep their order of first appearance in the input file
    # (stat_tests already sorts by panel / reference / method).
    method_b_order = list(dict.fromkeys(sub['method_b'].tolist()))
    sub['method_b'] = pd.Categorical(sub['method_b'], categories=method_b_order, ordered=True)
    sub = sub.sort_values(['method_b'])

    show_panel = ('plot' in sub.columns and sub['plot'].astype(str).nunique() > 1)
    show_ref = ('method_a' in sub.columns and sub['method_a'].astype(str).nunique() > 1)
    if n_col is None:
        # [FIX] keep the documented four-statistic schema (n, p, r, 95% CI of r)
        # even when the file has no n column at all: the cell then renders as '--'.
        # The previous code created this column but dropped it from the header and
        # every row (show_n was computed before the fallback), so the table silently
        # lost the sample size the caption promises.
        sub['n_effective'] = float('nan')
        n_col = 'n_effective'
    show_n = True

    # ---- header: identity columns + exactly (n, p, r, 95% CI) ----
    header_cells = []
    if show_panel:
        header_cells.append('Panel')
    if show_ref:
        header_cells.append('Method $a$')
    header_cells.append('Method $b$')
    if show_n:
        header_cells.append('$n$')
    header_cells += ['$p$', '$r$', '95\\% CI']
    n_id = (1 if show_panel else 0) + (1 if show_ref else 0) + 1   # identity columns
    # identity columns incl. the n column are left/integer-style, p/r/CI numeric
    colspec = 'l' * (n_id + (1 if show_n else 0)) + 'r' * 3

    legend = _ploidy_legend(_format_reference(reference),
                            show_panel_col=show_panel, show_ref_col=show_ref)
    if caption_note:
        legend = F'{caption_note} {legend}'
    lines = [
        F'% LaTeX table generated by {os.path.basename(sys.argv[0])} '
        '(requires \\usepackage{booktabs})',
        '\\begin{table}[htbp]',
        '  \\centering',
        F'  \\caption{{{legend}}}',
        F'  \\label{{{table_label}}}',
        F'  \\begin{{tabular}}{{{colspec}}}',
        '    \\toprule',
        F'    {" & ".join(header_cells)} \\\\',
        '    \\midrule',
    ]
    for rec in sub.itertuples(index=False):
        cells = []
        if show_panel:
            cells.append(_tex_escape(rec.plot))
        if show_ref:
            cells.append(_tex_escape(rec.method_a))
        cells.append(_tex_escape(rec.method_b))
        if show_n:
            cells.append(_fmt_signed_float(getattr(rec, n_col), 0))
        cells.append(_fmt_pvalue_holm(rec.pvalue_holm, alpha=alpha))
        cells.append(_fmt_effect(rec.rank_biserial_r))
        cells.append(_fmt_ci(rec.ci95_r_low, rec.ci95_r_high))
        lines.append(F'    {" & ".join(cells)} \\\\')
    lines += [
        '    \\bottomrule',
        '  \\end{tabular}',
        '\\end{table}',
    ]
    return lines


def emit_latex_stats_table(tsv_path, reference='ginkgo|10', legend_kind='ploidy',
                           table_label='tab:pairwise', alpha=0.05, caption_note=None):
    """Print a copy-paste-ready booktabs LaTeX table from a pooled *.stats.pairwise.tsv file.

    See _ploidy_latex_lines for the table contents (``legend_kind`` is kept for
    call-compatibility; the only schema this script writes is the ploidy one).
    Returns 0 on success, 1 on any problem.
    """
    lines = _ploidy_latex_lines(tsv_path, reference=reference,
                                table_label=table_label, alpha=alpha,
                                caption_note=caption_note)
    if lines is None:
        return 1
    print('\n'.join(lines))
    return 0


def write_latex_stats_table(tsv_path, tex_path, reference='ginkgo|10',
                            table_label='tab:pairwise', alpha=0.05,
                            caption_note=None):
    """[REV] Write the booktabs LaTeX table to `tex_path` (default pipeline output).

    Same table as emit_latex_stats_table, but into a file instead of stdout.
    Returns 0 on success, 1 on any problem (logged; never raises -- a cosmetics
    failure must not break the benchmark run).
    """
    try:
        lines = _ploidy_latex_lines(tsv_path, reference=reference,
                                    table_label=table_label, alpha=alpha,
                                    caption_note=caption_note)
        if lines is None:
            return 1
        with open(tex_path, 'w') as fh:
            fh.write('\n'.join(lines) + '\n')
        logging.info('LaTeX pairwise table written to %s', tex_path)
        return 0
    except Exception as exc:
        logging.warning('could not write the LaTeX table %s: %s', tex_path, exc)
        return 1


# ======================================================================================
# Donor-level pooling of the ploidy long table (all panels -> one test family)
# ======================================================================================
# [REV v4] These helpers implement the pooling requested for statistical power:
# every donor of every panel becomes exactly one independent sample.  The
# germline donors are reused by the three emulated cell-line panels (their
# dataset labels are identical there), so pooling first aggregates all of a
# donor's rows -- across panels and across the donor's sampleType / avgSpotLen
# datasets -- into ONE observation per method (median, non-failed rows
# preferred).  The pooled table therefore has exactly one row per (donor,
# method): 'dataset' holds the donor id and 'plot' is the constant 'POOLED',
# which restores the per-sample pivot invariant inside stat_tests and makes
# every Holm family span the entire donor pool.

def _usable_value(v):
    """True when a cluster-key cell holds a usable (non-null, non-empty) value."""
    import pandas as pd
    try:
        if pd.isna(v):
            return False
    except (TypeError, ValueError):
        pass
    return str(v).strip().lower() not in ('', 'none', 'nan', 'null')


def _to_bool_series(s):
    """Coerce a 'failed'-like column to a robust boolean Series."""
    import pandas as pd
    if s is None:
        return pd.Series(dtype=bool)
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    txt = s.astype(str).str.strip().str.lower()
    return txt.isin(('true', '1', 'yes', 't', 'y'))


def _donor_from_dataset_label(label, split=True):
    """Best-effort independent-unit (donor) id parsed from a row label.

    The base engine composes the germline-derived row labels as
    ``donor · sampleType · avgSpotLen`` (the emulated cell line is
    deliberately omitted, so identical labels appear in the three cell-line
    panels); with ``split=True`` the donor is therefore the first
    middle-dot-separated segment, which pools the sample-type / read-length
    variants of the same donor.  Labels without a middle dot (e.g. the ACT
    sample labels) are returned unchanged, i.e. each such dataset is treated
    as its own donor unless a donor column is available.  With
    ``split=False`` the full label is the unit (naive per-dataset mode).
    """
    s = str(label).strip()
    if not s or s.lower() in ('none', 'nan', 'nat', 'null'):
        return None
    if split:
        for sep in ('·', '•'):
            if sep in s:
                head = s.split(sep, 1)[0].strip()
                return head if head else s
    return s


def _resolve_unit_ids(tab, cluster_key_cols, naive=False):
    """Independent-unit (donor) id for every row of the ploidy long table.

    Resolution order, per row:
    1. the comma-separated columns named by ``--stats-cluster-key``
       (default ``donor``), when the column exists and the row's value is
       usable;
    2. otherwise the donor parsed from the dataset label (first
       middle-dot-separated segment; the whole label in naive mode).

    Because the germline-derived dataset labels -- and the donor column
    values, when present -- are identical in the three cell-line panels, this
    key pools each donor's rows across panels while keeping different donors
    apart.
    """
    import pandas as pd
    cluster_cols = [c for c in (cluster_key_cols or []) if c in tab.columns]
    if cluster_key_cols:
        missing = [c for c in cluster_key_cols if c not in cluster_cols]
        if missing:
            logging.info('pooled stats: cluster column(s) %s absent from the long '
                         'table; the affected rows fall back to the donor parsed '
                         'from the dataset label', ', '.join(map(str, missing)))
    units, n_fallback = [], 0
    for _idx, row in tab.iterrows():
        u = None
        if (not naive) and cluster_cols:
            vals = [row[c] for c in cluster_cols]
            if all(_usable_value(v) for v in vals):
                u = ' | '.join(str(v).strip() for v in vals)
        if u is None:
            u = _donor_from_dataset_label(row['dataset'], split=not naive)
            if cluster_cols and not naive:
                n_fallback += 1
        units.append(u)
    return pd.Series(units, index=tab.index), cluster_cols, n_fallback


def _pool_units(tab, cluster_key_cols, naive=False):
    """Pool all panels of the ploidy long table to ONE row per (unit, method).

    ``tab`` must already be cleaned (placeholder rows dropped, per-panel
    not-applicable caller rows dropped).  Every independent unit's rows --
    which may span several panels (the germline donors are reused by the
    three emulated cell-line panels) and several datasets (sampleType /
    avgSpotLen variants of the same donor) -- are reduced to a single
    observation per method:

    * non-failed rows are preferred, so a method that succeeded on a donor in
      at least one panel still contributes a valid observation;
    * every numeric column (the benchmark value(s), including any
      scenario-specific value columns, and the cell counts) becomes the
      median over the unit's usable rows, matching the cluster aggregation
      ('median') of the previous per-panel tests;
    * ``dataset`` becomes the unit id and ``plot`` becomes the constant
      'POOLED', so the pooled table satisfies the one-row-per-(dataset,
      method) pivot invariant required by stat_tests and every Holm family
      spans the whole donor pool.

    Returns ``(pooled, provenance, cluster_cols, n_fallback)``: ``pooled``
    carries the same columns as ``tab``; ``provenance`` records, per pooled
    row, how many source rows were aggregated and from which panels and
    datasets.
    """
    import pandas as pd
    tab = tab.copy()
    units, cluster_cols, n_fallback = _resolve_unit_ids(tab, cluster_key_cols, naive=naive)
    tab['_unit'] = units
    n_unresolved = int(tab['_unit'].isna().sum())
    if n_unresolved:
        logging.warning('pooled stats: %d row(s) without a resolvable independent unit '
                        'were dropped from the pooled tests', n_unresolved)
        tab = tab[tab['_unit'].notna()]

    if 'method' in tab.columns:
        method_cols = ['method']
    else:
        method_cols = [c for c in ('tool', 'max_cn') if c in tab.columns]
    if tab.empty:
        return tab.drop(columns=['_unit']), pd.DataFrame(), cluster_cols, n_fallback
    if not method_cols:
        logging.warning('pooled stats: the long table has neither "method" nor '
                        '"tool" columns; the panels cannot be pooled')
        return None, None, cluster_cols, n_fallback

    # Defensive: if the long table happens to be in a true long format with
    # scenario / metric row dimensions, keep those dimensions apart while
    # pooling so different metrics are never median-mixed.
    dim_cols = [c for c in ('scenario', 'metric') if c in tab.columns]
    group_cols = ['_unit'] + method_cols + dim_cols

    agg_rows, prov_rows = [], []
    for key, g_all in tab.groupby(group_cols, dropna=False, sort=False):
        keymap = dict(zip(group_cols, key if isinstance(key, tuple) else (key,)))
        unit = keymap['_unit']
        g = g_all
        if 'failed' in g.columns:
            ok = ~_to_bool_series(g['failed'])
            if ok.any():
                g = g[ok]
        rec = {}
        for col in tab.columns:
            if col == '_unit':
                continue
            if col in keymap:
                rec[col] = keymap[col]
            elif col == 'dataset':
                rec[col] = unit
            elif col == 'plot':
                rec[col] = 'POOLED'
            elif col == 'failed':
                rec[col] = (bool(_to_bool_series(g['failed']).all())
                            if 'failed' in g.columns else False)
            elif pd.api.types.is_numeric_dtype(g[col]):
                vals = pd.to_numeric(g[col], errors='coerce')
                vals = vals[vals.notna()]
                rec[col] = float(vals.median()) if len(vals) else float('nan')
            else:
                nn = g[col][g[col].notna()]
                rec[col] = nn.iloc[0] if len(nn) else None
        agg_rows.append(rec)
        prov_rows.append({
            'pooled_unit': unit,
            'n_source_rows': int(len(g_all)),
            'n_failed_source_rows': (int(_to_bool_series(g_all['failed']).sum())
                                     if 'failed' in g_all.columns else 0),
            'source_plots': ('+'.join(sorted({str(p) for p in g_all['plot']}))
                             if 'plot' in g_all.columns else ''),
            'source_datasets': '+'.join(sorted({str(d) for d in g_all['dataset']})),
        })
    cols = [c for c in tab.columns if c != '_unit']
    pooled = pd.DataFrame(agg_rows, columns=cols).reset_index(drop=True)
    prov = pd.DataFrame(prov_rows)
    return pooled, prov, cluster_cols, n_fallback


# ======================================================================================
# Entrypoint
# ======================================================================================
def _stat_arg_parser():
    """[REV] Flags consumed by this wrapper; everything else is forwarded to the
    base engine parser, build_parser() (the same -o default is re-injected
    when delegating)."""
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument('-o', '--output', default='scWGS-ploidy-performances')
    p.add_argument('--no-stats', dest='stats', action='store_false', default=True,
                   help='Skip the statistical tests (default: run them after the figures).')
    p.add_argument('--stats-reference', default='ginkgo|10', metavar='METHOD',
                   help="Reference method for the pairwise tests, as 'tool|max_cn' "
                        "(default: ginkgo|10, i.e. Ginkgo capped at 10).")
    p.add_argument('--stats-boot', type=int, default=10000, metavar='N',
                   help='Bootstrap resamples for the BCa CIs (default: 10000).')
    p.add_argument('--stats-seed', type=int, default=1, metavar='SEED',
                   help='Seed of the bootstrap RNG (default: 1).')
    p.add_argument('--stats-alpha', type=float, default=0.05, metavar='ALPHA',
                   help='Family-wise alpha for Holm rejection (default: 0.05).')
    # [REV v2/v4] cluster (independent-unit) option.  In the POOLED tests the
    # cluster key defines the independent unit: all of a unit's rows -- across
    # the panels and datasets in which it occurs -- are aggregated to ONE
    # observation per method, and the unit (the donor) is the independent
    # sample of the single pooled test family.
    p.add_argument('--stats-cluster-key', default='donor', metavar='COLS',
                   help='Comma-separated columns defining the independent unit '
                        '(donor) of the pooled tests. Default: "donor" -- the '
                        '"donor" column when usable, otherwise the donor parsed '
                        'from the dataset label (its first middle-dot-separated '
                        'segment). Pass "none" to use the dataset label itself as '
                        'the unit (discouraged: datasets of the same donor then '
                        'count as separate samples).')
    # [REV v4] iterate on the pooled tests without redrawing the figures.
    p.add_argument('--stats-only', action='store_true', default=False,
                   help='Run only the pooled ploidy statistical tests (all donors '
                        'from all panels, one test family) on the existing '
                        '<output>_pct_within_long.tsv and exit, without redrawing '
                        'the figures.')
    # [REV] The LaTeX pairwise table is generated BY DEFAULT (written to
    # <output>.pooled.stats.pairwise.tex after every successful pooled-stats
    # run); --no-latex-table opts out and --latex-table prints it to stdout.
    p.add_argument('--no-latex-table', dest='latex_table_auto',
                   action='store_false', default=True,
                   help='Do not write <output>.pooled.stats.pairwise.tex after '
                        'the pooled tests (default: write it).')
    p.add_argument('--latex-table', action='store_true', default=False,
                   help='Print one booktabs LaTeX table built from the existing '
                        '<output>.pooled.stats.pairwise.tsv file (columns method_a, '
                        'method_b, pvalue_holm; pooled over all panels, donors as '
                        'independent samples) to stdout and exit, without running '
                        'the pipeline.')
    return p


def _run_ploidy_stats(known):
    """[REV v4] ONE pooled set of statistical tests over all donors of all panels.

    Motivation: the previous per-panel families split the donors across four
    small Friedman/Wilcoxon tests, which clearly reduces statistical power.
    Design: all donors of all four panels (COLO-829, HCC1395, HeLa, ACT) are
    pooled into a single family, with each donor denoting one independent
    sample.  The germline-derived row labels (and the underlying donors) are
    reused across the three cell-line panels, so the same donor legitimately
    occurs in three panels: a naive concatenation would (a) duplicate
    (dataset, method) keys and break the per-sample pivot inside stat_tests
    (ValueError: duplicate entries) and (b) count one donor three times,
    violating the independence of the samples.  This function therefore first
    aggregates, for every donor, all of that donor's per-dataset evaluations
    -- across the panels, sample types and read lengths in which the donor
    occurs -- into ONE observation per method (the median, matching the
    cluster aggregation used previously; non-failed rows are preferred).  The
    resulting pooled table keeps exactly one row per (donor, method), so
    ``stat_tests.run_ploidy_benchmark_stats`` is called ONCE on the whole
    pool: Friedman omnibus + two-sided Wilcoxon signed-rank post-hoc paired
    by donor, Holm-corrected within each (scenario, metric) family over the
    entire donor pool, with effect sizes and BCa bootstrap CIs.  A single set
    of ``<output>.pooled.stats.*.tsv`` files is written; the pooled
    donor-level table and its provenance are exported as
    ``<output>.pooled.long.tsv`` for inspection.

    Returns True when the pooled tests ran and wrote their stats files, and
    False when they were skipped for any reason (missing stat_tests module,
    missing/invalid long table, empty pool, too few units/methods, or a
    failed run).  [REV] The caller uses this to decide whether the default
    LaTeX table can be generated.
    """
    import pandas as pd
    here = os.path.dirname(os.path.abspath(__file__))
    if here not in sys.path:
        sys.path.insert(0, here)
    try:
        import stat_tests
    except ImportError:
        logging.warning('stat_tests.py not found next to this script: statistical tests skipped')
        return False
    tsv_path = F'{known.output}_pct_within_long.tsv'
    if not os.path.isfile(tsv_path):
        logging.warning('ploidy statistics skipped: %s not found', tsv_path)
        return False
    tab = pd.read_csv(tsv_path, sep='\t')
    for col in ('plot', 'dataset'):
        if col not in tab.columns:
            logging.warning('ploidy statistics skipped: %s has no %r column', tsv_path, col)
            return False
    if 'method' not in tab.columns and 'tool' not in tab.columns:
        logging.warning('ploidy statistics skipped: %s has neither "method" nor "tool" columns',
                        tsv_path)
        return False

    # Independent-unit (donor) resolution used to pool the panels.
    key_txt = (known.stats_cluster_key or '').strip()
    if key_txt.lower() in ('none', 'naive', 'off'):
        cluster_key_cols, naive = [], True
        logging.warning('pooled stats: naive per-dataset units requested '
                        '(--stats-cluster-key none): datasets of the same donor '
                        'will count as separate samples (discouraged)')
    elif key_txt:
        cluster_key_cols, naive = [c.strip() for c in key_txt.split(',') if c.strip()], False
    else:
        cluster_key_cols, naive = ['donor'], False

    # Drop figure-only placeholder rows (failed runs with no dataset label).
    ds = tab['dataset']
    keep = (ds.notna() & (ds.astype(str).str.strip() != '')
            & (ds.astype(str) != 'None'))
    tab = tab[keep]
    if tab.empty:
        logging.warning('ploidy statistics skipped: no evaluable (dataset, method) rows in %s',
                        tsv_path)
        return False

    # Mirror the figure's "not applicable" encoding (CHISEL in ACT): the
    # caller is dropped from that panel only, so it keeps its valid rows in
    # the other panels of the pool.
    if 'tool' in tab.columns:
        for plot, na_tools in NOT_APPLICABLE_BY_PLOT.items():
            for na_tool in na_tools:
                is_na = ((tab['plot'].astype(str) == str(plot))
                         & (tab['tool'].astype(str).str.lower() == str(na_tool).lower()))
                n0 = len(tab)
                tab = tab[~is_na]
                if len(tab) < n0:
                    logging.info('pooled stats: dropped %d not-applicable %s row(s) from panel %s',
                                 n0 - len(tab), na_tool, plot)
    else:
        logging.info('pooled stats: no "tool" column; per-panel not-applicable rules skipped')

    # Report genuine (plot, dataset, method) duplicates from the input
    # summaries; they are absorbed by the per-donor median below.
    if 'method' in tab.columns:
        method_cols = ['method']
    else:
        method_cols = [c for c in ('tool', 'max_cn') if c in tab.columns]
    dup_key = ['plot', 'dataset'] + method_cols
    dup = tab.duplicated(subset=dup_key, keep=False)
    if dup.any():
        logging.warning('pooled stats: %d rows are duplicated (plot, dataset, %s) entries '
                        'from the input summaries; they are aggregated into the same '
                        'donor-method median (check the input summaries for duplicate '
                        'evaluations of the same dataset)',
                        int(dup.sum()), ', '.join(method_cols))

    # Pool: all panels -> one row per (independent unit, method).
    pooled, prov, cluster_cols, n_fallback = _pool_units(tab, cluster_key_cols, naive=naive)
    if pooled is None or pooled.empty:
        logging.warning('ploidy statistics skipped: the pooled table is empty')
        return False

    # Export the pooled donor-level table (+ provenance) for inspection.
    export = pooled.copy()
    for col in ('n_source_rows', 'n_failed_source_rows', 'source_plots', 'source_datasets'):
        if col in prov.columns:
            export[col] = prov[col].to_numpy()
    pooled_tsv = F'{known.output}.pooled.long.tsv'
    export.to_csv(pooled_tsv, sep='\t', index=False)
    logging.info('pooled donor-level table (one row per independent unit x method, with '
                 'provenance) written to %s', pooled_tsv)

    # The pivot invariant required by stat_tests: exactly one row per
    # (unit, method) in the pooled table.
    key = (['dataset', 'method'] if 'method' in pooled.columns
           else ['dataset'] + [c for c in ('tool', 'max_cn') if c in pooled.columns])
    dup_pool = pooled.duplicated(subset=key, keep=False)
    if dup_pool.any():
        logging.error('pooled stats: the pooled table still contains %d rows with duplicate '
                      '(unit, method) keys (%s); the long table carries additional per-row '
                      'dimensions that stat_tests cannot pivot -- please inspect %s',
                      int(dup_pool.sum()), ', '.join(key), pooled_tsv)
        return False

    n_units = int(pooled['dataset'].nunique())
    if 'method' in pooled.columns:
        n_m = int(pooled['method'].nunique())
    else:
        n_m = int(pooled.groupby([c for c in ('tool', 'max_cn') if c in pooled.columns],
                                 dropna=False).ngroups)
    if 'source_plots' in prov.columns:
        n_multi = int(prov['source_plots'].astype(str).str.contains('+', regex=False).sum())
    else:
        n_multi = 0

    if naive:
        unit_desc = 'dataset label (naive per-dataset units)'
    elif cluster_cols:
        unit_desc = '+'.join(cluster_cols)
        if n_fallback:
            unit_desc += F' (+{n_fallback} rows via the dataset-label donor fallback)'
    else:
        unit_desc = 'donor parsed from the dataset label'

    logging.info('pooling all panels for the ploidy statistics: %d long-table rows from %d '
                 'panel(s) -> %d independent units x %d methods; %d unit(s) contribute '
                 'evaluations from more than one panel (aggregated as the median across '
                 'panels); unit key: %s',
                 len(tab), int(tab['plot'].nunique()), n_units, n_m, n_multi, unit_desc)
    if n_multi == 0:
        logging.info('pooled stats: no independent unit contributes rows from more than one '
                     'panel; the pool is a plain concatenation of the per-panel donors')

    if n_units < 2 or n_m < 3:
        logging.warning('pooled tests skipped: only %d independent units x %d methods; '
                        'Friedman/pairwise tests are not meaningful here', n_units, n_m)
        return False

    out_prefix = F'{known.output}.pooled'
    logging.info('running POOLED ploidy statistical tests over all donors of all panels '
                 '(%d independent units x %d methods; reference: %s; unit key: %s) -> %s.*.tsv',
                 n_units, n_m, known.stats_reference, unit_desc, out_prefix)
    # cluster_key_cols=[] (naive per-dataset mode) is intentional: the pooled
    # rows ARE the independent units (donors) already, so no further clustering
    # must be applied inside stat_tests.
    settings = stat_tests.run_ploidy_benchmark_stats(
        pooled, out_prefix,
        reference=known.stats_reference,
        cluster_key_cols=[],
        cluster_agg='median',
        n_resamples=known.stats_boot,
        seed=known.stats_seed,
        alpha=known.stats_alpha)
    return settings is not None


def _auto_write_latex_table(known, ran_stats):
    """[REV] Write <output>.pooled.stats.pairwise.tex BY DEFAULT after a successful
    pooled-stats run, so the manuscript table can never drift from the tested
    numbers.  Silently skipped when the tests did not run, when the user passed
    --no-latex-table, or when the pairwise TSV is missing (a stale table from an
    earlier run must never be presented as the result of this run)."""
    if not ran_stats or not known.latex_table_auto:
        return
    pooled_tsv = F'{known.output}.pooled.stats.pairwise.tsv'
    if not os.path.isfile(pooled_tsv):
        logging.warning('default LaTeX table skipped: %s not found', pooled_tsv)
        return
    write_latex_stats_table(
        pooled_tsv,
        F'{known.output}.pooled.stats.pairwise.tex',
        reference=known.stats_reference,
        table_label='tab:scwgs-ploidy-pairwise-pooled',
        alpha=known.stats_alpha,
        caption_note='(All panels pooled; each donor is one independent sample.)')


def main(argv=None):
    stat_parser = _stat_arg_parser()
    known, rest = stat_parser.parse_known_args(argv if argv is not None else sys.argv[1:])
    if known.latex_table:
        # [REV v4] ONE table built from the single pooled test family.
        pooled_tsv = F'{known.output}.pooled.stats.pairwise.tsv'
        if not os.path.isfile(pooled_tsv):
            logging.error('pooled pairwise stats file not found: %s', pooled_tsv)
            logging.error('run the script with the statistical tests enabled (default) to '
                          'generate it first (the per-panel stats files are no longer '
                          'written)')
            return 1
        print(F'\n% ==== pooled (all panels, donors as independent samples): {pooled_tsv} ====')
        return emit_latex_stats_table(
            pooled_tsv,
            reference=known.stats_reference,
            legend_kind='ploidy',
            table_label='tab:scwgs-ploidy-pairwise-pooled',
            alpha=known.stats_alpha,
            caption_note='(All panels pooled; each donor is one independent sample.)')
    if known.stats_only:
        # [REV v4] iterate on the pooled tests without redrawing the figures.
        ran_stats = _run_ploidy_stats(known)
        _auto_write_latex_table(known, ran_stats)
        return 0
    # [REV] split off the statistical-test flags, then delegate the rest
    ret = _base_main(rest + ['-o', known.output])
    if known.stats and (ret is None or int(ret) == 0):
        ran_stats = _run_ploidy_stats(known)
        _auto_write_latex_table(known, ran_stats)
    return ret


if __name__ == "__main__":
    sys.exit(main())
