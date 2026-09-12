import argparse
import logging
import os
import sys

from multiprocessing import Pool

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')  # Set non-interactive backend before importing pyplot

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
import seaborn as sns
from matplotlib.patches import Patch, Rectangle, PathPatch

# [REV] Statistical tests (Friedman omnibus + two-sided Wilcoxon signed-rank
# post-hoc with Holm correction, effect sizes, BCa bootstrap CIs).  The sibling
# module is found because Python puts this script's directory on sys.path.
try:
    import stat_tests
except ImportError:  # pragma: no cover - only if the file tree is broken
    stat_tests = None

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(filename)s %(levelname)s %(message)s')


# --------------------------------------------------------------------------- #
# [REV] LaTeX export of the pairwise statistical table                        #
# --------------------------------------------------------------------------- #
# The statistical tests write <output>.stats.pairwise.tsv (see README.md).
# This block turns the important part of that file -- the pairwise rows of the
# two ground-truth scenarios (scenario in {Hap_0, Hap_1}) restricted to the
# columns (scenario, metric, caller_a, caller_b, pvalue_holm) -- into a
# copy-paste-ready booktabs LaTeX table.  [REV v3] The table is pivoted: each
# row is one (metric, caller_b) pair and the two numeric columns hold the
# Holm-adjusted P values for Hap_0 and Hap_1 respectively, so the two
# ground-truth scenarios can be compared side by side.  With --stats-all-pairs
# a Caller-$a$ column is added to identify the pair.
# [REV] By default the same table is ALSO written to
# <output>.stats.pairwise.tex after every successful stats run
# (--no-latex-table disables this); --latex-table additionally prints it to
# stdout and exits without reading stdin.

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


def _perf_legend(ref_desc, ref_is_constant_column=True):
    """Table legend for the CNV-calling performance pairwise table."""
    ref_clause = (F'({ref_desc})'
                  if ref_is_constant_column else
                  F'($a$; see the Caller $a$ column, {ref_desc} by default)')
    return (
        'Pairwise comparison of CNV-calling performance between the reference caller '
        F'{ref_clause} and each other caller ($b$). '
        'Hap\\_0 (haploidy-assumed): the ground-truth CNs of the near-haploid cells are '
        'assumed to be one-valued vectors (CN = 1 across the whole genome); Hap\\_1 '
        '(aneuploidy-aware): the ground-truth CNs are the CNs called by the same caller '
        'from the pre-simulated data (Fig.~1a). $p$ values are two-sided Wilcoxon '
        'signed-rank tests on per-donor medians of the paired per-cell differences '
        '(reference vs.\\ caller $b$), Holm--Bonferroni-adjusted within each (scenario, '
        'metric) family; bold values are significant at the 0.05 family-wise level. '
        'CN, copy number.'
    )


def _tex_escape(s):
    """Escape LaTeX special characters in a table text cell."""
    return (str(s)
            .replace('\\', r'\textbackslash{}')
            .replace('&', r'\&')
            .replace('%', r'\%')
            .replace('_', r'\_')
            .replace('#', r'\#')
            .replace('$', r'\$')
            .replace('^', r'\^{}'))


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


def _latex_table_lines(tsv_path, reference='ginkgo', legend_kind='perf',
                       table_label='tab:pairwise', alpha=0.05):
    """Build the booktabs LaTeX table from a *.stats.pairwise.tsv file.

    [REV v3] Only the pairwise rows of the two ground-truth scenarios
    (scenario in {Hap_0, Hap_1}) are kept.  The table is pivoted so that each
    row is one (metric, caller_b) pair and the two numeric columns hold the
    Holm-adjusted P values for Hap_0 and Hap_1 respectively.  With
    --stats-all-pairs the compared pairs are not all referenced to one caller,
    so a Caller-$a$ column is added to identify pairs.  Returns the list of
    table lines, or None on any problem (a reason is logged).
    """
    if not os.path.isfile(tsv_path):
        logging.error('pairwise stats file not found: %s', tsv_path)
        logging.error('run the script with the statistical tests enabled (default) to generate it first')
        return None
    tab = pd.read_csv(tsv_path, sep='\t')
    for col in ('scenario', 'metric', 'caller_b', 'pvalue_holm'):
        if col not in tab.columns:
            logging.error('column %r missing from %s (available: %s)',
                          col, tsv_path, ', '.join(map(str, tab.columns)))
            return None
    sub = tab.loc[tab['scenario'].isin(('Hap_0', 'Hap_1'))].copy()
    if sub.empty:
        logging.error('no rows with scenario in {Hap_0, Hap_1} in %s', tsv_path)
        return None
    sub = sub.dropna(subset=['scenario', 'metric', 'caller_b'])
    # Row order: metrics and compared callers keep their order of first
    # appearance in the Hap_0 rows (as in the main figures); fall back to the
    # full file if a scenario is missing.
    hap0 = sub.loc[sub['scenario'] == 'Hap_0']
    metric_order = list(dict.fromkeys(hap0['metric'].tolist())) or list(dict.fromkeys(sub['metric'].tolist()))
    caller_order = list(dict.fromkeys(hap0['caller_b'].tolist())) or list(dict.fromkeys(sub['caller_b'].tolist()))
    # With one fixed reference the Caller $a$ column is redundant (it goes
    # into the caption); with --stats-all-pairs it is needed to identify pairs.
    ref_callers = list(dict.fromkeys(sub['caller_a'].tolist())) if 'caller_a' in sub.columns else []
    show_caller_a = len(ref_callers) > 1
    # ---- pivot: one row per (metric, [caller_a,] caller_b), cols = scenario ----
    index_cols = ['metric', 'caller_a', 'caller_b'] if show_caller_a else ['metric', 'caller_b']
    piv = sub.pivot_table(index=index_cols, columns='scenario',
                          values='pvalue_holm', aggfunc='first')
    for sc in ('Hap_0', 'Hap_1'):
        if sc not in piv.columns:
            piv[sc] = float('nan')
    piv = piv[['Hap_0', 'Hap_1']].reset_index()
    piv['metric'] = pd.Categorical(piv['metric'], categories=metric_order, ordered=True)
    piv['caller_b'] = pd.Categorical(piv['caller_b'], categories=caller_order, ordered=True)
    if show_caller_a:
        caller_a_order = list(dict.fromkeys(sub['caller_a'].tolist()))
        piv['caller_a'] = pd.Categorical(piv['caller_a'], categories=caller_a_order, ordered=True)
        piv = piv.sort_values(['metric', 'caller_a', 'caller_b'])
    else:
        piv = piv.sort_values(['metric', 'caller_b'])
    ref_desc = _format_reference(reference)
    legend = _perf_legend(ref_desc, ref_is_constant_column=not show_caller_a)
    if show_caller_a:
        header = ('Metric & Caller $a$ & Caller $b$ '
                  '& Hap\\_0 Holm $p$ & Hap\\_1 Holm $p$ \\\\')
        colspec = 'llrrr'
    else:
        header = ('Metric & Caller $b$ '
                  '& Hap\\_0 Holm $p$ & Hap\\_1 Holm $p$ \\\\')
        colspec = 'llrr'
    lines = [
        F'% LaTeX table generated by {os.path.basename(sys.argv[0])} '
        '(requires \\usepackage{booktabs})',
        '\\begin{table}[htbp]',
        '  \\centering',
        F'  \\caption{{{legend}}}',
        F'  \\label{{{table_label}}}',
        F'  \\begin{{tabular}}{{{colspec}}}',
        '    \\toprule',
        F'    {header}',
        '    \\midrule',
    ]
    for rec in piv.itertuples(index=False):
        p_hap0 = _fmt_pvalue_holm(rec.Hap_0, alpha=alpha)
        p_hap1 = _fmt_pvalue_holm(rec.Hap_1, alpha=alpha)
        # the cell type a metric applies to is part of its display name (e.g.
        # 'PCC_intCN (aneuploid cells)'), matching the figure labels
        metric_label = F'{rec.metric} ({metric_cell_type_label(rec.metric)})'
        if show_caller_a:
            lines.append(F'    {_tex_escape(metric_label)} '
                         F'& {_tex_escape(caller_display_name_tex(rec.caller_a))} '
                         F'& {_tex_escape(caller_display_name_tex(rec.caller_b))} '
                         F'& {p_hap0} & {p_hap1} \\\\')
        else:
            lines.append(F'    {_tex_escape(metric_label)} '
                         F'& {_tex_escape(caller_display_name_tex(rec.caller_b))} '
                         F'& {p_hap0} & {p_hap1} \\\\')
    lines += [
        '    \\bottomrule',
        '  \\end{tabular}',
        '\\end{table}',
    ]
    return lines


def emit_latex_stats_table(tsv_path, reference='ginkgo', legend_kind='perf',
                           table_label='tab:pairwise', alpha=0.05):
    """Print a copy-paste-ready booktabs LaTeX table from a *.stats.pairwise.tsv file.

    See _latex_table_lines for the table contents.  Returns 0 on success,
    1 on any problem.
    """
    lines = _latex_table_lines(tsv_path, reference=reference,
                               legend_kind=legend_kind,
                               table_label=table_label, alpha=alpha)
    if lines is None:
        return 1
    print('\n'.join(lines))
    return 0


def write_latex_stats_table(tsv_path, tex_path, reference='ginkgo', legend_kind='perf',
                            table_label='tab:pairwise', alpha=0.05):
    """[REV] Write the booktabs LaTeX table to `tex_path` (default pipeline output).

    Same table as emit_latex_stats_table, but into a file instead of stdout.
    Returns 0 on success, 1 on any problem (logged; never raises).
    """
    try:
        lines = _latex_table_lines(tsv_path, reference=reference,
                                   legend_kind=legend_kind,
                                   table_label=table_label, alpha=alpha)
        if lines is None:
            return 1
        with open(tex_path, 'w') as fh:
            fh.write('\n'.join(lines) + '\n')
        logging.info('LaTeX pairwise table written to %s', tex_path)
        return 0
    except Exception as exc:  # pragma: no cover - never break the pipeline on cosmetics
        logging.warning('could not write the LaTeX table %s: %s', tex_path, exc)
        return 1


parser1 = argparse.ArgumentParser()
parser1.add_argument('-t', '--type', type=int, default=0, help='Output type. 0: all features. 1: testing features. 2: only plot the main fig. ')
parser1.add_argument('-o', '--output', default='scWGS-performances')
# [REV] statistical-test options (see bench_results/stat_tests.py)
parser1.add_argument('--no-stats', dest='stats', action='store_false', default=True,
                    help='Skip the statistical tests (default: run them and write '
                         '<output>.stats.{pairwise,friedman,concordance}.tsv + .stats.json).')
parser1.add_argument('--stats-only', action='store_true', default=False,
                    help='Run only the statistical tests; exit before any figure is drawn.')
parser1.add_argument('--stats-reference', default='ginkgo', metavar='CALLER',
                    help='Reference caller for the pairwise tests (default: ginkgo).')
parser1.add_argument('--stats-all-pairs', action='store_true', default=False,
                    help='Test all caller pairs instead of reference vs. every other caller.')
parser1.add_argument('--stats-boot', type=int, default=10000, metavar='N',
                    help='Bootstrap resamples for the BCa CIs (default: 10000).')
parser1.add_argument('--stats-seed', type=int, default=1, metavar='SEED',
                    help='Seed of the bootstrap RNG (default: 1).')
parser1.add_argument('--stats-alpha', type=float, default=0.05, metavar='ALPHA',
                    help='Family-wise alpha for Holm rejection (default: 0.05).')
parser1.add_argument('--stats-pair-key', default=None, metavar='COLS',
                    help='Comma-separated columns identifying one simulated cell '
                         '(default: accession_1,accession_2,cellLine,overall_ploidy,CNA_percent).')
# [REV v3] cluster (independent-experimental-unit) options: per-cell results of
# the same caller are correlated within donor groups (all ~1,989 simulated cells
# derive from only nine donors), so inference is aggregated to the per-donor
# level by default: each donor is one effective sample, and the per-cell
# differences within a donor are aggregated to a median before the Wilcoxon
# signed-rank test.  This avoids the pseudoreplication that arises from treating
# every cell as independent (see stat_tests.py for the full rationale and the
# ICC / design-effect diagnostics).
# Finer keys (e.g. accession_1,accession_2,cellLine) are available as a
# sensitivity analysis; pass "none" to revert to the naive per-cell tests.
parser1.add_argument('--stats-cluster-key', default='donor', metavar='COLS',
                    help='Comma-separated columns defining the independent experimental '
                         'unit (cluster) for the statistical tests. Default: donor '
                         '(each human donor is one effective sample; per-cell results '
                         'of the same donor are aggregated to a median before testing). '
                         'Finer sensitivity analysis: '
                         '--stats-cluster-key accession_1,accession_2,cellLine '
                         '(shared-material cells; in the current long TSV this '
                         'key is finer than the donor-level default). '
                         'Pass "none" to revert to the naive per-cell tests that treat '
                         'every cell as independent (discouraged: pseudoreplication).')
parser1.add_argument('--stats-no-cluster', action='store_true', default=False,
                    help='Alias of --stats-cluster-key none (discouraged: '
                         'pseudoreplication; per-cell results of the same caller are '
                         'correlated within shared-material groups).')
# [REV] The LaTeX pairwise table is now generated BY DEFAULT (written to
# <output>.stats.pairwise.tex after every successful stats run);
# --no-latex-table opts out and --latex-table prints it to stdout and exits.
parser1.add_argument('--no-latex-table', dest='latex_table_auto',
                    action='store_false', default=True,
                    help='Do not write <output>.stats.pairwise.tex after the '
                         'statistical tests (default: write it).')
parser1.add_argument('--latex-table', action='store_true', default=False,
                    help='Print the booktabs LaTeX table built from the existing '
                         '<output>.stats.pairwise.tsv (rows with scenario in '
                         '{Hap_0, Hap_1}; pivoted to one row per (metric, caller_b) '
                         'with two Holm-adjusted $p$ columns, Hap_0 and Hap_1) '
                         'to stdout and exit, before reading stdin.')

args = parser1.parse_args()

if args.latex_table:
    sys.exit(emit_latex_stats_table(
        args.output + '.stats.pairwise.tsv',
        reference=args.stats_reference,
        legend_kind='perf',
        table_label='tab:scwgs-perf-pairwise',
        alpha=args.stats_alpha))

# Column (caller) metadata: the manuscript name, the publication year shown under the
# name, and the exact publication date used to sort the columns (journal issue date
# where the journal assigns one, otherwise the online first-publication date; Europe
# PMC, verified 2026-09-12).  The DOI of each paper is listed for traceability:
#   HMMcopy    10.1093/bioinformatics/btl238   2006-07-01
#   Copynumber 10.1186/1471-2164-13-591        2012-11-04
#   Ginkgo     10.1038/nmeth.3578              2015-09-07
#   AneuFinder 10.1186/s13059-016-0971-7       2016-05-31
#   SCCNV      10.3389/fgene.2020.505441       2020-11-16
#   CHISEL     10.1038/s41587-020-0661-6       2021-02-01 (issue; online 2020-09-02)
#   SCYN       10.1186/s12864-021-07941-3      2021-11-16
#   SeCNV      10.1093/bib/bbac264             2022-07-01
#   FLCNA      10.1101/gr.278098.123           2024-02-07
CALLER_PUBLICATION = {
    'hmmcopy'   : ('HMMcopy',    2006, '2006-07-01'),
    'copynumber': ('Copynumber', 2012, '2012-11-04'),
    'ginkgo'    : ('Ginkgo',     2015, '2015-09-07'),
    'aneufinder': ('AneuFinder', 2016, '2016-05-31'),
    'sccnv'     : ('SCCNV',      2020, '2020-11-16'),
    'chisel'    : ('CHISEL',     2021, '2021-02-01'),
    'scyn'      : ('SCYN',       2021, '2021-11-16'),
    'secnv'     : ('SeCNV',      2022, '2022-07-01'),
    'flcna'     : ('FLCNA',      2024, '2024-02-07'),
}


def caller_display_name(caller):
    """Caller column label: manuscript name plus its exact publication year, e.g.
    'Ginkgo\\n2015'."""
    name, year, _ = CALLER_PUBLICATION.get(str(caller), (str(caller), None, '9999-99-99'))
    return F'{name}\n{year}' if year else name


def caller_sort_key(caller):
    """Chronological column order (exact publication date); unknown callers last."""
    name, _, date = CALLER_PUBLICATION.get(str(caller), (str(caller), None, '9999-99-99'))
    return (date, name)


def caller_display_name_tex(caller):
    """LaTeX-table label: manuscript name plus publication year, e.g. 'Ginkgo 2015'."""
    name, year, _ = CALLER_PUBLICATION.get(str(caller), (str(caller), None, ''))
    return F'{name} {year}'.strip()

logscale_features = [
        'with_aneuploidy_aware_gametes.obs2exp_ploidy_ratio',
        'with_haploidy_assumed_gametes.obs2exp_ploidy_ratio',
        'with_aneuploidy_aware_gametes.expected_ploidy',
        'with_haploidy_assumed_gametes.expected_ploidy',
]

CONTINUOUS_FEATURES_NAME2DESC = {
    # benchmarking-strategy-dependent
    'with_aneuploidy_aware_gametes.obs2exp_ploidy_ratio': 'Ratio of the observed (called) ploidy to the expected (ground-truth) ploidy of each simulated cell (Fig. 1a, path Hap_1). '  # [REV]
            '\nThe ground-truth CNs of the near-haploid cells are the CNs called by the same caller from the pre-simulated data. ',
    'with_haploidy_assumed_gametes.obs2exp_ploidy_ratio': 'Ratio of the observed (called) ploidy to the expected (ground-truth) ploidy of each simulated cell (Fig. 1a, path Hap_0). '  # [REV]
            '\nThe ground-truth CNs of the near-haploid cells are assumed to be one-valued vectors (i.e., CN = 1 across the whole genome). ',
    'with_aneuploidy_aware_gametes.expected_ploidy': 'Expected (ground-truth) ploidy of each simulated cell (Fig. 1a, path Hap_1). '  # [REV]
            '\nThe ground-truth CNs of the near-haploid cells are the CNs called by the same caller from the pre-simulated data. ',
    'with_haploidy_assumed_gametes.expected_ploidy': 'Expected (ground-truth) ploidy of each simulated cell (Fig. 1a, path Hap_0). '  # [REV]
            '\nThe ground-truth CNs of the near-haploid cells are assumed to be one-valued vectors (i.e., CN = 1 across the whole genome). ',
    # sample-dependent
    'average_seq_depth': 'Average sequencing depth of each simulated cell',  # [REV]
    'raw_total_sequences': 'Total number of sequenced reads of each simulated cell',  # [REV]
    'bases_mapped_cigar': 'Total number of bases mapped to the reference genome (CIGAR-aware) for each simulated cell',  # [REV]
    'reads_mapped': 'Total number of reads mapped to the reference genome for each simulated cell',  # [REV]
    # simulation dependent
    'CNA_percent': 'Percentage of the genome affected by the simulated copy-number alterations (CNAs)',  # [REV]
    # result-dependent
    'observed_ploidy': 'Ploidy of each simulated cell, as observed (called) by the copy-number caller',  # [REV]
    'bed_1_cn0_genome_size': 'Number of base pairs with copy number (CN) = 0 called by the same caller from the first of the two merged near-haploid samples',  # [REV]
    'bed_1_cn1_genome_size': 'Number of base pairs with copy number (CN) = 1 called by the same caller from the first of the two merged near-haploid samples',  # [REV]
    'bed_1_cn2plus_genome_size': 'Number of base pairs with copy number (CN) > 1 called by the same caller from the first of the two merged near-haploid samples',  # [REV]
    'bed_2_cn0_genome_size': 'Number of base pairs with copy number (CN) = 0 called by the same caller from the second of the two merged near-haploid samples',  # [REV]
    'bed_2_cn1_genome_size': 'Number of base pairs with copy number (CN) = 1 called by the same caller from the second of the two merged near-haploid samples',  # [REV]
    'bed_2_cn2plus_genome_size': 'Number of base pairs with copy number (CN) > 1 called by the same caller from the second of the two merged near-haploid samples',  # [REV]
}

continuous_features = list(CONTINUOUS_FEATURES_NAME2DESC.keys())

CATEGORICAL_FEATURES_NAME2DESC = {
    # sample-dependent
    'donor': 'The human donor from whom the near-haploid cells were derived',  # [REV]
    'sampleType' : 'Cell type of the near-haploid cells (e.g., sperm, polar body, or female pronucleus)',  # [REV]
    'avgSpotLen' : 'Average spot length (i.e., sequencing read length), which depends on the single-cell sequencing technology', # a few unique values  # [REV]
    # simulation-dependent
    'overall_ploidy' : 'Overall ploidy of the simulated cells, either diploid or aneuploid', # diploid or aneuploid  # [REV]
    'cellLine' : 'Cancer cell line (e.g., COLO-829, HCC1395, or HeLa) whose copy-number profile is emulated by the simulation',  # [REV]
    'n_samples_mixed' : 'Number of near-haploid samples merged to simulate each cell (a technical detail)', # 1 or 2, as computed below from accession_1 vs accession_2  # [REV]
}

categorical_features = list(CATEGORICAL_FEATURES_NAME2DESC.keys())

FEATURES_NAME2DESC = (CONTINUOUS_FEATURES_NAME2DESC | CATEGORICAL_FEATURES_NAME2DESC)

# --------------------------------------------------------------------------- #
# Metric IDs, scopes and order                                                #
# --------------------------------------------------------------------------- #
# The pipeline's historical metric IDs are renamed to the manuscript IDs when the
# long TSV is read, so an already-generated TSV keeps working (and a TSV that already
# uses the new IDs is accepted unchanged).  Refresh the mapping in
# cnv_bedset_to_consistency.py / cnv_gather_results.py when the data are regenerated.
METRIC_ID_RENAMES = {
    'accuracy'       : 'intCN_accuracy',
    'PCC_intCN'      : 'intCN_PCC',
    'PCC_nonintCN'   : 'nonintCN_PCC',
    'frac_cov_genome': 'CN_genome_cov_frac',
    'frac_modal_cn'  : 'intCN_modal_frac',
    'frac_mod_cn^N'  : 'intCN_modal_frac',   # transitional name from an earlier revision
}


def rename_metric_columns(df):
    """Apply METRIC_ID_RENAMES to every '<scenario>.<metric>' column."""
    renames = {}
    for col in df.columns:
        if '.' not in col:
            continue
        prefix, metric = col.rsplit('.', 1)
        new_metric = METRIC_ID_RENAMES.get(metric)
        if new_metric and new_metric not in df.columns:
            renames[col] = F'{prefix}.{new_metric}'
    if renames:
        shown = ', '.join(F'{k} -> {v}' for k, v in sorted(renames.items())[:4])
        logging.info('renaming %d metric column(s): %s%s', len(renames), shown,
                     ' ...' if len(renames) > 4 else '')
    return df.rename(columns=renames)


# Rows are grouped by what they measure and, within a group, by importance:
# integer-CN agreement first, then the correlation metrics, then the breakpoint
# metrics, and finally the two "call-quality" fractions (modal-CN fraction, then
# covered-genome fraction), which are the last two rows.
METRIC_ORDER = [
    'intCN_accuracy',
    'intCN_PCC',
    'nonintCN_PCC',
    'breakpoint_f1score',
    'breakpoint_precision',
    'breakpoint_recall',
    'intCN_modal_frac',
    'CN_genome_cov_frac',
]

# Cell type(s) each metric is evaluated on.  intCN_accuracy and CN_genome_cov_frac
# are meaningful for every simulated cell, and their figures/statistics therefore keep
# all cells.  The correlation and breakpoint metrics are restricted to the ANEUPLOID
# simulations: the normal (diploid) simulations have a constant ground-truth total CN
# (CN=2), so they have no ground-truth CN transitions (the Hap_0 breakpoint ratios are
# undefined there) and their CN profile is dominated by allelic/copy-neutral structure
# that these metrics cannot score.  intCN_modal_frac is by definition a property of the
# normal (non-tumor, diploid) simulations.
METRIC_CELL_TYPES = {
    'intCN_accuracy'      : 'all',
    'intCN_PCC'           : 'aneuploid',
    'nonintCN_PCC'        : 'aneuploid',
    'CN_genome_cov_frac'  : 'all',
    'intCN_modal_frac'    : 'diploid',
    'breakpoint_f1score'  : 'aneuploid',
    'breakpoint_precision': 'aneuploid',
    'breakpoint_recall'   : 'aneuploid',
}
METRIC_CELL_TYPE_LABELS = {'all': 'all cells', 'diploid': 'diploid cells', 'aneuploid': 'aneuploid cells'}

# intCN_modal_frac uses only the observed (called) CN profile, so it does not depend on
# the ground-truth scenario; it is displayed once instead of once per scenario (the
# legend carries an 'N/A' key for it).
SCENARIO_INDEPENDENT_PERF_METRICS = frozenset({'intCN_modal_frac'})
SCENARIO_INDEPENDENT_TAG = 'observed calls'


def metric_cell_type(metric):
    return METRIC_CELL_TYPES.get(metric, 'all')


def metric_cell_type_label(metric):
    return METRIC_CELL_TYPE_LABELS[metric_cell_type(metric)]


def metric_display_name(metric, sep='\n'):
    """Metric ID plus the cell type it is evaluated on, e.g.
    'breakpoint_recall\\n(aneuploid cells)'."""
    return F'{metric}{sep}({metric_cell_type_label(metric)})'


def metric_scenario_order(metric, scenario_order):
    """Scenarios plotted for a metric: both ground-truth scenarios, or the single
    'observed calls' series for scenario-independent metrics."""
    return [SCENARIO_INDEPENDENT_TAG] if metric in SCENARIO_INDEPENDENT_PERF_METRICS else scenario_order


def restrict_metrics_to_cell_types(df):
    """NaN out every metric/scenario column outside the cell type it applies to, so
    that both the figures and the statistical tests use that metric's own cell subset
    (all / aneuploid / diploid)."""
    if 'overall_ploidy' not in df.columns:
        logging.warning('column overall_ploidy is missing: metric cell-type restrictions are skipped')
        return df
    ploidy = df['overall_ploidy'].astype(str)
    for metric, cell_type in METRIC_CELL_TYPES.items():
        if cell_type == 'all':
            continue
        mask = ploidy != cell_type
        for col in [c for c in df.columns if c.endswith('.' + metric)]:
            df.loc[mask, col] = np.nan
    return df


if (args.type & 0x1):
    continuous_features = [continuous_features[0]]
    categorical_features = [categorical_features[0]]

df = pd.read_csv(sys.stdin, sep='\t')
df = rename_metric_columns(df)
sortby_columns = (['Caller'] + [x for x in (categorical_features + continuous_features) if x in df.columns])
df = df.sort_values(by=sortby_columns)
df = restrict_metrics_to_cell_types(df)
df['n_samples_mixed'] = np.where(df['accession_1'] == df['accession_2'], 1, 2)
the_df = df.copy()
the_callers = set(df['Caller'].unique())
# columns (callers) are ordered by the exact publication date of each tool
caller_and_its_df_iterable = sorted(df.groupby('Caller'), key=lambda kv: caller_sort_key(kv[0]))

# --------------------------------------------------------------------------- #
# Plot configuration                                                          #
# --------------------------------------------------------------------------- #

SHOW_MEDIAN_LABELS = True    # print the median of every box in a row just above the maximal-performance line
SHOW_ABBREV_FOOTNOTE = True  # one-line abbreviation note at the very bottom of each figure
SKIP_KDEPLOTS = True         # [FIX] the 2-D KDE layer carries so much vector graphics that the
                             # supplementary PDF becomes enormous; keep the CODE_v2 behaviour
                             # (KDE plots skipped) but make it an explicit, documented knob.

# Verbose definitions: NOT drawn inside the figures anymore. Paste them into the figure caption.
THE_PERF_METRIC_NAME2DESC = {
    'intCN_accuracy': 'Accuracy (Acc) of the observed (called) versus expected (ground-truth) integer copy numbers (CNs)',
    'CN_genome_cov_frac': 'Fraction of the human reference genome hg19 covered by the observed (called) copy-number profile',
    'intCN_modal_frac': 'Fraction of the CNV-call-covered genome (base-pair weighted) that is assigned to the modal (most frequent) observed integer copy-number state; in the normal (non-tumor, diploid) simulations the mode is CN=2 unless the caller\'s ploidy estimate is off. The metric uses the observed calls only, so it is identical under the Hap_0 and Hap_1 ground-truth derivations and is shown once (legend key N/A)',
    'intCN_PCC': 'Pearson correlation coefficient (PCC) of the observed (called) versus expected (ground-truth) integer copy numbers (CNs)',
    'nonintCN_PCC': 'Pearson correlation coefficient (PCC) of the observed (called) non-integer copy numbers versus the expected (ground-truth) integer copy numbers (CNs)',
    'breakpoint_f1score': 'F1-score of detecting copy-number changes (breakpoints), balancing breakpoint precision and recall',
    'breakpoint_precision': 'Breakpoint precision: a copy-number transition called by the caller is a true positive if it is matched one-to-one (closest pairs first) with a ground-truth copy-number transition within 200 kb',  # [REV]
    'breakpoint_recall': 'Breakpoint recall: a ground-truth copy-number transition is a true positive if it is matched one-to-one (closest pairs first) with a copy-number transition called by the caller within 200 kb',  # [REV]
}
assert set(METRIC_ORDER) == set(THE_PERF_METRIC_NAME2DESC), \
    'METRIC_ORDER and THE_PERF_METRIC_NAME2DESC disagree'
the_perf_metrics = list(METRIC_ORDER)   # rows are grouped and importance-sorted
THE_PERF_METRIC_NAME2SHORT = {  # short names for the caption; the figures show the IDs only
    'intCN_accuracy': 'Accuracy (Acc)',
    'CN_genome_cov_frac': 'Genome coverage',
    'intCN_modal_frac': 'Fraction at the modal CN',
    'intCN_PCC': 'PCC of integer CNs',
    'nonintCN_PCC': 'PCC of non-integer CNs',
    'breakpoint_f1score': 'Breakpoint F1-score',
    'breakpoint_precision': 'Breakpoint precision',
    'breakpoint_recall': 'Breakpoint recall',
}

aneu_gametes_perf_metrics = [f'with_aneuploidy_aware_gametes.{m}' for m in the_perf_metrics]
hapl_gametes_perf_metrics = [f'with_haploidy_assumed_gametes.{m}' for m in the_perf_metrics]

gamete_type_to_perf_metrics = {
    'aneuploidy_aware_gametes': aneu_gametes_perf_metrics,
    'haploidy_assumed_gametes': hapl_gametes_perf_metrics,
}

gamete_type2desc = {  # verbose: caption only
    'with_haploidy_assumed_gametes': 'Haploidy-assumed (path Hap_0): the ground-truth CNs of the near-haploid cells are assumed to be one-valued vectors (i.e., CN = 1 across the whole genome)',  # [REV]
    'with_aneuploidy_aware_gametes': 'Aneuploidy-aware (path Hap_1): the ground-truth CNs of the near-haploid cells are the CNs called by the same caller from the pre-simulated data',  # [REV]
}

# [FIX] The figures carry only these short tags; their meaning is stated exactly once,
# in the one-line note at the bottom of each figure (and in the caption).
gamete_type2short = {
    'with_haploidy_assumed_gametes': 'Hap_0',
    'with_aneuploidy_aware_gametes': 'Hap_1',
}
the_gamete_legend_title = 'Ground-truth derivation'
THE_ABBREV_NOTE = (
    'Hap_0: haploidy-assumed ground truth; Hap_1: aneuploidy-aware ground truth (Fig. 1a). '
    'CN, copy number; PCC, Pearson correlation coefficient.'
)
THE_ABBREV_NOTE_MEDIANS = THE_ABBREV_NOTE + '\n' + (
    ' Coloured numbers just above the dashed maximal-performance line (typically 1) give each box\'s median, coloured by ground-truth scenario.'
)
THE_ABBREV_NOTE_MEDIANS_SUPP = THE_ABBREV_NOTE + '\n' + (
    ' Dashed horizontal lines mark the median performance of each ground-truth scenario.'
)

# Scenario colours: the boxes use seaborn's 'colorblind' palette in gamete_type2short
# order; the median numbers/lines reuse the same hues, darkened so they stay readable.
SCENARIO_ORDER = list(gamete_type2short.values())
_scenario_palette = sns.color_palette('colorblind', len(SCENARIO_ORDER))
SCENARIO_TEXT_COLORS = {s: tuple(0.55 * ch for ch in _scenario_palette[i]) for i, s in enumerate(SCENARIO_ORDER)}
SCENARIO_LINE_COLORS = {s: tuple(0.75 * ch for ch in _scenario_palette[i]) for i, s in enumerate(SCENARIO_ORDER)}
# Neutral colour and legend key for the scenario-independent metrics
# (intCN_modal_frac): they use the observed calls only and are shown once.
SCENARIO_NA_COLOR = '0.55'
SCENARIO_NA_LABEL = 'N/A'

# [STYLE] Bold, centered overall figure title, mirroring the reference
# metric-by-method grid figure ('scRNA-seq CNV caller performance across
# datasets, methods, and metrics').
THE_FIG_TITLE = 'scWGS CNV caller performance across simulated cells, callers, and metrics'


# --------------------------------------------------------------------------- #
# [REV v3] Statistical tests (donor-level, cluster-robust)                    #
# --------------------------------------------------------------------------- #
# Every caller is evaluated on the SAME simulated cells, so per-cell metrics are
# paired (blocked) by cell: this handles the correlation ACROSS CALLERS within a
# cell.  What a block design additionally assumes is that the BLOCKS (cells) are
# mutually independent - i.e. that the many per-cell results produced by the SAME
# caller are independent.  That assumption is questionable here: all cells are
# downsamplings of the haplotype-normalized BAMs of only NINE donors scored
# against three COSMIC templates, so per-cell differences are correlated within
# donor groups (and the ITH deletions are nested across CNA_percent).  Treating
# ~1,989 cells as ~1,989 independent observations is pseudoreplication: the
# design effect 1 + (m-1)*ICC inflates the test statistics and shrinks the P
# values (see bench_results/test_stat_tests.py, demo_independence_failure, for a
# measured demonstration).
#
# [REV v3] The tests therefore run at the DONOR level by default (each donor is
# one effective sample): per-donor medians of the per-cell differences are the
# units of the Wilcoxon signed-rank and exact sign tests; Friedman on per-donor
# caller medians; Holm on donor-level P; cluster bootstrap CIs; ICC /
# design-effect / effective-n diagnostics per comparison.  The per-cell
# quantities remain in the outputs as descriptive statistics and as flagged
# naive comparisons (pvalue_cell_naive).  Finer cluster keys
# (--stats-cluster-key accession_1,accession_2,cellLine) are available as a
# sensitivity analysis.
# Exact P values, effect sizes and CIs land in:
#   <output>.stats.pairwise.tsv, <output>.stats.friedman.tsv,
#   <output>.stats.concordance.tsv, <output>.stats.json
if args.stats:
    if stat_tests is None:
        logging.warning('stat_tests.py not found next to this script: statistical tests skipped')
    else:
        _stats_prefix = args.output
        if args.stats_no_cluster or (args.stats_cluster_key or '').strip().lower() in ('none', 'naive', 'off'):
            _stats_cluster_key = []          # naive per-cell mode (discouraged)
        elif args.stats_cluster_key:
            _stats_cluster_key = [c.strip() for c in args.stats_cluster_key.split(',') if c.strip()]
        else:
            _stats_cluster_key = None        # module default: donor
        logging.info('running statistical tests (reference caller: %s, cluster key: %s) ...',
                     args.stats_reference,
                     'none (NAIVE per-cell)' if _stats_cluster_key == []
                     else (_stats_cluster_key or stat_tests.DEFAULT_CLUSTER_KEY))
        _stats_settings = stat_tests.run_caller_benchmark_stats(
            the_df, _stats_prefix,
            perf_metrics=the_perf_metrics,
            gamete_type2short=gamete_type2short,
            reference=args.stats_reference,
            all_pairs=args.stats_all_pairs,
            pair_key_cols=(args.stats_pair_key.split(',') if args.stats_pair_key else None),
            cluster_key_cols=_stats_cluster_key,
            cluster_agg='median',
            n_resamples=args.stats_boot,
            seed=args.stats_seed,
            alpha=args.stats_alpha)
        # [REV] The LaTeX pairwise table is generated BY DEFAULT after every
        # successful stats run (also under --stats-only), so the numbers behind
        # the figures and the manuscript table can never drift apart.
        if _stats_settings is not None and args.latex_table_auto:
            write_latex_stats_table(
                args.output + '.stats.pairwise.tsv',
                args.output + '.stats.pairwise.tex',
                reference=args.stats_reference,
                legend_kind='perf',
                table_label='tab:scwgs-perf-pairwise',
                alpha=args.stats_alpha)
        if args.stats_only:
            logging.info('--stats-only: exiting before the figures')
            sys.exit(0)


# --------------------------------------------------------------------------- #
# Shared plotting helpers                                                     #
# --------------------------------------------------------------------------- #

def make_tidy_perf_df():
    """Tidy long dataframe: one row per (caller, cell, metric, ground-truth scenario).
    Scenario-independent metrics contribute a single 'observed calls' series."""
    callers, metrics, scenarios, vals = [], [], [], []
    for metric in the_perf_metrics:
        scenario_independent = metric in SCENARIO_INDEPENDENT_PERF_METRICS
        gamete_types = ['with_haploidy_assumed_gametes'] if scenario_independent else list(gamete_type2short)
        for gamete_type in gamete_types:
            scenario = SCENARIO_INDEPENDENT_TAG if scenario_independent else gamete_type2short[gamete_type]
            callers.extend(list(the_df['Caller']))
            metrics.extend([metric] * len(the_df))
            scenarios.extend([scenario] * len(the_df))
            vals.extend(list(the_df[(gamete_type + '.' + metric)]))
    tidy = pd.DataFrame({
        'Caller': [caller_display_name(c) for c in callers],
        'Metric': metrics,
        'Scenario': scenarios,
        'Performance': vals,
    })
    caller_order = [caller_display_name(c) for c in sorted(dict.fromkeys(callers), key=caller_sort_key)]
    scenario_order = list(gamete_type2short.values())
    return tidy, caller_order, scenario_order


def boxplot_with_style(ax, **kw):
    """sns.boxplot; degrades gracefully if the installed seaborn rejects flierprops.

    NOTE: do NOT set rasterized=True on the fliers. Seaborn >= 0.13 draws every flier
    as its own Line2D artist, and the PDF backend turns each rasterized artist into a
    separate 300-dpi image XObject -- with hundreds of fliers this silently eats
    gigabytes of RAM (observed OOM kills at the PDF-save step).
    """
    kw.setdefault('saturation', 1)  # keep the box colors identical to the legend swatches
    try:
        return sns.boxplot(ax=ax, flierprops=dict(markersize=2, alpha=0.5), **kw)
    except TypeError:
        return sns.boxplot(ax=ax, **kw)


def metric_row_ylim(metric, vals):
    """[STYLE] Fixed, canonical y window shared by every panel of a metric row,
    mirroring the reference metric-by-method grid figure: correlation metrics
    (PCC) are pinned to (-1, 1), bounded-score metrics to (0, 1), and anything
    else falls back to a data-driven window with 10% padding. A shared window
    per row keeps panels visually comparable across callers."""
    vals = pd.Series(vals).dropna()
    if 'PCC' in metric:  # correlation metrics -> canonical full range
        return (-1.0, 1.0)
    if metric in BOUNDED_PERF_METRIC_RANGES:
        return BOUNDED_PERF_METRIC_RANGES[metric]
    if vals.empty:
        return (0.0, 1.0)
    lo, hi = float(vals.min()), float(vals.max())
    pad = max(0.02, (hi - lo) * 0.1)
    return (lo - pad, hi + pad)


# Fixed, canonical ranges of the bounded performance metrics (used by metric_row_ylim).
BOUNDED_PERF_METRIC_RANGES = {
    'intCN_accuracy': (0.0, 1.0),
    'CN_genome_cov_frac': (0.0, 1.0),
    'intCN_modal_frac': (0.0, 1.0),
    'breakpoint_f1score': (0.0, 1.0),
    'breakpoint_precision': (0.0, 1.0),
    'breakpoint_recall': (0.0, 1.0),
}


def swarm_grid_canvas(n_rows, n_cols, panel_w=1.1, min_panel_h=1.0, max_panel_h=4.0,
                      legend_h=1.3, label_pad_left=1.8, label_pad_top=0.62):
    """[STYLE] Square canvas with fixed absolute-inch zones, reproducing the
    layout machinery of the reference metric-by-method grid figure: row labels
    on the left (label_pad_left), column labels on top (label_pad_top), a
    shared legend strip at the bottom (legend_h), a FIXED panel width, and a
    panel height that stretches (within [min_panel_h, max_panel_h]) to fill
    whatever vertical room is left.

    [FIX] Zone heights are sized to SNUGLY fit their contents so no large
    blank bands appear between the grid and the figure title above it, or
    between the grid and the boxed legend below it: label_pad_top fits the
    column titles plus the row of scenario-median numbers that sit just
    above the first-row panels, and legend_h fits the boxed titled legend
    plus the two-line abbreviation note.  [FIX] legend_h was raised from 1.0
    in to 1.3 in: the boxed legend has to sit high enough for the note under
    it to clear its lower frame edge, while still leaving a gap to the
    bottom-row caller labels of the grid.  Returns (fig, gridspec,
    canvas_side_inch).

    NOTE: figures built on this canvas must be saved WITHOUT
    bbox_inches='tight' -- tight-bbox cropping would re-fit the output to the
    content's own (non-square) bounding box and defeat the explicit layout."""
    grid_w = panel_w * n_cols
    width_driven_side = label_pad_left + grid_w
    avail_h = width_driven_side - label_pad_top - legend_h
    panel_h = avail_h / n_rows if n_rows > 0 else min_panel_h
    panel_h = max(min_panel_h, min(max_panel_h, panel_h))
    grid_h = panel_h * n_rows
    side = max(width_driven_side, label_pad_top + grid_h + legend_h)
    fig = plt.figure(figsize=(side, side))
    gs = fig.add_gridspec(n_rows, n_cols,
                          left=label_pad_left / side,
                          right=min((label_pad_left + grid_w) / side, 0.995),
                          top=min((legend_h + grid_h) / side, 0.98),
                          bottom=legend_h / side,
                          wspace=0.15, hspace=0.25)
    return fig, gs, side


def annotate_box_medians(ax, plot_df, x_order, hue_order, fontsize=6, decimals=3, y_ref=1.0):
    """[tie-breaking] Thicken the median of every box and print all medians in a single
    row just above the dashed maximal-performance line (y_ref, typically 1), colour-coded
    by ground-truth scenario, so that almost-tied performances stay legible and directly
    comparable (the numbers are aligned instead of floating at unequal box heights).
    Implemented on top of the box patches, so it works with both the old and the new
    seaborn boxplot back-ends."""
    if not SHOW_MEDIAN_LABELS:
        return
    meds = plot_df.groupby(['Caller', 'Scenario'])['Performance'].median().to_dict()
    ax.axhline(y_ref, color='0.55', linewidth=0.8, linestyle=(0, (5, 3)), zorder=1)
    # seaborn >= 0.13 draws every box as a PathPatch (and adds invisible full-axes
    # Rectangles that merely carry the hue labels for the legend); older versions
    # draw plain Rectangles instead. Use whichever of the two exists -- never mix.
    path_boxes = [p for p in ax.patches if isinstance(p, PathPatch)]
    if path_boxes:
        box_x0w = []
        for p in path_boxes:
            vs = np.asarray(p.get_path().vertices)
            if vs.size == 0:
                continue
            x0, x1 = float(vs[:, 0].min()), float(vs[:, 0].max())
            box_x0w.append((x0, x1 - x0))
    else:
        box_x0w = [(p.get_x(), p.get_width()) for p in ax.patches if isinstance(p, Rectangle)]
    # group the (dodged) boxes per x category and map them onto hue_order left->right
    by_x = {}
    for x0, w in box_x0w:
        by_x.setdefault(int(round(x0 + w / 2.0)), []).append((x0, w))
    lo, hi = ax.get_ylim()
    off = 0.015 * (hi - lo)
    for xi, box_list in by_x.items():
        if not (0 <= xi < len(x_order)) or len(box_list) != len(hue_order):
            continue
        for (x0, w), hue_label in zip(sorted(box_list), hue_order):
            med = meds.get((x_order[xi], hue_label))
            if med is None or pd.isna(med):
                continue
            ax.plot([x0, x0 + w], [med, med], color='0.1', linewidth=2.0,
                    solid_capstyle='butt', zorder=4)
            ax.text(x0 + w / 2.0, y_ref + off, F'{med:.{decimals}f}',
                    ha='center', va='bottom', fontsize=fontsize,
                    color=SCENARIO_TEXT_COLORS.get(hue_label, '0.25'), zorder=5, clip_on=False)


def add_fig_title(fig):
    """[STYLE] Bold, horizontally centered overall title at the very top of the
    figure, as in the reference metric-by-method grid figure."""
    fig.suptitle(THE_FIG_TITLE, fontsize=12, fontweight='bold', y=0.995)


def add_bottom_legend(fig, handles, side, ncol, fontsize=8.5):
    """[STYLE] Shared legend strip, horizontally centered and anchored a small
    margin above the absolute figure bottom edge (inside the legend_h zone
    reserved by swarm_grid_canvas). [FIX] The legend now carries its name
    ('Ground-truth derivation') as the legend title and is drawn in a thin
    box, with the title centred over the entries.
    [FIX] The anchor was raised from 0.27 in to 0.46 in: the abbreviation note
    printed below the legend (add_bottom_note) keeps its place, and the legend box
    used to touch -- and slightly overlap -- that note. The legend zone reserved by
    swarm_grid_canvas (1.3 in) has room for the raised legend plus the note without
    touching the bottom-row caller labels of the grid."""
    leg = fig.legend(handles=handles, loc='lower center', ncol=ncol,
                     fontsize=fontsize, title=the_gamete_legend_title,
                     title_fontsize=fontsize, frameon=True, fancybox=False,
                     framealpha=1.0, edgecolor='0.65', borderpad=0.5,
                     bbox_to_anchor=(0.5, 0.46 / side),
                     columnspacing=1.6, handlelength=1.4, handletextpad=0.5)
    try:  # centre the title over the entries (matplotlib keeps this private)
        leg._legend_box.align = 'center'
    except (AttributeError, TypeError):
        pass
    return leg


def add_bottom_note(fig, note, side, fontsize=7.5):
    """[STYLE] One-line abbreviation note tucked under the legend strip, still
    inside the reserved legend zone."""
    if SHOW_ABBREV_FOOTNOTE:
        fig.text(0.5, 0.10 / side, note, ha='center', va='bottom',
                 fontsize=fontsize, color='0.25')


def scenario_na_handle():
    """Legend key for scenario-independent metrics (intCN_modal_frac): the metric does
    not depend on the ground-truth derivation, so it is shown once in neutral grey."""
    return Patch(facecolor=SCENARIO_NA_COLOR, edgecolor='0.2', linewidth=0.8, label=SCENARIO_NA_LABEL)


def gamete_legend_handles(include_na=True):
    palette = sns.color_palette('colorblind', len(gamete_type2short))
    handles = [Patch(facecolor=c, edgecolor='0.2', linewidth=0.8, label=lab)
               for c, lab in zip(palette, gamete_type2short.values())]
    if include_na:
        handles.append(scenario_na_handle())
    return handles


def caller_legend_handles(caller_order):
    palette = sns.color_palette('colorblind', len(caller_order))
    return [Patch(facecolor=c, edgecolor='0.2', linewidth=0.8, label=lab)
            for c, lab in zip(palette, caller_order)]


# --------------------------------------------------------------------------- #
# Main figures                                                                #
# --------------------------------------------------------------------------- #

def plot_multirow_main():
    """[STYLE] Single-column variant of the reference metric-by-method grid
    figure: one row per metric, all callers along x, and one boxplot pair
    (Hap_0 vs Hap_1) per caller instead of the reference's swarmplot. Metric
    IDs become left-side row labels (first column only), caller names sit
    below the bottom row, and the shared legend + abbreviation note move to a
    bottom strip, all inside a square reference-style canvas."""
    tidy, caller_order, scenario_order = make_tidy_perf_df()
    n_metrics = len(the_perf_metrics)
    # One wide panel per metric row: the panel width accommodates all callers.
    fig, gs, side = swarm_grid_canvas(n_metrics, 1, panel_w=1.1 * len(caller_order))
    for r, metric in enumerate(the_perf_metrics):
        ax = fig.add_subplot(gs[r, 0])
        sub_df = tidy[tidy['Metric'] == metric]
        hue_order = metric_scenario_order(metric, scenario_order)
        if sub_df['Performance'].notna().any():
            if metric in SCENARIO_INDEPENDENT_PERF_METRICS:
                # single series: no ground-truth scenario to compare for this metric
                boxplot_with_style(ax, data=sub_df, x='Caller', y='Performance',
                                   order=caller_order, color=SCENARIO_NA_COLOR, linewidth=1.0)
            else:
                boxplot_with_style(ax, data=sub_df, x='Caller', y='Performance',
                                   hue='Scenario', order=caller_order, hue_order=hue_order,
                                   palette='colorblind', linewidth=1.0)
        if ax.legend_ is not None:
            ax.legend_.remove()
        row_ylim = metric_row_ylim(metric, sub_df['Performance'])
        ax.set_ylim(*row_ylim)
        # Medians stay (tie-breaking); the dashed maximal-performance line sits at 1.
        annotate_box_medians(ax, sub_df, caller_order, hue_order, fontsize=6, y_ref=1.0)
        # [STYLE] metric ID + its cell-type scope as the left-side row label
        ax.set_ylabel(metric_display_name(metric), fontsize=8, rotation=0, ha='right', va='center', labelpad=8)
        ax.set_xlabel('')  # kill seaborn's auto 'Caller' x label on every row
        ax.tick_params(axis='y', labelsize=6)
        if r == n_metrics - 1:
            ax.tick_params(axis='x', labelsize=8.5)
            ax.set_xticklabels([c.replace('_', '\n') for c in caller_order])
        else:
            ax.tick_params(axis='x', labelbottom=False, length=0)
        for spine in ('top', 'right'):
            ax.spines[spine].set_visible(False)
    add_fig_title(fig)
    add_bottom_legend(fig, gamete_legend_handles(), side, ncol=3)
    add_bottom_note(fig, THE_ABBREV_NOTE_MEDIANS, side)
    try:
        fig.tight_layout()  # same call sequence as the reference figure
    except (ValueError, TypeError):
        pass
    plt.savefig(args.output + '_multirow_main.pdf', dpi=300)
    plt.savefig(args.output + '_multirow_main.png', dpi=300)
    plt.close()


def plot_grid_main():
    """[STYLE] The reference metric-by-method grid figure, 1:1: a square canvas
    with fixed-inch label/legend zones, one narrow panel per (metric row,
    caller column), panels of a row sharing its canonical y window, caller
    names above the first row AND below the last row, metric IDs as left-side
    row labels on the first column only, a bold overall title, and a shared
    boxed, titled legend at the bottom. Each panel holds one boxplot pair
    (Hap_0 vs Hap_1) instead of the reference's swarmplot."""
    tidy, caller_order, scenario_order = make_tidy_perf_df()
    n_metrics, n_callers = len(the_perf_metrics), len(caller_order)
    fig, gs, side = swarm_grid_canvas(n_metrics, n_callers)
    for r, metric in enumerate(the_perf_metrics):
        row_df = tidy[tidy['Metric'] == metric]
        row_ylim = metric_row_ylim(metric, row_df['Performance'])  # canonical window, shared per row
        hue_order = metric_scenario_order(metric, scenario_order)
        scenario_independent = metric in SCENARIO_INDEPENDENT_PERF_METRICS
        for c, caller in enumerate(caller_order):
            ax = fig.add_subplot(gs[r, c])
            sub_df = row_df[row_df['Caller'] == caller]
            if not sub_df.empty and sub_df['Performance'].notna().any():
                if scenario_independent:
                    boxplot_with_style(ax, data=sub_df, x='Caller', y='Performance',
                                       order=[caller], color=SCENARIO_NA_COLOR, linewidth=0.9)
                else:
                    boxplot_with_style(ax, data=sub_df, x='Caller', y='Performance',
                                       hue='Scenario', order=[caller], hue_order=hue_order,
                                       palette='colorblind', linewidth=0.9)
                if ax.legend_ is not None:
                    ax.legend_.remove()
            ax.set_ylim(*row_ylim)
            if not sub_df.empty:
                annotate_box_medians(ax, sub_df, [caller], hue_order,
                                     fontsize=5.5, y_ref=1.0)
            # [STYLE] panel geometry exactly as in the reference grid
            ax.set_xlim(-0.6, 0.6)
            ax.set_xticks([])
            ax.set_xlabel('')
            if c == 0:
                ax.set_ylabel(metric_display_name(metric), fontsize=8, rotation=0, ha='right', va='center', labelpad=8)
            else:
                ax.set_ylabel('')
                ax.set_yticklabels([])
            ax.tick_params(axis='y', labelsize=6)
            if r == 0:
                # [FIX] the scenario-median numbers sit just above the top edge
                # of the first-row panels (between them and the dashed
                # maximal-performance line at y_ref); pad the titles so the
                # two never overlap.
                ax.set_title(caller.replace('_', '\n'), fontsize=8.5, pad=10)
            elif r == n_metrics - 1:
                ax.set_xlabel(caller.replace('_', '\n'), fontsize=8.5)
            for spine in ('top', 'right'):
                ax.spines[spine].set_visible(False)
    add_fig_title(fig)
    add_bottom_legend(fig, gamete_legend_handles(), side, ncol=3)
    add_bottom_note(fig, THE_ABBREV_NOTE_MEDIANS, side)
    try:
        fig.tight_layout()  # same call sequence as the reference figure
    except (ValueError, TypeError):
        pass
    plt.savefig(args.output + '_final_grid.pdf', dpi=300)
    plt.savefig(args.output + '_final_grid.png', dpi=300)
    plt.close()


def plot_main():
    # [STYLE] the metric column labels carry the metric ID, the cell-type scope and the
    # ground-truth scenario, so the canvas is slightly wider and the labels slightly
    # smaller than before to keep 15 adjacent categories from overlapping.
    fig, ax = plt.subplots(figsize=(17.5, 7.5))
    xcats, callers, vals = [], [], []
    for metric in the_perf_metrics:
        scenario_independent = metric in SCENARIO_INDEPENDENT_PERF_METRICS
        gamete_types = ['with_haploidy_assumed_gametes'] if scenario_independent else list(gamete_type2short)
        for gamete_type in gamete_types:
            # [STYLE] horizontal category label: metric + cell-type scope (+ scenario),
            # like the reference figure's bottom column labels ('method\nvariant')
            xcat = (metric_display_name(metric) if scenario_independent
                    else F'{metric_display_name(metric)}\n({gamete_type2short[gamete_type]})')
            xcats.extend([xcat] * len(the_df))
            callers.extend(list(the_df['Caller']))
            vals.extend(list(the_df[(gamete_type + '.' + metric)]))
    dfm = pd.DataFrame({
        'Metric': xcats,
        'Caller': [caller_display_name(c) for c in callers],
        'Performance': vals,
    })
    xcat_order = []
    for m in the_perf_metrics:
        if m in SCENARIO_INDEPENDENT_PERF_METRICS:
            xcat_order.append(metric_display_name(m))
        else:
            xcat_order.extend(F'{metric_display_name(m)}\n({gamete_type2short[gt]})' for gt in gamete_type2short)
    dfm['Metric'] = pd.Categorical(dfm['Metric'], categories=xcat_order, ordered=True)
    dfm = dfm.dropna(subset=['Performance'])  # metrics restricted to one cell type have no point elsewhere
    caller_order = [caller_display_name(c) for c in sorted(dict.fromkeys(callers), key=caller_sort_key)]
    boxplot_with_style(ax, data=dfm, x='Metric', y='Performance', hue='Caller',
                       hue_order=caller_order, palette='colorblind', linewidth=0.9)
    if ax.legend_ is not None:
        ax.legend_.remove()
    # NOTE: per-box median numbers are intentionally omitted here (14 x 9 boxes is too dense);
    # the multirow and grid main figures carry all the medians instead.
    ax.set_xlabel('')
    ax.set_ylabel('Performances', fontsize=10)
    ax.tick_params(axis='y', labelsize=8)
    ax.tick_params(axis='x', labelsize=6.5)
    for spine in ('top', 'right'):
        ax.spines[spine].set_visible(False)
    # [STYLE] bold overall title + shared legend strip at the BOTTOM, as in the
    # reference figure (callers laid out over two rows, 5 + 4); margins are set
    # manually so the legend and the abbreviation note get a dedicated strip
    # at the bottom, mirroring the reference figure's fixed-zone layout.
    # [FIX] bottom margin raised from 0.15 to 0.20: the legend moved up (see below), so
    # the panel x tick labels need to stay clear of its upper frame edge as well.
    # [FIX] bottom margin raised from 0.20 to 0.24: the x categories now carry three
    # lines (metric ID / cell-type scope / ground-truth scenario), so the labels need
    # to stay clear of the boxed caller legend below them.
    fig.subplots_adjust(left=0.07, right=0.985, top=0.915, bottom=0.24)
    add_fig_title(fig)
    caller_handles = caller_legend_handles(caller_order)
    fig.legend(handles=caller_handles, loc='lower center', ncol=5,
               fontsize=8.5, frameon=True, fancybox=False,
               framealpha=1.0, edgecolor='0.65', borderpad=0.5,
               # [FIX] raised from 0.022 so the abbreviation note below cannot touch
               # the legend frame (same overlap as in the final-grid figure).
               bbox_to_anchor=(0.5, 0.045),
               columnspacing=1.8, handlelength=1.4, handletextpad=0.5)
    if SHOW_ABBREV_FOOTNOTE:
        fig.text(0.5, 0.008, THE_ABBREV_NOTE, ha='center', va='bottom',
                 fontsize=7.5, color='0.25')
    plt.savefig(args.output + '_main.pdf')
    plt.savefig(args.output + '_main.png', dpi=300)
    plt.close()


# --------------------------------------------------------------------------- #
# Supplementary figures                                                       #
# --------------------------------------------------------------------------- #

def plot_onepage(page_args): # (continuous_features, categorical_features, the_perf_metrics):
    # [FIX] the argument used to be called `args`, shadowing the global argparse
    # namespace of the same name -- a maintenance trap.  Renamed to page_args.
    feature, page_num = page_args
    if feature in logscale_features: feat_scale = 'log'
    else: feat_scale = ''
    logging.info(F'START plotting {feature} with scale={feat_scale}')
    if feature in continuous_features: assert feature not in categorical_features, F'The feature {feature} cannot be both continous and categorical'
    fig1 = plt.figure(figsize=(30, 20), constrained_layout=True)
    # top strip: centered legend + one-line abbreviation note. [FIX] the strip
    # is tall enough that the BOXED legend (its title + entries) clears the
    # note lines below it -- with the old frameless legend the overlap existed
    # but was invisible; an opaque box would hide the note text.
    gs = gridspec.GridSpec(1+len(the_perf_metrics), len(the_callers),
            height_ratios=[8.5]+[10]*len(the_perf_metrics), figure=fig1, wspace=0, hspace=0.1)
    legend_ax = fig1.add_subplot(gs[0,:])
    legend_ax.set_axis_off()
    if feature in continuous_features:
        # [FIX] NaN-safe data-driven windows: plain min()/max() return
        # order-dependent results in the presence of NaNs.
        feat_vals = pd.to_numeric(df[feature], errors='coerce').dropna()
        feat_min = float(feat_vals.min()) if len(feat_vals) else 0.0
        feat_max = float(feat_vals.max()) if len(feat_vals) else 1.0
        plot_feat_min = feat_min - (feat_max - feat_min) * 0.05
        plot_feat_max = feat_max + (feat_max - feat_min) * 0.05
        if feat_min > 0: plot_feat_minmax = feat_min * (1-0.05)
        else: plot_feat_minmax = -1e99
        if plot_feat_min < plot_feat_minmax: plot_feat_min = plot_feat_minmax
    handles, labels = [], []
    for rowidx, perf_metric in enumerate(the_perf_metrics):
        scenario_independent = perf_metric in SCENARIO_INDEPENDENT_PERF_METRICS
        gamete_types_here = ['with_haploidy_assumed_gametes'] if scenario_independent else list(gamete_type2short)
        feature_all_perf_vals = []
        for gamete_type in gamete_types_here:
            feature_all_perf_vals += pd.to_numeric(
                df[gamete_type + '.' + perf_metric], errors='coerce').dropna().tolist()
        if feature_all_perf_vals:
            min_perf_val = min(feature_all_perf_vals)
            max_perf_val = max(feature_all_perf_vals)
        else:
            min_perf_val, max_perf_val = 0.0, 1.0
        plot_perf_min = min_perf_val - (max_perf_val - min_perf_val) * 0.05
        plot_perf_max = max_perf_val + (max_perf_val - min_perf_val) * 0.05
        for colidx, (caller, caller_df) in enumerate(caller_and_its_df_iterable):
            plot_dfs = []
            for gamete_type in gamete_types_here:
                gamete_perf_metric = gamete_type+'.'+perf_metric
                # x: feature; y: performance metric
                plot_df = caller_df[[feature]].copy()
                plot_df[perf_metric] = caller_df[gamete_perf_metric]
                plot_df['gamete_type'] = (SCENARIO_INDEPENDENT_TAG if scenario_independent
                                          else gamete_type2short[gamete_type])  # short tag only
                plot_dfs.append(plot_df)
            plot_df = pd.concat(plot_dfs).reset_index(drop=True)
            ax2 = fig1.add_subplot(gs[rowidx+1, colidx])
            # a metric restricted to one cell type has no data for the other cells,
            # so panels can legitimately be empty; seaborn's boxplot/stripplot/
            # scatterplot crash on an all-NaN panel, so skip drawing it
            if plot_df[perf_metric].notna().any():
                if feature in categorical_features:
                    if feature == 'donor':  # [REV] shorten the *values*
                        plot_df[feature] = plot_df[feature].replace('345HS1', 'HS1') # prevent cluttering of words for the donor categorical variable
                    if scenario_independent:
                        plot_ret = sns.stripplot(data=plot_df, x=feature, y=perf_metric, ax=ax2,
                                                 color=SCENARIO_NA_COLOR, alpha=0.125, rasterized=True)
                    else:
                        plot_ret = sns.stripplot (data=plot_df, x=feature, y=perf_metric, hue='gamete_type', ax=ax2, palette='colorblind', alpha=0.125, rasterized=True)
                else:
                    logging.info(F'plotting {perf_metric} versus {feature} for {caller}')
                    if scenario_independent:
                        plot_ret = sns.scatterplot(data=plot_df, x=feature, y=perf_metric, ax=ax2,
                                                   color=SCENARIO_NA_COLOR, alpha=0.125, rasterized=True)
                    else:
                        plot_ret = sns.scatterplot(data=plot_df, x=feature, y=perf_metric, hue='gamete_type', style='gamete_type', ax=ax2, palette='colorblind', alpha=0.125, markers=['x', '+'], rasterized=True)
                    # [FIX] explicit, documented skip: the KDE layer carries so much
                    # vector graphics that the multipage PDF becomes enormous.
                    # (CODE_v2 initialised this flag to 1, which silently skipped
                    # *every* KDE plot; the behaviour is kept, but now it is named,
                    # and skipping it for data reasons is still reported.)
                    skip_kdeplot = SKIP_KDEPLOTS
                    if not skip_kdeplot:
                        for _gt, plot_df_2 in plot_df.groupby('gamete_type'):
                            if len(set(plot_df_2[perf_metric])) == 1:
                                skip_kdeplot = True
                                logging.warning(F'KDEplot of {perf_metric} versus {feature} for {caller} skipped: the {_gt} values are constant')
                    if skip_kdeplot:
                        logging.debug(F'Skip the KDEplot of {perf_metric} versus {feature} for {caller}')
                    else:
                        if scenario_independent:
                            sns.kdeplot(data=plot_df, x=feature, y=perf_metric, ax=ax2, color=SCENARIO_NA_COLOR,
                                        levels=10, fill=True, alpha=0.5, legend=False)
                        else:
                            sns.kdeplot(data=plot_df, x=feature, y=perf_metric, hue='gamete_type', ax=ax2, palette='colorblind', levels=10, fill=True, alpha=0.5, legend=False)
                # [NEW] one dashed reference line per scenario at that scenario's median performance,
                # so that almost-tied callers can still be told apart in the supplementary pages
                for scenario_short in ([SCENARIO_INDEPENDENT_TAG] if scenario_independent else SCENARIO_ORDER):
                    the_median = plot_df.loc[plot_df['gamete_type'] == scenario_short, perf_metric].median()
                    if pd.notna(the_median):
                        ax2.axhline(the_median, color=SCENARIO_LINE_COLORS.get(scenario_short, '0.55'), linewidth=1.2,
                                    linestyle=(0, (5, 3)), alpha=0.9, zorder=3)
                _handles, _labels = plot_ret.get_legend_handles_labels()
                if _labels:
                    handles, labels = _handles, _labels
                if plot_ret.legend_ is not None:
                    plot_ret.legend_.remove()
            if feat_scale:
                ax2.set_xscale(feat_scale)
            if feature in continuous_features:
                ax2.set_xlim(plot_feat_min, plot_feat_max)
            ax2.set_ylim(plot_perf_min, plot_perf_max)
            if rowidx == 0:
                # [FIX] single-line title: the caller name only. Caller names
                # appear exactly ONCE per page (top row); the columns stay
                # readable and the panels below carry no repeated names.
                ax2.set_title(caller_display_name(caller), fontsize=18)
            ax2.set_xlabel('')
            if colidx == 0:
                # [FIX] metric ID + its cell-type scope (same label as in the main figures)
                ax2.set_ylabel(metric_display_name(perf_metric), fontsize=15)
            else:
                ax2.set_ylabel('')
                ax2.tick_params(labelleft=False)
            sns.despine(ax=ax2)
    # [FIX] legend centered at the top of the strip; entries side by side below
    # the title. [STYLE] boxed, like every other legend in the script.  The 'N/A'
    # key marks scenario-independent metrics (intCN_modal_frac), shown once in grey.
    if labels:
        handles, labels = list(handles) + [scenario_na_handle()], list(labels) + [SCENARIO_NA_LABEL]
    else:
        handles = gamete_legend_handles()
        labels = [h.get_label() for h in handles]
    leg = legend_ax.legend(handles, labels,
            title=the_gamete_legend_title,
            loc='upper center', ncol=3, frameon=True, fancybox=False,
            framealpha=1.0, edgecolor='0.65', borderpad=0.6,
            fontsize=16, title_fontsize=16,
            markerscale=3, columnspacing=2.0)
    leg_handles = getattr(leg, 'legend_handles', None)
    if leg_handles is None:  # older matplotlib attribute name
        leg_handles = getattr(leg, 'legendHandles', [])
    for handle in leg_handles:
        handle.set_alpha(1.0)
    if SHOW_ABBREV_FOOTNOTE:  # all the definitions that used to clutter the legend, now one line
        legend_ax.text(0.5, 0.04, THE_ABBREV_NOTE_MEDIANS_SUPP, ha='center', va='bottom',
                       fontsize=11, color='0.25')
    factor = 'Factor: ' + FEATURES_NAME2DESC.get(feature, 'TODO')
    if len(factor) > 200*5:
        fig1.supxlabel(factor, fontsize=16)
    else:
        fig1.supxlabel(factor, fontsize=20)
    fig1.supylabel('Performances', fontsize=24)
    A2Z = [chr(i) for i in range(ord('a'), ord('z') + 1)]  # [REV] lowercase panel letters, matching the a-f style of Fig. 1
    # [FIX] more than 26 supplementary pages would previously raise an IndexError
    panel_letter = A2Z[page_num] if page_num < len(A2Z) else F's{page_num + 1}'
    sublabel_ax = fig1.add_subplot(gs[0,0])
    sublabel_ax.set_axis_off()
    sublabel_ax.set_title(panel_letter, fontsize=30, ha='left', fontweight='bold')

    logging.info(F'END: plotting {feature} with scale={feat_scale}')
    # Detach the figure from this worker process's pyplot registry *before* it is
    # pickled back to the parent. Without this, every finished page would stay alive
    # in the worker for the whole run and the pool would slowly eat gigabytes of RAM.
    plt.close(fig1)
    return fig1


plot_multirow_main()
plot_grid_main()
plot_main()
if (args.type & 0x2): sys.exit(0)

# Create a PDF file to save the pages
the_labels = None
with PdfPages(args.output + '-all.pdf') as pdf:
    n_cores = min([os.cpu_count(), 32])
    with Pool(processes=n_cores) as pool:
        my_map = pool.imap  # imap preserves the page order of the submitted tasks
        for fig1 in my_map(plot_onepage,
                [(feature_withscale, page_num)
                for page_num, feature_withscale in enumerate(continuous_features + categorical_features)]):
            pdf.savefig(fig1, bbox_inches='tight', dpi=75)
            plt.close(fig1)
