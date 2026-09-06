import argparse
import logging
import os
import sys

from multiprocessing import Pool
from functools import partial

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')  # Set non-interactive backend before importing pyplot

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
import seaborn as sns
from matplotlib.patches import Patch, Rectangle, PathPatch

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(filename)s %(levelname)s %(message)s')

parser1 = argparse.ArgumentParser()
parser1.add_argument('-t', '--type', type=int, default=0, help='Output type. 0: all features. 1: testing features. 2: only plot the main fig. ')
parser1.add_argument('-o', '--output', default='scWGS-performances')

args = parser1.parse_args()

# The triple-quoted string below maps each caller to its journal and publication year
'''
    'aneufinder': 'AneuFinder  Genome Biology               2016',
    'flcna'     : 'FLCNA       Genome Research              2024',
    'chisel'    : 'Chisel      Nature Biotechnology         2021',
    'copynumber': 'CopyNumber  BMC Genomics                 2012',
    'ginkgo'    : 'Ginkgo      Nature Methods               2015',
    'hmmcopy'   : 'HMMcopy     Bioinformatics               2006',
    'secnv'     : 'SeCNV       Briefings in Bioinformatics  2022',
    'sccnv'     : 'SCCNV       Frontiers in Genetics        2020',
    'scyn'      : 'SCYN/SCOPE  Cell Systems                 2020',
'''

caller2desc = {
    'aneufinder': 'AneuFinder',
    'flcna'     : 'FLCNA',
    'chisel'    : 'Chisel',
    'copynumber': 'Copynumber',
    'ginkgo'    : 'Ginkgo',
    'hmmcopy'   : 'HMMcopy',
    'secnv'     : 'SeCNV',
    'sccnv'     : 'SCCNV',
    'scyn'      : 'SCYN',
}

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
    'n_samples_mixed' : 'Number of near-haploid samples merged to simulate each cell (a technical detail)', # 0 or 1  # [REV]
}

categorical_features = list(CATEGORICAL_FEATURES_NAME2DESC.keys())

FEATURES_NAME2DESC = (CONTINUOUS_FEATURES_NAME2DESC | CATEGORICAL_FEATURES_NAME2DESC)

if (args.type & 0x1):
    continuous_features = [continuous_features[0]]
    categorical_features = [categorical_features[0]]

df = pd.read_csv(sys.stdin, sep='\t')
sortby_columns = (['Caller'] + [x for x in (categorical_features + continuous_features) if x in df.columns])
df = df.sort_values(by=sortby_columns)
df['n_samples_mixed'] = np.where(df['accession_1'] == df['accession_2'], 1, 2)
the_df = df.copy()
the_callers = set(df['Caller'].unique())
caller_and_its_df_iterable = df.groupby('Caller')

# --------------------------------------------------------------------------- #
# Plot configuration                                                          #
# --------------------------------------------------------------------------- #

SHOW_MEDIAN_LABELS = True    # print the median of every box in a row just above the maximal-performance line
SHOW_ABBREV_FOOTNOTE = True  # one-line abbreviation note at the very bottom of each figure

# Verbose definitions: NOT drawn inside the figures anymore. Paste them into the figure caption.
THE_PERF_METRIC_NAME2DESC = {
    'accuracy': 'Accuracy (Acc) of the observed (called) versus expected (ground-truth) integer copy numbers (CNs)',
    'PCC_intCN': 'Pearson correlation coefficient (PCC) of the observed (called) versus expected (ground-truth) integer copy numbers (CNs)',
    'PCC_nonintCN': 'Pearson correlation coefficient (PCC) of the observed (called) non-integer copy numbers versus the expected (ground-truth) integer copy numbers (CNs)',
    'frac_cov_genome': 'Fraction of the human reference genome hg19 covered by the observed (called) copy-number profile',
    'breakpoint_f1score': 'F1-score of detecting copy-number changes (breakpoints), balancing breakpoint precision and recall',
    'breakpoint_precision': 'Breakpoint precision: an observed (called) breakpoint is a true positive if at least one expected (ground-truth) breakpoint is within 200 kb',  # [REV]
    'breakpoint_recall': 'Breakpoint recall: an expected (ground-truth) breakpoint is a true positive if at least one observed (called) breakpoint is within 200 kb',  # [REV]
}
the_perf_metrics = list(THE_PERF_METRIC_NAME2DESC.keys())
THE_PERF_METRIC_NAME2SHORT = {  # short names for the caption; the figures show the IDs only
    'accuracy': 'Accuracy (Acc)',
    'PCC_intCN': 'PCC of integer CNs',
    'PCC_nonintCN': 'PCC of non-integer CNs',
    'frac_cov_genome': 'Genome coverage',
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

# [STYLE] Bold, centered overall figure title, mirroring the reference
# metric-by-method grid figure ('scRNA-seq CNV caller performance across
# datasets, methods, and metrics').
THE_FIG_TITLE = 'scWGS CNV caller performance across simulated cells, callers, and metrics'


# --------------------------------------------------------------------------- #
# Shared plotting helpers                                                     #
# --------------------------------------------------------------------------- #

def make_tidy_perf_df():
    """Tidy long dataframe: one row per (caller, cell, metric, ground-truth scenario)."""
    callers, metrics, scenarios, vals = [], [], [], []
    for metric in the_perf_metrics:
        for gamete_type in gamete_type2short:
            callers.extend(list(the_df['Caller']))
            metrics.extend([metric] * len(the_df))
            scenarios.extend([gamete_type2short[gamete_type]] * len(the_df))
            vals.extend(list(the_df[(gamete_type + '.' + metric)]))
    tidy = pd.DataFrame({
        'Caller': [caller2desc.get(c, c) for c in callers],
        'Metric': metrics,
        'Scenario': scenarios,
        'Performance': vals,
    })
    caller_order = list(dict.fromkeys(tidy['Caller']))   # appearance order (== sorted)
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
    'accuracy': (0.0, 1.0),
    'frac_cov_genome': (0.0, 1.0),
    'breakpoint_f1score': (0.0, 1.0),
    'breakpoint_precision': (0.0, 1.0),
    'breakpoint_recall': (0.0, 1.0),
}


def swarm_grid_canvas(n_rows, n_cols, panel_w=1.1, min_panel_h=1.0, max_panel_h=4.0,
                      legend_h=1.0, label_pad_left=1.8, label_pad_top=0.62):
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
    plus the one-line abbreviation note. Returns (fig, gridspec,
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
    box, with the title centred over the entries."""
    leg = fig.legend(handles=handles, loc='lower center', ncol=ncol,
                     fontsize=fontsize, title=the_gamete_legend_title,
                     title_fontsize=fontsize, frameon=True, fancybox=False,
                     framealpha=1.0, edgecolor='0.65', borderpad=0.5,
                     bbox_to_anchor=(0.5, 0.27 / side),
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


def gamete_legend_handles():
    palette = sns.color_palette('colorblind', len(gamete_type2short))
    return [Patch(facecolor=c, edgecolor='0.2', linewidth=0.8, label=lab)
            for c, lab in zip(palette, gamete_type2short.values())]


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
        boxplot_with_style(ax, data=sub_df, x='Caller', y='Performance',
                           hue='Scenario', order=caller_order, hue_order=scenario_order,
                           palette='colorblind', linewidth=1.0)
        if ax.legend_ is not None:
            ax.legend_.remove()
        row_ylim = metric_row_ylim(metric, sub_df['Performance'])
        ax.set_ylim(*row_ylim)
        # Medians stay (tie-breaking); the dashed maximal-performance line sits at 1.
        annotate_box_medians(ax, sub_df, caller_order, scenario_order, fontsize=6, y_ref=1.0)
        # [STYLE] metric ID as the left-side row label, like the reference rows
        ax.set_ylabel(metric, fontsize=8, rotation=0, ha='right', va='center', labelpad=8)
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
    add_bottom_legend(fig, gamete_legend_handles(), side, ncol=2)
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
        for c, caller in enumerate(caller_order):
            ax = fig.add_subplot(gs[r, c])
            sub_df = row_df[row_df['Caller'] == caller]
            if not sub_df.empty:
                boxplot_with_style(ax, data=sub_df, x='Caller', y='Performance',
                                   hue='Scenario', order=[caller], hue_order=scenario_order,
                                   palette='colorblind', linewidth=0.9)
                if ax.legend_ is not None:
                    ax.legend_.remove()
            ax.set_ylim(*row_ylim)
            if not sub_df.empty:
                annotate_box_medians(ax, sub_df, [caller], scenario_order,
                                     fontsize=5.5, y_ref=1.0)
            # [STYLE] panel geometry exactly as in the reference grid
            ax.set_xlim(-0.6, 0.6)
            ax.set_xticks([])
            ax.set_xlabel('')
            if c == 0:
                ax.set_ylabel(metric, fontsize=8, rotation=0, ha='right', va='center', labelpad=8)
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
    add_bottom_legend(fig, gamete_legend_handles(), side, ncol=2)
    add_bottom_note(fig, THE_ABBREV_NOTE_MEDIANS, side)
    try:
        fig.tight_layout()  # same call sequence as the reference figure
    except (ValueError, TypeError):
        pass
    plt.savefig(args.output + '_final_grid.pdf', dpi=300)
    plt.savefig(args.output + '_final_grid.png', dpi=300)
    plt.close()


def plot_main():
    fig, ax = plt.subplots(figsize=(15.5, 7.5))
    xcats, callers, vals = [], [], []
    for metric in the_perf_metrics:
        for gamete_type in gamete_type2short:
            # [STYLE] two-line horizontal category label, like the reference
            # figure's bottom column labels ('method\nvariant')
            xcat = F'{metric}\n({gamete_type2short[gamete_type]})'
            xcats.extend([xcat] * len(the_df))
            callers.extend(list(the_df['Caller']))
            vals.extend(list(the_df[(gamete_type + '.' + metric)]))
    dfm = pd.DataFrame({
        'Metric': xcats,
        'Caller': [caller2desc.get(c, c) for c in callers],
        'Performance': vals,
    })
    xcat_order = [F'{m}\n({gamete_type2short[gt]})'
                  for m in the_perf_metrics for gt in gamete_type2short]
    dfm['Metric'] = pd.Categorical(dfm['Metric'], categories=xcat_order, ordered=True)
    caller_order = list(dict.fromkeys(dfm['Caller']))
    boxplot_with_style(ax, data=dfm, x='Metric', y='Performance', hue='Caller',
                       hue_order=caller_order, palette='colorblind', linewidth=0.9)
    if ax.legend_ is not None:
        ax.legend_.remove()
    # NOTE: per-box median numbers are intentionally omitted here (14 x 9 boxes is too dense);
    # the multirow and grid main figures carry all the medians instead.
    ax.set_xlabel('')
    ax.set_ylabel('Performances', fontsize=10)
    ax.tick_params(axis='y', labelsize=8)
    ax.tick_params(axis='x', labelsize=8)
    for spine in ('top', 'right'):
        ax.spines[spine].set_visible(False)
    # [STYLE] bold overall title + shared legend strip at the BOTTOM, as in the
    # reference figure (callers laid out over two rows, 5 + 4); margins are set
    # manually so the legend and the abbreviation note get a dedicated strip
    # at the bottom, mirroring the reference figure's fixed-zone layout.
    fig.subplots_adjust(left=0.07, right=0.985, top=0.915, bottom=0.15)
    add_fig_title(fig)
    caller_handles = caller_legend_handles(caller_order)
    fig.legend(handles=caller_handles, loc='lower center', ncol=5,
               fontsize=8.5, frameon=True, fancybox=False,
               framealpha=1.0, edgecolor='0.65', borderpad=0.5,
               bbox_to_anchor=(0.5, 0.022),
               columnspacing=1.8, handlelength=1.4, handletextpad=0.5)
    if SHOW_ABBREV_FOOTNOTE:
        fig.text(0.5, 0.006, THE_ABBREV_NOTE, ha='center', va='bottom',
                 fontsize=7.5, color='0.25')
    plt.savefig(args.output + '_main.pdf')
    plt.savefig(args.output + '_main.png', dpi=300)
    plt.close()


# --------------------------------------------------------------------------- #
# Supplementary figures                                                       #
# --------------------------------------------------------------------------- #

def plot_onepage(args): # (continuous_features, categorical_features, the_perf_metrics):
    feature, page_num = args
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
        feat_min = min(df[feature])
        feat_max = max(df[feature])
        plot_feat_min = feat_min - (feat_max - feat_min) * 0.05
        plot_feat_max = feat_max + (feat_max - feat_min) * 0.05
        if feat_min > 0: plot_feat_minmax = feat_min * (1-0.05)
        else: plot_feat_minmax = -1e99
        if plot_feat_min < plot_feat_minmax: plot_feat_min = plot_feat_minmax
    for rowidx, perf_metric in enumerate(the_perf_metrics):
        feature_all_perf_vals = list(df[('with_aneuploidy_aware_gametes.'+perf_metric)]) + list(df[('with_haploidy_assumed_gametes.'+perf_metric)])
        min_perf_val = min(feature_all_perf_vals)
        max_perf_val = max(feature_all_perf_vals)
        plot_perf_min = min_perf_val - (max_perf_val - min_perf_val) * 0.05
        plot_perf_max = max_perf_val + (max_perf_val - min_perf_val) * 0.05
        for colidx, (caller, caller_df) in enumerate(caller_and_its_df_iterable):
            plot_dfs = []
            for gamete_type in gamete_type2short:
                gamete_perf_metric = gamete_type+'.'+perf_metric
                # x: feature; y: performance metric
                plot_df = caller_df[[feature]].copy()
                plot_df[perf_metric] = caller_df[gamete_perf_metric]
                plot_df['gamete_type'] = gamete_type2short[gamete_type]  # short tag only
                plot_dfs.append(plot_df)
            plot_df = pd.concat(plot_dfs).reset_index(drop=True)
            ax2 = fig1.add_subplot(gs[rowidx+1, colidx])
            if feature in categorical_features:
                if feature == 'donor':  # [REV] shorten the *values*
                    plot_df[feature] = plot_df[feature].replace('345HS1', 'HS1') # prevent cluttering of words for the donor categorical variable
                plot_ret = sns.stripplot (data=plot_df, x=feature, y=perf_metric, hue='gamete_type', ax=ax2, palette='colorblind', alpha=0.125, rasterized=True)
            else:
                logging.info(F'plotting {perf_metric} versus  {feature} for {caller}')
                plot_ret = sns.scatterplot(data=plot_df, x=feature, y=perf_metric, hue='gamete_type', style='gamete_type', ax=ax2, palette='colorblind', alpha=0.125, markers=['x', '+'], rasterized=True)
                skip_kdeplot = 1  # KDE has too much vector graphics in it, resulting in very big PDF # [FIX] CODE_v2 initialised this to 1, which skipped *every* KDE plot
                for gt, plot_df_2 in plot_df.groupby('gamete_type'):
                    if len(set(plot_df_2[perf_metric])) == 1:
                        skip_kdeplot += 1
                if skip_kdeplot:
                    logging.warning(f'Skip the KDEplot of {perf_metric} versus {feature} for {caller}')
                else:
                    plot_ret2= sns.kdeplot(data=plot_df, x=feature, y=perf_metric, hue='gamete_type', ax=ax2, palette='colorblind', levels=10, fill=True, alpha=0.5, legend=False)
            # [NEW] one dashed reference line per scenario at that scenario's median performance,
            # so that almost-tied callers can still be told apart in the supplementary pages
            for scenario_short in SCENARIO_ORDER:
                the_median = plot_df.loc[plot_df['gamete_type'] == scenario_short, perf_metric].median()
                if pd.notna(the_median):
                    ax2.axhline(the_median, color=SCENARIO_LINE_COLORS[scenario_short], linewidth=1.2,
                                linestyle=(0, (5, 3)), alpha=0.9, zorder=3)
            handles, labels = plot_ret.get_legend_handles_labels()
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
                ax2.set_title(caller2desc.get(caller, caller), fontsize=18)
            ax2.set_xlabel('')
            if colidx == 0:
                # [FIX] metric ID only (same ID as in the main figure)
                ax2.set_ylabel(perf_metric, fontsize=15)
            else:
                ax2.set_ylabel('')
                ax2.tick_params(labelleft=False)
            sns.despine(ax=ax2)
    # [FIX] legend centered at the top of the strip; entries side by side below
    # the title. [STYLE] boxed, like every other legend in the script.
    leg = legend_ax.legend(handles, labels,
            title=the_gamete_legend_title,
            loc='upper center', ncol=2, frameon=True, fancybox=False,
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
    sublabel_ax = fig1.add_subplot(gs[0,0])
    sublabel_ax.set_axis_off()
    sublabel_ax.set_title(A2Z[page_num], fontsize=30, ha='left', fontweight='bold')

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
        my_map = pool.imap # map # pool.imap
        for fig1 in my_map(plot_onepage,
                [(feature_withscale, page_num)
                for page_num, feature_withscale in enumerate(continuous_features + categorical_features)]):
            pdf.savefig(fig1, bbox_inches='tight', dpi=75)
            plt.close(fig1)
  
