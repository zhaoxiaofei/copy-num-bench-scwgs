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

THE_PERF_METRIC_NAME2DESC = {
    'accuracy': 'Accuracy (Acc) of the observed (called) versus expected (ground-truth)\ninteger copy numbers (CNs)',
    'PCC_intCN': 'Pearson correlation coefficient (PCC) of the observed (called) versus expected\n(ground-truth) integer copy numbers (CNs)',
    'PCC_nonintCN': 'Pearson correlation coefficient (PCC) of the observed (called) non-integer copy numbers\nversus the expected (ground-truth) integer copy numbers (CNs)',
    'frac_cov_genome': 'Fraction of the human reference genome hg19 covered by the observed (called)\ncopy-number profile',
    'breakpoint_f1score': 'F1-score of detecting copy-number changes (breakpoints), balancing breakpoint\nprecision and recall',
    'breakpoint_precision': 'Breakpoint precision: an observed (called) breakpoint is a true positive if at least one\nexpected (ground-truth) breakpoint is within 200 kb',  # [REV]
    'breakpoint_recall': 'Breakpoint recall: an expected (ground-truth) breakpoint is a true positive if at least one\nobserved (called) breakpoint is within 200 kb',  # [REV]
}
the_perf_metrics = list(THE_PERF_METRIC_NAME2DESC.keys())
THE_PERF_METRIC_NAME2SHORT = {
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

gamete_type2desc = {
    'with_haploidy_assumed_gametes': 'Haploidy-assumed (path Hap_0): the ground-truth CNs of the near-haploid cells are assumed to be one-valued vectors (i.e., CN = 1 across the whole genome)',  # [REV]
    'with_aneuploidy_aware_gametes': 'Aneuploidy-aware (path Hap_1): the ground-truth CNs of the near-haploid cells are the CNs called by the same caller from the pre-simulated data',  # [REV]
}
gamete_type2short = {
    'with_haploidy_assumed_gametes': 'Haploidy-assumed (path Hap_0)',  # [REV]
    'with_aneuploidy_aware_gametes': 'Aneuploidy-aware (path Hap_1)',  # [REV]
}
the_gamete_legend_title = 'Ground-truth assumption about the copy numbers (CNs) of the near-haploid cells (Fig. 1a)'  # [REV]

# pivot: 'Caller'

def plot_main():
    fig1 = plt.figure(figsize=(2*7, 1*5), constrained_layout=True)
    callers = []
    perf_names = []
    perf_vals = []
    for rowidx, perf_metric in enumerate(the_perf_metrics):
        for gamete_type in gamete_type2desc:
            callers.extend(list(the_df['Caller']))
            perf_names.extend([F'{THE_PERF_METRIC_NAME2SHORT[perf_metric]}\n{gamete_type2short[gamete_type]}'] * len(the_df))  # [REV] parentheses dropped
            perf_vals.extend(list(the_df[(gamete_type+'.'+perf_metric)]))
    df = pd.DataFrame({'Copy-number callers': [caller2desc.get(c, c) for c in callers], 'Performance metrics': perf_names, 'Performance': perf_vals})  # [REV] 'Performances' -> 'Performance'
    plot_ret = sns.boxplot(data=df, x='Performance metrics', y='Performance', hue='Copy-number callers')  # [REV]
    # Find the outlier artists and rasterize them
    #for artist in plot_ret.get_lines():
    #    # In matplotlib, outliers are often 'line' objects with 0 line width
    #    if artist.get_linestyle() == 'None':
    #        pass #artist.set_rasterized(True)
    #plot_ret = sns.stripplot(data=df, x='Performance metrics', y='Performances', hue='Copy-number callers', alpha=0.1, jitter=True,
    #          rasterized=True)

    handles, labels = plot_ret.get_legend_handles_labels()
    plt.tick_params(axis='both', which='major', labelsize=10)

    plot_ret.set_xticklabels(plot_ret.get_xticklabels(), rotation=20, ha='right')

    plt.savefig(args.output + '_main.pdf')
    plt.savefig(args.output + '_main.png', dpi=300)
    plt.close()

# https://www.doubao.com/chat/38421257704964354
def plot_grid_main():
    # --------------------------
    # 1. Data preparation (original logic)
    # --------------------------
    callers = []
    perf_names = []
    perf_vals = []
    scen_names = []
    for rowidx, perf_metric in enumerate(the_perf_metrics):
        for gamete_type in gamete_type2desc:
            callers.extend(list(the_df['Caller']))
            perf_names.extend([perf_metric] * len(the_df))
            scen_names.extend([gamete_type2short[gamete_type]] * len(the_df))
            perf_vals.extend(list(the_df[(gamete_type+'.'+perf_metric)]))

    # Tidy dataframe
    df = pd.DataFrame({
        'Caller': [caller2desc.get(c, c) for c in callers],
        'Metric': perf_names,
        'Scenario': scen_names,
        'Performance': perf_vals
    })

    # --------------------------
    # 2. Create grid: Rows=Metrics, Cols=Callers, one boxplot per gamete scenario in each cell
    # --------------------------
    g = sns.FacetGrid(
        data=df,
        row='Metric',        # Rows = performance metrics (Y)
        col='Caller',        # Columns = callers (X)
        sharex=True,
        sharey=True,
        height=2.2,
        aspect=0.7
    )

    # Plot one boxplot per ground-truth scenario in every grid cell
    g.map_dataframe(
        sns.boxplot,
        x='Caller',
        y='Performance',
        hue='Scenario',
        hue_order=list(gamete_type2short.values()),
        palette='colorblind',
        linewidth=1.0,
        flierprops=dict(markersize=1.5, alpha=0.4)
    )
    g.add_legend(title=the_gamete_legend_title)

    # --------------------------
    # 3. CRITICAL: Show REAL labels (not literals)
    # --------------------------
    # Remove ALL automatic titles and inner labels
    g.set_titles('')
    g.set_axis_labels('', '')

    # --------------------------
    # FIRST COLUMN ONLY: Show REAL METRIC names on Y-axis
    # --------------------------
    for ax, metric in zip(g.axes[:, 0], the_perf_metrics):
        ax.set_ylabel(THE_PERF_METRIC_NAME2SHORT[metric], fontsize=10, labelpad=8)  # REAL metric label

    # --------------------------
    # LAST ROW ONLY: Show REAL CALLER names on X-axis
    # --------------------------
    unique_callers = df['Caller'].unique()
    for ax, caller in zip(g.axes[-1, :], unique_callers):
        caller_name = caller2desc.get(caller, caller)
        ax.set_xlabel(caller_name, fontsize=10, labelpad=8)  # REAL caller label
        # Rotate long caller names if needed
        # ax.tick_params(axis='x', labelrotation=30)
        # Get the existing tick labels and update their properties
        for label in ax.get_xticklabels():
            label.set_rotation(25)       # Tilt at 25 degrees
            label.set_ha('right')        # Anchor at the right edge
            label.set_va('top')          # Anchor at the top edge
    # Hide all tick marks and the redundant per-cell tick labels (clean look)
    for ax in g.axes.flat:
        ax.tick_params(bottom=False, left=False, labelbottom=False)

    # --------------------------
    # 4. Save final figure
    # --------------------------
    plt.tight_layout()
    plt.savefig(args.output + '_final_grid.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(args.output + '_final_grid.png', dpi=300, bbox_inches='tight')
    plt.close()

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

def plot_multirow_main():
    # --------------------------
    # 1. Data preparation (same as original)
    # --------------------------
    callers = []
    perf_names = []
    perf_vals = []
    scen_names = []
    for rowidx, perf_metric in enumerate(the_perf_metrics):
        for gamete_type in gamete_type2desc:
            callers.extend(list(the_df['Caller']))
            perf_names.extend([perf_metric] * len(the_df))
            scen_names.extend([gamete_type2short[gamete_type]] * len(the_df))
            perf_vals.extend(list(the_df[(gamete_type+'.'+perf_metric)]))

    # Tidy dataframe
    df = pd.DataFrame({
        'Caller': [caller2desc.get(c, c) for c in callers],
        'Metric': perf_names,
        'Scenario': scen_names,
        'Performance': perf_vals
    })

    # --------------------------
    # 2. Create figure: 7 rows, 1 column (one per metric)
    # --------------------------
    n_metrics = len(the_perf_metrics)
    fig, axes = plt.subplots(nrows=n_metrics, ncols=1, figsize=(170*0.06, 225*0.06), constrained_layout=True)
    # https://link.springer.com/journal/13073/submission-guidelines

    # If only 1 metric (avoid axis list error)
    if n_metrics == 1:
        axes = [axes]

    # --------------------------
    # 3. Plot one boxplot per metric (each row = one metric, hue = gamete scenario)
    # --------------------------
    for ax, metric in zip(axes, the_perf_metrics):
        # Subset data for THIS metric only
        sub_df = df[df['Metric'] == metric]

        # Plot boxplot: X = Callers, Y = Performance, Hue = Scenario
        sns.boxplot(
            data=sub_df,
            x='Caller',
            y='Performance',
            hue='Scenario',
            hue_order=list(gamete_type2short.values()),
            palette='colorblind',
            ax=ax,
            linewidth=1.2,
            flierprops=dict(markersize=2, alpha=0.5)
        )
        handles, labels = ax.get_legend_handles_labels()
        if ax.legend_ is not None: ax.legend_.remove()  # one figure-level legend below instead

        # --------------------------
        # Style each subplot
        # --------------------------
        ax.set_title(THE_PERF_METRIC_NAME2DESC[metric], fontsize=12, weight='bold')  # TITLE = metric definition
        ax.set_xlabel('')                                 # No repeated x-label
        ax.tick_params(axis='x', labelsize=10, rotation=30)
        ax.tick_params(axis='y', labelsize=10)
        ax.grid(axis='y', alpha=0.3)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=15, ha='right')

    # --------------------------
    # One figure-level legend for the two ground-truth scenarios
    # --------------------------
    fig.legend(handles, labels, title=the_gamete_legend_title, loc='outside lower right', fontsize=10, title_fontsize=10)

    # --------------------------
    # 4. Save final figure
    # --------------------------
    plt.savefig(args.output + '_multirow_main.pdf', dpi=300)
    plt.savefig(args.output + '_multirow_main.png', dpi=300)
    plt.close()

def plot_onepage(args): # (continuous_features, categorical_features, the_perf_metrics):
    feature, page_num = args
    if feature in logscale_features: feat_scale = 'log'
    else: feat_scale = ''
    logging.info(F'START plotting {feature} with scale={feat_scale}')
    if feature in continuous_features: assert feature not in categorical_features, F'The feature {feature} cannot be both continous and categorical'
    fig1 = plt.figure(figsize=(25, 20), constrained_layout=True)
    gs = gridspec.GridSpec(1+len(the_perf_metrics), len(the_callers), height_ratios=[3]+[10]*len(the_perf_metrics), figure=fig1, wspace=0, hspace=0.1)
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
            for gamete_type in gamete_type2desc:
                gamete_perf_metric = gamete_type+'.'+perf_metric
                # x: feature; y: performance metric
                plot_df = caller_df[[feature]].copy()
                plot_df[perf_metric] = caller_df[gamete_perf_metric]
                plot_df['gamete_type'] = gamete_type2desc[gamete_type]
                plot_dfs.append(plot_df)
            plot_df = pd.concat(plot_dfs).reset_index(drop=True)
            #print(plot_df)
            ax2 = fig1.add_subplot(gs[rowidx+1, colidx])
            if feature in categorical_features:
                #plot_ret = sns.boxplot    (data=plot_df, x=feature, y=perf_metric, hue='gamete_type', ax=ax2, whis=np.inf)
                if feature == 'donor':  # [REV] shorten the *values* (the old line renamed column names, so 'HS1' never rendered)
                    plot_df[feature] = plot_df[feature].replace('345HS1', 'HS1') # prevent cluttering of words for the donor categorical variable
                plot_ret = sns.stripplot (data=plot_df, x=feature, y=perf_metric, hue='gamete_type', ax=ax2, alpha=0.125, rasterized=True)
                handles, labels = plot_ret.get_legend_handles_labels()
            else:
                logging.info(F'plotting {perf_metric} versus  {feature} for {caller}')
                #plot_ret = sns.scatterplot(data=plot_df.sample(n=100, random_state=0), x=feature, y=perf_metric, hue='gamete_type', style='gamete_type', ax=ax2)
                plot_ret = sns.scatterplot(data=plot_df, x=feature, y=perf_metric, hue='gamete_type', style='gamete_type', ax=ax2, alpha=0.125, markers=['x', '+'], rasterized=True)
                skip_kdeplot = 1
                for gt, plot_df_2 in plot_df.groupby('gamete_type'):
                    if len(set(plot_df_2[perf_metric])) == 1:
                        skip_kdeplot += 1
                if skip_kdeplot:
                    logging.warning(f'Skip the KDEplot of {perf_metric} versus {feature} for {caller}')
                else:
                    plot_ret2= sns.kdeplot(data=plot_df, x=feature, y=perf_metric, hue='gamete_type', ax=ax2, levels=10, fill=True, alpha=0.5, legend=False)
            handles, labels = plot_ret.get_legend_handles_labels()
            plot_ret.legend_.remove()
            if feat_scale:
                ax2.set_xscale(feat_scale)
            if feature in continuous_features:
                ax2.set_xlim(plot_feat_min, plot_feat_max)
            ax2.set_ylim(plot_perf_min, plot_perf_max)

            if rowidx == 0:
                ax2.set_title('Copy-number caller\n' + caller2desc.get(caller, caller), fontsize=20)  # [REV] was 'Method\n'

            if rowidx == len(the_perf_metrics) - 1:
                #ax2.set_xlabel(feature)
                ax2.set_xlabel('')
            else:
                ax2.set_xlabel('')
                #ax2.set_xticklabels('')
            if colidx == 0:
                ax2.set_ylabel(THE_PERF_METRIC_NAME2SHORT[perf_metric], fontsize=16)
            else:
                ax2.set_ylabel('')
                ax2.set_yticklabels('')
                #sns           .kdeplot    (data=df, x=feature, y=perf_metric, hue='Caller', ax=ax2, palette='husl', alpha=0.5, level=5)
                #ax_histx = fig1.add_subplot(gs[rowidx, colidx].get_gridspec()[rowidx, colidx], sharex=ax2)
                #ax_histy = fig1.add_subplot(gs[rowidx, colidx].get_gridspec()[rowidx, colidx], sharey=ax2)
                #sns.histplot(data=df, x=feature    , ax=ax_histx, kde=True, hue='Caller', palette='husl')
                #sns.histplot(data=df, y=perf_metric, ax=ax_histy, kde=True, hue='Caller', palette='husl')
                #grid = sns.JointGrid(x=feature, y=perf_metric, hue='Caller', data=df, ax=ax2, palette='husl')
                #grid.plot_joint(sns.scatterplot, alpha=0.25)
                #grid.plot_joint(sns.kdeplot, level=5, alpha=0.75)
                #grid.plot_marginals(sns.kdeplot, fill=True, alpha=0.5)
                #handles, labels = grid.ax_joint.get_legend_handles_labels()
            #labels = df['Caller'].unique()
            #plot_ret.ax_joint.legend_.remove()
            #if the_labels:
            #    assert (list(labels) == list(the_labels)), F'{labels} == {the_labels} failed for feature {feature}'
            #else: the_labels = labels
    leg = legend_ax.legend(handles, labels,
            title=the_gamete_legend_title
            #'\n' 'maybe_aneuploid: let the cells CNs be called from the original pre-simulated data; always_haploid: let the cell CNs be set to ones'
            , loc='center', fontsize=20, title_fontsize=20, markerscale=3,
            ncol=1 # ncol=len(labels)
            )
    for handle in leg.legend_handles:
        handle.set_alpha(1.0)
    #print(F'Hanle_label_list={handles},{labels}')
    #legend_ax.legend(handles, labels,
    #        title='Copy-number callers',
    #        loc='center', fontsize=12, title_fontsize=14, ncol=len(labels))
    factor = 'Factor: ' + FEATURES_NAME2DESC.get(feature, 'TODO')
    if len(factor) > 200*5:
        fig1.supxlabel(factor, fontsize=16)
    else:
        fig1.supxlabel(factor, fontsize=20)
    fig1.supylabel('Performance metrics', fontsize=24)
    A2Z = [chr(i) for i in range(ord('a'), ord('z') + 1)]  # [REV] lowercase panel letters, matching the a-f style of Fig. 1
    sublabel_ax = fig1.add_subplot(gs[0,0])
    sublabel_ax.set_axis_off()
    sublabel_ax.set_title(A2Z[page_num], fontsize=30, ha='left', fontweight='bold')
    #fig1.suptitle(A2Z[page_num], fontsize=30, ha='left', fontweight='bold')

    #plt.tight_layout()
    logging.info(F'END: plotting {feature} with scale={feat_scale}')
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
