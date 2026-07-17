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

caller2desc = {
    'chisel':     'Chisel      Nature Biotechnology         2021',
    'copynumber': 'CopyNumber  BMC Genomics                 2012',
    'ginkgo':     'Ginkgo      Nature Methods               2015',
    'hmmcopy':    'HMMcopy     Bioinformatics               2006',
    'secnv':      'SeCNV       Briefings in Bioinformatics  2022',
    'sccnv':      'SCCNV       Frontiers in Genetics        2020',
    'scyn':       'SCYN/SCOPE  Cell Systems                 2020',
}

logscale_features = [
        'with_aneuploidy_aware_gametes.obs2exp_ploidy_ratio',
        'with_haploidy_assumed_gametes.obs2exp_ploidy_ratio',
        'with_aneuploidy_aware_gametes.expected_ploidy',
        'with_haploidy_assumed_gametes.expected_ploidy',
]

CONTINUOUS_FEATURES_NAME2DESC = {
    # benchmarking-strategy-dependent
    'with_aneuploidy_aware_gametes.obs2exp_ploidy_ratio': 'Observed (called) to expected (ground-truth) ploidy ratio. '
            '\nThe ground-truth CNs of near-haploid cells are called by the same caller from the pre-simulated data. ',
    'with_haploidy_assumed_gametes.obs2exp_ploidy_ratio': 'Observed (called) to expected (ground-truth) ploidy ratio. '
            '\nThe ground-truth CNs of near-haploid cells are one-valued vectors (i.e., CN=1 everywhere in the genome). ',
    'with_aneuploidy_aware_gametes.expected_ploidy': 'Expected (called) ploidy of the simulated sequencing data. '
            '\nThe ground-truth CNs of near-haploid cells are called by the same caller from the pre-simulated data. ',
    'with_haploidy_assumed_gametes.expected_ploidy': 'Expected (called) ploidy of the simulated sequencing data. '
            '\nThe ground-truth CNs of near-haploid cells are one-valued vectors (i.e., CN=1 everywhere in the genome). ',
    # sample-dependent
    'average_seq_depth': 'average sequencing depth',
    'raw_total_sequences': 'total number of sequencing reads',
    'bases_mapped_cigar': 'total number of cigar-aware bases that are mapped',
    'reads_mapped': 'total number of sequencing reads that are mapped',
    # simulation dependent 
    'CNA_percent': 'percentage of copy-number alterations (CNAs) that are simulated',
    # result-dependent
    'observed_ploidy': 'observed ploidy inferred by the copy-number caller from the simulated sequencing data',
    'bed_1_cn0_genome_size': 'number of base pairs in the genome of the first sample with CN=0 (CN: copy number)',
    'bed_1_cn1_genome_size': 'number of base pairs in the genome of the first sample with CN=1 (CN: copy number)',
    'bed_1_cn2plus_genome_size': 'number of base pairs in the genome of the first sample with CN>1 (CN: copy number)',
    'bed_2_cn0_genome_size': 'number of base pairs in the genome of the second sample with CN=0 (CN: copy number)',
    'bed_2_cn1_genome_size': 'number of base pairs in the genome of the second sample with CN=1 (CN: copy number)',
    'bed_2_cn2plus_genome_size': 'number of base pairs in the genome of the second sample with CN>1 (CN: copy number)',
}

continuous_features = list(CONTINUOUS_FEATURES_NAME2DESC.keys())

CATEGORICAL_FEATURES_NAME2DESC = {
    # sample-dependent
    'donor': 'The human subject from which is the near-haploid cells are derived',
    'sampleType' : 'near-haploid cell type',
    'avgSpotLen' : 'average sequening read length' , # a few unique values
    # simulation-dependent
    'overall_ploidy' : 'either diploid or aneuploid', # diploid or aneuploid
    'cellLine' : 'cancer cell-line to simulate',
    'n_samples_mixed' : 'number of sequencing samples that are mixed for simulation (this is a technical detail)', # 0 or 1
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
    'accuracy': 'The accuracy of observed (called) versus expected (groundtruth) integer copy numbers',
    'PCC_intCN': 'The Pearson correlation coefficient of observed (called) versus expected (groundtruth) integer copy numbers',
    'PCC_nonintCN': 'The Pearson correlation coefficient of observed (called) non-integer copy number versus expected (groundtruth) integer copy number',
    'frac_cov_genome': 'The fraction of the human reference genome hg19 that is covered by the observed (called) copy-number profile. ',
    'breakpoint_f1score': 'The F1-score of detecting copy-number changes (breakpoints), representing precision-recall tradeoff' , 
    'breakpoint_precision': 'An observed (called) breakpoint is precise (true positive) if at least one expected (grountruth) breakpoint is within 200,000 base pairs', 
    'breakpoint_recall': 'An expected (groundtruth) breakpoint is recalled (true positive) if at least one observed (called) breakpoint is within 200,000 base pairs',    
}
the_perf_metrics = list(THE_PERF_METRIC_NAME2DESC.keys())
aneu_gametes_perf_metrics = [f'with_aneuploidy_aware_gametes.{m}' for m in the_perf_metrics]
hapl_gametes_perf_metrics = [f'with_haploidy_assumed_gametes.{m}' for m in the_perf_metrics]

gamete_type_to_perf_metrics = {
    'aneuploidy_aware_gametes': aneu_gametes_perf_metrics,
    'haploidy_assumed_gametes': hapl_gametes_perf_metrics,
}

# pivot: 'Caller'

def plot_main():
    fig1 = plt.figure(figsize=(1*7, 1*5), constrained_layout=True)
    callers = []
    perf_names = []
    perf_vals = []
    for rowidx, perf_metric in enumerate(the_perf_metrics):
        callers.extend(list(the_df['Caller']))
        perf_names.extend([perf_metric] * len(the_df))
        perf_vals.extend(list(the_df[('with_aneuploidy_aware_gametes.'+perf_metric)]))
    df = pd.DataFrame({'Copy-number callers': callers, 'Performance metrics': perf_names, 'Performances': perf_vals})
    plot_ret = sns.boxplot(data=df, x='Performance metrics', y='Performances', hue='Copy-number callers')
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
    plt.savefig(args.output + '_main.png', dpi=600)
    plt.close()

# https://www.doubao.com/chat/38421257704964354
def plot_grid_main():
    # --------------------------
    # 1. Data preparation (original logic)
    # --------------------------
    callers = []
    perf_names = []
    perf_vals = []
    for rowidx, perf_metric in enumerate(the_perf_metrics):
        callers.extend(list(the_df['Caller']))
        perf_names.extend([perf_metric] * len(the_df))
        perf_vals.extend(list(the_df[('with_aneuploidy_aware_gametes.'+perf_metric)]))
    
    # Tidy dataframe
    df = pd.DataFrame({
        'Caller': callers,
        'Metric': perf_names,
        'Performance': perf_vals
    })

    # --------------------------
    # 2. Create grid: Rows=Metrics, Cols=Callers
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

    # Plot boxplot in every grid cell
    g.map(
        sns.boxplot,
        'Caller',
        'Performance',
        color='#457B9D',
        linewidth=1.0,
        flierprops=dict(markersize=1.5, alpha=0.4)
    )

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
        ax.set_ylabel(metric, fontsize=10, labelpad=8)  # REAL metric label

    # --------------------------
    # LAST ROW ONLY: Show REAL CALLER names on X-axis
    # --------------------------
    unique_callers = df['Caller'].unique()
    for ax, caller in zip(g.axes[-1, :], unique_callers):
        ax.set_xlabel(caller, fontsize=10, labelpad=8)  # REAL caller label
        # Rotate long caller names if needed
        ax.tick_params(axis='x', labelrotation=30)

    # Hide all tick marks (clean look)
    for ax in g.axes.flat:
        ax.tick_params(bottom=False, left=False)

    # --------------------------
    # 4. Save final figure
    # --------------------------
    plt.tight_layout()
    plt.savefig(args.output + '_final_grid.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(args.output + '_final_grid.png', dpi=600, bbox_inches='tight')
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
    for rowidx, perf_metric in enumerate(the_perf_metrics):
        callers.extend(list(the_df['Caller']))
        perf_names.extend([perf_metric] * len(the_df))
        perf_vals.extend(list(the_df[('with_aneuploidy_aware_gametes.'+perf_metric)]))

    # Tidy dataframe
    df = pd.DataFrame({
        'Caller': callers,
        'Metric': perf_names,
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
    # 3. Plot one boxplot per metric (each row = one metric)
    # --------------------------
    for ax, metric in zip(axes, the_perf_metrics):
        # Subset data for THIS metric only
        sub_df = df[df['Metric'] == metric]
        
        # Plot boxplot: X = Callers, Y = Performance
        sns.boxplot(
            data=sub_df,
            x='Caller',
            y='Performance',
            ax=ax,
            color='steelblue',
            linewidth=1.2,
            flierprops=dict(markersize=2, alpha=0.5)
        )
        
        # --------------------------
        # Style each subplot
        # --------------------------
        ax.set_title(metric, fontsize=12, weight='bold')  # TITLE = metric name
        ax.set_xlabel('')                                 # No repeated x-label
        ax.tick_params(axis='x', labelsize=10, rotation=30)
        ax.tick_params(axis='y', labelsize=10)
        ax.grid(axis='y', alpha=0.3)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=15, ha='right')

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
            for gamete_type in ['with_haploidy_assumed_gametes', 'with_aneuploidy_aware_gametes']:
                gamete_perf_metric = gamete_type+'.'+perf_metric
                # x: feature; y: performance metric
                plot_df = caller_df[[feature]].copy()
                plot_df[perf_metric] = caller_df[gamete_perf_metric]
                # maybe_aneuploid
                # always_haploid
                # maybe_aneuploid: let the cells CNs be called from the original pre-simulated data; always_haploid: let the cell CNs be set to ones'
                plot_df['gamete_type'] = (
                        'maybe_aneuploid: Let CNs be called from the original pre-simulated data' 
                        if gamete_type == 'with_aneuploidy_aware_gametes' else 
                        'always_haploid: Let CNs be set to one-valued vectors')
                plot_dfs.append(plot_df)
            plot_df = pd.concat(plot_dfs).reset_index(drop=True)
            #print(plot_df)
            ax2 = fig1.add_subplot(gs[rowidx+1, colidx])
            if feature in categorical_features:
                #plot_ret = sns.boxplot    (data=plot_df, x=feature, y=perf_metric, hue='gamete_type', ax=ax2, whis=np.inf)
                plot_df.columns = [('HS1' if x == '345HS1' else x) for x in plot_df.columns] # prevent cluttering of words for the donor categorical variable
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
                ax2.set_title('Method\n' + caller, fontsize=20)

            if rowidx == len(the_perf_metrics) - 1:
                #ax2.set_xlabel(feature)
                ax2.set_xlabel('')
            else:
                ax2.set_xlabel('')
                #ax2.set_xticklabels('')
            if colidx == 0:
                ax2.set_ylabel(perf_metric, fontsize=16)
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
            title='Assumptions about near-haploid cell copy numbers (CNs)'
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
    factor = 'Factor ' + feature + ': ' + FEATURES_NAME2DESC.get(feature, 'TODO')
    if len(factor) > 200*5:
        fig1.supxlabel(factor, fontsize=16)
    else:
        fig1.supxlabel(factor, fontsize=20)
    fig1.supylabel('Performance metrics', fontsize=24)
    A2Z = [chr(i) for i in range(ord('A'), ord('Z') + 1)]
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
            pdf.savefig(fig1, bbox_inches='tight', dpi=100)
            plt.close(fig1)

