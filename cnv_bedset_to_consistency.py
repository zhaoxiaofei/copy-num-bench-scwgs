#!/usr/bin/env python
import argparse, collections, functools, json, logging, os, subprocess, sys
import numpy as np
import pandas as pd
#from functools import reduce

def nandiv(a, b, nan_val=np.nan):
    if b != 0: return a/b
    else: return nan_val

def intersect_intervals(set1, set2):
    result = []
    i, j = 0, 0
    while i < len(set1) and j < len(set2):
        # Find the overlap
        start = max(set1[i][0], set2[j][0])
        end = min(set1[i][1], set2[j][1])

        if start <= end:  # Valid overlap
            result.append((start, end))

        # Move the pointer that points to the earlier ending interval
        if set1[i][1] < set2[j][1]:
            i += 1
        else:
            j += 1

    return result

def find_intersection(intervals):
    interval = intervals[0]
    for i in range(1, len(intervals)):
        interval = intersect_intervals(interval, intervals[i])
    return interval

def bedfile2dict(bed_filename):
    chrom2intervals = collections.defaultdict(list)
    with open(bed_filename) as file:
        for line in file:
            if line.startswith('#') or line.lower().startswith('track'): continue
            tokens = line.split()
            chrom, start, end = str(tokens[0]), int(tokens[1]), int(tokens[2])
            chrom2intervals[chrom].append((start, end))
    return {chrom: sorted(intervals) for chrom, intervals in sorted(chrom2intervals.items())}
    
def change_file_ext(file_path, new_extension, old_extension=''):
    base_name, old_ext = os.path.splitext(file_path)
    if old_extension: assert '.'+old_extension == old_ext, F'File extension check: {old_extension} == {old_ext} failed!'
    return base_name + "." + new_extension

def pre_sim_df_to_obsCN_to_genome_size(pre_sim_df, start_colname, end_colname):
    obsCN_to_genome_size = [0, 0, 0]
    for obsCN, genome_size in zip(pre_sim_df['obsCN'], pre_sim_df[end_colname] - pre_sim_df[start_colname]):
        obsCN = min((obsCN, 2))
        obsCN_to_genome_size[obsCN] += genome_size
    return obsCN_to_genome_size

def weighted_lin_corr_coef(X, Y, weights):
    # Example weighted arrays
    #X = np.array([1, 2, 3, 4, 5])  # Array X
    #Y = np.array([0.5, 1.5, 2.5, 3.5, 4.5])  # Array Y
    #weights = np.array([0.1, 0.1, 0.1, 0.1, 0.1])  # Array of weights

    # Compute weighted means
    mean_X = np.average(X, weights=weights)
    mean_Y = np.average(Y, weights=weights)

    # Compute weighted covariance
    cov_xy = nandiv(np.sum(weights * (X - mean_X) * (Y - mean_Y)), np.sum(weights))

    # Compute weighted standard deviations
    std_X = np.sqrt(nandiv(np.sum(weights * (X - mean_X)**2), np.sum(weights)))
    std_Y = np.sqrt(nandiv(np.sum(weights * (Y - mean_Y)**2), np.sum(weights)))

    # Compute weighted correlation coefficient
    weighted_corr = nandiv(cov_xy, (std_X * std_Y))
    return weighted_corr
    #print("Weighted Pearson correlation coefficient:", weighted_corr)

def cmat_to_genome_size_and_accuracy(confusion_matrix_int):
    cmat_size = len(confusion_matrix_int)
    for row in confusion_matrix_int: assert len(row) == cmat_size, F'{len(row)} == {cmat_size}'
    genome_size = sum([col for row in confusion_matrix_int for col in row])
    accuracy = nandiv(sum([confusion_matrix_int[i][i] for i in range(cmat_size)]), float(genome_size))
    return genome_size, accuracy

def cns2ploidy(cns, sizes):
    above = sum([(int(cn) * int(size)) for cn, size in zip(cns, sizes)])
    below = sum([(      1 * int(size)) for cn, size in zip(cns, sizes)])
    return nandiv(float(above), float(below))

# cat /stor/zxf/cnv/refs/hg19.fa.dict | tail -n+2 | awk '{print $3}' | sed 's/LN://g'  | awk '{s += $1} END {print s}'
# 3095677412
N_REF_BASES = 3095677412

def statsfile_to_avgdp(sim_bam_stats, n_ref_bases=N_REF_BASES):
    #SN      bases mapped (cigar):   568 526 831       # more accurate
    ret = -1
    with open(sim_bam_stats) as file:
        for line in file:
            if line.startswith('SN\tbases mapped (cigar):\t'):
                ret = int(line.split('\t')[2]) / float(n_ref_bases)
                break
    return ret

def merge_bed_with_cn(df, cn_col_name):
    df = df.copy()
    #df = pd.read_csv(bed_path, sep="\t", header=None)
    if df.shape[1] < 4:
        raise ValueError("BED file must have at least 4 columns (chr, start, end, cn).")
    df.columns = [f'{col}_{i}' if col in df.columns[:i] else col 
             for i, col in enumerate(df.columns)]
    cn_col = next((col for col in df.columns if col.strip().upper() == cn_col_name.upper()), None)
    if cn_col is None: raise ValueError(f"Column '{cn_col_name}' not found in BED file.")
    df = df.sort_values(by=list(df.columns))# .reset_index(drop=True)
    merged = []
    current_chr, current_start, current_end, current_cn = None, None, None, None
    for _, row in df.iterrows():
        chr_, start, end, cn = row[0], row[1], row[2], row[cn_col]
        if (chr_ == current_chr) and (cn == current_cn) and (start <= current_end):
            current_end = max(current_end, end)
        else:
            if current_chr is not None:
                merged.append([current_chr, current_start, current_end, current_cn])
            current_chr, current_start, current_end, current_cn = chr_, start, end, cn
    if current_chr is not None:
        merged.append([current_chr, current_start, current_end, current_cn])
    return pd.DataFrame(merged, columns=["chr", "start", "end", "cn"])

def evaluate_breakpoints(obs_df, obs_col, exp_df, exp_col, window_size=200000):
    obs_df_1 = merge_bed_with_cn(obs_df, obs_col) #pd.read_csv(obs_df_1_bed, sep="\t", header=None, names=["chr", "start", "end"])
    exp_df_1 = merge_bed_with_cn(exp_df, exp_col) #pd.read_csv(exp_df_1_bed, sep="\t", header=None, names=["chr", "start", "end"])

    exp_df_1_dict = {}
    for chr_group, group in exp_df_1.groupby("chr"):
        exp_df_1_dict[chr_group] = group[["start", "end"]].values

    TP = 0
    FP = 0
    FN = 0

    exp_df_1_matched = set()

    for _, obs_df_1_row in obs_df_1.iterrows():
        chr_obs_df_1 = obs_df_1_row["chr"]
        start_obs_df_1 = obs_df_1_row["start"]
        matched = False

        if chr_obs_df_1 in exp_df_1_dict:
            for exp_df_1_entry in exp_df_1_dict[chr_obs_df_1]:
                start_exp_df_1 = exp_df_1_entry[0]
                distance = abs(start_obs_df_1 - start_exp_df_1)
                if distance <= window_size:
                    TP += 1
                    matched = True
                    exp_df_1_matched.add((chr_obs_df_1, start_exp_df_1))
                    break  # Match to closest (optional)

        if not matched:
            FP += 1

    for chr_exp_df_1, group in exp_df_1.groupby("chr"):
        for _, row in group.iterrows():
            start_exp_df_1 = row["start"]
            if (chr_exp_df_1, start_exp_df_1) not in exp_df_1_matched:
                FN += 1

    precision = TP / (TP + FP) if (TP + FP) > 0 else 0
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

    return {
        "TP": TP,
        "FP": FP,
        "FN": FN,
        "precision": precision,
        "recall": recall,
        "f1score": f1
    }

# Example usage
#metrics = evaluate_breakpoints("inferCNV_breakpoints.bed", "ground_truth_dna_breakpoints.bed")
#print(metrics)

# /mnt/e/cnv/data/S04.data2to4/SRR926960_12_sort_markdup_mut_ginkgo_intcns.bed
# /mnt/e/cnv/data/S04_FPN_202_COLO-829.data3/SRR926953_SRR926953_12_COLO-829.approx_truth.bed
# /mnt/e/cnv/data/S04_FPN_202_COLO-829.data3to4/SRR926953_SRR926958_12_COLO-829_SRR926953_SRR926958_12_COLO-829_ginkgo_intcns.bed
def bedset_to_consistency(pre_sim_call_bed_1_fname, pre_sim_call_bed_2_fname, approx_truth_bed_fname, post_sim_call_bed_int_fname, post_sim_call_bed_dep_fname, sim_bam_stats,
        chrom_colname, start_colname, end_colname):
    if not sim_bam_stats:
        sim_bam_stats = approx_truth_bed_fname.replace('_simtruth.bed', '.bam.stats')
    avgDP = statsfile_to_avgdp(sim_bam_stats)
    #pre_sim_call_bed_1_dict    = bedfile2dict(pre_sim_call_bed_1_fname)
    #pre_sim_call_bed_2_dict    = bedfile2dict(pre_sim_call_bed_2_fname)
    #approx_truth_bed_dict      = bedfile2dict(approx_truth_bed_fname)
    #post_sim_call_bed_int_dict = bedfile2dict(post_sim_call_int_fname)
    #post_sim_call_bed_dep_dict = bedfile2dict(post_sim_call_dep_fname)
    #bed_dicts = [pre_sim_call_bed_1_dict, pre_sim_call_bed_2_dict, approx_truth_bed_dict, post_sim_call_bed_int_dict, post_sim_call_bed_dep_dict]
    #chroms = set.intersection(*[set(bed_dict.keys()) for bed_dict in bed_dicts])
    #for chrom in chroms: ...

    post_sim_call_bed_multiinter = change_file_ext(post_sim_call_bed_int_fname, 'multiinter.bed',             'bed')
    pre_sim_call_bed_1_inter     = change_file_ext(post_sim_call_bed_int_fname, 'intersect_pre_sim_1.bed',    'bed')
    pre_sim_call_bed_2_inter     = change_file_ext(post_sim_call_bed_int_fname, 'intersect_pre_sim_2.bed',    'bed')
    approx_truth_bed_inter       = change_file_ext(post_sim_call_bed_int_fname, 'intersect_approx_truth.bed', 'bed')
    post_sim_call_bed_inter      = change_file_ext(post_sim_call_bed_int_fname, 'intersect_post_sim_CN.bed',  'bed')
    post_sim_call_bed_by_DP_inter= change_file_ext(post_sim_call_bed_int_fname, 'intersect_post_sim_DP.bed',  'bed')
    post_sim_call_bed_multiinter_cmd = (
        F'''bedtools intersect -a {pre_sim_call_bed_1_fname} -b {pre_sim_call_bed_2_fname} '''
        F'''| bedtools intersect -a {approx_truth_bed_fname} -b - | bedtools intersect -a {post_sim_call_bed_int_fname} -b - > {post_sim_call_bed_multiinter}'''
        #F''' bedtools multiinter -i <(cat {pre_sim_call_bed_1_fname} | grep -v ^chr_37) <(cat {pre_sim_call_bed_2_fname} | grep -v ^chr_37) '''
        #F''' <(cat {approx_truth_bed_fname} | grep -v ^chr_37) <(cat {post_sim_call_bed_int_fname} | grep -v ^chr_37) | awk '$4==4' > {post_sim_call_bed_multiinter} '''
    )
    cmd1 = F''' bedtools intersect -header -a {pre_sim_call_bed_1_fname     } -b {post_sim_call_bed_multiinter} > {pre_sim_call_bed_1_inter} '''
    cmd2 = F''' bedtools intersect -header -a {pre_sim_call_bed_2_fname     } -b {post_sim_call_bed_multiinter} > {pre_sim_call_bed_2_inter} '''
    cmd3 = F''' bedtools intersect -header -a {approx_truth_bed_fname       } -b {post_sim_call_bed_multiinter} > {approx_truth_bed_inter  } '''
    cmd4 = F''' bedtools intersect -header -a {post_sim_call_bed_int_fname  } -b {post_sim_call_bed_multiinter} > {post_sim_call_bed_inter } '''    
    cmd5 = F''' bedtools intersect -header -a {post_sim_call_bed_dep_fname  } -b {post_sim_call_bed_multiinter} > {post_sim_call_bed_by_DP_inter} '''
    if not post_sim_call_bed_dep_fname: cmd5 = F'printf "Skip generating {post_sim_call_bed_by_DP_inter}\\n"'
    
    for cmd in [post_sim_call_bed_multiinter_cmd, cmd1, cmd2, cmd3, cmd4, cmd5]:
        logging.info('Executing: ' + cmd)
        subprocess.run(cmd, shell=True, check=True, executable='/usr/bin/bash')

    pre_sim_df_1_raw = pd.read_csv(pre_sim_call_bed_1_fname, sep='\t', header=0)
    pre_sim_df_2_raw = pd.read_csv(pre_sim_call_bed_2_fname, sep='\t', header=0)

    obsCN_to_genome_size_1 = pre_sim_df_to_obsCN_to_genome_size(pre_sim_df_1_raw, start_colname, end_colname)
    obsCN_to_genome_size_2 = pre_sim_df_to_obsCN_to_genome_size(pre_sim_df_2_raw, start_colname, end_colname)

    pre_sim_df_1     = pd.read_csv(pre_sim_call_bed_1_inter, sep='\t', header=0)
    pre_sim_df_2     = pd.read_csv(pre_sim_call_bed_2_inter, sep='\t', header=0)

    approx_truth_df  = pd.read_csv(approx_truth_bed_inter  , sep='\t', header=0)
    post_sim_call_df = pd.read_csv(post_sim_call_bed_inter , sep='\t', header=0)

    assert len(post_sim_call_df) == len(pre_sim_df_1),    F'{len(post_sim_call_df)} == {len(pre_sim_df_1)} failed!'
    assert len(post_sim_call_df) == len(pre_sim_df_2),    F'{len(post_sim_call_df)} == {len(pre_sim_df_2)} failed!'
    assert len(post_sim_call_df) == len(approx_truth_df), F'{len(post_sim_call_df)} == {len(approx_truth_df)} failed!'
    
    #merged_df = functools.reduce(lambda left, right: pd.merge(left, right, on=['#chr_37', start_colname, end_colname], how='inner'), [pre_sim_df_1, pre_sim_df_2, approx_truth_df, post_sim_call_df])
    merged_df = pd.concat([pre_sim_df_1, pre_sim_df_2, approx_truth_df, post_sim_call_df], axis=1)
    merged_df['expMajorCN'] = (pre_sim_df_1['obsCN'] * approx_truth_df['majorCN'])
    merged_df['expMinorCN'] = (pre_sim_df_2['obsCN'] * approx_truth_df['minorCN'])
    merged_df['expCN']      = merged_df['expMajorCN'] + merged_df['expMinorCN']
    merged_df['approx_expCN'] = approx_truth_df['majorCN'] + approx_truth_df['minorCN']
    #print(merged_df)
    expCN_ploidy =        cns2ploidy(merged_df['expCN'],        post_sim_call_df[end_colname] - post_sim_call_df[start_colname])
    approx_expCN_ploidy = cns2ploidy(merged_df['approx_expCN'], post_sim_call_df[end_colname] - post_sim_call_df[start_colname])
    obsCN_ploidy =        cns2ploidy(post_sim_call_df['obsCN'], post_sim_call_df[end_colname] - post_sim_call_df[start_colname])

    obsCN_bed1_ploidy =   cns2ploidy(pre_sim_df_1['obsCN'], pre_sim_df_1[end_colname] - pre_sim_df_1[start_colname])
    obsCN_bed2_ploidy =   cns2ploidy(pre_sim_df_2['obsCN'], pre_sim_df_2[end_colname] - pre_sim_df_2[start_colname])

    exact_breakpoint_metrics = evaluate_breakpoints(merged_df, 'obsCN', merged_df, 'expCN')
    approx_breakpoint_metrics = evaluate_breakpoints(merged_df, 'obsCN', merged_df, 'approx_expCN')
 
    expCN_to_genome_size_accuracy_w_lin_corr_coef = {}
    for expCN_colname in ['expCN', 'approx_expCN']:
        confusion_matrix_int = [([0]*(8+1)) for _ in range(8+1)]
        obs_exp_to_cn = collections.defaultdict()
        for obsCN, expCN, genomesize in zip(post_sim_call_df['obsCN'], merged_df[expCN_colname], post_sim_call_df[end_colname] - post_sim_call_df[start_colname]):
            obsCN, expCN = min((8, obsCN)), min((8, expCN))
            confusion_matrix_int[obsCN][expCN] += genomesize
        genome_size, accuracy = cmat_to_genome_size_and_accuracy(confusion_matrix_int)
        if post_sim_call_bed_dep_fname:
            post_sim_call_df_by_DP = pd.read_csv(post_sim_call_bed_by_DP_inter, sep='\t', header=0)
            xyw = [(obsRD, expCN, genomesize) for 
                    obsRD, expCN, genomesize in 
                    zip(post_sim_call_df_by_DP['obsDP'], merged_df[expCN_colname], post_sim_call_df[end_colname] - post_sim_call_df[start_colname])]
            x, y, w = zip(*xyw)
            w_lin_corr_coef_byDP = weighted_lin_corr_coef(x, y, w)
        else:
            w_lin_corr_coef_byDP = np.nan
        if True:
            xyw = [(obsRD, expCN, genomesize) for 
                    obsRD, expCN, genomesize in 
                    zip(post_sim_call_df['obsCN'], merged_df[expCN_colname], post_sim_call_df[end_colname] - post_sim_call_df[start_colname])]
            if xyw:
                x, y, w = zip(*xyw)
                w_lin_corr_coef_byCN = weighted_lin_corr_coef(x, y, w)
            else:
                w_lin_corr_coef_byCN = np.nan
        expCN_to_genome_size_accuracy_w_lin_corr_coef[expCN_colname] = (genome_size, accuracy, w_lin_corr_coef_byCN, w_lin_corr_coef_byDP)
        post_sim_call_perf_cmat = change_file_ext(post_sim_call_bed_int_fname, F'perf.{expCN_colname}.confusion_matrix', 'bed')
        pd.DataFrame(confusion_matrix_int, index=[F'obsCN={i}' for i in range(8+1)], columns=[F'expCN={i}' for i in range(8+1)]).to_csv(post_sim_call_perf_cmat)
    
    #print(F'expCN_to_genome_size_accuracy_w_lin_corr_coef={expCN_to_genome_size_accuracy_w_lin_corr_coef}')
    data = {
        'pre_sim_call_bed_1'        : pre_sim_call_bed_1_fname,
        'pre_sim_call_bed_2'        : pre_sim_call_bed_2_fname,
        'approx_truth_bed'          : approx_truth_bed_fname,
        'post_sim_call_CN_bed'      : post_sim_call_bed_int_fname,
        'post_sim_call_DP_bed'      : post_sim_call_bed_dep_fname,
        'bed_1_cn0_genome_size'     : obsCN_to_genome_size_1[0],
        'bed_1_cn1_genome_size'     : obsCN_to_genome_size_1[1],
        'bed_1_cn2plus_genome_size' : obsCN_to_genome_size_1[2],
        'bed_2_cn0_genome_size'     : obsCN_to_genome_size_2[0],
        'bed_2_cn1_genome_size'     : obsCN_to_genome_size_2[1],
        'bed_2_cn2plus_genome_size' : obsCN_to_genome_size_2[2],
        
        'average_seq_depth'         : avgDP,
        'observed_bed1_ploidy'                          : obsCN_bed1_ploidy,
        'observed_bed2_ploidy'                          : obsCN_bed2_ploidy,
        'observed_ploidy'                               : obsCN_ploidy,
        # gini?
        
        'with_aneuploidy_aware_gametes.expected_ploidy' : expCN_ploidy, # observed_ploidy
        'with_aneuploidy_aware_gametes.obs2exp_ploidy_ratio' : obsCN_ploidy / expCN_ploidy,
        'with_aneuploidy_aware_gametes.genome_size'     : expCN_to_genome_size_accuracy_w_lin_corr_coef['expCN'][0],
        
        'with_aneuploidy_aware_gametes.accuracy'        : expCN_to_genome_size_accuracy_w_lin_corr_coef['expCN'][1],
        'with_aneuploidy_aware_gametes.PCC_intCN'       : expCN_to_genome_size_accuracy_w_lin_corr_coef['expCN'][2],
        'with_aneuploidy_aware_gametes.PCC_nonintCN'    : expCN_to_genome_size_accuracy_w_lin_corr_coef['expCN'][3],
        'with_aneuploidy_aware_gametes.frac_cov_genome' : expCN_to_genome_size_accuracy_w_lin_corr_coef['expCN'][0] / float(N_REF_BASES),
        
        'with_aneuploidy_aware_gametes.breakpoint_precision' : exact_breakpoint_metrics['precision'],
        'with_aneuploidy_aware_gametes.breakpoint_recall'    : exact_breakpoint_metrics['recall'],
        'with_aneuploidy_aware_gametes.breakpoint_f1score'   : exact_breakpoint_metrics['f1score'],
 
        'with_haploidy_assumed_gametes.expected_ploidy' : approx_expCN_ploidy,
        'with_haploidy_assumed_gametes.obs2exp_ploidy_ratio' : obsCN_ploidy / approx_expCN_ploidy,
        'with_haploidy_assumed_gametes.genome_size'     : expCN_to_genome_size_accuracy_w_lin_corr_coef['approx_expCN'][0],
        
        'with_haploidy_assumed_gametes.accuracy'        : expCN_to_genome_size_accuracy_w_lin_corr_coef['approx_expCN'][1],
        'with_haploidy_assumed_gametes.PCC_intCN'       : expCN_to_genome_size_accuracy_w_lin_corr_coef['approx_expCN'][2],
        'with_haploidy_assumed_gametes.PCC_nonintCN'    : expCN_to_genome_size_accuracy_w_lin_corr_coef['approx_expCN'][3],
        'with_haploidy_assumed_gametes.frac_cov_genome' : expCN_to_genome_size_accuracy_w_lin_corr_coef['approx_expCN'][0] / float(N_REF_BASES),
    
        'with_haploidy_assumed_gametes.breakpoint_precision' : approx_breakpoint_metrics['precision'],
        'with_haploidy_assumed_gametes.breakpoint_recall'    : approx_breakpoint_metrics['recall'],
        'with_haploidy_assumed_gametes.breakpoint_f1score'   : approx_breakpoint_metrics['f1score'],
    }

    post_sim_call_perf_json = change_file_ext(post_sim_call_bed_int_fname, 'perf.json', 'bed')
    with open(post_sim_call_perf_json, 'w') as file: json.dump(data, file, indent=2)
    return data

def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(pathname)s:%(lineno)d %(levelname)s - %(message)s')
    parser = argparse.ArgumentParser(description='Compute the copy-number (CN) consistency between BED files. ', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--pre-sim-call-bed-1', required=True, help='BED file corresponding to the BAM file used as input for simulating major CN. ')
    parser.add_argument('--pre-sim-call-bed-2', required=True, help='BED file corresponding to the BAM file used as input for simulating minor CN. ')
    parser.add_argument('--approx-truth-bed'  , required=True, help='BED file with simulated ground truth copy numbers assuming that the major and minor CNs are both one-valued vectors (i.e., from haploid cells). ')
    parser.add_argument('--post-sim-call-bed' , required=True, help='BED file with final integer (i.e., absolute ploidy) CN calling results. ')
    parser.add_argument('--post-sim-call-bed-by-DP', default='', help='BED file with final real-number (i.e., relative fragment depth) CN calling results. ')
    parser.add_argument('--sim-bam-stats'     , required=False,help='The output of running `samtools stats` on the BAM file taht serves as the input of the benchmarked CN caller. ')
    parser.add_argument('--chrom', required=False, default='#chr_37', help='Chromosome column name in BED files. ')
    parser.add_argument('--start', required=False, default='start_37', help='Chromosome start position column name in BED files. ')
    parser.add_argument('--end'  , required=False, default='end_37', help='Chromosome end position column name in BED files. ')

    args = parser.parse_args()
    consistency = bedset_to_consistency(args.pre_sim_call_bed_1, args.pre_sim_call_bed_2, args.approx_truth_bed, args.post_sim_call_bed, args.post_sim_call_bed_by_DP, args.sim_bam_stats,
            args.chrom, args.start, args.end)

if __name__ == '__main__': main()

