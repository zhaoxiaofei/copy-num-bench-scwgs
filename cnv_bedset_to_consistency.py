#!/usr/bin/env python
# Patched version of cnv_bedset_to_consistency.py
# (original: https://github.com/zhaoxiaofei/copy-num-bench-scwgs/blob/main/cnv_bedset_to_consistency.py)
#
# Fixes applied (each is marked with a [FIX n] comment at its site):
#  [FIX 1] (critical, breakpoint P/R) get_breakpoints() emitted a "breakpoint" between
#          every pair of consecutive segments without checking that the copy number
#          actually changes across the boundary.  Same-CN segments separated by a GAP
#          (routinely manufactured by the `bedtools intersect` truncation whenever any
#          input BED has an uncovered sub-interval) are not merged by merge_bed_with_cn(),
#          so both the observed and the expected side emitted a phantom breakpoint at the
#          identical gap midpoint; the pair matched at distance 0 and was counted as a TP.
#          This inflated precision AND recall -- including crediting breakpoints that the
#          benchmarked caller never called (a missed truth breakpoint hidden in a gap used
#          to come out as TP instead of FN).  Now only true CN transitions are emitted.
#  [FIX 2] Removed the dead evaluate_breakpoints_BUGGY(): it matched each observed
#          breakpoint to the FIRST expected one within the window (not the closest) and
#          allowed many-to-one matches (precision inflated to 1.0 by spamming).
#  [FIX 3] merge_bed_with_cn(): row[0]/row[1]/row[2] are integer keys on a label-indexed
#          Series (deprecated; FutureWarning on pandas 2.x, hard KeyError on pandas 3.x).
#          Now row.iloc[0]/iloc[1]/iloc[2].
#  [FIX 4] bedset_to_consistency(): the four `bedtools intersect` outputs are glued with
#          pd.concat(axis=1), which is positional; row-count asserts do not guarantee that
#          the (chr, start, end) columns are identical AND in the same order.  Now the
#          interval columns are asserted identical before concatenation, and the depth BED
#          row count is asserted before the zip() that used to truncate silently.
#  [FIX 5] The chained bedtools pipeline runs with `set -o pipefail`, so a failing
#          intermediate command is no longer masked by the last command's exit code.
#  [FIX 6] statsfile_to_avgdp(): returns None (rendered as JSON null) with a warning
#          instead of a -1 sentinel that silently poisoned average_seq_depth.
#  [FIX 7] --n-ref-bases CLI option (default hg19 = 3095677412); the covered-genome
#          fraction now uses the same reference size as the depth computation.
#  [FIX 8] obs2exp ploidy ratios use NaN-safe division (no ZeroDivisionError when the
#          expected ploidy is 0) and the ploidy assert no longer rejects NaN from empty
#          inputs.
#  [FIX 9] Copy-number values are coerced through a strict validator that rejects NaN and
#          non-integer CN values with a clear message (previously: obscure TypeError from
#          list indexing, or silent int() truncation of floats).
#  [FIX 10] perf.json is now strict JSON: non-finite floats (NaN/Inf) are written as null
#          and numpy scalars are converted through a json default= hook.
#  [FIX 11] --bp-window CLI option (default 200000, previously hardcoded).
#  [FIX 12] breakpoint TP/FP/FN/n_obs/n_exp are reported in perf.json next to
#          precision/recall, and evaluate_breakpoints() documents that both sides are
#          derived from the shared caller segmentation (CN-transition consistency, not a
#          raw-truth breakpoint benchmark).

import argparse, collections, functools, json, logging, math, os, subprocess, sys
import numpy as np
import pandas as pd
#from functools import reduce

def nandiv(a, b, nan_val=np.nan):
    if b != 0: return a/b
    else: return nan_val

# https://chat.z.ai/c/f51ab3db-ea55-4e03-a137-dd7c60278b65
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

def is_intlike(v, tol=1e-9):
    """True if v is a finite number arbitrarily close to an integer. [FIX 9]"""
    try:
        fv = float(v)
    except (TypeError, ValueError):
        return False
    return math.isfinite(fv) and abs(fv - round(fv)) <= tol

def as_int_cn(v, name='CN'):
    """Strictly validate and coerce a copy-number value to int. [FIX 9]"""
    if not is_intlike(v):
        raise ValueError(f'{name}={v!r} is NaN/non-finite or not an integer; invalid copy-number value. ')
    return int(round(float(v)))

def pre_sim_df_to_obsCN_to_genome_size(pre_sim_df, start_colname, end_colname):
    obsCN_to_genome_size = [0, 0, 0]
    for obsCN, genome_size in zip(pre_sim_df['obsCN'], pre_sim_df[end_colname] - pre_sim_df[start_colname]):
        obsCN = as_int_cn(obsCN, 'obsCN')                       # [FIX 9]
        assert obsCN >= 0, f"The obsCN={obsCN} is invalid!"
        obsCN = min((obsCN, 2))
        obsCN_to_genome_size[obsCN] += int(genome_size)         # [FIX 9] explicit cast for dtype safety
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
    above = sum([(as_int_cn(cn, 'cns2ploidy CN') * int(size)) for cn, size in zip(cns, sizes)])   # [FIX 9]
    below = sum([(      1 * int(size)) for cn, size in zip(cns, sizes)])
    return nandiv(float(above), float(below))

# cat /stor/zxf/cnv/refs/hg19.fa.dict | tail -n+2 | awk '{print $3}' | sed 's/LN://g'  | awk '{s += $1} END {print s}'
# 3095677412
# [FIX 7] hg19 default; override with --n-ref-bases.
N_REF_BASES = 3095677412

def statsfile_to_avgdp(sim_bam_stats, n_ref_bases=N_REF_BASES):
    #SN      bases mapped (cigar):   568 526 831       # more accurate
    ret = None                                                                                    # [FIX 6] was: -1 sentinel
    with open(sim_bam_stats) as file:
        for line in file:
            if line.startswith('SN\tbases mapped (cigar):\t'):
                ret = int(line.split('\t')[2]) / float(n_ref_bases)
                break
    if ret is None:
        logging.warning(f'No "SN<TAB>bases mapped (cigar):" line found in {sim_bam_stats}; average_seq_depth will be null. ')
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
        # [FIX 3] row[0]/row[1]/row[2] treated integer keys as positions on a label-indexed
        # Series (FutureWarning on pandas 2.x, hard KeyError on pandas 3.x).
        chr_, start, end, cn = row.iloc[0], row.iloc[1], row.iloc[2], row[cn_col]
        if (chr_ == current_chr) and (cn == current_cn) and (start <= current_end):
            current_end = max(current_end, end)
        else:
            if current_chr is not None:
                merged.append([current_chr, current_start, current_end, current_cn])
            current_chr, current_start, current_end, current_cn = chr_, start, end, cn
    if current_chr is not None:
        merged.append([current_chr, current_start, current_end, current_cn])
    return pd.DataFrame(merged, columns=["chr", "start", "end", "cn"])

# [FIX 2] The old evaluate_breakpoints_BUGGY() was removed: it matched each observed
# breakpoint to the FIRST expected one within the window (not the closest) and allowed
# many-to-one matches, which inflated precision.

# Example usage
#metrics = evaluate_breakpoints("inferCNV_breakpoints.bed", "ground_truth_dna_breakpoints.bed")
#print(metrics)

def get_breakpoints(merged_df):
    """
    Given a merged CN-segment bed dataframe (columns: chr, start, end, cn),
    sorted within each chromosome, return the set of breakpoints.

    A breakpoint is defined as the midpoint between the end of one region
    and the start of the immediately following region on the same
    chromosome:   bp = (prev_region.end + next_region.start) / 2

    The very first region's start and the very last region's end are NOT
    breakpoints (they are chromosome edges, not CN transitions).

    [FIX 1] A breakpoint is emitted ONLY where the copy number actually changes
    between the two adjacent regions.  Previously one was emitted between every
    pair of consecutive segments: same-CN segments separated by a gap (e.g. one
    created by `bedtools intersect` truncation where another BED has an uncovered
    sub-interval) are not merged by merge_bed_with_cn(), so both the observed and
    the expected side emitted a phantom breakpoint at the identical gap midpoint,
    which then matched at distance 0 and was counted as a TP -- inflating both
    precision and recall (and even crediting breakpoints the caller never called).

    Returns: dict {chr: sorted 1D numpy array of breakpoint positions}
    """
    bp_dict = {}
    for chr_name, group in merged_df.groupby("chr"):
        group = group.sort_values("start")
        if len(group) < 2:
            bp_dict[chr_name] = np.array([])
            continue
        ends    = group["end"].values[:-1]    # end of region i
        starts  = group["start"].values[1:]   # start of region i+1
        cns     = group["cn"].values
        changed = cns[:-1] != cns[1:]         # [FIX 1] keep only true CN transitions
        bp_dict[chr_name] = np.sort(((ends + starts) / 2.0)[changed])
    return bp_dict


def evaluate_breakpoints(obs_df, obs_col, exp_df, exp_col, window_size=200_000):
    """
    Match observed to expected breakpoints within window_size, enforcing one-to-one
    matching (greedy, closest pair first).

    [FIX 12] Semantics note: both the observed and the expected breakpoints are derived
    from the SAME (post-sim caller) segmentation restricted to the genome covered by all
    input BEDs; the expected side simply marks the boundaries at which the projected
    expected CN changes.  The resulting precision/recall therefore measure CN-transition
    CONSISTENCY on the shared segmentation (with the approximate truth projected onto it),
    not a classic breakpoint benchmark against the raw truth interval segmentation.
    """
    obs_df_1 = merge_bed_with_cn(obs_df, obs_col)
    exp_df_1 = merge_bed_with_cn(exp_df, exp_col)

    obs_bp_dict = get_breakpoints(obs_df_1)
    exp_bp_dict = get_breakpoints(exp_df_1)

    n_obs = sum(len(v) for v in obs_bp_dict.values())
    n_exp = sum(len(v) for v in exp_bp_dict.values())

    TP = 0

    for chr_name, obs_bps in obs_bp_dict.items():
        exp_bps = exp_bp_dict.get(chr_name, np.array([]))
        if len(exp_bps) == 0 or len(obs_bps) == 0:
            continue

        # Build all candidate (obs_idx, exp_idx, distance) pairs within window,
        # then greedily match closest pairs first, enforcing ONE-TO-ONE matching
        # on both sides. This prevents many observed breakpoints from all being
        # credited against a single expected breakpoint (which would let
        # precision be inflated for free by spamming observed breakpoints
        # around one true positive), and likewise prevents one observed
        # breakpoint from being double-counted against multiple expected ones.
        candidates = []
        for oi, obs_bp in enumerate(obs_bps):
            distances = np.abs(exp_bps - obs_bp)
            within = np.where(distances <= window_size)[0]
            for ei in within:
                candidates.append((distances[ei], oi, ei))

        candidates.sort(key=lambda x: x[0])  # closest pairs first

        matched_obs_idx = set()
        matched_exp_idx = set()
        for dist, oi, ei in candidates:
            if oi in matched_obs_idx or ei in matched_exp_idx:
                continue
            matched_obs_idx.add(oi)
            matched_exp_idx.add(ei)

        TP += len(matched_obs_idx)

    FP = n_obs - TP
    FN = n_exp - TP

    precision = TP / n_obs if n_obs else 0
    recall = TP / n_exp if n_exp else 0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

    return {
        "TP": TP,
        "FP": FP,
        "FN": FN,
        "n_obs": n_obs,    # [FIX 12]
        "n_exp": n_exp,    # [FIX 12]
        "precision": precision,
        "recall": recall,
        "f1score": f1
    }

# /mnt/e/cnv/data/S04.data2to4/SRR926960_12_sort_markdup_mut_ginkgo_intcns.bed
# /mnt/e/cnv/data/S04_FPN_202_COLO-829.data3/SRR926953_SRR926953_12_COLO-829.approx_truth.bed
# /mnt/e/cnv/data/S04_FPN_202_COLO-829.data3to4/SRR926953_SRR926958_12_COLO-829_SRR926953_SRR926958_12_COLO-829_ginkgo_intcns.bed
def bedset_to_consistency(pre_sim_call_bed_1_fname, pre_sim_call_bed_2_fname, approx_truth_bed_fname, post_sim_call_bed_int_fname, post_sim_call_bed_dep_fname, sim_bam_stats,
        chrom_colname, start_colname, end_colname, n_ref_bases=N_REF_BASES, bp_window=200_000):
    if not sim_bam_stats:
        sim_bam_stats = approx_truth_bed_fname.replace('_simtruth.bed', '.bam.stats')
    if os.path.exists(sim_bam_stats):                                        # [FIX 6] missing file -> null, not -1
        avgDP = statsfile_to_avgdp(sim_bam_stats, n_ref_bases=n_ref_bases)   # [FIX 7]
    else:
        logging.warning(f'sim_bam_stats file not found: {sim_bam_stats}; average_seq_depth will be null. ')
        avgDP = None

    post_sim_call_bed_multiinter = change_file_ext(post_sim_call_bed_int_fname, 'multiinter.bed',             'bed')
    pre_sim_call_bed_1_inter     = change_file_ext(post_sim_call_bed_int_fname, 'intersect_pre_sim_1.bed',    'bed')
    pre_sim_call_bed_2_inter     = change_file_ext(post_sim_call_bed_int_fname, 'intersect_pre_sim_2.bed',    'bed')
    approx_truth_bed_inter       = change_file_ext(post_sim_call_bed_int_fname, 'intersect_approx_truth.bed', 'bed')
    post_sim_call_bed_inter      = change_file_ext(post_sim_call_bed_int_fname, 'intersect_post_sim_CN.bed',  'bed')
    post_sim_call_bed_by_DP_inter= change_file_ext(post_sim_call_bed_int_fname, 'intersect_post_sim_DP.bed',  'bed')
    post_sim_call_bed_multiinter_cmd = (
        # [FIX 5] `set -o pipefail` so that a failure of an intermediate bedtools command
        # is not masked by the exit code of the last command in the pipeline.
        F'''set -o pipefail; bedtools intersect -a {pre_sim_call_bed_1_fname} -b {pre_sim_call_bed_2_fname} '''
        F'''| bedtools intersect -a {approx_truth_bed_fname} -b - | bedtools intersect -a {post_sim_call_bed_int_fname} -b - > {post_sim_call_bed_multiinter}'''
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

    # [FIX 4] The four intersect outputs are glued position-wise by pd.concat(axis=1)
    # below, so their interval columns must be identical in BOTH values and row order;
    # equal row counts alone do not guarantee that, and a violation would silently
    # misassign every CN column to the wrong interval.
    interval_cols = [chrom_colname, start_colname, end_colname]
    for df_name, other_df in (
            ('pre_sim_call_bed_1', pre_sim_df_1),
            ('pre_sim_call_bed_2', pre_sim_df_2),
            ('approx_truth_bed'  , approx_truth_df)):
        assert other_df[interval_cols].reset_index(drop=True).equals(post_sim_call_df[interval_cols].reset_index(drop=True)), (
            F'Interval columns of {df_name} differ from {post_sim_call_bed_int_fname} (values or row order): '
            F'the four "bedtools intersect" outputs are not 1:1, and pd.concat(axis=1) would silently misassign the CN columns. ')

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

    bp_cols   = [chrom_colname, start_colname, end_colname]
    obs_bp_df = post_sim_call_df[bp_cols + ['obsCN']].copy()
    exp_bp_df = post_sim_call_df[bp_cols].copy()
    exp_bp_df['expCN']        = merged_df['expCN'].values
    exp_bp_df['approx_expCN'] = merged_df['approx_expCN'].values

    exact_breakpoint_metrics  = evaluate_breakpoints(obs_bp_df, 'obsCN', exp_bp_df, 'expCN',        window_size=bp_window)   # [FIX 11]
    approx_breakpoint_metrics = evaluate_breakpoints(obs_bp_df, 'obsCN', exp_bp_df, 'approx_expCN', window_size=bp_window)   # [FIX 11]

    expCN_to_genome_size_accuracy_w_lin_corr_coef = {}
    for expCN_colname in ['expCN', 'approx_expCN']:
        confusion_matrix_int = [([0]*(8+1)) for _ in range(8+1)]
        obs_exp_to_cn = collections.defaultdict()
        for obsCN, expCN, genomesize in zip(post_sim_call_df['obsCN'], merged_df[expCN_colname], post_sim_call_df[end_colname] - post_sim_call_df[start_colname]):
            obsCN = as_int_cn(obsCN, 'obsCN')                   # [FIX 9]
            expCN = as_int_cn(expCN, 'expCN')                   # [FIX 9]
            assert obsCN >= 0, f"The obsCN={obsCN} is invalid!"
            assert expCN >= 0, f"The expCN={expCN} is invalid!"
            obsCN, expCN = min((8, obsCN)), min((8, expCN))
            confusion_matrix_int[obsCN][expCN] += int(genomesize)
        genome_size, accuracy = cmat_to_genome_size_and_accuracy(confusion_matrix_int)
        if post_sim_call_bed_dep_fname:
            post_sim_call_df_by_DP = pd.read_csv(post_sim_call_bed_by_DP_inter, sep='\t', header=0)
            # [FIX 4] the zip() below used to truncate silently if the DP file had fewer rows.
            assert len(post_sim_call_df_by_DP) == len(post_sim_call_df), (
                F'{len(post_sim_call_df_by_DP)} == {len(post_sim_call_df)} failed: '
                F'the depth-based BED is not 1:1 with the CN-based BED. ')
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

    assert pd.isna(expCN_ploidy) or expCN_ploidy >= 0, f"The expCN_ploidy={expCN_ploidy} is invalid!"   # [FIX 8] NaN-safe
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
        'with_aneuploidy_aware_gametes.obs2exp_ploidy_ratio' : nandiv(float(obsCN_ploidy), float(expCN_ploidy)),   # [FIX 8] NaN-safe division
        'with_aneuploidy_aware_gametes.genome_size'     : expCN_to_genome_size_accuracy_w_lin_corr_coef['expCN'][0],

        'with_aneuploidy_aware_gametes.accuracy'        : expCN_to_genome_size_accuracy_w_lin_corr_coef['expCN'][1],
        'with_aneuploidy_aware_gametes.PCC_intCN'       : expCN_to_genome_size_accuracy_w_lin_corr_coef['expCN'][2],
        'with_aneuploidy_aware_gametes.PCC_nonintCN'    : expCN_to_genome_size_accuracy_w_lin_corr_coef['expCN'][3],
        'with_aneuploidy_aware_gametes.frac_cov_genome' : expCN_to_genome_size_accuracy_w_lin_corr_coef['expCN'][0] / float(n_ref_bases),   # [FIX 7]

        'with_aneuploidy_aware_gametes.breakpoint_n_obs'     : exact_breakpoint_metrics['n_obs'],     # [FIX 12]
        'with_aneuploidy_aware_gametes.breakpoint_n_exp'     : exact_breakpoint_metrics['n_exp'],     # [FIX 12]
        'with_aneuploidy_aware_gametes.breakpoint_TP'        : exact_breakpoint_metrics['TP'],        # [FIX 12]
        'with_aneuploidy_aware_gametes.breakpoint_FP'        : exact_breakpoint_metrics['FP'],        # [FIX 12]
        'with_aneuploidy_aware_gametes.breakpoint_FN'        : exact_breakpoint_metrics['FN'],        # [FIX 12]
        'with_aneuploidy_aware_gametes.breakpoint_precision' : exact_breakpoint_metrics['precision'],
        'with_aneuploidy_aware_gametes.breakpoint_recall'    : exact_breakpoint_metrics['recall'],
        'with_aneuploidy_aware_gametes.breakpoint_f1score'   : exact_breakpoint_metrics['f1score'],

        'with_haploidy_assumed_gametes.expected_ploidy' : approx_expCN_ploidy,
        'with_haploidy_assumed_gametes.obs2exp_ploidy_ratio' : nandiv(float(obsCN_ploidy), float(approx_expCN_ploidy)),   # [FIX 8] NaN-safe division
        'with_haploidy_assumed_gametes.genome_size'     : expCN_to_genome_size_accuracy_w_lin_corr_coef['approx_expCN'][0],

        'with_haploidy_assumed_gametes.accuracy'        : expCN_to_genome_size_accuracy_w_lin_corr_coef['approx_expCN'][1],
        'with_haploidy_assumed_gametes.PCC_intCN'       : expCN_to_genome_size_accuracy_w_lin_corr_coef['approx_expCN'][2],
        'with_haploidy_assumed_gametes.PCC_nonintCN'    : expCN_to_genome_size_accuracy_w_lin_corr_coef['approx_expCN'][3],
        'with_haploidy_assumed_gametes.frac_cov_genome' : expCN_to_genome_size_accuracy_w_lin_corr_coef['approx_expCN'][0] / float(n_ref_bases),   # [FIX 7]

        'with_haploidy_assumed_gametes.breakpoint_n_obs'     : approx_breakpoint_metrics['n_obs'],    # [FIX 12]
        'with_haploidy_assumed_gametes.breakpoint_n_exp'     : approx_breakpoint_metrics['n_exp'],    # [FIX 12]
        'with_haploidy_assumed_gametes.breakpoint_TP'        : approx_breakpoint_metrics['TP'],       # [FIX 12]
        'with_haploidy_assumed_gametes.breakpoint_FP'        : approx_breakpoint_metrics['FP'],       # [FIX 12]
        'with_haploidy_assumed_gametes.breakpoint_FN'        : approx_breakpoint_metrics['FN'],       # [FIX 12]
        'with_haploidy_assumed_gametes.breakpoint_precision' : approx_breakpoint_metrics['precision'],
        'with_haploidy_assumed_gametes.breakpoint_recall'    : approx_breakpoint_metrics['recall'],
        'with_haploidy_assumed_gametes.breakpoint_f1score'   : approx_breakpoint_metrics['f1score'],
    }

    # [FIX 10] strict-JSON hardening: non-finite floats (NaN/Inf) -> null; numpy scalars
    # converted through a default= hook.
    def _to_builtin(o):
        if hasattr(o, 'item'):
            try: return o.item()
            except Exception: pass
        return str(o)
    def _nonfinite_to_null(o):
        if isinstance(o, dict):  return {k: _nonfinite_to_null(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)): return [_nonfinite_to_null(v) for v in o]
        if isinstance(o, float) and not math.isfinite(o): return None
        return o
    post_sim_call_perf_json = change_file_ext(post_sim_call_bed_int_fname, 'perf.json', 'bed')
    with open(post_sim_call_perf_json, 'w') as file: json.dump(_nonfinite_to_null(data), file, indent=2, default=_to_builtin)
    return data

def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(pathname)s:%(lineno)d %(levelname)s - %(message)s')
    parser = argparse.ArgumentParser(description='Compute the copy-number (CN) consistency between BED files. ', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--pre-sim-call-bed-1', required=True, help='BED file corresponding to the BAM file used as input for simulating major CN. ')
    parser.add_argument('--pre-sim-call-bed-2', required=True, help='BED file corresponding to the BAM file used as input for simulating minor CN. ')
    parser.add_argument('--approx-truth-bed'  , required=True, help='BED file with simulated ground truth copy numbers assuming that the major and minor CNs are both one-valued vectors (i.e., from haploid cells). ')
    parser.add_argument('--post-sim-call-bed' , required=True, help='BED file with final integer (i.e., absolute ploidy) CN calling results. ')
    parser.add_argument('--post-sim-call-bed-by-DP', default='', help='BED file with final real-number (i.e., relative fragment depth) CN calling results. ')
    parser.add_argument('--sim-bam-stats'     , required=False, help='The output of running `samtools stats` on the BAM file that serves as the input of the benchmarked CN caller. ')
    parser.add_argument('--chrom', required=False, default='#chr_37', help='Chromosome column name in BED files. ')
    parser.add_argument('--start', required=False, default='start_37', help='Chromosome start position column name in BED files. ')
    parser.add_argument('--end'  , required=False, default='end_37', help='Chromosome end position column name in BED files. ')
    parser.add_argument('--bp-window'  , required=False, type=int, default=200_000     , help='Max distance (bp) for matching an observed breakpoint to an expected one. ')             # [FIX 11]
    parser.add_argument('--n-ref-bases', required=False, type=int, default=N_REF_BASES , help='Total number of reference bases (denominator of the average depth, and the genome size used by frac_cov_genome; default is hg19). ')   # [FIX 7]

    args = parser.parse_args()
    consistency = bedset_to_consistency(args.pre_sim_call_bed_1, args.pre_sim_call_bed_2, args.approx_truth_bed, args.post_sim_call_bed, args.post_sim_call_bed_by_DP, args.sim_bam_stats,
            args.chrom, args.start, args.end, n_ref_bases=args.n_ref_bases, bp_window=args.bp_window)

if __name__ == '__main__': main()

