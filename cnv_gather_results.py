import argparse
import collections
import json
import logging
import os
import statistics
import sys

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.DEBUG)

def main():
    PERF_KEYS = [
        "bed_1_cn0_genome_size",
        "bed_1_cn0_genome_size",
        "bed_1_cn1_genome_size",
        "bed_1_cn2plus_genome_size",
        "bed_2_cn0_genome_size",
        "bed_2_cn1_genome_size",
        "bed_2_cn2plus_genome_size",
        "observed_ploidy",
        
        "with_aneuploidy_aware_gametes.expected_ploidy",
        "with_aneuploidy_aware_gametes.genome_size",
        'with_aneuploidy_aware_gametes.obs2exp_ploidy_ratio',

        "with_aneuploidy_aware_gametes.accuracy",
        "with_aneuploidy_aware_gametes.PCC_intCN",
        "with_aneuploidy_aware_gametes.PCC_nonintCN",
        "with_aneuploidy_aware_gametes.frac_cov_genome",
        
        "with_aneuploidy_aware_gametes.breakpoint_precision",
        "with_aneuploidy_aware_gametes.breakpoint_recall",
        "with_aneuploidy_aware_gametes.breakpoint_f1score",

        "with_haploidy_assumed_gametes.expected_ploidy",
        "with_haploidy_assumed_gametes.genome_size",
        'with_haploidy_assumed_gametes.obs2exp_ploidy_ratio',

        "with_haploidy_assumed_gametes.accuracy",
        "with_haploidy_assumed_gametes.PCC_intCN",
        "with_haploidy_assumed_gametes.PCC_nonintCN",
        "with_haploidy_assumed_gametes.frac_cov_genome",
        
        "with_haploidy_assumed_gametes.breakpoint_precision",
        "with_haploidy_assumed_gametes.breakpoint_recall",
        "with_haploidy_assumed_gametes.breakpoint_f1score",
    ]
    nan_keys = [
            "with_aneuploidy_aware_gametes.PCC_intCN", "with_aneuploidy_aware_gametes.PCC_nonintCN",
            "with_haploidy_assumed_gametes.PCC_intCN", "with_aneuploidy_aware_gametes.PCC_nonintCN",
    ]
    parser = argparse.ArgumentParser(description='Compute summary statistics related to mean (avg, sd) and median (min, Q1, Q2, Q3, max). ')
    #parser.add_argument('--caller', type=str, required=True, help='Single-cell copy-number caller (can be set to hmmcopy, ginkgo, etc.). ')
    parser.add_argument('-i', '--infiles', nargs='+', required=True, help='Input files with the .perf.json extension. ')
    parser.add_argument('-o', '--outprefix', required=True, help='Prefix of the output files, ending with the ``.short.tsv` and ``long.tsv`` extensions. ')
    args = parser.parse_args()
    
    long_dfs = []
    listof_dicts = []
    for caller in ['hmmcopy', 'ginkgo', 'copynumber', 'secnv', 'sccnv', 'scyn', 'chisel', 'aneufinder', 'flcna']:
        filenames = [filename for filename in args.infiles if F'_{caller}_' in filename.split('/')[-1]]
        key2vals = collections.defaultdict(list)
        for filename in filenames:
            key2vals['Caller'].append(caller)
            with open(filename) as file:
                js = json.load(file)
                for key in PERF_KEYS: assert key in js, F'The key {key} is not found in the json at {filename}'
                for key in js:
                    key2vals[key].append((js[key]))
                infojson_completed = False
                bamstats_completed = False
                infojson_fname, bamstats_fname = 'NA', 'NA'
                if 'approx_truth_bed' in js:
                    # /stor/zxf/cnv/data/S04/3from2_2_S04_3_FPN_202_COLO-829.datdir/2_S04_3_FPN_202_COLO-829_3from2_SRR926953_SRR926953_simtruth.bed
                    infojson_fname = js['approx_truth_bed'].replace('_simtruth.bed', '_info.json')
                    if os.path.exists(infojson_fname):
                        with open(infojson_fname) as file:
                            infojson = json.load(file)
                            for jsonkey in ['donor', 'sampleType', 'avgSpotLen', 'cellLine', 'accession_1', 'accession_2']:
                                key2vals[jsonkey].append((infojson[jsonkey]))
                            random_thres = infojson['random_thres_str'].split('_')
                            key2vals['overall_ploidy'].append(random_thres[0])
                            key2vals['CNA_percent'].append(100.0 - 0.01*float(random_thres[1]))
                            infojson_completed = True
                    bamstats_fname = js['approx_truth_bed'].replace('_simtruth.bed', '.bam.stats')
                    if os.path.exists(bamstats_fname):
                        raw_total_sequences = reads_mapped = bases_mapped_cigar = -1
                        with open(bamstats_fname) as file:
                            for line in file:
                                if line.startswith('SN\traw total sequences:\t'): raw_total_sequences = int(line.split('\t')[2])
                                if line.startswith('SN\treads mapped:\t'): reads_mapped = int(line.split('\t')[2])
                                if line.startswith('SN\tbases mapped (cigar):\t'): bases_mapped_cigar = int(line.split('\t')[2])
                        assert raw_total_sequences != -1, F'The file {bamstats_fname} does not contain a line starting with ``SN\traw total sequences:\t``!'
                        assert reads_mapped != -1, F'The file {bamstats_fname} does not contain a line starting with ``SN\treads mapped:\t``!'
                        assert bases_mapped_cigar != -1, F'The file {bamstats_fname} does not contain a line starting with ``SN\tbases mapped (cigar):\t``!'
                        key2vals['raw_total_sequences'].append(raw_total_sequences)
                        key2vals['reads_mapped'].append(reads_mapped)
                        key2vals['bases_mapped_cigar'].append(bases_mapped_cigar)
                        bamstats_completed = True
                if not infojson_completed:
                    logging.warning(F'Was unable to follow the link: {filename} > {infojson_fname}')
                    for jsonkey in ['donor', 'sampleType', 'avgSpotLen', 'cellLine', 'accession_1', 'accession_2', 'overall_ploidy', 'CNA_percent']:
                        key2vals[jsonkey].append(np.nan)
                if not bamstats_completed:
                    logging.warning(F'Was unable to follow the link: {filename} > {bamstats_fname}')
                    for jsonkey in ['raw_total_sequences', 'reads_mapped', 'bases_mapped_cigar']:
                        key2vals[jsonkey].append(np.nan)
        maxlen = max([len(vs) for k,vs in (key2vals.items())])
        for k,vs in sorted(key2vals.items()):
            if len(vs) != maxlen: logging.info(F'{caller}: Skipping the key {k} because the key is only present in {len(vs)} rows out of {maxlen} rows')
        key2vals = {k : vs for k,vs in sorted(key2vals.items()) if (len(vs) == maxlen)}
        caller_specific_df = pd.DataFrame(key2vals)
        long_dfs.append(caller_specific_df)
        for key in PERF_KEYS:
            vals = key2vals[key]
            q1, q2, q3 = statistics.quantiles(vals, n=4)
            if key in nan_keys:
                avg = np.nanmean(vals)
                sd = np.nanstd(vals)
            else:
                avg = np.mean(vals)
                sd = np.std(vals)
            listof_dicts.append({'Caller': caller, 'metric': key, 'mean': avg,'std': sd, 'Q1': q1, 'Q2': q2, 'Q3': q3})
            #print(F'{key}:\tQ1={q1}\tQ2={q2}\tQ3={q3}')
    pd.concat(long_dfs).to_csv(args.outprefix + '.long.tsv', sep='\t', index=False, na_rep='NA')
    df = pd.DataFrame(listof_dicts)
    df = df.sort_values(['metric', 'Caller'])
    df.to_csv(args.outprefix + '.short.tsv', sep='\t', index=False, na_rep='NA')
    with open(args.outprefix + '.cmd.sh', 'w') as file: file.write('\\\n\t'.join(sys.argv))

if __name__ == '__main__': main()
