#!/usr/bin/env python

import argparse, io, copy, json, logging, os, random, sys
import pandas as pd

from types import SimpleNamespace

import common as cm
from common import write2file, find_replace_all

## Requirements to run the script generated to stdout
##   bwa, samtools, hg19
##   bcftools, eagleimp with eagleimp's database in the directory EAGLE_IMP_DB_DIR (if running chisel is required)

DOWNSAMPLE_METHODS = ['bases', 'flagstat']
cosmic_cell_lines = ['COLO-829', 'HCC1395', 'HeLa']    

def is_male  (sex): return sex and sex.lower() in ['m', 'male', 'man', 'guy', 'boy']
def is_female(sex): return sex and sex.lower() in ['f', 'female', 'woman', 'girl', 'w']

def GRCH37_CONST_PLOID_BEDTABLE(cn=1, sex=None):
    sexstr = ''
    if is_male(sex):
        sexstr = F'''
chrX,   chrX,   0,      155270560,   {cn},    0
chrY,   chrY,   0,      59373566,    {cn},    0
'''.strip()
    if is_female(sex):
        sexstr = F'''
chrX,   chrX,   0,      155270560,   {cn},    {cn}
chrY,   chrY,   0,      59373566,    0,       0
'''.strip()

    return  pd.read_csv(io.StringIO(F'''
,       chr_37, start_37,  end_37,  majorCN,    minorCN
chr1,   chr1,   0,      249250621,  {cn},       {cn}
chr2,   chr2,   0,      243199373,  {cn},       {cn}
chr3,   chr3,   0,      198022430,  {cn},       {cn}
chr4,   chr4,   0,      191154276,  {cn},       {cn}
chr5,   chr5,   0,      180915260,  {cn},       {cn}
chr6,   chr6,   0,      171115067,  {cn},       {cn}
chr7,   chr7,   0,      159138663,  {cn},       {cn}
chr8,   chr8,   0,      146364022,  {cn},       {cn}
chr9,   chr9,   0,      141213431,  {cn},       {cn}
chr10,  chr10,  0,      135534747,  {cn},       {cn}
chr11,  chr11,  0,      135006516,  {cn},       {cn}
chr12,  chr12,  0,      133851895,  {cn},       {cn}
chr13,  chr13,  0,      115169878,  {cn},       {cn}
chr14,  chr14,  0,      107349540,  {cn},       {cn}
chr15,  chr15,  0,      102531392,  {cn},       {cn}
chr16,  chr16,  0,      90354753,   {cn},       {cn}
chr17,  chr17,  0,      81195210,   {cn},       {cn}
chr18,  chr18,  0,      78077248,   {cn},       {cn}
chr19,  chr19,  0,      59128983,   {cn},       {cn}
chr20,  chr20,  0,      63025520,   {cn},       {cn}
chr21,  chr21,  0,      48129895,   {cn},       {cn}
chr22,  chr22,  0,      51304566,   {cn},       {cn}
'''.strip() + '\n' + sexstr), skipinitialspace=True, index_col=0, header=0)

#GRCH37_DIPLOID_BEDTABLE = GRCH37_CONST_PLOID_BEDTABLE(1, sex=None)
#GRCH37_DIPLOID_BEDTABLE = GRCH37_CONST_PLOID_BEDTABLE(1, sex=None)

def read_flagstat(flagstat_filename, flagstat2json):    
    if flagstat_filename in flagstat2json: jsonfile = flagstat2json[flagstat_filename]
    else:
        try:
            with open(flagstat_filename) as file:
                jsonfile = json.load(file)
                logging.info(F'The json file {flagstat_filename} is loaded')
        except FileNotFoundError as err:
            logging.error(err)
            jsonfile = {'QC-passed reads': {'primary mapped': 1}}
        except ValueError as err:
            logging.error(err)
            jsonfile = {'QC-passed reads': {'primary mapped': 1}}    
        flagstat2json[flagstat_filename] = jsonfile
    ret = int(jsonfile['QC-passed reads']['primary mapped'])
    return ret

def samtoolsfrac(r, frac):
    assert frac >= 0, F'{frac} >=0 failed!'
    if frac >= 1.0 - sys.float_info.epsilon:
        return ''
    else:
        rint = r.randint(0, 2**16-1)
        return F'-s {rint+frac}'
def flagstats2downfracs(flagstat1, flagstat2, flagstat2json, are_bases):
    n_prim_mapped_1 = (flagstat1 if are_bases else read_flagstat(flagstat1, flagstat2json))
    n_prim_mapped_2 = (flagstat2 if are_bases else read_flagstat(flagstat2, flagstat2json))
    n_prim_mapped_max = max((n_prim_mapped_1, n_prim_mapped_2))
    r = random.Random()
    r.seed(n_prim_mapped_1 + n_prim_mapped_2)
    return (samtoolsfrac(r, float(n_prim_mapped_2) / n_prim_mapped_max), 
            samtoolsfrac(r, float(n_prim_mapped_1) / n_prim_mapped_max))

# proc_haploCN(alltable, chr_tuple5, 'majorCN', major_cn_bams, major_bams, majorfrac, majorCN_order)
def proc_haploCN(cosmic_df, chr_tuple5, haploCN_colname, haplo_cn_bams, haplo_bams, haplo_frac, haploCN_order, max_haploCN):
    assert haploCN_colname in ['majorCN', 'minorCN'], F'The haploCN_colname {haploCN_colname} is not valid!'
    chr_37_colidx, chr_37_start_colidx, chr_37_end_colidx, chr_37_prev_end_colidx, chr_37_next_start_colidx = chr_tuple5
    spike_cmds = []
    haplo_gen_bams = []
    for haploCN in range(1, max_haploCN+1, 1):
        haplo_df = cosmic_df.loc[(cosmic_df[haploCN_colname]>=haploCN-0.5), :]
        regions = []
        for rowidx in range(haplo_df.shape[0]):
            norm_start = (haplo_df.iat[rowidx,chr_37_prev_end_colidx] + haplo_df.iat[rowidx,chr_37_start_colidx]) // 2
            norm_end = (haplo_df.iat[rowidx,chr_37_next_start_colidx] + haplo_df.iat[rowidx,chr_37_end_colidx]) // 2
            region = F'{haplo_df.iat[rowidx,chr_37_colidx]}:{int(norm_start)}-{int(norm_end)}'
            regions.append(region)
        if regions:
            region_str = ' '.join(regions)
            haplo_bam_idx = haploCN_order[haploCN-1]
            haplo_cmd = F'samtools view {haplo_frac} -M -bh -o {haplo_cn_bams[haplo_bam_idx]} {haplo_bams[haplo_bam_idx]} {region_str} && echo processed {haploCN_colname}={haploCN} '
            haplo_gen_bams.append(haplo_cn_bams[haplo_bam_idx])
        else:
            haplo_cmd = F'echo {haploCN_colname}={haploCN} is not applicable because the ground-truth does not contain such haplotype-specific ploidy. '
        spike_cmds.append(haplo_cmd)
    return spike_cmds, haplo_gen_bams

def simulate(infodict,
        cell_line_df, cell_line,
        lib_1, acc_1, lib_2, acc_2,
        inst2into3script, inst2into3tmpdir, inst3from2simbam, inst3from2simbed, inst3from2dedupb,
        flagstat2json, bases1, bases2, downsample_method, writing_mode,
        MAX_HAPLO_CN=cm.MAX_HAPLO_CN):
    
    subclonal_frac_idxs = []
    infodict = copy.deepcopy(infodict)
    infodict['accession'] = acc_1
    infodict['GT'] = 'A'
    inst2from1mutbam1, inst2from1cnbams1, inst2from1flagjs1 = find_replace_all([
    cm.t2from1mutbam , cm.t2from1cnbams , cm.t2from1flagjs], infodict)
    infodict['accession'] = acc_2
    infodict['GT'] = 'B'
    inst2from1mutbam2, inst2from1cnbams2, inst2from1flagjs2 = find_replace_all([
    cm.t2from1mutbam , cm.t2from1cnbams , cm.t2from1flagjs], infodict)
    
    random2 = random.Random()
    random2.seed(acc_1 + ' ' + acc_2)
    aneuploid_test_0to1rv = random2.random()
    subclonal_frac_0to1thres_rv = random2.random() * 0.5 # PMC8054914 : Figs 3 and S4
    majorCN_order = random2.sample(range(MAX_HAPLO_CN), MAX_HAPLO_CN)
    minorCN_order = random2.sample(range(MAX_HAPLO_CN), MAX_HAPLO_CN)
    random3 = random.Random()
    random3.seed(' '.join((infodict['donor'], infodict['sampleType'], infodict['avgSpotLen'], infodict['cellLine'])))
    if lib_1 == lib_2 or acc_1 == acc_2:
        assert acc_1 == acc_2 and lib_1 == lib_2, F'{acc_1} == {acc_2} and {lib_1} == {lib_2} failed!'
        assert MAX_HAPLO_CN >= 5, F'{MAX_HAPLO_CN} >= 5 failed!'
        logging.info(F'Setting MAX_MAJOR_CN={MAX_HAPLO_CN}-2 and MAX_MINOR_CN=2 because {acc_1}=={acc_2} and {lib_1}=={lib_2}. ')
        minorCN_order = majorCN_order[-2:]
        majorCN_order = majorCN_order[:-2]
        MAX_MAJOR_CN = MAX_HAPLO_CN - 2
        MAX_MINOR_CN = 2
    else:
        MAX_MAJOR_CN = MAX_MINOR_CN = MAX_HAPLO_CN
    is_diploid = (cell_line.split('-')[0].lower() == 'normal' or aneuploid_test_0to1rv < 0.5)
    aneuploid_frac_int_repr = int(subclonal_frac_0to1thres_rv*1000*10)
    random_thres_str = (F'diploid_{aneuploid_frac_int_repr:0>4}' if is_diploid else F'aneuploid_{aneuploid_frac_int_repr:0>4}')
    
    assert 'bases' == DOWNSAMPLE_METHODS[0]
    majorfrac, minorfrac = flagstats2downfracs(bases1, bases2, flagstat2json, (downsample_method in DOWNSAMPLE_METHODS[0]))
    #logging.info(F'{inprefix_1} majorfrac={majorfrac}, {inprefix_2} minorfrac={minorfrac}')
    major_bams, minor_bams = inst2from1cnbams1, inst2from1cnbams2 # [F'{inprefix_1}.cn_{i}.bam' for i in range(1,1+MAX_HAPLO_CN,1)]
    major_cn_bams = [F'{inst2into3tmpdir}/{acc_1}_{acc_2}_{cell_line}_cn{i}_{random_thres_str}_major.bam' for i in range(1,1+MAX_HAPLO_CN,1)]
    minor_cn_bams = [F'{inst2into3tmpdir}/{acc_1}_{acc_2}_{cell_line}_cn{i}_{random_thres_str}_minor.bam' for i in range(1,1+MAX_HAPLO_CN,1)]
    
    subtable = cell_line_df.copy(deep=True)
    subtable['chr_37'] = [F'chr{chrname}' for chrname in subtable['chr_37']]
    if      len(cell_line.split('-')) >= 2 and cell_line.split('-')[0].lower() == 'normal' and is_male  (cell_line.split('-')[1]):
        GRCH37_DIPLOID_BEDTABLE = GRCH37_CONST_PLOID_BEDTABLE(1, sex='M')
    elif len(cell_line.split('-')) >= 2 and cell_line.split('-')[0].lower() == 'normal' and is_female(cell_line.split('-')[1]):
        GRCH37_DIPLOID_BEDTABLE = GRCH37_CONST_PLOID_BEDTABLE(1, sex='F')
    elif 'chrY' in subtable['chr_37'] or 'Y' in subtable['chr_37']:
        GRCH37_DIPLOID_BEDTABLE = GRCH37_CONST_PLOID_BEDTABLE(1, sex='M')
    else:
        GRCH37_DIPLOID_BEDTABLE = GRCH37_CONST_PLOID_BEDTABLE(1, sex='F')

    with cm.myopen(inst2into3script, writing_mode) as shfile:        
        if is_diploid:
            spike_cmds = []
            majorCN = minorCN = 1
            spike_cmds.append(F'echo The bam file {major_cn_bams[majorCN-1]} is substituted by {major_bams[majorCN_order[majorCN-1]]}')
            spike_cmds.append(F'echo The bam file {minor_cn_bams[minorCN-1]} is substituted by {minor_bams[minorCN_order[minorCN-1]]}')
            major_gen_bams = [major_bams[majorCN_order[majorCN-1]]]
            minor_gen_bams = [minor_bams[minorCN_order[minorCN-1]]]
            for i in range(1, MAX_MAJOR_CN + 1, 1):
                spike_cmds.append(F'echo majorCN={i} is not applicable because the ground-truth is diploid-normal. ')
            for i in range(1, MAX_MINOR_CN + 1, 1):
                spike_cmds.append(F'echo minorCN={i} is not applicable because the ground-truth is diploid-normal. ')            
            retdf = GRCH37_DIPLOID_BEDTABLE.loc[~GRCH37_DIPLOID_BEDTABLE['chr_37'].isin(['chrX', 'chrY']),:]
            retdf.columns = retdf.columns.str.replace('chr_37','#chr_37')
            retdf.to_csv(inst3from2simbed, header=True, index=False, sep='\t')
        else:
            prev_ends = []
            next_starts = []            
            subtable['start_37'] = subtable['start_37'].astype('int32')
            subtable['end_37'] = subtable['end_37'].astype('int32')
            subtable['minorCN'] = subtable['minorCN'].astype('int32')
            subtable['totalCN'] = subtable['totalCN'].astype('int32')

            subtable_chrom_colidx = subtable.columns.get_loc('chr_37')
            subtable_start_colidx = subtable.columns.get_loc('start_37')
            subtable_end_colidx = subtable.columns.get_loc('end_37')
            subtable_minorCN_colidx = subtable.columns.get_loc('minorCN')
            subtable_totalCN_colidx = subtable.columns.get_loc('totalCN')
            prev_chrom = ''
            prev_rowidx = -1
            prev_region_size = -1
            rowidxs = list(range(subtable.shape[0]))
            for rowidx in (rowidxs + rowidxs[::-1][1:]):
                assert subtable.iat[rowidx, subtable_minorCN_colidx] * 2 <= subtable.iat[rowidx, subtable_totalCN_colidx], F'minorCN * 2 <= totalCN failed for the line {subtable.iat[rowidx,:]}'
                chrom = subtable.iat[rowidx, subtable_chrom_colidx]
                region_size = subtable.iat[rowidx, subtable_end_colidx] - subtable.iat[rowidx, subtable_start_colidx]
                subclonal_frac_0to1rv = random3.random()
                tie_break_rv = random3.randint(0, 2)
                if (chrom == prev_chrom and subclonal_frac_0to1rv < subclonal_frac_0to1thres_rv and tie_break_rv 
                        # and (region_size < prev_region_size or (region_size == prev_region_size and tie_break_rv)):
                        ):
                    subtable.iat[rowidx, subtable_minorCN_colidx] = subtable.iat[rowidx-1, subtable_minorCN_colidx]
                    subtable.iat[rowidx, subtable_totalCN_colidx] = subtable.iat[rowidx-1, subtable_totalCN_colidx]
                    subclonal_frac_idxs.append(rowidx if (rowidx > prev_rowidx) else -rowidx)
                prev_chrom = chrom
                prev_rowidx = rowidx
                prev_region_size = region_size
            bedtable = subtable[['chr_37', 'start_37', 'end_37', 'minorCN', 'totalCN']]
            bedtable = bedtable.astype({'start_37': int, 'end_37': int, 'minorCN': int, 'totalCN': int})
            minorCN_colidx = bedtable.columns.get_loc('minorCN')
            totalCN_colidx = bedtable.columns.get_loc('totalCN')
            
            # csvformat -T /mnt/d/oocyte/cosmic-v97/cell_lines_copy_number.csv | awk -F"\t" '{ if ($12-$11 <= 5) {s+=$14} else {t += $14} } END {print (t)/(s+t)}'
            # 0.00417797
            # csvformat -T /mnt/d/oocyte/cosmic-v97/cell_lines_copy_number.csv | awk -F"\t" '{ if ($11 <= 3) {s+=$14} else {t += $14} } END {print (t)/(s+t)}'
            # 0.00681828
            def rowidx2cn(i):
                totalCN = bedtable.iat[i,totalCN_colidx]
                minorCN = bedtable.iat[i,minorCN_colidx]
                majorCN = totalCN - minorCN
                return min((MAX_MAJOR_CN, majorCN)), min((MAX_MINOR_CN, minorCN))
                #return min((4, majorCN)) + min((2, minorCN))
                #return int(2*min(2,bedtable.iat[i,minorCN_colidx]) + min((4,bedtable.iat[i,totalCN_colidx] - 2*bedtable.iat[i,minorCN_colidx])))
            cns = [rowidx2cn(i) for i in range(bedtable.shape[0])]
            bedtable['CN'] = [(a+b) for a,b in cns]
            bedtable['majorCN'] = [a for a,b in cns]
            bedtable['minorCN'] = [b for a,b in cns]
            retdf = bedtable.loc[~bedtable['chr_37'].isin(['chrX', 'chrY']),:]
            retdf.columns = retdf.columns.str.replace('chr_37','#chr_37')
            retdf.to_csv(inst3from2simbed, header=True, index=False, sep='\t')
            
            chr_37_colidx = subtable.columns.get_loc('chr_37')
            chr_37_start_colidx = subtable.columns.get_loc('start_37')
            chr_37_end_colidx = subtable.columns.get_loc('end_37')
            for rowidx in range(subtable.shape[0]):
                if (rowidx == 0) or (subtable.iat[rowidx-1,chr_37_colidx] != subtable.iat[rowidx,chr_37_colidx]):
                    prev_chr = subtable.iat[rowidx,chr_37_colidx]
                    prev_start = -subtable.iat[rowidx,chr_37_start_colidx]
                    prev_end = -subtable.iat[rowidx,chr_37_start_colidx]
                else:
                    prev_chr = subtable.iat[rowidx-1,chr_37_colidx]
                    prev_start = subtable.iat[rowidx-1,chr_37_start_colidx]
                    prev_end = subtable.iat[rowidx-1,chr_37_end_colidx]
                if (rowidx == subtable.shape[0]-1) or (subtable.iat[rowidx+1,chr_37_colidx] != subtable.iat[rowidx,chr_37_colidx]):
                    next_chr = subtable.iat[rowidx,chr_37_colidx]
                    next_start = GRCH37_DIPLOID_BEDTABLE.loc[next_chr,'end_37']*2 - subtable.iat[rowidx,chr_37_end_colidx] # 1000*1000*1000 #infinity
                    next_end = next_start #1000*1000*1000
                else:
                    next_chr = subtable.iat[rowidx+1,chr_37_colidx]
                    next_start = subtable.iat[rowidx+1,chr_37_start_colidx]
                    next_end = subtable.iat[rowidx+1,chr_37_end_colidx]
                prev_ends.append(prev_end)
                next_starts.append(next_start)
            alltable = subtable
            alltable['prev_end'] = prev_ends
            alltable['next_start'] = next_starts
            assert alltable.shape[0] == subtable.shape[0]
            chr_37_colidx = alltable.columns.get_loc('chr_37')
            chr_37_start_colidx = alltable.columns.get_loc('start_37')
            chr_37_end_colidx = alltable.columns.get_loc('end_37')
            chr_37_prev_end_colidx = alltable.columns.get_loc('prev_end')
            chr_37_next_start_colidx = alltable.columns.get_loc('next_start') 
            alltable['majorCN'] = alltable['totalCN']-alltable['minorCN']
            
            chr_tuple5 = chr_37_colidx, chr_37_start_colidx, chr_37_end_colidx, chr_37_prev_end_colidx, chr_37_next_start_colidx
            major_spike_cmds, major_gen_bams = proc_haploCN(alltable, chr_tuple5, 'majorCN', major_cn_bams, major_bams, majorfrac, majorCN_order, MAX_MAJOR_CN)
            minor_spike_cmds, minor_gen_bams = proc_haploCN(alltable, chr_tuple5, 'minorCN', minor_cn_bams, minor_bams, minorfrac, minorCN_order, MAX_MINOR_CN)
            spike_cmds = major_spike_cmds + minor_spike_cmds
        # enf of if (is_diploid) else ...
        
        for cmd in spike_cmds: write2file(cmd, shfile, inst2into3script)
        indices_cmd = ' && '.join([F'samtools index {bam}' for bam in (minor_gen_bams + major_gen_bams)])
        write2file(indices_cmd, shfile, inst2into3script)
        #chrs='chr1,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr2,chr20,chr21,chr22,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chrX,chrY'
        write2file(F'samtools merge -f - ' + ' '.join(minor_gen_bams + major_gen_bams)
               +F' | samtools sort -n -o - - | samtools fixmate -u - - | samtools view -bhu -e \'flag.paired\' -F 12 '
               +F' | samtools sort -o {inst3from2simbam} - && samtools index {inst3from2simbam} && samtools stats {inst3from2simbam} > {inst3from2simbam}.stats', shfile, inst2into3script)
        #write2file(F'samtools view -F 0x400 -o {inst3from2dedupb} {inst3from2simbam} && samtools index {inst3from2dedupb}', shfile, inst2into3script)
        if is_diploid:
            write2file(F'echo approx_truth_bed={inst3from2simbed} for cell_line={cell_line} ploidy=diploid   cell_type=normal simulating {inst3from2simbam}', shfile, inst2into3script)
        else:
            write2file(F'echo approx_truth_bed={inst3from2simbed} for cell_line={cell_line} ploidy=aneuploid cell_type=tumor  simulating {inst3from2simbam}', shfile, inst2into3script)
    return random_thres_str, subclonal_frac_idxs

def main(args1=None):
    ret = []
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(pathname)s:%(lineno)d %(levelname)s - %(message)s')
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    datadir = os.path.abspath(os.path.sep.join([script_dir, '..', 'data']))
    data0to1dir, data1to2dir, data2to3dir, data2to4dir, data3to4dir, data4to5dir = cm.get_varnames(datadir)
    
    defaultSraRunTable = os.path.sep.join([script_dir, 'scDNAaccessions.tsv'])
    cosmic_cn_filename = os.path.sep.join([script_dir, 'cosmic-v97', 'cell_lines_copy_number.csv'])
    
    parser = argparse.ArgumentParser(description='Generate bash commands to in-silico mix the single-cell sequencing data with copy numbers simulated from the COSMIC cell line database. ',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--SraRunTable', type=str, default=defaultSraRunTable, help=(
            'SraRunTable in TSV format containing the columns '
            '#Run, AvgSpotLen, Library~Name, Sample~Name, sample-type, Oocyte_ID, Donor, and SRA~Study'))
    parser.add_argument('--cosmic',      type=str, default=cosmic_cn_filename, help='Copy-number profile TSV file downloaded from cancer.sanger.ac.uk/cosmic/download/cell-lines-project/v97')
    parser.add_argument('--cell-lines',  nargs='+',default=cosmic_cell_lines, help='Cell-lines to be used in --cosmic')
    parser.add_argument('--downsample-method', choices=DOWNSAMPLE_METHODS, default=DOWNSAMPLE_METHODS[0], help='Downsampling method')
    parser.add_argument('-w', '--writing-mode', type=str, default=cm.DEFAULT_WRITING_MODE,
        help='File open mode for writing commands to shell script, pass any of {cm.OVERWRITING_PREVENTION_MODES} to prevent overwriting existing scripts (or w to do not prevent such thing). ')
    args = (args1 if args1 else parser.parse_args())
    
    cosmic_df = pd.read_csv(args.cosmic)
    df0 = pd.read_csv(args.SraRunTable, sep='\t', header=0)
    df0['sample-type'] = cm.norm_sample_type(df0)
    grouped = df0.groupby(['AvgSpotLen', 'sample-type', 'Donor'])
    partitioned_dfs = {partkey: df1 for partkey, df1 in grouped}
    flagstat2json = {}
    for cell_line in args.cell_lines:
        cell_line_df = cosmic_df.loc[(cosmic_df['#sample_name'] == cell_line),:].copy(deep=True)
        for (avgSpotLen, sampleType, donor), df1 in sorted(partitioned_dfs.items()):
            infodict = {
                'data0to1dir': data0to1dir,
                'data1to2dir': data1to2dir,
                'data2to3dir': data2to3dir,
                'data2to4dir': data2to4dir,
                'donor'      : str(donor),
                'sampleType' : str(sampleType),
                'avgSpotLen' : str(avgSpotLen),
                'cellLine'   : str(cell_line),
            }
            inst2into3end, inst2into3logdir, inst2into3tmpdir, inst3from2datdir = find_replace_all([
            cm.t2into3end, cm.t2into3logdir, cm.t2into3tmpdir, cm.t3from2datdir], infodict)

            cm.makedirs((inst2into3logdir, inst2into3tmpdir, inst3from2datdir))
            with cm.myopen(inst2into3end, args.writing_mode) as file: write2file(F'echo {inst2into3end} is done', file, inst2into3end)

            for     rowidx_1, (acc_1, lib_1, sample_1, bases1) in enumerate(zip(df1['#Run'], df1['Library~Name'], df1['Sample~Name'], df1['Bases'])):
                for rowidx_2, (acc_2, lib_2, sample_2, bases2) in enumerate(zip(df1['#Run'], df1['Library~Name'], df1['Sample~Name'], df1['Bases'])):
                    if not cm.circular_dist_below(rowidx_1, rowidx_2, len(df1)): continue
                    infodict.update({
                        'accession_1': str(acc_1),
                        'accession_2': str(acc_2),
                    })
                    
                    inst2into3script, inst3from2simbam, inst3from2simbed, inst3from2dedupb, inst3from2infojs = find_replace_all([
                    cm.t2into3script, cm.t3from2simbam, cm.t3from2simbed, cm.t3from2dedupb, cm.t3from2infojs], infodict)
                    random_thres_str, subclonal_frac_idxs = simulate(infodict,
                            cell_line_df, cell_line, 
                            lib_1, acc_1, lib_2, acc_2,
                            inst2into3script, inst2into3tmpdir, inst3from2simbam, inst3from2simbed, inst3from2dedupb,
                            flagstat2json, bases1, bases2, args.downsample_method, args.writing_mode)
                    infodict.update({
                        'inst2into3logdir': inst2into3logdir,
                        'inst2into3script': inst2into3script,
                        'inst2into3tmpdir': inst2into3tmpdir,
                        'inst3from2datdir': inst3from2datdir,
                        'inst3from2simbam': inst3from2simbam,
                        'inst3from2simbed': inst3from2simbed,
                        'inst3from2dedupb': inst3from2dedupb,
                        'inst3from2infojs': inst3from2infojs,
                        'random_thres_str': str(random_thres_str),
                        'subclonal_frac_idxs': ','.join([str(x) for x in subclonal_frac_idxs]) # for debugging purpose
                    })
                    with cm.myopen(inst3from2infojs, args.writing_mode) as jsfile: json.dump(infodict, jsfile, indent=2)
                    for GT, acc in [('A', acc_1), ('B', acc_2)]:
                        infodict['GT'] = GT
                        infodict['accession'] = acc
                        inst1into2sh5, = find_replace_all([
                        cm.t1into2sh5], infodict)
                        ret.append((inst1into2sh5, inst2into3script))
                        ret.append((inst2into3script, inst2into3end))                    
                    ret.append((inst2into3end, F'data3from2_DSA_{infodict["donor"]}_{infodict["sampleType"]}_{infodict["avgSpotLen"]}.rule'))
                    ret.append((inst2into3end, F'data3from2_all.rule'))
                    
    return ret
if __name__ == '__main__': print(cm.list2snakemake(main()))
