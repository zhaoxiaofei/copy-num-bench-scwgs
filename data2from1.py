#!/usr/bin/env python

import argparse, logging, os, sys
import pandas as pd

from types import SimpleNamespace

import common as cm
from common import find_replace_all, write2file
 
NUM_CPUS = 4
MAX_HAPLO_CN = 6

## Requirements to run the script generated to stdout
##   bwa, samtools, hg19
##   bcftools, eagleimp with eagleimp's database in the directory EAGLE_IMP_DB_DIR (if running chisel is required)

# >>> accessions=$(csvformat -T < /mnt/d/sperm/data/PMC7923680-SRP194057.SraRunTable.txt | tr ' ' '~' | grep single~cell | grep Sperm~Cell | awk '{print $1}')
# >>> ACC_LB_SM=$(csvformat -T < /mnt/d/sperm/data/PMC7923680-SRP194057.SraRunTable.txt | tr ' ' '~' | grep single~cell | grep Sperm~Cell | awk '{print $1, $19, $29}')
# >>> # accessions="SRR8981905 SRR8981918 SRR8981919 SRR8981920 SRR8981921 SRR8981922 SRR8981923 SRR8981924 SRR8981925 SRR8981928 SRR8981953 SRR8981954 SRR8981957 SRR8981958 SRR8981959 SRR8981960 SRR8981961 SRR8981962 SRR8981963 SRR8981964 SRR8981985 SRR8981986 SRR8981988 SRR8981989 SRR8981990 SRR8981991 SRR8981992 SRR8982030 SRR8982031 SRR8982034 SRR8982035 SRR8982036 SRR8982037 SRR8982038 SRR8982039 SRR8982040 SRR8982041 SRR8982064"

def main(args1=None):
    ret = []
    
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(pathname)s:%(lineno)d %(levelname)s - %(message)s')

    script_dir = os.path.dirname(os.path.abspath(__file__))
    datadir = os.path.abspath(os.path.sep.join([script_dir, '..', 'data']))
    data0to1dir, data1to2dir, data2to3dir, data2to4dir, data3to4dir, data4to5dir  = cm.get_varnames(datadir)

    defaultSraRunTable = os.path.sep.join([script_dir, 'scDNAaccessions.tsv'])

    root = os.path.abspath(F'{script_dir}/../')
    root = os.getenv('cnvguiderRoot', root)
    ref = F'{root}/refs/hg19.fa'
    ref = os.getenv('cnvguiderRef', ref)
    
    parser = argparse.ArgumentParser(description='Generate bash commands to in-silico mix the single-cell sequencing data with copy numbers simulated from the COSMIC cell line database. ',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--SraRunTable', type=str, default=defaultSraRunTable, help=(
            'SraRunTable in TSV format containing the columns '
            '#Run, AvgSpotLen, Library~Name, Sample~Name, sample-type, Oocyte_ID, Donor, and SRA~Study'))
    parser.add_argument('--bwa-ncpus', type=int, default=NUM_CPUS, help='Number of CPUs used by BWA MEM ')
    parser.add_argument('-w', '--writing-mode', type=str, default=cm.DEFAULT_WRITING_MODE,
        help='File open mode for writing commands to shell script, pass any of {cm.OVERWRITING_PREVENTION_MODES} to prevent overwriting existing scripts (or w to do not prevent such thing). ')
    args = (args1 if args1 else parser.parse_args())
    
    df0 = pd.read_csv(args.SraRunTable, sep='\t')
    grouped = df0.groupby('Donor')
    partitioned_dfs = {donor: df1 for donor, df1 in grouped}
    
    for donor, df1 in sorted(partitioned_dfs.items()):
        infodict = {'data0to1dir': data0to1dir, 'data1to2dir': data1to2dir, 'donor': donor}
        inst1into2sh0, inst1into2sh2, inst1into2end, inst1into2log, inst1into2tmp, inst2from1datdir = find_replace_all([
        cm.t1into2sh0, cm.t1into2sh2, cm.t1into2end, cm.t1into2log, cm.t1into2tmp, cm.t2from1datdir], infodict)
        
        cm.makedirs((inst1into2log, inst1into2tmp, inst2from1datdir))
        with cm.myopen(inst1into2sh0, args.writing_mode) as shfile0: write2file(F'echo {inst1into2sh0} is started', shfile0, inst1into2sh0)

        for acc, LB, SM in zip(df1['#Run'], df1['Library~Name'], df1['Sample~Name']):
            infodict['accession'] = acc
            inst1into2sh1, inst0into1fq1, inst0into1fq2 = find_replace_all([
            cm.t1into2sh1, cm.t0into1fq1, cm.t0into1fq2], infodict)
            
            ret.append((inst1into2sh0, inst1into2sh1, ['resources: mem_mb = 9000']))
            ret.append((inst1into2sh1, inst1into2sh2, ['resources: mem_mb = 9000']))
            with cm.myopen(inst1into2sh1, args.writing_mode) as shfile1:
                bam2 = inst1into2tmp + F'/{acc}_12_sort_markdup.bam'
                cmd = (F'''(bwa mem -R "@RG\\tID:{acc}\\tSM:{SM}\\tLB:{LB}\\tPU:L001\\tPL:ILLUMINA" -t {args.bwa_ncpus} {ref} {inst0into1fq1} {inst0into1fq2} '''
                      +F''' | samtools fixmate -m - - | samtools sort -o - - | samtools markdup - {bam2} && samtools index {bam2}) #parallel=bwa/''')
                write2file(cmd, shfile1, inst1into2sh1)
        
        germvcf1, germvcf2, germvcf2txt, germvcf3a, germvcf3b = find_replace_all(
            [cm.t2from1vcf010, cm.t2from1vcf020, cm.t2from1vcf021, cm.t2from1vcf03A, cm.t2from1vcf03B], 
            {'data1to2dir': data1to2dir, 'donor': donor}
        )
        cmd_callsnps = (
                F''' samtools merge --threads 4 -u - {inst1into2tmp}/*_12_sort_markdup.bam '''
               +F''' | bcftools mpileup --threads 4 -a "-FORMAT/AD,INFO/AD,INFO/ADF,INFO/ADR" --ignore-RG -Ou -f {ref} - | bcftools call --ploidy GRCh37 -mv -Oz -o {germvcf1}'''
               +F''' && bcftools index -ft {germvcf1}'''
               +F''' && bcftools norm -m -any {germvcf1} '''
               +F''' | bcftools filter -s q10       -e "QUAL<10" '''
               +F''' | bcftools filter -s lowADtoDP -e "(INFO/AD[1])/(INFO/DP) < 0.05" ''' 
               +F''' | bcftools filter -s lowAD     -e "INFO/AD[1] < 3 " '''
               +F''' | awk 'BEGIN {{ sample1="HaplotypeA"; sample2="HaplotypeB"; srand(0); }} '''
               +F''' /^#CHROM/ {{print $0"\t"sample1"\t"sample2"\tDiplotype"; next}} '''
               +F''' /^#/ {{print $0; next}} '''
               +F''' {{ gt1 = int(rand() * 2) ; gt2 = int(rand() * 2); print $0"\t"gt1"\t"gt2"\t"gt1"|"gt2 }}' '''
               +F''' | bcftools view -Oz -o {germvcf2} && bcftools view -s "Diplotype" -Ov -o {germvcf2txt} --targets {cm.chrsNS} {germvcf2} ''')
        for GT, germvcf3 in  [('A', germvcf3a), ('B', germvcf3b)]:
            cmd = F'''bcftools view -fPASS {germvcf2} -s Haplotype{GT} -Oz -o {germvcf3} '''
            cmd_callsnps += ' && ' + cmd
        with cm.myopen(inst1into2sh2, args.writing_mode) as shfile2: write2file(cmd_callsnps, shfile2, inst1into2sh2)

        for GT, germvcf3 in [('A', germvcf3a), ('B', germvcf3b)]:
            infodict['GT'] = GT
            for acc, LB, SM in zip(df1['#Run'], df1['Library~Name'], df1['Sample~Name']):
                infodict['accession'] = acc
                inst1into2sh3, inst1into2sh4, inst1into2sh5, inst2from1mutbam, inst2from1cnbams, inst2from1flagjs = find_replace_all([
                cm.t1into2sh3, cm.t1into2sh4, cm.t1into2sh5, cm.t2from1mutbam, cm.t2from1cnbams, cm.t2from1flagjs], infodict)
                
                ret.extend([(inst1into2sh2, inst1into2sh3), (inst1into2sh3, inst1into2sh4, ['resources: mem_mb = 9000']), (inst1into2sh4, inst1into2sh5, ['resources: mem_mb = 9000'])])
                with    cm.myopen(inst1into2sh3, args.writing_mode) as shfile3, \
                        cm.myopen(inst1into2sh4, args.writing_mode) as shfile4, \
                        cm.myopen(inst1into2sh5, args.writing_mode) as shfile5:
                    bam2 = F'{inst1into2tmp}/{acc}_12_sort_markdup.bam'
                    fq20 = F'{inst1into2tmp}/{acc}_gt{GT}_0.fastq.gz'
                    fq21 = F'{inst1into2tmp}/{acc}_gt{GT}_1.fastq.gz'
                    fq22 = F'{inst1into2tmp}/{acc}_gt{GT}_2.fastq.gz'
                    fq21sorted = F'{inst1into2tmp}/{acc}_gt{GT}_1_sorted.fastq.gz'
                    fq22sorted = F'{inst1into2tmp}/{acc}_gt{GT}_2_sorted.fastq.gz'
                    cmd3 = F'{root}/copy-num-bench-scwgs/data1to2code/safesim/bin/safemut -b {bam2} -v {germvcf3} -1 {fq21} -2 {fq22} -0 {fq20} #parallel=safemut/'
                    write2file(cmd3, shfile3, inst1into2sh3)
                    for fqidx, (fq, fqsorted) in enumerate([(fq21, fq21sorted), (fq22, fq22sorted)]):
                        tmpdir = F'/tmp/{acc}_gt{GT}_{fqidx+1}.sorted.tmpdir'
                        cmd4 = F'mkdir -p {tmpdir} && {root}/copy-num-bench-scwgs/data1to2code/safesim/bin/fastq-sort -n -S 4G --temporary-directory={tmpdir} <(zcat {fq}) | gzip --fast > {fqsorted} #parallel=fastq-sort/'
                        write2file(cmd4, shfile4, inst1into2sh4)
                    cn_bams_str = ' '.join(inst2from1cnbams)
                    cmd5= (F'''(bwa mem -R "@RG\\tID:{acc}\\tSM:{SM}\\tLB:{LB}\\tPU:L001\\tPL:ILLUMINA" -t {args.bwa_ncpus} {ref} {fq21sorted} {fq22sorted} ''' 
                          +F'''| samtools fixmate -m - - | samtools sort -o - - | samtools markdup - {inst2from1mutbam} && samtools index {inst2from1mutbam})'''
                          +F''' && samtools flagstat -O json {inst2from1mutbam} > {inst2from1flagjs} '''
                          +F''' && {root}/copy-num-bench-scwgs/data1to2code/splitbam.out {inst2from1mutbam} {MAX_HAPLO_CN} {cn_bams_str} && for b in {cn_bams_str}; do samtools index $b ; done #parallel=safemut.bwa/''')
                    write2file(cmd5, shfile5, inst1into2sh5)
                    ret.append((inst1into2sh5, inst1into2end))
        with cm.myopen(inst1into2end, args.writing_mode) as file_end: write2file(F'echo {inst1into2end} is done', file_end, inst1into2end)
        ret.append((inst1into2end, F'data2from1_donor_{donor}.rule'))
        ret.append((inst1into2end, F'data2from1_all.rule'))
    return ret

if __name__ == '__main__': print(cm.list2snakemake(main()))
