import argparse, os

import common as cm
import data2from1, data3from2, data4from2and3, data_tumor

from data2from1 import NUM_CPUS
from data3from2 import cosmic_cell_lines, DOWNSAMPLE_METHODS
from data4from2and3 import SC_CN_TOOLS

EVAL_STEPS = ['2from1', '3from2', '4from2and3']

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    defaultSraRunTable = os.path.sep.join([script_dir, 'scDNAaccessions.tsv'])
    cosmic_cn_filename = os.path.sep.join([script_dir, 'cosmic-v97', 'cell_lines_copy_number.csv'])    

    parser = argparse.ArgumentParser(description='Generate bash commands to evaluate ',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('--SraRunTable', type=str, default=defaultSraRunTable, help=(
            'SraRunTable in TSV format containing the columns '
            '#Run, AvgSpotLen, Library~Name, Sample~Name, sample-type, Oocyte_ID, Donor, and SRA~Study'))
    parser.add_argument('--tumor-fastq', action='store_true', help=(
        'Treat --SraRunTable as a list of real tumor FASTQ samples (instead of near-haploid germline samples). '
        'When set, only alignment + CNV calling + clustermap are run; the haplotype-mixing simulation is skipped.'))
    parser.add_argument('-w', '--writing-mode', type=str, default=cm.DEFAULT_WRITING_MODE,
        help='File open mode for writing commands to shell script, pass any of {cm.OVERWRITING_PREVENTION_MODES} to prevent overwriting existing scripts (or w to do not prevent such thing). ')
    # 2from1
    parser.add_argument('--bwa-ncpus', type=int, default=NUM_CPUS, help='Number of CPUs used by BWA MEM ')
    # 3from2
    parser.add_argument('--cosmic',      type=str, default=cosmic_cn_filename, help='Copy-number profile TSV file downloaded from cancer.sanger.ac.uk/cosmic/download/cell-lines-project/v97')
    parser.add_argument('--cell-lines',  nargs='+',default=cosmic_cell_lines, help='Cell-lines to be used in --cosmic')
    parser.add_argument('--downsample-method', choices=DOWNSAMPLE_METHODS, default=DOWNSAMPLE_METHODS[0], help='Downsampling method') # This should not be changed
    # 4from2and3
    parser.add_argument('--tools', nargs='+', default=SC_CN_TOOLS, choices=SC_CN_TOOLS, help='Software tools calling cell-specific copy numbers from from single-cell DNA-seq data')

    parser.add_argument('--steps', nargs='+', default=EVAL_STEPS, choices=EVAL_STEPS, help='Main steps')

    args = parser.parse_args()

    ret = []
    if args.tumor_fastq:
        ret.extend(data_tumor.main(args))
        return ret

    if '2from1'     in args.steps: ret.extend(data2from1.main    (args))
    if '3from2'     in args.steps: ret.extend(data3from2.main    (args))
    if '4from2and3' in args.steps: ret.extend(data4from2and3.main(args))
    return ret

if __name__ == '__main__': print(cm.list2snakemake(main()))

