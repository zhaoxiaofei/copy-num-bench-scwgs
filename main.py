import argparse, os

import common as cm
import data2from1, data3from2, data4from2and3, data_tumor
import gink_custom_binning
import ploidy_tools

from data2from1 import NUM_CPUS
from data3from2 import cosmic_cell_lines, DOWNSAMPLE_METHODS
from data4from2and3 import SC_CN_TOOLS
from ploidy_eval import DEFAULT_PLOIDY_WINDOW

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
    parser.add_argument('--fqs', nargs='+', default=None)
    parser.add_argument('--fastq-layout', choices=cm.FASTQ_LAYOUTS, default=cm.FASTQ_LAYOUT_DEFAULT, help=(
        'How tumor-mode input FASTQs under 1from0.datdir/ are located. '
        '"flat" (default) is the historical <accession>_1.fastq.gz with no filesystem lookup. '
        '"auto" searches with cm.fastq_pair: the sharded path if that file exists, else the '
        'flat path if that exists, else flat. "sharded" always uses the shard sub-path.'))
    parser.add_argument('--fastq-shard-template', default=cm.FASTQ_SHARD_TEMPLATE,
        help='Shard sub-path built from <study>, <accprefix> and <accession>')
    parser.add_argument('--fastq-shard-prefix-len', type=int, default=cm.FASTQ_SHARD_PREFIX_LEN,
        help='Number of leading accession characters forming <accprefix>')
    parser.add_argument('--infer-library-layout', action='store_true', help=(
        'Trust the FASTQ files on disk over the LibraryLayout column. Many SRA/ENA runs are '
        'annotated PAIRED although only single-end reads were submitted, which otherwise '
        'fails with a missing <accession>_2.fastq.gz. When set, any run whose read 2 is '
        'absent or holds no reads is aligned as SINGLE. Off by default: it stats each read 2.'))
    parser.add_argument('--tumor-datdir', default=os.path.abspath(os.path.sep.join([script_dir, '..', 'real_tumor_data'])))
    parser.add_argument('--ploidy-file', default=None, help=data_tumor.PLOIDY_FILE_HELP)
    parser.add_argument('--ploidy-window', type=float, default=DEFAULT_PLOIDY_WINDOW,
        help=data_tumor.PLOIDY_WINDOW_HELP)
    parser.add_argument('--ploidy-chroms', choices=['autosomes', 'all'], default='autosomes',
        help=data_tumor.PLOIDY_CHROMS_HELP)
    parser.add_argument('--ploidy-tools', nargs='*', default=[], choices=ploidy_tools.TOOLS,
        help=ploidy_tools.PLOIDY_TOOLS_HELP)
    parser.add_argument('--ploidy-facs', action='store_true', help=ploidy_tools.PLOIDY_FACS_HELP)

    parser.add_argument('--donor', default='tumor')
    parser.add_argument('--sampleType', default='tumor')
    parser.add_argument('--avgSpotLen', type=int, default=0)
    parser.add_argument('--phased-vcf', default=None, help='Phased VCF file required by haplotype-aware CNV callers such as Chisel.')
    parser.add_argument('-w', '--writing-mode', type=str, default=cm.DEFAULT_WRITING_MODE,
        help=F'File open mode for writing commands to shell script, pass any of {cm.OVERWRITING_PREVENTION_MODES} to prevent overwriting existing scripts (or w to do not prevent such thing). ')
    # 2from1
    parser.add_argument('--bwa-ncpus', type=int, default=NUM_CPUS, help='Number of CPUs used by BWA MEM ')
    # 3from2
    parser.add_argument('--cosmic',      type=str, default=cosmic_cn_filename, help='Copy-number profile TSV file downloaded from cancer.sanger.ac.uk/cosmic/download/cell-lines-project/v97')
    parser.add_argument('--cell-lines',  nargs='+',default=cosmic_cell_lines, help='Cell-lines to be used in --cosmic')
    parser.add_argument('--downsample-method', choices=DOWNSAMPLE_METHODS, default=DOWNSAMPLE_METHODS[0], help='Downsampling method') # This should not be changed
    # 4from2and3
    parser.add_argument('--tools', nargs='+', default=SC_CN_TOOLS, choices=SC_CN_TOOLS + [gink_custom_binning.TOOL], help='Software tools calling cell-specific copy numbers from from single-cell DNA-seq data')
    parser.add_argument('--excluded-tools', nargs='+', default=[], choices=SC_CN_TOOLS,
        help='Tools from SC_CN_TOOLS to exclude, together with any tools that depend on them')
    parser.add_argument('--binnings', nargs='+', default=[], help='Additional Ginkgo binning options')
    parser.add_argument('--steps', nargs='+', default=EVAL_STEPS, choices=EVAL_STEPS, help='Main steps')

    args = parser.parse_args()
    if args.fqs: args.tumor_fastq = True
    # Drop --excluded-tools plus every tool that (transitively) depends on them.
    excluded = set(args.excluded_tools)
    while True:
        dependents = {dependent for dep, dependents in data4from2and3.SC_CN_TOOL_DEPENDENCY_TO_DEPENDENT.items()
                      if dep in excluded for dependent in dependents}
        if dependents <= excluded: break
        excluded |= dependents
    args.tools = [tool for tool in args.tools if tool not in excluded]
    try: args.tools = gink_custom_binning.setup(data4from2and3, args.tools, args.binnings)
    except ValueError as exc: parser.error(str(exc))
    # Opt-in: without --ploidy-tools this returns args.tools as it is and leaves both
    # data4from2and3 and the generated workflow completely unchanged.
    try: args.tools = ploidy_tools.setup(data4from2and3, args.tools, args.ploidy_tools, args)
    except ValueError as exc: parser.error(str(exc))

    ret = []
    if args.tumor_fastq:
        ret.extend(data_tumor.main(args))
        return ret

    if '2from1'     in args.steps: ret.extend(data2from1.main    (args))
    if '3from2'     in args.steps: ret.extend(data3from2.main    (args))
    if '4from2and3' in args.steps: ret.extend(data4from2and3.main(args))
    return ret
if __name__ == '__main__': print(cm.list2snakemake(main()))
