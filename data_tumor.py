#!/usr/bin/env python
"""
Pipeline for running scWGS-based CNV callers on real tumor FASTQ files.

This module is a thin tumor-mode counterpart to (data2from1 + data3from2 + data4from2and3):
  * data2from1: only the BWA-MEM alignment step is performed (no germline VCF, no haplotype
    splitting, no copy-number simulation).
  * data3from2: skipped entirely (no in-silico CN simulation).
  * data4from2and3: the per-tool CNV calling commands are reused via
    `data4from2and3.run_tool_1`, but the post-simulation evaluation step (which compares
    against simulated ground truth) is replaced by clustermap generation.

The flow per donor / sample-type / spot-length group is:
    align FASTQ -> BAM   (one BAM per tumor cell, in a persistent datdir)
    -> run each requested CNV caller on the set of BAMs
    -> aggregate per-cell CNV BEDs into a clustered heatmap (PDF + PNG) per caller.

Activated when main.py is invoked with --tumor-fastq.
"""
import argparse, logging, os
import pandas as pd

import common as cm
from common import find_replace_all, write2file
import data2from1
import data4from2and3 as d4
from data4from2and3 import (
    SC_CN_TOOLS, SC_CN_EVAL_TOOLS, SC_CN_TOOL_TO_RUN_ORDER,
    SC_CN_TOOL_DEPENDENCY_TO_DEPENDENT, ResultCN, bamfilename2samplename,
)


# Path templates for the tumor-mode files we add. We keep them under data1to2dir/<donor>
# (the same place data2from1 already writes its alignment outputs) so the rest of the
# pipeline can pick them up without changing the directory layout.
t_tumor_bam      = '<data1to2dir>/<donor>/2from1_2_<donor>.datdir/2_<donor>_2from1_tumor_<accession>_sort_markdup.bam'
t_tumor_dedupbam = '<data1to2dir>/<donor>/2from1_2_<donor>.datdir/2_<donor>_2from1_tumor_<accession>_sort_markdup_dedup.bam'

t_clmap_logdir = '<data2to4dir>/<donor>/2into4_2_<donor>_3_<sampleType>_<avgSpotLen>_4_step<tool_order>_<tool>.logdir/'
t_clmap_script = '<data2to4dir>/<donor>/2into4_2_<donor>_3_<sampleType>_<avgSpotLen>_4_step<tool_order>_<tool>.logdir/2_<donor>_3_<sampleType>_<avgSpotLen>_4_step<tool_order>_<tool>_tumor_clustermap.sh'
t_clmap_prefix = '<data2to4dir>/<donor>/4from2_2_<donor>_3_<sampleType>_<avgSpotLen>_4_step<tool_order>_<tool>_clustermap'

def _gen_alignment_rules(args, root, ref, df0, data0to1dir, data1to2dir):
    """Generate align-only rules (no germline VCF, no haplotype splitting).

    Returns (snakemake-deps list, donor_to_bam_dict) where donor_to_bam_dict maps
    donor -> list of (acc, lib, sample_name, bam_path).
    """
    deps = []
    donor_to_bams = {}
    for donor, df1 in sorted(df0.groupby('Donor')):
        infodict = {'data0to1dir': data0to1dir, 'data1to2dir': data1to2dir, 'donor': str(donor)}
        inst_sh0, inst_end, inst_log, inst_datdir = find_replace_all(
            [cm.t1into2sh0, cm.t1into2end, cm.t1into2log, cm.t2from1datdir], infodict)
        cm.makedirs((inst_log, inst_datdir))
        with cm.myopen(inst_sh0, args.writing_mode) as f0:
            write2file(F'echo {inst_sh0} is started', f0, inst_sh0)
        donor_bams = []
        for acc, LB, SM, platform, LL in zip(df1['#Run'], df1['Library~Name'], df1['Sample~Name'], df1['Platform'], df1['LibraryLayout']):
            infodict_acc = dict(infodict, accession=str(acc))
            inst_sh1, inst_fq1, inst_fq2, tumor_bam, tumor_dedupbam = find_replace_all(
                [cm.t1into2sh1, cm.t0into1fq1, cm.t0into1fq2, t_tumor_bam, t_tumor_dedupbam],
                infodict_acc)
            if LL == 'SINGLE': inst_fq2 = ''
            else: assert LL == 'PAIRED'
            with cm.myopen(inst_sh1, args.writing_mode) as f1:
                cmd_align = (
                    F'(bwa mem -R "@RG\\tID:{acc}\\tSM:{SM}\\tLB:{LB}\\tPU:L001\\tPL:ILLUMINA" '
                    F'-t {args.bwa_ncpus} {ref} {inst_fq1} {inst_fq2} '
                    F'| samtools fixmate -m - - | samtools sort -o - - '
                    F'| samtools markdup - {tumor_bam} && samtools index {tumor_bam}) #parallel=tumor.bwa/'
                )
                cmd_dedup = (
                    F'samtools view -F 0x400 -o {tumor_dedupbam} {tumor_bam} '
                    F'&& samtools index {tumor_dedupbam} #parallel=tumor.dedup/'
                )
                write2file(cmd_align + ' && ' + cmd_dedup, f1, inst_sh1)
            deps.append((inst_sh0, inst_sh1, ['resources: mem_mb = 9000']))
            deps.append((inst_sh1, inst_end))
            donor_bams.append((str(acc), str(LB), str(SM), tumor_bam, tumor_dedupbam))
        with cm.myopen(inst_end, args.writing_mode) as fe:
            write2file(F'echo {inst_end} is done', fe, inst_end)
        deps.append((inst_end, F'data2from1_donor_{donor}.rule'))
        deps.append((inst_end, F'data2from1_all.rule'))
        donor_to_bams[str(donor)] = donor_bams
    return deps, donor_to_bams


def _gen_caller_and_clustermap_rules(args, root, df0, data2to4dir, donor_to_bams, visited_scripts, phased_vcf, ref):
    """For each (avgSpotLen, sample-type, donor) group and each requested CNV tool,
    generate the run-tool snakemake rules (delegated to data4from2and3.run_tool_1) and a
    clustermap rule that aggregates the resulting per-cell CNV BEDs into a heatmap.
    """
    deps = []
    df0 = df0.copy()
    logging.info(df0[['AvgSpotLen', 'sample-type', 'Donor']])
    grouped = df0.groupby(['AvgSpotLen', 'sample-type', 'Donor', 'Platform'])
    print(f'# Start iterating over df0 with columns={list(df0.columns)}')

    for (avgSpotLen, sampleType1, donor, platform), df1 in sorted(grouped):
        sampleType = (sampleType1 + '_' + platform)
        donor_bams = donor_to_bams.get(str(donor), [])
        print(f'# Start iterating over (avgSpotLen, sampleType, donor)={(avgSpotLen, sampleType, donor)}')
        if not donor_bams:
            print(f'# Skipping (avgSpotLen, sampleType, donor)={(avgSpotLen, sampleType, donor)}')
            continue
        # Tool ordering matches the germline-mode: dependencies first, then dependents.
        tool2script_dict = {}
        for tool in sorted(args.tools, key=lambda x: (-SC_CN_TOOL_TO_RUN_ORDER[x], x)):
            infodict = {
                'donor': str(donor), 'sampleType': str(sampleType),
                'avgSpotLen': str(avgSpotLen), 'cellLine': 'tumor',
                'tool': str(tool),
                'tool_order': str(SC_CN_TOOL_TO_RUN_ORDER[tool]),
                'tool_ord_1': str(SC_CN_TOOL_TO_RUN_ORDER[tool] + 1),
                'tool_ord_2': str(SC_CN_TOOL_TO_RUN_ORDER[tool] + 2),
                'data2to4dir': data2to4dir,
            }
            (logdir, script, script2, tmpdir, datdir, clmap_logdir, clmap_script, clmap_prefix
             ) = find_replace_all(
                [cm.t2into4logdir, cm.t2into4script, cm.t2into4scrip2, cm.t2into4tmpdir,
                 cm.t4from2datdir, t_clmap_logdir, t_clmap_script, t_clmap_prefix],
                infodict)
            cm.makedirs((logdir, tmpdir, datdir, clmap_logdir))

            # Build the inbam2call mapping for this group: tumor BAMs go in, per-cell CNV
            # BEDs come out (under datdir).
            inbam2call = {}
            for acc, LB, SM, tumor_bam, tumor_dedupbam in donor_bams:
                infodict_acc = dict(infodict, accession=str(acc), samplename=None)
                # Use the standard 4from2 BED naming; only the input BAMs differ.
                depcns_tpl, intcns_tpl = find_replace_all(
                    [cm.t4from2depcns, cm.t4from2intcns], infodict_acc)
                samplename = bamfilename2samplename(tumor_bam)
                infodict_acc['samplename'] = samplename
                depcns, intcns = find_replace_all([depcns_tpl, intcns_tpl], infodict_acc)
                inbam2call[tumor_bam] = ResultCN(
                    input_bam=tumor_bam, dedup_bam=tumor_dedupbam,
                    simul_bed='', info_json='', depCN_bed=depcns, intCN_bed=intcns)

            # Sentinel start: the per-donor data2from1 end script.
            start_script = find_replace_all([cm.t1into2end], {'data1to2dir': cm.get_varnames(
                os.path.abspath(os.path.sep.join([os.path.dirname(os.path.abspath(__file__)),
                '..', 'real_tumor_data'])))[1], 'donor': str(donor)})[0]

            # Reuse the existing per-tool run/normalize logic.
            run_deps, _, _, lib2bed = d4.run_tool_1(
                infodict, tool, inbam2call, tmpdir, script, script2,
                clmap_script, root, vcf=phased_vcf, tool2script_dict=tool2script_dict,
                start_script=start_script, is_overall_haploid=False,
                writing_mode=args.writing_mode, visited_scripts=visited_scripts)
            tool2script_dict[tool] = script
            deps.extend(run_deps)

            # Add the clustermap step for evaluation tools (the ones that actually produce CN BEDs).
            if tool in SC_CN_EVAL_TOOLS:
                bed_glob = F'{datdir}*intcns.bed'
                title = F'{tool} | donor={donor} sampleType={sampleType} avgSpotLen={avgSpotLen}'
                with cm.myopen(clmap_script, args.writing_mode) as cf:
                    cmd = (F'python {root}/copy-num-bench-scwgs/cnv_clustermap.py '
                           F'-i {bed_glob} -o {clmap_prefix} --fai {ref}.fai --bin-size 50000 '
                           F'--title "{title}" #sequential=clustermap.{tool}/')
                    write2file(cmd, cf, clmap_script)
                deps.append((script2, clmap_script))
                deps.append((clmap_script, F'data4from2and3_3_clustermap_DSA_{donor}_{sampleType}_{avgSpotLen}.rule'))
                deps.append((clmap_script, F'data4from2and3_3_clustermap_tool_{tool}.rule'))
                deps.append((clmap_script, F'data4from2and3_3_clustermap_all.rule'))
    return deps


def main(args1=None):
    """Entry point used by main.py when --tumor-fastq is set."""
    logging.basicConfig(level=logging.INFO,
        format='%(asctime)s %(pathname)s:%(lineno)d %(levelname)s - %(message)s')
    script_dir = os.path.dirname(os.path.abspath(__file__))
    datadir = os.path.abspath(os.path.sep.join([script_dir, '..', 'real_tumor_data']))
    data0to1dir, data1to2dir, data2to3dir, data2to4dir, data3to4dir, data4to5dir = cm.get_varnames(datadir)
    root = os.path.abspath(F'{script_dir}/../')
    root = os.getenv('cnvguiderRoot', root)
    ref = F'{root}/refs/hg19.fa'
    ref = os.getenv('cnvguiderRef', ref)
    phased_vcf = f'{datadir}/HG008-N-P.phased.hg19.pos.tsv'

    if args1 is None:
        # Standalone use: minimal argparse stub. main.py normally fills args1 in.
        parser = argparse.ArgumentParser(description='Tumor-mode pipeline (alignment + CNV calling + clustermap)',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('--SraRunTable', type=str, required=True)
        parser.add_argument('--bwa-ncpus', type=int, default=data2from1.NUM_CPUS)
        parser.add_argument('--tools', nargs='+', default=SC_CN_TOOLS, choices=SC_CN_TOOLS)
        parser.add_argument('-w', '--writing-mode', type=str, default=cm.DEFAULT_WRITING_MODE)
        parser.add_argument('--tumor-fastq', action='store_true')
        parser.add_argument('--phased-vcf', default=phased_vcf)
        args = parser.parse_args()
    else:
        args = args1
        if not hasattr(args, 'phased_vcf') or not args.phased_vcf:
            logging.info(f"Setting phased_vcf={phased_vcf} by default.")
            args.phased_vcf = phased_vcf

    df0 = pd.read_csv(args.SraRunTable, sep='\t', header=0)

    df0.columns = df0.columns.str.replace(' ', '~')
    df0 = df0.astype(str).apply(lambda x: x.str.replace(' ', '-'))
    if '#Run' not in df0.columns: df0['#Run'] = df0['Run']
    if 'AvgSpotLen' not in df0.columns: df0['AvgSpotLen'] = 0
    df0['AvgSpotLen'] = df0['AvgSpotLen'].replace('', 0).fillna('0')
    if 'sample-type' not in df0.columns: df0['sample-type'] = df0['sample_type']
    df0['sample-type'] = cm.norm_sample_type(df0)
    if 'Donor' not in df0.columns: df0['Donor'] = df0['isolate'].str.replace(' ', '-')

    deps_align, donor_to_bams = _gen_alignment_rules(args, root, ref, df0, data0to1dir, data1to2dir)
    visited_scripts = set()
    print(f"#Going over df_with_shape={df0.shape}")
    deps_calls = _gen_caller_and_clustermap_rules(args, root, df0, data2to4dir, donor_to_bams, visited_scripts, args.phased_vcf, ref)
    print(f"#num_deps_align={len(deps_align)} num_deps_calls={len(deps_calls)}")
    return deps_align + deps_calls


if __name__ == '__main__':
    print(cm.list2snakemake(main()))
