#!/usr/bin/env python

import argparse, io, json, logging, os, random, sys
import pandas as pd

from collections.abc import Iterable
from types import SimpleNamespace

import common as cm
from common import find_replace_all, get_varnames, change_file_ext, write2file, norm_sample_type

SC_CN_TOOL_DEPENDENCY_TO_DEPENDENT = {
    'readcounter'   : {'hmmcopy': '', 'copynumber': ''},
    'bam2bed'       : {'ginkgo' : ''},
    'mq30bedcov500k': {'sccnv'  : ''},
    'dedup'         : {'secnv'  : ''},
    'nop' : {
    'scyn'          : '',
    'chisel'        : '', 
    'flcna'         : '',
    'aneufinder'    : '',
    # 'alleloscope'   : '',
    }
}
# The following CNV callers were affected by runtime errors that were described in the github links. 
# 'alleloscope', # https://github.com/seasoncloud/Alleloscope/issues/16 # multiple types of errors
# https://github.com/EngeLab/ASCENT/issues/2
# https://github.com/elkebir-group/CNRein/issues/9 # very slow runtime too
# https://github.com/parklab/HiScanner/issues/11   # very slow runtime too
# https://github.com/zhyu-lab/cot/issues/1         # compiling error
# https://github.com/zhyu-lab/rccae/issues/1       # compiling error
# AFAIK, there are no other scWGS-based CNV callers.

def dict_append(d, to_append):
    if isinstance(d, dict): return {k: dict_append(v, to_append) for k, v in sorted(d.items())}
    elif isinstance(d, list): return [dict_append(subd, to_append) for subd in d]
    else: return {d: to_append}
def dict_rev_recursively(d):
    if isinstance(d, dict): return [dict_append(dict_rev_recursively(v), k) for k, v in sorted(d.items())]
    else: return d

def dict_simplify(d):
    if isinstance(d, list) and len(d) == 1:            return dict_simplify(d[0])
    if isinstance(d, list):                            return [dict_simplify(subd) for subd in d]
    if isinstance(d, dict) and [''] == list(d.keys()): return dict_simplify(list(d.values())[0])
    if isinstance(d, dict):                            return {k: dict_simplify(v) for k, v in sorted(d.items())}
    return d


SC_CN_EVAL_TOOLS = set(['hmmcopy', 'copynumber', 'ginkgo', 'sccnv', 'secnv', 'scyn', 'chisel', 'flcna', 'aneufinder'])

SC_CN_TOOL_TO_RUN_MODE = {
    'readcounter' : 'parallel',
    'hmmcopy'     : 'parallel',
    'bam2bed'     : 'parallel',
    'ginkgo'      : 'sequential',
    'copynumber'  : 'parallel',
    'sccnv'       : 'sequential',
    'secnv'       : 'sequential',
    'scyn'        : 'sequential',
    'chisel'      : 'sequential',
    'flcna'       : 'sequential',
    'aneufinder'  : 'sequential',
    'alleloscope' : 'sequential',
}

chrs = 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY'

class ResultCN:
    def __init__(self, input_bam, dedup_bam, simul_bed, info_json, depCN_bed, intCN_bed):
        self.input_bam = input_bam
        self.dedup_bam = dedup_bam
        self.simul_bed = simul_bed
        self.info_json = info_json
        self.depCN_bed = depCN_bed
        self.intCN_bed = intCN_bed
        self.lib = bamfilename2samplename(input_bam)

def from_DEPENDENCY_TO_DEPENDENT_to_TOOL_TO_RUN_ORDER(d, order=1, avoided_str=''):
    if isinstance(d, str): return {d: order} if (d != avoided_str) else {}
    else:
        ret = {}
        for k,v in sorted(d.items()):
            if k != avoided_str: ret[k] = max([ret.get(k, -1), order])
            ret.update(from_DEPENDENCY_TO_DEPENDENT_to_TOOL_TO_RUN_ORDER(v, order+1, avoided_str))
        return ret

def norm_RUN_ORDER(d): return {k : 
    #(max((v, 2 + (1 if k in [' chisel '] else 0))) if k in SC_CN_EVAL_TOOLS else v) 
    v for (k, v) in sorted(d.items())}

SC_CN_TOOL_TO_RUN_ORDER = norm_RUN_ORDER(from_DEPENDENCY_TO_DEPENDENT_to_TOOL_TO_RUN_ORDER(SC_CN_TOOL_DEPENDENCY_TO_DEPENDENT))
SC_CN_TOOLS = sorted(SC_CN_TOOL_TO_RUN_ORDER.keys())

## Requirements to run the script generated to stdout
##   bwa, samtools, hg19
##   bcftools, eagleimp with eagleimp's database in the directory EAGLE_IMP_DB_DIR (if running chisel is required)

def norm_tool_result(call, result_filename_1, tool, lib_1, pyscripts, result_filename_2=None, lib_from_tool=None, ref=None):
    # (infodict, t4fromXdepcns, t4fromXintcns, result_filename_1, tool, lib_1, tobed, result_filename_2=None, lib_from_tool=None):
    if not result_filename_2: result_filename_2 = result_filename_1
    if not lib_from_tool: lib_from_tool = lib_1
    toBed, cnv_norm_with_prior = pyscripts
    # infodict['bamfilename'] = lib_1
    # depcns, intcns = find_replace_all([t4fromXdepcns, t4fromXintcns], infodict)
    depcns, intcns = (call.depCN_bed, call.intCN_bed)
    if cnv_norm_with_prior:
        cmdnorm = F'cat {depcns} | python {cnv_norm_with_prior} > {intcns}.overall_haploid'
    else:
        cmdnorm = F'echo Skip generating {intcns}.overall_haploid'
    add_args = (f'--fai {ref}.fai ' if ref else '')
    cmd1 = F'cat {result_filename_1} | python {toBed} --caller {tool} --sample {lib_from_tool} {add_args}--dp 1> {depcns} && {cmdnorm} #parallel=eval.{tool}.{lib_1}.by.dep/'
    cmd2 = F'cat {result_filename_2} | python {toBed} --caller {tool} --sample {lib_from_tool} {add_args}     1> {intcns} #parallel=eval.{tool}.{lib_1}.by.int/'
    
    return [cmd1] + [cmd2], [depcns, intcns]

def bamfilename2samplename(bam): return bam.split(os.path.sep)[-1].split('.')[0]

def get_cleanup(script): return change_file_ext(script, 'cleanup.sh')

def run_tool_1(infodict, tool, inbam2call, tmpdir, script, script2, script_eval, rootdir, vcf, tool2script_dict, start_script, is_overall_haploid, writing_mode, visited_scripts, normal_bams_dir=None):
    cellLine = ('hap' if is_overall_haploid else infodict['cellLine'])
    
    ref=F'{rootdir}/refs/hg19.fa'
    bigwig=F'{rootdir}/refs/wgEncodeCrgMapabilityAlign36mer.bigWig'
    bigwig100=F'{rootdir}/refs/wgEncodeCrgMapabilityAlign100mer.bigWig'
    window_size = 200*1000
    readcounter_bin       = F'{rootdir}/cnvguider/data3to4code/hmmcopy_utils/bin/readCounter'
    correct_read_count_py = F'{rootdir}/cnvguider/data3to4code/single_cell_pipeline-0.8.26/single_cell/workflows/hmmcopy/scripts/correct_read_count.py'
    hmmcopy_R             = F'{rootdir}/cnvguider/data3to4code/single_cell_pipeline-0.8.26/single_cell/workflows/hmmcopy/scripts/hmmcopy_single_cell.R'
    ginkgo_sh             = F'{rootdir}/cnvguider/data3to4code/ginkgo/cli/ginkgo.sh'
    tobed                 =(F'{rootdir}/cnvguider/cnv_raw_to_bed.py', (F'{rootdir}/cnvguider/cnv_norm_with_prior.py' if is_overall_haploid else ''))
    
    bams       = sorted(inbam2call.keys())
    subscripts = [change_file_ext(bam, F'{tool}.sh') for bam in bams]
    bais       = [change_file_ext(bam, 'bam.bai') for bam in bams]
    beds       = [change_file_ext(bam, 'bed'    ) for bam in bams]
    mapq30bams = [change_file_ext(bam, 'mq30.bam')        for bam in bams]
    mapq30bais = [change_file_ext(bam, 'mq30.bam.bai')    for bam in bams]
    mapq30beds = [change_file_ext(bam, 'mq30.depth.bin500000.bed') for bam in bams] # This extension has to be hard coded for sccnv_v1.0.2.py
    dedup_bams = [inbam2call[bam].dedup_bam       for bam in bams]       
    dedup_bais = [change_file_ext(bam, 'bam.bai') for bam in dedup_bams]
    libs       = [bamfilename2samplename(bam)     for bam in bams]
    cnvs       = libs
    depcn_beds = [inbam2call[bam].depCN_bed       for bam in bams]
    intcn_beds = [inbam2call[bam].intCN_bed       for bam in bams]
    cmds, cmds2 = [], []
    bam2bed = {}
    lib2bed = {}
    
    sccnv_py = F'{rootdir}/cnvguider/data3to4code/SCCNV/sccnv_v1.0.2.py'
    sccnv_dir = os.path.dirname(os.path.abspath(sccnv_py))
    
    hg19bed = F'{sccnv_dir}/resource/hg19.bin500000.bed'

    deps = []
    # This should be already done in the setup
    #if tool == 'setup':
    #    cmds.append(F'time -p mapCounter -w {window_size} {bigwig} -c {chrs} > {ref}.mp.seg #parallel=setup.ref/')
    #    cmds.append(F'time -p gcCounter -w {window_size} {ref} -c {chrs} > {ref}.gc.seg #parallel=setup.ref/')
    #    cmds.append(F'time -p bedtools makewindows -g {ref}.fai -w 250000 > {ref}.250000.bed #parallel=setup.ref/')
    #    cmds.append(F'time -p bedtools makewindows -g {ref}.fai -w 500000 > {ref}.500000.bed #parallel=setup.ref/')
    if tool == 'nop':
        deps.append((start_script, script))
        cmds.append(F'echo performed no-operation')
        for tool_next in SC_CN_TOOL_DEPENDENCY_TO_DEPENDENT[tool]:
            script_next = tool2script_dict[tool_next]
            # These mem_mb's are empirically measured with "command time -v"
            if tool_next == 'scyn':
                deps.append((script, script_next, ['resources: mem_mb = 5000']))
            elif tool_next == 'chisel':
                deps.append((script, script_next, ['resources: mem_mb = 18000']))
            elif tool_next == 'alleloscope':
                deps.append((script, script_next, ['resources: mem_mb = 10000']))
            elif tool_next == 'flcna':
                deps.append((script, script_next, ['resources: mem_mb = 40000']))
            elif tool_next == 'aneufinder':
                deps.append((script, script_next, ['resources: mem_mb = 8000']))
            else:
                deps.append((script, script_next))
    elif tool in SC_CN_TOOL_DEPENDENCY_TO_DEPENDENT:
        assert len(subscripts) == len(set(subscripts)), F'The subscript files {sorted(subscripts)} are duplicated!'
        for bam, bed, lib, mq30bam, mq30bed, dedup_bam, subscript in zip(bams, beds, libs, mapq30bams, mapq30beds, dedup_bams, subscripts):
            assert bed != bam, F'{bam} != {bed} failed'
            if tool == 'readcounter':
                cmd = (F'time -p {readcounter_bin} -w {window_size} {bam} -c {chrs} > {bam}.wig ' 
                       F' && time -p conda run -n single_cell_pipeline-0.8.26 python {correct_read_count_py} {ref}.gc.seg {ref}.mp.seg {bam}.wig {bam}.wig.csv '
                       F' #parallel=prerun.{tool}/')
                cmd_cleanup = F'echo Keeping {bam}.wig and {bam}.wig.csv #parallel=cleanup.{tool}'
            elif tool == 'mq30bedcov500k':
                cmd = (F"samtools view -bh -q 30 {bam} { chrs.replace(',', ' ') } > {mq30bam} && samtools index {mq30bam} "
                       F" && samtools bedcov {hg19bed} {mq30bam} > {mq30bed} #parallel=prerun.{tool}/")
                cmd_cleanup = F'rm {mq30bam} && echo Keeping {mq30bed} #parallel=cleanup.{tool}/'
            elif tool == 'bam2bed':                
                cmd = F' time -p bedtools bamtobed -i {bam} > {bed} #parallel=prerun.{tool}/'
                cmd_cleanup = F'echo Keeping {bed} #parallel=cleanup.{tool}/'
            elif tool == 'dedup':
                cmd = F'samtools view -F 0x400 -o {dedup_bam} {bam} && samtools index {dedup_bam} #parallel=prerun.{tool}/'
                cmd_cleanup = F'rm {dedup_bam} #parallel=cleanup.{tool}/'
            else: 
                sys.stderr.write(F'The tool {tool} is unknown, aborting!')
                sys.exit(1)
            if subscript in visited_scripts:
                logging.info(F'  Skip generating the script {subscript} because it has already been generated. ')
            else:
                with cm.myopen(subscript, writing_mode) as file: write2file(cmd, file, subscript)
                visited_scripts.add(subscript)
                with cm.myopen(get_cleanup(subscript), writing_mode) as file: write2file(cmd_cleanup, file, get_cleanup(subscript))
            deps.append((start_script, subscript))
            deps.append((subscript, script))
            deps.append((get_cleanup(subscript), F'data4from2and3_4_cleanup_DSA_{infodict["donor"]}_{infodict["sampleType"]}_{infodict["avgSpotLen"]}.rule'))
            deps.append((get_cleanup(subscript), F'data4from2and3_4_cleanup_cellLine_{cellLine}.rule'))
            deps.append((get_cleanup(subscript), F'data4from2and3_4_cleanup_tool_{tool}.rule'))
            deps.append((get_cleanup(subscript), F'data4from2and3_4_cleanup_all.rule'))

        cmds.append(F'echo performed {tool}')
        for tool_next in SC_CN_TOOL_DEPENDENCY_TO_DEPENDENT[tool]:
            script_next = tool2script_dict[tool_next]
            deps.append((script, script_next))
            for subscript in subscripts:
                deps.append((script_next, get_cleanup(subscript)))

    if tool == 'hmmcopy':
        for bam, lib, cnv in zip(bams, libs, cnvs):
            cmd = (
            F's=1000 && e=0.999999 && nu=2.1 && min_mapqual=20 '
            F'&& time -p conda run -n hmmcopy Rscript {hmmcopy_R} --corrected_data={bam}.wig.csv --outdir={tmpdir} --sample_id={lib} '
            F'--param_multiplier="1,2,3,4,5,6" --param_mu="0,1,2,3,4,5,6,7,8,9,10,11" --param_k="100,100,700,100,25,25,25,25,25,25,25,25" --param_m="0,1,2,3,4,5,6,7,8,9,10,11" '
            F'--param_str=$s --param_e=$e --param_l=20 --param_nu=$nu --param_eta=50000 --param_g=3 --param_s=1 #parallel=run.{tool}/'
            )
            cmds.append(cmd)
            norm_cmds, norm_beds = norm_tool_result(inbam2call[bam], F'{tmpdir}/0/segs.csv', tool, lib, tobed) # infodict, t4fromXdepcns, t4fromXintcns, result_filename_1, t
            cmds2.extend(norm_cmds)
            bam2bed[bam] = norm_beds[1]
            lib2bed[lib] = norm_beds[1]
    if tool == 'ginkgo':
        cmd = (
            F'rm {tmpdir}/*.bed || true && cp -s {" ".join(beds)} {tmpdir}/ && ls {tmpdir}/*.bed > {tmpdir}/cells.list'
            F'&& time -p bash -evx {ginkgo_sh} --input {tmpdir} --genome hg19 --binning variable_175000_48_bwa --cells {tmpdir}/cells.list #sequential=run.{tool}/'
        )
        cmds.append(cmd)
        for bam, lib, cnv in zip(bams, libs, cnvs):
            norm_cmds, norm_beds = norm_tool_result(inbam2call[bam], F'{tmpdir}/SegFixed', tool, lib, tobed, result_filename_2=F'{tmpdir}/SegCopy', lib_from_tool=lib.replace('-', '.'))
            cmds2.extend(norm_cmds)
            bam2bed[bam] = norm_beds[1]
            lib2bed[lib] = norm_beds[1]
    if tool == 'copynumber':
        gamma = 40
        collect_rc_4copynumber_py = F'{rootdir}/cnvguider/data3to4code/collect_rc_4copynumber.py'
        CopyNumber_R = F'{rootdir}/cnvguider/data3to4code/CopyNumber.R'
        cmd = ' && '.join([F'cp {bam}.wig.csv {tmpdir}/{lib}.leaf.wig.csv' for bam, lib in zip(bams, libs)])
        # This code maintains compatibilty with the originally documented code from SingleCellCNABenchmark
        #for bam, lib in zip(bams, libs): 
        #    cmd = F'cp {bam}.wig.csv {tmpdir}/{lib}.leaf.wig.csv #parallel=prerun.{tool}/'
        #    cmds.append(cmd)
        cmd += F' && time -p python {collect_rc_4copynumber_py} {tmpdir} && time -p Rscript {CopyNumber_R} {tmpdir}/copynumber.input.csv {tmpdir}/copynumber.output {gamma} #sequential=run.{tool}/'
        cmds.append(cmd)
        for lib, cnv, bam in zip(libs, cnvs, bams):        
            norm_cmds, norm_beds = norm_tool_result(inbam2call[bam], F'{tmpdir}/copynumber.output.csv', tool, lib, tobed, lib_from_tool=lib.replace('-', '.'))
            cmds2.extend(norm_cmds)
            bam2bed[bam] = norm_beds[1]
            lib2bed[lib] = norm_beds[1]    
    if tool == 'sccnv':
        # The following cmd was already run during the installation of ScCNV
        # cmd = F'bedtools makewindows -g {sccnv_dir}/resource/hg19.chrlength.txt  -w 500000 > {hg19bed}'
        #for bam, lib in zip(bams, libs):
        #    cmd = (F'rm -r {tmpdir} || true && mkdir -p {tmpdir}/bam_mapq30/cells'
        #           F" && samtools view -bh -q 30 {bam} { chrs.replace(',', ' ') } > {tmpdir}/bam_mapq30/{lib}.bam && samtools index {tmpdir}/bam_mapq30/{lib}.bam "
        #           F" && samtools bedcov {hg19bed} {tmpdir}/bam_mapq30/{lib}.bam > {tmpdir}/bam_mapq30/cells/{lib}.depth.bin500000.bed #parallel=prerun.{tool}/")
        #    cmds.append(cmd)
        cmd = (F' rm -r {tmpdir} || true && mkdir -p {tmpdir}/bam_mapq30/cells && cp -s {" ".join(mapq30bams + mapq30bais)} {tmpdir}/bam_mapq30/ && cp -s {" ".join(mapq30beds)} {tmpdir}/bam_mapq30/cells/'
               F' && ls {tmpdir}/bam_mapq30/*.bam > {tmpdir}/bam_mapq30.bamlist.txt '
               F' && pushd {sccnv_dir} && time -p python {sccnv_py} -i {tmpdir}/bam_mapq30.bamlist.txt -o {tmpdir}/bam_mapq30/cells -k False -r True && popd #sequential=run.{tool}/')
        cmds.append(cmd)
        for lib, cnv, bam in zip(libs, cnvs, bams):        
            norm_cmds, norm_beds = norm_tool_result(inbam2call[bam], F'{tmpdir}/bam_mapq30/cells/result.dat6_cnvsmooth.txt', tool, lib, tobed)
            cmds2.extend(norm_cmds)
            bam2bed[bam] = norm_beds[1]
            lib2bed[lib] = norm_beds[1]    
    if tool == 'secnv':
        secnv_py = F'{rootdir}/cnvguider/data3to4code/SeCNV/Scripts/SeCNV.py'
        secnv_dir = os.path.dirname(os.path.abspath(secnv_py))
        cmd = (F'''rm -r {tmpdir}/bam_dedup/ || true && mkdir -p {tmpdir}/bam_dedup/ && cp -s {' '.join(dedup_bams + dedup_bais)} {tmpdir}/bam_dedup/ '''
               F' && pushd {secnv_dir} && mkdir -p {tmpdir}/output/ && export PATH="{secnv_dir}:${{PATH}}" && time -p python {secnv_py} {tmpdir}/bam_dedup {tmpdir}/output {ref} #sequential=run.{tool}/')
        cmds.append(cmd)
        for lib, cnv, bam in zip(libs, cnvs, bams):        
            norm_cmds, norm_beds = norm_tool_result(inbam2call[bam], F'{tmpdir}/output/genome_cov.bed', tool, lib, tobed, result_filename_2=F'{tmpdir}/output/cnv_matrix.csv')
            cmds2.extend(norm_cmds)
            bam2bed[bam] = norm_beds[1]
            lib2bed[lib] = norm_beds[1]
    # SCOPE1 has poorly written docs and is superseded by SCOPE2, so SCOPE1 is not run.
    # SCOPE2 is integrated into scyn, so scyn is run instead of SCOPE2
    if tool == 'scyn':
        cmd = (F'rm -r {tmpdir}/scyn_input/ || true && mkdir -p {tmpdir}/scyn_input/ {tmpdir}/scyn_output/ && cp -s {" ".join(bams + bais)} {tmpdir}/scyn_input/ '
               F' && time -p python {rootdir}/cnvguider/data3to4code/simplerun_scyn.py {tmpdir}/scyn_input/ {tmpdir}/scyn_output/ {tmpdir}/scyn_output.csv #sequential=run.{tool}/')
        cmds.append(cmd)
        for lib, cnv, bam in zip(libs, cnvs, bams):
            norm_cmds, norm_beds = norm_tool_result(inbam2call[bam], F'{tmpdir}/scyn_output.csv', tool, lib, tobed)
            cmds2.extend(norm_cmds)
            bam2bed[bam] = norm_beds[1]
            lib2bed[lib] = norm_beds[1]
    if tool == 'chisel':
        cmd1 = F'''rm -r {tmpdir} || true && mkdir -p {tmpdir}/chisel_input/ {tmpdir}/chisel_output/'''
        cmd2 = F'''cp -s { ' '.join(bams + bais) } {tmpdir}/chisel_input/'''
        cmd3 = F'''cd {tmpdir}/chisel_output/ && conda run -n chisel chisel_prep {tmpdir}/chisel_input/*.bam'''
        cmd4 = F'''conda run -n chisel chisel_nonormal -t barcodedcells.bam -r {ref} -l {vcf} || true'''
        # chisel_pseudonormal mislabeled all the normal samples as tumor samples with -e 0.9
        # chisel_pseudonormal+chisel seems to run forever
        #cmd4 =(F'''conda run -n chisel chisel_pseudonormal -e 0.6 -n pseudonormal.bam -r {ref} barcodedcells.bam && '''
        #       F'''conda run -n chisel chisel -n pseudonormal.bam -t barcodedcells.bam -r {ref} -l {vcf} || true''')
        cmd5 = F'''cat {tmpdir}/chisel_output/barcodedcells.info.tsv {tmpdir}/chisel_output/calls/calls.tsv > {tmpdir}/chisel_output/barcodedcells.info.calls.tsv'''
        cmds.append(' && '.join([cmd1, cmd2, cmd3, cmd4, cmd5]))
        for lib, cnv, bam in zip(libs, cnvs, bams):
            norm_cmds, norm_beds = norm_tool_result(inbam2call[bam], F'{tmpdir}/chisel_output/barcodedcells.info.calls.tsv', tool, lib, tobed)
            cmds2.extend(norm_cmds)
            bam2bed[bam] = norm_beds[1]
            lib2bed[lib] = norm_beds[1]
    
    # revised by https://sorryios.ai/chat/1eb9e83b-aed2-4007-a65f-c9cd4e61b989
    if tool == 'alleloscope':
         alleloscope_R = F'{rootdir}/cnvguider/data3to4code/simplerun_alleloscope_v21.R'
         allelo_tumor_dir  = F'{tmpdir}/allelo_tumor/'
         allelo_normal_dir = F'{tmpdir}/allelo_normal/'
         allelo_out_dir    = F'{tmpdir}/allelo_out/'
         allelo_out_csv    = F'{tmpdir}/allelo_out.csv'
         # Prepare: symlink tumor BAMs into allelo_tumor/
         cmd_parts = [F'rm -r {tmpdir} || true',
                      F'mkdir -p {allelo_tumor_dir} {allelo_normal_dir} {allelo_out_dir}']
         for bam in bams:
             bai = change_file_ext(bam, 'bam.bai')
             cmd_parts.append(F'ln -sf {bam} {allelo_tumor_dir}/')
             cmd_parts.append(F'ln -sf {bai} {allelo_tumor_dir}/')
         # Prepare: symlink normal BAMs into allelo_normal/
         # If normal_bams_dir is provided (post-sim case), use those;
         # otherwise (pre-sim case), reuse the same BAMs as both tumor and normal.
         if normal_bams_dir:
             cmd_parts.append(F'ln -sf {normal_bams_dir}/*.bam {allelo_normal_dir}/')
             cmd_parts.append(F'ln -sf {normal_bams_dir}/*.bam.bai {allelo_normal_dir}/')
         else:
             for bam in bams:
                 bai = change_file_ext(bam, 'bam.bai')
                 cmd_parts.append(F'ln -sf {bam} {allelo_normal_dir}/')
                 cmd_parts.append(F'ln -sf {bai} {allelo_normal_dir}/')
         # Run Alleloscope
         vcf_clean = F'{allelo_out_csv}_clean.vcf.gz'
         cmd_parts.append(f'bcftools view -m2 -M2 -v snps -O z -o {vcf_clean} {vcf}')
         cmd_parts.append(f'bcftools index {vcf_clean}')
         cmd_parts.append(
             F'time -p conda run -n alleloscope Rscript {alleloscope_R}'
             F' {allelo_tumor_dir} {allelo_normal_dir} {ref} {vcf_clean}'
             F' {allelo_out_dir} {allelo_out_csv}'
         )
         cmds.append(' && '.join(cmd_parts) + F' #sequential=run.{tool}/')
         for bam, lib, cnv in zip(bams, libs, cnvs):
             norm_cmds, norm_beds = norm_tool_result(
                 inbam2call[bam], allelo_out_csv, tool, lib, tobed)
             cmds2.extend(norm_cmds)
             bam2bed[bam] = norm_beds[1]
             lib2bed[lib] = norm_beds[1]
    # revised by https://sorryios.ai/chat/e62de0ac-e570-4d62-9672-1ba233c8db90
    if tool == 'flcna':
        flcna_R = F'{rootdir}/cnvguider/data3to4code/simplerun_flcna.R'
        cmd = (F'rm -r {tmpdir}/flcna_input/ || true && mkdir -p {tmpdir}/flcna_input/ {tmpdir}/flcna_output/ && cp -s {" ".join(bams + bais)} {tmpdir}/flcna_input/ '
               F' && time -p Rscript {flcna_R} {tmpdir}/flcna_input/ {ref} {tmpdir}/flcna_output/flcna_result.csv {bigwig100} #sequential=run.{tool}/')
        cmds.append(cmd)
        for bam, lib, cnv in zip(bams, libs, cnvs):
            norm_cmds, norm_beds = norm_tool_result(inbam2call[bam], F'{tmpdir}/flcna_output/flcna_result.csv', tool, lib, tobed, ref=ref)
            cmds2.extend(norm_cmds)
            bam2bed[bam] = norm_beds[1]
            lib2bed[lib] = norm_beds[1]
    if tool == 'aneufinder':
        aneufinder_R = F'{rootdir}/cnvguider/data3to4code/simplerun_aneufinder.R'
        cmd = (F'rm -r {tmpdir}/aneufinder_input/ || true && mkdir -p {tmpdir}/aneufinder_input/ && cp -s {" ".join(bams + bais)} {tmpdir}/aneufinder_input/ '
               F' && time -p Rscript {aneufinder_R} {tmpdir}/aneufinder_input/ {tmpdir}/aneufinder_input/ {tmpdir}/aneufinder_output.bed #sequential=run.{tool}/')
        cmds.append(cmd)
        for bam, lib, cnv in zip(bams, libs, cnvs):
            norm_cmds, norm_beds = norm_tool_result(inbam2call[bam], F'{tmpdir}/aneufinder_output.bed', tool, lib, tobed)
            cmds2.extend(norm_cmds)
            bam2bed[bam] = norm_beds[1]
            lib2bed[lib] = norm_beds[1]

    if cmds and script in visited_scripts:
        logging.info(F'  Skip generating the script {script} because it has already been generated. ')
    elif cmds:
        with cm.myopen(script, writing_mode) as shfile:
            for cmd in cmds:
                write2file(cmd, shfile, script)
        visited_scripts.add(script)
        deps.append((script2, F'data4from2and3_1_run_DSA_{infodict["donor"]}_{infodict["sampleType"]}_{infodict["avgSpotLen"]}.rule'))
        deps.append((script2, F'data4from2and3_1_run_cellLine_{cellLine}.rule'))
        deps.append((script2, F'data4from2and3_1_run_tool_{tool}.rule'))
        deps.append((script2, F'data4from2and3_1_run_all.rule'))
    if cmds2 and script2 in visited_scripts:
        logging.info(F'  Skip generating the script {script2} because it has already been generated. ')
    elif cmds2:
        with cm.myopen(script2, writing_mode) as shfile:
            for cmd in cmds2:
                write2file(cmd, shfile, script2)
        visited_scripts.add(script2)   
        deps.append((script, script2))
        deps.append((script2, script_eval))
        deps.append((script2, F'data4from2and3_2_norm_DSA_{infodict["donor"]}_{infodict["sampleType"]}_{infodict["avgSpotLen"]}.rule'))
        deps.append((script2, F'data4from2and3_2_norm_cellLine_{cellLine}.rule'))
        deps.append((script2, F'data4from2and3_2_norm_tool_{tool}.rule'))
        deps.append((script2, F'data4from2and3_2_norm_all.rule'))
    return deps, cmds, bam2bed, lib2bed

def run_tool(infodict, df1, tool,
    rootdir, avgSpotLen, sampleType, donor, tool2script_dicts, writing_mode, visited_scripts):
    
    inst1into2end, inst2into3end, inst2from1vcf021, inst2into4logdir, inst2into4script, inst2into4scrip2, inst2into4tmpdir, inst4from2datdir, inst3into4logdir, inst3into4script, inst3into4scrip2, inst3into4tmpdir, inst4from3datdir, inst4into5logdir, inst4into5script, = find_replace_all([
    cm.t1into2end, cm.t2into3end, cm.t2from1vcf021, cm.t2into4logdir, cm.t2into4script, cm.t2into4scrip2, cm.t2into4tmpdir, cm.t4from2datdir, cm.t3into4logdir, cm.t3into4script, cm.t3into4scrip2, cm.t3into4tmpdir, cm.t4from3datdir, cm.t4into5logdir, cm.t4into5script], infodict)
    
    cm.makedirs((inst2into4logdir, inst2into4tmpdir, inst4from2datdir, inst3into4logdir, inst3into4tmpdir, inst4from3datdir, inst4into5logdir))
    post2pre_simbam = {}
    assert len(pd.unique(df1['Donor'])) == 1, F'The dataframe {df1} has multiple donors'
    donor = list(pd.unique(df1['Donor']))[0]
    inbam2call1, inbam2call2 = {}, {}
    for rowidx_1, (acc_1, lib_1, sample_1) in enumerate(zip(df1['#Run'], df1['Library~Name'], df1['Sample~Name'])):
        for GT in ['B', 'A']:
            infodict['GT'] = GT
            infodict['accession'] = acc_1
            infodict['samplename'] = None
            inst2from1mutbam1, inst2from1dedupb1, inst4from2depcns1, inst4from2intcns1 = find_replace_all([
            cm.t2from1mutbam , cm.t2from1dedupb , cm.t4from2depcns , cm.t4from2intcns ], infodict)
            
            infodict['samplename'] = bamfilename2samplename(inst2from1mutbam1)
            inst4from2depcns1, inst4from2intcns1 = find_replace_all([
            inst4from2depcns1, inst4from2intcns1], infodict)
            
            inbam2call1[inst2from1mutbam1] = ResultCN(input_bam=inst2from1mutbam1, dedup_bam=inst2from1dedupb1, simul_bed='', info_json='', depCN_bed=inst4from2depcns1, intCN_bed=inst4from2intcns1)
        
        for rowidx_2, (acc_2, lib_2, sample_2) in enumerate(zip(df1['#Run'], df1['Library~Name'], df1['Sample~Name'])):
            if not cm.circular_dist_below(rowidx_1, rowidx_2, len(df1)): continue
            
            infodict['GT'] = 'B'
            infodict['accession'] = acc_2
            infodict['samplename'] = None
            inst2from1mutbam2, inst4from2depcns2, inst4from2intcns2 = find_replace_all([
            cm.t2from1mutbam , cm.t4from2depcns , cm.t4from2intcns ], infodict)
            
            infodict['samplename'] = bamfilename2samplename(inst2from1mutbam2)
            inst4from2depcns2, inst4from2intcns2 = find_replace_all([
            inst4from2depcns2, inst4from2intcns2], infodict)
            
            infodict['accession_1'] = acc_1
            infodict['accession_2'] = acc_2
            infodict['samplename'] = None
            inst3from2simbam, inst3from2simbed, inst3from2dedupb, inst3from2infojs, inst4from3depcns, inst4from3intcns = find_replace_all([
            cm.t3from2simbam, cm.t3from2simbed, cm.t3from2dedupb, cm.t3from2infojs, cm.t4from3depcns, cm.t4from3intcns], infodict)
            
            infodict['samplename'] = bamfilename2samplename(inst3from2simbam)
            inst4from3depcns, inst4from3intcns = find_replace_all([
            inst4from3depcns, inst4from3intcns], infodict)
            
            inbam2call2[inst3from2simbam] = ResultCN(
                    input_bam=inst3from2simbam, 
                    dedup_bam=inst3from2dedupb, 
                    simul_bed=inst3from2simbed, 
                    info_json=inst3from2infojs, 
                    depCN_bed=inst4from3depcns, 
                    intCN_bed=inst4from3intcns)
            post2pre_simbam[inst3from2simbam] = (inst2from1mutbam1, inst2from1mutbam2)
           
    outeval_fname = inst4into5script # os.path.dirname(data3prefix1) + F'to4.step{tool_order}_{donor}_{sampleType}_{avgSpotLen}_{tool}_eval.sh'

    # For alleloscope, pre-sim BAMs serve as matched normals for post-sim runs
    normal_bams_dir_for_postsim = None
    if tool == 'alleloscope':
        normal_bams_dir_for_postsim = inst2into4tmpdir + '/allelo_tumor'

    deps1, cmds1, bam2bed1, lib2bed1 = run_tool_1(infodict, tool, inbam2call1, inst2into4tmpdir, inst2into4script, inst2into4scrip2, 
            outeval_fname, rootdir, inst2from1vcf021, tool2script_dicts[0], inst1into2end, 
            is_overall_haploid=True , writing_mode=writing_mode, visited_scripts=visited_scripts)
    tool2script_dicts[0][tool] = inst2into4script
    deps2, cmds2, bam2bed2, lib2bed2 = run_tool_1(infodict, tool, inbam2call2, inst3into4tmpdir, inst3into4script, inst3into4scrip2,
            outeval_fname, rootdir, inst2from1vcf021, tool2script_dicts[1], inst2into3end, 
            is_overall_haploid=False, writing_mode=writing_mode, visited_scripts=set([]), normal_bams_dir=normal_bams_dir_for_postsim)

    tool2script_dicts[1][tool] = inst3into4script
    deps3, cmds3 = [], []
    if tool in SC_CN_EVAL_TOOLS:
        with cm.myopen(outeval_fname, writing_mode) as outeval_file:
            for postsim, (presim1, presim2) in post2pre_simbam.items():
                presim_bed1          = inbam2call1[presim1].intCN_bed + '.overall_haploid'
                presim_bed2          = inbam2call1[presim2].intCN_bed + '.overall_haploid'
                postsim_int_bed      = inbam2call2[postsim].intCN_bed
                postsim_dep_bed      = inbam2call2[postsim].depCN_bed
                sim_approx_truth_bed = inbam2call2[postsim].simul_bed
                cmd3 = F'python {rootdir}/cnvguider/cnv_bedset_to_consistency.py --pre-sim-call-bed-1 {presim_bed1} --pre-sim-call-bed-2 {presim_bed2} --approx-truth-bed {sim_approx_truth_bed} --post-sim-call-bed {postsim_int_bed} --post-sim-call-bed-by-DP {postsim_dep_bed} '
                cmds3.append(cmd3)
                write2file(cmd3, outeval_file, outeval_fname)
        deps3.append((inst2into4scrip2, outeval_fname))
        deps3.append((inst3into4scrip2, outeval_fname))
        deps3.append((outeval_fname, F'data4from2and3_3_eval_DSA_{infodict["donor"]}_{infodict["sampleType"]}_{infodict["avgSpotLen"]}.rule'))
        deps3.append((outeval_fname, F'data4from2and3_3_eval_cellLine_{infodict["cellLine"]}.rule'))
        deps3.append((outeval_fname, F'data4from2and3_3_eval_tool_{tool}.rule'))
        deps3.append((outeval_fname, F'data4from2and3_3_eval_all.rule'))
    return deps1 + deps2 + deps3

def main(args1=None):
    ret = []
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(pathname)s:%(lineno)d %(levelname)s - %(message)s')

    script_dir = os.path.dirname(os.path.abspath(__file__))
    datadir = os.path.abspath(os.path.sep.join([script_dir, '..', 'data']))
    data0to1dir, data1to2dir , data2to3dir, data2to4dir, data3to4dir, data4to5dir = cm.get_varnames(datadir)
    root = os.path.abspath(os.path.sep.join([script_dir, '..']))

    defaultSraRunTable = os.path.sep.join([script_dir, 'scDNAaccessions.tsv'])
    defaultSraRunTable = os.path.sep.join([script_dir, 'scDNAaccessions.tsv'])
    cosmic_cn_filename = os.path.sep.join([script_dir, 'cosmic-v97', 'cell_lines_copy_number.csv'])
    cosmic_cell_lines = ['COLO-829', 'HCC1395', 'HeLa']

    parser = argparse.ArgumentParser(description='Generate bash commands to in-silico mix the single-cell sequencing data with copy numbers simulated from the COSMIC cell line database. ',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--SraRunTable',type=str,default=defaultSraRunTable, help=(
        'SraRunTable in TSV format containing the columns '
        '#Run, AvgSpotLen, Library~Name, Sample~Name, sample-type, Oocyte_ID, Donor, and SRA~Study'))
    parser.add_argument('--tools', nargs='+', default=SC_CN_TOOLS, choices=SC_CN_TOOLS,
        help='Software tools calling cell-specific copy numbers from from single-cell DNA-seq data')
    parser.add_argument('--cell-lines',nargs='+',default=cosmic_cell_lines, help='Cell-lines')
    parser.add_argument('-w', '--writing-mode', type=str, default=cm.DEFAULT_WRITING_MODE,
        help='File open mode for writing commands to shell script, pass any of {cm.OVERWRITING_PREVENTION_MODES} to prevent overwriting existing scripts (or w to do not prevent such thing). ')
    
    args = (args1 if args1 else parser.parse_args())
    
    visited_scripts = set([])
    df0 = pd.read_csv(args.SraRunTable, sep='\t', header=0)
    df0['sample-type'] = norm_sample_type(df0)
    grouped = df0.groupby(['AvgSpotLen', 'sample-type', 'Donor'])
    partitioned_dfs = {partkey: df1 for partkey, df1 in grouped}
    for cell_line in args.cell_lines:
        for (avgSpotLen, sampleType, donor), df1 in sorted(partitioned_dfs.items()):
            tool2script_dicts = {}, {}
            for tool in sorted(args.tools, key=lambda x:(-SC_CN_TOOL_TO_RUN_ORDER[x],x)):
                logging.info(F'cell_line={cell_line} avgSpotLen={avgSpotLen} sampleType={sampleType} donor={donor} tool={tool}')
                infodict = {
                        'data0to1dir': data0to1dir,
                        'data1to2dir': data1to2dir,
                        'data2to3dir': data2to3dir,
                        'data2to4dir': data2to4dir,
                        'data3to4dir': data3to4dir,
                        'data4to5dir': data4to5dir,
                        'donor'      : str(donor),
                        'sampleType' : str(sampleType),
                        'avgSpotLen' : str(avgSpotLen),
                        'cellLine'   : str(cell_line),
                        'tool'       : str(tool),
                        'tool_order' : str(SC_CN_TOOL_TO_RUN_ORDER[tool]),
                        'tool_ord_1' : str(SC_CN_TOOL_TO_RUN_ORDER[tool]+1),
                        'tool_ord_2' : str(SC_CN_TOOL_TO_RUN_ORDER[tool]+2),
                }
                deps = run_tool(infodict, df1, tool, root, avgSpotLen, sampleType, donor, tool2script_dicts, args.writing_mode, visited_scripts)
                ret.extend(deps)
    return ret
if __name__ == '__main__': print(cm.list2snakemake(main()))

