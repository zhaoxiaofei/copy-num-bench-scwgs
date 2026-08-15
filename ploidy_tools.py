#!/usr/bin/env python
# revised by: https://claude.ai/chat/1154badf-2241-4cce-9e52-7c5565f9abe3
# revised by: https://claude.ai/chat/3874e270-ba5a-433d-ab66-cefaff6fb550

"""
Benchmark ploidy inference, for every tool in this repository: scAbsolute, whose primary output
is a ploidy estimate rather than a per-cell copy-number profile, and each of the CNV callers of
data4from2and3.py, whose ploidy is implied by the copy numbers it calls. The reference is
    Schneider MP, Cullen AE, Pangonyte J, et al.
    "scAbsolute: measuring single-cell ploidy and replication status."
    Genome Biology 2024;25:62.  https://doi.org/10.1186/s13059-024-03204-y
Why a separate module
---------------------
Every caller benchmarked by data4from2and3.py emits per-cell CN segments, so its ploidy is
*derived*: ploidy_eval.py reads the caller's `*intcns.bed` and takes the length-weighted mean
copy number (scAbsolute Eq. 1). A ploidy-inference tool instead *reports* ploidy directly, so
there is no BED to read and the existing evaluation step cannot be pointed at it. This module
supplies the missing piece -- running the tool, normalising its output, and scoring it with the
very same metrics -- and nothing else.
A CNV caller named in --ploidy-tools takes the same path, with one difference: it has already
been run by run_tool_1, so there is nothing to run and `eval` simply reads its ploidy out of the
`*intcns.bed` files it writes, with the same bed_to_ploidy() that ploidy_eval.py uses. The point
is not a second derivation -- it is that the caller then gets the whole report a ploidy tool
gets, including the per-sample point estimate of `consensus_ploidy` that ploidy_eval.py does not
compute, so that scAbsolute and the callers are ranked on identical numbers.
Where the ground truth comes from
--------------------------------
  * Real tumors (`--tumor-fastq`): the experimental FACS/DAPI ploidy of `--ploidy-file`, one
    number per sample.
  * Simulated data (no `--tumor-fastq`): the simulator's own per-cell CN profile, the
    `*_simtruth.bed` that data3from2.py writes next to each simulated BAM. The expected ploidy
    is then known exactly, per cell rather than per sample, and is read with the same
    bed_to_ploidy() that turns a caller's BED into an observed ploidy -- so truth and estimate
    share the `--chroms` and `--max-cn` conventions and the comparison is like for like.
    Only the post-simulation cells are scored; the pre-simulation germline cells are the raw
    material of the simulation and the pipeline normalises their calls to an overall-haploid
    scale, so no absolute ploidy is defined for them. The ploidy tool still runs on them, as
    every caller in this repository does, and its calls are written out for inspection.
In the simulated mode `--ploidy-tools` also switches on a ploidy evaluation of the ordinary
CNV callers against the same simulated truth, which is what makes the benchmark a comparison:
scAbsolute, which reports ploidy, against nine callers whose ploidy is implied by their
segments -- on data whose answer is known by construction. Naming those callers in
`--ploidy-tools` as well adds the per-sample point estimate for them too.
Strict backward compatibility
-----------------------------
The feature is inert until `setup()` is called, which main.py does only when `--ploidy-tools`
is passed. Concretely:
  * no existing function is edited; the tool-specific work is done in a wrapper around
    `data4from2and3.run_tool_1` installed by `setup()`, exactly as gink_custom_binning.py does,
    and the wrapper delegates untouched for every tool it does not own;
  * the metrics and the plot come from ploidy_eval.py by *import*, never by modification;
  * a ploidy tool is deliberately kept out of `SC_CN_EVAL_TOOLS`, so the clustermap and
    ploidy-evaluation branches of data_tumor.py stay switched off for it;
  * a CNV caller named in --ploidy-tools is registered nowhere at all: its run order, run mode,
    dependencies and `SC_CN_EVAL_TOOLS` membership stay exactly as they are, and it gains one
    evaluation script and no run rule, so the CNV benchmark is bit-for-bit unaffected;
  * new results are written to new files (`*_ploidy_calls.tsv`, `*_ploidy_tool_eval_*`) under
    the existing directory layout; no default output file is read, moved or overwritten.
Running the pipeline without `--ploidy-tools` therefore reproduces the previous Snakefile
verbatim (rule for rule, byte for byte, apart from the git-commit banner).
What is measured
----------------
1. Ploidy accuracy, per cell, with the scAbsolute metrics already implemented in
   ploidy_eval.py: the percentage of cells outside the +/- `--ploidy-window` experimental
   window, the mean absolute ploidy distance, and the 2x / 0.5x scaling-error diagnostics.
   The `_percell.tsv` / `_persample.tsv` / `_summary.json` triple has the same columns as the
   CNV-caller output of ploidy_eval.py, so the two can simply be concatenated for a joint
   ranking of "tools that infer ploidy" against "callers that imply ploidy".
2. Ploidy accuracy of the single per-sample point estimate that a ploidy-inference tool is
   usually asked for. It is derived from the per-cell values the way scAbsolute's own
   scripts/estimatePloidy.R does it (see `consensus_ploidy`) and reported in extra
   `sample_*` columns, appended to -- not substituted into -- the per-sample table.
3. Optionally (`--ploidy-facs`), the downstream effect of feeding the inferred ploidy into a
   CN caller: the calls are converted to a Ginkgo FACS file and a second Ginkgo pass named
   `ginkgo_facs_<tool>` is generated. That pass is a normal evaluation tool, so it picks up
   the existing clustermap and ploidy-evaluation steps for free and can be compared with the
   untouched `ginkgo` run in the usual result tables.
Command-line use (each subcommand also runs standalone, outside snakemake)
-------------------------------------------------------------------------
    # score a ploidy-calls TSV against the experimental ploidy
    python ploidy_tools.py eval -i '<...>_ploidy_calls.tsv' -o <out_prefix> \
        --ploidy-file ploidy.PRJNA629885.tsv --metadata-tsv SraRunTable.tsv --plot
    # turn the same TSV into a Ginkgo --facs file
    python ploidy_tools.py facs -i '<...>_ploidy_calls.tsv' -o <out.txt>
A ploidy-calls TSV is the canonical hand-off format between a tool and this module: a header
line plus one row per cell, with the required columns `cell` and `ploidy` and any number of
tool-specific extras (scAbsolute adds `rpc`, `used_reads` and `failure_reason`), which are
carried through to `_percell.tsv` untouched.
"""
import argparse, glob, json, logging, os, sys

import numpy as np
import pandas as pd
import common as cm
from common import change_file_ext, find_replace_all, write2file
import data_tumor
import ploidy_eval as pe
# ---------------------------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------------------------
TOOLS = ['scabsolute'] + sorted(data_tumor.SC_CN_EVAL_TOOLS)
# A CNV caller infers a ploidy too, it just does not report one: its ploidy is the length-weighted
# mean of the integer copy numbers it called (scAbsolute Eq. 1). Naming a caller in --ploidy-tools
# scores that ploidy with the very same metrics, and adds the per-sample point estimate that
# ploidy_eval.py does not compute, so that every tool in the benchmark is measured the same way.
# It costs no extra run: `eval` reads the ploidy out of the *intcns.bed files the caller has
# already written, which is why a caller is registered nowhere below and only gains a script.
CALLER_TOOLS = TOOLS[1:]
FACS_PREFIX = 'ginkgo_facs_'
# A ploidy tool runs before the CN callers (whose order is 2) so that the optional
# ginkgo_facs_<tool> pass, which consumes its output, is generated after it.
RUN_ORDER = 3
# scAbsolute's working point: its scaling step is only defined for bins of at most 1 Mb, and
# 500 kb is the bin size used throughout the paper.
DEFAULT_BIN_SIZE = 500
DEFAULT_GENOME = 'hg19'
PLOIDY_TOOLS_HELP = (
    'Ploidy-inference tools to benchmark in addition to the CNV callers. These report a ploidy '
    'estimate instead of a per-cell CN profile, so they are scored by ploidy_tools.py against '
    'the same --ploidy-file (which is therefore required) with the scAbsolute metrics. Leaving '
    'this empty keeps the pipeline exactly as it is.')
PLOIDY_FACS_HELP = (
    'Additionally re-run Ginkgo with each inferred ploidy supplied through its --facs option '
    '(as tool ginkgo_facs_<ploidy tool>), to measure how the ploidy estimate propagates into '
    'the downstream integer copy-number states.')
# Same naming scheme as the CNV callers, so the new files land beside the existing ones.
t_calls = ('<data2to4dir>/<donor>/4from2_2_<donor>_3_<sampleType>_<avgSpotLen>_4_step<tool_order>_<tool>.datdir/'
           '2_<donor>_3_<sampleType>_<avgSpotLen>_4_step<tool_order>_<tool>_4from2_ploidy_calls.tsv')
t_facs = ('<data2to4dir>/<donor>/4from2_2_<donor>_3_<sampleType>_<avgSpotLen>_4_step<tool_order>_<tool>.datdir/'
          '2_<donor>_3_<sampleType>_<avgSpotLen>_4_step<tool_order>_<tool>_4from2_ginkgo_facs.txt')
t_prefix = ('<data2to4dir>/<donor>/4from2_2_<donor>_3_<sampleType>_<avgSpotLen>_4_step<tool_order>_<tool>'
            '_ploidy_tool_eval')
# Post-simulation counterparts, derived from the canonical templates so that they cannot drift
# away from where data4from2and3.py puts the simulated cells and their per-cell truth BEDs.
g_datdir = cm.t4from3datdir                                             # .../4from3_..._<tool>.datdir/
g_calls = g_datdir + '2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>_4_step<tool_order>_<tool>_4from3_ploidy_calls.tsv'
g_facs = g_datdir + '2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>_4_step<tool_order>_<tool>_4from3_ginkgo_facs.txt'
g_prefix_tool = g_datdir[:-len('.datdir/')] + '_ploidy_tool_eval'       # ploidy-inference tools
g_prefix_call = g_datdir[:-len('.datdir/')] + '_ploidy_eval'            # CNV callers, as in data_tumor.py
g_script_call = cm.t3into4script.replace('_3into4_call.sh', '_3into4_ploidy_eval.sh')
g_script_tool = cm.t3into4script.replace('_3into4_call.sh', '_3into4_ploidy_tool_eval.sh')
g_truth_glob = os.path.dirname(cm.t3from2simbed) + '/*' + pe.TRUTH_BED_SUFFIX
# Both stay empty unless setup() runs, which is what makes the feature inert by default.
_CFG = {}
_CALLS = {}   # (donor, sampleType, avgSpotLen) -> dict(calls=, facs=, script=, tool=)
# Positional parameters of data4from2and3.run_tool_1 after (infodict, tool), so that the
# wrapper can read them by name however the caller chose to pass them.
_POS = ['inbam2call', 'tmpdir', 'script', 'script2', 'script_eval', 'rootdir', 'vcf',
        'tool2script_dict', 'start_script', 'is_overall_haploid', 'writing_mode',
        'visited_scripts', 'normal_bams_dir']

def setup(d4, tools, ploidy_tools, args):
    """Register the requested ploidy-inference tools and return the extended tool list.
    Returns `tools` unchanged, and leaves `d4` untouched, when no ploidy tool is requested.
    """
    ploidy_tools = list(dict.fromkeys(ploidy_tools or []))
    tools = list(tools)
    if not ploidy_tools: return tools
    unknown = [tool for tool in ploidy_tools if tool not in TOOLS]
    if unknown: raise ValueError(F'Unknown ploidy-inference tool(s) {unknown}; known: {TOOLS}')
    # Real tumors have no simulated truth, so there the experimental ploidy file is mandatory;
    # simulated data carries its own per-cell truth BEDs and needs nothing extra.
    simulated = not getattr(args, 'tumor_fastq', False)
    if not simulated and not getattr(args, 'ploidy_file', None):
        raise ValueError('--ploidy-tools requires --ploidy-file when used with --tumor-fastq')
    facs = bool(getattr(args, 'ploidy_facs', False))
    _CFG.update(tools=ploidy_tools, facs=facs, simulated=simulated, d4=d4,
                ploidy_file=(os.path.abspath(args.ploidy_file) if getattr(args, 'ploidy_file', None) else ''),
                window=getattr(args, 'ploidy_window', pe.DEFAULT_PLOIDY_WINDOW),
                chroms=getattr(args, 'ploidy_chroms', 'autosomes'),
                metadata_tsv=(os.path.abspath(args.SraRunTable) if getattr(args, 'SraRunTable', None) else ''),
                bin_size=DEFAULT_BIN_SIZE, genome=DEFAULT_GENOME)
    for tool in ploidy_tools:
        if tool in CALLER_TOOLS:
            # Registered already, as the caller it is. Its run order, run mode, dependencies and
            # SC_CN_EVAL_TOOLS membership all stay exactly as they are, so that the CNV benchmark
            # is bit-for-bit unaffected; only an evaluation of its ploidy is added.
            if tool not in tools:
                raise ValueError(F'--ploidy-tools {tool} scores the ploidy of the CNV caller '
                                 F'{tool}, so {tool} has to be in --tools as well')
            continue
        if tool not in tools: tools.append(tool)
        # 'nop' is the sentinel that waits for the alignments, like it does for every caller.
        d4.SC_CN_TOOL_DEPENDENCY_TO_DEPENDENT['nop'][tool] = ''
        d4.SC_CN_TOOL_TO_RUN_ORDER[tool] = RUN_ORDER
        d4.SC_CN_TOOL_TO_RUN_MODE[tool] = 'sequential'
        # NOT added to d4.SC_CN_EVAL_TOOLS: no per-cell CN BED is produced, so the clustermap
        # and ploidy_eval branches of data_tumor.py must not fire for this tool.
        if not facs: continue
        ftool = FACS_PREFIX + tool
        if ftool not in tools: tools.append(ftool)
        d4.SC_CN_TOOL_DEPENDENCY_TO_DEPENDENT['bam2bed'][ftool] = ''
        d4.SC_CN_TOOL_TO_RUN_ORDER[ftool] = d4.SC_CN_TOOL_TO_RUN_ORDER['ginkgo']
        d4.SC_CN_TOOL_TO_RUN_MODE[ftool] = d4.SC_CN_TOOL_TO_RUN_MODE['ginkgo']
        d4.SC_CN_EVAL_TOOLS.add(ftool)   # it *is* a CN caller, so it is evaluated like one
    if set(ploidy_tools) - set(CALLER_TOOLS) and 'nop' not in tools: tools.append('nop')
    if facs and 'bam2bed' not in tools: tools.append('bam2bed')
    if not getattr(d4.run_tool_1, '_ploidy_tools', False): d4.run_tool_1 = _wrap(d4.run_tool_1)
    return tools

def _wrap(original):
    def run_tool_1(infodict, tool, *args, **kwargs):
        params = dict(zip(_POS, args)); params.update(kwargs)
        if tool in _CFG.get('tools', ()) and tool not in CALLER_TOOLS:
            return _gen_ploidy_tool(infodict, tool, params)
        if tool.startswith(FACS_PREFIX) and _CFG.get('facs'):
            return _gen_facs_ginkgo(original, infodict, tool, args, kwargs, params)
        out = original(infodict, tool, *args, **kwargs)
        deps = list(out[0])
        # A caller named in --ploidy-tools has just been run by the untouched run_tool_1 above, so
        # all that is added is an evaluation of the ploidy implied by the BEDs it writes.
        if tool in _CFG.get('tools', ()):
            deps += _gen_ploidy_tool(infodict, tool, params)[0]
        # On simulated data every CNV caller is scored on ploidy too, against the same truth,
        # which is what turns the ploidy tool's numbers into a comparison. Post-simulation
        # cells only: the pre-simulation calls are normalised to an overall-haploid scale.
        if (_CFG.get('simulated') and not params.get('is_overall_haploid', True)
                and tool in _CFG['d4'].SC_CN_EVAL_TOOLS):
            deps += _gen_caller_ploidy_eval(infodict, tool, params)
        return (deps,) + tuple(out[1:])
    run_tool_1._ploidy_tools = True
    return run_tool_1

def _sim_eval(infodict, params, script_tpl, prefix_tpl, cmd_fn, rule_stem, src_script):
    """One evaluation script per --max-cn setting, scored against the simulated per-cell truth.
    The cap changes what a ploidy means rather than just its precision (see ploidy_eval.py), and
    it is applied to the truth as well as to the estimate, so both settings have to be evaluated
    the same way here as data_tumor.py does for real tumors.
    """
    deps = []
    truth_glob, = find_replace_all([g_truth_glob], infodict)
    for max_cn in data_tumor.PLOIDY_MAX_CNS:
        sfx = '' if float(max_cn) == float(pe.DEFAULT_MAX_CN) else F'_maxcn_{max_cn}'
        script, prefix = find_replace_all([script_tpl, prefix_tpl], infodict)
        script, prefix = script.removesuffix('.sh') + sfx + '.sh', prefix + sfx
        cm.makedirs((script, prefix))
        with cm.myopen(script, params['writing_mode']) as fh:
            write2file(cmd_fn(truth_glob, prefix, max_cn, sfx), fh, script)
        deps.append((src_script, script))
        for rule in (F'data4from2and3_{rule_stem}_DSA_{infodict["donor"]}_{infodict["sampleType"]}_{infodict["avgSpotLen"]}.rule',
                     F'data4from2and3_{rule_stem}_cellLine_{infodict["cellLine"]}.rule',
                     F'data4from2and3_{rule_stem}_tool_{infodict["tool"]}.rule',
                     F'data4from2and3_{rule_stem}_all.rule'):
            deps.append((script, rule))
    return deps

def _gen_caller_ploidy_eval(infodict, tool, params):
    """ploidy_eval.py on one caller's simulated per-cell BEDs, in the layout of data_tumor.py."""
    datdir, = find_replace_all([cm.t4from3datdir], infodict)
    label = infodict['cellLine']
    def cmd_fn(truth_glob, prefix, max_cn, sfx):
        title = (F'{tool} | cellLine={label} donor={infodict["donor"]} '
                 F'sampleType={infodict["sampleType"]} avgSpotLen={infodict["avgSpotLen"]}'
                 + (F' | max-cn={max_cn}' if sfx else ''))
        return (F'python {params["rootdir"]}/copy-num-bench-scwgs/ploidy_eval.py '
                F'-i {datdir}*intcns.bed -o {prefix} --truth-bed "{truth_glob}" '
                F'--sample-label {label} --ploidy-window {_CFG["window"]} '
                F'--chroms {_CFG["chroms"]} --max-cn {max_cn} '
                F'--tool {tool} --donor {infodict["donor"]} --sample-type {infodict["sampleType"]} '
                F'--avg-spot-len {infodict["avgSpotLen"]} --plot --title "{title}" '
                F'#sequential=ploidy_eval.{tool}{sfx}/')
    return _sim_eval(infodict, params, g_script_call, g_prefix_call, cmd_fn,
                     '4_ploidy_eval', params['script2'])

def _run_cmd(tool, rootdir, indir, calls):
    """Shell command running one ploidy-inference tool over a directory of BAM files."""
    if tool == 'scabsolute':
        rscript = F'{rootdir}/copy-num-bench-scwgs/data3to4code/simplerun_scabsolute.R'
        basedir = os.getenv('scAbsoluteRoot', F'{rootdir}/copy-num-bench-scwgs/data3to4code/scAbsolute')
        return (F'time -p conda run -n scabsolute Rscript {rscript} '
                F'{indir} {calls} {basedir} {_CFG["bin_size"]} {_CFG["genome"]}')
    raise ValueError(F'No run command is defined for the ploidy-inference tool {tool}')

def _gen_ploidy_tool(infodict, tool, params):
    """Run script + evaluation script for one ploidy-inference tool, in the layout of run_tool_1.
    Mirrors the (deps, cmds, bam2bed, lib2bed) contract of data4from2and3.run_tool_1; the two
    BED dictionaries stay empty because a ploidy tool produces no CN segments.
    """
    script, script2, tmpdir = params['script'], params['script2'], params['tmpdir']
    rootdir, wmode, visited = params['rootdir'], params['writing_mode'], params['visited_scripts']
    donor, sampleType, avgSpotLen = infodict['donor'], infodict['sampleType'], infodict['avgSpotLen']
    # Simulated data has two runs per tool: the pre-simulation germline cells (whose calls the
    # pipeline normalises to an overall-haploid scale, so no absolute ploidy is defined for
    # them) and the post-simulation cells, which carry a per-cell truth BED.
    presim = bool(_CFG.get('simulated') and params.get('is_overall_haploid', False))
    postsim = bool(_CFG.get('simulated') and not presim)
    calls, facs, prefix = find_replace_all(
        [g_calls, g_facs, g_prefix_tool] if postsim else [t_calls, t_facs, t_prefix], infodict)
    caller = tool in CALLER_TOOLS
    if caller:
        # The caller was run by run_tool_1 itself, and `eval` reads its ploidy straight out of the
        # per-cell BEDs it writes, so there is no run command, no calls file and no run rule here.
        # Its evaluation script sits beside its normalisation script, which is already taken.
        datdir, = find_replace_all([cm.t4from3datdir if postsim else cm.t4from2datdir], infodict)
        calls, script, script2 = (F'"{datdir}*intcns.bed"', script2,
                                  script2.replace('_norm.sh', '_ploidy_tool_eval.sh'))
        if presim: return [], [], {}, {}
        cmds, deps = [], []
    else:
        cm.makedirs((calls,))
        _CALLS[(donor, sampleType, avgSpotLen, infodict.get('cellLine', ''), presim)] = dict(
            calls=calls, facs=facs, script=script, tool=tool)
        bams = sorted(params['inbam2call'].keys())
        bais = [change_file_ext(bam, 'bam.bai') for bam in bams]
        indir = F'{tmpdir}/{tool}_input'
        cmd = (F'rm -r {indir} || true && mkdir -p {indir} && cp -s {" ".join(bams + bais)} {indir}/ '
               F'&& {_run_cmd(tool, rootdir, indir, calls)}')
        if _CFG.get('facs'):
            cmd += F' && python {rootdir}/copy-num-bench-scwgs/ploidy_tools.py facs -i {calls} -o {facs}'
        cmd += F' #sequential=run.{tool}/'
        cmds, deps = [cmd], []
    scripts = [(script, cmd)] if not caller else []
    if not _CFG.get('simulated'):
        metadata_arg = (F'--metadata-tsv "{_CFG["metadata_tsv"]}" ' if _CFG.get('metadata_tsv') else '')
        title = F'{tool} | donor={donor} sampleType={sampleType} avgSpotLen={avgSpotLen}'
        scripts.append((script2,
            F'python {rootdir}/copy-num-bench-scwgs/ploidy_tools.py eval '
            F'-i {calls} -o {prefix} --ploidy-file "{_CFG["ploidy_file"]}" {metadata_arg}'
            F'--ploidy-window {_CFG["window"]} --chroms {_CFG["chroms"]} '
            F'--tool {tool} --donor {donor} --sample-type {sampleType} '
            F'--avg-spot-len {avgSpotLen} --plot --title "{title}" #sequential=ploidy_tool_eval.{tool}/'))
    for fname, text in scripts:
        if fname in visited:
            logging.info(F'  Skip generating the script {fname} because it has already been generated. ')
            continue
        with cm.myopen(fname, wmode) as file: write2file(text, file, fname)
        visited.add(fname)
    for rule in ([] if caller else
                 (F'data4from2and3_1_run_DSA_{donor}_{sampleType}_{avgSpotLen}.rule',
                  F'data4from2and3_1_run_tool_{tool}.rule', 'data4from2and3_1_run_all.rule')):
        deps.append((script, rule))
    if postsim:
        # Scored exactly like the CNV callers above, against the same per-cell truth BEDs.
        def cmd_fn(truth_glob, prefix_1, max_cn, sfx):
            title = (F'{tool} | cellLine={infodict["cellLine"]} donor={donor} '
                     F'sampleType={sampleType} avgSpotLen={avgSpotLen}'
                     + (F' | max-cn={max_cn}' if sfx else ''))
            return (F'python {rootdir}/copy-num-bench-scwgs/ploidy_tools.py eval '
                    F'-i {calls} -o {prefix_1} --truth-bed "{truth_glob}" '
                    F'--sample-label {infodict["cellLine"]} --ploidy-window {_CFG["window"]} '
                    F'--chroms {_CFG["chroms"]} --max-cn {max_cn} '
                    F'--tool {tool} --donor {donor} --sample-type {sampleType} '
                    F'--avg-spot-len {avgSpotLen} --plot --title "{title}" '
                    F'#sequential=ploidy_tool_eval.{tool}{sfx}/')
        deps.extend(_sim_eval(infodict, params, g_script_tool, g_prefix_tool, cmd_fn,
                              '5_ploidy_tool_eval', script))
    elif not _CFG.get('simulated'):
        deps.append((script, script2))
        for rule in (F'data4from2and3_5_ploidy_tool_eval_DSA_{donor}_{sampleType}_{avgSpotLen}.rule',
                     F'data4from2and3_5_ploidy_tool_eval_tool_{tool}.rule',
                     'data4from2and3_5_ploidy_tool_eval_all.rule'):
            deps.append((script2, rule))
    return deps, cmds, {}, {}

def _gen_facs_ginkgo(original, infodict, tool, args, kwargs, params):
    """A second Ginkgo pass whose per-cell ploidy is fixed to the inferred one, via --facs.
    Ginkgo is run through the unmodified run_tool_1, and only the two anchors that the FACS
    file has to reach are rewritten afterwards, the same trick gink_custom_binning.py uses.
    """
    key = (infodict['donor'], infodict['sampleType'], infodict['avgSpotLen'],
           infodict.get('cellLine', ''),
           bool(_CFG.get('simulated') and params.get('is_overall_haploid', False)))
    called = _CALLS.get(key)
    if not called: raise ValueError(F'{tool}: no ploidy calls were generated for {key}')
    deps, cmds, bam2bed, lib2bed = original(infodict, 'ginkgo', *args, **kwargs)
    tmpdir, script = params['tmpdir'], params['script']
    # Ginkgo edits its FACS file in place (it strips '.bed' and re-sorts it), so it gets a copy.
    local = F'{tmpdir}/ginkgo_facs.txt'
    patches = [(F'rm {tmpdir}/*.bed || true', F'cp {called["facs"]} {local} && (rm {tmpdir}/*.bed || true)'),
               ('--genome hg19 --binning', F'--genome hg19 --facs {local} --binning'),
               ('#sequential=run.ginkgo/', F'#sequential=run.{tool}/')]
    if os.path.exists(script):
        text = open(script).read()
        for old, new in patches: text = text.replace(old, new)
        open(script, 'w').write(text)
    for old, new in patches: cmds = [cmd.replace(old, new) for cmd in cmds]
    deps = [tuple((d.replace('_tool_ginkgo.rule', F'_tool_{tool}.rule') if isinstance(d, str) else d)
                  for d in dep) for dep in deps]
    deps.append((called['script'], script))   # the ploidy calls must exist before Ginkgo starts
    return deps, cmds, bam2bed, lib2bed
# ---------------------------------------------------------------------------------------------
# Ploidy calls
# ---------------------------------------------------------------------------------------------

def read_calls(patterns, chroms=None, max_cn=pe.DEFAULT_MAX_CN):
    """Read one or more ploidy-calls TSVs (columns `cell` and `ploidy`, plus any extras).
    A `*.bed` input is a CNV caller's per-cell integer copy numbers instead. A caller states no
    ploidy, so its ploidy is the length-weighted mean of those copy numbers (scAbsolute Eq. 1),
    read here with the very same bed_to_ploidy(), `chroms` and `max_cn` that ploidy_eval.py uses
    -- so a caller and a ploidy-inference tool end up being scored on the same quantity.
    """
    files = []
    for pat in patterns:
        matched = sorted(glob.glob(pat)) if any(c in pat for c in '*?[') else [pat]
        for path in matched:
            logging.info('ploidy call matched path: %s', os.path.abspath(path))
        files.extend(matched)
    frames = []
    for path in files:
        if not (os.path.isfile(path) and os.path.getsize(path)): continue
        if path.endswith('.bed'):
            # The full path is the cell name, which is what pe.cell_candidates() and
            # pe.truth_lookup() already know how to resolve back to a sample.
            frames.append(pd.DataFrame([{'cell': path, 'ploidy': pe.bed_to_ploidy(
                path, chroms=chroms, weight='length', max_cn=max_cn)[0]}]))
            continue
        df = pd.read_csv(path, sep='\t', comment='#')
        missing = [c for c in ('cell', 'ploidy') if c not in df.columns]
        if missing: raise ValueError(F'{path}: ploidy-calls file lacks the column(s) {missing}')
        frames.append(df.assign(calls_file=os.path.abspath(path)))
    if not frames: raise ValueError(F'No usable ploidy-calls file among {patterns}')
    out = pd.concat(frames, ignore_index=True)
    out['ploidy'] = pd.to_numeric(out['ploidy'], errors='coerce')
    return out

def consensus_ploidy(values):
    """The single per-sample ploidy of a set of per-cell estimates.

    Replicates scAbsolute's own scripts/estimatePloidy.R, which rounds the median over the
    cells of a sample and floors the result at 2 (their sample-level estimate is an integer).
    """
    finite = np.asarray([v for v in np.asarray(values, dtype=float) if np.isfinite(v)])
    if not len(finite): return float('nan')
    return float(max(2.0, round(float(np.median(finite)))))

def _facs_main(args):
    """Ginkgo FACS file: `<cell> <TAB> <ploidy>`, headerless, non-positive ploidies dropped."""
    calls = read_calls(args.input)
    calls = calls[np.isfinite(calls['ploidy']) & (calls['ploidy'] > 0)]
    if calls.empty:
        sys.stderr.write('ploidy_tools facs: no cell has a usable ploidy.\n')
        return 1
    if args.per_sample: calls = calls.assign(ploidy=consensus_ploidy(calls['ploidy']))
    out_dir = os.path.dirname(os.path.abspath(args.output))
    if out_dir: os.makedirs(out_dir, exist_ok=True)
    calls[['cell', 'ploidy']].to_csv(args.output, sep='\t', header=False, index=False)
    sys.stderr.write(F'Wrote {len(calls)} ploidies to the Ginkgo FACS file {args.output}\n')
    return 0

def _eval_main(args):
    """Score reported ploidies with the scAbsolute metrics implemented in ploidy_eval.py."""
    chroms = pe.AUTOSOMES if args.chroms == 'autosomes' else pe.ALL_CHROMS
    calls = read_calls(args.input, chroms=chroms, max_cn=args.max_cn)
    table = pe.load_ploidy_table(args.ploidy_file) if args.ploidy_file else None
    # Simulated data: the expected ploidy of each cell is its own simulated CN profile, read the
    # same way a caller's BED is read, so that truth and estimate are on the same footing.
    truth = (pe.load_truth_ploidies(args.truth_bed, chroms=chroms, weight='length',
                                    max_cn=args.max_cn) if args.truth_bed else None)
    run2labels = pe.load_run_labels(args.metadata_tsv, args.sample_key_columns)
    logging.info('truth: %s (%d entries); metadata: %d runs; %d ploidy calls',
                 (args.truth_bed if truth is not None else args.ploidy_file),
                 (len(truth) if truth is not None else len(table)), len(run2labels), len(calls))
    records, unresolved = [], []
    extras = [c for c in calls.columns if c not in ('cell', 'ploidy')]
    for row in calls.to_dict('records'):
        cands, run = pe.cell_candidates(str(row['cell']), run2labels, args.sample_regex or None)
        if truth is None:
            sample, expected, matched = table.lookup(cands)
        else:
            matched = pe.truth_lookup(str(row['cell']), truth)
            expected = truth[matched] if matched else float('nan')
            sample = ((args.sample_label or matched) if matched and np.isfinite(expected) else None)
        if sample is None:
            unresolved.append((row['cell'], cands[:4]))
            continue
        rec = {'cell': row['cell'], 'run': run, 'sample': sample, 'matched_label': matched}
        rec.update(pe.per_cell_metrics(float(row['ploidy']), expected, args.ploidy_window))
        rec.update({k: row[k] for k in extras})
        for k, v in (('tool', args.tool), ('donor', args.donor),
                     ('sampleType', args.sample_type), ('avgSpotLen', args.avg_spot_len)):
            if v: rec[k] = v
        records.append(rec)
    if unresolved:
        logging.warning('%d/%d cells could not be matched to a ploidy-file sample; '
                        'first unmatched labels: %s', len(unresolved), len(calls), unresolved[0][1])
        if args.strict:
            for cell, cands in unresolved[:20]: sys.stderr.write(F'unresolved: {cell} tried {cands}\n')
            return 2
    if not records:
        sys.stderr.write('ploidy_tools eval: no cell could be matched to the ploidy file.\n')
        return 1
    percell = pd.DataFrame(records)
    per_sample, overall = pe.summarize_by_sample(percell, args.ploidy_window)
    # The per-sample point estimate, appended as new `sample_*` columns so that the per-sample
    # table keeps the exact column set that ploidy_eval.py produces for the CNV callers.
    consensus = {s: consensus_ploidy(sub['observed_ploidy'])
                 for s, sub in percell.groupby('sample', sort=True)}
    per_sample['sample_ploidy'] = per_sample['sample'].map(consensus)
    expected = {s: (float(sub['expected_ploidy'].iloc[0]) if sub['expected_ploidy'].nunique() == 1
                    else consensus_ploidy(sub['expected_ploidy']))
                for s, sub in percell.groupby('sample', sort=True)}
    per_sample['sample_ploidy_error'] = per_sample['sample_ploidy'] - per_sample['sample'].map(expected)
    per_sample['sample_abs_ploidy_distance'] = per_sample['sample_ploidy_error'].abs()
    per_sample['sample_within_window'] = per_sample['sample_abs_ploidy_distance'] <= args.ploidy_window
    overall.update({
        'n_calls': int(len(calls)), 'n_cells_unresolved': len(unresolved),
        'ploidy_file': (os.path.abspath(args.ploidy_file) if args.ploidy_file else ''),
        'truth_bed': (args.truth_bed or ''), 'ploidy_calls': args.input,
        'max_cn': ('inf' if pe.is_uncapped(args.max_cn) else float(args.max_cn)),
        'mean_abs_sample_ploidy_distance': float(per_sample['sample_abs_ploidy_distance'].mean()),
        'pct_samples_within_window': float(100.0 * per_sample['sample_within_window'].mean()),
        'tool': args.tool, 'donor': args.donor,
        'sampleType': args.sample_type, 'avgSpotLen': args.avg_spot_len})
    out_dir = os.path.dirname(os.path.abspath(args.output_prefix))
    if out_dir: os.makedirs(out_dir, exist_ok=True)
    percell.sort_values(['sample', 'cell']).to_csv(args.output_prefix + '_percell.tsv', sep='\t', index=False)
    per_sample.to_csv(args.output_prefix + '_persample.tsv', sep='\t', index=False)
    with open(args.output_prefix + '_summary.json', 'w') as fh: json.dump(overall, fh, indent=2, sort_keys=True)
    if args.plot: pe.plot_ploidy(percell, per_sample, args.output_prefix, args.ploidy_window, args.title)
    sys.stderr.write(per_sample.to_string(index=False) + '\n')
    sys.stderr.write(F'mean %outliers across samples = {overall["mean_pct_outliers_across_samples"]:.1f}; '
                     F'mean |per-sample ploidy distance| = {overall["mean_abs_sample_ploidy_distance"]:.3f}\n')
    sys.stderr.write(F'Wrote {args.output_prefix}_percell.tsv, _persample.tsv, _summary.json\n')
    return 0
# ---------------------------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------------------------

def build_parser():
    p = argparse.ArgumentParser(
        description=('Benchmark ploidy-inference tools (scAbsolute) against the experimental '
                     'ploidy, with the metrics of Genome Biol 2024;25:62.'),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    sub = p.add_subparsers(dest='subcommand', required=True)
    ev = sub.add_parser('eval', help='Score a ploidy-calls TSV against a ploidy file',
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ev.add_argument('-i', '--input', nargs='+', required=True, help=(
                    'Ploidy-calls TSV files or globs, or the *intcns.bed files of a CNV caller '
                    'whose ploidy is to be derived from its own integer copy numbers'))
    ev.add_argument('-o', '--output-prefix', required=True)
    ev.add_argument('--ploidy-file', default='',
                    help='TSV mapping sample -> expected ploidy (columns: sample, ploidy, [aliases])')
    ev.add_argument('--truth-bed', nargs='+', default=None, help=(
                    F'Per-cell truth BED files or globs (typically *{pe.TRUTH_BED_SUFFIX}) written '
                    'by the CN simulator; the alternative to --ploidy-file on simulated data'))
    ev.add_argument('--sample-label', default='',
                    help='Group every cell under this sample name (used with --truth-bed)')
    ev.add_argument('--chroms', choices=['autosomes', 'all'], default='autosomes',
                    help='Chromosomes over which a truth BED is averaged into an expected ploidy')
    ev.add_argument('--max-cn', type=pe.parse_max_cn, default=pe.DEFAULT_MAX_CN,
                    help='Cap applied to the truth copy numbers before averaging; inf to not cap')
    ev.add_argument('--metadata-tsv', default='',
                    help='SraRunTable used to map run accessions in the cell names to sample labels')
    ev.add_argument('--sample-key-columns', nargs='+', default=None,
                    help=F'Metadata columns holding sample labels (default: {pe.DEFAULT_SAMPLE_KEY_COLUMNS[:4]} ...)')
    ev.add_argument('--sample-regex', default='',
                    help='Regex applied to the cell name to extract the sample label directly')
    ev.add_argument('--ploidy-window', type=float, default=pe.DEFAULT_PLOIDY_WINDOW,
                    help='Half-width of the experimental ploidy window; a cell outside it is an outlier')
    ev.add_argument('--strict', action='store_true',
                    help='Exit non-zero if any cell cannot be resolved to a ploidy-file sample')
    ev.add_argument('--plot', action='store_true', help='Also write a scAbsolute Fig. 4 style plot')
    ev.add_argument('--title', default='')
    ev.add_argument('--tool', default='')
    ev.add_argument('--donor', default='')
    ev.add_argument('--sample-type', default='')
    ev.add_argument('--avg-spot-len', default='')
    fa = sub.add_parser('facs', help='Convert a ploidy-calls TSV into a Ginkgo --facs file',
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    fa.add_argument('-i', '--input', nargs='+', required=True, help='Ploidy-calls TSV files or globs')
    fa.add_argument('-o', '--output', required=True)
    fa.add_argument('--per-sample', action='store_true',
                    help='Give every cell the single per-sample ploidy instead of its own estimate')
    return p

def main(argv=None):
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(pathname)s:%(lineno)d %(levelname)s - %(message)s')
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.subcommand == 'eval' and not (args.ploidy_file or args.truth_bed):
        parser.error('one of --ploidy-file (experimental ploidy) or --truth-bed (simulated truth) is required')
    return _eval_main(args) if args.subcommand == 'eval' else _facs_main(args)

if __name__ == '__main__':
    sys.exit(main())
