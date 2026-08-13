#!/usr/bin/env python
# Revised by: https://claude.ai/chat/1154badf-2241-4cce-9e52-7c5565f9abe3

"""
Benchmark *ploidy-inference* tools -- tools whose primary output is a ploidy estimate rather
than a per-cell copy-number profile -- starting with scAbsolute:

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
import ploidy_eval as pe

# ---------------------------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------------------------

TOOLS = ['scabsolute']
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
    # Experimental ploidy is what a ploidy tool is scored against, and only the tumor mode
    # (data_tumor.py) knows about it, so fail here instead of generating unusable rules.
    if not getattr(args, 'tumor_fastq', False): raise ValueError('--ploidy-tools requires --tumor-fastq')
    if not getattr(args, 'ploidy_file', None): raise ValueError('--ploidy-tools requires --ploidy-file')
    facs = bool(getattr(args, 'ploidy_facs', False))
    _CFG.update(tools=ploidy_tools, facs=facs, ploidy_file=os.path.abspath(args.ploidy_file),
                window=getattr(args, 'ploidy_window', pe.DEFAULT_PLOIDY_WINDOW),
                chroms=getattr(args, 'ploidy_chroms', 'autosomes'),
                metadata_tsv=(os.path.abspath(args.SraRunTable) if getattr(args, 'SraRunTable', None) else ''),
                bin_size=DEFAULT_BIN_SIZE, genome=DEFAULT_GENOME)
    for tool in ploidy_tools:
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
    if 'nop' not in tools: tools.append('nop')
    if facs and 'bam2bed' not in tools: tools.append('bam2bed')
    if not getattr(d4.run_tool_1, '_ploidy_tools', False): d4.run_tool_1 = _wrap(d4.run_tool_1)
    return tools


def _wrap(original):
    def run_tool_1(infodict, tool, *args, **kwargs):
        params = dict(zip(_POS, args)); params.update(kwargs)
        if tool in _CFG.get('tools', ()): return _gen_ploidy_tool(infodict, tool, params)
        if tool.startswith(FACS_PREFIX) and _CFG.get('facs'):
            return _gen_facs_ginkgo(original, infodict, tool, args, kwargs, params)
        return original(infodict, tool, *args, **kwargs)
    run_tool_1._ploidy_tools = True
    return run_tool_1


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
    calls, facs, prefix = find_replace_all([t_calls, t_facs, t_prefix], infodict)
    cm.makedirs((calls,))
    _CALLS[(donor, sampleType, avgSpotLen)] = dict(calls=calls, facs=facs, script=script, tool=tool)

    bams = sorted(params['inbam2call'].keys())
    bais = [change_file_ext(bam, 'bam.bai') for bam in bams]
    indir = F'{tmpdir}/{tool}_input'
    cmd = (F'rm -r {indir} || true && mkdir -p {indir} && cp -s {" ".join(bams + bais)} {indir}/ '
           F'&& {_run_cmd(tool, rootdir, indir, calls)}')
    if _CFG.get('facs'):
        cmd += F' && python {rootdir}/copy-num-bench-scwgs/ploidy_tools.py facs -i {calls} -o {facs}'
    cmd += F' #sequential=run.{tool}/'
    metadata_arg = (F'--metadata-tsv "{_CFG["metadata_tsv"]}" ' if _CFG.get('metadata_tsv') else '')
    title = F'{tool} | donor={donor} sampleType={sampleType} avgSpotLen={avgSpotLen}'
    cmd2 = (F'python {rootdir}/copy-num-bench-scwgs/ploidy_tools.py eval '
            F'-i {calls} -o {prefix} --ploidy-file "{_CFG["ploidy_file"]}" {metadata_arg}'
            F'--ploidy-window {_CFG["window"]} '
            F'--tool {tool} --donor {donor} --sample-type {sampleType} '
            F'--avg-spot-len {avgSpotLen} --plot --title "{title}" #sequential=ploidy_tool_eval.{tool}/')

    deps = []
    for fname, text in ((script, cmd), (script2, cmd2)):
        if fname in visited:
            logging.info(F'  Skip generating the script {fname} because it has already been generated. ')
            continue
        with cm.myopen(fname, wmode) as file: write2file(text, file, fname)
        visited.add(fname)
    deps.append((script, script2))
    for rule in (F'data4from2and3_1_run_DSA_{donor}_{sampleType}_{avgSpotLen}.rule',
                 F'data4from2and3_1_run_tool_{tool}.rule', 'data4from2and3_1_run_all.rule'):
        deps.append((script, rule))
    for rule in (F'data4from2and3_5_ploidy_tool_eval_DSA_{donor}_{sampleType}_{avgSpotLen}.rule',
                 F'data4from2and3_5_ploidy_tool_eval_tool_{tool}.rule',
                 'data4from2and3_5_ploidy_tool_eval_all.rule'):
        deps.append((script2, rule))
    return deps, [cmd], {}, {}


def _gen_facs_ginkgo(original, infodict, tool, args, kwargs, params):
    """A second Ginkgo pass whose per-cell ploidy is fixed to the inferred one, via --facs.

    Ginkgo is run through the unmodified run_tool_1, and only the two anchors that the FACS
    file has to reach are rewritten afterwards, the same trick gink_custom_binning.py uses.
    """
    key = (infodict['donor'], infodict['sampleType'], infodict['avgSpotLen'])
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


def read_calls(patterns):
    """Read one or more ploidy-calls TSVs (columns `cell` and `ploidy`, plus any extras)."""
    files = []
    for pat in patterns:
        files.extend(sorted(glob.glob(pat)) if any(c in pat for c in '*?[') else [pat])
    frames = []
    for path in files:
        if not (os.path.isfile(path) and os.path.getsize(path)): continue
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
    calls = read_calls(args.input)
    table = pe.load_ploidy_table(args.ploidy_file)
    run2labels = pe.load_run_labels(args.metadata_tsv, args.sample_key_columns)
    logging.info('ploidy file %s: %d samples; metadata: %d runs; %d ploidy calls',
                 args.ploidy_file, len(table), len(run2labels), len(calls))

    records, unresolved = [], []
    extras = [c for c in calls.columns if c not in ('cell', 'ploidy')]
    for row in calls.to_dict('records'):
        cands, run = pe.cell_candidates(str(row['cell']), run2labels, args.sample_regex or None)
        sample, expected, matched = table.lookup(cands)
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
    per_sample['sample_ploidy_error'] = per_sample['sample_ploidy'] - per_sample['expected_ploidy']
    per_sample['sample_abs_ploidy_distance'] = per_sample['sample_ploidy_error'].abs()
    per_sample['sample_within_window'] = per_sample['sample_abs_ploidy_distance'] <= args.ploidy_window
    overall.update({
        'n_calls': int(len(calls)), 'n_cells_unresolved': len(unresolved),
        'ploidy_file': os.path.abspath(args.ploidy_file), 'ploidy_calls': args.input,
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
    ev.add_argument('-i', '--input', nargs='+', required=True, help='Ploidy-calls TSV files or globs')
    ev.add_argument('-o', '--output-prefix', required=True)
    ev.add_argument('--ploidy-file', required=True,
                    help='TSV mapping sample -> expected ploidy (columns: sample, ploidy, [aliases])')
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
    args = build_parser().parse_args(argv)
    return _eval_main(args) if args.subcommand == 'eval' else _facs_main(args)


if __name__ == '__main__':
    sys.exit(main())
