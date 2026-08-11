import collections, gzip, logging, os, re, subprocess

### data-related variables and methods

# data1
t0into1fq1 = '<data0to1dir>/1from0.datdir/<accession>_1.fastq.gz'
t0into1fq2 = '<data0to1dir>/1from0.datdir/<accession>_2.fastq.gz'

# data1, sharded. Only the files under 1from0.datdir move; every other path template in
# this file is unchanged. A flat 1from0.datdir does not survive 10x-scale tumor data
# (thousands of cells per BAM x several BAMs = 10^4..10^5 entries in one directory), so
# tumor-mode FASTQs may instead live under a shard sub-path:
#   1from0.datdir/SRP114962/SRR633/SRR6337882_1.fastq.gz
#   1from0.datdir/PRJNA498809/colo82/colo829_10x_G1_1k_AACGTGATCCTTGCAA-1_1.fastq.gz
# t0into1fq1/t0into1fq2 above are untouched, so every existing caller keeps its behaviour.
FASTQ_SHARD_TEMPLATE   = '<study>/<accprefix>'
FASTQ_SHARD_PREFIX_LEN = 6
# 'flat' is first, and is the default everywhere: without an explicit --fastq-layout the
# resolver below is never called, so no os.path.exists probing happens at all.
FASTQ_LAYOUTS = ['flat', 'auto', 'sharded']
FASTQ_LAYOUT_DEFAULT = FASTQ_LAYOUTS[0]
t0into1fq1_sharded = '<data0to1dir>/1from0.datdir/<shard>/<accession>_1.fastq.gz'
t0into1fq2_sharded = '<data0to1dir>/1from0.datdir/<shard>/<accession>_2.fastq.gz'

# A gzip stream containing nothing is 20 bytes, and bgzip's EOF block is 28, while a
# genuine one-read FASTQ gzips to about 33 -- far too narrow a gap to judge by size. So
# only files below this bound are actually decompressed; anything larger plainly has
# content. Used to catch runs annotated PAIRED whose read 2 is absent or is a zero-read
# placeholder.
AMBIGUOUS_FASTQ_GZ_BYTES = 1024

def fastq_has_reads(path, ambiguous_bytes=AMBIGUOUS_FASTQ_GZ_BYTES):
    """True when path exists and yields at least one byte of FASTQ."""
    try: size = (os.path.getsize(path) if path else 0)
    except OSError: return False           # missing, or not a file we can stat
    if size <= 0: return False
    if size > ambiguous_bytes: return True # too large to be an empty stream
    if not str(path).endswith(('.gz', '.gzip')): return True
    try:
        with gzip.open(path, 'rb') as fh: return bool(fh.read(1))
    # BadGzipFile subclasses OSError, but a truncated stream raises EOFError instead.
    except (OSError, EOFError): return False

def sanitize_path_part(part, default='NA'):
    cleaned = re.sub(r'[^A-Za-z0-9._+-]', '-', str(part if part is not None else '').strip())
    return cleaned if cleaned and cleaned.lower() != 'nan' else default

def fastq_shard(accession, study=None, template=None, prefix_len=None):
    """Render the shard sub-path (no leading/trailing separator) for one accession."""
    template   = (template   if template   else FASTQ_SHARD_TEMPLATE)
    prefix_len = (prefix_len if prefix_len else FASTQ_SHARD_PREFIX_LEN)
    acc = sanitize_path_part(accession)
    rendered = template
    for key, val in (('study', sanitize_path_part(study)), ('accprefix', acc[:prefix_len] or acc),
                     ('accession', acc)):
        rendered = rendered.replace('<' + key + '>', val)
    return os.path.sep.join([p for p in rendered.replace('\\', '/').split('/') if p])

def fastq_pair(data0to1dir, accession, study=None, layout='auto', template=None, prefix_len=None):
    """Resolve (read1, read2) for one accession.

    layout='auto' (the default) prefers the sharded path when that file is already on
    disk, then the flat path when that one is, and otherwise falls back to flat. So a
    tree populated the old way yields exactly the strings t0into1fq1/t0into1fq2 would
    have yielded, and a tree holding both old and new files resolves per accession.
    """
    infodict = {'data0to1dir': data0to1dir, 'accession': str(accession),
                'shard': fastq_shard(accession, study, template, prefix_len)}
    sharded = find_replace_all([t0into1fq1_sharded, t0into1fq2_sharded], infodict)
    flat    = find_replace_all([t0into1fq1,         t0into1fq2        ], infodict)
    if layout == 'sharded': return sharded
    if layout == 'flat':    return flat
    assert layout == 'auto', F'Unknown FASTQ layout {layout}; expected one of {FASTQ_LAYOUTS}'
    if os.path.exists(sharded[0]): return sharded
    if os.path.exists(flat[0]):    return flat
    return flat

# data1to2
t1into2log = '<data1to2dir>/<donor>/1into2_2_<donor>.logdir/'
t1into2sh0 = '<data1to2dir>/<donor>/1into2_2_<donor>.logdir/2_<donor>_1into2_a.sentinel_begin.sh'
t1into2sh1 = '<data1to2dir>/<donor>/1into2_2_<donor>.logdir/2_<donor>_1into2_step1-gen-bam_<accession>.sh'
t1into2sh2 = '<data1to2dir>/<donor>/1into2_2_<donor>.logdir/2_<donor>_1into2_step2-gen-vcf.sh'
t1into2sh3 = '<data1to2dir>/<donor>/1into2_2_<donor>.logdir/2_<donor>_1into2_step3-gen-mut-unsort-fqs_<accession>_<GT>.sh'
t1into2sh4 = '<data1to2dir>/<donor>/1into2_2_<donor>.logdir/2_<donor>_1into2_step4-gen-mut-sorted-fqs_<accession>_<GT>.sh'
t1into2sh5 = '<data1to2dir>/<donor>/1into2_2_<donor>.logdir/2_<donor>_1into2_step5-gen-mut-bam_<accession>_<GT>.sh'
t1into2end = '<data1to2dir>/<donor>/1into2_2_<donor>.logdir/2_<donor>_1into2_z.sentinel_end.sh'
t1into2tmp = '<data1to2dir>/<donor>/1into2_2_<donor>.tmpdir/'

# data2
t2from1datdir = '<data1to2dir>/<donor>/2from1_2_<donor>.datdir/'
t2from1vcf010 = '<data1to2dir>/<donor>/2from1_2_<donor>.datdir/2_<donor>_2from1_step2out_10_init.vcf.gz'
t2from1vcf020 = '<data1to2dir>/<donor>/2from1_2_<donor>.datdir/2_<donor>_2from1_step2out_20_randomHaploPair.vcf.gz'
t2from1vcf021 = '<data1to2dir>/<donor>/2from1_2_<donor>.datdir/2_<donor>_2from1_step2out_21_randomHaploPair_diplotype.vcf'
t2from1vcf03A = '<data1to2dir>/<donor>/2from1_2_<donor>.datdir/2_<donor>_2from1_step2out_3A_haploA.vcf.gz'
t2from1vcf03B = '<data1to2dir>/<donor>/2from1_2_<donor>.datdir/2_<donor>_2from1_step2out_3B_haploB.vcf.gz'
t2from1mutbam = '<data1to2dir>/<donor>/2from1_2_<donor>.datdir/2_<donor>_2from1_step5out_<accession>_<GT>_sort_markdup_mut.bam'
t2from1dedupb = '<data1to2dir>/<donor>/2from1_2_<donor>.datdir/2_<donor>_2from1_step5out_<accession>_<GT>_sort_markdup_mut_dedup.bam'
t2from1cnbams=[F'<data1to2dir>/<donor>/2from1_2_<donor>.datdir/2_<donor>_2from1_step5out_<accession>_<GT>_sort_markdup_mut_cn{i}.bam' for i in range(1, 1+6, 1)]
t2from1flagjs = '<data1to2dir>/<donor>/2from1_2_<donor>.datdir/2_<donor>_2from1_step5out_<accession>_<GT>_sort_markdup_mut_flagstat.json'

# run the above first if samtools flagstat is used for downsampling sequencing reads

# data2to3
t2into3logdir = '<data2to3dir>/<donor>/2into3_2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>.logdir/'
t2into3script = '<data2to3dir>/<donor>/2into3_2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>.logdir/2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>_2into3_<accession_1>_<accession_2>.sh'
t2into3end    = '<data2to3dir>/<donor>/2into3_2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>.logdir/2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>_2into3_z.sentinel_end.sh'
t2into3tmpdir = '<data2to3dir>/<donor>/2into3_2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>.tmpdir/'

# data3
t3from2datdir = '<data2to3dir>/<donor>/3from2_2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>.datdir/'
t3from2simbam = '<data2to3dir>/<donor>/3from2_2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>.datdir/2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>_3from2_<accession_1>_<accession_2>.bam'
t3from2simbed = '<data2to3dir>/<donor>/3from2_2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>.datdir/2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>_3from2_<accession_1>_<accession_2>_simtruth.bed'
t3from2dedupb = '<data2to3dir>/<donor>/3from2_2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>.datdir/2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>_3from2_<accession_1>_<accession_2>_dedup.bam'
t3from2infojs = '<data2to3dir>/<donor>/3from2_2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>.datdir/2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>_3from2_<accession_1>_<accession_2>_info.json'

# data2to4
t2into4logdir = '<data2to4dir>/<donor>/2into4_2_<donor>_3_<sampleType>_<avgSpotLen>_4_step<tool_order>_<tool>.logdir/'
t2into4script = '<data2to4dir>/<donor>/2into4_2_<donor>_3_<sampleType>_<avgSpotLen>_4_step<tool_order>_<tool>.logdir/2_<donor>_3_<sampleType>_<avgSpotLen>_4_step<tool_order>_<tool>_2into4_call.sh'
t2into4scrip2 = '<data2to4dir>/<donor>/2into4_2_<donor>_3_<sampleType>_<avgSpotLen>_4_step<tool_order>_<tool>.logdir/2_<donor>_3_<sampleType>_<avgSpotLen>_4_step<tool_ord_1>_<tool>_2into4_norm.sh'
t2into4tmpdir = '<data2to4dir>/<donor>/2into4_2_<donor>_3_<sampleType>_<avgSpotLen>_4_step<tool_order>_<tool>.tmpdir/'

# data3to4
t3into4logdir = '<data3to4dir>/<donor>/3into4_2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>_4_step<tool_order>_<tool>.logdir/'
t3into4script = '<data3to4dir>/<donor>/3into4_2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>_4_step<tool_order>_<tool>.logdir/2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>_4_step<tool_order>_<tool>_3into4_call.sh'
t3into4scrip2 = '<data3to4dir>/<donor>/3into4_2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>_4_step<tool_order>_<tool>.logdir/2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>_4_step<tool_ord_1>_<tool>_3into4_norm.sh'

t3into4tmpdir = '<data3to4dir>/<donor>/3into4_2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>_4_step<tool_order>_<tool>.tmpdir/'

# data4

t4from2datdir = '<data2to4dir>/<donor>/4from2_2_<donor>_3_<sampleType>_<avgSpotLen>_4_step<tool_order>_<tool>.datdir/'
t4from2depcns = '<data2to4dir>/<donor>/4from2_2_<donor>_3_<sampleType>_<avgSpotLen>_4_step<tool_order>_<tool>.datdir/2_<donor>_3_<sampleType>_<avgSpotLen>_4_step<tool_order>_<tool>_4from2_<samplename>_depcns.bed'
t4from2intcns = '<data2to4dir>/<donor>/4from2_2_<donor>_3_<sampleType>_<avgSpotLen>_4_step<tool_order>_<tool>.datdir/2_<donor>_3_<sampleType>_<avgSpotLen>_4_step<tool_order>_<tool>_4from2_<samplename>_intcns.bed'

t4from3datdir = '<data2to4dir>/<donor>/4from3_2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>_4_step<tool_order>_<tool>.datdir/'
t4from3depcns = '<data2to4dir>/<donor>/4from3_2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>_4_step<tool_order>_<tool>.datdir/2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>_4_step<tool_order>_<tool>_4from3_<samplename>_depcns.bed'
t4from3intcns = '<data2to4dir>/<donor>/4from3_2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>_4_step<tool_order>_<tool>.datdir/2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>_4_step<tool_order>_<tool>_4from3_<samplename>_intcns.bed'

t4into5logdir = '<data3to4dir>/<donor>/4into5_2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>_4_step<tool_order>_<tool>.logdir/'
t4into5script = '<data3to4dir>/<donor>/4into5_2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>_4_step<tool_order>_<tool>.logdir/2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>_4_step<tool_ord_2>_<tool>_4into5_eval.sh'

def find_replace_all(args, old2new, prefix='<', suffix='>'):
    ret = []
    for arg1 in args:
        assert isinstance(arg1, (str, list, tuple, set, dict)), F'The variable {arg1} is neither a string nor an iterable!'
        if isinstance(arg1, str):
            arg2 = arg1
            for old, new in sorted(old2new.items()):
                if new:
                    arg2 = arg2.replace(prefix+old+suffix, new)
            ret.append(arg2)
        else:
            ret.append(list(find_replace_all(arg1, old2new, prefix, suffix)))
    return tuple(ret)

def get_varnames(default_value, varnames=['data0to1dir', 'data1to2dir', 'data2to3dir', 'data2to4dir', 'data3to4dir', 'data4to5dir']):
    ret = []
    for varname in varnames:
        v = os.getenv(varname, default_value)
        ret.append(v)
    return ret

### OS-related variables and methods

OVERWRITING_PREVENTION_MODES = ['no_overwriting', 'no']
DEFAULT_WRITING_MODE = OVERWRITING_PREVENTION_MODES[0]
def myopen(filename, mode):
    if mode in OVERWRITING_PREVENTION_MODES and os.path.exists(filename):
        logging.info(F'Redirect {filename} to {os.devnull} to prevent overwriting')
        return open(os.devnull, 'w')
    elif mode in OVERWRITING_PREVENTION_MODES:
        return open(filename, 'w')
    else:
        return open(filename, mode)

def change_file_ext(file_path, new_extension, old_extension=''):
    base_name, old_ext = os.path.splitext(file_path)
    # Bug found and fixed by: https://sorryios.ai/chat/8d3002cb-af21-4302-9c97-1db112d059f7
    # But this fix does not require any rerun since old_extension='' everywhere in the past and as of now.
    if old_extension: assert '.'+old_extension == old_ext, F'File extension check: {old_extension} == {old_ext} failed!'
    return base_name + "." + new_extension

def makedirs(args):
    for arg in args: os.makedirs(os.path.dirname(arg), exist_ok=True)

def write2file(cmd, file, filename):
    assert '<' not in filename and '>' not in filename, F'The filename {filename} is invalid!'
    ret = file.write(cmd + '\n')
    #print(cmd + F' # written to file {filename}')
    assert ret > 0, F'Writing ({cmd}) to {filename} failed!'
    return ret

# Apache airflow is better than Snakemake in terms of engineering, but it requires heavy setup, so Snakemake is used by default
def gen_airflow_content(tasks, dependencies):
    return F'''
from airflow import DAG
from airflow.operators.bash_operator import BashOperator
from airflow.operators.empty import EmptyOperator
from datetime import datetime

default_args = {{
    'owner': 'zxf',
    'start_date': datetime(2025, 1, 17),
}}

dag = DAG(
    dag_id="copy-num-bench-scwgs",
    schedule="0 0 * * *",
    catchup=False,
    default_args=default_args,
    #schedule_interval=None,  # If you want to trigger manually
)

# Define your bash scripts as tasks
{tasks}

run_this_last = EmptyOperator(
        task_id="run_this_last",
    )

#task1 = BashOperator(
#    task_id='task1',
#    bash_command='bash /path/to/script1.sh',
#    dag=dag,
#)

#task2 = BashOperator(
#    task_id='task2',
#    bash_command='bash /path/to/script2.sh',
#    dag=dag,
#)

# Set dependencies
{dependencies}

#task1 >> task2  # task1 must complete before task2 starts
'''
def script2task(script): return 'task_' + re.sub(r'[^a-zA-Z0-9]', '_', script.split(os.path.sep)[-1])
def list2airflow(listof_src_dst):
    script2string = {}
    dependencies = []
    visited_tasks = set([])
    for src_dst in sorted(listof_src_dst):
        src_script, dst_script = src_dst[0], src_dst[1]
        src_task, dst_task = script2task(src_script), script2task(dst_script)
        for script, task in zip([src_script, dst_script], [src_task, dst_task]):
            if script not in script2string:
                script2string[script] = F'''{task} = BashOperator(task_id='{task}', bash_command=" bash -vx {script} " , dag=dag)'''
        dependencies.append(F'''{src_task} >> {dst_task}''')
        if dst_task not in visited_tasks:
            dependencies.append(F'''{dst_task} >> run_this_last''')
            visited_tasks.add(dst_task)
    return gen_airflow_content('\n'.join(sorted(script2string.values())), '\n'.join(dependencies))


def list2snakemake(listof_src_dst, allrules='all'):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    init_version_info = subprocess.check_output(F'cd {script_dir} && git rev-parse HEAD && git diff HEAD', shell=True, text=True)
    init_version_comment = '\n'.join([F'#{line}' for line in init_version_info.split('\n')])
    rules = []
    all_scripts = set([])
    dst_to_srcs_dict = collections.defaultdict(set)
    rule_to_srcs_dict = collections.defaultdict(set)
    dst_to_params_dict = {}
    dst_dones = []   
    for src_dst in sorted(listof_src_dst):
        src_script, dst_script = src_dst[0], src_dst[1]
        if dst_script.endswith('.rule'):
            rule_to_srcs_dict[dst_script].add(src_script)
        else:
            dst_to_srcs_dict[dst_script].add(src_script)
            if len(src_dst) > 2:
                dst_to_params_dict[dst_script] = src_dst[2]
            all_scripts.add(src_script)
            all_scripts.add(dst_script)
    
    for dst_script in sorted(all_scripts):
        src_scripts = sorted(list(dst_to_srcs_dict[dst_script]))
        dst_task = script2task(dst_script)
        dst_done = F'"{dst_script}.done"'
        dst_dones.append(dst_done)
        params = '\n'.join(dst_to_params_dict.get(dst_script, []))
        if len(src_scripts) == 0:
            rule = F'''
# 1. Beginning
rule {dst_task}:
    input: "{dst_script}"
    output: {dst_done}
    {params}
    shell: "command time -v bash -evx {dst_script} 2> {dst_script}.stderr && (pushd {script_dir} && git rev-parse HEAD && git diff HEAD) > {dst_script}.done"'''
        else:    
            src_dones = [F'"{s}.done"' for s in src_scripts]
            rule = F'''
# 2. Middle
rule {dst_task}:
    input: {', '.join(src_dones)}
    output: {dst_done}
    {params}
    shell: "command time -v bash -evx {dst_script} 2> {dst_script}.stderr && (pushd {script_dir} && git rev-parse HEAD && git diff HEAD) > {dst_script}.done"'''
        rules.append(rule)
    for dst_rulename in sorted(rule_to_srcs_dict.keys()):
        dst_task = script2task(dst_rulename)
        src_scripts = sorted(list(rule_to_srcs_dict[dst_rulename]))
        rule = F'''
# 3. End
rule {dst_task}:
    input: "{', '.join(src_scripts)}"'''
        rules.append(rule)

    rule_all = F'''
### GIT_COMMIT_INFO
{init_version_comment}
rule {allrules}:
    input: {', '.join(dst_dones)}
'''
    return '\n'.join([rule_all] + sorted(rules))

### genomic-related variables and methods

MAX_HAPLO_CN = 6

chrsNS = 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22'
chrs   = 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY'

def norm_sample_type(df0):
    return (df0['sample-type']
            .str.replace('single~cell', 'sperm')
            .str.replace('Oocyte~product:~Second~poloar~body', 'PB2')
            .str.replace('Oocyte~product:~First~poloar~body' , 'PB1')
            .str.replace('Oocyte~product:~female~pro-nucleus', 'FPN'))


### math-related variables and methods

def circular_dist(n1, n2, nmax):
    assert n1 < nmax, F'{n1} < {nmax} failed'
    assert n2 < nmax, F'{n2} < {nmax} failed'
    n1a = n1 + nmax
    n2a = n2 + nmax
    return min([abs(n1 - n2), (n1a - n2),  (n2a - n1)])

def circular_dist_below(n1, n2, nmax, dist=1): return circular_dist(n1, n2, nmax) <= dist

