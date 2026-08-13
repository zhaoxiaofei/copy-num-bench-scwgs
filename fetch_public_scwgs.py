#!/usr/bin/env python3
"""
Download public scWGS data (ENA/SRA runs or 10x Cell-Ranger-DNA BAMs), convert it to
per-cell FASTQ files, and emit a metadata TSV that main.py --tumor-fastq accepts.

Output layout (matches common.t0into1fq1 / common.t0into1fq2):

    <outdir>/1from0.datdir/<accession>_1.fastq.gz
    <outdir>/1from0.datdir/<accession>_2.fastq.gz     # paired-end only
    <outdir>/<dataset_id>.SraRunTable.tsv            # feed this to main.py

The metadata TSV uses the same column names as real_tumor/GIAB_PRJNA200694_*.tsv, so it
is a drop-in for the existing tumor-mode entry point:

    python fetch_public_scwgs.py --dataset tnbc_kim2018 --outdir ../real_tumor_data
    python main.py --tumor-fastq \
        --SraRunTable ../real_tumor_data/tnbc_kim2018.SraRunTable.tsv > Snakefile
    snakemake --cores 200

Modes
  ena       One ENA/SRA run == one cell. Pulls fastq_ftp straight from the ENA portal
            API (works for SRP/PRJNA/PRJEB/E-MTAB/GSE-backed studies).
  sra       Same, but via prefetch + a dump tool. Fallback when ENA has no FASTQ
            mirror (some submitter-supplied BAM-only studies).
  sra_file  Download the raw .sra archive per run (lftp, parallel + resumable), then
            fastq-dump/fasterq-dump it locally. Same result as `sra`, but prefetch is a
            single stream and is usually the slowest part of a large fetch; pulling the
            .sra with `pget -n` instead is typically several times faster, and the
            download is separately resumable and md5-verifiable.
  tenx_bam  One 10x position-sorted BAM == thousands of cells. Splits on the CB tag
            with `samtools split -d CB`, filters background barcodes, then converts
            each per-cell BAM to a name-collated FASTQ pair.
  bam       A single BAM == a single cell/sample.

Everything is resumable: each output gets a `.done` sentinel and is skipped on rerun.

Requires: samtools >= 1.14 (for `samtools split -d`), and sra-tools for --mode sra.
"""

import argparse
import collections
import concurrent.futures
import csv
import gzip
import hashlib
import logging
import os
import re
import shutil
import subprocess
import sys
import urllib.parse
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as cm   # noqa: E402  (shard template + layout resolver live here)

ENA_PORTAL = 'https://www.ebi.ac.uk/ena/portal/api/filereport'
ENA_FIELDS = ('run_accession,experiment_accession,sample_accession,study_accession,'
              'library_name,sample_title,sample_alias,instrument_platform,'
              'library_layout,read_count,base_count,fastq_ftp,fastq_md5,'
              'submitted_ftp,scientific_name,sra_ftp,sra_md5,sra_bytes,fastq_bytes,'
              'secondary_study_accession')

# Fallbacks when ENA has no sra_ftp entry. <run> is substituted with the run accession.
# The AWS Open Data mirror needs no credentials and is what prefetch usually resolves to.
SRA_URL_TEMPLATES = ('https://sra-pub-run-odp.s3.amazonaws.com/sra/<run>/<run>',)

# cm.fastq_has_reads arrived with the LibraryLayout fix; degrade gracefully without it.
_has_reads = getattr(cm, 'fastq_has_reads', lambda p: bool(p) and os.path.exists(p))

# Column order of real_tumor/GIAB_PRJNA200694_SraRunTable_bioskryb_ILM.tsv. data_tumor.py
# rewrites ' ' -> '~' in headers, so 'SRA Study' becomes the 'SRA~Study' it looks for.
# 'Fastq1'/'Fastq2' are appended so the run table records the sharded paths explicitly;
# data_tumor.py already prefers those columns over any path convention, which makes the
# emitted table independent of --fastq-layout.
META_COLUMNS = ['Run', 'AvgSpotLen', 'Library Name', 'Sample Name', 'sample_type',
                'SRA Study', 'Bases', 'Bytes', 'SEQUENCING_CENTER (run)', 'isolate',
                'Platform', 'LibraryLayout', 'Fastq1', 'Fastq2']

log = logging.getLogger('fetch_public_scwgs')


# ----------------------------------------------------------------------------- helpers

def item_name(item):
    """A short label for an item: ENA rows are whole dicts and unreadable in a log line."""
    if isinstance(item, dict):
        return item.get('run_accession') or item.get('Run') or '<row>'
    return str(item)


def run_parallel(fn, items, jobs, label='item'):
    """Map fn over items with a thread pool, keeping the input order in the results.

    Fetching is I/O bound -- an lftp/prefetch subprocess, or a urllib socket -- so threads
    are the right tool here: no pickling of arguments, and the GIL is released for the
    whole of both. Returns (results, failures); a slot is None where fn raised, so one
    dead run out of a thousand does not throw the other 999 away.
    """
    results, failures = [None] * len(items), []
    if jobs <= 1:
        for idx, item in enumerate(items):
            try:
                results[idx] = fn(item)
            except (Exception, SystemExit) as exc:                  # noqa: BLE001
                log.error('%s failed: %s (%s)', label, item_name(item), exc)
                failures.append((item, exc))
        return results, failures

    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = {pool.submit(fn, item): idx for idx, item in enumerate(items)}
        done = 0
        for future in concurrent.futures.as_completed(futures):
            idx = futures[future]
            done += 1
            try:
                results[idx] = future.result()
            except (Exception, SystemExit) as exc:                  # noqa: BLE001
                log.error('%s failed: %s (%s)', label, item_name(items[idx]), exc)
                failures.append((items[idx], exc))
            if done % 50 == 0 or done == len(items):
                log.info('%d/%d %ss done (%d failed)', done, len(items), label, len(failures))
    return results, failures


def report_failures(failures, outdir, dataset_id, total):
    """Log a summary and leave a retry list, so a partial fetch is never silent."""
    if not failures:
        return None
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, F'{dataset_id}.failed.tsv')
    log.warning('%d/%d failed; the run table below covers only the %d that succeeded',
                len(failures), total, total - len(failures))
    try:
        with open(path, 'w') as fh:
            fh.write('#item\terror\n')
            for item, exc in failures:
                fh.write(F'{item_name(item)}\t{str(exc)[:300]}\n')
        log.warning('retry list written to %s', path)
        return path
    except OSError as exc:
        log.warning('could not write the retry list (%s)', exc)
        return None


def run_cmd(cmd, dry_run=False, **kwargs):
    log.info('+ %s', cmd if isinstance(cmd, str) else ' '.join(cmd))
    if dry_run:
        return 0
    return subprocess.run(cmd, shell=isinstance(cmd, str), check=True, **kwargs).returncode


def require_tool(name, hint=''):
    if shutil.which(name) is None:
        raise SystemExit(F'Required executable not found on PATH: {name}. {hint}')
    return shutil.which(name)


def samtools_supports_split_by_tag():
    """`samtools split -d TAG` landed in 1.14; older builds silently ignore -d."""
    try:
        res = subprocess.run(['samtools', 'split', '--help'], capture_output=True, text=True)
        out = (res.stdout or '') + (res.stderr or '')
    except OSError:
        return False
    return '-d TAG' in out


def is_done(path):
    return os.path.exists(path + '.done')


def mark_done(path, dry_run=False):
    if not dry_run:
        with open(path + '.done', 'w') as fh:
            fh.write('ok\n')


def md5_of(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(chunk), b''):
            h.update(block)
    return h.hexdigest()


def lftp_available():
    return shutil.which('lftp') is not None


def lftp_download(url, dest, conns=1000-1, dry_run=False):
    """Fetch one URL with lftp pget. lftp -O takes a *directory* and keeps the URL
    basename, so download there first and rename afterwards when the desired filename
    differs (ENA single-end runs arrive as <run>.fastq.gz but we want <run>_1.fastq.gz).
    """
    dest_dir = os.path.dirname(os.path.abspath(dest))
    if not dry_run: os.makedirs(dest_dir, exist_ok=True)
    landed = os.path.join(dest_dir, os.path.basename(urllib.parse.urlparse(url).path))
    run_cmd(F'lftp -c "set net:reconnect-interval-base 5; set net:max-retries 10; '
            F'pget -c -n {conns} -O {dest_dir} {url}"', dry_run)
    if not dry_run and landed != os.path.abspath(dest):
        os.replace(landed, dest)
    return dest


def download(url, dest, md5=None, dry_run=False, retries=3, downloader='lftp', conns=1000-1):
    """Download to dest.part then rename, so a killed job never leaves a truncated file."""
    if is_done(dest):
        log.info('skip (already done): %s', dest)
        return dest
    if not url.startswith(('http://', 'https://', 'ftp://')):
        url = 'https://' + url          # ENA returns bare ftp.sbi... host paths
    if not dry_run: os.makedirs(os.path.dirname(dest), exist_ok=True)
    log.info('GET %s -> %s', url, dest)
    if downloader == 'lftp' and lftp_available():
        lftp_download(url, dest, conns=conns, dry_run=dry_run)
        if md5 and md5 != 'NA' and not dry_run:
            got = md5_of(dest)
            if got != md5:
                os.remove(dest)
                raise RuntimeError(F'md5 mismatch for {url}: expected {md5}, got {got}')
        mark_done(dest, dry_run)
        return dest
    if downloader == 'lftp':
        log.warning('lftp not on PATH; falling back to the pure-python downloader')
    if dry_run:
        return dest
    part = dest + '.part'
    last = None
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=120) as resp, open(part, 'wb') as out:
                shutil.copyfileobj(resp, out, length=1 << 22)
            break
        except Exception as exc:                                   # noqa: BLE001
            last = exc
            log.warning('attempt %d/%d failed for %s: %s', attempt, retries, url, exc)
    else:
        raise RuntimeError(F'Download failed after {retries} attempts: {url} ({last})')
    if md5 and md5 != 'NA':
        got = md5_of(part)
        if got != md5:
            os.remove(part)
            raise RuntimeError(F'md5 mismatch for {url}: expected {md5}, got {got}')
    os.replace(part, dest)
    mark_done(dest)
    return dest


def url_exists(url):
    req = urllib.request.Request(url, method='HEAD')
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            return 200 <= resp.status < 400
    except Exception:                                              # noqa: BLE001
        return False


# ------------------------------------------------------------------------ ENA metadata

def ena_filereport(accession, fields=ENA_FIELDS):
    """Query the ENA portal API. Accepts SRP/ERP/DRP, PRJNA/PRJEB, SRR/ERR, E-MTAB-*."""
    query = urllib.parse.urlencode({'accession': accession, 'result': 'read_run',
                                    'fields': fields, 'format': 'tsv', 'limit': 0})
    url = F'{ENA_PORTAL}?{query}'
    log.info('ENA filereport: %s', url)
    with urllib.request.urlopen(url, timeout=180) as resp:
        text = resp.read().decode('utf-8')
    rows = list(csv.DictReader(text.splitlines(), delimiter='\t'))
    if not rows:
        raise SystemExit(F'ENA returned no runs for {accession}. Check the accession, or '
                         F'use --mode sra if the study is BAM-only.')
    log.info('ENA returned %d runs for %s', len(rows), accession)
    return rows


# Only these look like INSDC study identifiers that the ENA portal can resolve. Anything
# else in the accession column (a 10x dataset slug, a synthetic study id) has no second
# accession by construction and is marked 'none' rather than left as an unresolved 'NA'.
INSDC_STUDY_RE = re.compile(r'^(PRJ[EDN][A-Z]\d+|[SED]RP\d+|E-\w+-\d+|SRA\d+)$')
GEO_SERIES_RE = re.compile(r'^GS[EM]\d+$')


def fill_alt_accessions(args):
    """Resolve every manifest row's paired accession from ENA and rewrite the TSV.

    Each INSDC study carries both a BioProject (PRJNA/PRJEB) and a secondary study
    (SRP/ERP) accession -- PRJNA629885 and SRP259526 are the same Minussi 2021 data.
    Which one a paper cites is arbitrary and ENA does not always index both, so recording
    the pair makes each row resolvable from either. This asks ENA rather than guessing.
    """
    with open(args.manifest) as fh:
        lines = fh.read().splitlines()
    header = lines[0].lstrip('#').split('\t')
    if 'alt_accession' not in header:
        raise SystemExit(F'{args.manifest} has no alt_accession column.')
    i_acc, i_alt = header.index('accession'), header.index('alt_accession')

    out, changed, unresolved = [lines[0]], 0, []
    for line in lines[1:]:
        if not line.strip() or line.startswith('#'):
            out.append(line)
            continue
        cols = line.split('\t')
        acc, alt = cols[i_acc], cols[i_alt]
        if alt and alt.upper() not in ('NA', '') and not args.force:
            out.append(line)
            continue

        if GEO_SERIES_RE.match(acc):
            # ENA does not index GEO series ids; the linked SRP is on the GEO page.
            unresolved.append(F'{acc} (GEO series: take the SRP from its GEO page)')
            out.append(line)
            continue
        if not INSDC_STUDY_RE.match(acc):
            cols[i_alt] = 'none'                 # 10x-hosted or synthetic: no pair exists
            out.append('\t'.join(cols))
            changed += 1
            continue

        try:
            rows = ena_filereport_any([acc], fields='study_accession,secondary_study_accession')
        except (Exception, SystemExit) as exc:                      # noqa: BLE001
            unresolved.append(F'{acc} ({exc})')
            out.append(line)
            continue
        pair = {rows[0].get('study_accession') or '', rows[0].get('secondary_study_accession') or ''}
        other = sorted(a for a in pair if a and a != acc)
        if not other:
            unresolved.append(F'{acc} (ENA reported no second accession)')
            out.append(line)
            continue
        cols[i_alt] = other[0]
        log.info('%s -> %s', acc, other[0])
        out.append('\t'.join(cols))
        changed += 1

    if args.dry_run:
        print('\n'.join(out))
    else:
        with open(args.manifest, 'w') as fh:
            fh.write('\n'.join(out) + '\n')
    log.info('filled %d row(s)%s', changed, '' if not args.dry_run else ' (dry-run, not written)')
    for note in unresolved:
        log.warning('still unresolved: %s', note)
    return 0


def inspect_accession(args):
    """Summarise an ENA/SRA study without downloading anything.

    Answers the questions that decide whether a study belongs in the manifest: how many
    runs, what platform and species, single or paired, how big, and -- the one that bites
    downstream -- how many runs are annotated PAIRED while only one FASTQ was submitted.
    """
    rows = ena_filereport_any([args.accession, args.alt_accession])

    def tally(fn):
        return dict(collections.Counter(fn(r) for r in rows))

    def observed(row):
        ftp = [u for u in (row.get('fastq_ftp') or '').split(';') if u]
        if not ftp:
            return 'NO-FASTQ'
        mates = sum(u.endswith(('_1.fastq.gz', '_2.fastq.gz')) for u in ftp)
        return 'PAIRED' if mates >= 2 else 'SINGLE'

    def total(field):
        return sum(int(r.get(field) or 0) for r in rows)

    mismatched = [r['run_accession'] for r in rows
                  if (r.get('library_layout') or '').upper() == 'PAIRED'
                  and observed(r) == 'SINGLE']
    # avg_spot_len returns strings; min/max on those compares lexicographically ('150'
    # sorts below '75'), so coerce to int before summarising.
    spots = sorted(int(avg_spot_len(r) or 0) for r in rows)
    has_fastq = any((r.get('fastq_ftp') or '') for r in rows)

    out = [
        F'accession         : {args.accession}',
        F'runs              : {len(rows)}',
        F'study (BioProject): ' + ', '.join(sorted(tally(lambda r: r.get('study_accession') or '?'))),
        F'study (SRA)       : ' + ', '.join(sorted(tally(lambda r: r.get('secondary_study_accession') or '?'))),
        F'species           : {tally(lambda r: r.get("scientific_name") or "?")}',
        F'platform          : {tally(lambda r: r.get("instrument_platform") or "?")}',
        F'layout, annotated : {tally(lambda r: (r.get("library_layout") or "?").upper())}',
        F'layout, on ENA    : {tally(observed)}',
        F'read length       : min {min(spots, default="?")}, max {max(spots, default="?")}',
        F'total bases       : {total("base_count"):,}',
        F'fastq size        : {total("fastq_bytes") / 1e9:.1f} GB '
        F'(.sra: {total("sra_bytes") / 1e9:.1f} GB)',
        F'PAIRED-but-SINGLE : {len(mismatched)} run(s)'
        + (' -> pass --infer-library-layout to main.py' if mismatched else ''),
        F'example runs      : ' + ', '.join(r['run_accession'] for r in rows[:4]),
        '',
        F'Candidate manifest row (fill in TODO, then append to {os.path.basename(args.manifest)}):',
        '\t'.join([args.accession.lower(), ('ena' if has_fastq else 'sra_file'),
                   rows[0].get('study_accession') or args.accession,
                   rows[0].get('secondary_study_accession') or 'NA',
                   'TODO-donor', args.sample_type,
                   (rows[0].get('instrument_platform') or 'ILLUMINA').upper(),
                   'open', 'NA', 'TODO-notes']),
    ]
    print('\n' + '\n'.join(out))
    if mismatched:
        log.warning('annotated PAIRED but only one FASTQ on ENA: %s',
                    ', '.join(mismatched[:10]) + (' ...' if len(mismatched) > 10 else ''))
    return 0


def sra_url_candidates(run, ena_row=None, extra_template=None):
    """Ordered .sra URLs to try for one run: explicit template, then ENA, then AWS."""
    urls = []
    if extra_template:
        urls.append(extra_template.replace('<run>', run))
    for url in ((ena_row or {}).get('sra_ftp') or '').split(';'):
        if url.strip():
            urls.append(url.strip())            # ENA mirror; carries an md5 we can check
    urls += [t.replace('<run>', run) for t in SRA_URL_TEMPLATES]
    seen, ordered = set(), []
    for url in urls:
        if url not in seen:
            seen.add(url)
            ordered.append(url)
    return ordered


def dump_sra_to_fastq(args, source, run, celldir):
    """Dump `source` (a local .sra path, or a bare accession the toolkit resolves) into
    <celldir>/<run>_{1,2}.fastq.gz. Returns (fq1, fq2, layout); fq2 is '' for single-end.

    --split-3 writes <run>_1/<run>_2 when mates exist and <run> when they do not, so the
    files it produces settle the layout no matter what the run was annotated as.
    """
    fq1 = os.path.join(celldir, F'{run}_1.fastq.gz')
    fq2 = os.path.join(celldir, F'{run}_2.fastq.gz')
    if is_done(fq1):
        log.info('skip (already done): %s', fq1)
        return (fq1, fq2, 'PAIRED') if _has_reads(fq2) else (fq1, '', 'SINGLE')

    tmp = os.path.join(args.tmpdir, F'dump.{run}')
    if not args.dry_run:
        os.makedirs(tmp, exist_ok=True)
    if args.dump_tool == 'fastq-dump':
        # fastq-dump gzips as it writes: no second pass, and a far smaller disk peak.
        require_tool('fastq-dump', 'Install sra-tools.')
        run_cmd(['fastq-dump', '--split-3', '--gzip', '--outdir', tmp,
                 *args.dump_args, str(source)], args.dry_run)
        suffix = '.fastq.gz'
    else:
        require_tool('fasterq-dump', 'Install sra-tools.')
        run_cmd(['fasterq-dump', '--split-3', '--threads', str(args.threads),
                 '--outdir', tmp, '-t', tmp, *args.dump_args, str(source)], args.dry_run)
        suffix = '.fastq'

    def collect(name, dest):
        src = os.path.join(tmp, name)
        if args.dry_run:
            log.info('+ (dry-run) %s -> %s', src, dest)
            return True
        if not os.path.exists(src):
            return False
        if suffix.endswith('.gz'):
            os.replace(src, dest)                          # already gzipped
        else:
            run_cmd(F'{args.gzip_cmd} -c "{src}" > "{dest}"')
            os.remove(src)
        return True

    if collect(F'{run}_2{suffix}', fq2):
        collect(F'{run}_1{suffix}', fq1)
        layout = 'PAIRED'
        if os.path.exists(os.path.join(tmp, F'{run}{suffix}')):
            log.info('%s: discarding --split-3 orphan reads', run)
    else:
        # Single-end, or a run annotated PAIRED whose mates were never submitted.
        if not collect(F'{run}_1{suffix}', fq1):
            collect(F'{run}{suffix}', fq1)
        fq2, layout = '', 'SINGLE'
    if not args.dry_run:
        shutil.rmtree(tmp, ignore_errors=True)
        mark_done(fq1)
    return fq1, fq2, layout


def ena_filereport_any(accessions, fields=ENA_FIELDS):
    """Try each accession in turn and return the first that resolves.

    Every study carries both a BioProject (PRJNA/PRJEB) and a secondary SRA study
    (SRP/ERP) accession -- e.g. PRJNA629885 and SRP259526 are the same Minussi 2021 data.
    The ENA portal does not always index both, and legacy SRA* submission accessions
    often resolve under neither, so keeping the alternate around avoids a dead end.
    """
    tried = [a for a in accessions if a and str(a).upper() != 'NA']
    if not tried:
        raise SystemExit('No accession to look up: pass --accession.')
    problems = []
    for pos, acc in enumerate(tried):
        try:
            return ena_filereport(acc, fields)
        # ena_filereport signals "no such study" with SystemExit, which is a
        # BaseException and would sail straight past a bare `except Exception`.
        except (Exception, SystemExit) as exc:                     # noqa: BLE001
            problems.append(F'{acc}: {exc}')
            log.warning('ENA lookup failed for %s (%s)', acc, exc)
            if pos < len(tried) - 1:
                log.info('retrying with the alternate accession %s', tried[pos + 1])
    raise SystemExit('No ENA study resolved. Tried -> ' + '; '.join(problems))


def shard_dir(args, fqroot, accession, study):
    """Directory holding <accession>_{1,2}.fastq.gz, created on demand."""
    if args.fastq_layout == 'flat':
        outdir = fqroot
    else:
        outdir = os.path.join(fqroot, cm.fastq_shard(
            accession, study, args.fastq_shard_template, args.fastq_shard_prefix_len))
    if not args.dry_run:
        os.makedirs(outdir, exist_ok=True)
    return outdir


# ------------------------------------------------------------------------------- modes

def mode_ena(args, fqdir):
    """One run == one cell. Rename ENA FASTQs into <accession>_{1,2}.fastq.gz."""
    rows = ena_filereport_any([args.accession, args.alt_accession])
    if args.max_cells:
        rows = rows[:args.max_cells]
    def fetch_one(row):
        run = row['run_accession']
        ftp = [u for u in (row.get('fastq_ftp') or '').split(';') if u]
        md5 = [m for m in (row.get('fastq_md5') or '').split(';') if m]
        md5 += ['NA'] * (len(ftp) - len(md5))
        if not ftp:
            log.warning('%s has no fastq_ftp on ENA; skipping (try --mode sra_file)', run)
            return None
        # ENA emits [_1,_2] for paired runs, sometimes with a leading unpaired file.
        if len(ftp) >= 2:
            pairs = [(u, m) for u, m in zip(ftp, md5)
                     if u.endswith(('_1.fastq.gz', '_2.fastq.gz'))] or list(zip(ftp, md5))[:2]
            layout = 'PAIRED'
        else:
            pairs, layout = list(zip(ftp, md5)), 'SINGLE'
        if layout == 'PAIRED' and len(pairs) < 2:
            layout = 'SINGLE'          # declared paired, only one file actually offered
        study = row.get('study_accession') or args.accession
        celldir = shard_dir(args, fqdir, run, study)
        fqpaths = ['', '']
        for idx, (url, m) in enumerate(pairs, start=1):
            fqpaths[idx - 1] = download(url, os.path.join(celldir, F'{run}_{idx}.fastq.gz'),
                                        md5=m, dry_run=args.dry_run,
                                        downloader=args.downloader, conns=args.lftp_conns)
        return make_meta_row(
            fastq1=fqpaths[0], fastq2=fqpaths[1],
            run=run, library=row.get('library_name') or run,
            sample=row.get('sample_title') or row.get('sample_alias') or run,
            study=row.get('study_accession') or args.accession,
            donor=args.donor, sample_type=args.sample_type,
            platform=(row.get('instrument_platform') or 'ILLUMINA').upper(),
            layout=layout, bases=row.get('base_count', ''),
            avg_spot_len=avg_spot_len(row))

    results, failures = run_parallel(fetch_one, rows, args.jobs, 'run')
    args.failures = report_failures(failures, args.outdir, args.dataset_id, len(rows))
    return [m for m in results if m]


def avg_spot_len(row):
    try:
        reads, bases = int(row.get('read_count') or 0), int(row.get('base_count') or 0)
        return str(bases // reads) if reads else '0'
    except (TypeError, ValueError):
        return '0'


def _runs_and_ena_rows(args):
    """Run list plus whatever ENA knows about each run (may be empty)."""
    ena_rows = {}
    if args.accession:
        try:
            ena_rows = {r['run_accession']: r
                        for r in ena_filereport_any([args.accession, args.alt_accession])}
        except (Exception, SystemExit) as exc:                      # noqa: BLE001
            log.warning('ENA lookup failed for %s (%s); continuing without its metadata',
                        args.accession, exc)
    runs = list(args.runs or ena_rows)
    if not runs:
        raise SystemExit(F'No runs for --mode {args.mode}: pass --runs, or an --accession '
                         'that the ENA portal can resolve.')
    return (runs[:args.max_cells] if args.max_cells else runs), ena_rows


def mode_sra(args, fqdir):
    """prefetch + dump fallback for studies ENA does not mirror as FASTQ."""
    require_tool('prefetch', 'Install sra-tools (conda install -c bioconda sra-tools).')
    runs, ena_rows = _runs_and_ena_rows(args)
    def fetch_one(run):
        row = ena_rows.get(run, {})
        study = row.get('study_accession') or args.accession
        celldir = shard_dir(args, fqdir, run, study)
        # prefetch lays runs out as <outdir>/<run>/<run>.sra.
        cached = os.path.join(args.sradir, run, F'{run}.sra')
        if not os.path.exists(cached):
            run_cmd(['prefetch', '--max-size', 'u', '-O', args.sradir, run], args.dry_run)
        source = cached if os.path.exists(cached) else run
        fq1, fq2, layout = dump_sra_to_fastq(args, source, run, celldir)
        drop_sra(args, cached)
        return make_meta_row(run=run, library=row.get('library_name') or run,
                                  sample=row.get('sample_title') or run, study=study,
                                  donor=args.donor, sample_type=args.sample_type,
                                  platform=(row.get('instrument_platform') or 'ILLUMINA').upper(),
                                  layout=layout, fastq1=fq1, fastq2=fq2,
                                  bases=row.get('base_count', ''),
                                  avg_spot_len=avg_spot_len(row) if row else '0')

    results, failures = run_parallel(fetch_one, runs, args.jobs, 'run')
    args.failures = report_failures(failures, args.outdir, args.dataset_id, len(runs))
    return [m for m in results if m]


def drop_sra(args, path):
    """Delete a converted .sra unless --keep-sra. These are the bulk of the disk use."""
    if args.keep_sra or args.dry_run or not path or not os.path.exists(path):
        return
    os.remove(path)
    if os.path.exists(path + '.done'):
        os.remove(path + '.done')
    parent = os.path.dirname(path)          # prefetch's per-run directory, if now empty
    if parent != args.sradir and os.path.isdir(parent) and not os.listdir(parent):
        os.rmdir(parent)


def mode_sra_file(args, fqdir):
    """Download each run's raw .sra archive, then dump it locally.

    Same output as --mode sra, but the transfer goes through lftp `pget -c -n`, which is
    parallel, resumable and md5-checkable, instead of prefetch's single stream.
    """
    runs, ena_rows = _runs_and_ena_rows(args)

    def fetch_one(run):
        row = ena_rows.get(run, {})
        study = row.get('study_accession') or args.accession
        celldir = shard_dir(args, fqdir, run, study)
        sra = os.path.join(args.sradir, F'{run}.sra')
        if is_done(sra) or os.path.exists(sra):
            log.info('reusing local archive: %s', sra)
        else:
            md5s = [m for m in (row.get('sra_md5') or '').split(';') if m.strip()]
            candidates = sra_url_candidates(run, row, args.sra_url_template)
            for pos, url in enumerate(candidates):
                # ENA-derived URLs are the ones the md5 belongs to; templates are not.
                md5 = md5s[0] if (md5s and url in (row.get('sra_ftp') or '')) else None
                # No HEAD pre-check: some mirrors reject HEAD, and a failed transfer is
                # just as cheap to detect from the downloader's own exit status.
                try:
                    download(url, sra, md5=md5, dry_run=args.dry_run,
                             downloader=args.downloader, conns=args.lftp_conns)
                    break
                except Exception as exc:                           # noqa: BLE001
                    log.warning('failed from %s (%s)', url, exc)
                    if pos == len(candidates) - 1:
                        raise SystemExit(
                            F'Could not download the .sra for {run}. Tried: '
                            + ', '.join(candidates)
                            + '. Pass --sra-url-template with a <run> placeholder for a '
                              'mirror that works from this host, or use --mode sra.')
            else:
                raise SystemExit(F'No reachable .sra URL for {run}; tried '
                                 + ', '.join(candidates))
        fq1, fq2, layout = dump_sra_to_fastq(args, sra, run, celldir)
        drop_sra(args, sra)
        return make_meta_row(run=run, library=row.get('library_name') or run,
                             sample=row.get('sample_title') or run, study=study,
                             donor=args.donor, sample_type=args.sample_type,
                             platform=(row.get('instrument_platform') or 'ILLUMINA').upper(),
                             layout=layout, fastq1=fq1, fastq2=fq2,
                             bases=row.get('base_count', ''),
                             avg_spot_len=avg_spot_len(row) if row else '0')

    results, failures = run_parallel(fetch_one, runs, args.jobs, 'run')
    args.failures = report_failures(failures, args.outdir, args.dataset_id, len(runs))
    return [m for m in results if m]


def read_barcode_whitelist(path):
    """Accept a bare barcode list or a 10x per_cell_summary_metrics.csv."""
    barcodes = []
    with open(path) as fh:
        sniff = fh.readline()
        fh.seek(0)
        if ',' in sniff and 'barcode' in sniff.lower():
            for row in csv.DictReader(fh):
                key = next((k for k in row if k and k.lower().strip() == 'barcode'), None)
                if key:
                    barcodes.append(row[key].strip())
        else:
            barcodes = [ln.strip() for ln in fh if ln.strip() and not ln.startswith('#')]
    log.info('read %d barcodes from %s', len(barcodes), path)
    return barcodes


def mode_tenx_bam(args, fqdir):
    """Split a 10x position-sorted BAM on CB, then emit one FASTQ pair per cell."""
    require_tool('samtools', 'Install samtools >= 1.14.')
    if not samtools_supports_split_by_tag():
        raise SystemExit('samtools split -d TAG is unavailable (needs samtools >= 1.14). '
                         'Upgrade samtools, or pre-split the BAM and use --mode bam.')
    bam = args.bam
    if not bam:
        if not args.bam_url:
            raise SystemExit('--mode tenx_bam needs --bam (local) or --bam-url (remote).')
        if not args.dry_run and not url_exists(args.bam_url):
            raise SystemExit(
                F'BAM URL is not reachable: {args.bam_url}\n'
                'The 10x CDN layout changes between releases. Open the dataset landing '
                'page, copy the "Possorted BAM" link, and pass it via --bam-url.')
        bam = download(args.bam_url,
                       os.path.join(args.bamdir, os.path.basename(
                           urllib.parse.urlparse(args.bam_url).path)),
                       dry_run=args.dry_run, downloader=args.downloader,
                       conns=args.lftp_conns)

    splitdir = os.path.join(args.tmpdir, F'{args.dataset_id}.cellbams')
    if not args.dry_run: os.makedirs(splitdir, exist_ok=True)
    sentinel = os.path.join(splitdir, 'split')
    if is_done(sentinel):
        log.info('skip CB split (already done): %s', splitdir)
    else:
        # One pass over the BAM; samtools manages the per-tag output handles itself.
        # -d defaults to --max-split 100; a 10x run has thousands of barcodes, so the
        # default would silently drop most cells. samtools holds one file handle per
        # output, so raise `ulimit -n` above --max-split before running this.
        run_cmd(['samtools', 'split', '-d', 'CB', '-M', str(args.max_split),
                 '-@', str(args.threads),
                 '-u', os.path.join(splitdir, 'unassigned.bam'),
                 '-f', os.path.join(splitdir, '%!.bam'), bam], args.dry_run)
        mark_done(sentinel, args.dry_run)

    whitelist = set(read_barcode_whitelist(args.barcodes)) if args.barcodes else None
    cellbams = sorted(f for f in (os.listdir(splitdir) if os.path.isdir(splitdir) else [])
                      if f.endswith('.bam') and f != 'unassigned.bam')
    log.info('%d barcode BAMs in %s', len(cellbams), splitdir)

    selected = [f for f in cellbams if whitelist is None or f[:-4] in whitelist]
    if args.max_cells:
        # Bounded before conversion so a pilot run does not convert every barcode first.
        selected = selected[:args.max_cells * 4]

    def convert_one(fname):
        barcode = fname[:-4]
        cellbam = os.path.join(splitdir, fname)
        # Background barcodes carry a handful of reads and just add noise + scheduler load.
        if not args.dry_run and args.min_reads:
            n = int(subprocess.run(['samtools', 'view', '-c', cellbam],
                                   capture_output=True, text=True, check=True).stdout.strip())
            if n < args.min_reads:
                return None
        acc = F'{args.cell_prefix}{barcode}'.replace('.', '_')   # '.' breaks samplename parsing
        celldir = shard_dir(args, fqdir, acc, args.accession)
        fq1 = os.path.join(celldir, F'{acc}_1.fastq.gz')
        fq2 = os.path.join(celldir, F'{acc}_2.fastq.gz')
        if is_done(fq1):
            log.info('skip (already done): %s', fq1)
        else:
            # collate is mandatory: samtools fastq needs mates adjacent, and a 10x BAM
            # is position-sorted. -n keeps the original read names.
            run_cmd(F'samtools collate -u -O -@ {args.threads} "{cellbam}" '
                    F'{args.tmpdir}/collate.{acc} | '
                    F'samtools fastq -n -1 "{fq1}" -2 "{fq2}" -0 /dev/null -s /dev/null -',
                    args.dry_run)
            mark_done(fq1, args.dry_run)
        return make_meta_row(run=acc, library=args.dataset_id, sample=acc,
                             study=args.accession, donor=args.donor,
                             sample_type=args.sample_type, platform='ILLUMINA',
                             layout='PAIRED', fastq1=fq1, fastq2=fq2)

    # samtools already gets -@ threads, so keep this pool modest relative to --jobs.
    results, failures = run_parallel(convert_one, selected,
                                     max(1, min(args.jobs, args.convert_jobs)), 'cell')
    args.failures = report_failures(failures, args.outdir, args.dataset_id, len(selected))
    meta = [m for m in results if m]
    if args.max_cells:
        meta = meta[:args.max_cells]
    log.info('kept %d cells (min_reads=%s, whitelist=%s)', len(meta), args.min_reads,
             bool(whitelist))
    return meta


def mode_bam(args, fqdir):
    """A single BAM -> a single cell."""
    require_tool('samtools')
    if not args.bam:
        raise SystemExit('--mode bam needs --bam.')
    acc = args.cell_prefix or os.path.basename(args.bam).split('.')[0]
    celldir = shard_dir(args, fqdir, acc, args.accession)
    fq1, fq2 = (os.path.join(celldir, F'{acc}_{i}.fastq.gz') for i in (1, 2))
    if is_done(fq1):
        log.info('skip (already done): %s', fq1)
    else:
        run_cmd(F'samtools collate -u -O -@ {args.threads} "{args.bam}" '
                F'{args.tmpdir}/collate.{acc} | '
                F'samtools fastq -n -1 "{fq1}" -2 "{fq2}" -0 /dev/null -s /dev/null -',
                args.dry_run)
        mark_done(fq1, args.dry_run)
    return [make_meta_row(run=acc, library=acc, sample=acc, study=args.accession,
                          donor=args.donor, sample_type=args.sample_type,
                          platform='ILLUMINA', layout='PAIRED', fastq1=fq1, fastq2=fq2)]


# ---------------------------------------------------------------------------- metadata

def make_meta_row(run, library, sample, study, donor, sample_type, platform, layout,
                  bases='', avg_spot_len='0', center='NA', fastq1='', fastq2=''):
    """One row of the SraRunTable that data_tumor.py consumes.

    data_tumor derives Donor from the 'isolate' column (there is no 'tissue' column
    here) and then prefixes it with 'SRA~Study', so the effective donor directory
    becomes '<study>_<donor>'. Values must be space-free: data_tumor replaces ' '
    with '-' anyway, but doing it here keeps the TSV readable.
    """
    def clean(v):
        return str(v if v not in (None, '') else 'NA').replace(' ', '-').replace('\t', '-')
    return {'Run': clean(run), 'AvgSpotLen': avg_spot_len or '0',
            'Library Name': clean(library), 'Sample Name': clean(sample),
            'sample_type': clean(sample_type), 'SRA Study': clean(study),
            'Bases': bases or '', 'Bytes': '', 'SEQUENCING_CENTER (run)': clean(center),
            'isolate': clean(donor), 'Platform': clean(platform),
            'LibraryLayout': clean(layout),
            'Fastq1': (os.path.abspath(fastq1) if fastq1 else ''),
            'Fastq2': (os.path.abspath(fastq2) if fastq2 else '')}


def verify_layouts(rows, dry_run=False):
    """Trust the bytes on disk over the submitter's LibraryLayout annotation.

    ENA/SRA routinely annotate a run PAIRED when only single-end reads were submitted;
    downstream that becomes a missing <accession>_2.fastq.gz at bwa time. Anything we
    just wrote is on disk right now, so settle the question here instead.
    """
    if dry_run:
        return rows
    fixed = 0
    for row in rows:
        if row.get('LibraryLayout') != 'PAIRED':
            continue
        if cm.fastq_has_reads(row.get('Fastq2')):
            continue
        log.warning('%s is annotated PAIRED but read 2 has no reads; recording SINGLE',
                    row.get('Run'))
        row['LibraryLayout'], row['Fastq2'] = 'SINGLE', ''
        fixed += 1
    if fixed:
        log.warning('corrected LibraryLayout for %d/%d runs', fixed, len(rows))
    return rows


def write_metadata_tsv(rows, path, dry_run=False):
    if not rows:
        if dry_run:
            log.warning('dry-run produced no rows (expected: nothing was downloaded)')
            return path
        raise SystemExit('No cells were produced; refusing to write an empty run table.')
    log.info('writing %d rows -> %s', len(rows), path)
    if dry_run:
        return path
    os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
    with open(path, 'w', newline='') as fh:
        # lineterminator='\n': the csv default is '\r\n', which would leave a stray
        # '\r' on the last column (LibraryLayout) and break data_tumor's PAIRED assert.
        writer = csv.DictWriter(fh, fieldnames=META_COLUMNS, delimiter='\t',
                                lineterminator='\n', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    return path


def load_manifest(path):
    with open(path) as fh:
        header = fh.readline().lstrip('#').rstrip('\n').split('\t')
        return {r['dataset_id']: r
                for r in csv.DictReader(fh, fieldnames=header, delimiter='\t')
                if r.get('dataset_id') and not r['dataset_id'].startswith('#')}


# -------------------------------------------------------------------------------- main

MODES = {'ena': mode_ena, 'sra': mode_sra, 'sra_file': mode_sra_file,
         'tenx_bam': mode_tenx_bam, 'bam': mode_bam}


def build_parser(script_dir):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--manifest', default=os.path.join(script_dir, 'public_scwgs_datasets.tsv'))
    p.add_argument('--dataset', help='dataset_id from --manifest (fills the options below)')
    p.add_argument('--list', action='store_true', help='list datasets in the manifest and exit')
    p.add_argument('--fill-alt-accessions', action='store_true', help=(
        'ask ENA for each manifest row\'s other study accession (BioProject vs SRA '
        'study), write the pairs back into the alt_accession column and exit. Rows whose '
        'accession is not an INSDC study id are marked "none". Combine with --dry-run to '
        'preview, or --force to overwrite entries that are already filled.'))
    p.add_argument('--force', action='store_true',
                   help='with --fill-alt-accessions, re-resolve rows that are already filled')
    p.add_argument('--inspect', action='store_true', help=(
        'summarise --accession from the ENA portal (runs, platform, species, annotated '
        'vs actual layout, size) and print a candidate manifest row, then exit. '
        'Downloads nothing: use it to vet an accession before adding it to the TSV.'))
    p.add_argument('--mode', choices=sorted(MODES))
    p.add_argument('--accession', help='ENA/SRA study accession, e.g. SRP114962')
    p.add_argument('--alt-accession', default=None, help=(
        'the same study under its other accession form (BioProject vs SRA study), tried '
        'when the primary is not indexed by ENA'))
    p.add_argument('--runs', nargs='+', help='explicit run accessions (overrides --accession lookup)')
    p.add_argument('--bam', help='local BAM (tenx_bam / bam modes)')
    p.add_argument('--bam-url', help='remote BAM URL (tenx_bam mode)')
    p.add_argument('--barcodes', help='barcode whitelist or 10x per_cell_summary_metrics.csv')
    p.add_argument('--max-split', type=int, default=20000,
                   help='samtools split --max-split; must exceed the barcode count, and '
                        'ulimit -n must exceed it too')
    p.add_argument('--min-reads', type=int, default=20000,
                   help='drop per-barcode BAMs below this read count (0 disables)')
    p.add_argument('--cell-prefix', default='', help='prefix for synthesised per-cell accessions')
    p.add_argument('--donor', default='', help='donor label; defaults to the dataset_id')
    p.add_argument('--sample-type', default='tumor')
    p.add_argument('--outdir', default=os.path.abspath(
        os.path.join(script_dir, '..', 'real_tumor_data')),
        help='must match data0to1dir / --tumor-datdir used by main.py')
    p.add_argument('--tmpdir', default=None)
    p.add_argument('--bamdir', default=None,
                   help='where downloaded BAMs land (default <outdir>/bam.datdir)')
    p.add_argument('--sradir', default=None,
                   help='where .sra archives land (default <outdir>/sra.datdir)')
    p.add_argument('--sra-url-template', default=None, help=(
        'tried before ENA and the built-in mirrors; <run> is replaced with the run '
        'accession, e.g. https://my.mirror/sra/<run>/<run>.sra'))
    p.add_argument('--dump-tool', choices=['fasterq-dump', 'fastq-dump'],
                   default='fasterq-dump', help=(
        'fasterq-dump is multi-threaded but writes plain FASTQ that must then be gzipped; '
        'fastq-dump is slower yet gzips inline, which roughly halves the disk peak'))
    p.add_argument('--dump-args', nargs='*', default=[],
                   help='extra flags passed through to the dump tool, e.g. --skip-technical')
    p.add_argument('--keep-sra', action='store_true',
                   help='keep the .sra archives after conversion (they dominate disk use)')
    p.add_argument('--fastq-layout', choices=cm.FASTQ_LAYOUTS, default='sharded', help=(
        'sharded (default) writes <outdir>/1from0.datdir/<shard>/<accession>_1.fastq.gz; '
        'flat writes the historical <outdir>/1from0.datdir/<accession>_1.fastq.gz'))
    p.add_argument('--fastq-shard-template', default=cm.FASTQ_SHARD_TEMPLATE,
                   help='shard sub-path from <study>, <accprefix>, <accession>')
    p.add_argument('--fastq-shard-prefix-len', type=int, default=cm.FASTQ_SHARD_PREFIX_LEN,
                   help='leading accession characters forming <accprefix>')
    p.add_argument('--downloader', choices=['lftp', 'python'], default='lftp',
                   help='lftp uses parallel pget with retries; falls back to python if absent')
    p.add_argument('--lftp-conns', type=int, default=1000-1, help='lftp pget -n')
    p.add_argument('-j', '--jobs', type=int, default=16, help=(
        'how many runs to fetch concurrently: parallel prefetch/.sra downloads in the sra '
        'and sra_file modes, and parallel FASTQ downloads from ENA in the ena mode. '
        'Downloads are I/O bound, so this is worth far more than --threads. Use 1 for '
        'strictly serial behaviour'))
    p.add_argument('--convert-jobs', type=int, default=8, help=(
        'concurrent per-cell BAM->FASTQ conversions in tenx_bam mode; capped by --jobs. '
        'Each conversion already gets --threads samtools threads'))
    p.add_argument('--threads', type=int, default=8,
                   help='threads *per* subprocess (samtools, fasterq-dump); multiplies with --jobs')
    p.add_argument('--max-cells', type=int, default=0, help='0 = no limit')
    p.add_argument('--gzip-cmd', default='gzip -1', help='e.g. "pigz -p 8" if available')
    p.add_argument('--allow-controlled', action='store_true',
                   help='permit datasets flagged CONTROLLED in the manifest')
    p.add_argument('--dry-run', action='store_true')
    p.add_argument('-v', '--verbose', action='store_true')
    return p


def main(argv=None):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    args = build_parser(script_dir).parse_args(argv)
    logging.basicConfig(level=(logging.DEBUG if args.verbose else logging.INFO),
                        format='%(asctime)s %(levelname)s %(message)s')

    manifest = load_manifest(args.manifest) if os.path.exists(args.manifest) else {}
    if args.list:
        for did, row in manifest.items():
            print(F"{did}\t{row['mode']}\t{row['accession']}\t{row['access']}\t{row['notes']}")
        return 0

    if args.fill_alt_accessions:
        return fill_alt_accessions(args)

    if args.inspect:
        if not args.accession and args.dataset in manifest:
            args.accession = manifest[args.dataset]['accession']
        if not args.accession:
            raise SystemExit('--inspect needs --accession (or a --dataset from the manifest).')
        return inspect_accession(args)

    if args.dataset:
        entry = manifest.get(args.dataset)
        if entry is None:
            raise SystemExit(F'Unknown dataset {args.dataset}. Try --list.')
        if entry.get('access', '').upper() == 'CONTROLLED' and not args.allow_controlled:
            raise SystemExit(
                F"{args.dataset} is controlled-access ({entry['accession']}). Obtain a DAC "
                'approval first, then rerun with --allow-controlled and a credentialed client.')
        args.mode = args.mode or entry['mode']
        args.accession = args.accession or entry['accession']
        args.alt_accession = args.alt_accession or entry.get('alt_accession')
        args.donor = args.donor or (entry.get('donor') or args.dataset)
        args.sample_type = entry.get('sample_type') or args.sample_type
        if entry.get('url', 'NA') not in ('NA', '') and not args.bam_url and not args.bam:
            args.bam_url = entry['url']
        args.dataset_id = args.dataset
    else:
        if not args.mode:
            raise SystemExit('Pass --dataset or --mode.')
        args.dataset_id = args.accession or (args.bam and os.path.basename(args.bam)) or 'adhoc'
    args.donor = args.donor or args.dataset_id
    args.cell_prefix = args.cell_prefix or (F'{args.dataset_id}_' if args.mode == 'tenx_bam' else '')
    args.tmpdir = args.tmpdir or os.path.join(args.outdir, 'fetch.tmpdir')
    args.bamdir = args.bamdir or os.path.join(args.outdir, 'bam.datdir')
    args.sradir = args.sradir or os.path.join(args.outdir, 'sra.datdir')

    args.failures = None
    # 16 jobs x `pget -n 999` is ~16k sockets against one host: ENA and NCBI will throttle
    # or ban that. Warn rather than override, since --lftp-conns is a deliberate choice.
    if args.downloader == 'lftp' and args.jobs * args.lftp_conns > 512:
        log.warning('--jobs %d x --lftp-conns %d = %d concurrent connections; ENA/NCBI will '
                    'throttle or refuse. Consider --lftp-conns %d, or --jobs 1 for one big '
                    'file (a BAM) where all the connections should go to that transfer.',
                    args.jobs, args.lftp_conns, args.jobs * args.lftp_conns,
                    max(1, 512 // max(1, args.jobs)))

    fqdir = os.path.join(args.outdir, '1from0.datdir')
    if not args.dry_run:
        for d in (fqdir, args.tmpdir, args.bamdir, args.sradir): os.makedirs(d, exist_ok=True)

    rows = MODES[args.mode](args, fqdir)
    rows = verify_layouts(rows, args.dry_run)
    tsv = write_metadata_tsv(rows, os.path.join(args.outdir,
                                                F'{args.dataset_id}.SraRunTable.tsv'),
                             args.dry_run)
    shard_hint = ('<study>/<accprefix>/' if args.fastq_layout != 'flat' else '')
    print(F'''
FASTQ   : {fqdir}/{shard_hint}<accession>_{{1,2}}.fastq.gz   ({len(rows)} cells)
Metadata: {tsv}   (carries absolute Fastq1/Fastq2 columns){F'{chr(10)}FAILED  : {args.failures}   (rerun to retry: every success is skipped)' if args.failures else ''}

Next:
  python main.py --tumor-fastq --SraRunTable {tsv} --tumor-datdir {os.path.abspath(args.outdir)} > Snakefile
  snakemake --cores ${{NUM_CPUS}}
'''.rstrip())
    return 0


if __name__ == '__main__':
    sys.exit(main())
