#!/usr/bin/env python
# https://claude.ai/chat/1154badf-2241-4cce-9e52-7c5565f9abe3

"""
Compare the observed per-cell ploidy (inferred by a benchmarked scWGS CNV caller) against the
expected per-cell ploidy (an orthogonal experimental estimate, e.g. FACS/DAPI or karyotype)
supplied through a ploidy file, using the evaluation metrics defined by scAbsolute:

    Schneider MP, Cullen AE, Pangonyte J, et al.
    "scAbsolute: measuring single-cell ploidy and replication status."
    Genome Biology 2024;25:62.  https://doi.org/10.1186/s13059-024-03204-y

Definitions taken from that paper
--------------------------------
* Ploidy (their Eq. 1) is the mean absolute copy number over the genomic bins of a cell,
      p = (1/M) * sum_j c_j ,
  i.e. it is proportional to the amount of DNA in the cell and is *not* a chromosome count.
  Our per-cell CNV calls are variable-length BED segments rather than fixed-size bins, so the
  direct generalisation is the segment-length-weighted mean copy number (identical to Eq. 1
  when the segments are the caller's fixed-size bins).  `--weight segment` gives the
  unweighted per-segment mean instead.

* Primary metric 1 -- percentage of ploidy outliers.  scAbsolute assess "the percentage of
  cells outside an experimental ploidy window of +/- 0.5 around the peak of the DAPI
  distribution"; the window absorbs segmentation and FACS uncertainty while still excluding
  genuine ploidy changes.  Reported here as `pct_outliers` (window size is `--ploidy-window`).

* Primary metric 2 -- mean absolute ploidy distance.  "the mean absolute distance across all
  cells in a sample from the experimental ploidy estimate".  Reported as
  `mean_abs_ploidy_distance`.  In scAbsolute's Table 1 this is the value in parentheses.

* Table 1 aggregate -- "the mean % of outliers across all samples per method", i.e. the
  unweighted mean of `pct_outliers` over samples.  Reported in the ALL row as
  `mean_pct_outliers_across_samples` (the cell-weighted pooled value is also given).

* Auxiliary diagnostics (from the Fig. 4 annotation, where "blue asterisks indicate ploidy
  levels of 1/2 or 2 times the experimental ploidy estimate"): the fraction of cells that
  land within the same +/- window around 2 x p_exp or 0.5 x p_exp.  These are the
  characteristic non-identifiability failures (whole-genome-doubling / G2 confusion and
  halving), and separating them from unstructured error is what makes the outlier rate
  interpretable.  Reported as `pct_near_2x`, `pct_near_half`, `pct_scaling_error`.

* Copy-number cap (`--max-cn`).  A ploidy is a genome-wide mean, so a few focal high-level
  amplifications can move it by more than the whole evaluation window.  Capping each segment
  at `--max-cn` (default 10) keeps such segments counted but bounded; `--max-cn inf` averages
  the copy numbers exactly as the caller reported them.  The two settings answer different
  questions -- how well the caller places the bulk of the genome, versus what its raw output
  literally implies -- so the pipeline runs both and writes them to separate output prefixes.

Usage
-----
    python ploidy_eval.py \
        -i '/path/to/4from2_..._*intcns.bed' \
        -o /path/to/out_prefix \
        --ploidy-file ploidy.PRJNA629885.tsv \
        --metadata-tsv SraRunTable.PRJNA629885.tsv \
        --plot

Outputs `<prefix>_percell.tsv`, `<prefix>_persample.tsv`, `<prefix>_summary.json` and, with
`--plot`, `<prefix>.pdf` / `<prefix>.png` in the style of scAbsolute Fig. 4.

This module is also importable: `load_ploidy_table`, `bed_to_ploidy`, `per_cell_metrics` and
`summarize_by_sample` are the pieces reused by data_tumor.py.
"""
import argparse, collections, glob, json, logging, os, re, sys

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------------------------

AUTOSOMES = [F'chr{i}' for i in range(1, 23)]
SEX_CHROMS = ['chrX', 'chrY']
ALL_CHROMS = AUTOSOMES + SEX_CHROMS

# scAbsolute's window around the experimental point estimate (their Table 1 / Fig. 4).
DEFAULT_PLOIDY_WINDOW = 0.5

# Cap applied to the per-segment copy number before it is averaged into a ploidy. Ploidy is a
# genome-wide mean, so without a cap a few focal high-level amplifications (which some callers
# report as CN in the hundreds) move it by more than the whole +/- 0.5 window; with a cap the
# amplified segments still count, only not unboundedly. `inf` disables the cap, which is the
# right setting when the question is what the caller literally reported.
DEFAULT_MAX_CN = 10.0
UNCAPPED_MAX_CN_ALIASES = ('inf', 'infinity', 'none', 'no', 'nan', '')


def is_uncapped(max_cn):
    """True when `max_cn` asks for the copy numbers to be averaged exactly as they were called."""
    return max_cn is None or not np.isfinite(float(max_cn))


def parse_max_cn(value):
    """--max-cn: a number, or any of `inf`/`none` to leave the copy numbers uncapped."""
    if value is None:
        return float('inf')
    if str(value).strip().lower() in UNCAPPED_MAX_CN_ALIASES:
        return float('inf')
    max_cn = float(value)                    # raises ValueError, which argparse turns into a message
    if not (max_cn > 0):
        raise ValueError(F'--max-cn {value} must be positive')
    return max_cn

# Columns of an SraRunTable that may carry a human-readable sample label, most specific first.
# `TN2_S6_C359`-style library names collapse to `TN2` through the prefix backoff in
# PloidyTable.lookup(), so listing the per-cell columns first costs nothing.
DEFAULT_SAMPLE_KEY_COLUMNS = [
    'Sample~Name', 'Sample Name', 'Sample_Name', 'SampleName', 'sample_name',
    'Library~Name', 'Library Name', 'Library_Name', 'LibraryName', 'library_name',
    'sample_title', 'Title', 'isolate', 'cell_line', 'Cell_Line', 'cell~line',
    'source_name', 'source~name', 'tissue', 'Donor', 'donor', 'submitted_subject_id',
]
RUN_COLUMNS = ['#Run', 'Run', 'run_accession', 'Run~accession']
RUN_ACCESSION_RE = re.compile(r'((?:[SED]RR|SAMN|GSM)\d+)')

# ---------------------------------------------------------------------------------------------
# Ploidy file
# ---------------------------------------------------------------------------------------------


def _norm_key(name):
    """Fold the many spellings of one sample id onto a single key.

    `MDA-MB-231`, `MDA_MB_231`, `MDAMB231` and `mdamb231` all become `mdamb231`, which is what
    lets one ploidy file serve run tables that were curated by different people.
    """
    return re.sub(r'[^a-z0-9]', '', str(name).strip().lower())


class PloidyTable:
    """sample -> expected ploidy, with tolerant name resolution.

    Resolution order for a candidate label (first hit wins):
      1. exact (normalised) sample id or alias
      2. user-supplied regex alias (`re:<pattern>`)
      3. token-prefix backoff: `TN2_S6_C359` -> `TN2_S6` -> `TN2`
      4. token-bounded substring, longest key first (so `MDA-MB-231-EX1` beats `MDA-MB-231`)
    """

    def __init__(self, records):
        self.records = list(records)                 # list of dict rows, in file order
        self.key2sample = {}                         # normalised key -> canonical sample id
        self.sample2ploidy = {}
        self.regexes = []                            # (compiled pattern, sample id)
        for rec in self.records:
            sample = str(rec['sample']).strip()
            self.sample2ploidy[sample] = float(rec['ploidy'])
            keys = [sample] + list(rec.get('alias_list', []))
            for key in keys:
                key = str(key).strip()
                if not key:
                    continue
                if key.startswith('re:'):
                    self.regexes.append((re.compile(key[3:]), sample))
                    continue
                nkey = _norm_key(key)
                if not nkey:
                    continue
                prev = self.key2sample.get(nkey)
                if prev is not None and prev != sample:
                    logging.warning('ploidy file: key %r maps to both %r and %r; keeping %r',
                                    key, prev, sample, prev)
                    continue
                self.key2sample[nkey] = sample
        # Longest keys first so that the substring pass prefers the most specific sample.
        self._keys_by_len = sorted(self.key2sample, key=len, reverse=True)

    def __len__(self):
        return len(self.sample2ploidy)

    def __contains__(self, sample):
        return sample in self.sample2ploidy

    def ploidy_of(self, sample):
        return self.sample2ploidy[sample]

    def _resolve_one(self, candidate):
        if candidate is None:
            return None
        cand = str(candidate).strip()
        if not cand:
            return None
        nkey = _norm_key(cand)
        if nkey in self.key2sample:
            return self.key2sample[nkey]
        for pattern, sample in self.regexes:
            if pattern.search(cand):
                return sample
        # Token-prefix backoff: drop trailing `_S6`, `_C359`, ... one token at a time.
        tokens = [t for t in re.split(r'[^A-Za-z0-9]+', cand) if t]
        for stop in range(len(tokens) - 1, 0, -1):
            sub = _norm_key(''.join(tokens[:stop]))
            if sub in self.key2sample:
                return self.key2sample[sub]
        # Token-bounded substring, longest key first.
        for key in self._keys_by_len:
            if len(key) >= 3 and key in nkey:
                return self.key2sample[key]
        return None

    def lookup(self, candidates):
        """Return (sample_id, expected_ploidy, matched_candidate) or (None, None, None)."""
        for cand in candidates:
            sample = self._resolve_one(cand)
            if sample is not None:
                return sample, self.sample2ploidy[sample], cand
        return None, None, None


def load_ploidy_table(path):
    """Read a ploidy TSV.

    Required columns: `sample` and `ploidy` (a `#` in front of the header is tolerated, as is a
    headerless two-column file).  Optional `aliases` holds `|`-separated alternative spellings;
    an alias of the form `re:<regex>` is matched as a regular expression.  Lines starting with
    `#` are comments.  Any other columns are carried through to the outputs untouched.
    """
    rows = []
    header = None
    with open(path) as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.rstrip('\n').rstrip('\r')
            if not line.strip():
                continue
            if line.startswith('#'):
                if header is None:                       # `#sample<TAB>ploidy<TAB>...`
                    cand = [c.strip() for c in line.lstrip('#').split('\t')]
                    if len(cand) >= 2 and _norm_key(cand[0]) in ('sample', 'sampleid', 'name'):
                        header = cand
                continue
            toks = [t.strip() for t in line.split('\t')]
            if header is None:
                if len(toks) >= 2 and _norm_key(toks[0]) in ('sample', 'sampleid', 'name'):
                    header = toks
                    continue
                header = ['sample', 'ploidy'] + [F'col{i}' for i in range(3, len(toks) + 1)]
            rec = dict(zip(header, toks))
            if 'sample' not in rec or 'ploidy' not in rec:
                raise ValueError(F'{path}:{lineno}: ploidy file needs `sample` and `ploidy` columns')
            if not rec['sample'] or not rec['ploidy']:
                continue
            try:
                rec['ploidy'] = float(rec['ploidy'])
            except ValueError:
                raise ValueError(F'{path}:{lineno}: ploidy {rec["ploidy"]!r} is not a number')
            if not (rec['ploidy'] > 0):
                raise ValueError(F'{path}:{lineno}: ploidy {rec["ploidy"]} must be positive')
            rec['alias_list'] = [a for a in str(rec.get('aliases', '')).split('|') if a.strip()]
            rows.append(rec)
    if not rows:
        raise ValueError(F'{path}: no ploidy records found')
    return PloidyTable(rows)


# ---------------------------------------------------------------------------------------------
# Observed ploidy from per-cell CNV BED files
# ---------------------------------------------------------------------------------------------


def load_bed(path):
    """Read a 4-column (chrom, start, end, copy-number) BED; malformed lines are skipped."""
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line or line.startswith(('#', 'track', 'browser')):
                continue
            toks = line.split('\t')
            if len(toks) < 4:
                continue
            try:
                rows.append((toks[0], int(round(float(toks[1]))),
                             int(round(float(toks[2]))), float(toks[3])))
            except ValueError:
                continue
    return pd.DataFrame(rows, columns=['chrom', 'start', 'end', 'cn'])


def _norm_chrom(chrom):
    c = str(chrom).strip()
    return c if c.lower().startswith('chr') else F'chr{c}'


def bed_to_ploidy(path_or_df, chroms=None, weight='length', max_cn=10):
    """Observed ploidy of one cell = mean absolute copy number (scAbsolute Eq. 1).

    Returns (ploidy, covered_bases, n_segments); ploidy is NaN when nothing is usable.
    `chroms` restricts the calculation (default: autosomes -- sex-chromosome calls are the
    least reliable output of every caller and scAbsolute likewise excludes them when reasoning
    about the diploid baseline).
    `max_cn` caps the per-segment copy number before averaging, so that a handful of focal
    high-level amplifications cannot dominate a mean that is meant to describe the whole
    genome; pass None or a non-finite value to average the copy numbers as they were called.
    """
    df = load_bed(path_or_df) if isinstance(path_or_df, str) else path_or_df.copy()
    if df.empty:
        return float('nan'), 0, 0
    df = df.assign(chrom=df['chrom'].map(_norm_chrom))
    if chroms:
        df = df[df['chrom'].isin(set(chroms))]
    df = df[np.isfinite(df['cn'])]
    if df.empty:
        return float('nan'), 0, 0
    lengths = (df['end'] - df['start']).astype(float).clip(lower=0)
    if weight == 'segment':
        w = pd.Series(np.ones(len(df)), index=df.index)
    else:
        w = lengths
    total = float(w.sum())
    if total <= 0:
        return float('nan'), 0, int(len(df))
    cn_series = df['cn'] if is_uncapped(max_cn) else df['cn'].clip(upper=float(max_cn))
    return float((cn_series * w).sum() / total), int(lengths.sum()), int(len(df))


# ---------------------------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------------------------


def per_cell_metrics(observed, expected, window=DEFAULT_PLOIDY_WINDOW):
    """scAbsolute per-cell quantities for one cell."""
    err = observed - expected
    abs_err = abs(err)
    return {
        'expected_ploidy': expected,
        'observed_ploidy': observed,
        'ploidy_error': err,
        'abs_ploidy_distance': abs_err,
        'obs2exp_ploidy_ratio': (observed / expected) if expected else float('nan'),
        'is_outlier': bool(abs_err > window),
        'near_2x': bool(abs(observed - 2.0 * expected) <= window),
        'near_half': bool(abs(observed - 0.5 * expected) <= window),
    }


def _pct(mask):
    return float(100.0 * np.mean(mask)) if len(mask) else float('nan')


def summarize_by_sample(percell_df, window=DEFAULT_PLOIDY_WINDOW):
    """Per-sample table in the layout of scAbsolute Table 1, plus an ALL row.

    The ALL row carries both `mean_pct_outliers_across_samples` (scAbsolute's own last row:
    an unweighted mean over samples, so a 40-cell sample counts as much as a 1300-cell one)
    and the cell-weighted pooled figures.
    """
    rows = []
    for sample, sub in percell_df.groupby('sample', sort=True):
        ok = sub[np.isfinite(sub['observed_ploidy'])]
        if ok.empty:
            continue
        abs_dist = ok['abs_ploidy_distance'].to_numpy(dtype=float)
        err = ok['ploidy_error'].to_numpy(dtype=float)
        rows.append({
            'sample': sample,
            'n_cells': int(len(ok)),
            'expected_ploidy': float(ok['expected_ploidy'].iloc[0]),
            # --- scAbsolute primary metrics ---
            'pct_outliers': _pct(ok['is_outlier'].to_numpy(dtype=bool)),
            'mean_abs_ploidy_distance': float(np.mean(abs_dist)),
            # --- auxiliary ---
            'median_abs_ploidy_distance': float(np.median(abs_dist)),
            'rmse_ploidy': float(np.sqrt(np.mean(err ** 2))),
            'mean_signed_ploidy_error': float(np.mean(err)),
            'pct_within_window': _pct(~ok['is_outlier'].to_numpy(dtype=bool)),
            'pct_near_2x': _pct(ok['near_2x'].to_numpy(dtype=bool)),
            'pct_near_half': _pct(ok['near_half'].to_numpy(dtype=bool)),
            'pct_scaling_error': _pct((ok['near_2x'] | ok['near_half']).to_numpy(dtype=bool)),
            'mean_observed_ploidy': float(np.mean(ok['observed_ploidy'])),
            'median_observed_ploidy': float(np.median(ok['observed_ploidy'])),
            'sd_observed_ploidy': float(np.std(ok['observed_ploidy'], ddof=1)) if len(ok) > 1 else 0.0,
        })
    per_sample = pd.DataFrame(rows)
    if per_sample.empty:
        return per_sample, {}
    ok_all = percell_df[np.isfinite(percell_df['observed_ploidy'])]
    abs_all = ok_all['abs_ploidy_distance'].to_numpy(dtype=float)
    err_all = ok_all['ploidy_error'].to_numpy(dtype=float)
    overall = {
        'n_samples': int(len(per_sample)),
        'n_cells': int(len(ok_all)),
        'ploidy_window': float(window),
        # scAbsolute Table 1, last row.
        'mean_pct_outliers_across_samples': float(per_sample['pct_outliers'].mean()),
        'mean_abs_ploidy_distance_across_samples': float(per_sample['mean_abs_ploidy_distance'].mean()),
        # Cell-weighted equivalents.
        'pooled_pct_outliers': _pct(ok_all['is_outlier'].to_numpy(dtype=bool)),
        'pooled_mean_abs_ploidy_distance': float(np.mean(abs_all)) if len(abs_all) else float('nan'),
        'pooled_rmse_ploidy': float(np.sqrt(np.mean(err_all ** 2))) if len(err_all) else float('nan'),
        'pooled_mean_signed_ploidy_error': float(np.mean(err_all)) if len(err_all) else float('nan'),
        'pooled_pct_near_2x': _pct(ok_all['near_2x'].to_numpy(dtype=bool)),
        'pooled_pct_near_half': _pct(ok_all['near_half'].to_numpy(dtype=bool)),
        'pooled_pct_scaling_error': _pct((ok_all['near_2x'] | ok_all['near_half']).to_numpy(dtype=bool)),
    }
    return per_sample, overall


# ---------------------------------------------------------------------------------------------
# Cell -> sample resolution
# ---------------------------------------------------------------------------------------------


def load_run_labels(metadata_tsv, key_columns=None):
    """run accession -> [candidate sample labels], most specific first."""
    if not metadata_tsv or not os.path.isfile(metadata_tsv):
        return {}
    sep = ',' if metadata_tsv.lower().endswith('.csv') else '\t'
    meta = pd.read_csv(metadata_tsv, sep=sep, dtype=str).fillna('')
    meta.columns = [c.strip() for c in meta.columns]
    run_col = next((c for c in RUN_COLUMNS if c in meta.columns), None)
    if run_col is None:
        # Tolerate the `~`-for-space rewriting that data_tumor.py applies in memory.
        alt = {c.replace('~', ' '): c for c in meta.columns}
        run_col = next((alt[c] for c in RUN_COLUMNS if c in alt), None)
    if run_col is None:
        logging.warning('%s: no run column among %s; run->sample mapping disabled',
                        metadata_tsv, RUN_COLUMNS)
        return {}
    cols = list(key_columns) if key_columns else DEFAULT_SAMPLE_KEY_COLUMNS
    present, seen = [], set()
    for c in cols:
        for actual in (c, c.replace('~', ' '), c.replace(' ', '~')):
            if actual in meta.columns and actual not in seen:
                present.append(actual)
                seen.add(actual)
                break
    out = {}
    for _, row in meta.iterrows():
        run = str(row[run_col]).strip()
        if not run:
            continue
        labels = [str(row[c]).strip() for c in present if str(row[c]).strip()]
        out[run] = labels
    return out


def cell_candidates(bed_path, run2labels, sample_regex=None):
    """Ordered list of labels to try when resolving one per-cell BED to a ploidy-file sample."""
    base = re.sub(r'\.bed$', '', os.path.basename(bed_path))
    base = re.sub(r'_(int|dep)cns$', '', base)
    cands = []
    if sample_regex:
        m = re.search(sample_regex, base)
        if m:
            cands.append(m.group(1) if m.groups() else m.group(0))
    runs = RUN_ACCESSION_RE.findall(base)
    for run in runs:
        cands.extend(run2labels.get(run, []))
        cands.append(run)
    cands.append(base)
    # De-duplicate, keeping order.
    seen, out = set(), []
    for c in cands:
        if c and c not in seen:
            seen.add(c)
            out.append(c)
    return out, (runs[0] if runs else '')


# ---------------------------------------------------------------------------------------------
# Plot (scAbsolute Fig. 4 style)
# ---------------------------------------------------------------------------------------------


def plot_ploidy(percell_df, per_sample_df, out_prefix, window, title=''):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    samples = list(per_sample_df['sample'])
    fig, ax = plt.subplots(figsize=(max(5.0, 1.35 * len(samples) + 2.0), 4.6))
    rng = np.random.default_rng(0)
    for i, sample in enumerate(samples):
        exp = float(per_sample_df.loc[per_sample_df['sample'] == sample, 'expected_ploidy'].iloc[0])
        obs = percell_df.loc[percell_df['sample'] == sample, 'observed_ploidy'].to_numpy(dtype=float)
        obs = obs[np.isfinite(obs)]
        ax.add_patch(plt.Rectangle((i - 0.42, exp - window), 0.84, 2 * window,
                                   facecolor='0.82', edgecolor='none', zorder=1))
        ax.scatter(i + rng.uniform(-0.28, 0.28, size=len(obs)), obs,
                   s=5, alpha=0.45, color='#31688e', linewidths=0, zorder=2)
        ax.scatter([i], [exp], marker='x', s=90, color='crimson', linewidths=2.0, zorder=4)
        for mult in (0.5, 2.0):
            ax.scatter([i], [mult * exp], marker='*', s=55, color='#1f77b4', zorder=3)
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=30, ha='right', fontsize=8)
    ax.set_ylabel('ploidy (mean absolute copy number)')
    ax.set_xlim(-0.6, len(samples) - 0.4)
    ax.grid(axis='y', linestyle=':', linewidth=0.6, alpha=0.6)
    if title:
        ax.set_title(title, fontsize=9)
    fig.tight_layout()
    fig.savefig(out_prefix + '.pdf', bbox_inches='tight')
    fig.savefig(out_prefix + '.png', dpi=150, bbox_inches='tight')
    plt.close(fig)


# ---------------------------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------------------------


def build_parser():
    p = argparse.ArgumentParser(
        description=('Compare caller-inferred per-cell ploidy with the expected ploidy from a '
                     'ploidy file, using the scAbsolute metrics (Genome Biol 2024;25:62).'),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('-i', '--input', nargs='+', required=True,
                   help='Per-cell integer-CN BED files or globs (typically *_intcns.bed)')
    p.add_argument('-o', '--output-prefix', required=True)
    p.add_argument('--ploidy-file', required=True,
                   help='TSV mapping sample -> expected ploidy (columns: sample, ploidy, [aliases])')
    p.add_argument('--metadata-tsv', default='',
                   help='SraRunTable used to map run accessions in the BED names to sample labels')
    p.add_argument('--sample-key-columns', nargs='+', default=None,
                   help=F'Metadata columns holding sample labels (default: {DEFAULT_SAMPLE_KEY_COLUMNS[:4]} ...)')
    p.add_argument('--sample-regex', default='',
                   help='Regex applied to the BED basename to extract the sample label directly')
    p.add_argument('--ploidy-window', type=float, default=DEFAULT_PLOIDY_WINDOW,
                   help='Half-width of the experimental ploidy window; a cell outside it is an outlier')
    p.add_argument('--chroms', default='autosomes', choices=['autosomes', 'all'],
                   help='Chromosomes used to compute observed ploidy')
    p.add_argument('--weight', default='length', choices=['length', 'segment'],
                   help='Weighting of the mean copy number: by segment length or per segment')
    p.add_argument('--min-covered-bases', type=float, default=0.0,
                   help='Drop cells whose called segments cover fewer bases than this')
    p.add_argument('--strict', action='store_true',
                   help='Exit non-zero if any cell cannot be resolved to a ploidy-file sample')
    p.add_argument('--max-cn', type=parse_max_cn, default=DEFAULT_MAX_CN, help=(
                   'Cap applied to each segment copy number before averaging; pass inf (or none) '
                   'to average the copy numbers exactly as the caller reported them'))
    p.add_argument('--plot', action='store_true', help='Also write a scAbsolute Fig. 4 style plot')
    p.add_argument('--title', default='')
    # Free-form provenance columns, filled in by data_tumor.py so that per-tool result files
    # from different groups can simply be concatenated.
    p.add_argument('--tool', default='')
    p.add_argument('--donor', default='')
    p.add_argument('--sample-type', default='')
    p.add_argument('--avg-spot-len', default='')
    return p


def main(argv=None):
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(pathname)s:%(lineno)d %(levelname)s - %(message)s')
    args = build_parser().parse_args(argv)

    files = []
    for pat in args.input:
        files.extend(sorted(glob.glob(pat)) if any(c in pat for c in '*?[') else [pat])
    files = [f for f in files if os.path.isfile(f) and os.path.getsize(f) > 0]
    if not files:
        sys.stderr.write('ploidy_eval: no usable input BED files.\n')
        return 1

    table = load_ploidy_table(args.ploidy_file)
    run2labels = load_run_labels(args.metadata_tsv, args.sample_key_columns)
    chroms = AUTOSOMES if args.chroms == 'autosomes' else ALL_CHROMS
    logging.info('ploidy file %s: %d samples; metadata: %d runs; %d BED files',
                 args.ploidy_file, len(table), len(run2labels), len(files))

    records, unresolved = [], []
    for path in files:
        cands, run = cell_candidates(path, run2labels, args.sample_regex or None)
        sample, expected, matched = table.lookup(cands)
        if sample is None:
            unresolved.append((path, cands[:4]))
            continue
        observed, covered, n_seg = bed_to_ploidy(path, chroms=chroms, weight=args.weight, max_cn=args.max_cn)
        if args.min_covered_bases and covered < args.min_covered_bases:
            logging.warning('skipping %s: only %d covered bases', os.path.basename(path), covered)
            continue
        rec = {'cell': re.sub(r'\.bed$', '', os.path.basename(path)), 'run': run,
               'sample': sample, 'matched_label': matched, 'bed': os.path.abspath(path),
               'n_segments': n_seg, 'covered_bases': covered}
        rec.update(per_cell_metrics(observed, expected, args.ploidy_window))
        for k, v in (('tool', args.tool), ('donor', args.donor),
                     ('sampleType', args.sample_type), ('avgSpotLen', args.avg_spot_len)):
            if v:
                rec[k] = v
        records.append(rec)

    if unresolved:
        logging.warning('%d/%d cells could not be matched to a ploidy-file sample; '
                        'first unmatched labels: %s', len(unresolved), len(files),
                        unresolved[0][1] if unresolved else '')
        if args.strict:
            for path, cands in unresolved[:20]:
                sys.stderr.write(F'unresolved: {os.path.basename(path)} tried {cands}\n')
            return 2
    if not records:
        sys.stderr.write('ploidy_eval: no cell could be matched to the ploidy file.\n')
        return 1

    percell = pd.DataFrame(records)
    per_sample, overall = summarize_by_sample(percell, args.ploidy_window)
    overall.update({'n_cells_unresolved': len(unresolved), 'n_bed_files': len(files),
                    'ploidy_file': os.path.abspath(args.ploidy_file),
                    'chroms': args.chroms, 'weight': args.weight,
                    # A string when uncapped, so that the file stays valid JSON for readers
                    # that reject the Infinity literal.
                    'max_cn': ('inf' if is_uncapped(args.max_cn) else float(args.max_cn)),
                    'tool': args.tool, 'donor': args.donor,
                    'sampleType': args.sample_type, 'avgSpotLen': args.avg_spot_len})

    out_dir = os.path.dirname(os.path.abspath(args.output_prefix))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    percell.sort_values(['sample', 'cell']).to_csv(args.output_prefix + '_percell.tsv',
                                                   sep='\t', index=False)
    per_sample.to_csv(args.output_prefix + '_persample.tsv', sep='\t', index=False)
    with open(args.output_prefix + '_summary.json', 'w') as fh:
        json.dump(overall, fh, indent=2, sort_keys=True)
    if args.plot:
        plot_ploidy(percell, per_sample, args.output_prefix, args.ploidy_window, args.title)

    sys.stderr.write(per_sample.to_string(index=False) + '\n')
    sys.stderr.write(F'mean %outliers across samples = '
                     F'{overall["mean_pct_outliers_across_samples"]:.1f}; '
                     F'pooled mean |ploidy distance| = '
                     F'{overall["pooled_mean_abs_ploidy_distance"]:.3f}\n')
    sys.stderr.write(F'Wrote {args.output_prefix}_percell.tsv, _persample.tsv, _summary.json\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
