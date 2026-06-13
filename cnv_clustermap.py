#!/usr/bin/env python
"""
Generate a clustered heatmap (clustermap) of CNV calls from a set of per-sample BED files.

Each input BED file has the format produced by cnv_raw_to_bed.py:
    #chr_37  start_37  end_37  obsCN
The sample name is taken from the file basename (strip directory and any .bed suffix), or
optionally from a regex extracted with --sample-regex.

Rows of the clustermap are samples, columns are genomic bins (one column per chromosome
slice when --by-chrom is used, otherwise one per fixed-size bin).

Hierarchical clustering is applied to rows (samples) only; columns are kept in genomic
order to preserve interpretability.

The figure is written to PDF and PNG (same basename as -o).
"""
import argparse, glob, os, re, sys
import numpy as np
import pandas as pd

# Use a non-interactive backend so the script runs on headless cluster nodes.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns


CHROM_ORDER = [F'chr{i}' for i in list(range(1, 23)) + ['X', 'Y']]


def chrom_sort_key(chrom):
    c = str(chrom).replace('chr', '')
    if c == 'X': return 23
    if c == 'Y': return 24
    try: return int(c)
    except Exception: return 99


def load_bed(path):
    """Load a single BED file produced by cnv_raw_to_bed.py."""
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line: continue
            if line.startswith('#'): continue
            toks = line.split('\t')
            if len(toks) < 4: continue
            chrom = toks[0]
            try:
                start = int(round(float(toks[1])))
                end = int(round(float(toks[2])))
                cn = float(toks[3])
            except ValueError:
                continue
            rows.append((chrom, start, end, cn))
    return pd.DataFrame(rows, columns=['chrom', 'start', 'end', 'cn'])


def sample_name_from_path(path, regex=None):
    base = os.path.basename(path)
    base = re.sub(r'\.bed$', '', base)
    if regex:
        m = re.search(regex, base)
        if m: return m.group(1) if m.groups() else m.group(0)
    return base


def bed_to_chrom_means(df, chroms=CHROM_ORDER):
    """Aggregate CNV per chromosome by length-weighted mean."""
    out = {}
    for chrom in chroms:
        sub = df[df['chrom'] == chrom]
        if sub.empty:
            out[chrom] = np.nan
        else:
            widths = (sub['end'] - sub['start']).clip(lower=1).astype(float)
            out[chrom] = float((sub['cn'] * widths).sum() / widths.sum())
    return out


def bed_to_fixed_bins(df, bin_size, chrom_sizes):
    """Aggregate CNV onto fixed-size genomic bins; return ordered (label, value) list."""
    cols = []
    for chrom in sorted(chrom_sizes, key=chrom_sort_key):
        chrom_len = chrom_sizes[chrom]
        sub = df[df['chrom'] == chrom]
        n_bins = max(1, int(np.ceil(chrom_len / bin_size)))
        for b in range(n_bins):
            bs = b * bin_size
            be = min((b + 1) * bin_size, chrom_len)
            label = F'{chrom}:{bs}-{be}'
            if sub.empty:
                cols.append((label, np.nan))
                continue
            overlap_start = np.maximum(sub['start'], bs)
            overlap_end = np.minimum(sub['end'], be)
            ov = (overlap_end - overlap_start).clip(lower=0).astype(float)
            tot = ov.sum()
            if tot <= 0:
                cols.append((label, np.nan))
            else:
                cols.append((label, float((sub['cn'] * ov).sum() / tot)))
    return cols


def parse_chrom_sizes(fai_path):
    sizes = {}
    with open(fai_path) as fh:
        for line in fh:
            toks = line.rstrip('\n').split('\t')
            if len(toks) >= 2:
                try: sizes[toks[0]] = int(toks[1])
                except ValueError: pass
    return sizes


def main():
    parser = argparse.ArgumentParser(description='Generate a clustered CNV heatmap from per-sample BED files.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', nargs='+', required=True,
        help='BED files (or glob patterns) produced by cnv_raw_to_bed.py')
    parser.add_argument('-o', '--output-prefix', required=True,
        help='Output filename prefix; files <prefix>.pdf and <prefix>.png will be written')
    parser.add_argument('--title', default='', help='Title shown above the heatmap')
    parser.add_argument('--by-chrom', action='store_true',
        help='Use one column per chromosome (length-weighted mean CN); default mode')
    parser.add_argument('--bin-size', type=int, default=0,
        help='If >0, use fixed-size genomic bins of this many bp (requires --fai)')
    parser.add_argument('--fai', type=str, default='',
        help='Reference fasta index (.fai); needed when --bin-size is given')
    parser.add_argument('--sample-regex', type=str, default='',
        help='Regex applied to BED basename to extract sample name (first capture group)')
    parser.add_argument('--cmap', type=str, default='RdBu_r', help='Matplotlib colormap')
    parser.add_argument('--vmin', type=float, default=0.0)
    parser.add_argument('--vmax', type=float, default=4.0)
    args = parser.parse_args()

    # Expand globs.
    files = []
    for pat in args.input:
        if any(c in pat for c in '*?['):
            files.extend(sorted(glob.glob(pat)))
        else:
            files.append(pat)
    files = [f for f in files if os.path.isfile(f) and os.path.getsize(f) > 0]
    if not files:
        sys.stderr.write('No input BED files found.\n')
        sys.exit(1)

    use_fixed_bins = (args.bin_size > 0)
    if use_fixed_bins:
        if not args.fai:
            sys.stderr.write('--bin-size requires --fai\n'); sys.exit(1)
        chrom_sizes = parse_chrom_sizes(args.fai)
        chrom_sizes = {c: s for c, s in chrom_sizes.items() if c in CHROM_ORDER}

    rows = {}
    for path in files:
        name = sample_name_from_path(path, args.sample_regex or None)
        df = load_bed(path)
        if df.empty:
            sys.stderr.write(F'Warning: {path} has no usable rows; skipping.\n')
            continue
        if use_fixed_bins:
            cols = bed_to_fixed_bins(df, args.bin_size, chrom_sizes)
            rows[name] = pd.Series({lbl: v for lbl, v in cols})
        else:
            rows[name] = pd.Series(bed_to_chrom_means(df))

    if not rows:
        sys.stderr.write('No data to plot.\n'); sys.exit(1)

    mat = pd.DataFrame(rows).T
    # Preserve genomic column order.
    if not use_fixed_bins:
        mat = mat[[c for c in CHROM_ORDER if c in mat.columns]]
    # Drop columns that are all-NaN, then fill remaining NaNs with the row mean (for clustering).
    mat = mat.dropna(axis=1, how='all')
    fill = mat.apply(lambda r: r.fillna(r.mean()), axis=1)

    # Cluster only rows; preserve column order for biological interpretability.
    figsize = (max(8, min(0.25 * mat.shape[1] + 4, 24)),
               max(4, min(0.25 * mat.shape[0] + 2, 24)))
    g = sns.clustermap(
        fill,
        row_cluster=True, col_cluster=False,
        method='average', metric='euclidean',
        cmap=args.cmap, vmin=args.vmin, vmax=args.vmax, center=2.0,
        figsize=figsize,
        cbar_kws={'label': 'Copy number'},
        xticklabels=True, yticklabels=True,
    )
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=8)
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0,  fontsize=8)
    if args.title:
        g.fig.suptitle(args.title, y=1.02)

    out_dir = os.path.dirname(os.path.abspath(args.output_prefix))
    if out_dir: os.makedirs(out_dir, exist_ok=True)
    g.savefig(args.output_prefix + '.pdf', bbox_inches='tight')
    g.savefig(args.output_prefix + '.png', dpi=150, bbox_inches='tight')
    plt.close('all')
    sys.stderr.write(F'Wrote {args.output_prefix}.pdf and {args.output_prefix}.png\n')


if __name__ == '__main__':
    main()
