#!/usr/bin/env python

"""
From https://xsimplechat.com/chat?session=ssn_zkOZ9p4TrStJ&topic=tpc_QjGSOTAvuffT

Generate a clustered CNV heatmap from per-sample BED files.

* Rows = samples (hierarchically clustered).
* Columns = fixed-size genomic bins, kept in genomic order.
* Heatmap uses *integer* copy numbers with a *discrete* colourbar.
* Dashed vertical lines mark chromosome boundaries; chromosome names are
  drawn once, centred under each block.
"""
import argparse, glob, os, re, sys
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm    ### CHANGED
import seaborn as sns

CHROM_ORDER = [F'chr{i}' for i in list(range(1, 23)) + ['X', 'Y']]

def chrom_sort_key(chrom):
    c = str(chrom).replace('chr', '')
    if c == 'X': return 23
    if c == 'Y': return 24
    try: return int(c)
    except Exception: return 99

def load_bed(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line or line.startswith('#'): continue
            toks = line.split('\t')
            if len(toks) < 4: continue
            try:
                rows.append((toks[0],
                             int(round(float(toks[1]))),
                             int(round(float(toks[2]))),
                             float(toks[3])))
            except ValueError:
                continue
    return pd.DataFrame(rows, columns=['chrom', 'start', 'end', 'cn'])

def sample_name_from_path(path, regex=None):
    base = re.sub(r'\.bed$', '', os.path.basename(path))
    if regex:
        m = re.search(regex, base)
        if m: return m.group(1) if m.groups() else m.group(0)
    return base

def bed_to_chrom_means(df, chroms=CHROM_ORDER):
    out = {}
    for chrom in chroms:
        sub = df[df['chrom'] == chrom]
        if sub.empty:
            out[chrom] = np.nan
        else:
            w = (sub['end'] - sub['start']).clip(lower=1).astype(float)
            out[chrom] = float((sub['cn'] * w).sum() / w.sum())
    return out

def bed_to_fixed_bins(df, bin_size, chrom_sizes):
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
                cols.append((label, np.nan)); continue
            ov = (np.minimum(sub['end'], be) - np.maximum(sub['start'], bs)
                  ).clip(lower=0).astype(float)
            tot = ov.sum()
            cols.append((label, float((sub['cn'] * ov).sum() / tot) if tot > 0 else np.nan))
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
    p = argparse.ArgumentParser(
        description='Generate a clustered CNV heatmap from per-sample BED files.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('-i', '--input', nargs='+', required=True)
    p.add_argument('-o', '--output-prefix', required=True)
    p.add_argument('--title', default='')
    p.add_argument('--by-chrom', action='store_true')
    p.add_argument('--bin-size', type=int, default=0)
    p.add_argument('--fai',    type=str, default='')
    p.add_argument('--sample-regex', type=str, default='')
    p.add_argument('--cmap',   type=str, default='RdBu_r')
    p.add_argument('--vmin',   type=int, default=0)           ### CHANGED: int, not float
    p.add_argument('--vmax',   type=int, default=6)           ### CHANGED: int, not float
    p.add_argument('--center', type=int, default=2,           ### CHANGED: diploid value
                   help='Diploid copy number used to fill NaN bins')
    p.add_argument('--show-sample-labels', type=bool, default=1, # action='store_true',
                   help='Print sample names on y-axis (on by default; matches the reference figure)')
    args = p.parse_args()

    # ---- collect input files ------------------------------------------------
    files = []
    for pat in args.input:
        files.extend(sorted(glob.glob(pat)) if any(c in pat for c in '*?[') else [pat])
    files = [f for f in files if os.path.isfile(f) and os.path.getsize(f) > 0]
    if not files:
        sys.stderr.write('No input BED files found.\n'); sys.exit(1)

    use_fixed_bins = (args.bin_size > 0)
    if use_fixed_bins:
        if not args.fai:
            sys.stderr.write('--bin-size requires --fai\n'); sys.exit(1)
        chrom_sizes = {c: s for c, s in parse_chrom_sizes(args.fai).items()
                       if c in CHROM_ORDER}

    rows = {}
    for path in files:
        name = sample_name_from_path(path, args.sample_regex or None)
        df = load_bed(path)
        if df.empty:
            sys.stderr.write(F'Warning: {path} has no usable rows; skipping.\n'); continue
        if use_fixed_bins:
            rows[name] = pd.Series({lbl: v for lbl, v in bed_to_fixed_bins(
                df, args.bin_size, chrom_sizes)})
        else:
            rows[name] = pd.Series(bed_to_chrom_means(df))

    if not rows:
        sys.stderr.write('No data to plot.\n'); sys.exit(1)

    samp_names = list(rows.keys())
    common_prefix = os.path.commonprefix(samp_names)
    common_preadd = common_prefix.split('_')[-1]
    common_suffix = (os.path.commonprefix([x[::-1] for x in samp_names]))[::-1]
    def rm_suffix(x, suffix): return (x[0:-len(suffix)] if x.endswith(suffix) else x)
    samp_name_old2new = {x: rm_suffix(x.removeprefix(common_prefix), common_suffix) for x in samp_names}
    rows = {(common_preadd + samp_name_old2new[sn]) : val for (sn, val) in rows.items()}
    mat = pd.DataFrame(rows).T
    if not use_fixed_bins:
        mat = mat[[c for c in CHROM_ORDER if c in mat.columns]]

    # ============== CHANGED: integer copy numbers ===========================
    mat = mat.round().clip(lower=args.vmin, upper=args.vmax)
    mat = mat.dropna(axis=1, how='all')
    fill = mat.fillna(args.center)

    # ============== CHANGED: discrete colormap + boundary norm ==============
    n_levels = args.vmax - args.vmin + 1
    base = plt.get_cmap(args.cmap, n_levels)
    discrete_cmap = ListedColormap([base(i) for i in range(n_levels)])
    norm = BoundaryNorm(np.arange(args.vmin - 0.5, args.vmax + 1.5, 1.0),
                        discrete_cmap.N)

    figsize = (max(8, min(0.09 * mat.shape[1] + 6, 12)),
               max(8, min(0.18 * mat.shape[0] + 3, 12)))
    print(f"figsize={figsize}")
    g = sns.clustermap(
        fill,
        row_cluster=True, col_cluster=False,
        method='average', metric='euclidean',
        # cmap='coolwarm',
        cmap=discrete_cmap, norm=norm,                     ### CHANGED: norm replaces vmin/vmax/center
        figsize=figsize,
        cbar_kws={'label':  'Copy numbers',
                  'ticks':  np.arange(args.vmin, args.vmax + 1),
                  'spacing': 'proportional',                 ### CHANGED: discrete colorbar
                  'orientation': 'horizontal'},              ### CHANGED: horizontal color bar  
                  xticklabels=False,                         ### CHANGED: hide per-bin x labels
        yticklabels=args.show_sample_labels,
        dendrogram_ratio=(0.15, 0.055),                      ### https://chat.deepseek.com/a/chat/s/e2f22ab0-ceb2-47f1-8e8a-005ef1810b65
        # center = 2,
    )
    # Adjust position: [x0, y0, width, height]
    # This example places it below the heatmap
    g.ax_cbar.set_position([0.25, 0.98, 0.5, 0.01])
    # g.ax_cbar.set_title('Copy numbers')

    # Change tick length
    g.ax_cbar.tick_params(axis='x', length=3)

    ax = g.ax_heatmap
    ax.set_xlabel(''); ax.set_ylabel('')

    # ============== CHANGED: chromosome dividers + centred labels ===========
    col_chroms = [c.split(':')[0] for c in fill.columns]     # works for both modes

    prev = None
    for i, c in enumerate(col_chroms):
        if prev is not None and c != prev:
            ax.axvline(i, color='black', linestyle='--', linewidth=0.8)
        prev = c

    pos = {}
    for i, c in enumerate(col_chroms):
        pos.setdefault(c, []).append(i)
    centers, labels = [], []
    for c in CHROM_ORDER:
        if c in pos:
            centers.append((pos[c][0] + pos[c][-1] + 1) / 2.0)
            labels.append(c.replace('chr', ''))
    ax.set_xticks(centers)
    ax.set_xticklabels(labels, rotation=15, fontsize=8)
    ax.tick_params(axis='x', length=0)
    # ========================================================================

    if args.show_sample_labels:
        plt.setp(ax.get_yticklabels(), rotation=0, fontsize=6)

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
