#!/usr/bin/env python3
"""Synthetic end-to-end test for bench_results/stat_tests.py (v3, donor-level).

Part A  Fig. 2 mode (long TSV): 300 simulated cells built from 12 shared
        material units (accession_1 x accession_2 x cellLine) belonging to
        NINE donors, with donor-correlated caller effects - the same
        correlation structure that the real benchmark has (all cells of one
        donor share that donor's haploid material).  Verifies the default
        donor-level outputs, a fine-grained (shared-material) sensitivity
        run, the naive (--cluster-key none) fallback, and the diagnostics
        columns.
Part B  Fig. 3 mode (pct_within_long.tsv): datasets of one plot group share
        donors (germline groups) or are independent (ACT, empty donor).
Part C  demo_independence_failure: Monte-Carlo simulation under a TRUE null
        (no caller difference) with clustered per-cell differences; measures
        the empirical Type I error of the naive per-cell Wilcoxon vs. the
        cluster-level Wilcoxon.  This is the reproducible demonstration of
        what happens when the independence assumption fails.
"""
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
OUT = os.path.join(REPO, 'work_test_out')
os.makedirs(OUT, exist_ok=True)
STAT = os.path.join(HERE, 'stat_tests.py')
sys.path.insert(0, HERE)
import stat_tests  # noqa: E402  (for the Part C demo)

rng = np.random.default_rng(7)

# ------------------------------------------------------------------ Fig. 2 --
callers = ['aneufinder', 'chisel', 'copynumber', 'flcna', 'ginkgo', 'hmmcopy',
           'sccnv', 'scyn', 'secnv']
metrics = ['accuracy', 'PCC_intCN', 'PCC_nonintCN', 'frac_cov_genome',
           'breakpoint_precision', 'breakpoint_recall', 'breakpoint_f1score']
scenarios = ['with_haploidy_assumed_gametes', 'with_aneuploidy_aware_gametes']

# 12 shared-material units, mapped onto 9 donors: (acc1, acc2, cellLine)
units = []
for i in range(12):
    units.append((F'SRR{100000 + i}', F'SRR{200000 + i * 5}',
                  ['COLO-829', 'HCC1395', 'HeLa'][i % 3]))
donor_of_unit = [i % 9 for i in range(len(units))]   # nine donors, as in the benchmark
N_PER_UNIT = 25            # cells per material unit -> 300 cells total
cells = []
for u_idx, (a1, a2, cl) in enumerate(units):
    for j in range(N_PER_UNIT):
        # (overall_ploidy, CNA_percent) unique per j -> unique cell keys
        cells.append((a1, a2, cl,
                      ['diploid', 'aneuploid'][j % 2], 20.0 + 1.0 * j))

# per-cell baseline, caller effect (ginkgo best), and DONOR x caller
# interaction: the same-caller results of one donor are correlated, exactly
# like downsamplings of the same donor's haplotype BAMs.
cell_base = rng.beta(5, 2, len(cells))
caller_effect = {'ginkgo': 0.10, 'aneufinder': 0.02, 'chisel': -0.05,
                 'copynumber': -0.02, 'flcna': 0.0, 'hmmcopy': -0.06,
                 'sccnv': -0.01, 'scyn': -0.03, 'secnv': 0.01}
unit_of_cell = [i // N_PER_UNIT for i in range(len(cells))]
donor_caller_noise = {(donor_of_unit[u], c): rng.normal(0, 0.05)
                      for u in range(len(units)) for c in callers}

rows = []
for c in callers:
    for i, (a1, a2, cl, op, cna) in enumerate(cells):
        u = unit_of_cell[i]
        d = donor_of_unit[u]
        row = {'Caller': c, 'accession_1': a1, 'accession_2': a2, 'cellLine': cl,
               'overall_ploidy': op, 'CNA_percent': cna,
               'donor': F'donor{d}', 'sampleType': ['sperm', 'PB1', 'PB2'][u % 3],
               'avgSpotLen': 100, 'observed_ploidy': 2.0 + rng.normal(0, 0.3),
               'raw_total_sequences': int(1e6), 'reads_mapped': int(9e5),
               'bases_mapped_cigar': int(6e7)}
        for sc in scenarios:
            for m in metrics:
                base = (cell_base[i] + caller_effect[c]
                        + donor_caller_noise[(d, c)]     # donor-shared effect
                        + rng.normal(0, 0.04))
                if m in ('PCC_intCN', 'PCC_nonintCN'):
                    base = 0.7 * base
                row[F'{sc}.{m}'] = float(np.clip(base, -1, 1))
        rows.append(row)
long_df = pd.DataFrame(rows)
long_tsv = os.path.join(OUT, 'synthetic.long.tsv')
long_df.to_csv(long_tsv, sep='\t', index=False, na_rep='NA')
print(F'wrote {long_tsv} ({len(long_df)} rows, {len(units)} material units '
      F'x {N_PER_UNIT} cells across {len(set(donor_of_unit))} donors)')

with open(long_tsv) as fh:
    ret = subprocess.run([sys.executable, STAT, '-o', os.path.join(OUT, 'fig2stats'),
                          '--reference', 'ginkgo', '--boot', '500'],
                         stdin=fh, capture_output=True, text=True)
print('Fig. 2 mode (donor level, default) exit code:', ret.returncode)
if ret.returncode != 0:
    print(ret.stdout); print(ret.stderr)
    sys.exit(1)

pw = pd.read_csv(os.path.join(OUT, 'fig2stats.stats.pairwise.tsv'), sep='\t')
fr = pd.read_csv(os.path.join(OUT, 'fig2stats.stats.friedman.tsv'), sep='\t')
need_cols = ['inference_level', 'cluster_key', 'n_clusters', 'n_cells_paired',
             'pvalue_two_sided', 'pvalue_sign_test', 'pvalue_holm',
             'rank_biserial_r', 'cl_effect_paired', 'ci95_median_diff_low',
             'ci95_median_diff_high', 'pvalue_cell_naive',
             'rank_biserial_r_cell_naive', 'icc_within_cluster_d',
             'design_effect', 'n_effective_cells', 'p_inflation_ratio']
missing = [c for c in need_cols if c not in pw.columns]
assert not missing, F'pairwise TSV missing columns: {missing}'
assert set(pw['inference_level']) == {'cluster'}, set(pw['inference_level'])
assert pw['n_clusters'].eq(9).all(), pw['n_clusters'].unique()
assert pw['n_cells_paired'].eq(300).all()
assert fr['level'].eq('cluster').any() and fr['level'].eq('cell (naive)').any()
assert fr.loc[fr['level'] == 'cluster', 'n_blocks'].eq(9).all()
assert fr.loc[fr['level'] == 'cell (naive)', 'n_blocks'].eq(300).all()
caller_rows = pw[pw['caller_b'] != '<scenario>']       # caller-vs-caller rows
scen_rows = pw[pw['caller_b'] == '<scenario>']         # Hap_0-vs-Hap_1 rows
# caller-vs-caller differences carry a cluster-shared component in this
# synthetic design, so their ICC must be positive; the scenario-difference
# rows are built WITHOUT a cluster-shared component, so their ICC hovers
# around 0 (method-of-moments estimates may be slightly negative - design
# effect correctly stays 1).
assert pw['icc_within_cluster_d'].between(-1, 1).all()
assert caller_rows['icc_within_cluster_d'].gt(0).all()
assert pw['design_effect'].ge(1).all()
assert pw['n_effective_cells'].le(300).all()
assert caller_rows['p_inflation_ratio'].median() >= 1.0
with open(os.path.join(OUT, 'fig2stats.stats.json')) as fh:
    js = json.load(fh)
assert js['independence']['cluster_key'] == ['donor']
assert js['independence']['n_clusters'] == 9
print(F'OK fig2stats: {len(pw)} pairwise rows, icc median '
      F'{pw["icc_within_cluster_d"].median():.3f}, design effect median '
      F'{pw["design_effect"].median():.2f}, n_eff median '
      F'{pw["n_effective_cells"].median():.1f}, P-inflation median '
      F'{pw["p_inflation_ratio"].median():.1f}x')

# fine-grained sensitivity key: the 12 shared-material units (accession pair
# x cellLine), i.e. NOT the default donor-level analysis
with open(long_tsv) as fh:
    ret = subprocess.run([sys.executable, STAT, '-o', os.path.join(OUT, 'fig2fine'),
                          '--reference', 'ginkgo', '--boot', '300',
                          '--cluster-key', 'accession_1,accession_2,cellLine'],
                         stdin=fh, capture_output=True, text=True)
print('Fig. 2 mode (--cluster-key accession_1,accession_2,cellLine) exit code:', ret.returncode)
if ret.returncode != 0:
    print(ret.stdout); print(ret.stderr)
    sys.exit(1)
pw_fine = pd.read_csv(os.path.join(OUT, 'fig2fine.stats.pairwise.tsv'), sep='\t')
assert pw_fine['cluster_key'].eq('accession_1|accession_2|cellLine').all()
print(F'OK fig2fine: material-unit clusters = {sorted(pw_fine["n_clusters"].unique())}')

# naive mode (--cluster-key none): v1 behaviour, explicitly flagged
with open(long_tsv) as fh:
    ret = subprocess.run([sys.executable, STAT, '-o', os.path.join(OUT, 'fig2naive'),
                          '--reference', 'ginkgo', '--boot', '300',
                          '--cluster-key', 'none'],
                         stdin=fh, capture_output=True, text=True)
print('Fig. 2 mode (--cluster-key none, naive) exit code:', ret.returncode)
if ret.returncode != 0:
    print(ret.stdout); print(ret.stderr)
    sys.exit(1)
pw_naive = pd.read_csv(os.path.join(OUT, 'fig2naive.stats.pairwise.tsv'), sep='\t')
assert set(pw_naive['inference_level']) == {'cell (naive)'}
assert pw_naive['note'].str.contains('NAIVE').any()
with open(os.path.join(OUT, 'fig2naive.stats.json')) as fh:
    js_n = json.load(fh)
assert 'NAIVE' in js_n['independence']['inference_level']
print('OK fig2naive: naive mode flagged and complete')

# ------------------------------------------------------------------ Fig. 3 --
tools = ['aneufinder', 'chisel', 'copynumber', 'flcna', 'ginkgo', 'hmmcopy',
         'sccnv', 'scyn', 'secnv', 'scabsolute']
plots = {
    'COLO-829': ([F'COLO-{i}' for i in range(6)],
                 ['donorA', 'donorA', 'donorB', 'donorB', 'donorC', 'donorC']),
    'HCC1395': ([F'HCC-{i}' for i in range(5)],
                ['donorD', 'donorD', 'donorE', 'donorE', 'donorF']),
    'HeLa': ([F'HeLa-{i}' for i in range(4)],
             ['donorG', 'donorG', 'donorH', 'donorH']),
    'ACT': ([F'TN{i}' for i in range(1, 9)] + ['MDA-MB-453', 'SK-BR-3', 'T47D', 'ZR-75-1'],
            [''] * 12),          # real tumour samples: no donor metadata
}
ploidy_effect = {'ginkgo': 6.0, 'aneufinder': 2.0, 'chisel': -3.0,
                 'copynumber': -1.0, 'flcna': 0.0, 'hmmcopy': -4.0,
                 'sccnv': -0.5, 'scyn': -2.0, 'secnv': 1.0, 'scabsolute': 0.5}
rows = []
for plot, (datasets, donors) in plots.items():
    donor_shock = {d: rng.normal(0, 4.0) for d in set(donors) if d}
    for ds, dn in zip(datasets, donors):
        for t in tools:
            caps = ['10', 'inf'] if t != 'scabsolute' else ['inf']
            for cap in caps:
                if t == 'scyn' and plot == 'ACT':
                    continue  # SCYN produced no output on ACT
                if t == 'chisel' and plot == 'ACT':
                    continue  # not applicable on ACT
                shared = (donor_shock.get(dn, 0.0)            # donor-shared shock
                          + ploidy_effect.get(t, 0.0))
                pct = np.clip(50 + shared + rng.normal(0, 6), 0, 100)
                rows.append({'plot': plot, 'dataset': ds, 'tool': t, 'max_cn': cap,
                             'method': F'{t}|{cap}', 'window': 0.5,
                             'n_cells': 300, 'n_cells_finite': 250 + int(rng.integers(0, 50)),
                             'pct_within': pct, 'mean_abs_ploidy_error': abs(rng.normal(0, 0.5)),
                             'failed': False, 'donor': dn, 'sampleType': '',
                             'avgSpotLen': '', 'cellLine': plot if plot == 'ACT' else ''})
ploidy_df = pd.DataFrame(rows)
ploidy_tsv = os.path.join(OUT, 'synthetic_pct_within_long.tsv')
ploidy_df.to_csv(ploidy_tsv, sep='\t', index=False, na_rep='NA')
print(F'wrote {ploidy_tsv} ({len(ploidy_df)} rows)')

ret = subprocess.run([sys.executable, STAT, '-i', ploidy_tsv,
                      '-o', os.path.join(OUT, 'fig3stats'), '--reference', 'ginkgo|10',
                      '--boot', '500'], capture_output=True, text=True)
print('Fig. 3 mode (donor-level clusters, default) exit code:', ret.returncode)
if ret.returncode != 0:
    print(ret.stdout); print(ret.stderr)
    sys.exit(1)

pw3 = pd.read_csv(os.path.join(OUT, 'fig3stats.stats.pairwise.tsv'), sep='\t')
fr3 = pd.read_csv(os.path.join(OUT, 'fig3stats.stats.friedman.tsv'), sep='\t')
with open(os.path.join(OUT, 'fig3stats.stats.json')) as fh:
    js3 = json.load(fh)
per_group = js3['independence']['per_plot_group']
assert per_group['COLO-829']['cluster_column'] == 'donor'
assert per_group['COLO-829']['n_clusters'] == 3
assert per_group['HCC1395']['n_clusters'] == 3
assert per_group['HeLa']['n_clusters'] == 2
assert per_group['ACT']['cluster_column'] == 'dataset', per_group['ACT']
germline = pw3[pw3['plot'] != 'ACT']
assert germline['cluster_key'].eq('donor').all()
assert fr3['level'].eq('cluster').any() and fr3['level'].eq('dataset (naive)').any()
print(F'OK fig3stats: cluster columns per group '
      F'{ {k: v.get("cluster_column") for k, v in per_group.items()} }, '
      F'n_clusters { {k: v.get("n_clusters") for k, v in per_group.items()} }')

for f in ['fig2stats.stats.pairwise.tsv', 'fig2stats.stats.friedman.tsv',
          'fig2stats.stats.concordance.tsv', 'fig2stats.stats.json',
          'fig2fine.stats.pairwise.tsv', 'fig2naive.stats.pairwise.tsv',
          'fig3stats.stats.pairwise.tsv', 'fig3stats.stats.friedman.tsv',
          'fig3stats.stats.json']:
    p = os.path.join(OUT, f)
    assert os.path.isfile(p), F'MISSING {p}'
    print(F'OK {f}: {os.path.getsize(p)} bytes')

# ------------------------------------------------------------------ Part C --
def demo_independence_failure(n_reps=300, n_clusters=45, n_total=1989,
                              icc=0.3, alpha=0.05, seed=0):
    """What happens if per-cell results of the same caller are treated as
    independent although they are clustered: empirical Type I error of the
    naive per-cell Wilcoxon signed-rank test vs. the cluster-level test,
    under a TRUE null (median difference exactly 0).

    The cluster structure mirrors the benchmark: ~45 (donor-pair x template)
    material units of ~44 cells each (1,989 cells in total), paired
    differences d_i = shared cluster effect + cell-level noise, with the
    between/within variance ratio set to hit the target ICC."""
    rs = np.random.default_rng(seed)
    sizes = np.full(n_clusters, n_total // n_clusters)
    sizes[:n_total - sizes.sum()] += 1
    var_b = icc / (1.0 - icc)          # within-cluster variance = 1
    sd_b = np.sqrt(var_b)
    rej_naive = rej_cluster = rej_sign = 0
    pvs_naive, pvs_cluster = [], []
    for _ in range(n_reps):
        g = rs.normal(0.0, sd_b, n_clusters)
        d = np.concatenate([np.full(s, gv) for s, gv in zip(sizes, g)])
        d = d + rs.normal(0.0, 1.0, n_total)      # cell-level noise
        cl = np.repeat(np.arange(n_clusters), sizes)
        # naive per-cell test (v1 behaviour)
        w_cell = stat_tests.wilcoxon_signed_rank(d, np.zeros_like(d))
        rej_naive += int(w_cell['pvalue'] <= alpha)
        pvs_naive.append(w_cell['pvalue'])
        # donor/cluster-level test (default): per-cluster medians
        cm = pd.Series(d).groupby(cl).median().to_numpy()
        w_cl = stat_tests.wilcoxon_signed_rank(cm, np.zeros_like(cm))
        rej_cluster += int(w_cl['pvalue'] <= alpha)
        pvs_cluster.append(w_cl['pvalue'])
        st = stat_tests.sign_test_two_sided(cm)
        rej_sign += int(st['pvalue'] <= alpha)
    rate_naive = rej_naive / n_reps
    rate_cluster = rej_cluster / n_reps
    rate_sign = rej_sign / n_reps
    print('\n--- demo_independence_failure (TRUE null, no caller difference) ---')
    print(F'  simulated structure : {n_clusters} clusters x ~{n_total // n_clusters} cells '
          F'= {n_total} per-cell differences, ICC target {icc}')
    print(F'  design effect at this ICC ~ 1 + (m-1)*ICC ~ {1 + (n_total // n_clusters - 1) * icc:.0f} '
          F'(effective n ~ {n_total / (1 + (n_total // n_clusters - 1) * icc):.0f} cells)')
    print(F'  naive per-cell Wilcoxon  : rejects H0 in {rej_naive:3d}/{n_reps} runs '
          F'({100 * rate_naive:.0f}% Type I error; nominal {100 * alpha:.0f}%)')
    print(F'  cluster-level Wilcoxon   : rejects H0 in {rej_cluster:3d}/{n_reps} runs '
          F'({100 * rate_cluster:.0f}% Type I error)')
    print(F'  cluster-level sign test  : rejects H0 in {rej_sign:3d}/{n_reps} runs '
          F'({100 * rate_sign:.0f}% Type I error)')
    print(F'  median naive P = {np.median(pvs_naive):.2e}; '
          F' median cluster P = {np.median(pvs_cluster):.3f}')
    # the cluster-level test must stay near the nominal level
    lo, hi = alpha - 3 * np.sqrt(alpha * (1 - alpha) / n_reps), alpha + 3 * np.sqrt(alpha * (1 - alpha) / n_reps)
    assert lo <= rate_cluster <= hi, F'cluster-level Type I error {rate_cluster:.3f} outside [{lo:.3f}, {hi:.3f}]'
    # while the naive per-cell test must be grossly anticonservative
    assert rate_naive >= 3 * alpha, (F'naive Type I error {rate_naive:.3f} not clearly '
                                     F'anticonservative; the demo lost its point')
    return rate_naive, rate_cluster


demo_independence_failure()
print('\nAll stat_tests end-to-end tests PASSED '
      '(including the independence-failure demonstration)')
