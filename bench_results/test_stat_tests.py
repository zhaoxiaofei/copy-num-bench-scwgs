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
import re
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
metrics = ['intCN_accuracy', 'intCN_PCC', 'nonintCN_PCC', 'CN_genome_cov_frac',
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
                if m in ('intCN_PCC', 'nonintCN_PCC'):
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
             'rank_biserial_r', 'ci95_r_low', 'ci95_r_high',
             'cl_effect_paired', 'ci95_median_diff_low',
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
# 95% CI of the effect size r: plausible bounds and it brackets the estimate
assert pw['ci95_r_low'].between(-1.01, 1.01).all() and pw['ci95_r_high'].between(-1.01, 1.01).all()
assert (pw['ci95_r_low'] <= pw['ci95_r_high']).all()
_brackets = ((pw['ci95_r_low'] <= pw['rank_biserial_r'])
             & (pw['rank_biserial_r'] <= pw['ci95_r_high']))
assert _brackets.mean() >= 0.9, F'bootstrap CI of r brackets the estimate only in {_brackets.mean() * 100:.0f}% of rows'
print(F'OK fig2stats: 95% CI of r present ({_brackets.mean() * 100:.0f}% of rows bracketed)')
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
assert pw3['ci95_r_low'].notna().any() and pw3['ci95_r_high'].notna().any()
assert (pw3['ci95_r_low'] <= pw3['ci95_r_high']).all()
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

# ------------------------------------------------- Part B2: merged runs ---
# Repeated runs of ONE biological sample must count as ONE effective sample
# in every statistical consumer (pairwise, Friedman, effect sizes and MCB):
# MDA-MB-231 enters the ACT arm as MDAMB231c28 / MDAMB231c8 / MDAMB231_popp31
# and TN6 as a 36 bp + a 152 bp run.  The per-run values below have easy
# medians (MDA-MB-231: median(90, 80, 100) = 90; TN6: median(60, 40) = 50).
assert stat_tests.canonical_ploidy_sample('MDAMB231c28') == 'MDAMB231'
assert stat_tests.canonical_ploidy_sample('SRP259526_MDAMB231c8') == 'MDAMB231'
assert stat_tests.canonical_ploidy_sample('SRP259526_MDAMB231_popp31') == 'MDAMB231'
assert stat_tests.canonical_ploidy_sample('MDAMB231c28 · 50 bp') == 'MDAMB231'
assert stat_tests.canonical_ploidy_sample('TN6 · 152 bp') == 'TN6'
assert stat_tests.canonical_ploidy_sample('SRP259526_TN6') == 'TN6'
assert stat_tests.canonical_ploidy_sample('S01') == 'S01'
assert stat_tests.canonical_ploidy_sample('345HS1') == '345HS1'
# a label not listed in the table stays its own unit (and is reported)
assert stat_tests.canonical_ploidy_sample('NEWSAMPLE · 50 bp') == 'NEWSAMPLE · 50 bp'
assert stat_tests.ploidy_sample_key('NEWSAMPLE · 50 bp') is None
assert stat_tests.PLOIDY_SAMPLE_DATASETS['MDAMB231'] == (
    'MDAMB231c28', 'MDAMB231c8', 'MDAMB231_popp31')
print('OK fig3merge: canonical effective-sample ids (MDA-MB-231 x3, TN6 runs)')

MERGE_METHODS = [('ginkgo', '10'), ('aneufinder', '10'), ('hmmcopy', '10')]
merge_rows = []
for ds, donor, spot, vals in [
        ('MDAMB231c28 · 50 bp', 'SRP259526_MDAMB231c28', '50', (90.0, 80.0, 70.0)),
        ('MDAMB231c8 · 50 bp', 'SRP259526_MDAMB231c8', '50', (80.0, 70.0, 60.0)),
        ('MDAMB231_popp31 · 50 bp', 'SRP259526_MDAMB231_popp31', '50', (100.0, 90.0, 50.0)),
        ('TN6 · 36 bp', 'SRP259526_TN6', '36', (60.0, 50.0, 40.0)),
        ('TN6 · 152 bp', 'SRP259526_TN6', '152', (40.0, 30.0, 20.0)),
        ('TN1 · 50 bp', 'SRP259526_TN1', '50', (70.0, 60.0, 30.0))]:
    for (tool, cap), pct in zip(MERGE_METHODS, vals):
        merge_rows.append({'plot': 'ACT', 'dataset': ds, 'tool': tool, 'max_cn': cap,
                           'method': F'{tool}|{cap}', 'window': 0.5,
                           'n_cells': 300, 'n_cells_finite': 280, 'pct_within': pct,
                           'mean_abs_ploidy_error': 0.2, 'failed': False, 'donor': donor,
                           'sampleType': 'NA_ILLUMINA', 'avgSpotLen': spot, 'cellLine': ''})
merge_df = pd.DataFrame(merge_rows)
merge_tsv = os.path.join(OUT, 'fig3merge_pct_within_long.tsv')
merge_df.to_csv(merge_tsv, sep='\t', index=False)

# B2.1 standalone stat_tests: the ACT group collapses 6 datasets into 3 units
ret = subprocess.run([sys.executable, STAT, '-i', merge_tsv,
                      '-o', os.path.join(OUT, 'fig3merge_raw'),
                      '--reference', 'ginkgo|10', '--boot', '300'],
                     capture_output=True, text=True)
print('Fig. 3 mode (correlated-run merge, standalone) exit code:', ret.returncode)
if ret.returncode != 0:
    print(ret.stdout); print(ret.stderr)
    sys.exit(1)
pw_m = pd.read_csv(os.path.join(OUT, 'fig3merge_raw.stats.pairwise.tsv'), sep='\t')
assert pw_m['n_pairs'].eq(3).all(), pw_m[['method_b', 'n_pairs', 'n_clusters']]
assert pw_m['n_clusters'].eq(3).all()
with open(os.path.join(OUT, 'fig3merge_raw.stats.json')) as fh:
    js_m = json.load(fh)
assert js_m['independence']['per_plot_group']['ACT']['n_clusters'] == 3, \
    js_m['independence']['per_plot_group']
assert js_m['effective_sample_merging']['primary_samples']['MDAMB231'] == [
    'MDAMB231c28', 'MDAMB231c8', 'MDAMB231_popp31']
print('OK fig3merge: standalone tests run on 3 effective samples, not 6 datasets')

# B2.2 MCB panel family uses the same merged units (medians verified)
sys.path.insert(0, HERE)
import mcb as _mcb  # noqa: E402
act_fam = next(f for f in _mcb._ploidy_panel_families(merge_df) if f['panel'] == 'ACT')
assert set(act_fam['matrix'].index) == {'MDAMB231', 'TN6', 'TN1'}, act_fam['matrix'].index
assert act_fam['matrix'].loc['MDAMB231', 'ginkgo|10'] == 90.0
assert act_fam['matrix'].loc['TN6', 'ginkgo|10'] == 50.0
assert act_fam['matrix'].loc['TN1', 'ginkgo|10'] == 70.0
print('OK fig3merge: MCB panel uses one row per effective sample (merged medians)')

# B2.3 pooled pipeline path: one pooled row per effective sample x method, the
#      per-run values merged by median, and the pooling provenance recorded.
EVAL_P = os.path.join(HERE, 'scWGS-ploidy-performances-eval.py')
ret = subprocess.run([sys.executable, EVAL_P, '--stats-only', '--no-mcb',
                      '-o', os.path.join(OUT, 'fig3merge'), '--stats-boot', '300'],
                     capture_output=True, text=True, cwd=HERE)
print('Fig. 3 pooled stats (correlated-run merge) exit code:', ret.returncode)
if ret.returncode != 0:
    print(ret.stdout); print(ret.stderr)
    sys.exit(1)
pool_m = pd.read_csv(os.path.join(OUT, 'fig3merge.pooled.long.tsv'), sep='\t')
gink = pool_m[pool_m['method'] == 'ginkgo|10'].set_index('dataset')
assert set(gink.index) == {'MDAMB231', 'TN6', 'TN1'}, gink.index
assert gink.loc['MDAMB231', 'pct_within'] == 90.0
assert gink.loc['TN6', 'pct_within'] == 50.0
assert gink.loc['TN1', 'pct_within'] == 70.0
mda_prov = pool_m[(pool_m['dataset'] == 'MDAMB231')
                  & (pool_m['method'] == 'ginkgo|10')].iloc[0]
assert mda_prov['n_source_rows'] == 3
assert 'MDAMB231c28' in mda_prov['source_datasets']
assert 'MDAMB231_popp31' in mda_prov['source_datasets']
pw_p = pd.read_csv(os.path.join(OUT, 'fig3merge.pooled.stats.pairwise.tsv'), sep='\t')
assert pw_p['n_pairs'].eq(3).all()
with open(os.path.join(OUT, 'fig3merge.pooled.stats.json')) as fh:
    js_p = json.load(fh)
merged_units = {g['unit']: g for g in js_p['pooling']['merged_units']}
assert merged_units['MDAMB231']['n_runs'] == 3 and merged_units['TN6']['n_runs'] == 2
print('OK fig3merge: pooled pipeline merges MDA-MB-231 x3 and TN6 x2 into one unit each')

# B2.4 an ACT sample that is NOT listed in PLOIDY_SAMPLE_DATASETS must be
#      reported and kept as its own unit - never guessed or merged silently.
import importlib.util as _ilu  # noqa: E402
import logging as _logging  # noqa: E402
_spec = _ilu.spec_from_file_location('ploidy_eval_for_merge', EVAL_P)
_mod = _ilu.module_from_spec(_spec)
_spec.loader.exec_module(_mod)
unknown_df = merge_df.copy()
unknown_df.loc[unknown_df['dataset'] == 'TN1 · 50 bp', 'dataset'] = 'NEWSAMPLE · 50 bp'
unknown_df.loc[unknown_df['dataset'] == 'NEWSAMPLE · 50 bp', 'donor'] = 'SRP259526_NEWSAMPLE'
_captured = []


class _Cap(_logging.Handler):
    def emit(self, record):
        _captured.append(record.getMessage())


_root = _logging.getLogger()
_root.addHandler(_Cap())
try:
    _pooled_u, _prov_u, _cols_u, _fb_u = _mod._pool_units(unknown_df, ['donor'])
finally:
    _root.handlers = [h for h in _root.handlers if not isinstance(h, _Cap)]
assert set(_pooled_u['dataset']) == {'MDAMB231', 'TN6', 'SRP259526_NEWSAMPLE'}
assert any('NEWSAMPLE' in m and 'PLOIDY_SAMPLE_DATASETS' in m for m in _captured), _captured
print('OK fig3merge: unlisted ACT sample warned about and kept as its own unit')

# ------------------------------------------------------------------ Part D --
# LaTeX table regression: the booktabs tables written by the two evaluation
# scripts must carry, per comparison, EXACTLY the four statistics
# (n, p, r, 95% CI of r) in this order - nothing more, nothing less.
import shutil

# D.1 - CNV-caller table (scWGS-performances-eval.py --stats-only on the
#       synthetic long TSV of Part A; the table is auto-written as .tex)
EVAL = os.path.join(HERE, 'scWGS-performances-eval.py')
tex1 = os.path.join(OUT, 'fig2tex.stats.pairwise.tex')
with open(long_tsv) as fh:
    ret = subprocess.run([sys.executable, EVAL, '-t', '1',
                          '-o', os.path.join(OUT, 'fig2tex'), '--stats-only',
                          '--stats-boot', '300'],
                         stdin=fh, capture_output=True, text=True, cwd=HERE)
assert ret.returncode == 0, ret.stderr[-2000:]
assert os.path.isfile(tex1), F'MISSING {tex1}'


def _check_tex_columns(text, ncols):
    """Assert that every table row occupies exactly `ncols` columns.

    Rows are the lines ending in the LaTeX row terminator; a
    \\multicolumn{N}{..}{..} cell counts as N columns.
    """
    for ln in text.splitlines():
        s = ln.strip()
        if not s or not s.endswith('\\\\'):
            continue
        if s.startswith('\\') and not s.startswith('\\multicolumn'):
            continue
        cells = 0
        for part in s.split('&'):
            m = re.match(r'\s*\\multicolumn\{(\d+)\}', part)
            cells += int(m.group(1)) if m else 1
        assert cells == ncols, \
            F'LaTeX row carries {cells} columns, expected {ncols}: {s[:100]}'


text = open(tex1).read()
assert text.count(r'\begin{landscape}') == 1 and text.count(r'\end{landscape}') == 1
assert text.count(r'\captionof{table}{') == 1
assert text.count(r'\begin{longtable}{llrrrrrrrr}') == 1
assert text.count(r'\endfirsthead') == 1 and text.count(r'\endhead') == 1
assert text.count(r'\endfoot') == 1
assert text.count(r'\addtocounter{table}{-1}') == 1
assert text.count(r'\multicolumn{4}{c}{Hap\_0}') == 2     # first head + head
assert text.count(r'\multicolumn{4}{c}{Hap\_1}') == 2
assert text.count(r'$n$ & $p$ & $r$ & \qty{95}{\percent} CI') == 4
assert r'\cmidrule(lr){3-6}' in text and r'\cmidrule(lr){7-10}' in text
assert 'Holm $p$' not in text and 'Median diff' not in text and '95\\%' not in text
_check_tex_columns(text, 10)
body = [ln for ln in text.splitlines()
        if ln.startswith('\t\t') and ln.rstrip().endswith('\\')
        and 'multicolumn' not in ln and '$n$ & $p$' not in ln
        and 'toprule' not in ln and 'midrule' not in ln and 'bottomrule' not in ln]
assert body and all(ln.count('&') == 9 for ln in body), \
    'caller table rows must have exactly 10 cells (metric, caller_b, 8 statistics)'
assert len(re.findall(r'\[-?\d+\.\d{2}, -?\d+\.\d{2}\]', text)) >= 2 * len(body) - 4
print(F'OK Part D.1: caller LaTeX table = n, p, r, CI per scenario ({len(body)} rows)')

# D.2 - pooled ploidy table (stat_tests.py ploidy mode -> synthetic pooled
#       donor table -> scWGS-ploidy-performances-eval.py --latex-table)
methods_p = ['aneufinder|10', 'chisel|10', 'ginkgo|10', 'hmmcopy|10', 'scabsolute|inf']
effect_p = {'ginkgo|10': 8.0, 'aneufinder|10': 2.0, 'chisel|10': -3.0,
            'hmmcopy|10': -5.0, 'scabsolute|inf': 0.5}
rows_p = []
for d in [F'donor{i}' for i in range(12)]:
    shock = rng.normal(0, 4.0)
    for m in methods_p:
        rows_p.append({'plot': 'POOLED', 'dataset': F'{d}_ds', 'tool': m.split('|')[0],
                       'max_cn': m.split('|')[1], 'method': m, 'window': 0.5,
                       'n_cells': 300, 'n_cells_finite': 280,
                       'pct_within': float(np.clip(60 + shock + effect_p[m]
                                                   + rng.normal(0, 5), 0, 100)),
                       'failed': False, 'donor': d})
pooled_in = os.path.join(OUT, 'synthetic_pooled_long.tsv')
pd.DataFrame(rows_p).to_csv(pooled_in, sep='\t', index=False)
ret = subprocess.run([sys.executable, STAT, '-i', pooled_in,
                      '-o', os.path.join(OUT, 'fig3pooled'),
                      '--reference', 'ginkgo|10', '--cluster-key', 'none',
                      '--boot', '500'], capture_output=True, text=True)
assert ret.returncode == 0, ret.stderr[-2000:]
pooled_tsv = os.path.join(OUT, 'fig3pooled.stats.pairwise.tsv')
assert os.path.isfile(pooled_tsv)
EVAL_P = os.path.join(HERE, 'scWGS-ploidy-performances-eval.py')
ploidy_out = os.path.join(OUT, 'fig3pooledtex')
os.makedirs(OUT, exist_ok=True)
shutil.copyfile(pooled_tsv, ploidy_out + '.pooled.stats.pairwise.tsv')
ret = subprocess.run([sys.executable, EVAL_P, '-o', ploidy_out, '--latex-table'],
                     capture_output=True, text=True, cwd=HERE)
assert ret.returncode == 0, ret.stderr[-2000:]
text_p = ret.stdout
assert text_p.count(r'\begin{table}[htbp]') == 1
assert r'\caption{' in text_p and r'\label{tab:scwgs-ploidy-pairwise-pooled}' in text_p
assert r'\begin{tabular}{lllrrr}' in text_p
assert (r'Method $b$ & CN cap & $n$ & $p$ & $r$ & \qty{95}{\percent} CI'
        in text_p)
assert 'AneuFinder 2016' in text_p and 'CHISEL 2021' in text_p
assert 'AneuFinder 2016 & 10 &' in text_p
assert '|10' not in text_p and '|inf' not in text_p     # raw ids are rendered
assert 'Median diff' not in text_p and 'Paired $n$' not in text_p
assert '95\\%' not in text_p
_check_tex_columns(text_p, 6)
body_p = [ln for ln in text_p.splitlines()
          if ln.startswith('\t\t') and ln.rstrip().endswith('\\')
          and '$n$ & $p$' not in ln and 'toprule' not in ln
          and 'midrule' not in ln and 'bottomrule' not in ln]
assert body_p and all(ln.count('&') == 5 for ln in body_p), \
    'ploidy table rows must have exactly 6 cells (method_b, CN cap, n, p, r, CI)'
assert any(' & none & ' in ln for ln in body_p), 'uncapped rows must show "none"'
print(F'OK Part D.2: pooled ploidy LaTeX table = method, CN cap, n, p, r, CI '
      F'({len(body_p)} rows)')

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
print()


# ------------------------------------------------------------------ Part F --
# [REV v5] Hsu's MCB (comparison with the best) - bench_results/mcb.py.
# F.1 core: planted matrices (unique winner / clear separation / null) and
#            the k=2 exact reduction to the paired-t interval;
# F.2 CLI: caller mode (Fig. 2) and ploidy mode (Fig. 3) on the synthetic
#            tables of Parts A/B;
# F.3 integration: scWGS-performances-eval.py --stats-only (Part D.1) must
#            also write .stats.mcb.* and the Fig. 2 source data by default;
# F.4 LaTeX regression: exactly n, p, r, CI per row; tectonic compile.
import mcb  # noqa: E402

# F.1a - planted unique winner / clear ranking (k = 5 methods, n = 20 units)
rngF = np.random.default_rng(23)
methF = ['M1', 'M2', 'M3', 'M4', 'M5']
effF = {'M1': 0.30, 'M2': 0.20, 'M3': 0.12, 'M4': 0.02, 'M5': 0.00}
rowsF = []
for u in range(20):
    shock = rngF.normal(0, 0.03)
    for m in methF:
        rowsF.append({'unit': F'u{u}', 'method': m,
                      'value': 0.5 + effF[m] + shock + rngF.normal(0, 0.05)})
matF = pd.DataFrame(rowsF).pivot(index='unit', columns='method', values='value')
resF = mcb.mcb_analyse(matF, alpha=0.05, n_resamples=500, seed=1)
assert resF is not None
vF = resF['verdict']
assert vF['sample_best'] == 'M1'
assert vF['unique_winner'] and vF['leading_group'] == ['M1'], vF
assert 'M5' in vF['inferior'] and 'M4' in vF['inferior'], vF
rF = {r['method']: r for r in resF['rows']}
assert rF['M1']['mcb_low'] > 0                       # the winner beats all others
assert rF['M5']['mcb_high'] < 0 and rF['M5']['pvalue_mcb_one_sided'] <= 0.05
assert rF['M1']['pvalue_mcb_one_sided'] >= 0.5       # the sample-best cannot be 'inferior'
assert rF['M1']['rank_biserial_r_vs_best'] > 0 > rF['M5']['rank_biserial_r_vs_best']
assert all(r['n_units'] == 20 for r in resF['rows'])
print(F'OK F.1a: planted ranking -> unique winner {vF["sample_best"]}, '
      F'g={vF["leading_group_size"]}, {len(vF["inferior"])} method(s) inferior')

# F.1b - planted null (no method difference): family-wise Type I error over
# replicated null families must stay near the nominal level.  Two null
# designs are probed: (i) the benchmark's correlation structure (a shared
# unit shock -> method values positively correlated within a unit), where
# the MCB calibration must be at or below the nominal level; (ii) fully
# independent within-unit values (the worst case for the max-t calibration;
# the bootstrap-t max is known to be mildly liberal at moderate n there, so
# the guard is a documented anticonservatism ceiling, not the nominal level).
n_null_reps, B_null = 60, 300
for label, shared, guard_hi in (('correlated (benchmark structure)', 0.10, None),
                                ('independent (worst case)', 0.0, 0.15)):
    rej_null = 0
    for rep in range(n_null_reps):
        rng_rep = np.random.default_rng(1000 + rep)
        rows_rep = []
        for u in range(20):
            shock = rng_rep.normal(0, shared)
            for m in methF:
                rows_rep.append({'unit': F'u{u}', 'method': m,
                                 'value': 0.5 + shock + rng_rep.normal(0, 0.08)})
        mat_rep = pd.DataFrame(rows_rep).pivot(index='unit', columns='method', values='value')
        res_rep = mcb.mcb_analyse(mat_rep, alpha=0.05, n_resamples=B_null, seed=rep)
        if res_rep is None:
            continue
        rej_null += int(any(r['significant_inferior'] for r in res_rep['rows']))
    rate_null = rej_null / n_null_reps
    lo_null = max(0.0, 0.05 - 3 * np.sqrt(0.05 * 0.95 / n_null_reps))
    hi_null = 0.05 + 3 * np.sqrt(0.05 * 0.95 / n_null_reps)
    if guard_hi is None:
        assert rate_null <= hi_null, \
            F'MCB family-wise Type I error {rate_null:.3f} > {hi_null:.3f} under the correlated null'
    else:
        assert rate_null <= guard_hi, \
            F'MCB family-wise Type I error {rate_null:.3f} > {guard_hi} under the independent null'
    print(F'OK F.1b [{label}]: family-wise rejection in {rej_null}/{n_null_reps} '
          F'families ({100 * rate_null:.0f}%; nominal 5%, guard [{lo_null:.3f}, {hi_null:.3f}])')

# F.1c - k = 2 reduces exactly to the paired-t interval (up to the bootstrap
# Monte Carlo error of the critical value)
d2 = rngF.normal(0.12, 0.05, 25)
mat2 = pd.DataFrame({'A': 0.5 + d2, 'B': 0.5}, index=[F'u{i}' for i in range(25)])
res2 = mcb.mcb_analyse(mat2, alpha=0.05, n_resamples=2000, seed=1)
from scipy import stats as _sps
t_ci = _sps.t.interval(0.95, len(d2) - 1, loc=np.mean(d2), scale=_sps.sem(d2))
rA = [r for r in res2['rows'] if r['method'] == 'A'][0]
tol = 0.10 * (t_ci[1] - t_ci[0])
assert abs(rA['mcb_low'] - t_ci[0]) <= tol and abs(rA['mcb_high'] - t_ci[1]) <= tol, \
    (rA['mcb_low'], rA['mcb_high'], t_ci)
assert rA['mcb_low'] > 0 and [r for r in res2['rows'] if r['method'] == 'B'][0]['mcb_high'] < 0
print(F'OK F.1c: k=2 interval [{rA["mcb_low"]:.4f}, {rA["mcb_high"]:.4f}] matches the '
      F'paired-t CI [{t_ci[0]:.4f}, {t_ci[1]:.4f}] (tol {tol:.4f}); A unique-best, B inferior')

# F.2 - CLI: caller mode on the Part A long TSV
mcb_out = os.path.join(OUT, 'figFmcb')
with open(long_tsv) as fh:
    ret = subprocess.run([sys.executable, os.path.join(HERE, 'mcb.py'),
                          '-t', 'caller', '-o', mcb_out, '--boot', '300'],
                         stdin=fh, capture_output=True, text=True, cwd=HERE)
print('F.2 caller-mode CLI exit code:', ret.returncode)
if ret.returncode != 0:
    print(ret.stdout); print(ret.stderr); sys.exit(1)
mc = pd.read_csv(mcb_out + '.tsv', sep='\t')
need_mcb = ['n_units', 'gap_to_best', 'se_gap', 'mcb_low', 'mcb_high',
            'pvalue_mcb_one_sided', 'rank_biserial_r_vs_best',
            'ci95_r_low', 'ci95_r_high', 'significant_inferior',
            'family_leading_group', 'family_leading_group_size']
missing = [c for c in need_mcb if c not in mc.columns]
assert not missing, F'MCB TSV missing columns: {missing}'
assert mc['n_units'].eq(9).all()          # nine donors, as the pairwise tests
assert set(mc['scenario']) == {'Hap_0', 'Hap_1'}
assert mc['metric'].nunique() == len(metrics)
with open(mcb_out + '.json') as fh:
    jm = json.load(fh)
ginkgo_fams = [v for v in jm['per_family_verdicts'] if v['sample_best'] == 'ginkgo']
assert len(ginkgo_fams) == len(jm['per_family_verdicts'])   # planted best everywhere
assert 'ginkgo' in jm['consensus_per_scenario']['Hap_0']['methods_never_significantly_inferior']
print(F'OK F.2a: caller MCB CLI -> {len(mc)} rows, {len(jm["per_family_verdicts"])} families, '
      F'ginkgo sample-best in every family')

# F.2b - CLI: ploidy mode on the Part B balloon long table
mcb_out_p = os.path.join(OUT, 'figFmcbPloidy')
ret = subprocess.run([sys.executable, os.path.join(HERE, 'mcb.py'),
                      '-t', 'ploidy', '-i', ploidy_tsv,
                      '-o', mcb_out_p, '--boot', '300'],
                     capture_output=True, text=True, cwd=HERE)
print('F.2b ploidy-mode CLI exit code:', ret.returncode)
if ret.returncode != 0:
    print(ret.stdout); print(ret.stderr); sys.exit(1)
for f in [mcb_out_p + '.pooled.stats.mcb.tsv', mcb_out_p + '.pooled.stats.mcb.json',
          mcb_out_p + '.pooled.stats.mcb.tex', mcb_out_p + '.stats.mcb_perpanel.tsv']:
    assert os.path.isfile(f), F'MISSING {f}'
mcp = pd.read_csv(mcb_out_p + '.pooled.stats.mcb.tsv', sep='\t')
mcpp = pd.read_csv(mcb_out_p + '.stats.mcb_perpanel.tsv', sep='\t')
assert set(mcpp['panel']) == {'COLO-829', 'HCC1395', 'HeLa', 'ACT'}
with open(mcb_out_p + '.pooled.stats.mcb.json') as fh:
    jpp = json.load(fh)
best_pooled = jpp['consensus']['pooled']['sample_best']
assert best_pooled.startswith('ginkgo'), best_pooled   # planted ploidy effect
print(F'OK F.2b: ploidy MCB CLI -> pooled best {best_pooled} (n={mcp["n_units"].iloc[0]}), '
      F'per-panel families {sorted(set(mcpp["panel"]))}')

# F.3 - integration: the D.1 --stats-only run must have written the MCB files
# and the Fig. 2 source data by default
for f in ['fig2tex.stats.mcb.tsv', 'fig2tex.stats.mcb.json', 'fig2tex.stats.mcb.tex',
          'fig2tex.fig2_source_data.tsv', 'fig2tex.fig2_source_data.meta.json']:
    p = os.path.join(OUT, f)
    assert os.path.isfile(p), F'MISSING {p} (the --stats-only run of D.1)'
src = pd.read_csv(os.path.join(OUT, 'fig2tex.fig2_source_data.tsv'), sep='\t')
assert len(src) == 300 * 9 * 2 * len(metrics)      # cells x callers x scenarios x metrics
assert src['donor'].nunique() == 9
print(F'OK F.3: eval script writes MCB + Fig. 2 source data by default '
      F'({len(src)} source-data rows, {src["donor"].nunique()} donors)')

# F.4 - LaTeX regression of the two MCB tables
text_m = open(os.path.join(OUT, 'fig2tex.stats.mcb.tex')).read()
assert text_m.count(r'\begin{landscape}') == 1
assert text_m.count(r'\begin{longtable}{llrrrrrrrr}') == 1
assert text_m.count(r'\multicolumn{4}{c}{Hap\_0}') == 2
assert text_m.count(r'\multicolumn{4}{c}{Hap\_1}') == 2
assert text_m.count(r'$n$ & $p$ & $r$ & \qty{95}{\percent} CI') == 4
assert 'Ginkgo' in text_m
assert '95\\%' not in text_m
_check_tex_columns(text_m, 10)
body_m = [ln for ln in text_m.splitlines()
          if ln.startswith('\t\t') and ln.rstrip().endswith('\\')
          and 'multicolumn' not in ln and '$n$ & $p$' not in ln
          and 'toprule' not in ln and 'midrule' not in ln and 'bottomrule' not in ln]
assert body_m and all(ln.count('&') == 9 for ln in body_m), \
    'caller MCB rows must have exactly 10 cells (metric, caller, 8 statistics)'
text_p = open(mcb_out_p + '.pooled.stats.mcb.tex').read()
assert text_p.count(r'\begin{landscape}') == 1
assert text_p.count(r'\begin{longtable}{llrrrr}') == 1
assert text_p.count(r'\multicolumn{6}{l}') == 5          # pooled + 4 panels
assert r'Method & CN cap & $n$ & $p$ & $r$ & \qty{95}{\percent} CI' in text_p
assert '95\\%' not in text_p
_check_tex_columns(text_p, 6)
body_p = [ln for ln in text_p.splitlines()
          if ln.startswith('\t\t') and ln.rstrip().endswith('\\')
          and 'multicolumn' not in ln and '$n$ & $p$' not in ln
          and 'toprule' not in ln and 'midrule' not in ln and 'bottomrule' not in ln]
assert body_p and all(ln.count('&') == 5 for ln in body_p), \
    'ploidy MCB rows must have exactly 6 cells (method, CN cap, n, p, r, CI)'
print(F'OK F.4: MCB LaTeX tables carry exactly n, p, r, CI '
      f'({len(body_m)} caller rows; {len(body_p)} ploidy rows)')

# tectonic compile (optional; skipped when tectonic is absent)
if shutil.which('tectonic'):
    tex_dir = os.path.join(OUT, 'texcheckF')
    os.makedirs(tex_dir, exist_ok=True)
    shutil.copyfile(os.path.join(OUT, 'fig2tex.stats.mcb.tex'),
                    os.path.join(tex_dir, 'caller.tex'))
    shutil.copyfile(mcb_out_p + '.pooled.stats.mcb.tex',
                    os.path.join(tex_dir, 'ploidy.tex'))
    wrap = os.path.join(tex_dir, 'wrap.tex')
    with open(wrap, 'w') as fh:
        fh.write('\\documentclass{article}\n'
                 '\\usepackage{booktabs}\n\\usepackage{longtable}\n'
                 '\\usepackage{pdflscape}\n\\usepackage{caption}\n'
                 '\\usepackage{siunitx}\n'
                 '\\begin{document}\n'
                 '\\input{caller.tex}\\input{ploidy.tex}\\end{document}\n')
    ret = subprocess.run(['tectonic', 'wrap.tex'], cwd=tex_dir,
                         capture_output=True, text=True)
    assert ret.returncode == 0, ret.stderr[-1500:]
    print('OK F.5: both MCB tables compile with tectonic')

print('\nAll stat_tests end-to-end tests PASSED '
      '(including the independence-failure demonstration and Hsu\'s MCB)')
