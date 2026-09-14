#!/usr/bin/env python3
"""stat_tests.py - Statistical tests for CopyNumBench benchmark results.

Design rationale (why these tests)
==================================
The scWGS benchmark is a randomized-complete-block design: every CNV caller is
evaluated on the SAME simulated cells, so per-cell performances are paired
(blocked) by cell. The performance metrics (accuracy, Pearson correlation,
coverage fraction, breakpoint F1-score) are bounded and non-normal, so
nonparametric tests are used throughout, and all comparisons are two-sided
(two-tailed):

1.  Omnibus per (ground-truth scenario, metric): Friedman test (the
    repeated-measures rank ANOVA) across the k callers, with Kendall's W as
    concordance effect size and the per-caller mean ranks.
2.  Post-hoc pairwise: two-sided Wilcoxon signed-rank tests, reference caller
    vs. every other caller, paired per cell (Pratt-style zero handling,
    zero_method='zsplit'), Holm-Bonferroni family-wise correction within each
    (scenario, metric) family of comparisons. All pairwise comparisons are
    available with --all-pairs.
3.  Effect sizes per comparison: matched-pairs rank-biserial correlation r
    (positive = reference performs better) with a 95% percentile-bootstrap
    confidence interval obtained by resampling the independent units (seeded,
    hence fully reproducible), the paired common-language effect
    size P(ref > other) + 0.5*P(ref == other), and the median performance
    difference with a 95% confidence interval.
4.  Hap_0 vs. Hap_1 agreement (the CNP/aneuploidy-fix check): Spearman's rho
    between the per-caller median-performance rankings under the two
    ground-truth scenarios, plus a two-sided Wilcoxon signed-rank test of the
    within-caller Hap_1 - Hap_0 difference.
5.  Ploidy benchmark (Fig. 3): per plot group (COLO-829, HCC1395, HeLa, ACT),
    Friedman test across methods on the per-sample percentage of cells within
    the +/-0.5 window, then two-sided Wilcoxon signed-rank tests paired by
    sample (reference method vs. every other), Holm-corrected, with the same
    effect sizes. McNemar's exact test is provided for per-cell paired binary
    outcomes (within/outside the window) when two tools are evaluated on
    identical evaluable cells.

Independence assumptions (v3 revision - READ THIS)
==================================================
WHAT IS MODELLED: every test is paired/blocked. All k callers are run on the
SAME cells (Fig. 2) and all methods on the SAME datasets (Fig. 3); the
cross-caller correlation within one cell/dataset is exactly what the pairing
(blocking) accounts for. No independent-samples test is used anywhere.

WHAT IS ASSUMED (and was the weak spot of v1): a randomized complete block
design further requires the BLOCKS to be mutually independent - i.e. that
per-cell (Fig. 2) / per-dataset (Fig. 3) results, including the MANY results
produced by the SAME caller across blocks, are independent draws. v1 tested
all ~1,989 per-cell differences as if they were 1,989 independent
observations. For this benchmark that assumption is NOT credible:

* The ~1,989 simulated cells are not ~1,989 independent genomes. They are
  built from the haploid material of only NINE donors (eight oocyte donors +
  one sperm donor) and THREE COSMIC copy-number templates (COLO-829,
  HCC1395, HeLa). The design deliberately preserves the donors' coverage
  fluctuations, amplification bias and mapping artifacts. A caller confused
  by donor-specific material repeats that error on every cell carrying it:
  per-cell results of the same caller - and the paired differences between
  two callers - are positively correlated within donors (with finer
  (accession_1, accession_2, cellLine) sub-groups inside each donor).
* The ITH simulation is nested by construction (CNVs deleted at a lower CNA
  percentage p stay deleted at higher p), so neighbouring p values are
  correlated too.
* Fig. 3: within a germline-derived plot group the "datasets" are
  (donor, sampleType, avgSpotLen) combinations - several datasets share one
  donor's biological material, and every dataset of a group shares the same
  emulated cell-line template.

WHAT HAPPENS IF THE ASSUMPTION FAILS (it does, to a measurable degree):
classic pseudoreplication. The null variance of the Wilcoxon signed-rank /
Friedman statistics assumes independent blocks; with intra-cluster
correlation rho and average cluster size m the variance is inflated by the
design effect DE = 1 + (m - 1) * rho, so the effective sample size is only
n_eff = n / DE. With ~45 shared-material clusters of ~44 cells and a modest
ICC of 0.3, DE ~ 14: every W/chi2 statistic is ~sqrt(14) ~ 3.7x too large, a
difference whose true two-sided P is 0.05 can be reported orders of magnitude
smaller, Holm no longer controls the family-wise error (it corrects
already-anticonservative P values), Kendall's W is inflated, i.i.d.
cell-level bootstrap CIs lose coverage, and in the extreme the benchmark
"winner" can be decided by ONE donor whose artifacts favour a caller (a
cluster-level Simpson-type reversal). bench_results/test_stat_tests.py
(demo_independence_failure) reproduces exactly this failure mode and measures
it: under a true null with ICC = 0.3 the per-cell Wilcoxon rejects in the
large majority of runs while the cluster-level test stays at the nominal
level.

THE v3 FIX (default behaviour)
------------------------------
Inference is moved to the level of the independent experimental unit - the
human DONOR for the Fig. 2 caller benchmark - while per-cell quantities are
kept as DESCRIPTIVE statistics and as flagged naive (non-inferential)
comparisons:

* Fig. 2 default cluster key: `donor`.  All ~1,989 simulated cells are
  downsamplings of the haploid material of only NINE donors, so per-cell
  results of one caller - and the paired differences between callers - are
  correlated within donors.  Each donor is therefore one effective sample,
  and all cells of a donor are aggregated into one observation per caller
  before testing.  Finer keys
  (`--cluster-key accession_1,accession_2,cellLine`) are available as a
  sensitivity analysis; `--cluster-key none` reverts to the naive per-cell
  tests.
* Per comparison, the per-cell differences d = x - y are aggregated to
  per-cluster medians d_g (one value per cluster); the two-sided Wilcoxon
  signed-rank test and the exact two-sided sign test run on the d_g, and
  Holm-Bonferroni is applied to the CLUSTER-level P values.
* The Friedman omnibus likewise runs on per-cluster caller medians (rows =
  clusters); the per-cell Friedman is kept as a 'cell (naive)' row.
* The 95% CI of the median per-cell difference comes from a CLUSTER bootstrap
  (clusters resampled with replacement, all of a cluster's cells kept
  together), not from an i.i.d. cell bootstrap.
* Diagnostics per comparison: ICC(1,1) (one-way ANOVA, method of moments) of
  the paired differences within clusters, the design effect, the effective
  sample size, and the naive-vs-cluster P-value ratio - the degree of
  dependence in the actual data is measured, not assumed away.
* Fig. 3: per plot group the first usable column of donor -> cellLine ->
  dataset defines the clusters (a real-tumor sample is its own unit; a
  missing donor label collapses to one shared '(missing)' cluster, which is
  the conservative choice).

Estimand note: the cluster-level tests target the cluster-population
generalisation ("on a NEW donor, does caller A beat caller B?"), each cluster
(donor) weighted equally; the cell-level medians and
common-language effect sizes describe the benchmarked cell population, with
cluster-robust uncertainty. Both are reported side by side.

Exact P values are written in full to the output TSVs. When the asymptotic
P value underflows double precision, the note column records it; the TSV
carries the raw float.

Outputs (prefix = -o/--output)
==============================
<prefix>.stats.pairwise.tsv    one row per comparison: inference level,
                               n_clusters / n_cells_paired, medians, median
                               difference, W, two-sided P (donor/cluster level
                               by default), sign-test P, Holm-adjusted P,
                               rank-biserial r + its 95% bootstrap CI
                               (ci95_r_low/high), CL effect size, 95% CI of
                               the median difference, naive per-cell P +
                               rank-biserial, ICC, design effect, effective
                               n, P-inflation ratio, notes
<prefix>.stats.friedman.tsv    omnibus Friedman chi2, df, P, Kendall's W,
                               per-caller mean ranks; rows at both levels:
                               'cluster' (primary) and 'cell (naive)'
<prefix>.stats.concordance.tsv (caller benchmark) Hap_0 vs Hap_1 agreement
<prefix>.stats.json            settings, exact n per analysis, the
                               independence/clustering record, versions, seed

Usage
=====
# Fig. 2 (CNV-caller benchmark) statistics, on the long TSV:
cat ${BENCHMARK_RESULT_FILE_PREFIX}.long.tsv | \
    python bench_results/stat_tests.py -o ${BENCHMARK_RESULT_FILE_PREFIX}.stats \
    --reference ginkgo
# ... with a finer, shared-material sensitivity analysis:
cat ${BENCHMARK_RESULT_FILE_PREFIX}.long.tsv | \
    python bench_results/stat_tests.py -o ${BENCHMARK_RESULT_FILE_PREFIX}.stats.fine \
    --reference ginkgo --cluster-key accession_1,accession_2,cellLine

# Fig. 3 (ploidy benchmark) statistics, on the balloon-plot table:
python bench_results/stat_tests.py -i '*_pct_within_long.tsv' \
    -o scWGS-ploidy-performances.stats --reference 'ginkgo|10'
"""
from __future__ import annotations

import argparse
import glob
import json
import logging
import os
import platform
import sys
import warnings

import numpy as np
import pandas as pd

try:
    import scipy
    from scipy import stats as sps
    _HAVE_SCIPY = True
except ImportError:  # pragma: no cover - the frozen env always has scipy
    _HAVE_SCIPY = False

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(filename)s %(levelname)s %(message)s')

# Columns that identify one simulated cell in the long TSV. (accession_1,
# accession_2, cellLine) already determine the simulated cell (data3from2.py
# seeds every simulation decision from exactly these inputs); overall_ploidy
# and CNA_percent are derived from the same seed and are kept as a safety net.
DEFAULT_CELL_KEY = ['accession_1', 'accession_2', 'cellLine',
                    'overall_ploidy', 'CNA_percent']

# Default CLUSTER key for Fig. 2 (v3): donor.  In the current benchmark all
# simulated cells of one donor are downsamplings of that donor's haploid
# material, so they are not independent draws and each human donor is one
# effective sample.  (accession_1, accession_2, cellLine) identifies a single
# simulated cell, not a reusable cluster, in the exported long TSV - each
# accession-pair x template combination occurs once.  Pass
# --cluster-key accession_1,accession_2,cellLine for a finer sensitivity
# analysis, or --cluster-key none for the (discouraged) naive per-cell tests.
DEFAULT_CLUSTER_KEY = ['donor']

# Fig. 3: columns tried, in order, to find the independent unit within each
# plot group (first column with >= 2 distinct non-missing values wins; a
# real-tumor sample with no donor metadata ends up clustered by dataset,
# which is the honest choice for samples from different patients).
PLOIDY_CLUSTER_FALLBACK_CHAIN = ['donor', 'cellLine', 'dataset']

# Cluster bootstrap runtime cap (the pairwise loop can run >100 comparisons).
MAX_CLUSTER_BOOTSTRAP = 5000

TINY = float(np.finfo(float).tiny)

# --------------------------------------------------------------------------- #
# Core, dependency-light statistics helpers                                   #
# --------------------------------------------------------------------------- #
def holm_bonferroni(pvals):
    """Holm-Bonferroni family-wise adjusted p-values (step-down)."""
    p = np.asarray(pvals, dtype=float)
    m = len(p)
    adj = np.empty(m, dtype=float)
    order = np.argsort(p, kind='stable')
    running_max = 0.0
    for rank, i in enumerate(order):
        val = (m - rank) * p[i]
        running_max = max(running_max, val)
        adj[i] = min(1.0, running_max)
    return adj


def rank_biserial_matched(d):
    """Matched-pairs rank-biserial correlation for paired differences d.

    r = (R+ - R-) / (n(n+1)/2), Pratt-consistent: zero differences contribute
    half of their rank to each side (matching zero_method='zsplit').
    r > 0 means x tends to be larger than y (reference better when d = x - y).
    """
    d = np.asarray(d, dtype=float)
    d = d[np.isfinite(d)]
    n = len(d)
    if n == 0:
        return float('nan')
    ranks = sps.rankdata(np.abs(d))
    zero = d == 0
    r_pos = float(np.sum(ranks[d > 0]) + 0.5 * np.sum(ranks[zero]))
    r_neg = float(np.sum(ranks[d < 0]) + 0.5 * np.sum(ranks[zero]))
    total = n * (n + 1.0) / 2.0
    return (r_pos - r_neg) / total


def common_language_paired(d):
    """P(d > 0) + 0.5*P(d == 0): probability that the reference scores higher
    than the competitor on a random cell (ties count half)."""
    d = np.asarray(d, dtype=float)
    d = d[np.isfinite(d)]
    if len(d) == 0:
        return float('nan')
    return float(np.mean(d > 0) + 0.5 * np.mean(d == 0))


def wilcoxon_signed_rank(x, y, alternative='two-sided', zero_method='zsplit'):
    """Two-sided Wilcoxon signed-rank test on paired samples.

    Returns a dict with n_pairs, n_zero, n_nonzero, statistic, pvalue and an
    optional note. Complete pairs only (both values finite).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.shape != y.shape:
        raise ValueError(F'wilcoxon_signed_rank: shape mismatch {x.shape} vs {y.shape}')
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    d = x - y
    out = {
        'n_pairs': int(len(d)),
        'n_zero': int(np.sum(d == 0)),
        'n_nonzero': int(np.sum(d != 0)),
        'statistic': float('nan'),
        'pvalue': float('nan'),
        'note': '',
    }
    if len(d) == 0:
        out['note'] = 'no complete pairs'
        return out
    if np.all(d == 0):
        out['statistic'] = 0.0
        out['pvalue'] = 1.0
        out['note'] = 'all paired differences are zero'
        return out
    try:
        with warnings.catch_warnings():
            # small cluster counts make scipy warn about the normal
            # approximation while still returning a valid exact/asymptotic P
            warnings.simplefilter('ignore', UserWarning)
            res = sps.wilcoxon(x, y, zero_method=zero_method, alternative=alternative,
                               correction=True, method='auto')
        out['statistic'] = float(res.statistic)
        out['pvalue'] = float(res.pvalue)
    except (ValueError, RuntimeWarning) as exc:  # degenerate tiny samples
        out['note'] = F'Wilcoxon undefined on this sample ({exc}); P set to 1'
        out['pvalue'] = 1.0
        return out
    if out['pvalue'] <= TINY:
        out['note'] = F'asymptotic two-sided P underflows double precision (P < {TINY:.1e})'
    return out


def sign_test_two_sided(d):
    """Exact two-sided sign test on paired differences (zeros dropped).

    Complements the Wilcoxon signed-rank test when the number of independent
    units is too small for the exact signed-rank null to reach P < 0.05.
    Returns dict(n, k_positive, pvalue)."""
    d = np.asarray(d, dtype=float)
    d = d[np.isfinite(d)]
    pos = int(np.sum(d > 0))
    neg = int(np.sum(d < 0))
    n = pos + neg
    if n == 0:
        return {'n': 0, 'k_positive': pos, 'pvalue': 1.0}
    p = min(1.0, 2.0 * sps.binom.cdf(min(pos, neg), n, 0.5))
    return {'n': n, 'k_positive': pos, 'pvalue': float(p)}


def bca_bootstrap_ci(d, statistic=np.median, n_resamples=10000,
                     confidence_level=0.95, seed=1):
    """BCa bootstrap CI of a statistic of paired differences (seeded).

    NAIVE-MODE ONLY (assumes i.i.d. observations): used when clustering is
    disabled (--cluster-key none). With clustered data use
    cluster_bootstrap_median_ci instead.
    Falls back to the percentile method for degenerate samples, then to the
    point estimate itself.
    """
    d = np.asarray(d, dtype=float)
    d = d[np.isfinite(d)]
    if len(d) == 0:
        return float('nan'), float('nan'), 'no data'
    point = float(statistic(d))
    if np.all(d == d[0]):
        return point, point, 'degenerate sample (all differences equal)'
    rng = np.random.default_rng(seed)
    for method in ('BCa', 'percentile'):
        try:
            # scipy < 1.15 names the seed argument `random_state`; it accepts a
            # np.random.Generator (frozen env: scipy 1.14.1) as well as an int.
            # Degenerate inputs (e.g. a constant metric) raise or warn in scipy;
            # we catch/silence that and fall through to the next method, ending
            # at the point estimate, so a degenerate comparison never aborts
            # the whole benchmark statistics run.
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                res = sps.bootstrap((d,), statistic, n_resamples=int(n_resamples),
                                    confidence_level=confidence_level, method=method,
                                    random_state=rng)
            lo, hi = float(res.confidence_interval.low), float(res.confidence_interval.high)
            if np.isfinite(lo) and np.isfinite(hi):
                return lo, hi, method
            logging.debug('bootstrap (%s) returned a non-finite interval', method)
        except Exception as exc:  # degenerate jackknife / too few distinct values
            logging.debug('bootstrap (%s) failed: %s', method, exc)
    return point, point, 'bootstrap degenerate; point estimate returned'


def cluster_bootstrap_median_ci(d, cluster_ids, statistic=np.median,
                                n_resamples=2000, confidence_level=0.95, seed=1):
    """Cluster (hierarchical) bootstrap CI of the median of clustered d.

    Clusters - the independent experimental units - are resampled with
    replacement and all observations of a drawn cluster are kept together;
    the median is computed over the pooled resampled observations. This
    targets the cell-population median while respecting within-cluster
    correlation, unlike an i.i.d. cell bootstrap, whose intervals are too
    narrow under clustering (coverage below the nominal level).
    Returns (low, high, method)."""
    d = np.asarray(d, dtype=float)
    cl = np.asarray(list(cluster_ids), dtype=object)
    if len(d) != len(cl):
        raise ValueError('cluster_bootstrap_median_ci: d and cluster_ids length mismatch')
    ok = np.isfinite(d)
    d, cl = d[ok], cl[ok]
    if len(d) == 0:
        return float('nan'), float('nan'), 'no data'
    point = float(statistic(d))
    uniq, inv = np.unique(cl, return_inverse=True)
    g = int(len(uniq))
    if g < 2:
        return point, point, 'single cluster: CI undefined'
    if np.all(d == d[0]):
        return point, point, 'degenerate sample (all differences equal)'
    order = np.argsort(inv, kind='stable')
    d_sorted, inv_sorted = d[order], inv[order]
    bounds = list(np.searchsorted(inv_sorted, np.arange(g), side='left')) + [len(inv_sorted)]
    groups = [d_sorted[bounds[i]:bounds[i + 1]] for i in range(g)]
    rng = np.random.default_rng(seed)
    b = int(n_resamples)
    meds = np.empty(b, dtype=float)
    for i in range(b):
        pick = rng.integers(0, g, g)
        meds[i] = statistic(np.concatenate([groups[j] for j in pick]))
    a = 1.0 - float(confidence_level)
    lo, hi = np.quantile(meds, [a / 2.0, 1.0 - a / 2.0])
    return float(lo), float(hi), 'cluster-percentile'


def bootstrap_r_ci(units, n_resamples=2000, confidence_level=0.95, seed=1):
    """Percentile-bootstrap CI of the matched-pairs rank-biserial effect size r.

    `units` are the i.i.d. paired differences AT THE INFERENCE LEVEL: the
    per-cluster (donor) medians in cluster mode, or the raw per-cell
    differences in naive mode. The independent units are resampled with
    replacement, r is recomputed on every resample, and the percentile
    interval of the r distribution is returned. This targets the same estimand
    as the reported rank_biserial_r, so (r, CI) describe ONE effect size.
    Degenerate samples fall back to the point estimate with an explanatory
    method string. Returns (low, high, method).
    """
    u = np.asarray(units, dtype=float)
    u = u[np.isfinite(u)]
    if len(u) == 0:
        return float('nan'), float('nan'), 'no data'
    point = rank_biserial_matched(u)
    if len(u) < 2:
        return point, point, 'single unit: CI undefined'
    if np.all(u == u[0]):
        return point, point, 'degenerate sample (all differences equal)'
    rng = np.random.default_rng(seed)
    b = int(n_resamples)
    rs = np.empty(b, dtype=float)
    for i in range(b):
        pick = rng.integers(0, len(u), len(u))
        rs[i] = rank_biserial_matched(u[pick])
    a = 1.0 - float(confidence_level)
    lo, hi = np.quantile(rs, [a / 2.0, 1.0 - a / 2.0])
    return float(lo), float(hi), 'unit-percentile'


def icc_design_effect(d, cluster_ids):
    """ICC(1,1) of clustered observations + design effect + effective n.

    One-way random-effects ANOVA, method of moments:
        ICC = (MSB - MSW) / (MSB + (m0 - 1) * MSW),
    with the unbalanced-size adjustment
        m0 = (N - sum(n_g^2)/N) / (G - 1).
    Design effect DE = 1 + (m0 - 1) * max(ICC, 0) (a negative ICC estimate is
    reported as-is but treated as 0 for DE); effective sample size
    n_eff = N / DE. This quantifies, on the actual data, how badly a per-cell
    (naive) test would overstate significance. Returns a dict."""
    d = np.asarray(d, dtype=float)
    cl = np.asarray(list(cluster_ids), dtype=object)
    if len(d) != len(cl):
        raise ValueError('icc_design_effect: d and cluster_ids length mismatch')
    ok = np.isfinite(d)
    d, cl = d[ok], cl[ok]
    n = int(len(d))
    out = {'n_obs': n, 'n_clusters': 0, 'icc': float('nan'), 'm0': float('nan'),
           'design_effect': float('nan'), 'n_effective': float('nan')}
    uniq, inv = np.unique(cl, return_inverse=True)
    g = int(len(uniq))
    out['n_clusters'] = g
    if g < 2 or n <= g or n < 3:
        return out
    if np.all(d == d[0]):
        return out  # zero variance: ICC undefined
    cnt = np.bincount(inv, minlength=g).astype(float)
    sums = np.bincount(inv, weights=d, minlength=g)
    means = sums / cnt
    grand = float(d.mean())
    ssb = float(np.sum(cnt * (means - grand) ** 2))
    ssw = float(np.sum(d ** 2) - np.sum(cnt * means ** 2))
    msb = ssb / (g - 1.0)
    msw = ssw / (n - g)
    m0 = (n - float(np.sum(cnt ** 2)) / n) / (g - 1.0)
    denom = msb + (m0 - 1.0) * msw
    icc = float('nan')
    if denom > 0:
        icc = float((msb - msw) / denom)
    de = float(1.0 + (m0 - 1.0) * max(icc, 0.0)) if np.isfinite(icc) else float('nan')
    out.update({'icc': icc, 'm0': float(m0), 'design_effect': de,
                'n_effective': float(n / de) if np.isfinite(de) and de > 0 else float('nan')})
    return out


def friedman_test(samples):
    """Friedman test on complete blocks. samples: list of k equal-length arrays.

    Returns dict(chi2, df, pvalue, n_blocks, k, kendalls_w, mean_ranks).
    NOTE: the blocks are assumed mutually independent - with clustered blocks
    (repeated results from the same caller/material) run this on per-cluster
    aggregated values, not on raw per-cell values.
    """
    mat = np.column_stack([np.asarray(s, dtype=float) for s in samples])
    n, k = mat.shape
    if k < 3 or n < 2:
        return None
    res = sps.friedmanchisquare(*[mat[:, j] for j in range(k)])
    mean_ranks = sps.rankdata(mat, axis=1).mean(axis=0)
    return {
        'chi2': float(res.statistic),
        'df': int(k - 1),
        'pvalue': float(res.pvalue),
        'n_blocks': int(n),
        'k': int(k),
        'kendalls_w': float(res.statistic / (n * (k - 1))),
        'mean_ranks': [float(r) for r in mean_ranks],
    }


def mcnemar_exact(within_a, within_b):
    """Exact McNemar test on paired binary outcomes (e.g. per-cell
    within/outside the ploidy window for two tools evaluated on the same
    cells). Two-sided exact binomial on the discordant pairs.
    Returns dict(b, c, n_discordant, pvalue).

    CAUTION (v3): like every per-cell test, this assumes independent cells;
    with cells clustered by donor/material, aggregate the binary outcome per
    cluster (e.g. the majority/median outcome per donor) or restrict the test
    to cells of a single dataset before using it."""
    a = np.asarray(within_a, dtype=bool)
    b = np.asarray(within_b, dtype=bool)
    if a.shape != b.shape:
        raise ValueError('mcnemar_exact: paired vectors must have equal length')
    b01 = int(np.sum(~a & b))   # only B inside
    c10 = int(np.sum(a & ~b))   # only A inside
    n = b01 + c10
    if n == 0:
        return {'b': b01, 'c': c10, 'n_discordant': 0, 'pvalue': 1.0}
    p = min(1.0, 2.0 * sps.binom.cdf(min(b01, c10), n, 0.5))
    return {'b': b01, 'c': c10, 'n_discordant': n, 'pvalue': float(p)}


# --------------------------------------------------------------------------- #
# Cluster bookkeeping                                                         #
# --------------------------------------------------------------------------- #
def _norm_missing(v):
    """Normalise a metadata value to '' when missing, else to its str form."""
    if v is None:
        return ''
    if isinstance(v, float) and np.isnan(v):
        return ''
    try:
        if pd.isna(v):
            return ''
    except (TypeError, ValueError):
        pass
    return str(v).strip()


def _cluster_ids(frame, cols):
    """Per-row cluster id string built from the given columns.

    Missing components render as '' (empty), so rows with missing metadata
    form identifiable groups instead of crashing."""
    if not cols:
        return pd.Series([''] * len(frame), index=frame.index, dtype=object)
    parts = [frame[c].map(_norm_missing) for c in cols]
    cid = parts[0].astype(object)
    for p in parts[1:]:
        cid = cid + '|' + p.astype(object)
    return cid


def _as_key(t):
    """Hashable key from a (possibly 1-element) index entry."""
    if isinstance(t, (tuple, list)):
        return tuple(t)
    return (t,)


def _cell2cluster_map(df, key_cols, cluster_cols):
    """dict: cell key (tuple of key_cols values) -> cluster id.

    The first occurrence wins; cell-level metadata is identical across the
    caller rows of one cell, so this is well defined."""
    cid = _cluster_ids(df, cluster_cols)
    m = {}
    for k, v in zip((_as_key(t) for t in df[key_cols].itertuples(index=False, name=None)),
                    cid.tolist()):
        m.setdefault(k, v)
    return m


def _pivot_cluster_labels(piv, cell2cluster):
    """Cluster id per row of a cell-keyed pivot, plus a duplicate-safe index."""
    labels = [cell2cluster.get(_as_key(t), '__UNMAPPED__') for t in piv.index]
    if any(l == '__UNMAPPED__' for l in labels):
        logging.warning('stat_tests: some pivot rows have no cluster label '
                        '(cell key not found in the metadata map); they form '
                        "their own '__UNMAPPED__' group")
    return pd.Series(labels, index=piv.index, dtype=object)


def software_versions():
    info = {
        'python': platform.python_version(),
        'numpy': np.__version__,
        'pandas': pd.__version__,
    }
    if _HAVE_SCIPY:
        info['scipy'] = scipy.__version__
    return info


# --------------------------------------------------------------------------- #
# Shared pairwise-record builder (cluster-level inference + naive comparison)  #
# --------------------------------------------------------------------------- #
def _pairwise_record(x, y, labels, base, cluster_agg='median', n_resamples=10000,
                     seed=1, cluster_key_str=''):
    """One pairwise-comparison record.

    Primary inference on the CLUSTER level (per-cluster/donor aggregation of
    the per-cell differences), naive per-cell comparison kept for
    transparency.
    x, y: paired per-cell (or per-dataset) values, equal length; labels:
    cluster id per element (None -> naive mode). base: dict with the identity
    columns (scenario/metric/caller_a/caller_b or plot/method_a/method_b).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if labels is not None and len(labels) != len(x):
        raise ValueError('_pairwise_record: labels length does not match x')
    ok = np.isfinite(x) & np.isfinite(y)
    if labels is not None:
        labels = [l for l, m in zip(labels, ok) if m]
    x, y = x[ok], y[ok]
    d = x - y
    rec = dict(base)
    w_cell = wilcoxon_signed_rank(x, y)
    notes = [n for n in (w_cell['note'],) if n]
    rec.update({
        'inference_level': 'cluster' if labels is not None else 'cell (naive)',
        'cluster_key': cluster_key_str if labels is not None else '',
        'n_cells_paired': w_cell['n_pairs'],
        'median_a': float(np.median(x)) if len(x) else float('nan'),
        'median_b': float(np.median(y)) if len(y) else float('nan'),
        'median_diff_a_minus_b': float(np.median(d)) if len(d) else float('nan'),
        'mean_diff_a_minus_b': float(np.mean(d)) if len(d) else float('nan'),
        'cl_effect_paired': common_language_paired(d),
        'rank_biserial_r_cell_naive': rank_biserial_matched(d),
        'wilcoxon_W_cell_naive': w_cell['statistic'],
        'pvalue_cell_naive': w_cell['pvalue'],
        'n_zero_diffs': w_cell['n_zero'],
    })
    if labels is not None:
        cl_ser = pd.Series(d, index=pd.Index(labels, dtype=object))
        cm = cl_ser.groupby(level=0).agg(cluster_agg).sort_index()
        cm_v = cm.to_numpy(dtype=float)
        w_cl = wilcoxon_signed_rank(cm_v, np.zeros_like(cm_v))
        st = sign_test_two_sided(cm_v)
        n_boot = int(min(n_resamples, MAX_CLUSTER_BOOTSTRAP))
        ci_lo, ci_hi, ci_method = cluster_bootstrap_median_ci(
            d, labels, np.median, n_resamples=max(n_boot, 200), seed=seed)
        r_lo, r_hi, r_ci_method = bootstrap_r_ci(
            cm_v, n_resamples=max(n_boot, 200), seed=seed)
        icc = icc_design_effect(d, labels)
        rec.update({
            'n_pairs': int(len(cm_v)),            # units the primary test runs on
            'n_clusters': int(len(cm_v)),
            'wilcoxon_W': w_cl['statistic'],
            'pvalue_two_sided': w_cl['pvalue'],
            'pvalue_sign_test': st['pvalue'],
            'rank_biserial_r': rank_biserial_matched(cm_v),
            'ci95_r_low': r_lo, 'ci95_r_high': r_hi, 'ci_r_method': r_ci_method,
            'ci95_median_diff_low': ci_lo, 'ci95_median_diff_high': ci_hi,
            'ci_method': ci_method,
            'icc_within_cluster_d': icc['icc'],
            'design_effect': icc['design_effect'],
            'n_effective_cells': icc['n_effective'],
            'mean_cluster_size_m0': icc['m0'],
        })
        if w_cl['note']:
            notes.append(w_cl['note'])
        if 0 < len(cm_v) < 6:
            notes.append(F'only {len(cm_v)} independent clusters: the exact '
                         'two-sided Wilcoxon cannot reach P < 0.05 below 6 '
                         'units; see pvalue_sign_test')
        if np.isfinite(w_cl['pvalue']) and np.isfinite(w_cell['pvalue']):
            if w_cl['pvalue'] < 1e-12 and w_cell['pvalue'] < 1e-12:
                rec['p_inflation_ratio'] = float('nan')  # both underflow
                notes.append('both naive and cluster P underflow; ratio undefined')
            else:
                # > 1 means the naive per-cell P overstates significance
                rec['p_inflation_ratio'] = float(
                    w_cl['pvalue'] / max(w_cell['pvalue'], TINY))
        else:
            rec['p_inflation_ratio'] = float('nan')
    else:
        st = sign_test_two_sided(d)
        ci_lo, ci_hi, ci_method = bca_bootstrap_ci(
            d, np.median, n_resamples=n_resamples, seed=seed)
        r_lo, r_hi, r_ci_method = bootstrap_r_ci(
            d, n_resamples=min(n_resamples, MAX_CLUSTER_BOOTSTRAP), seed=seed)
        rec.update({
            'n_pairs': w_cell['n_pairs'],
            'n_clusters': float('nan'),
            'wilcoxon_W': w_cell['statistic'],
            'pvalue_two_sided': w_cell['pvalue'],
            'pvalue_sign_test': st['pvalue'],
            'rank_biserial_r': rank_biserial_matched(d),
            'ci95_r_low': r_lo, 'ci95_r_high': r_hi, 'ci_r_method': r_ci_method,
            'ci95_median_diff_low': ci_lo, 'ci95_median_diff_high': ci_hi,
            'ci_method': ci_method,
            'icc_within_cluster_d': float('nan'),
            'design_effect': float('nan'),
            'n_effective_cells': float('nan'),
            'mean_cluster_size_m0': float('nan'),
            'p_inflation_ratio': 1.0,
        })
        notes.append('NAIVE per-cell level: per-cell results treated as '
                     'independent (pseudoreplication risk); rerun with a '
                     '--cluster-key for valid inference')
    rec['note'] = '; '.join(notes)
    return rec


# --------------------------------------------------------------------------- #
# Fig. 2: CNV-caller benchmark statistics (long TSV)                           #
# --------------------------------------------------------------------------- #
def run_caller_benchmark_stats(df, out_prefix, perf_metrics, gamete_type2short,
                               reference='ginkgo', all_pairs=False,
                               pair_key_cols=None, cluster_key_cols=None,
                               cluster_agg='median', n_resamples=10000,
                               seed=1, alpha=0.05):
    """Statistical tests for the main CNV-caller benchmark (Fig. 2).

    df: long TSV dataframe (one row per caller per simulated cell);
    perf_metrics: metric ids; gamete_type2short: scenario prefix -> short tag;
    cluster_key_cols: None -> DEFAULT_CLUSTER_KEY (donor, i.e. each human
    donor is one effective sample); [] -> naive per-cell mode; a list ->
    custom cluster key (missing columns are dropped with a warning).
    Writes <out_prefix>.stats.{pairwise,friedman,concordance}.{tsv,json}.
    Returns the dict that is also written to .stats.json.
    """
    if not _HAVE_SCIPY:
        logging.warning('scipy is not available: statistical tests skipped')
        return None
    key_cols = [c for c in (pair_key_cols or DEFAULT_CELL_KEY) if c in df.columns]
    if not key_cols:
        logging.warning('no cell-identity columns (%s) found in the input; '
                        'falling back to order-based pairing (callers must '
                        'have equal row counts)', ', '.join(pair_key_cols or DEFAULT_CELL_KEY))
        counts = df.groupby('Caller').size()
        if counts.nunique() != 1:
            raise SystemExit(F'stat_tests: cannot pair callers with unequal row counts {counts.to_dict()}')
    if key_cols:
        dup = df.duplicated(subset=key_cols + ['Caller']).sum()
        if dup:
            logging.warning('%d duplicate (cell, caller) rows: aggregating by median', dup)
    callers = sorted(df['Caller'].unique())
    if reference not in callers:
        raise SystemExit(F'stat_tests: reference caller {reference!r} not among {callers}')

    # ---- resolve the cluster key (the independent experimental unit) ----
    requested_cluster = (cluster_key_cols is None
                         or list(cluster_key_cols) != [])
    if cluster_key_cols is None:
        cluster_cols = list(DEFAULT_CLUSTER_KEY)
        cluster_key_source = 'default (donor = independent human donor)'
    else:
        cluster_cols = list(cluster_key_cols)
        cluster_key_source = 'user-specified'
    dropped = [c for c in cluster_cols if c not in df.columns]
    cluster_cols = [c for c in cluster_cols if c in df.columns]
    if dropped:
        logging.warning('cluster-key columns %s not present in the input; '
                        'clustering uses the remaining %s',
                        dropped, cluster_cols or 'nothing')
    can_cluster = bool(cluster_cols) and bool(key_cols)
    if requested_cluster and not can_cluster:
        logging.warning('cluster analysis impossible (%s): falling back to the '
                        'NAIVE per-cell level, whose independence assumption is '
                        'questionable for this benchmark (see the module docstring)',
                        'no cluster columns' if key_cols else 'no cell-key columns')
    cell2cluster = _cell2cluster_map(df, key_cols, cluster_cols) if can_cluster else None
    cluster_key_str = '|'.join(cluster_cols) if can_cluster else ''

    independence = {
        'pairing': 'all callers evaluated on the same simulated cells; '
                   'per-cell metrics paired (blocked) by cell',
        'remaining_assumption_of_blocks': 'blocks (cells) mutually independent - '
                   'i.e. per-cell results of the SAME caller treated as independent',
        'assumption_violation': 'all ~1,989 cells are downsamplings of the haploid '
                   'material of only NINE donors scored against three COSMIC templates; '
                   'cells of one donor share that donor\'s BAM-derived noise and artifacts, '
                   'and ITH deletions are nested across CNA_percent, so same-caller '
                   'results are positively correlated within donors',
        'cluster_key': cluster_cols,
        'cluster_key_source': cluster_key_source,
        'aggregation': F'per-cluster {cluster_agg} of the per-cell paired '
                       'differences (donor level by default)',
        'inference_level': ('cluster (independent experimental units)'
                            if can_cluster else
                            'cell (NAIVE - independence assumed; pseudoreplication risk)'),
        'primary_tests': ('two-sided Wilcoxon signed-rank + exact sign test on '
                          'per-independent-unit (donor by default) medians; Friedman '
                          'on the same per-unit caller medians; Holm-Bonferroni on '
                          'unit-level P; unit-cluster bootstrap 95% CI'
                          if can_cluster else
                          'two-sided Wilcoxon signed-rank per cell; Friedman per cell; '
                          'BCa 95% CI (all NAIVE)'),
        'naive_columns_kept_for_comparison': ['pvalue_cell_naive',
                                              'rank_biserial_r_cell_naive',
                                              'wilcoxon_W_cell_naive'],
        'diagnostics': 'ICC(1,1) of the paired differences within clusters (one-way '
                       'ANOVA, method of moments); design effect DE = 1 + (m0-1)*ICC; '
                       'n_effective = n_cells / DE; p_inflation_ratio = cluster P / naive P',
        'recommended_sensitivity': 'rerun with --cluster-key '
                                   'accession_1,accession_2,cellLine (finer, '
                                   'shared-material level) and compare the conclusions',
    }
    if can_cluster:
        cid_counts = pd.Series(list(cell2cluster.values()), dtype=object).value_counts()
        independence['n_clusters'] = int(len(cid_counts))
        independence['cluster_sizes'] = {
            'min': int(cid_counts.min()), 'median': float(np.median(cid_counts)),
            'max': int(cid_counts.max()), 'total_cells': int(cid_counts.sum())}
        if len(cid_counts) < 10:
            logging.warning('only %d independent clusters detected: statistical power '
                            'is limited and conclusions should be phrased cautiously',
                            len(cid_counts))

    settings = {
        'analysis': 'scWGS CNV-caller benchmark (Fig. 2)',
        'design': ('randomized complete block; blocks = independent clusters '
                   '(%s, donor by default); per-cell metrics paired by cell and '
                   'aggregated per cluster (donor) before testing' % cluster_key_str
                   if can_cluster else
                   'randomized complete block; per-cell performances paired by '
                   'cell (NAIVE: cells treated as independent)'),
        'pair_key_columns': key_cols or ['<row order>'],
        'independence': independence,
        'omnibus_test': 'Friedman per (scenario, metric) at the cluster level '
                        "('cluster' rows) and at the per-cell level "
                        "('cell (naive)' rows)",
        'posthoc_test': ('two-sided Wilcoxon signed-rank + exact sign test on '
                         'per-cluster (donor-level by default) medians of the '
                         'paired differences '
                         "(zero_method='zsplit', continuity correction)"
                         if can_cluster else
                         "two-sided Wilcoxon signed-rank paired per cell "
                         "(zero_method='zsplit', continuity correction)"),
        'multiple_comparison_correction': 'Holm-Bonferroni within each (scenario, metric) '
                        'family on the cluster-level (donor-level by default) P values '
                        '(scenario-difference rows: '
                        'Holm within each (metric) family across callers)',
        'effect_sizes': ['matched-pairs rank-biserial r (cluster/donor level) with '
                         'a 95% percentile-bootstrap CI over the independent units',
                         'paired common-language effect size (cell population)',
                         'median per-cell difference with cluster-bootstrap 95% CI'],
        'tail': 'two-sided (two-tailed) for all tests',
        'reference_caller': reference,
        'all_pairs': bool(all_pairs),
        'bootstrap_n_resamples': int(n_resamples),
        'cluster_bootstrap_n_resamples': int(min(n_resamples, MAX_CLUSTER_BOOTSTRAP)),
        'bootstrap_seed': int(seed),
        'alpha': float(alpha),
        'n_cells_per_caller': {c: int(n) for c, n in df.groupby('Caller').size().items()},
        'software': software_versions(),
    }

    pairwise_rows, friedman_rows, concordance_rows, scen_diff_rows = [], [], [], []

    def _fr_row(level, mat, n_cells, tag, metric):
        complete = mat.dropna(axis=0, how='any')
        cols = list(mat.columns)
        row = {'level': level, 'scenario': tag, 'metric': metric,
               'n_blocks': int(len(complete)), 'n_cells_aggregated': int(n_cells),
               'k_callers': int(len(cols))}
        fr = friedman_test([complete[c].to_numpy() for c in cols]) if len(complete) else None
        if fr is not None:
            row.update({'friedman_chi2': fr['chi2'], 'df': fr['df'],
                        'pvalue': fr['pvalue'], 'kendalls_w': fr['kendalls_w'],
                        **{F'mean_rank_{c}': r for c, r in zip(cols, fr['mean_ranks'])}})
        else:
            row['note'] = 'omnibus undefined (needs >= 3 callers and >= 2 complete blocks)'
        return row

    for gamete_type, tag in gamete_type2short.items():
        for metric in perf_metrics:
            col = F'{gamete_type}.{metric}'
            if col not in df.columns:
                continue
            # ---- per-cell pivot: rows = simulated cells, columns = callers ----
            if key_cols:
                series = df.groupby(key_cols + ['Caller'])[col].median()
                piv = series.unstack('Caller')
            else:
                piv = df.pivot(index=None, columns='Caller', values=col)
                piv.index = np.arange(len(piv))
            piv = piv.reindex(columns=callers)

            lab_series = _pivot_cluster_labels(piv, cell2cluster) if cell2cluster else None

            # ---- omnibus Friedman: cluster level (primary) + cell level (naive) ----
            if lab_series is not None:
                cmat = pd.DataFrame(piv.to_numpy(dtype=float),
                                    index=pd.Index(lab_series.to_numpy(), dtype=object),
                                    columns=piv.columns)
                cmat = cmat.groupby(level=0).agg(cluster_agg).sort_index()
                n_cells_behind = int(piv.notna().any(axis=1).sum())
                friedman_rows.append(_fr_row('cluster', cmat, n_cells_behind, tag, metric))
            friedman_rows.append(_fr_row('cell (naive)', piv, int(piv.notna().any(axis=1).sum()),
                                         tag, metric))

            # ---- pairwise post-hoc ----
            others = [c for c in callers if c != reference]
            pairs = ([(a, b) for a in callers for b in callers if a < b]
                     if all_pairs else [(reference, b) for b in others])
            records = []
            for a, b in pairs:
                sub = piv[[a, b]].dropna(axis=0, how='any')
                if lab_series is not None:
                    labels = lab_series.loc[sub.index].tolist()
                else:
                    labels = None
                records.append(_pairwise_record(
                    sub[a].to_numpy(dtype=float), sub[b].to_numpy(dtype=float),
                    labels, {'scenario': tag, 'metric': metric,
                             'caller_a': a, 'caller_b': b},
                    cluster_agg=cluster_agg, n_resamples=n_resamples, seed=seed,
                    cluster_key_str=cluster_key_str))
            # Holm family = all comparisons within (scenario, metric)
            if records:
                adj = holm_bonferroni([r['pvalue_two_sided'] for r in records])
                for r, p_adj in zip(records, adj):
                    r['pvalue_holm'] = float(p_adj)
                    r['reject_holm'] = bool(p_adj <= alpha)
                pairwise_rows.extend(records)

    # ---- Hap_0 vs Hap_1 concordance of the caller rankings ----
    # caller medians per scenario (descriptive ranking, units = callers)
    def caller_medians(gamete_prefix):
        out = {}
        for c in callers:
            sub = df[df['Caller'] == c]
            out[c] = [float(sub[F'{gamete_prefix}.{m}'].median())
                      if F'{gamete_prefix}.{m}' in sub.columns else float('nan')
                      for m in perf_metrics]
        return pd.DataFrame(out, index=perf_metrics)
    prefixes = list(gamete_type2short.keys())
    if len(prefixes) == 2:
        med_a, med_b = caller_medians(prefixes[0]), caller_medians(prefixes[1])
        for metric in perf_metrics:
            if metric not in med_a.index or metric not in med_b.index:
                continue
            va, vb = med_a.loc[metric], med_b.loc[metric]
            ok = np.isfinite(va) & np.isfinite(vb)
            if ok.sum() >= 3:
                rho = sps.spearmanr(va[ok], vb[ok])
                concordance_rows.append({
                    'metric': metric, 'n_callers': int(ok.sum()),
                    'spearman_rho': float(rho.statistic),
                    'pvalue_two_sided': float(rho.pvalue),
                    'scenario_a': gamete_type2short[prefixes[0]],
                    'scenario_b': gamete_type2short[prefixes[1]],
                    'note': 'descriptive rank concordance across callers '
                            '(callers are the units; algorithmic families are '
                            'not strictly independent)',
                })
        # within-caller scenario difference (the CNP/aneuploidy-fix effect)
        for c in callers:
            sub = df[df['Caller'] == c]
            for metric in perf_metrics:
                ca, cb = F'{prefixes[0]}.{metric}', F'{prefixes[1]}.{metric}'
                if ca not in sub.columns or cb not in sub.columns:
                    continue
                if key_cols:
                    subg = sub.groupby(key_cols)[[ca, cb]].median()
                    x = subg[ca].to_numpy(dtype=float)
                    y = subg[cb].to_numpy(dtype=float)
                    labels = ([cell2cluster.get(_as_key(t), '__UNMAPPED__')
                               for t in subg.index] if cell2cluster else None)
                else:
                    x = sub[ca].to_numpy(dtype=float)
                    y = sub[cb].to_numpy(dtype=float)
                    labels = None
                if len(x) == 0 or x.shape != y.shape:
                    continue
                scen_diff_rows.append(_pairwise_record(
                    x, y, labels,
                    {'scenario': F'{gamete_type2short[prefixes[0]]}_vs_{gamete_type2short[prefixes[1]]}',
                     'metric': metric, 'caller_a': c, 'caller_b': '<scenario>'},
                    cluster_agg=cluster_agg, n_resamples=n_resamples, seed=seed,
                    cluster_key_str=cluster_key_str))
        # Holm family for the scenario-difference rows: (metric) across callers
        by_metric = {}
        for r in scen_diff_rows:
            by_metric.setdefault(r['metric'], []).append(r)
        for metric, rows in by_metric.items():
            adj = holm_bonferroni([r['pvalue_two_sided'] for r in rows])
            for r, p_adj in zip(rows, adj):
                r['pvalue_holm'] = float(p_adj)
                r['reject_holm'] = bool(p_adj <= alpha)
        pairwise_rows.extend(scen_diff_rows)
    # clean leftovers
    for r in pairwise_rows:
        r.setdefault('pvalue_holm', float('nan'))
        r.setdefault('reject_holm', '')

    _write_tables(out_prefix, pairwise_rows, friedman_rows, concordance_rows,
                  settings, ('scenario', 'metric', 'caller_a', 'caller_b'))
    return settings


def _write_tables(out_prefix, pairwise_rows, friedman_rows, concordance_rows,
                  settings, pairwise_sort_cols):
    out_dir = os.path.dirname(os.path.abspath(out_prefix))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    pd.DataFrame(pairwise_rows).sort_values(list(pairwise_sort_cols)).to_csv(
        F'{out_prefix}.stats.pairwise.tsv', sep='\t', index=False, na_rep='NA')
    pd.DataFrame(friedman_rows).to_csv(
        F'{out_prefix}.stats.friedman.tsv', sep='\t', index=False, na_rep='NA')
    if concordance_rows:
        pd.DataFrame(concordance_rows).to_csv(
            F'{out_prefix}.stats.concordance.tsv', sep='\t', index=False, na_rep='NA')
    with open(F'{out_prefix}.stats.json', 'w') as fh:
        json.dump(settings, fh, indent=2)
    logging.info('stat_tests: wrote %s.stats.pairwise.tsv, .stats.friedman.tsv%s, .stats.json',
                 out_prefix, ', .stats.concordance.tsv' if concordance_rows else '')


# --------------------------------------------------------------------------- #
# Fig. 3: ploidy-benchmark statistics (balloon-plot table)                     #
# --------------------------------------------------------------------------- #
def run_ploidy_benchmark_stats(tab, out_prefix, reference='ginkgo|10',
                               cluster_key_cols=None, cluster_agg='median',
                               n_resamples=10000, seed=1, alpha=0.05):
    """Statistical tests for the ploidy-estimation benchmark (Fig. 3).

    tab: the *_pct_within_long.tsv table written by
    scWGS-ploidy-performances-eval(.v02).py with columns plot, dataset, tool,
    max_cn, method, pct_within, n_cells_finite, failed (plus donor /
    sampleType / avgSpotLen / cellLine when produced by v02). Pairing is by
    dataset (sample) within each plot group.

    cluster_key_cols: None -> per plot group the first usable column of
    donor -> cellLine -> dataset defines the independent units (datasets of
    one germline donor share that donor's material, so per-method results on
    them are correlated); [] -> naive per-dataset mode; a list -> custom chain.
    Writes the same .stats.* files.
    """
    if not _HAVE_SCIPY:
        logging.warning('scipy is not available: statistical tests skipped')
        return None
    needed = {'plot', 'dataset', 'method', 'pct_within'}
    missing = needed - set(tab.columns)
    if missing:
        raise SystemExit(F'stat_tests: ploidy table lacks columns {sorted(missing)}')
    tab = tab.copy()
    if 'failed' in tab.columns:
        tab = tab[~tab['failed'].astype(bool)]
    tab['pct_within'] = pd.to_numeric(tab['pct_within'], errors='coerce')
    methods = sorted(tab['method'].astype(str).unique())
    if reference not in methods:
        raise SystemExit(F'stat_tests: reference method {reference!r} not among {methods}')
    n_cells_col = 'n_cells_finite' if 'n_cells_finite' in tab.columns else None

    if cluster_key_cols is None:
        chain = list(PLOIDY_CLUSTER_FALLBACK_CHAIN)
        chain_source = F'default fallback chain: {" -> ".join(chain)}'
    elif list(cluster_key_cols) == []:
        chain = []
        chain_source = 'clustering disabled (--cluster-key none): NAIVE per-dataset level'
    else:
        chain = [c for c in cluster_key_cols if c in tab.columns]
        dropped = [c for c in cluster_key_cols if c not in tab.columns]
        if dropped:
            logging.warning('cluster-key columns %s not present in the ploidy table; '
                            'using the remaining %s', dropped, chain or 'nothing')
        chain_source = 'user-specified chain: ' + (' -> '.join(chain) or '(none)')

    per_group_info = {}

    def _pick_cluster_column(sub, plot):
        """First column of the chain with >= 2 distinct non-missing values."""
        for cand in chain:
            if cand not in sub.columns:
                continue
            vals = sub[cand].map(_norm_missing)
            uniq = sorted(set(v for v in vals.unique() if v != ''))
            if len(uniq) >= 2:
                return cand, vals
            logging.info('plot %s: cluster column %r is missing or constant; '
                         'trying the next column', plot, cand)
        return None, None

    settings = {
        'analysis': 'ploidy-estimation benchmark (Fig. 3)',
        'unit': 'per-sample percentage of cells with inferred ploidy within the window',
        'design': 'randomized complete block; percentages paired by dataset (sample) '
                  'within each plot group; inference aggregated to the independent '
                  'unit (cluster) selected per plot group',
        'independence': {
            'pairing': 'all methods evaluated on the same datasets; percentages '
                       'paired by dataset within each plot group',
            'remaining_assumption_of_blocks': 'datasets (blocks) mutually independent - '
                       'i.e. per-dataset results of the SAME method treated as independent',
            'assumption_violation': 'within a germline-derived plot group, datasets with '
                       'the same donor share that donor haplotype material, and every '
                       'dataset of the group shares the emulated cell-line template; '
                       'per-method results on shared-donor datasets are correlated',
            'cluster_selection': chain_source,
            'missing_donor_policy': "datasets without donor metadata collapse into one "
                                    "shared '(missing)' cluster per plot group (conservative)",
            'per_plot_group': per_group_info,
            'aggregation': F'per-cluster {cluster_agg} of the per-dataset paired differences',
            'diagnostics': 'ICC(1,1) of the paired differences within clusters; design '
                           'effect; n_effective; p_inflation_ratio = cluster P / naive P',
        },
        'omnibus_test': 'Friedman per plot group across methods at the cluster level '
                        "('cluster' rows) and per-dataset level ('dataset (naive)' rows)",
        'posthoc_test': 'two-sided Wilcoxon signed-rank + exact sign test on per-cluster '
                        'medians, paired by dataset; Holm-Bonferroni within each plot group',
        'effect_sizes': ['matched-pairs rank-biserial r (cluster level) with a '
                         '95% percentile-bootstrap CI over the independent units',
                         'paired common-language effect size (dataset population)',
                         'median per-dataset difference (percentage points) with '
                         'cluster-bootstrap 95% CI'],
        'tail': 'two-sided (two-tailed) for all tests',
        'reference_method': reference,
        'bootstrap_n_resamples': int(n_resamples),
        'cluster_bootstrap_n_resamples': int(min(n_resamples, MAX_CLUSTER_BOOTSTRAP)),
        'bootstrap_seed': int(seed),
        'alpha': float(alpha),
        'per_method': {},
        'software': software_versions(),
    }

    pairwise_rows, friedman_rows = [], []

    def _fr_row(level, mat, n_datasets, plot):
        complete = mat.dropna(axis=0, how='any')
        cols = list(mat.columns)
        row = {'level': level, 'plot': str(plot),
               'n_blocks': int(len(complete)),
               'n_datasets_aggregated': int(n_datasets),
               'k_methods': int(len(cols))}
        fr = friedman_test([complete[c].to_numpy() for c in cols]) if len(complete) else None
        if fr is not None:
            row.update({'friedman_chi2': fr['chi2'], 'df': fr['df'],
                        'pvalue': fr['pvalue'], 'kendalls_w': fr['kendalls_w'],
                        **{F'mean_rank_{m}': r for m, r in zip(cols, fr['mean_ranks'])}})
        else:
            row['note'] = 'omnibus undefined (needs >= 3 methods and >= 2 complete blocks)'
        return row

    for plot, sub in tab.groupby('plot', sort=True):
        chosen, donor_vals = _pick_cluster_column(sub, plot)
        piv = sub.pivot(index='dataset', columns='method', values='pct_within')
        piv = piv.reindex(columns=[m for m in methods if m in piv.columns])
        ok_methods = [m for m in piv.columns if piv[m].notna().any()]
        if chosen:
            # dataset -> cluster id (missing values collapse to '(missing)')
            d2c = {}
            for ds, v in zip(sub['dataset'], donor_vals):
                d2c.setdefault(ds, v if v != '' else '(missing)')
            lab_series = pd.Series([d2c.get(ds, '(missing)') for ds in piv.index],
                                   index=piv.index, dtype=object)
            n_clusters = int(lab_series.nunique())
            per_group_info[str(plot)] = {'cluster_column': chosen,
                                         'n_clusters': n_clusters}
            if n_clusters < 10:
                logging.warning('plot %s: only %d independent clusters (%s); '
                                'power is limited', plot, n_clusters, chosen)
        else:
            lab_series = None
            per_group_info[str(plot)] = {'cluster_column': None,
                                         'note': 'no usable cluster column; '
                                                 'per-dataset (naive) level only'}

        # ---- omnibus Friedman: cluster level (primary) + dataset level (naive) ----
        if lab_series is not None:
            cmat = pd.DataFrame(piv.to_numpy(dtype=float),
                                index=pd.Index(lab_series.to_numpy(), dtype=object),
                                columns=piv.columns)
            cmat = cmat.groupby(level=0).agg(cluster_agg).sort_index()
            friedman_rows.append(_fr_row('cluster', cmat, int(piv.notna().any(axis=1).sum()), plot))
        friedman_rows.append(_fr_row('dataset (naive)', piv, int(piv.notna().any(axis=1).sum()), plot))

        # ---- pairwise vs reference ----
        records = []
        for m in ok_methods:
            if m == reference:
                continue
            both = piv[[reference, m]].dropna(axis=0, how='any')
            if lab_series is not None:
                labels = lab_series.loc[both.index].tolist()
            else:
                labels = None
            rec = _pairwise_record(
                both[reference].to_numpy(dtype=float), both[m].to_numpy(dtype=float),
                labels, {'plot': str(plot), 'method_a': reference, 'method_b': m},
                cluster_agg=cluster_agg, n_resamples=n_resamples, seed=seed,
                cluster_key_str=chosen or '')
            if n_cells_col and n_cells_col in sub.columns:
                tot = sub[sub['method'].isin([reference, m])].groupby('method')[n_cells_col].sum()
                rec['n_cells_evaluated_a'] = int(tot.get(reference, 0))
                rec['n_cells_evaluated_b'] = int(tot.get(m, 0))
            records.append(rec)
        if records:
            adj = holm_bonferroni([r['pvalue_two_sided'] for r in records])
            for r, p_adj in zip(records, adj):
                r['pvalue_holm'] = float(p_adj)
                r['reject_holm'] = bool(p_adj <= alpha)
            pairwise_rows.extend(records)
        for m in ok_methods:
            msub = sub[sub['method'] == m]
            entry = {'n_datasets': int(len(msub))}
            if n_cells_col and n_cells_col in msub.columns:
                entry['n_cells_total'] = int(pd.to_numeric(
                    msub[n_cells_col], errors='coerce').fillna(0).sum())
            settings['per_method'][m] = entry

    for r in pairwise_rows:
        r.setdefault('pvalue_holm', float('nan'))
        r.setdefault('reject_holm', '')

    _write_tables(out_prefix, pairwise_rows, friedman_rows, [], settings,
                  ('plot', 'method_a', 'method_b'))
    return settings


# --------------------------------------------------------------------------- #
# Command-line interface                                                       #
# --------------------------------------------------------------------------- #
def _expand_inputs(patterns):
    paths, seen = [], set()
    for pat in patterns:
        matched = sorted(glob.glob(pat)) if any(c in pat for c in '*?[') else [pat]
        for p in matched:
            ap = os.path.abspath(p)
            if ap not in seen and os.path.isfile(ap):
                seen.add(ap)
                paths.append(p)
    return paths


def _parse_cluster_key(raw):
    """None | 'none' | 'a,b,c' -> None (default) | [] (naive) | list of cols."""
    if raw is None:
        return None
    s = raw.strip().lower()
    if s in ('none', 'naive', 'off', 'no'):
        return []
    return [c.strip() for c in raw.split(',') if c.strip()]


def main(argv=None):
    parser = argparse.ArgumentParser(
        description='Statistical tests for CopyNumBench benchmark results '
                    '(cluster-level Friedman omnibus; two-sided Wilcoxon '
                    'signed-rank + exact sign test post-hoc with Holm '
                    'correction; rank-biserial and common-language effect '
                    'sizes; cluster-bootstrap CIs; ICC / design-effect '
                    'diagnostics of the independence assumption).')
    parser.add_argument('-i', '--input', nargs='*', default=None,
                        help='Input table(s): either one long TSV (Fig. 2, one '
                             'row per caller per cell) or one *_pct_within_long.tsv '
                             '(Fig. 3). Globs allowed. Default: stdin (Fig. 2 mode).')
    parser.add_argument('-o', '--output', default='bench-results',
                        help='Output prefix; writes <prefix>.stats.{pairwise,friedman}.tsv '
                             'and <prefix>.stats.json (+ .stats.concordance.tsv in Fig. 2 mode).')
    parser.add_argument('--reference', default='ginkgo',
                        help="Reference caller (Fig. 2 mode, default 'ginkgo') or method "
                             "(Fig. 3 mode, e.g. 'ginkgo|10' for CapAt10, default in that mode).")
    parser.add_argument('--all-pairs', action='store_true',
                        help='Fig. 2 mode: test all caller pairs, not only reference vs rest.')
    parser.add_argument('--pair-key', default=None,
                        help='Comma-separated columns identifying one simulated cell '
                             '(Fig. 2 mode). Default: accession_1,accession_2,cellLine,'
                             'overall_ploidy,CNA_percent.')
    parser.add_argument('--cluster-key', default=None,
                        help='Comma-separated columns defining the independent '
                             'experimental unit (cluster). Fig. 2 default: '
                             'donor (each human donor is one effective sample; '
                             'per-cell results of the same donor are aggregated '
                             'to a median before testing); Fig. 3 default chain: '
                             'donor -> cellLine -> dataset. Pass "none" to disable '
                             'clustering and revert to the naive per-cell/per-dataset '
                             'tests (discouraged: pseudoreplication).')
    parser.add_argument('--cluster-agg', choices=['median', 'mean'], default='median',
                        help='Aggregation of per-cell/per-dataset differences within '
                             'each cluster (default median).')
    parser.add_argument('--boot', type=int, default=10000,
                        help='Bootstrap resamples for the CIs (default 10000; the '
                             'cluster bootstrap is capped at 5000).')
    parser.add_argument('--seed', type=int, default=1,
                        help='Seed of the bootstrap RNG (default 1).')
    parser.add_argument('--alpha', type=float, default=0.05,
                        help='Family-wise significance level for Holm rejection (default 0.05).')
    args = parser.parse_args(argv)

    if not _HAVE_SCIPY:
        raise SystemExit('stat_tests: scipy is required (pip install scipy)')

    if args.input:
        paths = _expand_inputs(args.input)
        if not paths:
            raise SystemExit(F'stat_tests: no input files matched {args.input}')
        dfs = [pd.read_csv(p, sep='\t') for p in paths]
        df = pd.concat(dfs, ignore_index=True) if len(dfs) > 1 else dfs[0]
        source = F'{len(paths)} file(s): ' + ', '.join(paths)
    else:
        df = pd.read_csv(sys.stdin, sep='\t')
        source = '<stdin>'

    if 'pct_within' in df.columns:
        reference = args.reference if args.reference != 'ginkgo' else 'ginkgo|10'
        settings = run_ploidy_benchmark_stats(
            df, args.output, reference=reference,
            cluster_key_cols=_parse_cluster_key(args.cluster_key),
            cluster_agg=args.cluster_agg,
            n_resamples=args.boot, seed=args.seed, alpha=args.alpha)
    elif 'Caller' in df.columns:
        settings = run_caller_benchmark_stats(
            df, args.output,
            perf_metrics=['intCN_accuracy', 'intCN_PCC', 'nonintCN_PCC', 'CN_genome_cov_frac',
                          'breakpoint_precision', 'breakpoint_recall', 'breakpoint_f1score'],
            gamete_type2short={'with_haploidy_assumed_gametes': 'Hap_0',
                               'with_aneuploidy_aware_gametes': 'Hap_1'},
            reference=args.reference, all_pairs=args.all_pairs,
            pair_key_cols=(args.pair_key.split(',') if args.pair_key else None),
            cluster_key_cols=_parse_cluster_key(args.cluster_key),
            cluster_agg=args.cluster_agg,
            n_resamples=args.boot, seed=args.seed, alpha=args.alpha)
    else:
        raise SystemExit('stat_tests: input has neither a Caller column (Fig. 2 mode) '
                         'nor a pct_within column (Fig. 3 mode); cannot detect the analysis type')
    if settings:
        settings['input'] = source
        with open(F'{args.output}.stats.json', 'w') as fh:
            json.dump(settings, fh, indent=2)
    return 0


if __name__ == '__main__':
    sys.exit(main())
