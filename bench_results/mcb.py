#!/usr/bin/env python3
"""mcb.py - Hsu's Multiple Comparison with the Best (MCB) for the scWGS
benchmark figures.

WHAT THIS MODULE ADDS
=====================
stat_tests.py answers, per (scenario, metric) or plot group, "does the
REFERENCE caller outperform caller b?" with exactly four statistics
(n, p, r, 95% CI).  The manuscript questions, however, are reference-free:
which caller is the best overall, is there a UNIQUE winner, is there a
leading (top-g) group, and by how much does every other caller trail the
best?  Hsu's Multiple Comparison with the Best (MCB) is the classical tool
for exactly those questions: for every method i it builds a SIMULTANEOUS
confidence interval for

    theta_i - max_{j != i} theta_j        ("method i versus the best of the
                                           others"; larger = better)

so that, with family-wise confidence 1 - alpha over the whole family:

    * interval entirely BELOW 0   ->  i is significantly INFERIOR to the best
                                      (it can be excluded from the leading
                                      group);
    * interval entirely ABOVE 0   ->  i is significantly BETTER than every
                                      other method: i is the UNIQUE WINNER
                                      (possible only for the sample-best
                                      method);
    * interval CONTAINS 0         ->  i is INDISTINGUISHABLE from the best;
                                      i belongs to the leading group.

The MCB family is defined PER TASK (the figure or subfigure the numbers
back), PER METRIC and PER CALLER, exactly mirroring how the benchmark is
reported:

    Fig. 2 (caller benchmark)   one MCB family per (ground-truth scenario,
                                performance metric); rows = the k callers.
    Fig. 3 (ploidy benchmark)   one POOLED family over all donors of all
                                panels (primary, matching the pooled pairwise
                                table) plus one family per balloon panel
                                a-d (per-subfigure sensitivity analysis);
                                rows = the tool|cap methods.
    Fig. 4 (HG008 karyotype)    descriptive by design (the manuscript states
                                the HG008 validation is interpreted as
                                descriptive, not inferential); no MCB is run,
                                but the heatmap source data are exported by
                                cnv_clustermap.py.
    Fig. 5 (scRNA-seq)           handled by the companion repository's own
                                mcb.py (one family per swarm-grid metric).

DESIGN (blocked / repeated-measures MCB, cluster-robust)
=======================================================
All k methods are evaluated on the same cells / datasets, so performance is
paired (blocked) by unit.  Inference is run at the level of the INDEPENDENT
experimental unit - the DONOR for Fig. 2 (the ~1,989 simulated cells are
downsamplings of only nine donors), the pooled donor for Fig. 3 - exactly as
in stat_tests.py: per-unit medians of the per-cell values first, MCB on the
resulting n x k complete-block matrix afterwards (units with any missing
method are dropped and counted; methods missing everywhere are dropped from
the family).

For the complete-block matrix X (n units x k methods, higher = better):

1.  Point estimates: mu_i = mean over units; the "best competitor" of i is
    j*(i) = argmax_{j != i} mu_j and the gap is D_i = mu_i - mu_{j*(i)}.
2.  Per pair (i, j): paired per-unit differences d_ij with mean D_ij and
    standard error SE_ij = sd(d_ij) / sqrt(n).
3.  Simultaneous MCB intervals (PRIMARY, bootstrap): resample the n units
    with replacement B times (seeded, reproducible); T*_ij = (D*_ij -
    D_ij) / SE_ij; c_boot = (1 - alpha) quantile of max_{i<j} |T*_ij|.
    The MCB interval of method i is the Tukey-style projection

        L_i = min_{j != i} (D_ij - c_boot * SE_ij)
        U_i = min_{j != i} (D_ij + c_boot * SE_ij)

    (because theta_i - max_{j != i} theta_j = min_{j != i} (theta_i -
    theta_j), simultaneous pairwise bands project onto the min).  For k = 2
    this reduces exactly to the paired-t interval of the single difference.
4.  MCB-adjusted one-sided P per method (the table's p): the single-step
    bootstrap max-|t| adjusted two-sided P of method i versus its best
    competitor, converted to the one-sided H1 "i is inferior to the best"
    (H0: theta_i >= max_{j != i} theta_j).  Small p = significantly
    inferior; p >= 0.5 for the sample-best method by construction.
5.  Effect size (the table's r): matched-pairs rank-biserial correlation
    of the per-unit differences d_{i,j*(i)} against the best competitor
    (positive = i tends to score higher), with a per-method percentile
    bootstrap 95% CI (descriptive; NOT simultaneous - the simultaneous
    quantity of this analysis is the MCB interval of the gap, which is the
    CI shown in the tables).  The rank-biserial is computed by
    stat_tests.rank_biserial_matched so this module can never disagree with
    the pairwise tables on a shared statistic.
6.  Parametric cross-check: c_par = q_{k, n-1}(1 - alpha) / sqrt(2), the
    studentized-range critical value scaled so that k = 2 reproduces the
    paired-t critical value exactly; intervals [min(D_ij - c_par * SE_ij),
    min(D_ij + c_par * SE_ij)] under normality of the per-unit differences.
    Reported in the TSVs as mcb_low_param / mcb_high_param.
7.  Omnibus: the Friedman test (via stat_tests.friedman_test) on the same
    complete-block matrix - "does any ordering of the methods exist?" -
    with Kendall's W; the MCB verdicts (unique winner / leading group of
    size g / inferior set) then describe the ordering's separability.

Verdicts (per family):
    unique_winner     exactly one method whose interval lies entirely
                      above 0 (equivalently: the leading group has size 1)
    leading_group     every method whose interval is not entirely below 0
                      (the methods that cannot be excluded as best); g is
                      its size - g = 2 is the "top-2" situation
    inferior          methods with interval entirely below 0
    no_separation     no method is significantly inferior (g = k)

LATeX DISCIPLINE (same as every other table of this project)
============================================================
Every table row carries, after its identity columns, EXACTLY the four
statistics n (effective sample size = number of independent units),
p (MCB-adjusted one-sided P), r (rank-biserial effect size versus the best
competitor) and the 95% CI (the simultaneous MCB interval of the gap to the
best, on the metric's own scale) - nothing else.

CALIBRATION NOTE (regression-tested in test_stat_tests.py, Part F.1b)
=====================================================================
The bootstrap max-|t| calibration is a TRUE studentized bootstrap (every
replicate recomputes its own per-pair standard error), so the family-wise
error of the MCB verdicts under a true null is: at or below the nominal
level for the benchmark's positively-correlated within-unit structure at
n = 9 donors (conservative), near-nominal at n ~ 40 materials, and mildly
liberal (up to ~10% at a nominal 5%) in the worst case of independent
within-unit values at moderate n - the price of the nonparametric,
selection-aware calibration.  For k = 2 the interval reduces exactly to the
paired-t interval (up to the bootstrap Monte Carlo error).

Outputs (prefix = -o/--output)
==============================
Caller mode (Fig. 2):
    <prefix>.stats.mcb.tsv          one row per (scenario, metric, caller)
    <prefix>.stats.mcb.tex          booktabs table (two scenario groups,
                                    mirroring the pairwise table layout)
    <prefix>.stats.mcb.json         settings, per-family verdicts, consensus
Ploidy mode (Fig. 3):
    <prefix>.pooled.stats.mcb.tsv   pooled family (one row per method)
    <prefix>.stats.mcb_perpanel.tsv one row per (panel, method)
    <prefix>.stats.mcb.tex          booktabs table: pooled section + one
                                    section per balloon panel
    <prefix>.stats.mcb.json         settings, verdicts, consensus

Usage
=====
# Fig. 2 MCB, standalone on the long TSV:
cat ${PREFIX}.long.tsv | python bench_results/mcb.py -t caller \
    -o ${PREFIX}.plots.stats.mcb
# Fig. 3 MCB, standalone on the balloon long table:
python bench_results/mcb.py -t ploidy -i '<output>_pct_within_long.tsv' \
    -o ${PLOIDY_PREFIX}.stats.mcb --reference-not-needed
The two evaluation scripts call the same runners (run_caller_mcb /
run_ploidy_mcb) by default, so the tables and the figures can never drift
apart.
"""
from __future__ import annotations

import argparse
import glob
import json
import logging
import os
import sys

import numpy as np
import pandas as pd

try:
    import scipy
    from scipy import stats as sps
    _HAVE_SCIPY = True
except ImportError:  # pragma: no cover - the frozen env always has scipy
    _HAVE_SCIPY = False

# Everything numerical is reused from the validated sibling module so the
# pairwise tables, the winner-style analyses and MCB can never disagree on a
# shared statistic (rank-biserial r, Friedman, Holm, ...).
import stat_tests as _ST

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(filename)s %(levelname)s %(message)s')

MAX_BOOTSTRAP = 5000          # runtime cap for the unit bootstrap
DEFAULT_N_RESAMPLES = 2000    # MCB bootstrap resamples (seeded, reproducible)


# --------------------------------------------------------------------------- #
# Formatting helpers (identical style to the eval scripts' tables)            #
# --------------------------------------------------------------------------- #
def _tex_escape(s):
    """Escape the LaTeX specials of a plain identifier."""
    return (str(s).replace('\\', r'\textbackslash{}').replace('&', r'\&')
            .replace('%', r'\%').replace('_', r'\_').replace('#', r'\#')
            .replace('{', r'\{').replace('}', r'\}').replace('$', r'\$'))


def _fmt_effect(x, decimals=2):
    """Signed effect size at fixed decimals; 'n/a' for NaN."""
    if x is None or not np.isfinite(x):
        return 'n/a'
    return F'{x:+.{decimals}f}'


def _fmt_ci(low, high, decimals=2):
    """[low, high] at fixed decimals; 'n/a' when undefined."""
    if low is None or high is None or not (np.isfinite(low) and np.isfinite(high)):
        return 'n/a'
    return F'[{low:.{decimals}f}, {high:.{decimals}f}]'


def _fmt_pvalue_mcb(x, alpha=0.05):
    """p column: full precision, '< 10^-4' guard, '<= alpha' bolded marker."""
    if x is None or not np.isfinite(x):
        return 'n/a'
    if x < 1e-4:
        return r'$< 10^{-4}$'
    return F'{x:.3f}'


def _fmt_signed_int(x):
    """n column: plain integer."""
    if x is None or not np.isfinite(x):
        return 'n/a'
    return F'{int(x)}'


# --------------------------------------------------------------------------- #
# Core: Hsu's MCB on a complete-block unit x method matrix                    #
# --------------------------------------------------------------------------- #
def mcb_analyse(mat, alpha=0.05, n_resamples=DEFAULT_N_RESAMPLES, seed=1):
    """MCB analysis of one family.

    mat: DataFrame, rows = independent units, columns = methods, values
    oriented so that LARGER IS BETTER.  Rows with any missing value are
    dropped (complete blocks); columns that are entirely missing are
    dropped.  Returns a dict:

        {'n_units', 'k_methods', 'units_dropped', 'methods_dropped',
         'friedman' (dict|None), 'c_boot', 'c_par',
         'rows': [per-method record dict], 'verdict': {...}, 'notes': [...]}

    or None when the family is not analysable (fewer than 2 units or
    fewer than 2 methods after the complete-block filter).
    """
    if not _HAVE_SCIPY:
        logging.warning('scipy is not available: MCB analysis skipped')
        return None
    X = mat.copy()
    methods = list(X.columns)
    notes = []
    all_missing = [m for m in methods if not X[m].notna().any()]
    if all_missing:
        X = X.drop(columns=all_missing)
        notes.append(F'methods with no values dropped: {", ".join(map(str, all_missing))}')
    methods = list(X.columns)
    before = len(X)
    X = X.dropna(axis=0, how='any')
    units_dropped = before - len(X)
    if units_dropped:
        notes.append(F'{units_dropped} of {before} units dropped (missing at least '
                     F'one method; complete blocks only)')
    if len(X) < 2 or len(methods) < 2:
        notes.append(F'not analysable: {len(X)} complete units x {len(methods)} methods')
        return None
    X = X.astype(float)
    n, k = X.shape
    vals = X.to_numpy(dtype=float)
    mu = vals.mean(axis=0)
    sd = vals.std(axis=0, ddof=1) if n > 1 else np.zeros(k)

    # ---- per-pair paired-difference statistics ------------------------------
    pairs = [(i, j) for i in range(k) for j in range(i + 1, k)]
    D = np.zeros((k, k))       # D[i, j] = mean(x_i - x_j); D[j, i] = -D[i, j]
    SE = np.zeros((k, k))
    DD = {}
    for i, j in pairs:
        d = vals[:, i] - vals[:, j]
        D[i, j] = d.mean()
        D[j, i] = -D[i, j]
        se = (d.std(ddof=1) / np.sqrt(n)) if (n > 1 and d.std(ddof=1) > 0) else 0.0
        SE[i, j] = SE[j, i] = se
        DD[(i, j)] = d

    # ---- best competitor, gap ------------------------------------------------
    jstar = np.zeros(k, dtype=int)
    for i in range(k):
        others = [j for j in range(k) if j != i]
        jstar[i] = max(others, key=lambda j: mu[j])
    gap = np.array([D[i, jstar[i]] for i in range(k)])

    # ---- bootstrap max-|t| over all pairs (simultaneous calibration) --------
    # TRUE bootstrap-t: every replicate recomputes its own per-pair standard
    # error S* from the resampled differences, so T* = (D* - D)/S* has the
    # same studentized form as the observed statistics and the (1-alpha)
    # quantile of max|T*| calibrates the family of all pairwise t's INCLUDING
    # the denominator's sampling variability and the cross-pair selection (a
    # plug-in fixed-SE bootstrap would under-estimate the critical value and
    # reject too often - the calibration is regression-tested in
    # test_stat_tests.py Part F.1b).
    # D*_ij is the difference of the replicate means (mean of paired
    # differences = difference of means).
    B = int(min(max(int(n_resamples), 200), MAX_BOOTSTRAP))
    # the per-method r CI is DESCRIPTIVE (not the simultaneous quantity), so
    # its resamples are capped tighter than the calibration bootstrap to keep
    # the runtime bounded on large families (k(k-1)/2 pairs x k methods).
    r_B = int(min(B, 1000))
    rng = np.random.default_rng(seed)
    max_abs_t = np.empty(B)
    r_boot = {i: np.empty(r_B) for i in range(k)}   # per-method r resamples
    pair_idx = np.array(pairs)                    # (n_pairs, 2)
    D_pairs = D[pair_idx[:, 0], pair_idx[:, 1]] if len(pairs) else np.zeros(0)
    for b in range(B):
        idx = rng.integers(0, n, n)
        vb = vals[idx]
        if len(pairs):
            dij = vb[:, pair_idx[:, 0]] - vb[:, pair_idx[:, 1]]   # n x n_pairs
            D_b = dij.mean(axis=0)
            SE_b = dij.std(axis=0, ddof=1) / np.sqrt(n) if n > 1 else np.zeros(len(pairs))
            ok = SE_b > 0
            t = np.zeros(len(pairs))
            t[ok] = (D_b[ok] - D_pairs[ok]) / SE_b[ok]
            max_abs_t[b] = float(np.max(np.abs(t)))
        else:
            max_abs_t[b] = 0.0
        if b < r_B:
            for i in range(k):
                r_boot[i][b] = _ST.rank_biserial_matched(vb[:, i] - vb[:, jstar[i]])
    c_boot = float(np.quantile(max_abs_t, 1.0 - alpha))

    # ---- MCB intervals: Tukey-style projection of the pairwise bands --------
    L = np.array([min(D[i, j] - c_boot * SE[i, j] for j in range(k) if j != i)
                  for i in range(k)])
    U = np.array([min(D[i, j] + c_boot * SE[i, j] for j in range(k) if j != i)
                  for i in range(k)])

    # ---- MCB-adjusted one-sided P per method ---------------------------------
    # single-step max-|t| adjusted two-sided P of i vs its best competitor,
    # converted to the one-sided H1 "i is inferior to the best of the others".
    p_one = np.ones(k)
    for i in range(k):
        se_i = SE[i, jstar[i]]
        t_i = (gap[i] / se_i) if se_i > 0 else 0.0
        p_two = float(np.mean(max_abs_t >= abs(t_i)))
        p_one[i] = (p_two / 2.0) if gap[i] < 0 else 1.0 - p_two / 2.0

    # ---- parametric cross-check (classical Hsu/Tukey, pooled scale) ----------
    # The classical construction assumes one common standard deviation of the
    # paired differences; the pooled estimate removes the per-pair-SE selection
    # that would otherwise make the cross-check anticonservative (the per-pair
    # maximum t is stochastically larger than the studentized range).  With
    # c_par = q_{k, n-1}/sqrt(2) the k = 2 case reproduces the paired-t
    # critical value exactly.
    if pairs:
        s_pool = float(np.sqrt(np.mean([np.var(DD[(i, j)], ddof=1) if n > 1
                                        else 0.0 for i, j in pairs])))
    else:
        s_pool = 0.0
    try:
        c_par = float(sps.studentized_range.ppf(1.0 - alpha, k, df=n - 1) / np.sqrt(2.0))
    except Exception:                                   # pragma: no cover
        c_par = float('nan')
    if np.isfinite(c_par) and s_pool > 0:
        se_pool = s_pool / np.sqrt(n)
        Lp = np.array([min(D[i, j] - c_par * se_pool for j in range(k) if j != i)
                       for i in range(k)])
        Up = np.array([min(D[i, j] + c_par * se_pool for j in range(k) if j != i)
                       for i in range(k)])
    else:
        Lp = np.full(k, np.nan)
        Up = np.full(k, np.nan)

    # ---- Friedman omnibus on the same complete blocks ------------------------
    fr = _ST.friedman_test([vals[:, j] for j in range(k)])
    mean_ranks = fr['mean_ranks'] if fr else [float('nan')] * k

    # ---- verdicts -------------------------------------------------------------
    inferior = [i for i in range(k) if U[i] < 0]
    leading = [i for i in range(k) if U[i] >= 0]
    order = sorted(range(k), key=lambda i: (-mu[i], methods[i]))
    best = order[0]
    unique_winner = bool(L[best] > 0 and len(leading) == 1)

    rows = []
    for i in range(k):
        d_vs_best = vals[:, i] - vals[:, jstar[i]]
        r_i = _ST.rank_biserial_matched(d_vs_best)
        r_lo, r_hi = (float(np.quantile(r_boot[i], alpha / 2.0)),
                      float(np.quantile(r_boot[i], 1.0 - alpha / 2.0)))
        rows.append({
            'method': str(methods[i]),
            'n_units': int(n),
            'mean': float(mu[i]),
            'sd_units': float(sd[i]),
            'friedman_mean_rank': float(mean_ranks[i]),
            'rank_position': int(sum(1 for j in range(k) if mu[j] > mu[i]) + 1),
            'best_competitor': str(methods[jstar[i]]),
            'gap_to_best': float(gap[i]),
            'se_gap': float(SE[i, jstar[i]]),
            'mcb_low': float(L[i]),
            'mcb_high': float(U[i]),
            'mcb_low_param': float(Lp[i]),
            'mcb_high_param': float(Up[i]),
            'c_boot': float(c_boot),
            'c_par': float(c_par),
            's_pool': float(s_pool),
            'pvalue_mcb_one_sided': float(p_one[i]),
            'pvalue_mcb_two_sided_adj': float(2.0 * min(p_one[i], 1.0 - p_one[i])),
            'rank_biserial_r_vs_best': float(r_i),
            'ci95_r_low': r_lo,
            'ci95_r_high': r_hi,
            'cl_effect_vs_best': _ST.common_language_paired(d_vs_best),
            'significant_inferior': bool(U[i] < 0),
            'unique_best': bool(L[i] > 0),
            'in_leading_group': bool(U[i] >= 0),
            'note': '',
        })
    for r in rows:
        if r['se_gap'] == 0.0:
            r['note'] = 'zero within-unit variance of the gap; interval and P are degenerate'

    verdict = {
        'family': None,     # filled by the caller
        'n_units': int(n),
        'k_methods': int(k),
        'friedman_pvalue': float(fr['pvalue']) if fr else float('nan'),
        'friedman_chi2': float(fr['chi2']) if fr else float('nan'),
        'kendalls_w': float(fr['kendalls_w']) if fr else float('nan'),
        'ordering_exists': bool(fr is not None and fr['pvalue'] <= alpha),
        'order_by_mean': [str(methods[i]) for i in order],
        'sample_best': str(methods[best]),
        'unique_winner': unique_winner,
        'leading_group': [str(methods[i]) for i in leading],
        'leading_group_size': int(len(leading)),
        'inferior': [str(methods[i]) for i in inferior],
        'top2_exists': bool(len(leading) == 2),
    }
    return {'n_units': int(n), 'k_methods': int(k), 'units_dropped': int(units_dropped),
            'methods_dropped': all_missing, 'friedman': fr, 'c_boot': c_boot,
            'c_par': c_par, 'rows': rows, 'verdict': verdict, 'notes': notes}


# --------------------------------------------------------------------------- #
# Family builders: the SAME unit x method matrices the pairwise tests use      #
# --------------------------------------------------------------------------- #
def _resolve_cluster_cols(df, cluster_key_cols):
    """Resolve the cluster key exactly like stat_tests.run_caller_benchmark_stats.

    Returns (cluster_cols, cluster_key_str, cell2cluster_needed).  None ->
    DEFAULT_CLUSTER_KEY (donor); [] -> naive per-cell mode; a list -> filtered
    by availability (missing columns warn)."""
    if cluster_key_cols is None:
        cluster_cols = list(_ST.DEFAULT_CLUSTER_KEY)
        cluster_key_source = 'default (donor = independent human donor)'
    else:
        cluster_cols = list(cluster_key_cols)
        cluster_key_source = 'user-specified'
    dropped = [c for c in cluster_cols if c not in df.columns]
    cluster_cols = [c for c in cluster_cols if c in df.columns]
    if dropped:
        logging.warning('mcb: cluster-key columns %s not present in the input; '
                        'clustering uses the remaining %s',
                        dropped, cluster_cols or 'nothing')
    return cluster_cols, cluster_key_source


def _iter_caller_families(df, perf_metrics, gamete_type2short,
                          pair_key_cols=None, cluster_key_cols=None,
                          cluster_agg='median'):
    """Yield the Fig. 2 MCB families: (scenario tag, metric, donor x caller
    matrix, cluster-key string), one per (ground-truth scenario, metric).

    The matrix is built EXACTLY like the Friedman 'cluster' rows of
    stat_tests.run_caller_benchmark_stats: cell-level pivot first (median over
    duplicated (cell, caller) rows), then aggregation to the independent unit
    (donor by default) with the same cluster_agg, so MCB and the pairwise
    tests operate on identical numbers."""
    key_cols = [c for c in (pair_key_cols or _ST.DEFAULT_CELL_KEY) if c in df.columns]
    if key_cols:
        dup = df.duplicated(subset=key_cols + ['Caller']).sum()
        if dup:
            logging.warning('mcb: %d duplicate (cell, caller) rows: aggregating by median', dup)
    callers = sorted(df['Caller'].unique())
    cluster_cols, _src = _resolve_cluster_cols(df, cluster_key_cols)
    can_cluster = bool(cluster_cols) and bool(key_cols)
    cell2cluster = _ST._cell2cluster_map(df, key_cols, cluster_cols) if can_cluster else None
    cluster_key_str = '|'.join(cluster_cols) if can_cluster else ''
    if not can_cluster:
        logging.warning('mcb: cluster analysis impossible; the families fall back to the '
                        'NAIVE per-cell level (pseudoreplication risk, as in stat_tests)')
    for gamete_type, tag in gamete_type2short.items():
        for metric in perf_metrics:
            col = F'{gamete_type}.{metric}'
            if col not in df.columns:
                continue
            if key_cols:
                series = df.groupby(key_cols + ['Caller'])[col].median()
                piv = series.unstack('Caller')
            else:
                piv = df.pivot(index=None, columns='Caller', values=col)
                piv.index = np.arange(len(piv))
            piv = piv.reindex(columns=callers)
            if cell2cluster is not None:
                lab_series = _ST._pivot_cluster_labels(piv, cell2cluster)
                cmat = pd.DataFrame(piv.to_numpy(dtype=float),
                                    index=pd.Index(lab_series.to_numpy(), dtype=object),
                                    columns=piv.columns)
                cmat = cmat.groupby(level=0).agg(cluster_agg).sort_index()
            else:
                cmat = piv
            yield {'scenario': tag, 'metric': metric, 'matrix': cmat,
                   'cluster_key': cluster_key_str}


def _usable_ploidy_value(v):
    """A metadata value is usable when it is not None/NaN/empty-ish."""
    return _ST._norm_missing(v) != ''


def _ploidy_unit_from_label(label, split=True):
    """Donor id parsed from a germline row label ``donor . sampleType . avgSpotLen``
    (first middle-dot segment); the whole label in naive mode.  Mirrors the
    evaluation script's _donor_from_dataset_label."""
    s = str(label).strip()
    if not s or s.lower() in ('none', 'nan', 'nat', 'null'):
        return None
    if split:
        for sep in ('·', '•'):
            if sep in s:
                head = s.split(sep, 1)[0].strip()
                return head if head else s
    return s


def _ploidy_clean(tab):
    """Drop failed and placeholder rows from a ploidy long table (defensive)."""
    tab = tab.copy()
    if 'failed' in tab.columns:
        keep = ~tab['failed'].astype(str).str.strip().str.lower().isin(
            ('true', '1', 'yes', 't', 'y'))
        tab = tab[keep]
    if 'dataset' in tab.columns:
        ds = tab['dataset']
        keep = (ds.notna() & (ds.astype(str).str.strip() != '')
                & (ds.astype(str) != 'None'))
        tab = tab[keep]
    if 'pct_within' in tab.columns:
        tab['pct_within'] = pd.to_numeric(tab['pct_within'], errors='coerce')
    return tab


def _ploidy_pool(tab, cluster_key_cols=None, cluster_agg='median', naive=False):
    """Pool all panels into one row per (independent unit, method).

    Mirrors the evaluation script's _pool_units semantics for the standalone
    CLI: unit = the requested cluster column(s) (default 'donor') when usable,
    otherwise the donor parsed from the dataset label; non-failed rows are
    preferred; the value is the median over the unit's usable rows.  Returns
    (pooled DataFrame with columns unit/method/pct_within, n_fallback)."""
    tab = _ploidy_clean(tab)
    cluster_cols = [c for c in (cluster_key_cols or ['donor']) if c in tab.columns]
    units = []
    n_fallback = 0
    for _idx, row in tab.iterrows():
        u = None
        if (not naive) and cluster_cols:
            vals = [row[c] for c in cluster_cols]
            if all(_usable_ploidy_value(v) for v in vals):
                u = ' | '.join(str(v).strip() for v in vals)
        if u is None:
            u = _ploidy_unit_from_label(row['dataset'], split=not naive)
            if cluster_cols and not naive:
                n_fallback += 1
        units.append(u)
    tab['_unit'] = units
    tab = tab[tab['_unit'].notna()]
    method_col = 'method' if 'method' in tab.columns else 'tool'
    rows = []
    for (unit, meth), g in tab.groupby(['_unit', method_col], sort=False):
        vals = pd.to_numeric(g['pct_within'], errors='coerce').dropna()
        rows.append({'unit': unit, 'method': meth,
                     'pct_within': float(vals.median()) if len(vals) else float('nan')})
    return pd.DataFrame(rows, columns=['unit', 'method', 'pct_within']), n_fallback


def _ploidy_panel_families(tab, cluster_agg='median'):
    """Yield the per-panel (subfigure a-d) MCB families of Fig. 3.

    Within each plot group the independent unit is chosen by the same
    donor -> cellLine -> dataset fallback chain as
    stat_tests.run_ploidy_benchmark_stats (a real-tumor sample is its own
    unit; missing donor metadata collapses to a shared '(missing)' cluster
    only when the whole chain is unusable, which the chain prevents by its
    dataset column)."""
    tab = _ploidy_clean(tab)
    chain = list(_ST.PLOIDY_CLUSTER_FALLBACK_CHAIN)
    for plot, sub in tab.groupby('plot', sort=True):
        chosen, unit_vals = None, None
        for cand in chain:
            if cand not in sub.columns:
                continue
            vals = sub[cand].map(_ST._norm_missing)
            uniq = set(v for v in vals.unique() if v != '')
            if len(uniq) >= 2:
                chosen, unit_vals = cand, vals
                break
        if chosen is None:
            logging.warning('mcb: panel %s has no usable cluster column; using the '
                            'dataset label as the unit', plot)
            chosen, unit_vals = 'dataset', sub['dataset'].map(_ST._norm_missing)
        sub = sub.copy()
        sub['_unit'] = unit_vals.to_numpy()
        rows = []
        for (unit, meth), g in sub.groupby(['_unit', 'method'], sort=False):
            vals = pd.to_numeric(g['pct_within'], errors='coerce').dropna()
            rows.append({'unit': unit, 'method': meth,
                         'pct_within': float(vals.median()) if len(vals) else float('nan')})
        long = pd.DataFrame(rows, columns=['unit', 'method', 'pct_within'])
        mat = long.pivot(index='unit', columns='method', values='pct_within')
        yield {'panel': str(plot), 'matrix': mat, 'cluster_key': chosen}


def _matrix_from_pooled(pooled):
    """Pooled long table (columns unit/method/pct_within, or the evaluation
    script's pooled table with dataset+method) -> unit x method matrix."""
    unit_col = 'unit' if 'unit' in pooled.columns else 'dataset'
    method_col = 'method' if 'method' in pooled.columns else 'tool'
    val_col = 'pct_within' if 'pct_within' in pooled.columns else 'pct_within'
    return pooled.pivot(index=unit_col, columns=method_col, values=val_col)


# --------------------------------------------------------------------------- #
# Runners                                                                      #
# --------------------------------------------------------------------------- #
def _write_outputs(out_prefix, tsv_rows, settings, tex_lines=None):
    out_dir = os.path.dirname(os.path.abspath(out_prefix))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    tsv_path = out_prefix + '.tsv'
    pd.DataFrame(tsv_rows).to_csv(tsv_path, sep='\t', index=False, na_rep='NA')
    with open(out_prefix + '.json', 'w') as fh:
        json.dump(settings, fh, indent=2)
    if tex_lines is not None:
        with open(out_prefix + '.tex', 'w') as fh:
            fh.write('\n'.join(tex_lines) + '\n')
        logging.info('mcb: wrote %s, .json, .tex', tsv_path)
    else:
        logging.info('mcb: wrote %s, .json', tsv_path)


def run_caller_mcb(df, out_prefix, perf_metrics, gamete_type2short,
                   pair_key_cols=None, cluster_key_cols=None, cluster_agg='median',
                   n_resamples=DEFAULT_N_RESAMPLES, seed=1, alpha=0.05,
                   metric_display=None, caller_display=None, write_tex=True,
                   table_label='tab:scwgs-perf-mcb'):
    """Fig. 2 (CNV-caller benchmark) MCB: one family per (scenario, metric).

    Writes <out_prefix>.tsv / .json and, by default, the booktabs table
    <out_prefix>.tex (layout mirrors the pairwise table: rows =
    (metric, caller), two scenario groups of the four statistics n, p, r,
    95% CI - the CI being the simultaneous MCB interval of the gap to the
    best).  Returns the settings dict, or None when nothing was analysable.
    """
    if not _HAVE_SCIPY:
        logging.warning('scipy is not available: MCB skipped')
        return None
    rows_out, verdicts, skipped = [], [], []
    n_families = 0
    for fam in _iter_caller_families(df, perf_metrics, gamete_type2short,
                                     pair_key_cols=pair_key_cols,
                                     cluster_key_cols=cluster_key_cols,
                                     cluster_agg=cluster_agg):
        n_families += 1
        res = mcb_analyse(fam['matrix'], alpha=alpha, n_resamples=n_resamples,
                          seed=seed)
        if res is None:
            skipped.append(F"{fam['scenario']}/{fam['metric']}")
            continue
        res['verdict']['family'] = F"{fam['scenario']}/{fam['metric']}"
        verdicts.append(res['verdict'])
        for rec in res['rows']:
            rec = dict(rec)
            rec.update({
                'scenario': fam['scenario'], 'metric': fam['metric'],
                'task': 'Fig2 (scWGS CNV-caller benchmark)',
                'cluster_key': fam['cluster_key'],
                'friedman_pvalue': res['verdict']['friedman_pvalue'],
                'kendalls_w': res['verdict']['kendalls_w'],
                'family_unique_winner': res['verdict']['unique_winner'],
                'family_leading_group': '|'.join(res['verdict']['leading_group']),
                'family_leading_group_size': res['verdict']['leading_group_size'],
                'family_notes': '; '.join(res['notes']),
            })
            rows_out.append(rec)
    if not rows_out and not verdicts:
        logging.warning('mcb (caller): no analysable (scenario, metric) family%s',
                        F' (skipped: {", ".join(skipped)})' if skipped else '')
        return None

    settings = {
        'analysis': "Hsu's MCB (comparison with the best) - scWGS CNV-caller "
                    'benchmark (Fig. 2)',
        'design': 'per (ground-truth scenario, metric) family: donor x caller '
                  'complete-block matrix (per-donor medians of the per-cell '
                  'values, cluster key as in stat_tests); MCB intervals for '
                  'theta_i - max_{j != i} theta_j',
        'inference_level': 'cluster (independent experimental units; donor by default)',
        'interval_method': 'studentized cluster bootstrap max-|t| (Tukey-style '
                           'projection; single-step, family = all caller pairs '
                           'within the (scenario, metric) family)',
        'p_definition': 'MCB-adjusted one-sided P of H0: theta_i >= max_{j != i} '
                        'theta_j against H1: inferior to the best (single-step '
                        'bootstrap max-|t|, halved two-sided)',
        'effect_size': 'matched-pairs rank-biserial r versus the best competitor '
                       '(positive = the caller scores higher), 95% percentile '
                       'bootstrap CI (descriptive, per method); the table CI is '
                       'the simultaneous MCB interval of the gap to the best',
        'parametric_cross_check': 'studentized range q(k, n-1)/sqrt(2) (k = 2 '
                                  'reproduces the paired-t critical value)',
        'omnibus': 'Friedman test on the same complete-block matrix',
        'tail': 'one-sided p (inferiority to the best); intervals two-sided',
        'verdict_rules': 'interval entirely below 0 = significantly inferior; '
                         'entirely above 0 = unique winner; contains 0 = member '
                         'of the leading group (g = group size; g = 2 is the '
                         'top-2 situation)',
        'n_families': n_families,
        'n_families_skipped': len(skipped),
        'skipped_families': skipped,
        'bootstrap_n_resamples': int(min(max(int(n_resamples), 200), MAX_BOOTSTRAP)),
        'bootstrap_seed': int(seed),
        'alpha': float(alpha),
        'per_family_verdicts': verdicts,
        'software': _ST.software_versions(),
    }
    # consensus across the metric families, per scenario: the methods that are
    # never significantly inferior in any family of that scenario, and how
    # often each method is the sample-best.
    cons = {}
    for tag in sorted({v['family'].split('/')[0] for v in verdicts}):
        fams = [v for v in verdicts if v['family'].split('/')[0] == tag]
        never_inferior = sorted(set.intersection(*[
            set(v['leading_group']) for v in fams])) if fams else []
        best_counts = {}
        for v in fams:
            best_counts[v['sample_best']] = best_counts.get(v['sample_best'], 0) + 1
        cons[tag] = {
            'n_metric_families': len(fams),
            'methods_never_significantly_inferior': never_inferior,
            'sample_best_counts': best_counts,
            'families_with_unique_winner': sum(1 for v in fams if v['unique_winner']),
            'families_with_top2': sum(1 for v in fams if v['top2_exists']),
        }
    settings['consensus_per_scenario'] = cons

    tex_lines = None
    if write_tex:
        tex_lines = caller_mcb_latex_lines(
            pd.DataFrame(rows_out), metric_display=metric_display,
            caller_display=caller_display, alpha=alpha, table_label=table_label)
    _write_outputs(out_prefix, rows_out, settings, tex_lines)
    return settings


def run_ploidy_mcb(tab, out_prefix, cluster_key_cols=None, cluster_agg='median',
                   n_resamples=DEFAULT_N_RESAMPLES, seed=1, alpha=0.05,
                   already_pooled=False, long_tab=None, per_panel=True,
                   write_tex=True, table_label='tab:scwgs-ploidy-mcb',
                   method_display=None):
    """Fig. 3 (ploidy benchmark) MCB.

    Primary family: ALL panels POOLED (one row per independent unit x method,
    the same pooled table the pairwise tests use when called from the
    evaluation script with already_pooled=True; the standalone CLI pools
    first).  Sensitivity families: one per balloon panel (subfigure a-d),
    taken from long_tab when given (the evaluation script passes its full
    long table alongside the pooled table), else from tab itself.
    Writes <out_prefix>.tsv (pooled rows; out_prefix conventionally ends with
    '.pooled.stats.mcb'), <root>_perpanel.tsv (per-panel rows) and the
    booktabs table <out_prefix>.tex with one section per family.  Returns the
    settings dict or None."""
    if not _HAVE_SCIPY:
        logging.warning('scipy is not available: MCB skipped')
        return None
    rows_pooled, rows_panels, verdicts = [], [], []
    if already_pooled:
        pooled_mat = _matrix_from_pooled(tab)
        pooled_desc = 'pooled by the calling evaluation script (_pool_units)'
    else:
        pooled, n_fallback = _ploidy_pool(tab, cluster_key_cols=cluster_key_cols,
                                          cluster_agg=cluster_agg)
        pooled_mat = _matrix_from_pooled(pooled)
        pooled_desc = (F'pooled by mcb.py (donor column / label fallback; '
                       F'{n_fallback} label-fallback row(s))')
    res = mcb_analyse(pooled_mat, alpha=alpha, n_resamples=n_resamples, seed=seed)
    pooled_verdict = None
    if res is not None:
        res['verdict']['family'] = 'POOLED (all panels)'
        pooled_verdict = res['verdict']
        verdicts.append(pooled_verdict)
        for rec in res['rows']:
            rec = dict(rec)
            rec.update({
                'family': 'POOLED (all panels)', 'panel': 'POOLED',
                'task': 'Fig3 (scWGS ploidy benchmark, all panels pooled)',
                'cluster_key': pooled_desc,
                'friedman_pvalue': res['verdict']['friedman_pvalue'],
                'kendalls_w': res['verdict']['kendalls_w'],
                'family_unique_winner': res['verdict']['unique_winner'],
                'family_leading_group': '|'.join(res['verdict']['leading_group']),
                'family_leading_group_size': res['verdict']['leading_group_size'],
                'family_notes': '; '.join(res['notes']),
            })
            rows_pooled.append(rec)

    panel_verdicts = []
    panel_source = long_tab if long_tab is not None else (None if already_pooled else tab)
    if per_panel and panel_source is not None:
        for fam in _ploidy_panel_families(panel_source, cluster_agg=cluster_agg):
            r = mcb_analyse(fam['matrix'], alpha=alpha, n_resamples=n_resamples,
                            seed=seed)
            if r is None:
                continue
            r['verdict']['family'] = F"panel {fam['panel']}"
            panel_verdicts.append(r['verdict'])
            for rec in r['rows']:
                rec = dict(rec)
                rec.update({
                    'family': F"panel {fam['panel']}", 'panel': fam['panel'],
                    'task': F'Fig3 (scWGS ploidy benchmark, subfigure {fam["panel"]})',
                    'cluster_key': fam['cluster_key'],
                    'friedman_pvalue': r['verdict']['friedman_pvalue'],
                    'kendalls_w': r['verdict']['kendalls_w'],
                    'family_unique_winner': r['verdict']['unique_winner'],
                    'family_leading_group': '|'.join(r['verdict']['leading_group']),
                    'family_leading_group_size': r['verdict']['leading_group_size'],
                    'family_notes': '; '.join(r['notes']),
                })
                rows_panels.append(rec)
    verdicts.extend(panel_verdicts)
    if not verdicts:
        logging.warning('mcb (ploidy): no analysable family')
        return None

    root = out_prefix.replace('.pooled', '', 1) if out_prefix.endswith('.pooled.stats.mcb') \
        else out_prefix
    if not root.endswith('.stats.mcb'):
        root = root + '.stats.mcb'
    settings = {
        'analysis': "Hsu's MCB (comparison with the best) - ploidy-estimation "
                    'benchmark (Fig. 3)',
        'design': 'pooled family (primary; every donor of all panels, one row per '
                  'unit x method) plus one sensitivity family per balloon panel',
        'unit': 'per-sample percentage of cells with inferred ploidy within the window',
        'interval_method': 'studentized cluster bootstrap max-|t| (Tukey-style '
                           'projection; single-step, family = all method pairs '
                           'within the MCB family)',
        'p_definition': 'MCB-adjusted one-sided P of H0: theta_i >= max_{j != i} '
                        'theta_j against H1: inferior to the best',
        'effect_size': 'matched-pairs rank-biserial r versus the best competitor; '
                       'the table CI is the simultaneous MCB interval of the gap '
                       'to the best (percentage points)',
        'parametric_cross_check': 'studentized range q(k, n-1)/sqrt(2)',
        'omnibus': 'Friedman test on the same complete-block matrix',
        'verdict_rules': 'interval entirely below 0 = significantly inferior; '
                         'entirely above 0 = unique winner; contains 0 = member '
                         'of the leading group',
        'pooled_family': pooled_desc,
        'bootstrap_n_resamples': int(min(max(int(n_resamples), 200), MAX_BOOTSTRAP)),
        'bootstrap_seed': int(seed),
        'alpha': float(alpha),
        'per_family_verdicts': verdicts,
        'consensus': {
            'pooled': (pooled_verdict or {}),
            'per_panel': [v for v in panel_verdicts],
        },
        'software': _ST.software_versions(),
    }
    tex_lines = None
    if write_tex:
        tex_lines = ploidy_mcb_latex_lines(rows_pooled, rows_panels,
                                           alpha=alpha, table_label=table_label,
                                           method_display=method_display)
    _write_outputs(out_prefix, rows_pooled, settings, tex_lines)
    if rows_panels:
        pd.DataFrame(rows_panels).to_csv(
            root + '_perpanel.tsv', sep='\t', index=False, na_rep='NA')
        logging.info('mcb: wrote %s_perpanel.tsv', root)
    return settings


# --------------------------------------------------------------------------- #
# LaTeX tables (exactly n, p, r, 95% CI after the identity columns)           #
# --------------------------------------------------------------------------- #
_MCB_CAPTION_CORE = (
    'Hsu\'s multiple comparison with the best (MCB): every row is one caller '
    'compared with the best of the others.  After the identity columns each '
    'row carries exactly $n$ (independent units), $p$ (MCB-adjusted one-sided '
    'P; H$_0$: the method is at least as good as the best, small $P$ = '
    'significantly inferior), $r$ (matched-pairs rank-biserial effect size '
    'versus the best competitor; positive = better) and the simultaneous 95\\% '
    'MCB interval of the performance gap to the best '
    '($\\theta_i - \\max_{j \\ne i} \\theta_j$; entirely below 0 = '
    'significantly inferior, entirely above 0 = unique best, brackets 0 = '
    'indistinguishable from the best, i.e.\\ member of the leading group).')


def _disp(mapping, key):
    """Display label from an optional {id: label} mapping, else the id."""
    if mapping and key in mapping:
        return str(mapping[key])
    return str(key)


def caller_mcb_latex_lines(rows_df, metric_display=None, caller_display=None,
                           alpha=0.05, table_label='tab:scwgs-perf-mcb'):
    """Booktabs table for the Fig. 2 caller MCB, mirroring the pairwise-table
    layout: rows = (metric, caller); each ground-truth scenario owns one group
    of the four statistics (n, p, r, 95% CI).  Returns the line list, or None
    on any problem (a reason is logged)."""
    if rows_df is None or rows_df.empty:
        logging.error('mcb latex: no rows')
        return None
    sub = rows_df.copy()
    sub['scenario'] = sub['scenario'].astype(str)
    metric_order = list(dict.fromkeys(sub['metric'].tolist()))
    caller_order = list(dict.fromkeys(sub.sort_values(['metric', 'rank_position'])
                                      ['method'].tolist()))
    stat_keys = [('n', 'n_units'), ('p', 'pvalue_mcb_one_sided'),
                 ('r', 'rank_biserial_r_vs_best'),
                 ('lo', 'mcb_low'), ('hi', 'mcb_high')]
    piv = None
    for tag, key in stat_keys:
        p = sub.pivot_table(index=['metric', 'method'], columns='scenario',
                            values=key, aggfunc='first')
        p.columns = [F'{sc}_{tag}' for sc in p.columns]
        piv = p if piv is None else piv.join(p)
    scenarios = sorted(set(sub['scenario']))
    for sc in scenarios:
        for tag, _key in stat_keys:
            col = F'{sc}_{tag}'
            if col not in piv.columns:
                piv[col] = float('nan')
    piv = piv[[F'{sc}_{tag}' for sc in scenarios for tag, _key in stat_keys]]
    piv = piv.reset_index()
    piv['metric'] = pd.Categorical(piv['metric'], categories=metric_order, ordered=True)
    piv['method'] = pd.Categorical(piv['method'], categories=caller_order, ordered=True)
    piv = piv.sort_values(['metric', 'method'])
    n_groups = len(scenarios)
    stat_header = '$n$ & $p$ & $r$ & 95\\% CI'
    mc_cells = ' & '.join(F'\\multicolumn{{4}}{{c}}{{{_tex_escape(sc)}}}' for sc in scenarios)
    n_id = 2          # metric & caller
    first_group = n_id + 1
    cmid = '    ' + ''.join(F'\\cmidrule(lr){{{first_group + 4 * g}-{first_group + 4 * g + 3}}}'
                            for g in range(n_groups))
    header = [
        F'    Metric & Caller & {mc_cells} \\\\',
        cmid,
        F'     & {" & ".join([stat_header] * n_groups)} \\\\',
    ]
    colspec = 'll' + 'rrrr' * n_groups
    caption = (F'Hsu\'s MCB (comparison with the best) for the scWGS CNV-caller '
               F'benchmark (Fig.~2), per ground-truth scenario and performance '
               F'metric: each row compares one caller with the best of the others '
               F'at the donor level (per-donor medians; complete blocks).  '
               + _MCB_CAPTION_CORE)
    lines = [
        F'% LaTeX table generated by mcb.py (requires \\usepackage{{booktabs}})',
        '\\begin{table}[htbp]',
        '  \\centering',
        '  \\small',
        F'  \\caption{{{caption}}}',
        F'  \\label{{{table_label}}}',
        F'  \\begin{{tabular}}{{{colspec}}}',
        '    \\toprule',
        *header,
        '    \\midrule',
    ]
    current_metric = None
    for rec in piv.itertuples(index=False):
        sc_cells = []
        for sc in scenarios:
            n_v = _fmt_signed_int(getattr(rec, F'{sc}_n', None))
            p_v = _fmt_pvalue_mcb(getattr(rec, F'{sc}_p', None), alpha=alpha)
            r_v = _fmt_effect(getattr(rec, F'{sc}_r', None))
            ci_v = _fmt_ci(getattr(rec, F'{sc}_lo', None), getattr(rec, F'{sc}_hi', None))
            sc_cells.append(F'{n_v} & {p_v} & {r_v} & {ci_v}')
        metric_cell = _tex_escape(_disp(metric_display, str(rec.metric))) \
            if str(rec.metric) != current_metric else ''
        current_metric = str(rec.metric)
        lines.append(F'    {metric_cell} & {_tex_escape(_disp(caller_display, str(rec.method)))} '
                     F'& {" & ".join(sc_cells)} \\\\')
    lines += ['    \\bottomrule', '  \\end{tabular}', '\\end{table}']
    return lines


def ploidy_mcb_latex_lines(rows_pooled, rows_panels, alpha=0.05,
                           table_label='tab:scwgs-ploidy-mcb', method_display=None):
    """Booktabs table for the Fig. 3 ploidy MCB: one section per family (the
    pooled family first, then one section per balloon panel), rows = methods,
    columns = exactly n, p, r, 95% CI.  Returns the line list or None."""
    sections = []
    if rows_pooled:
        sections.append(('All panels pooled (each donor one independent sample)',
                         rows_pooled))
    by_panel = {}
    for rec in rows_panels:
        by_panel.setdefault(rec['panel'], []).append(rec)
    for panel in sorted(by_panel):
        sections.append((F'Panel {panel}', by_panel[panel]))
    if not sections:
        logging.error('mcb latex (ploidy): no rows')
        return None
    caption = (F'Hsu\'s MCB (comparison with the best) for the ploidy-estimation '
               F'benchmark (Fig.~3): each row compares one method with the best of '
               F'the others on the per-sample percentage of cells within the '
               F'window.  The first section pools every donor of the four panels '
               F'(primary); the remaining sections repeat the analysis per '
               F'balloon panel (sensitivity).  ' + _MCB_CAPTION_CORE)
    lines = [
        F'% LaTeX table generated by mcb.py (requires \\usepackage{{booktabs}})',
        '\\begin{table}[htbp]',
        '  \\centering',
        '  \\small',
        F'  \\caption{{{caption}}}',
        F'  \\label{{{table_label}}}',
        '  \\begin{tabular}{lrrrr}',
        '    \\toprule',
        '    Method & $n$ & $p$ & $r$ & 95\\% CI \\\\',
        '    \\midrule',
    ]
    for title, rows in sections:
        lines.append(F'    \\multicolumn{{5}}{{l}}{{\\textit{{{_tex_escape(title)}}}}} \\\\')
        rows = sorted(rows, key=lambda r: (r.get('rank_position', 99), r['method']))
        for rec in rows:
            n_v = _fmt_signed_int(rec.get('n_units'))
            p_v = _fmt_pvalue_mcb(rec.get('pvalue_mcb_one_sided'), alpha=alpha)
            r_v = _fmt_effect(rec.get('rank_biserial_r_vs_best'))
            ci_v = _fmt_ci(rec.get('mcb_low'), rec.get('mcb_high'))
            lines.append(F'    {_tex_escape(_disp(method_display, rec["method"]))} '
                         F'& {n_v} & {p_v} & {r_v} & {ci_v} \\\\')
        lines.append('    \\midrule')
    if lines[-1] == '    \\midrule':
        lines[-1] = '    \\bottomrule'
    lines += ['  \\end{tabular}', '\\end{table}']
    return lines


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
    if raw is None:
        return None
    s = raw.strip().lower()
    if s in ('none', 'naive', 'off', 'no'):
        return []
    return [c.strip() for c in raw.split(',') if c.strip()]


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Hsu's MCB (comparison with the best) for the scWGS benchmark "
                    'figures: Fig. 2 (caller benchmark, from the long TSV on stdin '
                    'or -i) and Fig. 3 (ploidy benchmark, from the balloon long '
                    'table via -i).')
    ap.add_argument('-t', '--task', choices=['caller', 'ploidy'], default='caller',
                    help='Which figure task to analyse: caller = Fig. 2 (long TSV), '
                         'ploidy = Fig. 3 (pct_within long table).')
    ap.add_argument('-i', '--input', nargs='*', default=[],
                    help='Input file(s)/glob(s). Caller mode also reads the long TSV '
                         'from stdin when -i is empty.')
    ap.add_argument('-o', '--output', required=True,
                    help='Output prefix. Caller mode writes <prefix>.tsv/.json/.tex; '
                         'ploidy mode writes <prefix>.pooled... files (see the '
                         'module docstring).')
    ap.add_argument('--metrics', default=None, metavar='IDS',
                    help='Comma-separated metric ids to analyse (caller mode; '
                         'default: every metric column found).')
    ap.add_argument('--scenarios', default=None, metavar='PREFIXES',
                    help='Comma-separated gamete-type scenario prefixes (caller '
                         'mode; default: with_haploidy_assumed_gametes,'
                         'with_aneuploidy_aware_gametes).')
    ap.add_argument('--cluster-key', default=None, metavar='COLS',
                    help='Independent-unit key: "donor" columns (caller/pooled '
                         'ploidy) or "none" for the naive per-unit level. Default: '
                         'donor (caller), the donor->cellLine->dataset chain '
                         '(per-panel ploidy).')
    ap.add_argument('--boot', type=int, default=DEFAULT_N_RESAMPLES, metavar='N',
                    help=F'Bootstrap resamples (default: {DEFAULT_N_RESAMPLES}, '
                         F'capped at {MAX_BOOTSTRAP}).')
    ap.add_argument('--seed', type=int, default=1, metavar='SEED',
                    help='Seed of the bootstrap RNG (default: 1).')
    ap.add_argument('--alpha', type=float, default=0.05, metavar='ALPHA',
                    help='Family-wise alpha (default: 0.05).')
    ap.add_argument('--no-tex', dest='tex', action='store_false', default=True,
                    help='Do not write the booktabs LaTeX table.')
    ap.add_argument('--latex-table', action='store_true', default=False,
                    help='(ploidy mode) rebuild and print the LaTeX table from an '
                         'existing <prefix>.pooled.stats.mcb.tsv and exit.')
    args = ap.parse_args(argv)

    cluster_key = _parse_cluster_key(args.cluster_key)
    if args.task == 'caller':
        if args.input:
            frames = [pd.read_csv(p, sep='\t') for p in _expand_inputs(args.input)]
            df = pd.concat(frames, ignore_index=True)
        else:
            df = pd.read_csv(sys.stdin, sep='\t')
        if args.metrics:
            perf_metrics = [m.strip() for m in args.metrics.split(',') if m.strip()]
        else:
            prefixes = ([p.strip() for p in args.scenarios.split(',')]
                        if args.scenarios else
                        ['with_haploidy_assumed_gametes', 'with_aneuploidy_aware_gametes'])
            found = set()
            for p in prefixes:
                for c in df.columns:
                    if c.startswith(p + '.'):
                        found.add(c.split('.', 1)[1])
            perf_metrics = list(dict.fromkeys(found))
        gamete2short = {'with_haploidy_assumed_gametes': 'Hap_0',
                        'with_aneuploidy_aware_gametes': 'Hap_1'}
        if args.scenarios:
            prefixes = [p.strip() for p in args.scenarios.split(',')]
            gamete2short = {p: gamete2short.get(p, p) for p in prefixes}
        settings = run_caller_mcb(
            df, args.output, perf_metrics, gamete2short,
            cluster_key_cols=cluster_key, n_resamples=args.boot, seed=args.seed,
            alpha=args.alpha, write_tex=args.tex)
        return 0 if settings is not None else 1

    # ploidy mode
    if args.latex_table:
        tsv = args.output if args.output.endswith('.pooled.stats.mcb') \
            else args.output + '.pooled.stats.mcb'
        tsv = F'{tsv}.tsv' if not tsv.endswith('.tsv') else tsv
        rows = pd.read_csv(tsv, sep='\t').to_dict('records')
        per_tsv = tsv.replace('.pooled.stats.mcb.tsv', '.stats.mcb_perpanel.tsv')
        rows_panels = pd.read_csv(per_tsv, sep='\t').to_dict('records') \
            if os.path.isfile(per_tsv) else []
        lines = ploidy_mcb_latex_lines(rows, rows_panels, alpha=args.alpha)
        if lines is None:
            return 1
        print('\n'.join(lines))
        return 0
    paths = _expand_inputs(args.input)
    if not paths:
        logging.error('ploidy mode needs -i <pct_within_long.tsv> (or a glob)')
        return 1
    frames = [pd.read_csv(p, sep='\t') for p in paths]
    tab = pd.concat(frames, ignore_index=True)
    settings = run_ploidy_mcb(
        tab, args.output + '.pooled.stats.mcb', cluster_key_cols=cluster_key,
        n_resamples=args.boot, seed=args.seed, alpha=args.alpha,
        per_panel=True, write_tex=args.tex)
    return 0 if settings is not None else 1


if __name__ == '__main__':
    sys.exit(main())
