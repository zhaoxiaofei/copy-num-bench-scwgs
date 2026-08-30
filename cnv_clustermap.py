#!/usr/bin/env python3

# https://sorryios.ai/c/6a7be034-0f44-83ea-a3e4-94d128612c3f

"""
Fast CNV clustermap from per-sample BED files.

Main design choices
-------------------
1. Rows (cells/samples) are clustered from CNV profiles with a length-weighted
   Manhattan (L1) distance rather than Euclidean distance.
2. Missing CN values are NOT treated as diploid for clustering. Pairwise
   distances use only bins observed in both cells. Missing values are filled
   with --center only for heatmap display.
3. Sex chromosomes are excluded from clustering by default, but remain visible
   in the heatmap. Use --include-sex-chromosomes-in-clustering to include them.
4. Clustering can automatically coarsen very high-resolution bin matrices to at
   most --cluster-max-bins features. The full-resolution matrix is still plotted.
5. Complete matrices use scipy.spatial.distance.pdist('cityblock'), which is
   implemented in compiled code. Matrices with NaNs use a Numba-parallel exact
   implementation when Numba is available, otherwise a NumPy block fallback.
6. Fixed-bin BED conversion is implemented as a segment-to-bin sweep, avoiding
   the expensive "all segments x every bin" overlap calculation.

Exact hierarchical clustering still requires O(n_cells^2) pairwise distances.
The --cluster-max-bins option reduces the expensive dependence on the number of
bins, but it cannot remove the quadratic dependence on the number of cells.
"""

import argparse
import glob
import os
import re
import sys
sys.setrecursionlimit(20000)
import time
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap

import seaborn as sns
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist, squareform

try:
    from numba import njit, prange
    HAVE_NUMBA = True
except Exception:
    HAVE_NUMBA = False


CHROM_ORDER = [f"chr{i}" for i in list(range(1, 23)) + ["X", "Y"]]
CHROM_SET = set(CHROM_ORDER)


def eprint(msg: str) -> None:
    sys.stderr.write(str(msg) + "\n")
    sys.stderr.flush()


def chrom_sort_key(chrom: str) -> int:
    c = str(chrom).replace("chr", "")
    if c == "X":
        return 23
    if c == "Y":
        return 24
    try:
        return int(c)
    except Exception:
        return 99


def load_bed(path: str) -> pd.DataFrame:
    """Read the first four BED columns: chrom, start, end, copy number."""
    try:
        df = pd.read_csv(
            path,
            sep="\t",
            comment="#",
            header=None,
            usecols=[0, 1, 2, 3],
            names=["chrom", "start", "end", "cn"],
            dtype={"chrom": str},
            engine="c",
            on_bad_lines="skip",
        )
    except Exception:
        # More permissive fallback for unusual BED-like text files.
        df = pd.read_csv(
            path,
            sep="\t",
            comment="#",
            header=None,
            usecols=[0, 1, 2, 3],
            names=["chrom", "start", "end", "cn"],
            dtype={"chrom": str},
            engine="python",
            on_bad_lines="skip",
        )

    for c in ["start", "end", "cn"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["chrom", "start", "end", "cn"]).copy()

    if df.empty:
        return pd.DataFrame(columns=["chrom", "start", "end", "cn"])

    df["chrom"] = df["chrom"].astype(str)
    df["start"] = np.rint(df["start"].to_numpy()).astype(np.int64)
    df["end"] = np.rint(df["end"].to_numpy()).astype(np.int64)
    df["cn"] = df["cn"].astype(np.float32)
    df = df[df["end"] > df["start"]]
    return df


def sample_name_from_path(path: str, regex: str = "") -> str:
    base = re.sub(r"\.bed$", "", os.path.basename(path))
    if regex:
        m = re.search(regex, base)
        if m:
            return m.group(1) if m.groups() else m.group(0)
    return base


def parse_chrom_sizes(fai_path: str) -> Dict[str, int]:
    sizes: Dict[str, int] = {}
    with open(fai_path) as fh:
        for line in fh:
            toks = line.rstrip("\n").split("\t")
            if len(toks) >= 2:
                try:
                    sizes[toks[0]] = int(toks[1])
                except ValueError:
                    pass
    return sizes


def make_fixed_bin_layout(
    chrom_sizes: Dict[str, int], bin_size: int
) -> Tuple[List[str], np.ndarray, Dict[str, Tuple[int, int, int]]]:
    """
    Return column labels, physical bin lengths, and per-chromosome layout.

    layout[chrom] = (global_column_offset, number_of_bins, chromosome_length)
    """
    labels: List[str] = []
    lengths: List[float] = []
    layout: Dict[str, Tuple[int, int, int]] = {}

    offset = 0
    for chrom in sorted(chrom_sizes, key=chrom_sort_key):
        if chrom not in CHROM_SET:
            continue
        chrom_len = int(chrom_sizes[chrom])
        n_bins = max(1, int(np.ceil(chrom_len / float(bin_size))))
        layout[chrom] = (offset, n_bins, chrom_len)

        for b in range(n_bins):
            bs = b * bin_size
            be = min((b + 1) * bin_size, chrom_len)
            labels.append(f"{chrom}:{bs}-{be}")
            lengths.append(float(max(be - bs, 1)))

        offset += n_bins

    return labels, np.asarray(lengths, dtype=np.float64), layout


def bed_to_fixed_bins_fast(
    df: pd.DataFrame,
    bin_size: int,
    layout: Dict[str, Tuple[int, int, int]],
    n_columns: int,
) -> np.ndarray:
    """
    Convert segmented CN BED data to fixed bins in O(number of segments + bins).

    For each chromosome, segment overlap is accumulated into the first/last
    partial bins directly and into fully covered interior bins with difference
    arrays. This is substantially faster than testing every segment against
    every genomic bin.
    """
    out = np.full(n_columns, np.nan, dtype=np.float32)

    grouped = {chrom: sub for chrom, sub in df.groupby("chrom", sort=False)}

    for chrom, (offset, n_bins, chrom_len) in layout.items():
        sub = grouped.get(chrom)
        if sub is None or sub.empty:
            continue

        weighted = np.zeros(n_bins, dtype=np.float64)
        coverage = np.zeros(n_bins, dtype=np.float64)
        weighted_diff = np.zeros(n_bins + 1, dtype=np.float64)
        coverage_diff = np.zeros(n_bins + 1, dtype=np.float64)

        starts = sub["start"].to_numpy(dtype=np.int64, copy=False)
        ends = sub["end"].to_numpy(dtype=np.int64, copy=False)
        cns = sub["cn"].to_numpy(dtype=np.float64, copy=False)

        for start, end, cn in zip(starts, ends, cns):
            start = max(0, int(start))
            end = min(chrom_len, int(end))
            if end <= start or not np.isfinite(cn):
                continue

            b0 = min(start // bin_size, n_bins - 1)
            b1 = min((end - 1) // bin_size, n_bins - 1)

            if b0 == b1:
                ov = float(end - start)
                weighted[b0] += cn * ov
                coverage[b0] += ov
                continue

            first_end = min((b0 + 1) * bin_size, chrom_len)
            first_ov = float(max(first_end - start, 0))
            weighted[b0] += cn * first_ov
            coverage[b0] += first_ov

            last_start = b1 * bin_size
            last_ov = float(max(end - last_start, 0))
            weighted[b1] += cn * last_ov
            coverage[b1] += last_ov

            # Every bin strictly between b0 and b1 is fully covered.
            if b1 > b0 + 1:
                weighted_diff[b0 + 1] += cn * bin_size
                weighted_diff[b1] -= cn * bin_size
                coverage_diff[b0 + 1] += bin_size
                coverage_diff[b1] -= bin_size

        weighted += np.cumsum(weighted_diff[:-1])
        coverage += np.cumsum(coverage_diff[:-1])

        vals = np.full(n_bins, np.nan, dtype=np.float32)
        valid = coverage > 0
        vals[valid] = (weighted[valid] / coverage[valid]).astype(np.float32)
        out[offset : offset + n_bins] = vals

    return out


def bed_to_chrom_means_vector(
    df: pd.DataFrame, chroms: Sequence[str] = CHROM_ORDER
) -> np.ndarray:
    out = np.full(len(chroms), np.nan, dtype=np.float32)
    chrom_to_index = {c: i for i, c in enumerate(chroms)}

    for chrom, sub in df.groupby("chrom", sort=False):
        idx = chrom_to_index.get(str(chrom))
        if idx is None or sub.empty:
            continue
        w = (sub["end"] - sub["start"]).clip(lower=1).to_numpy(dtype=np.float64)
        cn = sub["cn"].to_numpy(dtype=np.float64)
        denom = w.sum()
        if denom > 0:
            out[idx] = float(np.dot(cn, w) / denom)
    return out


def clean_sample_names(names: Sequence[str]) -> List[str]:
    """Preserve the original script's common-prefix/suffix cleanup, safely."""
    names = list(names)
    if len(names) <= 1:
        return names

    common_prefix = os.path.commonprefix(names)
    common_preadd = common_prefix.split("_")[-1]
    common_suffix = os.path.commonprefix([x[::-1] for x in names])[::-1]

    cleaned: List[str] = []
    for x in names:
        y = x.removeprefix(common_prefix)
        if common_suffix and y.endswith(common_suffix):
            y = y[: -len(common_suffix)]
        cleaned.append(common_preadd + y)

    # Avoid accidental overwriting if cleanup creates duplicate names.
    seen: Dict[str, int] = {}
    unique: List[str] = []
    for x in cleaned:
        count = seen.get(x, 0) + 1
        seen[x] = count
        unique.append(x if count == 1 else f"{x}_dup{count}")
    return unique


def chromosome_of_column(col: str) -> str:
    return str(col).split(":", 1)[0]


def coarsen_cluster_matrix(
    X: np.ndarray,
    columns: Sequence[str],
    weights: np.ndarray,
    max_bins: int,
) -> Tuple[np.ndarray, List[str], np.ndarray]:
    """
    Chromosome-aware aggregation of adjacent bins for faster clustering.

    The weighted mean is computed independently for every cell, ignoring NaNs.
    Groups never cross chromosome boundaries. Group weights equal the genomic
    lengths represented by the original bins.
    """
    columns = list(columns)
    weights = np.asarray(weights, dtype=np.float64)
    n_cells, n_features = X.shape

    if max_bins <= 0 or n_features <= max_bins:
        return np.ascontiguousarray(X, dtype=np.float32), columns, weights

    chroms = [chromosome_of_column(c) for c in columns]

    # Columns are genomic-order contiguous by chromosome in this script.
    chrom_ranges: List[Tuple[str, int, int]] = []
    start = 0
    while start < n_features:
        chrom = chroms[start]
        end = start + 1
        while end < n_features and chroms[end] == chrom:
            end += 1
        chrom_ranges.append((chrom, start, end))
        start = end

    # Choose one global group width, then increase it until the total number of
    # chromosome-respecting groups is <= max_bins.
    group_width = max(1, int(np.ceil(n_features / float(max_bins))))

    def number_of_groups(width: int) -> int:
        return sum(int(np.ceil((end - start) / float(width)))
                   for _, start, end in chrom_ranges)

    while number_of_groups(group_width) > max_bins:
        group_width += 1

    coarse_parts: List[np.ndarray] = []
    coarse_cols: List[str] = []
    coarse_weights: List[np.ndarray] = []

    for chrom, start, end in chrom_ranges:
        sub = np.ascontiguousarray(X[:, start:end], dtype=np.float32)
        w = weights[start:end]
        n = end - start
        starts = np.arange(0, n, group_width, dtype=np.int64)

        valid = np.isfinite(sub)
        sub0 = np.nan_to_num(sub, nan=0.0, copy=True)

        # np.add.reduceat performs all groups for a chromosome at once.
        numerator = np.add.reduceat(sub0 * w[None, :], starts, axis=1)
        denominator = np.add.reduceat(valid * w[None, :], starts, axis=1)

        coarse = np.full(numerator.shape, np.nan, dtype=np.float32)
        np.divide(
            numerator,
            denominator,
            out=coarse,
            where=denominator > 0,
            casting="unsafe",
        )

        group_w = np.add.reduceat(w, starts)
        coarse_parts.append(coarse)
        coarse_weights.append(group_w)

        for k, s in enumerate(starts):
            e = min(int(s) + group_width, n)
            first_col = columns[start + int(s)]
            last_col = columns[start + e - 1]
            coarse_cols.append(f"{chrom}|{first_col}..{last_col}")

    Xc = np.ascontiguousarray(np.concatenate(coarse_parts, axis=1), dtype=np.float32)
    wc = np.concatenate(coarse_weights).astype(np.float64, copy=False)
    return Xc, coarse_cols, wc


def prepare_cluster_matrix(
    mat: pd.DataFrame,
    column_lengths: pd.Series,
    include_sex_chromosomes: bool,
    max_bins: int,
    min_feature_coverage: float,
    drop_invariant: bool,
) -> Tuple[np.ndarray, List[str], np.ndarray]:
    columns = list(mat.columns)

    if not include_sex_chromosomes:
        columns = [c for c in columns if chromosome_of_column(c) not in {"chrX", "chrY"}]

    if not columns:
        raise ValueError("No columns remain for clustering.")

    X = np.ascontiguousarray(mat[columns].to_numpy(dtype=np.float32), dtype=np.float32)
    weights = column_lengths.loc[columns].to_numpy(dtype=np.float64)

    X, columns, weights = coarsen_cluster_matrix(X, columns, weights, max_bins=max_bins)

    coverage = np.mean(np.isfinite(X), axis=0)
    keep = coverage >= min_feature_coverage
    if not np.any(keep):
        raise ValueError(
            "No clustering bins pass --cluster-min-feature-coverage. "
            "Lower that threshold or inspect missing CN calls."
        )

    X = X[:, keep]
    weights = weights[keep]
    columns = [c for c, k in zip(columns, keep) if k]

    if drop_invariant and X.shape[1] > 1:
        with np.errstate(all="ignore"):
            col_min = np.nanmin(X, axis=0)
            col_max = np.nanmax(X, axis=0)
        variable = np.isfinite(col_min) & np.isfinite(col_max) & ((col_max - col_min) > 1e-8)
        if np.any(variable):
            X = X[:, variable]
            weights = weights[variable]
            columns = [c for c, k in zip(columns, variable) if k]
        else:
            eprint("Warning: all clustering bins are invariant; retaining them so zero distances can be computed.")

    if X.shape[1] == 0:
        raise ValueError("No clustering bins remain after filtering.")

    if not np.all(np.isfinite(weights)) or np.any(weights <= 0):
        raise ValueError("Clustering bin weights must be finite and positive.")

    return np.ascontiguousarray(X, dtype=np.float32), columns, weights


if HAVE_NUMBA:
    @njit(parallel=True, cache=True)
    def _nan_weighted_cityblock_numba(
        X: np.ndarray,
        weights: np.ndarray,
        min_shared_weight: float,
    ) -> np.ndarray:
        n_cells, n_features = X.shape
        out = np.empty(n_cells * (n_cells - 1) // 2, dtype=np.float64)

        for i in prange(n_cells - 1):
            # First condensed-distance index for pair (i, i+1).
            base = n_cells * i - (i * (i + 1)) // 2
            for j in range(i + 1, n_cells):
                numerator = 0.0
                denominator = 0.0
                for k in range(n_features):
                    a = X[i, k]
                    b = X[j, k]
                    if np.isfinite(a) and np.isfinite(b):
                        w = weights[k]
                        d = a - b
                        if d < 0:
                            d = -d
                        numerator += w * d
                        denominator += w

                idx = base + (j - i - 1)
                if denominator >= min_shared_weight and denominator > 0.0:
                    out[idx] = numerator / denominator
                else:
                    out[idx] = np.nan
        return out


def _nan_weighted_cityblock_numpy(
    X: np.ndarray,
    weights: np.ndarray,
    min_shared_weight: float,
    block_size: int,
) -> np.ndarray:
    """Exact NaN-aware fallback using vectorized j-blocks."""
    n_cells, _ = X.shape
    out = np.empty(n_cells * (n_cells - 1) // 2, dtype=np.float64)
    weights = weights.astype(np.float64, copy=False)

    for i in range(n_cells - 1):
        xi = X[i].astype(np.float64, copy=False)
        finite_i = np.isfinite(xi)
        base = n_cells * i - (i * (i + 1)) // 2

        for j0 in range(i + 1, n_cells, block_size):
            j1 = min(j0 + block_size, n_cells)
            Y = X[j0:j1].astype(np.float64, copy=False)
            valid = np.isfinite(Y) & finite_i[None, :]

            denom = np.sum(valid * weights[None, :], axis=1)
            diff = np.abs(Y - xi[None, :])
            diff[~valid] = 0.0
            numer = np.sum(diff * weights[None, :], axis=1)

            d = np.full(j1 - j0, np.nan, dtype=np.float64)
            ok = (denom >= min_shared_weight) & (denom > 0)
            d[ok] = numer[ok] / denom[ok]

            out_start = base + (j0 - i - 1)
            out[out_start : out_start + (j1 - j0)] = d

    return out


def weighted_cityblock_condensed(
    X: np.ndarray,
    weights: np.ndarray,
    min_pair_overlap: float,
    backend: str,
    numpy_block_size: int,
) -> Tuple[np.ndarray, str]:
    """
    Length-weighted mean absolute CN difference for all cell pairs.

    d(i,j) = sum_k w_k * |CN_ik - CN_jk| / sum_k w_k,
    where sums include only bins observed in both cells.
    """
    X = np.ascontiguousarray(X, dtype=np.float32)
    weights = np.asarray(weights, dtype=np.float64)
    total_weight = float(weights.sum())
    if total_weight <= 0:
        raise ValueError("Total clustering weight is zero.")

    has_missing = bool(np.isnan(X).any())

    # Fastest exact path: compiled SciPy cityblock when the matrix is complete.
    if not has_missing and backend in {"auto", "scipy"}:
        normalized_w = weights / total_weight
        X_scaled = X.astype(np.float64, copy=False) * normalized_w[None, :]
        return pdist(X_scaled, metric="cityblock"), "scipy-pdist"

    if backend == "scipy" and has_missing:
        raise ValueError(
            "--distance-backend scipy cannot handle NaN-aware Manhattan distances. "
            "Use auto, numba, or numpy."
        )

    min_shared_weight = float(min_pair_overlap) * total_weight

    if backend in {"auto", "numba"} and HAVE_NUMBA:
        return (
            _nan_weighted_cityblock_numba(
                X,
                weights.astype(np.float64),
                min_shared_weight,
            ),
            "numba-parallel",
        )

    if backend == "numba" and not HAVE_NUMBA:
        raise RuntimeError(
            "--distance-backend numba was requested, but numba is not installed. "
            "Install numba or use --distance-backend numpy."
        )

    if has_missing and backend == "auto" and not HAVE_NUMBA:
        eprint(
            "Warning: clustering matrix contains NaNs and numba is unavailable; "
            "using the slower NumPy exact fallback. Installing numba will usually speed this up."
        )

    return (
        _nan_weighted_cityblock_numpy(
            X,
            weights,
            min_shared_weight,
            block_size=max(1, int(numpy_block_size)),
        ),
        "numpy-block",
    )


def choose_optimal_ordering(mode: str, n_cells: int, max_cells: int) -> bool:
    if mode == "on":
        return True
    if mode == "off":
        return False
    return n_cells <= max_cells


def estimate_condensed_gb(n_cells: int) -> float:
    n_pairs = n_cells * (n_cells - 1) // 2
    return n_pairs * 8.0 / (1024.0 ** 3)


def save_clustering_outputs(
    output_prefix: str,
    Z: np.ndarray,
    sample_names: Sequence[str],
    distances: np.ndarray,
    save_distance_matrix: bool,
) -> None:
    np.savetxt(
        output_prefix + ".linkage.tsv",
        Z,
        delimiter="\t",
        header="left_cluster\tright_cluster\tdistance\tn_members",
        comments="",
        fmt=["%.0f", "%.0f", "%.8g", "%.0f"],
    )

    order = leaves_list(Z).astype(int)
    order_df = pd.DataFrame(
        {
            "rank": np.arange(1, len(order) + 1, dtype=int),
            "sample": [sample_names[i] for i in order],
            "original_index": order,
        }
    )
    order_df.to_csv(output_prefix + ".cluster_order.tsv", sep="\t", index=False)

    if save_distance_matrix:
        D = squareform(distances)
        dist_df = pd.DataFrame(D, index=sample_names, columns=sample_names)
        dist_df.to_csv(output_prefix + ".cnv_distance.tsv.gz", sep="\t", compression="gzip")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate a fast, CNV-aware clustered heatmap from per-sample BED files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    p.add_argument("-i", "--input", nargs="+", required=True,
                   help="BED files or shell-style glob patterns")
    p.add_argument("-o", "--output-prefix", required=True)
    p.add_argument("--title", default="")

    p.add_argument(
        "--by-chrom",
        action="store_true",
        help="Use one length-weighted mean CN value per chromosome instead of fixed bins. "
             "If --bin-size is also supplied, --by-chrom takes precedence.",
    )
    p.add_argument(
        "--bin-size",
        type=int,
        default=0,
        help="Plotting bin size in bp. Values >0 use fixed bins and require --fai. "
             "With 0, chromosome means are used.",
    )
    p.add_argument("--fai", type=str, default="", help="Reference FASTA .fai file")
    p.add_argument("--sample-regex", type=str, default="")
    p.add_argument("--metadata-tsv", type=str, default="")

    p.add_argument(
        "--float-copy-numbers",
        action="store_true",
        help="Keep continuous CN values instead of rounding to integer copy numbers.",
    )
    p.add_argument("--cmap", type=str, default="RdBu_r")
    p.add_argument("--vmin", type=float, default=0)
    p.add_argument("--vmax", type=float, default=6)
    p.add_argument(
        "--center",
        type=float,
        default=2,
        help="Diploid display value used to fill NaNs in the plotted heatmap only.",
    )
    p.add_argument(
        "--show-sample-labels",
        type=int,
        choices=[0, 1],
        default=1,
        help="Show sample names on the heatmap y-axis.",
    )

    # Clustering controls.
    p.add_argument(
        "--cluster-max-bins",
        type=int,
        default=1000,
        help="Maximum number of genomic features used for clustering. If the plotting "
             "matrix has more bins, adjacent bins are chromosome-aware coarsened for "
             "clustering only. Use 0 for no coarsening.",
    )
    p.add_argument(
        "--cluster-min-feature-coverage",
        type=float,
        default=0.50,
        help="Minimum fraction of cells with an observed CN for a clustering feature.",
    )
    p.add_argument(
        "--cluster-min-pair-overlap",
        type=float,
        default=0.50,
        help="Minimum weighted fraction of clustering features jointly observed for a cell pair.",
    )
    p.add_argument(
        "--include-sex-chromosomes-in-clustering",
        action="store_true",
        help="Include chrX/chrY in clustering. They are excluded by default but still plotted.",
    )
    p.add_argument(
        "--keep-invariant-bins",
        action="store_true",
        help="Do not remove invariant clustering features. Dropping them is faster and does "
             "not change complete-data Manhattan clustering topology.",
    )
    p.add_argument(
        "--linkage-method",
        choices=["average", "complete", "weighted", "single"],
        default="average",
        help="Hierarchical linkage method used with the precomputed CNV distances.",
    )
    p.add_argument(
        "--optimal-ordering",
        choices=["auto", "on", "off"],
        default="auto",
        help="Optimal leaf ordering. 'auto' enables it only for smaller datasets because it can be slow.",
    )
    p.add_argument(
        "--optimal-ordering-max-cells",
        type=int,
        default=1200,
        help="Maximum number of cells for automatic optimal leaf ordering.",
    )
    p.add_argument(
        "--distance-backend",
        choices=["auto", "scipy", "numba", "numpy"],
        default="auto",
        help="Pairwise-distance backend. auto uses SciPy pdist for complete data and Numba "
             "for NaN-aware distances when available.",
    )
    p.add_argument(
        "--numpy-distance-block-size",
        type=int,
        default=128,
        help="Block size for the NumPy NaN-aware distance fallback.",
    )
    p.add_argument(
        "--max-distance-gb",
        type=float,
        default=8.0,
        help="Refuse exact hierarchical clustering if the condensed float64 distance vector "
             "alone would exceed this size. Increase only if sufficient RAM is available.",
    )
    p.add_argument(
        "--save-distance-matrix",
        action="store_true",
        help="Also save the full square cell-by-cell CNV distance matrix as .tsv.gz. "
             "This can be very large.",
    )
    p.add_argument("--png-dpi", type=int, default=150)

    return p.parse_args()


def main() -> None:
    args = parse_args()

    if args.vmax <= args.vmin:
        raise ValueError("--vmax must be greater than --vmin")
    if not (0 < args.cluster_min_feature_coverage <= 1):
        raise ValueError("--cluster-min-feature-coverage must be in (0, 1]")
    if not (0 < args.cluster_min_pair_overlap <= 1):
        raise ValueError("--cluster-min-pair-overlap must be in (0, 1]")
    if args.cluster_max_bins < 0:
        raise ValueError("--cluster-max-bins must be >= 0")

    # ------------------------------------------------------------------
    # Collect inputs.
    # ------------------------------------------------------------------
    files: List[str] = []
    for pat in args.input:
        files.extend(sorted(glob.glob(pat)) if any(c in pat for c in "*?[") else [pat])
    files = [f for f in files if os.path.isfile(f) and os.path.getsize(f) > 0]
    if not files:
        eprint("No input BED files found.")
        sys.exit(1)

    use_fixed_bins = (args.bin_size > 0) and (not args.by_chrom)

    chrom_sizes: Dict[str, int] = {}
    if args.fai:
        chrom_sizes = {
            c: s for c, s in parse_chrom_sizes(args.fai).items() if c in CHROM_SET
        }

    if use_fixed_bins:
        if not args.fai:
            eprint("--bin-size requires --fai unless --by-chrom is used.")
            sys.exit(1)
        if not chrom_sizes:
            eprint("No chr1-chr22/X/Y chromosome sizes were found in --fai.")
            sys.exit(1)
        columns, column_lengths_arr, layout = make_fixed_bin_layout(chrom_sizes, args.bin_size)
    else:
        if args.bin_size > 0 and args.by_chrom:
            eprint("Note: --by-chrom is set, so --bin-size is ignored.")
        columns = list(CHROM_ORDER)
        if chrom_sizes:
            column_lengths_arr = np.asarray(
                [float(chrom_sizes.get(c, 1)) for c in columns], dtype=np.float64
            )
        else:
            column_lengths_arr = np.ones(len(columns), dtype=np.float64)
        layout = {}
        eprint(
            "Note: chromosome-mean mode is enabled. Fixed bins are recommended for CNV clone clustering."
        )

    # ------------------------------------------------------------------
    # Build cell x genomic-bin matrix.
    # ------------------------------------------------------------------
    t0 = time.time()
    raw_names: List[str] = []
    vectors: List[np.ndarray] = []

    for idx, path in enumerate(files, start=1):
        name = sample_name_from_path(path, args.sample_regex)
        df = load_bed(path)
        if df.empty:
            eprint(f"Warning: {path} has no usable rows; skipping.")
            continue

        if use_fixed_bins:
            vec = bed_to_fixed_bins_fast(
                df,
                bin_size=args.bin_size,
                layout=layout,
                n_columns=len(columns),
            )
        else:
            vec = bed_to_chrom_means_vector(df)

        if not np.isfinite(vec).any():
            eprint(f"Warning: {path} produced no finite CN bins; skipping.")
            continue

        raw_names.append(name)
        vectors.append(vec)

        if idx % 100 == 0 or idx == len(files):
            eprint(f"Loaded {idx}/{len(files)} BED files")

    if not vectors:
        eprint("No data to plot.")
        sys.exit(1)

    sample_names = clean_sample_names(raw_names)

    if args.metadata_tsv:
        metadata = pd.read_csv(args.metadata_tsv, sep="\t", dtype=str).fillna("")
        run_col = next((c for c in ["#Run", "Run"] if c in metadata.columns), None)
        if run_col and "treatment" in metadata.columns:
            treatments = dict(
                zip(metadata[run_col].str.strip(), metadata["treatment"].str.strip())
            )
            sample_names = [
                f"{name}_{treatments[name]}" if treatments.get(name) else name
                for name in sample_names
            ]
        else:
            eprint(
                "Warning: --metadata-tsv was supplied but no (#Run or Run) + treatment columns were found."
            )

    matrix = np.vstack(vectors).astype(np.float32, copy=False)
    mat = pd.DataFrame(matrix, index=sample_names, columns=columns)

    if not args.float_copy_numbers:
        mat = mat.round()
    mat = mat.clip(lower=args.vmin, upper=args.vmax)
    mat = mat.dropna(axis=1, how="all")

    column_lengths = pd.Series(column_lengths_arr, index=columns, dtype=np.float64)
    column_lengths = column_lengths.loc[mat.columns]

    eprint(
        f"Built matrix: {mat.shape[0]} cells x {mat.shape[1]} plotted bins "
        f"in {time.time() - t0:.1f} s"
    )

    if mat.shape[0] < 2:
        raise ValueError("At least two cells are required for hierarchical clustering.")

    # ------------------------------------------------------------------
    # Build a smaller, biologically appropriate clustering matrix.
    # ------------------------------------------------------------------
    t1 = time.time()
    X_cluster, cluster_columns, cluster_weights = prepare_cluster_matrix(
        mat=mat,
        column_lengths=column_lengths,
        include_sex_chromosomes=args.include_sex_chromosomes_in_clustering,
        max_bins=args.cluster_max_bins,
        min_feature_coverage=args.cluster_min_feature_coverage,
        drop_invariant=not args.keep_invariant_bins,
    )

    missing_fraction = float(np.mean(~np.isfinite(X_cluster)))
    eprint(
        f"Clustering matrix: {X_cluster.shape[0]} cells x {X_cluster.shape[1]} features; "
        f"missing={missing_fraction:.3%}; preparation={time.time() - t1:.1f} s"
    )

    # ------------------------------------------------------------------
    # Pairwise CNV distances.
    # ------------------------------------------------------------------
    distance_gb = estimate_condensed_gb(mat.shape[0])
    eprint(f"Condensed pairwise distance vector: ~{distance_gb:.3f} GiB")
    if distance_gb > args.max_distance_gb:
        raise MemoryError(
            f"The condensed distance vector alone is ~{distance_gb:.2f} GiB, exceeding "
            f"--max-distance-gb={args.max_distance_gb}. Exact hierarchical clustering is "
            "quadratic in the number of cells. Increase the limit only if you have enough RAM."
        )

    t2 = time.time()
    distances, distance_backend_used = weighted_cityblock_condensed(
        X=X_cluster,
        weights=cluster_weights,
        min_pair_overlap=args.cluster_min_pair_overlap,
        backend=args.distance_backend,
        numpy_block_size=args.numpy_distance_block_size,
    )

    bad = ~np.isfinite(distances)
    if np.any(bad):
        n_bad = int(bad.sum())
        raise ValueError(
            f"{n_bad} cell pairs have insufficient jointly observed CNV sequence for "
            f"--cluster-min-pair-overlap={args.cluster_min_pair_overlap}. Lower the threshold "
            "or inspect missing CNV calls."
        )

    eprint(
        f"Computed {len(distances):,} pairwise distances with {distance_backend_used} "
        f"in {time.time() - t2:.1f} s"
    )

    # ------------------------------------------------------------------
    # Hierarchical clustering.
    # ------------------------------------------------------------------
    use_optimal_ordering = choose_optimal_ordering(
        args.optimal_ordering,
        n_cells=mat.shape[0],
        max_cells=args.optimal_ordering_max_cells,
    )
    eprint(
        f"Linkage: method={args.linkage_method}; optimal_ordering={use_optimal_ordering}"
    )

    t3 = time.time()
    row_Z = linkage(
        distances,
        method=args.linkage_method,
        optimal_ordering=use_optimal_ordering,
    )
    eprint(f"Hierarchical linkage completed in {time.time() - t3:.1f} s")

    out_dir = os.path.dirname(os.path.abspath(args.output_prefix))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    save_clustering_outputs(
        output_prefix=args.output_prefix,
        Z=row_Z,
        sample_names=list(mat.index),
        distances=distances,
        save_distance_matrix=args.save_distance_matrix,
    )

    # ------------------------------------------------------------------
    # Heatmap display. Missing CN is filled only here, after clustering.
    # ------------------------------------------------------------------
    fill = mat.fillna(args.center)

    if args.float_copy_numbers:
        heatmap_cmap = args.cmap
        heatmap_norm = None
        cbar_kws = {
            "label": "Relative copy-number intensity",
            "orientation": "horizontal",
        }
        heatmap_vmin = args.vmin
        heatmap_vmax = args.vmax
        heatmap_center = args.center
    else:
        # Integer CN palette; vmin/vmax must define integer endpoints here.
        int_vmin = int(round(args.vmin))
        int_vmax = int(round(args.vmax))
        if not np.isclose(args.vmin, int_vmin) or not np.isclose(args.vmax, int_vmax):
            raise ValueError("Integer copy-number mode requires integer --vmin and --vmax values.")

        n_levels = int_vmax - int_vmin + 1
        base = plt.get_cmap(args.cmap, n_levels)
        discrete_cmap = ListedColormap([base(i) for i in range(n_levels)])
        heatmap_norm = BoundaryNorm(
            np.arange(int_vmin - 0.5, int_vmax + 1.5, 1.0),
            discrete_cmap.N,
        )
        heatmap_cmap = discrete_cmap
        cbar_kws = {
            "label": "Copy numbers",
            "ticks": np.arange(int_vmin, int_vmax + 1),
            "spacing": "proportional",
            "orientation": "horizontal",
        }
        heatmap_vmin = None
        heatmap_vmax = None
        heatmap_center = None

    figsize = (
        max(8, min(0.09 * mat.shape[1] + 6, 12)),
        max(8, min(0.18 * mat.shape[0] + 3, 12)),
    )
    eprint(f"figsize={figsize}")

    if args.show_sample_labels and mat.shape[0] > 500:
        eprint(
            "Warning: >500 cells with sample labels may be visually crowded and slow to render. "
            "Use --show-sample-labels 0 for a cleaner/faster figure."
        )

    t4 = time.time()
    g = sns.clustermap(
        fill,
        row_cluster=True,
        col_cluster=False,
        row_linkage=row_Z,
        cmap=heatmap_cmap,
        norm=heatmap_norm,
        vmin=heatmap_vmin,
        vmax=heatmap_vmax,
        center=heatmap_center,
        figsize=figsize,
        cbar_kws=cbar_kws,
        xticklabels=False,
        yticklabels=bool(args.show_sample_labels),
        dendrogram_ratio=(0.15, 0.055),
    )

    g.ax_cbar.set_position([0.25, 0.98, 0.5, 0.01])
    g.ax_cbar.tick_params(axis="x", length=3)

    ax = g.ax_heatmap
    ax.set_xlabel("")
    ax.set_ylabel("")

    # Chromosome dividers and centered labels.
    col_chroms = [chromosome_of_column(c) for c in fill.columns]
    prev = None
    for i, c in enumerate(col_chroms):
        if prev is not None and c != prev:
            ax.axvline(i, color="black", linestyle="--", linewidth=0.8)
        prev = c

    pos: Dict[str, List[int]] = {}
    for i, c in enumerate(col_chroms):
        pos.setdefault(c, []).append(i)

    centers: List[float] = []
    labels: List[str] = []
    for c in CHROM_ORDER:
        if c in pos:
            centers.append((pos[c][0] + pos[c][-1] + 1) / 2.0)
            labels.append(c.replace("chr", ""))

    ax.set_xticks(centers)
    ax.set_xticklabels(labels, rotation=15, fontsize=8)
    ax.tick_params(axis="x", length=0)

    if args.show_sample_labels:
        plt.setp(ax.get_yticklabels(), rotation=0, fontsize=6)

    if args.title:
        g.fig.suptitle(args.title, y=1.02)

    g.savefig(args.output_prefix + ".pdf", bbox_inches="tight")
    g.savefig(args.output_prefix + ".png", dpi=args.png_dpi, bbox_inches="tight")
    plt.close("all")

    eprint(f"Rendering completed in {time.time() - t4:.1f} s")
    eprint(
        "Wrote:\n"
        f"  {args.output_prefix}.pdf\n"
        f"  {args.output_prefix}.png\n"
        f"  {args.output_prefix}.linkage.tsv\n"
        f"  {args.output_prefix}.cluster_order.tsv"
        + (f"\n  {args.output_prefix}.cnv_distance.tsv.gz" if args.save_distance_matrix else "")
    )


if __name__ == "__main__":
    main()

