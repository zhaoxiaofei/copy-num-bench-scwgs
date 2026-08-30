#!/usr/bin/env python3
"""
https://chatgpt.com/c/6a947653-f9f0-83ea-9bac-a97e08c5b168

Final-manuscript plotting front-end for scWGS-ploidy-performances-eval_v02.py.

Place this file beside ``scWGS-ploidy-performances-eval_v02.py`` in
``bench_results/`` and invoke it with the same CLI arguments as v02.  Data
loading, filtering, row ordering, benchmark calculations, and TSV export are
still performed by v02; this file replaces only the main-text plotting layer.

Visual encoding
---------------
* Balloon DIAMETER = percentage of cells whose ploidy estimate is within ±0.5
  of the ground truth.
* Balloon COLOUR INTENSITY = log10(sample size), by default log10 of the
  number of finite/evaluable cells contributing to that percentage
  (``n_cells_finite``).
* Red × = a result was expected but is missing / failed.
* Grey hatched block = method is not applicable to that panel.

Method-specific rules
---------------------
* scAbsolute: shown once, because copy-number cap vs no-cap does not apply; it
  produces a final ploidy estimate rather than a CNV profile.
* CHISEL in ACT: shown as not applicable because phased genotypes are not
  readily available for the ACT samples.

The combined 2×2 figure follows common manuscript conventions: lowercase bold
panel letters (a–d) at the upper-left outside each axes, no embedded figure
headline/caption, fully horizontal, vertically interleaved method labels, explicit CapAt10 labels, and one shared legend below the panels, and collision-free two-tier x-axis labels.
"""

from __future__ import annotations

import importlib.util
import os
import sys
from collections import OrderedDict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Patch, Rectangle


# ======================================================================================
# User-facing style/configuration knobs
# ======================================================================================
COUNT_FIELD = "n_cells_finite"        # use "n_cells" for total loaded cells instead
COUNT_LEGEND_LABEL = r"$\log_{10}$(evaluable cells, n)"
PLOIDY_TOLERANCE_TEXT = "±0.5"

# Academic figures usually place the explanatory title/caption outside the figure.
SHOW_COMBINED_SUPTITLE = False
SHOW_COMBINED_FOOTNOTE = False
SHOW_INDIVIDUAL_FOOTNOTE = False

# Perceptually ordered, colour-blind-safe-ish sequential blue scale.
COUNT_CMAP = LinearSegmentedColormap.from_list(
    "sample_size_blues",
    ["#eff6ff", "#c6dbef", "#9ecae1", "#6baed6", "#3182bd", "#08519c", "#08306b"],
)

TEXT = "#20262B"
MUTED = "#66717A"
GRID = "#E2E7EB"
GROUP_GRID = "#B8C1C8"
ALT_BAND = "#F7F9FA"
EDGE = "#21384E"
MISSING = "#D84A4A"
NA_FACE = "#F0F2F4"
NA_EDGE = "#A5AEB5"

# Methods for which cap/no-cap is conceptually not applicable.
SINGLE_COLUMN_TOOLS = {"scabsolute"}

# Entire caller is not applicable in these panels.
NOT_APPLICABLE_BY_PLOT = {
    "ACT": {"chisel"},
}

# Marker diameter (points).  We square these before passing to scatter because
# matplotlib's ``s`` is marker area in pt^2.  Thus the user-visible *diameter*
# is linear in the percentage, exactly as requested.
D_MIN_INDIV = 2.2
D_MAX_INDIV = 8.0
D_MIN_COMBINED = 1.5
D_MAX_COMBINED = 6.0

# The actual maximum diameter is also limited by the densest panel row pitch so
# 100% circles can never touch the rows above/below after manuscript reduction.
MAX_DIAMETER_ROW_FRACTION = 0.72

# Designed around a full-width/two-column manuscript figure.  Vector PDF output
# remains crisp if the publisher rescales it modestly.
COMBINED_FIGSIZE = (8.25, 9.25)
INDIVIDUAL_FIGSIZE = (8.1, 5.9)

# Cache display payloads because the original v02 passes only the percentage
# matrix to its combined-figure function.
_PANEL_CACHE = {}


# ======================================================================================
# Base-module loader and style
# ======================================================================================
def _load_base_module():
    here = os.path.dirname(os.path.abspath(__file__))
    base_path = os.path.join(here, "scWGS-ploidy-performances-eval_v02.py")
    if not os.path.isfile(base_path):
        raise SystemExit(
            "Could not find scWGS-ploidy-performances-eval_v02.py next to this file.\n"
            "Place this script in the same bench_results directory as v02."
        )
    spec = importlib.util.spec_from_file_location("ploidy_v02_base", base_path)
    if spec is None or spec.loader is None:
        raise SystemExit("Could not import the v02 script: " + base_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _apply_style():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "Liberation Sans"],
        "font.size": 7.0,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.55,
        "axes.edgecolor": GROUP_GRID,
        "axes.labelcolor": TEXT,
        "xtick.color": TEXT,
        "ytick.color": TEXT,
        "text.color": TEXT,
        "xtick.major.size": 0,
        "ytick.major.size": 0,
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
    })


# ======================================================================================
# Data/display helpers
# ======================================================================================
def _pct_to_area(pct, d_min, d_max):
    if not np.isfinite(pct):
        return 0.0
    p = np.clip(float(pct), 0.0, 100.0) / 100.0
    diameter = d_min + (d_max - d_min) * p
    return float(diameter ** 2)


def _count_to_log10(n):
    """Transform a positive cell count to log10(n); non-positive values are NaN."""
    try:
        x = float(n)
    except Exception:
        return float("nan")
    if not np.isfinite(x) or x <= 0:
        return float("nan")
    return float(np.log10(x))


def _count_norm(entries):
    """Shared linear normalization *after* log10-transforming cell counts."""
    if entries is None or COUNT_FIELD not in entries.columns:
        return Normalize(vmin=0.0, vmax=1.0)
    raw = np.asarray(entries[COUNT_FIELD], dtype=float)
    raw = raw[np.isfinite(raw) & (raw > 0)]
    if raw.size == 0:
        return Normalize(vmin=0.0, vmax=1.0)
    logx = np.log10(raw)
    vmax = float(np.max(logx))
    # Start at log10(1)=0 to keep the colour scale interpretable and stable.
    if not np.isfinite(vmax) or vmax <= 0:
        vmax = 1.0
    return Normalize(vmin=0.0, vmax=vmax)


def _nice_log_ticks(norm):
    """Integer log10 ticks (0,1,2,...) within the shared colour scale."""
    hi = float(norm.vmax)
    ticks = np.arange(0, int(np.floor(hi)) + 1, dtype=float)
    if ticks.size == 0:
        ticks = np.array([0.0])
    return ticks


def _limit_diameter_to_rows(fig, panel_height_frac, n_rows, requested_dmax):
    """Prevent circles from touching neighbouring rows while keeping one scale."""
    if n_rows <= 0:
        return requested_dmax
    row_pitch_pt = panel_height_frac * fig.get_figheight() * 72.0 / float(n_rows)
    return min(float(requested_dmax), MAX_DIAMETER_ROW_FRACTION * row_pitch_pt)


def _choose_single_cap(pairs):
    """Prefer cap=10 for scAbsolute if both duplicate variants exist."""
    if not pairs:
        return None
    for _tool, cap in pairs:
        try:
            if float(cap) == 10.0:
                return cap
        except Exception:
            pass
    return pairs[0][1]


def _cap_state_label(cap):
    """Publication label for a CN-calling variant.

    The repository defines the default copy-number cap as 10, so the capped
    variant is labelled ``CapAt10`` rather than the symbolic ``<=10``.
    """
    try:
        x = float(cap)
        if np.isfinite(x) and np.isclose(x, 10.0):
            return "CapAt10"
        if not np.isfinite(x):
            return "NC"
    except Exception:
        pass
    text = str(cap).strip().lower()
    if text in {"inf", "+inf", "infinity", "none", "no cap", "nocap", "nc"}:
        return "NC"
    # In normal v02 output only 10 and infinity occur.  Retain an informative
    # fallback for user-filtered/legacy inputs rather than silently mislabelling.
    try:
        return f"CapAt{float(cap):g}"
    except Exception:
        return str(cap)


def _display_colspecs(plot_id, raw_col_pairs, base):
    """Build the *displayed* columns from the raw tool×cap column pairs."""
    grouped = OrderedDict()
    for tool, cap in raw_col_pairs:
        grouped.setdefault(tool, []).append((tool, cap))

    displayed = []
    na_tools = set(NOT_APPLICABLE_BY_PLOT.get(plot_id, set()))

    for tool, pairs in grouped.items():
        if tool in SINGLE_COLUMN_TOOLS:
            displayed.append({
                "tool": tool,
                "cap": _choose_single_cap(pairs),
                "cap_label": "",
                "applicable": True,
            })
            continue

        if tool in na_tools:
            # Keep the pair width so ACT remains directly comparable with the other
            # panels.  The whole pair is rendered as one hatched NA block later.
            if pairs:
                for _t, cap in pairs:
                    displayed.append({
                        "tool": tool,
                        "cap": cap,
                        "cap_label": _cap_state_label(cap),
                        "applicable": False,
                    })
            else:
                displayed.append({
                    "tool": tool,
                    "cap": None,
                    "cap_label": "",
                    "applicable": False,
                })
            continue

        for _t, cap in pairs:
            displayed.append({
                "tool": tool,
                "cap": cap,
                "cap_label": _cap_state_label(cap),
                "applicable": True,
            })

    return displayed


def _column_groups(colspecs):
    """Return [(tool, left_col, right_col, applicable), ...]."""
    if not colspecs:
        return []
    groups = []
    start = 0
    tool = colspecs[0]["tool"]
    applicable = colspecs[0]["applicable"]
    for j in range(1, len(colspecs) + 1):
        if j == len(colspecs) or colspecs[j]["tool"] != tool:
            groups.append((tool, start, j - 1, applicable))
            if j < len(colspecs):
                start = j
                tool = colspecs[j]["tool"]
                applicable = colspecs[j]["applicable"]
    return groups


def _build_payload(entries, plot_id, row_labels, raw_col_pairs, base):
    colspecs = _display_colspecs(plot_id, raw_col_pairs, base)
    sub = entries[entries["plot"] == plot_id]

    lookup = {}
    for rec in sub.itertuples(index=False):
        if rec.dataset is None or (isinstance(rec.dataset, float) and np.isnan(rec.dataset)):
            continue
        lookup[(rec.dataset, rec.tool, rec.max_cn)] = rec

    nr, nc = len(row_labels), len(colspecs)
    pct = np.full((nr, nc), np.nan, dtype=float)
    count = np.full((nr, nc), np.nan, dtype=float)
    status = np.full((nr, nc), "missing", dtype=object)  # valid / missing / na

    for i, dataset in enumerate(row_labels):
        for j, spec in enumerate(colspecs):
            if not spec["applicable"]:
                status[i, j] = "na"
                continue

            rec = lookup.get((dataset, spec["tool"], spec["cap"]))
            if rec is None:
                status[i, j] = "missing"
                continue

            n = getattr(rec, COUNT_FIELD, np.nan)
            if np.isfinite(n):
                count[i, j] = float(n)

            if (not rec.failed) and np.isfinite(rec.pct_within):
                pct[i, j] = float(rec.pct_within)
                status[i, j] = "valid"
            else:
                status[i, j] = "missing"

    payload = {
        "plot_id": plot_id,
        "row_labels": list(row_labels),
        "raw_col_pairs": list(raw_col_pairs),
        "colspecs": colspecs,
        "pct": pct,
        "count": count,
        "status": status,
    }
    _PANEL_CACHE[plot_id] = payload
    return payload


# ======================================================================================
# Figure annotations / legends
# ======================================================================================
def _panel_letter(fig, ax, letter, fontsize=13.5):
    """Large bold lowercase panel label, just outside the upper-left axes corner."""
    if not letter:
        return
    bbox = ax.get_position()
    fig.text(
        bbox.x0 - 0.030,
        bbox.y1 + 0.008,
        str(letter).lower(),
        ha="left",
        va="bottom",
        fontsize=fontsize,
        fontweight="bold",
        color=TEXT,
    )


def _size_legend(ax, d_min, d_max, values=(0, 25, 50, 75, 100), fontsize=6.3):
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.0, 0.98, f"Cells within {PLOIDY_TOLERANCE_TEXT} of truth",
            ha="left", va="top", fontsize=fontsize + 0.3, fontweight="bold", color=TEXT)
    ys = np.linspace(0.80, 0.24, len(values))
    ref_colour = COUNT_CMAP(0.55)
    for y, v in zip(ys, values):
        ax.scatter([0.20], [y], s=_pct_to_area(v, d_min, d_max),
                   facecolor=ref_colour, edgecolor=EDGE, linewidth=0.32)
        ax.text(0.43, y, f"{v}%", ha="left", va="center", fontsize=fontsize, color=TEXT)


def _status_legend(ax, fontsize=6.3):
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.scatter([0.18], [0.68], marker="x", s=28, color=MISSING, linewidth=1.0)
    ax.text(0.38, 0.68, "Missing result", va="center", fontsize=fontsize, color=TEXT)

    rect = Rectangle((0.10, 0.20), 0.16, 0.22, facecolor=NA_FACE,
                     edgecolor=NA_EDGE, hatch="///", linewidth=0.6)
    ax.add_patch(rect)
    ax.text(0.38, 0.31, "Not applicable", va="center", fontsize=fontsize, color=TEXT)
    ax.text(0.10, -0.03, "CapAt10 = CN cap at 10;  NC = no cap",
            va="top", ha="left", fontsize=max(fontsize - 0.35, 4.8), color=MUTED)


def _legend_card(ax, title, fontsize=6.2):
    """Style a small boxed legend region with a clear section title."""
    ax.set_facecolor("#FCFDFE")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#CDD5DB")
        spine.set_linewidth(0.65)
    ax.text(0.045, 0.91, title, ha="left", va="top",
            fontsize=fontsize, fontweight="bold", color=TEXT)


def _colourbar_in_card(card_ax, norm, fontsize=5.8):
    cax = card_ax.inset_axes([0.08, 0.34, 0.84, 0.20])
    sm = ScalarMappable(norm=norm, cmap=COUNT_CMAP)
    sm.set_array([])
    cb = card_ax.figure.colorbar(sm, cax=cax, orientation="horizontal")
    cb.outline.set_linewidth(0.45)
    cb.outline.set_edgecolor(GROUP_GRID)
    cb.ax.tick_params(labelsize=fontsize, length=1.8, width=0.45, color=TEXT, pad=1.2)
    cb.set_ticks(_nice_log_ticks(norm))
    cb.set_label(COUNT_LEGEND_LABEL, fontsize=fontsize + 0.1, color=TEXT, labelpad=2.0)
    return cb


def _size_legend_horizontal(ax, d_min, d_max, values=(0, 25, 50, 75, 100), fontsize=6.0, embedded=False):
    """Compact horizontal size legend for the shared combined-figure legend."""
    if not embedded:
        ax.set_axis_off()
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.text(0.0, 0.95, f"Cells within {PLOIDY_TOLERANCE_TEXT} of truth",
                ha="left", va="top", fontsize=fontsize + 0.25, fontweight="bold", color=TEXT)
        y_circle, y_text = 0.53, 0.16
    else:
        ax.text(0.045, 0.72, f"Cells within {PLOIDY_TOLERANCE_TEXT} of truth",
                ha="left", va="top", fontsize=fontsize, color=MUTED)
        y_circle, y_text = 0.39, 0.13
    xs = np.linspace(0.09, 0.91, len(values))
    ref_colour = COUNT_CMAP(0.55)
    for x, v in zip(xs, values):
        ax.scatter([x], [y_circle], s=_pct_to_area(v, d_min, d_max),
                   facecolor=ref_colour, edgecolor=EDGE, linewidth=0.30)
        ax.text(x, y_text, f"{v}%", ha="center", va="center", fontsize=fontsize, color=TEXT)


def _status_legend_horizontal(ax, fontsize=6.0, embedded=False):
    if not embedded:
        ax.set_axis_off()
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        y = 0.68
    else:
        y = 0.56
    ax.scatter([0.09], [y], marker="x", s=25, color=MISSING, linewidth=0.95)
    ax.text(0.17, y, "Missing", va="center", fontsize=fontsize, color=TEXT)
    rect = Rectangle((0.51, y - 0.10), 0.08, 0.20, facecolor=NA_FACE,
                     edgecolor=NA_EDGE, hatch="///", linewidth=0.55)
    ax.add_patch(rect)
    ax.text(0.63, y, "Not applicable", va="center", fontsize=fontsize, color=TEXT)
    ax.text(0.05, 0.16, "CapAt10 = CN cap at 10    |    NC = no cap",
            va="center", ha="left", fontsize=max(fontsize - 0.25, 4.8), color=MUTED)


def _colourbar_horizontal(fig, rect, norm, fontsize=6.0):
    cax = fig.add_axes(rect)
    sm = ScalarMappable(norm=norm, cmap=COUNT_CMAP)
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.outline.set_linewidth(0.5)
    cb.outline.set_edgecolor(GROUP_GRID)
    cb.ax.tick_params(labelsize=fontsize, length=2.0, width=0.5, color=TEXT, pad=1.5)
    cb.set_ticks(_nice_log_ticks(norm))
    cb.set_label(COUNT_LEGEND_LABEL, fontsize=fontsize + 0.2, color=TEXT, labelpad=2.5)
    return cb


def _colourbar(fig, rect, norm, fontsize=6.2):
    cax = fig.add_axes(rect)
    sm = ScalarMappable(norm=norm, cmap=COUNT_CMAP)
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cax)
    cb.outline.set_linewidth(0.5)
    cb.outline.set_edgecolor(GROUP_GRID)
    cb.ax.tick_params(labelsize=fontsize, length=2.0, width=0.5, color=TEXT, pad=2)
    cb.set_ticks(_nice_log_ticks(norm))
    cb.set_label(COUNT_LEGEND_LABEL, fontsize=fontsize + 0.2, color=TEXT, labelpad=4)
    return cb


# ======================================================================================
# Plotting overrides
# ======================================================================================
def install_plotting_overrides(base):
    def build_matrix(entries, plot_id, row_labels, col_pairs):
        # base.main() only needs a matrix-like object to decide panel availability;
        # the actual displayed payload is cached and used by our custom plotters.
        payload = _build_payload(entries, plot_id, row_labels, col_pairs, base)
        return payload["pct"]

    def draw_panel(ax, payload, norm, d_min, d_max, *,
                   method_font=5.2, row_font=5.1, cap_font=4.35,
                   method_rotation=0.0, show_ylabels=True,
                   stagger_method_labels=False):
        pct = payload["pct"]
        count = payload["count"]
        status = payload["status"]
        rows = payload["row_labels"]
        cols = payload["colspecs"]
        groups = _column_groups(cols)
        nr, nc = pct.shape

        # Alternating caller bands improve pair grouping without a boxed-table look.
        for gi, (_tool, left, right, applicable) in enumerate(groups):
            if not applicable:
                ax.axvspan(left - 0.5, right + 0.5, facecolor=NA_FACE,
                           edgecolor=NA_EDGE, linewidth=0.0, hatch="///", zorder=0)
            elif gi % 2 == 1:
                ax.axvspan(left - 0.5, right + 0.5, facecolor=ALT_BAND,
                           edgecolor="none", zorder=0)

        # Quiet grid.
        for y in np.arange(-0.5, nr, 1.0):
            ax.axhline(y, color=GRID, linewidth=0.36, zorder=0.5)
        for x in np.arange(-0.5, nc, 1.0):
            ax.axvline(x, color=GRID, linewidth=0.28, zorder=0.5)
        for _tool, left, _right, _app in groups[1:]:
            ax.axvline(left - 0.5, color=GROUP_GRID, linewidth=0.65, zorder=1)

        xs, ys, sizes, colours = [], [], [], []
        mx, my = [], []
        for i in range(nr):
            for j in range(nc):
                if status[i, j] == "valid":
                    xs.append(j)
                    ys.append(i)
                    sizes.append(_pct_to_area(pct[i, j], d_min, d_max))
                    n = count[i, j]
                    logn = _count_to_log10(n)
                    colours.append(COUNT_CMAP(norm(0.0 if not np.isfinite(logn) else logn)))
                elif status[i, j] == "missing":
                    mx.append(j)
                    my.append(i)
                # status == 'na' is intentionally represented by the group-wide hatch.

        if xs:
            ax.scatter(xs, ys, s=sizes, c=colours, edgecolor=EDGE,
                       linewidth=0.28, zorder=3)
        if mx:
            ax.scatter(mx, my, marker="x", s=19, color=MISSING,
                       linewidth=0.85, zorder=4)

        # Not-applicable callers are indicated by the hatched band only; the
        # shared legend explains the hatch.  Avoiding text inside the data grid
        # prevents collisions and keeps the ACT panel visually quiet.

        ax.set_xlim(-0.5, nc - 0.5)
        ax.set_ylim(nr - 0.5, -0.5)

        ax.set_yticks(np.arange(nr))
        if show_ylabels:
            ax.set_yticklabels([base.wrap_label(str(r), width=34) for r in rows], fontsize=row_font)
            ax.tick_params(axis="y", pad=2.5)
        else:
            # Panels a-c share the same germline-derived row ordering, so the
            # upper-right panel can omit repeated labels.  This is standard
            # multi-panel practice and prevents labels from spilling into panel a.
            ax.set_yticklabels([])
            ax.tick_params(axis="y", pad=0.0)

        # Two-tier x axis.  Cap-state labels are vertical because ``CapAt10``
        # is necessarily much wider than a subcolumn.  Method names are the
        # semantic x-axis labels; they are fully horizontal and alternate
        # between two vertical levels so adjacent long names remain readable.
        ax.set_xticks(np.arange(nc))
        ax.set_xticklabels([])
        ax.tick_params(axis="x", pad=0.0)
        for j, spec in enumerate(cols):
            if not spec["cap_label"]:
                continue
            ax.text(j, -0.014, spec["cap_label"],
                    transform=ax.get_xaxis_transform(),
                    ha="center", va="top", rotation=90,
                    fontsize=cap_font, clip_on=False, color=MUTED)

        for gi, (tool, left, right, _applicable) in enumerate(groups):
            center = 0.5 * (left + right)
            # A shallow two-level stagger is used only in the combined figure.
            # It preserves the requested near-horizontal labels while ensuring
            # adjacent long names never collide after manuscript-scale reduction.
            if stagger_method_labels:
                y_method = -0.158 if gi % 2 == 0 else -0.205
            else:
                y_method = -0.165
            # scAbsolute occupies one column rather than a two-column pair; a
            # tiny left offset gives it visual breathing room from HMMcopy.
            if tool == "scabsolute":
                center -= 0.10
            ax.text(center, y_method, base.pretty_tool(tool),
                    transform=ax.get_xaxis_transform(),
                    ha="center", va="top", rotation=0,
                    fontsize=method_font,
                    clip_on=False, color=TEXT)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color(GROUP_GRID)
        ax.spines["bottom"].set_color(GROUP_GRID)
        ax.spines["left"].set_linewidth(0.55)
        ax.spines["bottom"].set_linewidth(0.55)
        ax.set_facecolor("white")
        return ax

    def plot_one_balloon(entries, plot_id, title, row_labels, col_pairs, args,
                         panel_letter=None, s_min=None, s_max=None):
        payload = _build_payload(entries, plot_id, row_labels, col_pairs, base)
        norm = _count_norm(entries)
        scale = np.sqrt(max(float(args.dot_scale), 1e-6))
        d_min, d_max = D_MIN_INDIV * scale, D_MAX_INDIV * scale

        fig = plt.figure(figsize=INDIVIDUAL_FIGSIZE)
        ax_rect = [0.155, 0.275, 0.675, 0.605]
        d_max = _limit_diameter_to_rows(fig, ax_rect[3], len(row_labels), d_max)
        d_min = min(d_min, 0.34 * d_max)
        ax = fig.add_axes(ax_rect)
        draw_panel(ax, payload, norm, d_min, d_max,
                   method_font=6.0, row_font=6.2, cap_font=5.0, method_rotation=0,
                   show_ylabels=True, stagger_method_labels=False)
        _panel_letter(fig, ax, panel_letter, fontsize=14.5)
        ax.set_title(title, loc="left", fontsize=8.8, fontweight="bold", pad=6, color=TEXT)
        ax.set_ylabel("Dataset", fontsize=7.0, labelpad=6)

        # Compact legend block at right.
        size_ax = fig.add_axes([0.848, 0.55, 0.14, 0.27])
        _size_legend(size_ax, d_min, d_max, fontsize=6.2)
        status_ax = fig.add_axes([0.848, 0.405, 0.14, 0.105])
        _status_legend(status_ax, fontsize=6.2)
        _colourbar(fig, [0.905, 0.225, 0.014, 0.15], norm, fontsize=6.1)

        if SHOW_INDIVIDUAL_FOOTNOTE:
            fig.text(0.155, 0.035,
                     "Balloon diameter: accuracy; colour: sample size; red ×: missing; hatched: not applicable.",
                     fontsize=6.0, color=MUTED)
        return fig

    def plot_combined_four(figs_spec, args):
        nonempty = [s for s in figs_spec if s[2] and s[3] is not None]
        if not nonempty:
            return None

        entries = getattr(base, "_LAST_ENTRIES_FOR_REVISED", None)
        norm = _count_norm(entries)
        scale = np.sqrt(max(float(args.dot_scale), 1e-6))
        d_min, d_max = D_MIN_COMBINED * scale, D_MAX_COMBINED * scale

        fig = plt.figure(figsize=COMBINED_FIGSIZE)

        # Use nearly the full page width for the panels.  The shared legend is
        # moved below the 2x2 grid instead of consuming a right-hand column.
        # This gives long method names (CopyNumber, AneuFinder, scAbsolute, ...)
        # enough horizontal room to remain only slightly tilted without overlap.
        left, right = 0.082, 0.988
        bottom, top = 0.245, 0.965
        wspace, hspace = 0.095, 0.155
        panel_w = (right - left - wspace) / 2.0
        panel_h = (top - bottom - hspace) / 2.0

        max_rows = max((len(s[2]) for s in nonempty), default=1)
        d_max = _limit_diameter_to_rows(fig, panel_h, max_rows, d_max)
        d_min = min(d_min, 0.34 * d_max)

        positions = [
            [left,                  bottom + panel_h + hspace, panel_w, panel_h],
            [left + panel_w+wspace, bottom + panel_h + hspace, panel_w, panel_h],
            [left,                  bottom,                    panel_w, panel_h],
            [left + panel_w+wspace, bottom,                    panel_w, panel_h],
        ]

        letters = "abcd"
        for k, (plot_id, title, row_labels, matrix, col_pairs) in enumerate(
                base.PLOT_SPECS_TO_SPEC(figs_spec)):
            ax = fig.add_axes(positions[k])

            payload = _PANEL_CACHE.get(plot_id)
            if payload is None and row_labels:
                payload = _build_payload(entries, plot_id, row_labels, col_pairs, base)

            if payload is None or not row_labels:
                ax.set_axis_off()
                _panel_letter(fig, ax, letters[k], fontsize=15.5)
                ax.text(0.0, 1.01, title, transform=ax.transAxes,
                        fontsize=7.3, fontweight="bold", ha="left", va="bottom", color=MUTED)
                ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                        ha="center", va="center", fontsize=6.5, color=MUTED)
                continue

            # b repeats exactly the germline-derived row labels shown in a;
            # suppressing them prevents center-gutter overlap.  d keeps its
            # distinct ACT sample labels, at a slightly smaller size.
            show_y = (k != 1)
            row_fs = 4.55 if k != 3 else 4.15
            draw_panel(ax, payload, norm, d_min, d_max,
                       method_font=5.05, row_font=row_fs, cap_font=3.8,
                       method_rotation=0, show_ylabels=show_y,
                       stagger_method_labels=True)
            _panel_letter(fig, ax, letters[k], fontsize=15.5)
            ax.set_title(title, loc="left", fontsize=7.35, fontweight="bold", pad=5, color=TEXT)

        # Shared legend strip: three visually separated cards.  This avoids the
        # previous "clustered" appearance and makes each visual encoding obvious.
        acc_ax = fig.add_axes([0.082, 0.040, 0.285, 0.125])
        _legend_card(acc_ax, "Balloon diameter")
        _size_legend_horizontal(acc_ax, d_min, d_max, fontsize=5.9, embedded=True)

        stat_ax = fig.add_axes([0.390, 0.040, 0.240, 0.125])
        _legend_card(stat_ax, "Result status")
        _status_legend_horizontal(stat_ax, fontsize=5.9, embedded=True)

        count_ax = fig.add_axes([0.653, 0.040, 0.335, 0.125])
        _legend_card(count_ax, "Balloon colour")
        _colourbar_in_card(count_ax, norm, fontsize=5.8)

        if SHOW_COMBINED_SUPTITLE:
            fig.suptitle("Ploidy-estimation accuracy", fontsize=9.2, fontweight="bold", y=0.993)
        if SHOW_COMBINED_FOOTNOTE:
            fig.text(left, 0.018,
                     "Diameter = % cells within ±0.5 of truth; colour = evaluable-cell count; "
                     "red × = missing; hatch = not applicable.",
                     fontsize=5.8, color=MUTED)
        return fig

    # Install our replacements into v02.
    base._apply_style = _apply_style
    base.build_matrix = build_matrix
    base.plot_one_balloon = plot_one_balloon
    base.plot_combined_four = plot_combined_four

    # Store entries so the combined figure can rebuild/cache payloads if needed.
    if not hasattr(base, "_LAST_ENTRIES_FOR_REVISED"):
        base._LAST_ENTRIES_FOR_REVISED = None
    original_expand = base.expand_runs_to_entries

    def expand_runs_to_entries_and_store(*args, **kwargs):
        out = original_expand(*args, **kwargs)
        base._LAST_ENTRIES_FOR_REVISED = out
        return out

    base.expand_runs_to_entries = expand_runs_to_entries_and_store


# ======================================================================================
# Entrypoint
# ======================================================================================
def main(argv=None):
    base = _load_base_module()
    install_plotting_overrides(base)
    return base.main(argv)


if __name__ == "__main__":
    sys.exit(main())

