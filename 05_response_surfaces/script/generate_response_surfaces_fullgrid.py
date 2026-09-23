"""
Generate the INFORM-HIV response-surface figures from the surrogate grid.

Reads  : 04_cost_mapping/output/response_surface_grid_full.csv
         (written by 04_cost_mapping/script/export_response_surface_grid.R)
Writes : the three PNGs INTO THIS DIRECTORY (the delivery folder), not the CWD.

Why this file exists
--------------------
The previous version of this script drew a TRANSPOSED surface: it sorted the CSV
by ['prep_red', 'art_red'] -- so the flat array varied fastest in ART -- and then
called z.reshape(ny, nx), which makes the rows PrEP, while contourf(X, Y, Z)
with meshgrid(prep, art) expects rows = ART. The drawn field was therefore the
truth mirrored about the diagonal (max 1.31 IRR; 67% of the plane off by >0.25).
This version PIVOTS the long table into the array (never reshapes it) and asserts
the result node-by-node before plotting.

Conventions (must match the CSV columns; state them on the figure)
-----------------------------------------------------------------
* time summary : TIME_SUMMARY = "horizon_mean"  -> mean over ticks of the
  per-draw ratio, the convention Figure 1 (forest) uses. The CSV also carries
  "tick10" for the scenario-table convention.
* insulation   : one `case` per run (the committed cost model by default). If the
  CSV carries several cases this script refuses to guess and stops.
* Panel A = intervention USE reduction (identity mapping); Panel B = GOVERNMENT
  FUNDING reduction (Delta mapped to coverage through the cost model, then
  queried). Panel B must differ from Panel A -- asserted below.

Run:  python3 generate_response_surfaces_fullgrid.py
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------------
CSV_PATH = Path(
    "/Users/r_vardavas/Documents/Projects_2026/HIV-UChicago/"
    "INFORM-HIV-Funding-Scenario-Exploration/04_cost_mapping/output/"
    "response_surface_grid_full.csv")
OUT_DIR = Path(__file__).resolve().parent        # the delivery folder
TIME_SUMMARY = "horizon_mean"
CASE = None            # None -> auto-detect; must be unambiguous

PANEL_A = "A. Intervention use reduction"
PANEL_B = "B. Government funding reduction"

# Slices are drawn only at policy-relevant levels: 101 lines with 101 legend
# entries is unreadable (the previous version produced a 3473 x 5894 px canvas
# for a figsize=(14, 5) figure). 0.0 is kept as the no-cut reference.
SLICE_LEVELS = [0.0, 0.10, 0.25, 0.40, 0.75]
POLICY_LEVELS = [0.10, 0.25, 0.40]      # "a few values that represent policies"

CMAP = plt.get_cmap("YlOrRd")
SLICE_CMAP = plt.get_cmap("viridis")
# ---------------------------------------------------------------------------
# LOAD + ORIENT (pivot, never reshape)
# ---------------------------------------------------------------------------
def load_surface():
    df = pd.read_csv(CSV_PATH)
    need = {"panel", "case", "time_summary", "art_red", "prep_red",
            "art_cov", "prep_cov", "mean", "lower", "upper"}
    missing = need - set(df.columns)
    if missing:
        raise SystemExit(f"{CSV_PATH} is missing columns {sorted(missing)}; "
                         "re-run export_response_surface_grid.R")

    cases = sorted(df["case"].unique())
    case = cases[0] if CASE is None else CASE
    if len(cases) > 1 and CASE is None:
        raise SystemExit(f"CSV holds several gamma cases {cases}; pass CASE=<name> "
                         "explicitly so denominators cannot be mixed silently.")
    ts = sorted(df["time_summary"].unique())
    if TIME_SUMMARY not in ts:
        raise SystemExit(f"time_summary {TIME_SUMMARY!r} not in CSV ({ts})")

    sub = df[(df["case"] == case) & (df["time_summary"] == TIME_SUMMARY)]
    art_nodes = np.sort(sub["art_red"].unique())
    prep_nodes = np.sort(sub["prep_red"].unique())

    out = {}
    for panel in (PANEL_A, PANEL_B):
        p = sub[sub["panel"] == panel]
        if p.empty:
            raise SystemExit(f"panel {panel!r} absent for case={case}")
        grids = {}
        for col in ("mean", "lower", "upper", "art_cov", "prep_cov"):
            z = p.pivot_table(index="art_red", columns="prep_red", values=col,
                              aggfunc="first")
            z = z.reindex(index=art_nodes, columns=prep_nodes)
            if z.isna().any().any():
                raise SystemExit(f"{panel}: {col} does not fill a complete "
                                 f"{len(art_nodes)}x{len(prep_nodes)} grid")
            grids[col] = z.values                      # rows = ART, cols = PrEP
        out[panel] = grids

    # --- orientation contract -------------------------------------------------
    # Z[i, j] must be the value at ART = art_nodes[i], PrEP = prep_nodes[j].
    # This is exactly the check the reshape version failed.
    for panel, g in out.items():
        for i in (0, len(art_nodes) // 2, len(art_nodes) - 1):
            for j in (0, len(prep_nodes) // 3, len(prep_nodes) - 1):
                row = sub[(sub["panel"] == panel)
                          & (sub["art_red"] == art_nodes[i])
                          & (sub["prep_red"] == prep_nodes[j])]
                if not np.isclose(g["mean"][i, j], row["mean"].iloc[0],
                                  rtol=0, atol=1e-12):
                    raise SystemExit(
                        f"orientation contract failed for {panel} at "
                        f"(ART={art_nodes[i]}, PrEP={prep_nodes[j]}): "
                        f"{g['mean'][i, j]} vs {row['mean'].iloc[0]}")

    # --- asymmetric-node check (equal-cut nodes cannot detect an axis swap) ---
    def at(g, a, p_):
        i = int(np.argmin(np.abs(art_nodes - a)))
        j = int(np.argmin(np.abs(prep_nodes - p_)))
        return g["mean"][i, j]

    a1, a2 = at(out[PANEL_A], 0.40, 0.10), at(out[PANEL_A], 0.10, 0.40)
    if abs(a1 - a2) < 0.05:
        raise SystemExit("surface looks symmetric at asymmetric nodes -- "
                         "suspect an axis swap")
    b1 = at(out[PANEL_B], 0.40, 0.40)
    a3 = at(out[PANEL_A], 0.40, 0.40)
    if np.isclose(b1, a3, rtol=0, atol=1e-9):
        raise SystemExit("Panel B equals Panel A at (0.40, 0.40) -- the cost "
                         "mapping is not applied")

    print(f"case={case}  time_summary={TIME_SUMMARY}")
    print(f"  nodes: {len(art_nodes)} x {len(prep_nodes)}  "
          f"Delta in [{art_nodes[0]:.4f}, {art_nodes[-1]:.4f}]")
    print(f"  A(ART .40, PrEP .10) = {a1:.6f} ; A(ART .10, PrEP .40) = {a2:.6f}")
    print(f"  at (0.40, 0.40): A = {a3:.6f}  B = {b1:.6f}  B/A = {b1 / a3:.6f}")
    return out, art_nodes, prep_nodes, case


SURF, ART_NODES, PREP_NODES, CASE_USED = load_surface()
X, Y = np.meshgrid(PREP_NODES, ART_NODES)      # X = PrEP, Y = ART (rows = ART)
VMIN = 1.0
VMAX = max(SURF[PANEL_A]["mean"].max(), SURF[PANEL_B]["mean"].max())
NORM = plt.Normalize(vmin=VMIN, vmax=VMAX)
# ---------------------------------------------------------------------------
# PLOT HELPERS
# ---------------------------------------------------------------------------
def slice_index(level):
    return int(np.argmin(np.abs(ART_NODES - level)))


def draw_heatmap(ax, g, title, xlabel, ylabel):
    """Z = g['mean'] with Z[i, j] at (ART_NODES[i], PREP_NODES[j]) -- pivoted."""
    # Explicit, IDENTICAL fill levels for every panel, so the shared colour scale
    # is legible from the colourbars (integer `levels=` lets matplotlib derive
    # different ranges per panel, which reads as two different scales).
    levels = np.linspace(VMIN, VMAX, 21)
    im = ax.contourf(X, Y, g["mean"], levels=levels, cmap=CMAP, norm=NORM)
    cs = ax.contour(X, Y, g["mean"], levels=np.arange(1.25, VMAX, 0.25),
                    colors="black", linewidths=0.4, alpha=0.35)
    ax.clabel(cs, inline=True, fontsize=7, fmt="%.1f", colors="black")
    for lv in POLICY_LEVELS:
        ax.axhline(lv, color="white", linewidth=0.8, linestyle=":", alpha=0.9, zorder=3)
        ax.axvline(lv, color="white", linewidth=0.8, linestyle=":", alpha=0.9, zorder=3)
        ax.plot([lv], [0.0], marker="o", color="black", ms=5, zorder=6)
        ax.plot([0.0], [lv], marker="o", color="black", ms=5, zorder=6)
    ax.set_xlabel(xlabel, fontsize=9, fontweight="bold")
    ax.set_ylabel(ylabel, fontsize=9, fontweight="bold")
    ax.set_title(title, fontsize=10, fontweight="bold", loc="left")
    ax.set_xlim(-0.01, 0.76)
    ax.set_ylim(-0.01, 0.76)
    ax.tick_params(labelsize=8)
    return im


def draw_slices(ax, g, title, xlabel, ylabel, legend_title):
    n = max(len(SLICE_LEVELS) - 1, 1)
    for k, lv in enumerate(SLICE_LEVELS):
        i = slice_index(lv)
        col = SLICE_CMAP(k / n)
        ax.plot(PREP_NODES, g["mean"][i, :], "-", linewidth=1.8, color=col,
                label=f"$\\Delta_{{ART}}$ = {lv:.2f}")
        ax.fill_between(PREP_NODES, g["lower"][i, :], g["upper"][i, :],
                        alpha=0.12, color=col, linewidth=0)
    ax.axhline(1.0, color="black", linestyle="--", linewidth=0.8, alpha=0.5)
    ax.set_xlabel(xlabel, fontsize=9, fontweight="bold")
    ax.set_ylabel(ylabel, fontsize=9, fontweight="bold")
    ax.set_title(title, fontsize=10, fontweight="bold", loc="left")
    ax.set_xlim(-0.01, 0.76)
    ax.grid(True, alpha=0.25)
    ax.tick_params(labelsize=8)
    leg = ax.legend(fontsize=8, ncol=min(len(SLICE_LEVELS), 3), loc="upper left",
                    title=legend_title, framealpha=0.9)
    leg.get_title().set_fontsize(8)
    return leg


def add_colorbar(fig, im, ax):
    cbar = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.046)
    cbar.set_label("IRR vs. status quo", fontsize=8, labelpad=4)
    cbar.ax.tick_params(labelsize=8)
    cbar.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))


def save(fig, name):
    out = OUT_DIR / name
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  wrote {out.name}  ({out.stat().st_size / 1024:.0f} KB)")


# ---------------------------------------------------------------------------
# FIGURE 1 -- Panel A: intervention USE reduction
# ---------------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5),
                               gridspec_kw={"width_ratios": [1, 1.05], "wspace": 0.42})

im = draw_heatmap(ax1, SURF[PANEL_A],
                  "(A) Intervention use reduction",
                  "$\\Delta_{PrEP}$ (proportional reduction)",
                  "$\\Delta_{ART}$ (proportional reduction)")
add_colorbar(fig, im, ax1)

draw_slices(ax2, SURF[PANEL_A],
            "(B) Slices at fixed $\\Delta_{ART}$",
            "$\\Delta_{PrEP}$ (proportional reduction)",
            "Incidence rate ratio vs. status quo",
            "ART use reduction")
fig.suptitle("Incidence rate ratio surface over intervention use reductions "
             f"($\\Delta_j \\in [0, 0.75]$; horiz.-mean IRR, CRN; case: {CASE_USED})",
             fontsize=9, y=1.02)
save(fig, "fig_response_surface_intervention.png")

# ---------------------------------------------------------------------------
# FIGURE 2 -- Panel B: GOVERNMENT FUNDING reduction (cost model applied)
# ---------------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5),
                               gridspec_kw={"width_ratios": [1, 1.05], "wspace": 0.42})

im = draw_heatmap(ax1, SURF[PANEL_B],
                  "(A) Government funding reduction",
                  "$\\Delta_{PrEP}$ (proportional funding reduction)",
                  "$\\Delta_{ART}$ (proportional funding reduction)")
add_colorbar(fig, im, ax1)

draw_slices(ax2, SURF[PANEL_B],
            "(B) Slices at fixed $\\Delta_{ART}$",
            "$\\Delta_{PrEP}$ (proportional funding reduction)",
            "Incidence rate ratio vs. status quo",
            "ART funding reduction")
fig.suptitle("Incidence rate ratio surface over government funding reductions "
             "(coverage mapped by the cost model; horiz.-mean IRR, CRN; "
             f"case: {CASE_USED})",
             fontsize=9, y=1.02)
save(fig, "fig_response_surface_funding.png")
# ---------------------------------------------------------------------------
# FIGURE 3 -- A vs B side by side (the insulation buffer)
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(11, 9))

imA = draw_heatmap(axes[0, 0], SURF[PANEL_A],
                   "(A) Intervention use reduction",
                   "$\\Delta_{PrEP}$ (proportional reduction)",
                   "$\\Delta_{ART}$ (proportional reduction)")
add_colorbar(fig, imA, axes[0, 0])
imB = draw_heatmap(axes[0, 1], SURF[PANEL_B],
                   "(B) Government funding reduction",
                   "$\\Delta_{PrEP}$ (proportional funding reduction)",
                   "$\\Delta_{ART}$ (proportional funding reduction)")
add_colorbar(fig, imB, axes[0, 1])

# Bottom row: the same cut under the two interpretations. Shared colour scale is
# intentional so the two surfaces can be read against each other; the slices
# therefore share a y-scale as well.
ymax = max(SURF[PANEL_A]["mean"].max(), SURF[PANEL_B]["mean"].max()) * 1.02

# (C) ART-only cut: PrEP held at 0, vary Delta_ART
i0, j0 = 0, 0
ax = axes[1, 0]
for g, lbl, col in ((SURF[PANEL_A], "Use reduction (Panel A)", "#597cbe"),
                    (SURF[PANEL_B], "Funding reduction (Panel B)", "#45aF84")):
    ax.plot(ART_NODES, g["mean"][:, j0], "-", linewidth=2.0, color=col, label=lbl)
    ax.fill_between(ART_NODES, g["lower"][:, j0], g["upper"][:, j0],
                    alpha=0.12, color=col, linewidth=0)
ax.axhline(1.0, color="black", linestyle="--", linewidth=0.8, alpha=0.5)
ax.set_xlabel("$\\Delta_{ART}$ (0 = no cut)", fontsize=9, fontweight="bold")
ax.set_ylabel("Incidence rate ratio", fontsize=9, fontweight="bold")
ax.set_title("(C) ART-only cut ($\\Delta_{PrEP} = 0$)", fontsize=10,
             fontweight="bold", loc="left")
ax.set_xlim(-0.01, 0.76); ax.set_ylim(0.99, ymax)
ax.grid(True, alpha=0.25); ax.tick_params(labelsize=8)
ax.legend(fontsize=8, loc="upper left", framealpha=0.9)

# (D) PrEP-only cut: ART held at 0, vary Delta_PrEP
ax = axes[1, 1]
for g, lbl, col in ((SURF[PANEL_A], "Use reduction (Panel A)", "#597cbe"),
                    (SURF[PANEL_B], "Funding reduction (Panel B)", "#45aF84")):
    ax.plot(PREP_NODES, g["mean"][i0, :], "-", linewidth=2.0, color=col, label=lbl)
    ax.fill_between(PREP_NODES, g["lower"][i0, :], g["upper"][i0, :],
                    alpha=0.12, color=col, linewidth=0)
ax.axhline(1.0, color="black", linestyle="--", linewidth=0.8, alpha=0.5)
ax.set_xlabel("$\\Delta_{PrEP}$ (0 = no cut)", fontsize=9, fontweight="bold")
ax.set_ylabel("Incidence rate ratio", fontsize=9, fontweight="bold")
ax.set_title("(D) PrEP-only cut ($\\Delta_{ART} = 0$)", fontsize=10,
             fontweight="bold", loc="left")
ax.set_xlim(-0.01, 0.76); ax.set_ylim(0.99, ymax)
ax.grid(True, alpha=0.25); ax.tick_params(labelsize=8)
ax.legend(fontsize=8, loc="upper left", framealpha=0.9)

fig.suptitle("Same cut, two interpretations: coverage mapped through the cost "
             f"model buffers the IRR (horiz.-mean IRR, CRN, shared colour scale; "
             f"case: {CASE_USED})", fontsize=9, y=0.995)
fig.tight_layout(rect=(0, 0, 1, 0.99))
save(fig, "fig_response_surfaces_comparison.png")

# ---------------------------------------------------------------------------
# CLOSING SELF-CHECK: every figure must have gone to the delivery folder
# ---------------------------------------------------------------------------
expected = ["fig_response_surface_intervention.png",
            "fig_response_surface_funding.png",
            "fig_response_surfaces_comparison.png"]
for name in expected:
    f = OUT_DIR / name
    if not f.exists() or f.stat().st_size < 20_000:
        raise SystemExit(f"expected output missing or suspiciously small: {f}")
print("all three figures are in the delivery folder:", OUT_DIR)
print("done.")



