"""
Greyscale variant of the "Mean Incidence Risk Ratio by Funding Reduction" plot.

Distinguishes the three scenarios WITHOUT color, using marker shape plus three
evenly-stepped greyscale fills: circle (light grey), square (mid grey),
triangle (near-black). Each marker keeps a black outline so its shape stays
readable even where the fill is close in value to the background band.
"""

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

# Notes:
# colors for previous plots that are not greyscale
colors = {
    "Vary PrEP (ART = 0)": "#A50D3F",   # deep rose/crimson
    "Vary ART (PrEP = 0)": "#2A9CED",   # bright sky blue
    "Vary both equally":   "#675496",   # violet = avg(rose, blue)
}


# ---- Data: reduction (%) -> scenario -> (mean, ci_low, ci_high) ----
data = {
    10: {
        "Vary PrEP (ART = 0)":   (1.05, 1.01, 1.10),
        "Vary ART (PrEP = 0)":   (1.19, 1.13, 1.26),
        "Vary both equally":     (1.24, 1.19, 1.30),
    },
    20: {
        "Vary PrEP (ART = 0)":   (1.09, 1.03, 1.16),
        "Vary ART (PrEP = 0)":   (1.39, 1.30, 1.48),
        "Vary both equally":     (1.51, 1.40, 1.63),
    },
    30: {
        "Vary PrEP (ART = 0)":   (1.14, 1.06, 1.23),
        "Vary ART (PrEP = 0)":   (1.59, 1.47, 1.72),
        "Vary both equally":     (1.78, 1.65, 1.92),
    },
    40: {
        "Vary PrEP (ART = 0)":   (1.18, 1.08, 1.29),
        "Vary ART (PrEP = 0)":   (1.77, 1.62, 1.94),
        "Vary both equally":     (2.05, 1.87, 2.24),
    },
}

scenarios = ["Vary PrEP (ART = 0)", "Vary ART (PrEP = 0)", "Vary both equally"]

# Distinguish scenarios by marker shape + three evenly-stepped greyscale fills
style = {
    "Vary PrEP (ART = 0)": dict(marker="o", facecolor="#D9D9D9", linestyle="-"),  # light grey
    "Vary ART (PrEP = 0)": dict(marker="s", facecolor="#737373", linestyle="-"),  # mid grey
    "Vary both equally":   dict(marker="^", facecolor="#0D0D0D", linestyle="-"),  # near-black
}
EDGE = "black"

reductions = sorted(data.keys())
n_scen = len(scenarios)
offset_max = 0.16
band_padding = 0.12
band_gap = 0.08
row_gap = 2 * (offset_max + band_padding) + band_gap
offsets = np.linspace(-offset_max, offset_max, n_scen)

fig, ax = plt.subplots(figsize=(10, 6.3))

band_colors = ["#ececec", "#dcdcdc", "#cacaca", "#b8b8b8"]
for i, red in enumerate(reductions):
    y_center = i * row_gap
    ax.axhspan(y_center - row_gap / 2 + band_gap / 2, y_center + row_gap / 2 - band_gap / 2,
               color=band_colors[i % len(band_colors)], zorder=0)

for i, red in enumerate(reductions):
    y_center = i * row_gap
    for j, scen in enumerate(scenarios):
        mean, lo, hi = data[red][scen]
        y = y_center + offsets[j]
        s = style[scen]

        # Draw whiskers separately so we can set a distinct linestyle per scenario
        _, caps, bars = ax.errorbar(
            mean, y,
            xerr=[[mean - lo], [hi - mean]],
            fmt="none", ecolor="black",
            elinewidth=1.6, capsize=4, zorder=3,
        )
        for b in bars:
            b.set_linestyle(s["linestyle"])

        ax.plot(mean, y, marker=s["marker"], markersize=9,
                 markerfacecolor=s["facecolor"], markeredgecolor=EDGE,
                 markeredgewidth=1.3, linestyle="none", zorder=4)

        ax.text(hi * 1.02, y, f"{mean:.2f}", va="center", ha="left",
                 fontsize=10, fontweight="bold", color="black", zorder=5)

ax.set_xscale("log")
ax.set_xlim(0.97, 2.4)
ax.set_xticks([1.0, 1.2, 1.5, 2.0])
ax.set_xticklabels(["1.0", "1.2", "1.5", "2.0"])
ax.xaxis.set_minor_locator(mticker.NullLocator())

ax.xaxis.grid(True, color="white", linewidth=1.2, linestyle="--", zorder=1)
ax.set_axisbelow(False)

ax.set_yticks([])
ax.set_ylim(-row_gap / 2, (len(reductions) - 1) * row_gap + row_gap / 2)

label_transform = ax.get_yaxis_transform()
for i, red in enumerate(reductions):
    y_center = i * row_gap
    ax.text(-0.01, y_center, f"{red}%", transform=label_transform,
             ha="right", va="center", fontsize=12, fontweight="bold", color="#333333")

ax.set_xlabel("Mean incidence risk ratio (scenario vs. baseline, natural-log scale)")
ax.text(-0.09, 0.5, "Funding\nReduction", transform=ax.transAxes, rotation=90,
        ha="center", va="center", fontsize=11, color="dimgray")

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)

fig.text(0.06, 0.97, "Mean Incidence Risk Ratio by Funding Reduction",
          fontsize=16, fontweight="bold", ha="left", va="top")
fig.text(0.06, 0.92, "Scenario vs. no-reduction baseline; points are means, whiskers are 90% CIs",
          fontsize=11, color="dimgray", ha="left", va="top")

# Legend: shape + fill + linestyle combo per scenario
# handles = [plt.Line2D([0], [0], marker="o", color=colors[s], linestyle="-",
#                        markersize=9, linewidth=2.5) for s in scenarios]
handles = [
    plt.Line2D([0], [0], marker=style[s]["marker"], color="black",
                markerfacecolor=style[s]["facecolor"], markeredgecolor=EDGE,
                markeredgewidth=1.3, linestyle=style[s]["linestyle"],
                markersize=9, linewidth=1.8)
    for s in scenarios
]
legend = fig.legend(handles, scenarios, title="Scenario", loc="lower center",
                     ncol=3, bbox_to_anchor=(0.5, 0.01), frameon=True,
                     fontsize=10, title_fontsize=11)
legend.get_frame().set_edgecolor("#cccccc")
legend.get_frame().set_facecolor("white")

plt.tight_layout(rect=[0, 0.1, 1, 0.87])
plt.savefig("incidence_risk_ratio_plot_greyscale.png", dpi=200)
plt.show()