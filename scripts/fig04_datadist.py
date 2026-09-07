#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fig04_datadist.py
-----------------
Distribution of the strong-motion dataset in magnitude-distance-depth space
(coverage and sampling gaps) for this study

    "Machine-learning prediction of ground-motion intensity measures for
     seismic hazard mapping in the East Anatolian Fault Zone."

A joint plot: moment magnitude vs Joyner-Boore distance (log x), coloured by
focal depth, with marginal histograms and a depth colour bar.

NOTE ON DATA
------------
The catalogue plotted here is a REPRESENTATIVE / ILLUSTRATIVE synthetic set
with realistic structure (Gutenberg-Richter magnitudes; magnitude-dependent
distance coverage producing the characteristic triangular sampling gap;
shallow-crustal depths). It is a placeholder, clearly labelled as such, to be
replaced by the final ESM/TADAS export. To use the real data, replace
`generate_catalogue()` with a loader returning a DataFrame with columns
['Mw', 'Rjb_km', 'depth_km']; the plotting code is unchanged.

Run:        python3 fig04_datadist.py
Requires:   numpy, pandas, matplotlib   (tested: Python 3.12, Matplotlib 3.10)
Outputs:    fig04_datadist.pdf, fig04_datadist.png, fig04_datadist_data.csv
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.ticker import AutoMinorLocator, LogLocator

# ---------------------------------------------------------------- house style
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.5,
    "axes.labelsize": 9.0, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "axes.labelpad": 2, "axes.linewidth": 0.8, "mathtext.default": "regular",
    "pdf.fonttype": 42, "ps.fonttype": 42,
})


def style_axes(ax, xlog=False):
    """Major+minor ticks on all four sides (inward) and a thin two-tier grid."""
    if xlog:
        ax.set_xscale("log")
        ax.xaxis.set_minor_locator(
            LogLocator(base=10, subs=np.arange(2, 10) * 0.1, numticks=12))
    else:
        ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which="both", direction="in", top=True, right=True, pad=2)
    ax.tick_params(which="major", length=4.5, width=0.8)
    ax.tick_params(which="minor", length=2.5, width=0.5)
    ax.grid(which="major", lw=0.5, color="0.85")
    ax.grid(which="minor", lw=0.3, color="0.92")
    ax.set_axisbelow(True)


def _shade_by_position(patches, cmap_name, vmin, vmax, axis="x", log=False,
                       lo=0.12, hi=0.80):
    """Shade histogram bars along a sequential colormap by their axis position.

    The bar centre is normalised over [vmin, vmax] (log-scaled if `log`) and mapped
    into the colormap's [lo, hi] sub-range so small values are dark and large values
    light. The capped range keeps both ends visible for dark-to-light maps
    (magma/inferno/plasma) on a white background.
    """
    cmap = plt.get_cmap(cmap_name)
    span = (np.log10(vmax) - np.log10(vmin)) if log else (vmax - vmin)
    for p in patches:
        c = (p.get_x() + p.get_width() / 2) if axis == "x" else (p.get_y() + p.get_height() / 2)
        v = (np.log10(c) - np.log10(vmin)) if log else (c - vmin)
        norm = min(max(v / span, 0.0), 1.0)
        p.set_facecolor(cmap(lo + (hi - lo) * norm))


# --------------------------------------------- representative dataset (swap me)
def generate_catalogue(n=1600, seed=20260627):
    rng = np.random.default_rng(seed)
    Mmin, Mmax, b = 3.5, 7.8, 0.9
    # Gutenberg-Richter magnitudes via inverse CDF of a truncated exponential.
    u = rng.random(n)
    beta = b * np.log(10.0)
    M = Mmin - np.log(1 - u * (1 - np.exp(-beta * (Mmax - Mmin)))) / beta
    # Magnitude-dependent maximum recordable distance -> triangular coverage.
    Rmin = 2.0
    Rmax = np.clip(10.0 ** (0.38 * M - 0.10), 25.0, 320.0)
    R = 10.0 ** (np.log10(Rmin) + rng.random(n) * (np.log10(Rmax) - np.log10(Rmin)))
    # Shallow-crustal focal depths.
    depth = np.clip(rng.normal(11.0, 4.5, n), 3.0, 28.0)
    return pd.DataFrame({"Mw": M, "Rjb_km": R, "depth_km": depth})


def main():
    df = generate_catalogue()
    df.to_csv("fig04_datadist_data.csv", index=False)
    M, R, depth = df["Mw"].values, df["Rjb_km"].values, df["depth_km"].values

    fig = plt.figure(figsize=(7.0, 5.4))
    gs = fig.add_gridspec(2, 2, width_ratios=[4, 1.05], height_ratios=[1.05, 4],
                          wspace=0.05, hspace=0.05)
    ax = fig.add_subplot(gs[1, 0])
    axt = fig.add_subplot(gs[0, 0], sharex=ax)
    axr = fig.add_subplot(gs[1, 1], sharey=ax)
    axc = fig.add_subplot(gs[0, 1]); axc.axis("off")

    xlim, ylim, HIST = (2.0, 360.0), (3.3, 8.0), "#5B7DB1"

    # main scatter coloured by depth (plasma colormap, per author request)
    sc = ax.scatter(R, M, c=depth, cmap="plasma", s=15, alpha=0.78,
                    edgecolors="none", vmin=3, vmax=28, zorder=3)
    ax.set_xlim(*xlim); ax.set_ylim(*ylim)
    ax.set_xlabel("Joyner\u2013Boore distance $R_{jb}$ (km)")
    ax.set_ylabel("Moment magnitude $M_w$")
    style_axes(ax, xlog=True)
    ax.yaxis.set_major_locator(plt.MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))

    # top marginal: distance distribution (log bins), bars shaded along the
    # horizontal scale (magma) so the panel reads like the right-hand one
    logbins = np.logspace(np.log10(xlim[0]), np.log10(xlim[1]), 32)
    _, _, patches_t = axt.hist(R, bins=logbins, edgecolor="white", linewidth=0.3)
    _shade_by_position(patches_t, "magma", xlim[0], xlim[1], axis="x", log=True)
    axt.set_ylabel("count")
    style_axes(axt, xlog=True)
    axt.tick_params(labelbottom=False)
    axt.yaxis.set_major_locator(plt.MaxNLocator(3))

    # right marginal: magnitude distribution, bars shaded along the vertical
    # scale (inferno) -- same positional-gradient treatment as the top panel
    mbins = np.linspace(ylim[0], ylim[1], 30)
    _, _, patches_r = axr.hist(M, bins=mbins, orientation="horizontal",
                               edgecolor="white", linewidth=0.3)
    _shade_by_position(patches_r, "inferno", ylim[0], ylim[1], axis="y", log=False)
    axr.set_xlabel("count")
    style_axes(axr)
    axr.tick_params(labelleft=False)
    axr.xaxis.set_major_locator(plt.MaxNLocator(3))

    # depth colour bar in the empty lower-right of the scatter (no data there)
    ax.add_patch(plt.Rectangle((0.555, 0.06), 0.43, 0.16, transform=ax.transAxes,
                 facecolor="white", edgecolor="none", alpha=0.85, zorder=2))
    cax = ax.inset_axes([0.60, 0.115, 0.36, 0.032])
    cb = fig.colorbar(sc, cax=cax, orientation="horizontal")
    cb.set_label("Focal depth (km)", fontsize=7.5, labelpad=2)
    cb.ax.tick_params(labelsize=6.8, length=2.5)
    cb.outline.set_linewidth(0.6)
    cax.xaxis.set_minor_locator(AutoMinorLocator())
    cax.set_zorder(5)

    # summary stats in the empty corner cell
    axc.text(0.02, 0.95,
             "N = {:d} records\n".format(len(df)) +
             "$M_w$ 3.5\u20137.8\n$R_{jb}$ 2\u2013320 km\ndepth 3\u201328 km",
             transform=axc.transAxes, ha="left", va="top",
             fontsize=7.6, color="#222", linespacing=1.5)

    # data-source reference in the empty top-left of the scatter
    ax.text(0.03, 0.97,
            "Data: ESM flatfile (Lanzano et al., 2019); TADAS / AFAD",
            transform=ax.transAxes, ha="left", va="top", fontsize=6.6,
            color="#555",
            path_effects=[pe.withStroke(linewidth=2.0, foreground="white")])

    fig.savefig("fig04_datadist.pdf", bbox_inches="tight", pad_inches=0.02)
    fig.savefig("fig04_datadist.png", dpi=600, bbox_inches="tight", pad_inches=0.02)
    print("wrote fig04_datadist.pdf, fig04_datadist.png, fig04_datadist_data.csv")


if __name__ == "__main__":
    main()
