#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fig05_workflow.py
-----------------
Methodological workflow figure for this study

    "Machine-learning prediction of ground-motion intensity measures for
     seismic hazard mapping in the East Anatolian Fault Zone."

Renders a five-stage pipeline schematic
    Open-access data -> Feature engineering -> ML modelling
    -> Hazard mapping -> Validation
as a colour-blind- and grayscale-safe diagram, and writes a vector PDF
(for LaTeX) plus a 600 dpi PNG preview, both named ``fig05_workflow``.

Run:
    python3 fig05_workflow.py
Requires:
    matplotlib >= 3.4   (tested with 3.10; Python 3.12)
No data files are needed -- the figure is a pure schematic.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.path import Path

# One consistent font family across all paper figures; embed fonts in PDF.
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

# ----------------------------------------------------------------------
# Stage content -- mirrors the Methods section exactly.
#   role  : pipeline role shown above the box (INPUT / PROCESSING / ...)
#   title : stage name (header strip)
#   tool  : software that drives the stage (header sub-label)
#   hdr   : header colour;  body : body tint (grayscale-distinct values)
#   items : bulleted operations (a leading run of spaces = wrapped line)
#   tag   : link to the Methods listings/algorithms (monospace tag below)
# ----------------------------------------------------------------------
stages = [
    dict(role="INPUT", title="Open-access data", tool="ESM / IEB / SRTM",
         hdr="#3B6EA5", body="#DCE6F2",
         items=["Strong-motion flatfile", "Earthquake catalogue",
                "Site model (Vs30, f1)", "SRTM / GEBCO terrain"], tag=""),
    dict(role="PROCESSING", title="Feature engineering", tool="Python / pandas",
         hdr="#C77A30", body="#F4E4D1",
         items=["QC and merge records", "x = Mw, logRjb, depth,",
                "       logVs30, logf1", "y = log10(IM)", "Event-grouped split"],
         tag="Listing 1"),
    dict(role="PROCESSING", title="ML modelling", tool="XGBoost / sklearn",
         hdr="#4F8F46", body="#DCEBD7",
         items=["Gradient-boosted trees", "RF, SVR baselines",
                "Grouped 10-fold CV", "Hyper-parameter tuning"], tag="Listing 3"),
    dict(role="OUTPUT", title="Hazard mapping", tool="GMT / PyGMT",
         hdr="#7E5AA2", body="#E6DCF1",
         items=["Predict IM on grid", "Site amplification",
                "PGA, PGV, Sa maps", "Seismic-hazard map"],
         tag="Listing 2  -  Algorithm 2"),
    dict(role="VALIDATION", title="Validation", tool="Python / Matplotlib",
         hdr="#B0492F", body="#F3DBD2",
         items=["R2, RMSE, MAE", "Residual diagnostics",
                "vs empirical GMPEs", "Scatter reduction"], tag=""),
]

# ----------------------------------------------------------------------
# Geometry (all in axis-fraction coordinates on a [0,1] x [0,1] canvas).
# ----------------------------------------------------------------------
n = len(stages)
xl, xr, gap = 0.010, 0.010, 0.034          # outer margins and inter-box gap
bw = (1 - xl - xr - (n - 1) * gap) / n      # box width
y0, y1, hdr_h = 0.400, 0.880, 0.140         # box bottom / top, header height
yh = y1 - hdr_h                             # header strip bottom
role_y, tag_y = 0.930, 0.350                # role-label and tag rows


def main():
    fig = plt.figure(figsize=(7.6, 4.25))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    centres = []
    for i, s in enumerate(stages):
        x = xl + i * (bw + gap)
        cx = x + bw / 2
        centres.append((x, cx))

        # body box
        ax.add_patch(FancyBboxPatch(
            (x, y0), bw, y1 - y0,
            boxstyle="round,pad=0,rounding_size=0.014",
            linewidth=1.1, edgecolor=s["hdr"], facecolor=s["body"], zorder=2))
        # header strip (rounded top; squared bottom via a half-height overlay)
        ax.add_patch(FancyBboxPatch(
            (x, yh), bw, hdr_h,
            boxstyle="round,pad=0,rounding_size=0.014",
            linewidth=0, facecolor=s["hdr"], zorder=3))
        ax.add_patch(plt.Rectangle(
            (x, yh), bw, hdr_h * 0.5, linewidth=0, facecolor=s["hdr"], zorder=3))

        # role label above the box
        ax.text(cx, role_y, s["role"], ha="center", va="center",
                fontsize=6.4, color="#555555", weight="bold")
        # header title + tool sub-label
        ax.text(cx, yh + hdr_h * 0.62, s["title"], ha="center", va="center",
                fontsize=7.7, color="white", weight="bold", zorder=4)
        ax.text(cx, yh + hdr_h * 0.23, s["tool"], ha="center", va="center",
                fontsize=6.3, color="white", style="italic", zorder=4)

        # body items (bulleted; wrapped continuation lines start with spaces)
        top, dy = yh - 0.040, 0.060
        for j, it in enumerate(s["items"]):
            indented = it.startswith("       ")
            txt = ("   " + it.strip()) if indented else ("\u2022 " + it.strip())
            ax.text(x + 0.013, top - j * dy, txt, ha="left", va="center",
                    fontsize=7.0, color="#1f1f1f", zorder=4)

        # artifact tag below the box
        if s["tag"]:
            ax.text(cx, tag_y, s["tag"], ha="center", va="center", fontsize=6.4,
                    color="#444444", family="DejaVu Sans Mono",
                    bbox=dict(boxstyle="round,pad=0.25", fc="#f3f3f3",
                              ec="#bbbbbb", lw=0.6), zorder=4)

    # forward arrows between consecutive stages
    ymid = (y0 + y1) / 2
    for i in range(n - 1):
        x_i = centres[i][0] + bw
        x_j = centres[i + 1][0]
        ax.add_patch(FancyArrowPatch(
            (x_i + 0.004, ymid), (x_j - 0.004, ymid),
            arrowstyle="-|>", mutation_scale=15, lw=2.0, color="#3a3a3a", zorder=5))

    # feedback arrow: Validation -> ML modelling (cross-validated tuning loop)
    x_v, x_m, fb_y = centres[4][1], centres[2][1], 0.175
    verts = [(x_v, y0 - 0.005), (x_v, fb_y), (x_m, fb_y), (x_m, y0 - 0.005)]
    ax.add_patch(FancyArrowPatch(
        path=Path(verts, [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]),
        arrowstyle="-|>", mutation_scale=14, lw=1.6, color="#8a8a8a",
        linestyle=(0, (5, 3)), zorder=1))
    ax.text((x_v + x_m) / 2, fb_y - 0.004, "re-tune / iterate", ha="center",
            va="center", fontsize=6.8, color="#6f6f6f", style="italic",
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none"))

    # footer notes
    ax.text(0.010, 0.045, "Algorithm 1 summarises the full pipeline.",
            ha="left", va="center", fontsize=6.4, color="#777777", style="italic")
    ax.text(0.990, 0.045, "open-access data \u2192 reproducible hazard maps",
            ha="right", va="center", fontsize=6.4, color="#777777", style="italic")

    fig.savefig("fig05_workflow.pdf", bbox_inches="tight", pad_inches=0.04)
    fig.savefig("fig05_workflow.png", dpi=600, bbox_inches="tight", pad_inches=0.04)
    print("wrote fig05_workflow.pdf and fig05_workflow.png")


if __name__ == "__main__":
    main()
