#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fig06_mlarch -- ML architecture & feature set (gradient-boosted tree ensemble)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle

plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "ps.fonttype": 42,
                     "mathtext.fontset": "dejavusans"})

BLUE_H, BLUE_B = "#3B6EA5", "#DCE6F2"     # input features
GREEN_H, GREEN_B = "#4F8F46", "#E3EFDE"   # model (ensemble)
ORANGE_H, ORANGE_B = "#C77A30", "#F4E4D1" # targets
NODE = "#2f6b2a"; LEAF = "#9ec48f"; LINE = "#9aa0a6"; ARR = "#3a3a3a"

fig = plt.figure(figsize=(8.0, 4.7))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

def rbox(x, y, w, h, fc, ec, lw=1.1, z=2, r=0.014):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
        boxstyle=f"round,pad=0,rounding_size={r}", linewidth=lw,
        edgecolor=ec, facecolor=fc, zorder=z))

def header(x, y, w, h, fc, title, z=3):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
        boxstyle="round,pad=0,rounding_size=0.014", linewidth=0, facecolor=fc, zorder=z))
    ax.add_patch(plt.Rectangle((x, y), w, h*0.5, linewidth=0, facecolor=fc, zorder=z))
    ax.text(x+w/2, y+h*0.5, title, ha="center", va="center",
            fontsize=7.6, color="white", weight="bold", zorder=z+1)

# ---------------------------------------------------------------- INPUT panel
ix, iw = 0.018, 0.170
iy0, iy1, ihdr = 0.150, 0.820, 0.082
rbox(ix, iy0, iw, iy1-iy0, BLUE_B, BLUE_H)
header(ix, iy1-ihdr, iw, ihdr, BLUE_H, "Input vector  x")
feats = [(r"$M_w$", "moment magnitude"),
         (r"$\log R_{jb}$", "distance (km)"),
         (r"depth", "focal depth (km)"),
         (r"$\log V_{s30}$", "site stiffness"),
         (r"$\log f_1$", "site frequency")]
btop, bbot = iy1-ihdr-0.035, iy0+0.030
fy = np.linspace(btop, bbot, len(feats))
fnode_h = 0.090
feat_y = []
for (sym, desc), yc in zip(feats, fy):
    rbox(ix+0.012, yc-fnode_h/2, iw-0.024, fnode_h, "white", BLUE_H, lw=0.9, z=4, r=0.012)
    ax.text(ix+0.026, yc+0.014, sym, ha="left", va="center", fontsize=8.2, color="#13314f", zorder=5)
    ax.text(ix+0.026, yc-0.022, desc, ha="left", va="center", fontsize=5.9, color="#444", zorder=5)
    feat_y.append(yc)
ax.text(ix+iw/2, 0.880, "FEATURES", ha="center", va="center", fontsize=6.6, color="#555", weight="bold")
ax.text(ix+iw/2, 0.095, "Listing 1", ha="center", va="center", fontsize=6.3, color="#444",
        family="DejaVu Sans Mono",
        bbox=dict(boxstyle="round,pad=0.25", fc="#f3f3f3", ec="#bbb", lw=0.6))

# ---------------------------------------------------------------- ENSEMBLE panel
ex, ew = 0.300, 0.430
ey0, ey1, ehdr = 0.150, 0.840, 0.082
rbox(ex, ey0, ew, ey1-ey0, GREEN_B, GREEN_H)
header(ex, ey1-ehdr, ew, ehdr, GREEN_H, "Gradient-boosted regression trees  (XGBoost)")
ax.text(ex+ew/2, 0.905, "MODEL", ha="center", va="center", fontsize=6.6, color="#555", weight="bold")

def draw_tree(cx, cy, w=0.052, h=0.150):
    root = (cx, cy+h/2)
    L, R = (cx-w/3, cy), (cx+w/3, cy)
    leaves = [(cx-w/2, cy-h/2), (cx-w/6, cy-h/2), (cx+w/6, cy-h/2), (cx+w/2, cy-h/2)]
    for a, b in [(root, L), (root, R), (L, leaves[0]), (L, leaves[1]),
                 (R, leaves[2]), (R, leaves[3])]:
        ax.plot([a[0], b[0]], [a[1], b[1]], color=LINE, lw=0.8, zorder=4)
    pts = [root, L, R]
    ax.scatter([p[0] for p in pts], [p[1] for p in pts], s=26, color=NODE,
               edgecolors="white", linewidths=0.5, zorder=5)
    ax.scatter([p[0] for p in leaves], [p[1] for p in leaves], s=20, marker="s",
               color=LEAF, edgecolors=NODE, linewidths=0.5, zorder=5)

cy_t = 0.520
xs = [ex+0.060, ex+0.160, ex+0.260, ex+0.372]
labels = ["Tree 1", "Tree 2", "Tree 3", "Tree M"]
for k, (xc, lab) in enumerate(zip(xs, labels)):
    draw_tree(xc, cy_t)
    ax.text(xc, cy_t-0.115, lab, ha="center", va="center", fontsize=6.2, color="#2f6b2a", weight="bold")
    if k < len(xs)-1:
        mid = (xc+xs[k+1])/2
        sym = r"$\cdots$" if k == 2 else "+"
        ax.text(mid, cy_t, sym, ha="center", va="center", fontsize=10 if sym=="+" else 12,
                color="#555", weight="bold")

# sequential residual-fitting arrow above the trees
ax.add_patch(FancyArrowPatch((xs[0]-0.01, cy_t+0.135), (xs[-1]+0.01, cy_t+0.135),
    arrowstyle="-|>", mutation_scale=12, lw=1.3, color="#8a8a8a",
    linestyle=(0, (4, 2)), zorder=4))
ax.text((xs[0]+xs[-1])/2, cy_t+0.170,
        r"sequential fit to pseudo-residuals  $r_m=-\partial\mathcal{L}/\partial F$",
        ha="center", va="center", fontsize=6.3, color="#6f6f6f", style="italic")

# additive model + regularisation notes below the trees
ax.text(ex+ew/2, 0.300, r"$F_M(\mathbf{x})=F_0+\eta\sum_{m=1}^{M}h_m(\mathbf{x})$",
        ha="center", va="center", fontsize=9.0, color="#13320f")
ax.text(ex+ew/2, 0.232, r"shrinkage $\eta$  $\cdot$  L2 penalty $\lambda$  $\cdot$  row/column subsampling",
        ha="center", va="center", fontsize=6.2, color="#444")
ax.text(ex+ew/2, 0.185, "baselines: random forest, SVR (identical protocol)",
        ha="center", va="center", fontsize=6.0, color="#666", style="italic")
ax.text(ex+ew-0.052, 0.726, "Listing 3", ha="center", va="center", fontsize=6.3, color="#444",
        family="DejaVu Sans Mono",
        bbox=dict(boxstyle="round,pad=0.25", fc="#eef4ea", ec="#bccfb5", lw=0.6), zorder=6)

# fan-in lines: each feature -> ensemble left edge
for yc in feat_y:
    ax.add_patch(FancyArrowPatch((ix+iw, yc), (ex, np.clip(yc, ey0+0.06, ey1-ehdr-0.04)),
        arrowstyle="-", lw=0.8, color=LINE, zorder=1,
        connectionstyle="arc3,rad=0.05"))

# ---------------------------------------------------------------- SUM + OUTPUT
sx, sy, sr = 0.762, cy_t, 0.020
ax.add_patch(FancyArrowPatch((ex+ew, cy_t), (sx-sr-0.004, sy),
    arrowstyle="-|>", mutation_scale=14, lw=2.0, color=ARR, zorder=5))
ax.add_patch(Circle((sx, sy), sr, facecolor="white", edgecolor=GREEN_H, lw=1.3, zorder=5))
ax.text(sx, sy, r"$\Sigma$", ha="center", va="center", fontsize=11, color="#13320f", zorder=6)
ax.text(sx, sy+0.075, r"$\hat{y}=\log_{10}(\mathrm{IM})$", ha="center", va="center",
        fontsize=7.2, color="#222")

ox, ow = 0.808, 0.176
oy0, oy1, ohdr = 0.150, 0.820, 0.082
rbox(ox, oy0, ow, oy1-oy0, ORANGE_B, ORANGE_H)
header(ox, oy1-ohdr, ow, ohdr, ORANGE_H, "Intensity measures")
ims = [(r"PGA", "peak accel."), (r"PGV", "peak velocity"),
       (r"$S_a(T)$", "spectral accel."), (r"$I_A$", "Arias intensity"),
       (r"$D_{5-95}$", "sig. duration")]
oy = np.linspace(oy1-ohdr-0.035-0.045, oy0+0.030+0.045, len(ims))
onode_h = 0.090
out_y = []
for (sym, desc), yc in zip(ims, oy):
    rbox(ox+0.012, yc-onode_h/2, ow-0.024, onode_h, "white", ORANGE_H, lw=0.9, z=4, r=0.012)
    ax.text(ox+0.026, yc+0.014, sym, ha="left", va="center", fontsize=8.0, color="#6e3c10", zorder=5)
    ax.text(ox+0.026, yc-0.022, desc, ha="left", va="center", fontsize=5.9, color="#444", zorder=5)
    out_y.append(yc)
ax.text(ox+ow/2, 0.880, "TARGETS", ha="center", va="center", fontsize=6.6, color="#555", weight="bold")
ax.text(ox+ow/2, 0.095, "one model per IM", ha="center", va="center", fontsize=6.1,
        color="#666", style="italic")

# fan-out: sum -> each IM node
for yc in out_y:
    ax.add_patch(FancyArrowPatch((sx+sr+0.004, sy), (ox, yc),
        arrowstyle="-|>", mutation_scale=9, lw=0.9, color=LINE, zorder=1,
        connectionstyle="arc3,rad=0.04"))

fig.savefig("fig06_mlarch.pdf", bbox_inches="tight", pad_inches=0.04)
fig.savefig("fig06_mlarch.png", dpi=600, bbox_inches="tight", pad_inches=0.04)
print("wrote fig06_mlarch.pdf and .png")
