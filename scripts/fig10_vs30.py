#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fig10_vs30.py  --  Site-condition (Vs30) model and resulting amplification.

Reproduces Figure 10 of this study: a two-panel map across the East
Anatolian Fault Zone --
    (a) Vs30 site model (m/s)        (b) site amplification, log10(A_site)
The amplification is derived directly from Vs30 through
    A_site = (Vref / Vs30) ** k ,  log10(A_site) = k*(log10(Vref) - log10(Vs30))
so the two panels are internally consistent: soft, low-Vs30 basins amplify
(log10 A > 0); stiff, high-Vs30 ranges de-amplify relative to the rock reference.

DATA HOOK
---------
VS30_FILE: a 3-column table  lon  lat  Vs30[m/s]  (None -> representative model
built from basin/range structure). The amplification panel is always recomputed
from whichever Vs30 field is used.

REQUIREMENTS:  GMT >= 6.5,  pygmt,  numpy        (pip install pygmt numpy)
USAGE:         python fig10_vs30.py
"""

import os
import numpy as np
import pygmt

# ----------------------------------------------------------------------------
# 0. Configuration
# ----------------------------------------------------------------------------
VS30_FILE = None                        # -> "vs30_model.xyz" for a real model
REGION    = [36.0, 41.0, 36.0, 40.0]
SPACING   = 0.02
VREF      = 760.0                        # reference-rock Vs30 (m/s)
KAMP      = 0.45                         # amplification exponent
OUTSTEM   = "fig10_vs30"
DPI       = 400

# ----------------------------------------------------------------------------
# 1. Faults and provincial centres
# ----------------------------------------------------------------------------
N_STRAND = [(40.58, 39.30), (40.40, 39.10), (40.15, 38.92), (39.92, 38.78),
            (39.60, 38.62), (39.35, 38.45)]
S_STRAND = [(39.35, 38.45), (39.00, 38.25), (38.62, 38.06), (38.20, 37.88),
            (37.78, 37.78), (37.55, 37.72), (37.30, 37.52), (37.05, 37.42),
            (36.80, 37.20), (36.55, 36.92), (36.35, 36.62), (36.22, 36.40)]
CARDAK   = [(38.00, 38.02), (37.70, 38.08), (37.40, 38.12), (37.10, 38.16),
            (36.85, 38.20)]
DST      = [(36.55, 36.92), (36.42, 36.55), (36.30, 36.22), (36.20, 36.02)]
FAULTS   = [N_STRAND + S_STRAND[1:], CARDAK, DST]
CITIES   = [("Gaziantep", 37.38, 37.07), ("Kahramanmaras", 36.92, 37.58),
            ("Malatya", 38.31, 38.36), ("Elazig", 39.22, 38.68),
            ("Antakya", 36.16, 36.20), ("Diyarbakir", 40.23, 37.91),
            ("Adiyaman", 38.28, 37.76)]

# ----------------------------------------------------------------------------
# 2. Vs30 field (representative model OR real file) and derived amplification
# ----------------------------------------------------------------------------
_lon = np.round(np.arange(REGION[0], REGION[1] + 1e-9, SPACING), 3)
_lat = np.round(np.arange(REGION[2], REGION[3] + 1e-9, SPACING), 3)
_LON, _LAT = np.meshgrid(_lon, _lat)

def representative_vs30():
    """Slope/structure-style proxy: smooth background, sedimentary basins lower
    Vs30, mountain ranges raise it. Returns lon, lat, Vs30 (1-D)."""
    logV = np.full(_LON.shape, np.log10(480.0))
    rng = np.random.default_rng(11)
    coarse = rng.normal(size=(14, 18))
    cy, cx = np.linspace(36, 40, 14), np.linspace(36, 41, 18)
    tmp = np.vstack([np.interp(_lon, cx, coarse[i]) for i in range(14)])
    logV += 0.10 * np.vstack([np.interp(_lat, cy, tmp[:, j])
                              for j in range(len(_lon))]).T
    g = lambda lo, la, s: np.exp(-(((_LON - lo) ** 2 + (_LAT - la) ** 2) / (2 * s ** 2)))
    basins = [(36.25, 36.35, 0.30, 0.42), (36.95, 37.55, 0.26, 0.30),
              (37.40, 37.10, 0.28, 0.22), (38.35, 38.35, 0.26, 0.28),
              (39.25, 38.62, 0.24, 0.25), (40.20, 37.90, 0.34, 0.22),
              (38.30, 37.75, 0.22, 0.20), (38.00, 37.40, 0.24, 0.24),
              (38.55, 37.65, 0.22, 0.20), (38.85, 37.95, 0.22, 0.18),
              (37.30, 37.50, 0.20, 0.18)]
    ranges = [(36.30, 36.90, 0.34, 0.18), (36.50, 38.20, 0.40, 0.20),
              (40.50, 39.20, 0.45, 0.22), (39.60, 39.30, 0.40, 0.20),
              (37.80, 38.65, 0.34, 0.15), (40.55, 37.55, 0.40, 0.12)]
    for lo, la, s, d in basins: logV -= d * g(lo, la, s)
    for lo, la, s, d in ranges: logV += d * g(lo, la, s)
    return _LON.ravel(), _LAT.ravel(), np.clip(10 ** logV, 150, 900).ravel()

if VS30_FILE and os.path.exists(VS30_FILE):
    d = np.loadtxt(VS30_FILE); lon, lat, vs = d[:, 0], d[:, 1], d[:, 2]
    print(f"loaded {len(vs)} Vs30 values from {VS30_FILE}")
else:
    lon, lat, vs = representative_vs30()
    print("VS30_FILE not set/found -> representative Vs30 model")

loga = np.clip(KAMP * (np.log10(VREF) - np.log10(vs)), -0.32, 0.32)  # log10 amplification

vs_grid  = pygmt.surface(x=lon, y=lat, z=vs,   region=REGION, spacing=SPACING, tension=0.25)
amp_grid = pygmt.surface(x=lon, y=lat, z=loga, region=REGION, spacing=SPACING, tension=0.25)

# ----------------------------------------------------------------------------
# 3. Render (two panels, shared layout/style with figs 8-9)
# ----------------------------------------------------------------------------
pygmt.config(FONT_ANNOT_PRIMARY="9p,Helvetica", FONT_LABEL="10p,Helvetica",
             FONT_TITLE="13p,Helvetica-Bold", MAP_FRAME_TYPE="plain",
             MAP_FRAME_PEN="0.9p,black", MAP_TITLE_OFFSET="4p",
             MAP_GRID_PEN_PRIMARY="0.25p,gray80", MAP_TICK_LENGTH_PRIMARY="3p",
             PS_CHAR_ENCODING="ISOLatin1+")

fig = pygmt.Figure()
with fig.subplot(nrows=1, ncols=2, subsize=("8.0c", "8.1c"),
                 margins=["0.8c", "2.0c"], sharey=True,
                 frame=["WSne", "xa1f0.5", "ya1f0.5"]):

    # ---- (a) Vs30 site model ----
    with fig.set_panel(0):
        pygmt.makecpt(cmap="roma", series=[150, 850, 25])
        fig.grdimage(grid=vs_grid, projection="M8.0c", region=REGION,
                     cmap=True, frame="+t(a) V@-s30@- site model")
        fig.grdcontour(grid=vs_grid, levels=100, pen="0.2p,gray25")
        fig.coast(water="173/205/226", shorelines="0.4p,gray40",
                  borders="1/0.9p,magenta", area_thresh=40, resolution="i")
        for seg in FAULTS:
            x, y = zip(*seg); fig.plot(x=x, y=y, pen="1.8p,white@15")
        for seg in FAULTS:
            x, y = zip(*seg); fig.plot(x=x, y=y, pen="0.8p,gray10")
        fig.plot(x=[c[1] for c in CITIES], y=[c[2] for c in CITIES],
                 style="c0.12c", fill="black", pen="0.3p,white")
        for name, lo, la in CITIES:
            fig.text(x=lo, y=la, text=name, font="6.5p,Helvetica,black",
                     justify="LB", offset="0.12c/0.09c", fill="white@25",
                     pen="0.25p,gray40", clearance="1.5p/1.5p+tO")
        fig.colorbar(position="JBC+o0/0.75c+w6.6c/0.28c+h",
                     frame=r'xa200f50+l"V@-s30@-  (m s@+-1@+)"')

    # ---- (b) site amplification ----
    with fig.set_panel(1):
        pygmt.makecpt(cmap="vik", series=[-0.32, 0.32, 0.02])
        fig.grdimage(grid=amp_grid, projection="M8.0c", region=REGION,
                     cmap=True, frame="+t(b) Site amplification")
        fig.grdcontour(grid=amp_grid, levels=0.05, pen="0.2p,gray25")
        fig.coast(water="173/205/226", shorelines="0.4p,gray40",
                  borders="1/0.9p,magenta", area_thresh=40, resolution="i")
        for seg in FAULTS:
            x, y = zip(*seg); fig.plot(x=x, y=y, pen="1.8p,white@15")
        for seg in FAULTS:
            x, y = zip(*seg); fig.plot(x=x, y=y, pen="0.8p,gray10")
        fig.plot(x=[c[1] for c in CITIES], y=[c[2] for c in CITIES],
                 style="c0.12c", fill="black", pen="0.3p,white")
        fig.colorbar(position="JBC+o0/0.75c+w6.6c/0.28c+h",
                     frame=r'xa0.1f0.05+l"log@-10@-(A@-site@-)  (amplification)"')

fig.savefig(f"{OUTSTEM}.png", dpi=DPI)
fig.savefig(f"{OUTSTEM}.pdf")
print(f"wrote {OUTSTEM}.png ({DPI} dpi) and {OUTSTEM}.pdf")
