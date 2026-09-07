#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fig09_sa.py  --  Predicted spectral acceleration (and PGV) maps at selected periods.

Reproduces Figure 9 of this study: a 2x2 panel of predicted ground-motion
intensity measures across the East Anatolian Fault Zone --
    (a) PGV            (b) Sa(0.3 s)
    (c) Sa(1.0 s)      (d) Sa(3.0 s)
Each panel carries its own colour bar (units/ranges differ), active-fault traces,
the two 2023 Kahramanmaras epicentres, the national border and PGA-style isolines.
Short-period Sa stays concentrated on the rupture; PGV and long-period Sa broaden.

DATA HOOK
---------
PRED_FILES maps each intensity measure to a 3-column table
        lon  lat  log10(IM)            # IM in g for Sa, cm/s for PGV
Entries left as None fall back to a representative attenuation field so the
figure renders out of the box. Drop in your trained-model grids to make it real.

REQUIREMENTS:  GMT >= 6.5,  pygmt,  numpy        (pip install pygmt numpy)
USAGE:         python fig09_sa.py
"""

import os
import numpy as np
import pygmt

# ----------------------------------------------------------------------------
# 0. Configuration
# ----------------------------------------------------------------------------
REGION     = [36.0, 41.0, 36.0, 40.0]
SPACING    = 0.02
OUTSTEM    = "fig09_sa"
DPI        = 400

# real-data hook: lon, lat, log10(IM) tables (None -> representative field)
PRED_FILES = {"pgv": None, "sa03": None, "sa10": None, "sa30": None}

# per-panel definition:
#   key, panel title, colour-bar label, main-rupture peak, b, h, zmin, zmax
PANELS = [
    ("pgv",  "(a) PGV",           r"log@-10@-(PGV),  cm s@+-1@+",
      2.10, 0.95,  7,  0.30,  2.15),
    ("sa03", "(b) S@-a@-(0.3 s)", r"log@-10@- S@-a@-(0.3 s),  g",
      0.40, 1.05,  6, -1.30,  0.45),
    ("sa10", "(c) S@-a@-(1.0 s)", r"log@-10@- S@-a@-(1.0 s),  g",
     -0.10, 0.92,  8, -1.70, -0.05),
    ("sa30", "(d) S@-a@-(3.0 s)", r"log@-10@- S@-a@-(3.0 s),  g",
     -0.60, 0.85, 10, -2.30, -0.55),
]

# ----------------------------------------------------------------------------
# 1. Active-fault geometry and 2023 epicentres
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
EVENTS = [   # lon, lat, label, star-size, label-justify, label-offset
    (37.22, 37.17, "2023 Pazarcik M7.8",        "0.42c", "TL",  "0.10c/-0.10c"),
    (37.20, 38.02, "2023 Elbistan M7.5",        "0.42c", "BC",  "0c/0.12c"),
    (39.07, 38.36, "2020 Elazig-Sivrice M6.7",  "0.32c", "TL",  "0.10c/-0.10c"),
    (40.42, 39.04, "2003 Bingol M6.4",          "0.32c", "RB", "-0.12c/0.06c"),
    (40.05, 38.80, "2010 Kovancilar M6.1",      "0.32c", "TC",  "0c/-0.11c"),
    (36.20, 36.20, "2023 Defne-Hatay M6.3",     "0.32c", "LB",  "0.10c/0.08c"),
]

# ----------------------------------------------------------------------------
# 2. Distance-to-fault helpers (computed once, reused by every IM)
# ----------------------------------------------------------------------------
LAT0 = 38.0
def _to_km(lon, lat):
    return ((np.asarray(lon) - 36.0) * 111.320 * np.cos(np.radians(LAT0)),
            (np.asarray(lat) - 36.0) * 110.574)

def _densify(poly, step=1.0):
    pts = []
    for (l0, a0), (l1, a1) in zip(poly[:-1], poly[1:]):
        x0, y0 = _to_km(l0, a0); x1, y1 = _to_km(l1, a1)
        n = max(2, int(np.hypot(x1 - x0, y1 - y0) / step) + 1)
        for t in np.linspace(0, 1, n):
            pts.append((x0 + (x1 - x0) * t, y0 + (y1 - y0) * t))
    return np.array(pts)

_lon = np.round(np.arange(REGION[0], REGION[1] + 1e-9, SPACING), 3)
_lat = np.round(np.arange(REGION[2], REGION[3] + 1e-9, SPACING), 3)
_LON, _LAT = np.meshgrid(_lon, _lat)
_GX, _GY = _to_km(_LON, _LAT)
_R = {}
for _name, _poly in {"main": S_STRAND + CARDAK, "north": N_STRAND, "dst": DST}.items():
    _r2 = np.full(_LON.shape, np.inf)
    for _px, _py in _densify(_poly):
        np.minimum(_r2, (_GX - _px) ** 2 + (_GY - _py) ** 2, out=_r2)
    _R[_name] = np.sqrt(_r2)
_rng = np.random.default_rng(7)
_coarse = _rng.normal(size=(11, 14))
_cy, _cx = np.linspace(36, 40, 11), np.linspace(36, 41, 14)
_tmp = np.vstack([np.interp(_lon, _cx, _coarse[i]) for i in range(11)])
_PERT = np.vstack([np.interp(_lat, _cy, _tmp[:, j]) for j in range(len(_lon))]).T

def representative_field(main_peak, b, h, zmin, zmax):
    """log10(IM) = max_sys[ peak_sys - b*log10(1 + R/h) ], enveloped over the
    fault systems with the 2023 rupture dominant."""
    peaks = {"main": main_peak, "north": main_peak - 0.45, "dst": main_peak - 0.30}
    f = np.full(_LON.shape, -9.9)
    for s in _R:
        np.maximum(f, peaks[s] - b * np.log10(1 + _R[s] / h), out=f)
    return _LON.ravel(), _LAT.ravel(), np.clip(f + 0.04 * _PERT, zmin, zmax).ravel()

def build_grid(key, main_peak, b, h, zmin, zmax):
    path = PRED_FILES.get(key)
    if path and os.path.exists(path):
        d = np.loadtxt(path); lon, lat, val = d[:, 0], d[:, 1], d[:, 2]
        print(f"{key}: loaded {len(val)} predictions from {path}")
    else:
        lon, lat, val = representative_field(main_peak, b, h, zmin, zmax)
        print(f"{key}: representative field")
    g = pygmt.surface(x=lon, y=lat, z=val, region=REGION, spacing=SPACING, tension=0.35)
    return pygmt.grdclip(g, below=[zmin, zmin], above=[zmax, zmax])

# ----------------------------------------------------------------------------
# 3. Render
# ----------------------------------------------------------------------------
pygmt.config(FONT_ANNOT_PRIMARY="9p,Helvetica", FONT_LABEL="9.5p,Helvetica",
             FONT_TITLE="12p,Helvetica-Bold", MAP_FRAME_TYPE="plain",
             MAP_FRAME_PEN="0.8p,black", MAP_TITLE_OFFSET="4p",
             MAP_GRID_PEN_PRIMARY="0.25p,gray80", MAP_TICK_LENGTH_PRIMARY="3p",
             PS_CHAR_ENCODING="ISOLatin1+")

fig = pygmt.Figure()
with fig.subplot(nrows=2, ncols=2, subsize=("7.0c", "7.1c"),
                 margins=["0.7c", "2.1c"], sharex=True, sharey=True,
                 frame=["WSne", "xa1f0.5", "ya1f0.5"]):
    for idx, (key, title, cblab, mp, b, h, zmin, zmax) in enumerate(PANELS):
        grid = build_grid(key, mp, b, h, zmin, zmax)
        with fig.set_panel(idx):
            pygmt.makecpt(cmap="lajolla", series=[zmin, zmax, 0.05], reverse=True)
            fig.grdimage(grid=grid, projection="M7.0c", region=REGION, cmap=True,
                         frame=["xa1f0.5", "ya1f0.5", f"+t{title}"])
            fig.grdcontour(grid=grid, levels=0.05, pen="0.15p,139/69/19")
            fig.coast(water="173/205/226", shorelines="0.3p,gray40",
                      borders="1/0.8p,magenta", area_thresh=40, resolution="i")
            for seg in FAULTS:
                x, y = zip(*seg); fig.plot(x=x, y=y, pen="1.6p,white@15")
            for seg in FAULTS:
                x, y = zip(*seg); fig.plot(x=x, y=y, pen="0.7p,gray10")
            for sz, pen in (("0.42c", "0.6p,black"), ("0.32c", "0.5p,black")):
                ev = [e for e in EVENTS if e[3] == sz]
                if ev:
                    fig.plot(x=[e[0] for e in ev], y=[e[1] for e in ev],
                             style=f"a{sz}", fill="yellow", pen=pen)
            for lo, la, txt, _sz, just, off in EVENTS:
                fig.text(x=lo, y=la, text=txt, font="5.5p,Helvetica-Bold,black",
                         justify=just, offset=off, fill="white@25",
                         pen="0.2p,gray40", clearance="1.2p/1.2p+tO")
            fig.colorbar(position="JBC+o0/0.7c+w5.9c/0.25c+h",
                         frame=[f"xa0.5f0.1+l{cblab}"])

fig.savefig(f"{OUTSTEM}.png", dpi=DPI)
fig.savefig(f"{OUTSTEM}.pdf")
print(f"wrote {OUTSTEM}.png ({DPI} dpi) and {OUTSTEM}.pdf")
