#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fig08_pga.py  --  Predicted peak ground acceleration (PGA) across the EAFZ.

Reproduces Figure 8 of this study (Methods Listing 2). Renders a PyGMT map
of predicted log10(PGA) [g] over the East Anatolian Fault Zone, with active-fault
traces, PGA isolines, the 2023 Kahramanmaras epicentres, 14 provincial centres,
country labels, and full cartographic furniture (graticule with 0.1-deg minor
ticks, scale bar, north arrow, colour bar with 0.1 minor ticks).

DATA HOOK
---------
Set PRED_FILE to a whitespace/CSV table with three columns:
        lon  lat  log10(PGA in g)          # one row per model-prediction point
If PRED_FILE is None or missing, the script builds a *representative*,
physically-based field from a distance-attenuation model anchored on the
2023 rupture geometry, so the figure renders without external data.

REQUIREMENTS:  GMT >= 6.5,  pygmt,  numpy        (pip install pygmt numpy)
USAGE:         python fig08_pga.py
"""

import os
import numpy as np
import pygmt

# ----------------------------------------------------------------------------
# 0. Configuration
# ----------------------------------------------------------------------------
PRED_FILE  = None                       # -> "predictions_pga.xyz" for real data
REGION     = [36.0, 41.0, 36.0, 40.0]   # EAFZ bounding box  [W, E, S, N]
PROJECTION = "T38.75/15c"               # Transverse Mercator, central mer. 38.75E
SPACING    = 0.02                       # output grid spacing (deg)
ZMIN, ZMAX = -1.55, 0.10                # colour range of log10(PGA) [g]
CINT       = 0.1                        # isoline / colour-bar minor interval
OUTSTEM    = "fig08_pga"
DPI        = 400

# ----------------------------------------------------------------------------
# 1. Active-fault geometry (lon, lat)
# ----------------------------------------------------------------------------
N_STRAND = [(40.58, 39.30), (40.40, 39.10), (40.15, 38.92), (39.92, 38.78),
            (39.60, 38.62), (39.35, 38.45)]
S_STRAND = [(39.35, 38.45), (39.00, 38.25), (38.62, 38.06), (38.20, 37.88),
            (37.78, 37.78), (37.55, 37.72), (37.30, 37.52), (37.05, 37.42),
            (36.80, 37.20), (36.55, 36.92), (36.35, 36.62), (36.22, 36.40)]
CARDAK   = [(38.00, 38.02), (37.70, 38.08), (37.40, 38.12), (37.10, 38.16),
            (36.85, 38.20)]
DST      = [(36.55, 36.92), (36.42, 36.55), (36.30, 36.22), (36.20, 36.02)]
EAFZ_FULL = N_STRAND + S_STRAND[1:]                 # continuous main strand
FAULTS = [EAFZ_FULL, CARDAK, DST]                   # plotted as cased lines

# epicentres of large events: lon, lat, label, label-justify, label-offset, star-size, font
EPICENTRES = [
    (37.22, 37.17, "Pazarcik M7.8",      "RM", "-0.24c/0.05c", "0.60c", "8p"),
    (37.20, 38.02, "Elbistan M7.5",      "RB", "-0.22c/0.10c", "0.60c", "8p"),
    (39.07, 38.36, "Elazig M6.7 (2020)", "TL",  "0.14c/-0.12c","0.46c", "7.5p"),
    (40.42, 39.04, "Bingol M6.4 (2003)", "BC",  "0c/0.16c",    "0.46c", "7.5p"),
]

# 14 provincial centres: name, lon, lat, label-justify, label-offset
#   (offsets chosen so labels do not overlap each other or the data)
CITIES = [
    ("Gaziantep",     37.38, 37.07, "LB",  "0.14c/0.10c", True),
    ("Kahramanmaras", 36.92, 37.58, "LB",  "0.14c/0.10c", False),
    ("Malatya",       38.31, 38.36, "LB",  "0.14c/0.10c", False),
    ("Adiyaman",      38.28, 37.76, "LB",  "0.14c/0.10c", False),
    ("Elazig",        39.22, 38.68, "LB",  "0.14c/0.10c", False),
    ("Sanliurfa",     38.79, 37.17, "LB",  "0.14c/0.10c", False),
    ("Antakya",       36.16, 36.20, "LB",  "0.14c/0.10c", False),
    ("Diyarbakir",    40.23, 37.91, "RB", "-0.14c/0.10c", False),
    ("Bingol",        40.50, 38.88, "RB", "-0.14c/0.10c", False),
    ("Tunceli",       39.55, 39.10, "RB", "-0.14c/0.10c", False),
    ("Mardin",        40.74, 37.31, "RB", "-0.14c/0.10c", False),
    ("Siverek",       39.32, 37.75, "LB",  "0.14c/0.10c", False),
    ("Osmaniye",      36.25, 37.07, "TC",  "0c/-0.12c",   False),  # moved below symbol
    ("Kilis",         37.12, 36.72, "RB", "-0.14c/0.10c", False),
    # --- 4 more along the EAFZ, placed in gaps ---
    ("Golbasi",       37.64, 37.78, "LB",  "0.14c/0.10c", False),
    ("Celikhan",      38.49, 38.03, "BC",  "0c/0.12c",    False),
    ("Palu",          39.93, 38.69, "TC",  "0c/-0.12c",   False),
    ("Turkoglu",      36.85, 37.39, "TC",  "0c/-0.12c",   False),
]

# country labels: text, lon, lat  (placed on each side of the border)
COUNTRIES = [("T\u00fcrkiye", 37.95, 39.62), ("Syria", 37.85, 36.36)]

# ----------------------------------------------------------------------------
# 2. Predicted-PGA point cloud  (real file  OR  representative model)
# ----------------------------------------------------------------------------
def representative_field():
    """log10PGA = a - b*log10(R+h), enveloped over the fault systems, with the
    2023 rupture dominant. Returns lon, lat, val (1-D arrays)."""
    LAT0, b, h = 38.0, 1.05, 6.0
    systems = {"main":  (S_STRAND + CARDAK, 0.838),   # 2023 rupture (dominant)
               "north": (N_STRAND,          0.428),   # northern EAFZ strand
               "dst":   (DST,               0.578)}   # Dead Sea / Amanos transform

    def to_km(lon, lat):
        return ((np.asarray(lon) - 36.0) * 111.320 * np.cos(np.radians(LAT0)),
                (np.asarray(lat) - 36.0) * 110.574)

    def densify(poly, step=1.0):
        pts = []
        for (l0, a0), (l1, a1) in zip(poly[:-1], poly[1:]):
            x0, y0 = to_km(l0, a0); x1, y1 = to_km(l1, a1)
            n = max(2, int(np.hypot(x1 - x0, y1 - y0) / step) + 1)
            for t in np.linspace(0, 1, n):
                pts.append((x0 + (x1 - x0) * t, y0 + (y1 - y0) * t))
        return np.array(pts)

    lon = np.round(np.arange(REGION[0], REGION[1] + 1e-9, SPACING), 3)
    lat = np.round(np.arange(REGION[2], REGION[3] + 1e-9, SPACING), 3)
    LON, LAT = np.meshgrid(lon, lat)
    GX, GY = to_km(LON, LAT)

    field = np.full(LON.shape, -9.9)
    for poly, a in systems.values():
        r2 = np.full(LON.shape, np.inf)
        for px, py in densify(poly):
            np.minimum(r2, (GX - px) ** 2 + (GY - py) ** 2, out=r2)
        np.maximum(field, a - b * np.log10(np.sqrt(r2) + h), out=field)

    rng = np.random.default_rng(7)
    coarse = rng.normal(size=(11, 14))
    cy, cx = np.linspace(36, 40, 11), np.linspace(36, 41, 14)
    tmp = np.vstack([np.interp(lon, cx, coarse[i]) for i in range(11)])
    pert = np.vstack([np.interp(lat, cy, tmp[:, j]) for j in range(len(lon))]).T
    field = np.clip(field + 0.05 * pert, ZMIN, ZMAX)
    return LON.ravel(), LAT.ravel(), field.ravel()


if PRED_FILE and os.path.exists(PRED_FILE):
    pred = np.loadtxt(PRED_FILE)
    lon, lat, val = pred[:, 0], pred[:, 1], pred[:, 2]
    print(f"Loaded {len(val)} predictions from {PRED_FILE}")
else:
    lon, lat, val = representative_field()
    print("PRED_FILE not set/found -> using representative attenuation field")

grid = pygmt.surface(x=lon, y=lat, z=val,
                     region=REGION, spacing=SPACING, tension=0.35)
grid = pygmt.grdclip(grid, below=[ZMIN, ZMIN], above=[ZMAX, ZMAX])

# ----------------------------------------------------------------------------
# 3. Render
# ----------------------------------------------------------------------------
pygmt.config(FONT_ANNOT_PRIMARY="10p,Helvetica", FONT_LABEL="11p,Helvetica",
             MAP_FRAME_TYPE="plain", MAP_FRAME_PEN="0.9p,black",
             MAP_GRID_PEN_PRIMARY="0.25p,gray75", MAP_TICK_LENGTH_PRIMARY="4p",
             PS_CHAR_ENCODING="ISOLatin1+")

fig = pygmt.Figure()
pygmt.makecpt(cmap="lajolla", series=[ZMIN, ZMAX, 0.05], reverse=True)

# raster + 0.1-deg minor ticks on the graticule (a1 f0.1 g1)
fig.grdimage(grid=grid, projection=PROJECTION, region=REGION,
             frame=["a1f0.1g1", "WESN"], cmap=True)

# PGA isolines every 0.1 (thin brown).  Older PyGMT: use interval=CINT
fig.grdcontour(grid=grid, levels=CINT, pen="0.3p,139/69/19")

fig.coast(water="173/205/226", shorelines="0.4p,gray40",
          borders="1/1.2p,magenta", area_thresh=30, resolution="i")

# cased fault traces (white halo + dark core)
for seg in FAULTS:
    x, y = zip(*seg); fig.plot(x=x, y=y, pen="2.6p,white@15")
for seg in FAULTS:
    x, y = zip(*seg); fig.plot(x=x, y=y, pen="1.1p,gray10")

# epicentres (stars) + provincial centres (circles; Gaziantep enlarged)
for size, pen in (("0.60c", "0.8p,black"), ("0.46c", "0.7p,black")):
    pts = [e for e in EPICENTRES if e[5] == size]
    if pts:
        fig.plot(x=[e[0] for e in pts], y=[e[1] for e in pts],
                 style=f"a{size}", fill="yellow", pen=pen)
fig.plot(x=[c[1] for c in CITIES], y=[c[2] for c in CITIES],
         style="c0.13c", fill="magenta", pen="0.4p,black")
fig.plot(x=[37.38], y=[37.07], style="c0.22c", fill="magenta", pen="0.7p,black")

# city labels (8 pt, boxed; per-city justify/offset to avoid overlaps)
for name, lo, la, just, off, _big in CITIES:
    fig.text(x=lo, y=la, text=name, font="8p,Helvetica,black",
             justify=just, offset=off, fill="white@25",
             pen="0.25p,gray40", clearance="2p/2p+tO")

# epicentre labels (font sized per event)
for lo, la, txt, just, off, _sz, fsz in EPICENTRES:
    fig.text(x=lo, y=la, text=txt, font=f"{fsz},Helvetica-Bold,black",
             justify=just, offset=off, fill="white@20",
             pen="0.3p,gray40", clearance="2p/2p+tO")

# fault label
fig.text(x=39.80, y=38.95, text="EAFZ", font="8p,Helvetica-BoldOblique,gray10",
         justify="LM", fill="white@25", pen="0.25p,gray40", clearance="2p/2p+tO")

# country names (placed by side of the border)
for txt, lo, la in COUNTRIES:
    fig.text(x=lo, y=la, text=txt, font="12p,Helvetica-Bold,gray20",
             justify="CM", fill="white@40", clearance="2p/2p+tO")

# colour bar with 0.1 minor ticks
fig.colorbar(position="JBC+w11c/0.35c+o0/1.05c+h",
             frame=r'xa0.5f0.1+l"Predicted log@-10@-(PGA),  PGA in g"')

# scale bar: bottom-right at 39 deg 30 min E, NO box, font 1 pt smaller
with pygmt.config(FONT_ANNOT_PRIMARY="9p,Helvetica"):
    fig.basemap(map_scale="g39.5/36.25+w100k+f+u+c38")

# north arrow: starting 36 deg 15 min E, pulled inside, NO box
fig.basemap(rose="g36.25/39.65+w0.9c+l,,,N")

fig.savefig(f"{OUTSTEM}.png", dpi=DPI)
fig.savefig(f"{OUTSTEM}.pdf")
print(f"wrote {OUTSTEM}.png ({DPI} dpi) and {OUTSTEM}.pdf")
