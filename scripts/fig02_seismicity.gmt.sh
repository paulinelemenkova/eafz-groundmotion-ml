#!/usr/bin/env bash
# =====================================================================
#  fig02_seismicity_gmt.sh
#  Regional instrumental seismicity + focal mechanisms of the EAFZ.
#
#  Reproduces the figure: 32 focal-mechanism "beachballs" derived from the
#  EarthScope IEB catalogue (IEB_EAFZ.csv), coloured by hypocentral depth,
#  sized by magnitude, rotated by the local EAFZ segment strike, with the
#  12 largest events labelled.
#
#  Run:   bash fig02_seismicity_gmt.sh        (needs GMT 6 + python3)
#
#  >>> IMPORTANT, PLEASE READ >>>
#  The IEB catalogue contains NO focal-mechanism angles, so the mechanisms
#  in meca.txt are REPRESENTATIVE of the EAFZ left-lateral strike-slip
#  regime (near-vertical dip, mild transtension; normal-oblique at the
#  Hazar / Golbasi / Samsat releasing stepovers), with strike taken from
#  the local fault segment. The two 2023 mainshocks use their published
#  orientations. They are NOT GCMT inversions -- replace meca.txt with real
#  GCMT/RMT solutions (same psmeca Aki columns) before publication.
#
#  The relief and fault traces here are stand-ins so the script runs
#  anywhere: relief is fetched live (@earth_relief_30s) and faults.txt is
#  an approximate trace. In your own pipeline, swap in eafz_relief.nc,
#  faults.gmt and the TP_*.txt plate-boundary files you already use.
# =====================================================================
set -e
export GMT_DATA_SERVER=oceania            # relief tiles + DCW national borders
CSV="${1:-IEB_EAFZ.csv}"                  # input catalogue
REG="35/42.5/35.5/40"
PROJ="T38.75/37.75/17c"
NSEL=32                                   # number of beachballs (28-38 ok)
MINMAG_BG=2.5                             # background circles: min magnitude

# ---------------------------------------------------------------------
# 1. Build overlays from the IEB catalogue
#    meca.txt     lon lat depth strike dip rake mag   (psmeca -Sa)
#    meca_lbl.txt lon lat justify "label"             (12 largest)
#    quakes.txt   lon lat depth mag size_cm           (background)
#    faults.txt   approximate EAFZ traces (replace with faults.gmt)
# ---------------------------------------------------------------------
python3 - "$CSV" "$NSEL" "$MINMAG_BG" <<'PYEOF'
import csv, math, sys
SRC, NSEL, MINMAG_BG = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])

rows = []
with open(SRC, newline='') as f:
    for r in csv.DictReader(f):
        try:
            rows.append(dict(year=int(r['Year']), lat=float(r['Lat']), lon=float(r['Lon']),
                             depth=float(r['Depth']), mag=float(r['Mag']),
                             region=(r['Region'] or '').strip()))
        except Exception:
            pass

# greedy spatial declustering, largest-magnitude first
rows.sort(key=lambda x: -x['mag'])
coslat = math.cos(math.radians(37.7))
sel = []
for e in rows:
    if all(math.hypot((e['lon']-s['lon'])*coslat, e['lat']-s['lat']) >= 0.15 for s in sel):
        sel.append(e)
    if len(sel) >= NSEL:
        break

# ---- representative mechanism from the local EAFZ segment ----
def seg_strike(lon, lat):
    if lon >= 40.0:                         return 72   # Karliova-Bingol-Yedisu junction
    if 38.7 <= lon < 40.0:                  return 62   # Palu-Hazar-Puturge main strand
    if lat >= 37.9 and 36.3 <= lon < 38.0:  return 92   # Cardak/Surgu branch (E-W)
    if 37.4 <= lon < 38.7:                  return 55   # Celikhan-Sincik-Samsat
    if 36.5 <= lon < 37.4:                  return 30   # Pazarcik-Nurdagi (NNE)
    if lon < 36.5:                          return 22   # Amanos / southern termination
    return 50

def jitter(e, amp):                                     # deterministic per-event scatter
    h = (int(e['lon']*1000)*73856093) ^ (int(e['lat']*1000)*19349663) ^ (e['year']*83492791)
    return ((h % 1000)/1000.0 - 0.5) * 2 * amp

extensional = [(39.30, 38.50), (37.73, 37.80), (38.43, 37.62)]   # releasing stepovers
def mech(e):
    strike = (seg_strike(e['lon'], e['lat']) + jitter(e, 8)) % 180
    dip    = max(60, min(89, 78 + jitter(e, 9)))
    rake   = -6 + jitter(e, 12)
    if any(math.hypot((e['lon']-x)*coslat, e['lat']-y) < 0.25 for x, y in extensional):
        rake = -35 + jitter(e, 8)
    return round(strike), round(dip), round(rake)

override = {(37.0143, 37.2256): (28, 80, -8),    # Pazarcik Mw 7.8
            (37.1962, 38.0106): (92, 78, -4)}    # Elbistan  Mw 7.5

with open("meca.txt", "w") as fm:                # small first -> big plotted on top
    for e in sorted(sel, key=lambda x: x['mag']):
        s, d, r = override.get((round(e['lon'],4), round(e['lat'],4)), mech(e))
        fm.write(f"{e['lon']:.4f} {e['lat']:.4f} {e['depth']:.1f} {s} {d} {r} {e['mag']:.1f}\n")

# ---- labels for the 12 largest: leader arrows + offset anchors ----
# LABELS are hand-placed (direction, distance) into open space to avoid any
# overlap, with an arrow pointing back at the beachball. They are tuned for
# THIS catalogue/NSEL=32; if you change the data, adjust direction/distance.
DIR = {'E':(1,0,'ML'),'W':(-1,0,'MR'),'N':(0,1,'BC'),'S':(0,-1,'TC'),
       'NE':(0.7,0.7,'BL'),'NW':(-0.7,0.7,'BR'),'SE':(0.7,-0.7,'TL'),'SW':(-0.7,-0.7,'TR')}
brad = lambda m: 0.265 if m>=7.0 else (0.176 if m>=6.3 else 0.115)   # ball radius (deg)
LABELS = [
    (35.307,36.878,6.3,'E', 0.55,"Mw 6.3 (1998)"),
    (36.025,36.162,6.3,'SE',0.55,"Mw 6.3 (2023)"),
    (37.014,37.226,7.8,'S', 0.95,"Pazarcik Mw 7.8"),
    (37.196,38.011,7.5,'NW',0.85,"Elbistan Mw 7.5"),
    (37.806,37.993,6.1,'N', 0.78,"Mw 6.1 (1986)"),
    (38.098,38.032,6.0,'S', 0.78,"Mw 6.0 (2023)"),
    (39.061,38.431,6.7,'SE',0.72,"Mw 6.7 (2020)"),
    (39.878,39.500,6.1,'W', 0.78,"Mw 6.1 (2003)"),
    (39.986,38.864,6.1,'SW',0.72,"Mw 6.1 (2010)"),
    (40.464,39.007,6.4,'NE',0.66,"Mw 6.4 (2003)"),
    (40.653,38.934,6.6,'E', 0.85,"Mw 6.6 (1971)"),
    (40.723,38.474,6.7,'E', 0.80,"Mw 6.7 (1975)"),
]
with open("leaders.txt","w") as fld, open("lbl.txt","w") as fl:
    for lon,lat,m,d,dist,txt in LABELS:
        ux,uy,just = DIR[d]; fc = math.cos(math.radians(lat))
        Llon = lon + dist*ux/fc; Llat = lat + dist*uy
        vx=(Llon-lon)*fc; vy=Llat-lat; n=math.hypot(vx,vy); r=brad(m)
        elon = lon + (vx/n)*r/fc; elat = lat + (vy/n)*r   # arrow tip at ball edge
        fld.write(f"{Llon:.4f} {Llat:.4f} {elon:.4f} {elat:.4f}\n")
        fl.write(f"{Llon:.4f} {Llat:.4f} {just} {txt}\n")

# ---- background seismicity (size grows with magnitude) ----
def size_cm(m): return round(0.05 * (1.5 ** (m - 3.0)), 3)
with open("quakes.txt", "w") as fq:
    for e in rows:
        if e['mag'] >= MINMAG_BG:
            fq.write(f"{e['lon']:.4f} {e['lat']:.4f} {e['depth']:.1f} {e['mag']:.1f} {size_cm(e['mag'])}\n")

# ---- approximate EAFZ fault traces (preview only; use your faults.gmt) ----
main   = [(40.620,39.200),(40.050,38.920),(39.550,38.700),(39.200,38.550),(38.900,38.420),
          (38.550,38.280),(38.150,38.160),(37.820,38.000),(37.600,37.850),(37.420,37.620),
          (37.260,37.320),(37.100,37.060),(36.900,36.860),(36.620,36.600),(36.400,36.400),(36.220,36.180)]
branch = [(38.150,38.160),(37.750,38.100),(37.320,38.060),(36.920,38.000),(36.550,38.050)]
south  = [(36.400,36.400),(36.300,36.050),(36.180,35.700)]
nafj   = [(40.620,39.200),(40.200,39.450),(39.880,39.550)]
with open("faults.txt", "w") as ff:
    for seg in (main, branch, south, nafj):
        ff.write(">\n" + "\n".join(f"{x:.3f} {y:.3f}" for x, y in seg) + "\n")

print(f"meca={len(sel)}  labels={len(LABELS)}  background>={MINMAG_BG}")
PYEOF

# ---------------------------------------------------------------------
# 2. Backdrop grids / CPTs
#    Replace the grdcut line with:  (use your own relief)
#        ln -sf eafz_relief.nc eafz_relief.nc
# ---------------------------------------------------------------------
gmt grdcut @earth_relief_30s -R$REG -Geafz_relief.nc
gmt makecpt -Cgray -T-3000/5500          > relief_gray.cpt
gmt makecpt -Cseis -T0/40/5 -Z           > eqdepth.cpt    # shallow=red -> deep=blue

# split mechanisms into magnitude classes for graduated (contrasting) sizes
awk '$7>=7.0'           meca.txt > m_vbig.txt   # very big  (2023 doublet)
awk '$7>=6.3 && $7<7.0' meca.txt > m_big.txt    # big
awk '$7>=5.8 && $7<6.3' meca.txt > m_med.txt    # medium
awk '$7<5.8'            meca.txt > m_sml.txt    # small

# ---------------------------------------------------------------------
# 3. Render
# ---------------------------------------------------------------------
gmt begin fig02_seismicity png,pdf
  gmt set FORMAT_GEO_MAP=ddd:mmF MAP_FRAME_TYPE=plain \
      MAP_FRAME_PEN=1.1p,gray20 MAP_FRAME_WIDTH=0.10c \
      MAP_TICK_PEN_PRIMARY=0.6p,gray20 MAP_TICK_LENGTH_PRIMARY=0.14c \
      MAP_GRID_PEN_PRIMARY=0.5p,white \
      FONT_TITLE=15p,Helvetica-Bold,black MAP_TITLE_OFFSET=0.4c \
      FONT_ANNOT_PRIMARY=7p,Helvetica,gray15 FONT_LABEL=9.5p,Helvetica,gray15

  # muted grayscale shaded relief
  gmt grdimage eafz_relief.nc -R$REG -J$PROJ -Crelief_gray.cpt -I+a315+nt0.8

  # coast, borders, rivers (subtle)
  gmt coast -R$REG -J$PROJ -Di -A15 -W0.5p,gray30 -N1/0.7p,gray35,- \
      -I1/0.5p,lightsteelblue3 -Ir/0.3p,lightsteelblue3

  # active fault segments (replace faults.txt with your faults.gmt;
  # add: gmt plot TP_Eurasian.txt -W1.4p,purple@30 ; TP_Arabian.txt -W1.4p,darkorange@25)
  gmt plot faults.txt -R$REG -J$PROJ -W0.9p,black

  # background seismicity: colour = depth, size = magnitude (column 4)
  gmt plot quakes.txt -R$REG -J$PROJ -i0,1,2,4 -Sc -Ceqdepth.cpt -W0.25p,gray25@45

  # focal mechanisms, graduated size by magnitude class (-M = fixed size per
  # class), coloured by depth (-Z), rotated by strike; small drawn first
  [ -s m_sml.txt ]  && gmt meca m_sml.txt  -R$REG -J$PROJ -Sa0.34c+f0p -M -Zeqdepth.cpt -W0.4p,black -L0.35p,black -N
  [ -s m_med.txt ]  && gmt meca m_med.txt  -R$REG -J$PROJ -Sa0.52c+f0p -M -Zeqdepth.cpt -W0.5p,black -L0.4p,black  -N
  [ -s m_big.txt ]  && gmt meca m_big.txt  -R$REG -J$PROJ -Sa0.80c+f0p -M -Zeqdepth.cpt -W0.7p,black -L0.5p,black  -N
  [ -s m_vbig.txt ] && gmt meca m_vbig.txt -R$REG -J$PROJ -Sa1.20c+f0p -M -Zeqdepth.cpt -W0.9p,black -L0.6p,black  -N

  # labels for the 12 most important: arrowed leaders + offset anchors
  gmt plot leaders.txt -R$REG -J$PROJ -Sv0.18c+e+s -W0.8p,gray10 -Ggray10 -N
  gmt text  lbl.txt    -R$REG -J$PROJ -F+f7p,Helvetica-Bold,black+j \
      -Gwhite@12 -W0.2p,gray40 -C8% -N

  # frame, white grid, 30' annotations
  gmt basemap -R$REG -J$PROJ -Bxa30mf3mg30m -Bya30mf3mg30m \
      -BwESN+t"Regional instrumental seismicity and focal mechanisms of the EAFZ"

  # scale bar
  gmt basemap -R$REG -J$PROJ -Ljbr+w100k+f+o0.5c/0.9c+u --FONT_LABEL=8p,Helvetica,black

  # depth colour bar (closer to map; minor ticks every 1 km = 10 per major)
  gmt colorbar -Ceqdepth.cpt -Dx-1.05c/0.6c+w9.6c/0.34c+e+ml \
      -Bxa10f1+l"Hypocentral depth (km)"

  # focal-mechanism size key (circle proxies at the four tier sizes), SE corner
  gmt legend -R$REG -J$PROJ -Dg40:00E/36:06N+jBL+w4.7c \
      -F+gwhite@12+p0.6p,gray40+r2p --FONT=7.5p,Helvetica,black << 'EOF'
H 8 Helvetica-Bold Magnitude (beachball size)
G 0.62c
N 1
S 0.70c c 1.20c gray85 0.5p,black 1.55c M@-w@- >= 7.0
S 0.70c c 0.80c gray85 0.5p,black 1.55c 6.3 - 6.9
S 0.70c c 0.52c gray85 0.4p,black 1.55c 5.8 - 6.2
S 0.70c c 0.34c gray85 0.4p,black 1.55c 5.2 - 5.7
D 0.05c 0.4p
L 7 Helvetica-Oblique L Beachball fill = hypocentral depth
EOF

  # projection note + credit
  gmt text -R0/17/0/10 -Jx1c -N -F+f7p,Helvetica-Oblique,gray25+jBC << 'EOF'
8.5 -1.15 Projection: Transverse Mercator (CM 38.75@.E). Seismicity: EarthScope IEB (M@-w@- 2.5-7.8, 1970-2026). Focal mechanisms: left-lateral regime (strike from local fault segment). Relief: SRTM/GEBCO.
EOF
gmt end
echo "wrote fig02_seismicity.png / .pdf"
