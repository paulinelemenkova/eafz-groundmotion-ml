#!/usr/bin/env bash
# =====================================================================
# fig03_stations -- EAFZ strong-motion station coverage map (GMT 6.5)
# Rectangular transverse Mercator (+r), central meridian 38.75E.
# DEM is cut from a LARGER box (34/43.5/35/41) so the bowed-out +r
# corners are filled with genuine relief.
# Stations coloured + sized by per-station record count (jet CPT);
# 200 m brown isolines; red faults; goldenrod1 country borders; gold
# country names; white 10p city labels (+ leader arrows for extras &
# Kahramanmaras); river names; 2-column legend; bar = map height; inset.
# Inputs (same dir): stations.csv, faults.gmt, studybox.txt.
# =====================================================================
set -e
export GMT_END_SHOW=off
IN=.; OUT=.
REGP="34/43.5/35/41"
REG="35/35.5/42.5/40+r"
PROJ="T38.75/37.75/17c"

awk -F',' 'NR>1{s=0.20+0.023*sqrt($3); printf "%s %s %s %.3f\n",$1,$2,$3,s}' "$IN/stations.csv" > stations_gmt.txt
cat > cities.txt << 'EOF'
38.31 38.36 LB Malatya
38.28 37.76 LM Adiyaman
39.22 38.67 LB Elazig
35.32 37.00 LM Adana
36.16 36.20 LM Antakya
40.23 37.91 LM Diyarbakir
36.25 37.07 LM Osmaniye
EOF
for J in LM RB LB RM; do
  awk -v j="$J" '$3==j{printf "%s %s %s\n",$1,$2,substr($0,index($0,$4))}' cities.txt > "cities_$J.txt" || true
done
cat > countries.txt << 'EOF'
39.75 37.25 CM T\374rkiye
38.70 36.20 CM Syria
42.10 36.25 RM Iraq
EOF
cat > new_cities.txt << 'EOF'
41.15 39.05 LT Karliova
40.73 38.70 LT Bingol
39.58 38.57 RT Palu
39.06 38.33 LB Puturge
37.64 38.50 RB Akcadag
36.93 38.43 RB Elbistan
36.22 38.20 RB Goksun
38.92 37.60 LT Kahta
38.06 36.82 LT Nizip
37.36 36.52 LT Kilis
35.78 36.62 RM Iskenderun
36.50 37.86 RB Kahramanmaras
EOF
cat > leaders.txt << 'EOF'
41.15 39.05 41.00 39.30
40.73 38.70 40.50 38.88
39.58 38.57 39.93 38.70
39.06 38.33 38.87 38.19
37.64 38.50 37.96 38.34
36.93 38.43 37.20 38.21
36.22 38.20 36.50 38.02
38.92 37.60 38.62 37.78
38.06 36.82 37.79 37.01
37.36 36.52 37.12 36.72
35.78 36.62 36.17 36.59
36.50 37.86 36.93 37.58
EOF

gmt grdcut @earth_relief_15s -R$REGP -Grelief.nc
gmt grdgradient relief.nc -A315 -Ne0.8 -Gint.nc
gmt makecpt -Cgray -T-2500/4500/100 -Z > relief.cpt
gmt makecpt -Cjet  -T0/420/20       -Z > counts.cpt
MAPH=$(gmt mapproject -R$REG -J$PROJ -Wh | tr -d ' '); MAPH=${MAPH:-12.86}

gmt begin fig03_stations png,pdf
  gmt set FORMAT_GEO_MAP=ddd:mmF MAP_FRAME_TYPE=plain MAP_FRAME_PEN=1.1p,gray20 MAP_FRAME_WIDTH=0.10c \
      MAP_TICK_PEN_PRIMARY=0.6p,gray20 MAP_TICK_LENGTH_PRIMARY=0.14c MAP_GRID_PEN_PRIMARY=0.5p,white \
      FONT_TITLE=15p,Helvetica-Bold,black MAP_TITLE_OFFSET=0.35c FONT_ANNOT_PRIMARY=7p,Helvetica,gray15 FONT_LABEL=9.5p,Helvetica,gray15
  gmt grdimage relief.nc -R$REG -J$PROJ -Crelief.cpt -Iint.nc
  gmt grdcontour relief.nc -R$REG -J$PROJ -C200 -W0.15p,sienna@55 -L0/4500
  gmt coast -R$REG -J$PROJ -Da -A15 -W0.6p,gray15 -N1/0.8p,goldenrod1 -I1/0.5p,steelblue2 -I2/0.35p,steelblue2 -Slightsteelblue1@55
  gmt plot "$IN/faults.gmt" -R$REG -J$PROJ -W0.9p,red
  gmt plot stations_gmt.txt -R$REG -J$PROJ -Si -Ccounts.cpt -W0.35p,black
  echo "38.55 38.52 Euphrates" | gmt text -R$REG -J$PROJ -F+f9p,Helvetica-BoldOblique,steelblue4+a34+jBC -Gwhite@35 -C6%
  echo "39.78 38.80 Murat"     | gmt text -R$REG -J$PROJ -F+f8p,Helvetica-BoldOblique,steelblue4+a-8+jBC -Gwhite@35 -C6%
  echo "40.55 37.72 Tigris"    | gmt text -R$REG -J$PROJ -F+f9p,Helvetica-BoldOblique,steelblue4+a0+jBC -Gwhite@35 -C6%
  gmt plot leaders.txt -R$REG -J$PROJ -Sv0.16c+s+e+a40 -W0.6p,gray20 -Ggray20
  for J in LM RB LB RM; do
    [ -s cities_$J.txt ] && gmt text cities_$J.txt -R$REG -J$PROJ -F+f10p,Helvetica,white+j$J -D0.15c/0c
  done
  gmt text new_cities.txt -R$REG -J$PROJ -F+f10p,Helvetica,white+j
  echo "38.79 37.17 Sanliurfa" | gmt text -R$REG -J$PROJ -F+f10p,Helvetica,white+jLM -D0.15c/-0.40c
  echo "37.38 37.07 Gaziantep" | gmt text -R$REG -J$PROJ -F+f10p,Helvetica,white+jLM -D0.15c/0.34c
  gmt text countries.txt -R$REG -J$PROJ -F+f13p,Helvetica-Bold,gold+j
  gmt basemap -R$REG -J$PROJ -Bxa30mf3mg30m -Bya30mf3mg30m -BwESN+t"Strong-motion station network and record coverage of the EAFZ"
  gmt basemap -R$REG -J$PROJ -Ljbr+w100k+f+o0.5c/0.9c+u --FONT_LABEL=8p,Helvetica,black
  gmt colorbar -Ccounts.cpt -DJML+w${MAPH}c/0.38c+o0.55c/0c+e+ml -Bxa50f5+l"Usable records per station (count)" --FONT_LABEL=10.5p,Helvetica,gray15
  gmt legend -R$REG -J$PROJ -DjTL+w6.0c+o0.2c/0.2c -F+gwhite@10+p0.6p,gray40+r2p --FONT=8p,Helvetica,black << 'EOF'
H 9 Helvetica-Bold Legend
D 0.08c 0.5p
N 2
S 0.34c i 0.30c gray80 0.3p,black 0.72c 25 records
S 0.34c - 0.45c - 1.2p,red 0.80c Active faults
G 0.16c
S 0.34c i 0.42c gray80 0.3p,black 0.72c 100 records
S 0.34c - 0.45c - 1.2p,goldenrod1 0.80c Country border
G 0.20c
S 0.34c i 0.54c gray80 0.3p,black 0.72c 250 records
S 0.34c - 0.45c - 1.2p,steelblue2 0.80c River
G 0.30c
S 0.34c i 0.66c gray80 0.3p,black 0.72c 400 records
G 0.30c
N 1
EOF
  gmt inset begin -DjTR+w2.2c+o0.2c/0.2c -F+gwhite@10+p0.8p,gray30+r2p
     gmt coast -Rg -JG38.75/37.75/2.05c -Gtan -Slightsteelblue1 -A8000 -Bg30 -W0.12p,gray40 --MAP_GRID_PEN_PRIMARY=0.25p,white
     gmt plot "$IN/studybox.txt" -Rg -JG38.75/37.75/2.05c -W1.1p,red -Gred@65 -L
  gmt inset end
  gmt text -R0/17/0/10 -Jx1c -N -F+f7.5p,Helvetica-Oblique,gray25+jBC << 'EOF'
8.5 -0.85 Projection: rectangular transverse Mercator (central meridian 38.75@.E). Stations: representative AFAD/ESM strong-motion network; per-station record counts illustrative. Relief: SRTM/GEBCO 15@+\042@+.
EOF
gmt end
