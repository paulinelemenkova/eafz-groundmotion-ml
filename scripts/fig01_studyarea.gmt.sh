#!/usr/bin/env bash
set -e
export GMT_END_SHOW=off
REG="35/42.5/35.5/40"
PROJ="T38.75/37.75/17c"          # Transverse Mercator, central meridian 38.75E

gmt makecpt -Cglobe -T-1600/3800/100 -Z > topo.cpt
gmt makecpt -Cseis -T0/90/10 -Z         > eqdepth.cpt   # shallow=red -> deep=blue

gmt begin fig01_studyarea png,pdf
  gmt set FORMAT_GEO_MAP=ddd:mmF MAP_FRAME_TYPE=plain \
      MAP_FRAME_PEN=1.1p,gray20 MAP_FRAME_WIDTH=0.10c \
      MAP_TICK_PEN_PRIMARY=0.6p,gray20 MAP_TICK_LENGTH_PRIMARY=0.14c \
      MAP_GRID_PEN_PRIMARY=0.5p,white \
      FONT_TITLE=15p,Helvetica-Bold,black MAP_TITLE_OFFSET=0.4c \
      FONT_ANNOT_PRIMARY=7p,Helvetica,gray15 FONT_LABEL=9.5p,Helvetica,gray15

  # 1. shaded relief ('topo' CPT)
  gmt grdimage eafz_relief.nc -R$REG -J$PROJ -Ctopo.cpt -I+a315+nt0.9

  # 2. topographic isolines every 500 m ('thin' pen)
  gmt grdcontour eafz_relief.nc -R$REG -J$PROJ -C500 -A1000+f5p,Helvetica,gray25 \
      -Wcthin,gray35@45 -Wathin,gray25@25 -L-1600/4000

  # 3. coastline, national borders, rivers (GMT/GSHHG free data)
  gmt coast -R$REG -J$PROJ -Df -A15 -W0.6p,gray15 -N1/0.8p,gray25 \
      -I1/0.8p,steelblue2 -I2/0.5p,steelblue2 -Ir/0.4p,steelblue2

  # 4. tectonic-plate boundaries (Eurasian / Arabian) -- thick coloured, under faults
  gmt plot TP_Eurasian.txt -R$REG -J$PROJ -W2.2p,purple@25
  gmt plot TP_Arabian.txt  -R$REG -J$PROJ -W2.2p,darkorange@20

  # 5. active faults
  gmt plot faults.gmt -R$REG -J$PROJ -W0.9p,black

  # 6. seismicity: colour = focal depth (blue shallow -> red deep), size = magnitude
  gmt plot quakes.txt -R$REG -J$PROJ -i0,1,2,4 -Sc -Ceqdepth.cpt -W0.25p,gray25@40

  # 7. eight major earthquakes: leaders, yellow stars, labels
  gmt plot eq8_leaders.txt -R$REG -J$PROJ -W0.5p,gray30
  gmt plot eq8_stars.txt   -R$REG -J$PROJ -Sa0.55c -Gyellow -W1.1p,black
  gmt text eq8_lab_LM.txt  -R$REG -J$PROJ -F+f7.5p,Helvetica-Bold,black+jLM -Gwhite@25 -W0.2p,gray40 -C16%
  gmt text eq8_lab_RM.txt  -R$REG -J$PROJ -F+f7.5p,Helvetica-Bold,black+jRM -Gwhite@25 -W0.2p,gray40 -C16%

  # 8. cities: Gaziantep (big red), nine others (small magenta)
  gmt plot cities.txt    -R$REG -J$PROJ -Sc0.13c -Gmagenta -W0.3p,black
  gmt plot bigcity.txt   -R$REG -J$PROJ -Sc0.34c -Gmagenta -W0.6p,black
  gmt text cities_lbl.txt  -R$REG -J$PROJ -F+f7p,Helvetica,gray10+jTC -D0/-0.14c -Gwhite@35 -C12%
  gmt text bigcity_lbl.txt -R$REG -J$PROJ -F+f9.5p,Helvetica-Bold,black+jBL -D0.18c/0.16c -Gwhite@20 -W0.3p,gray40 -C18%

  # 9. water-body label (navy)
  gmt text seas.txt -R$REG -J$PROJ -F+a60+f+jLM

  # 9b. country labels
  gmt text countries.txt -R$REG -J$PROJ -F+f13p,Helvetica-Bold,white

  # 10. frame + 30' annotations (10 minor ticks between), white grid, on w E S N
  gmt basemap -R$REG -J$PROJ -Bxa30mf3mg30m -Bya30mf3mg30m \
      -BwESN+t"Seismotectonic setting of the East Anatolian Fault Zone"

  # 11. scale bar (bottom-right)
  gmt basemap -R$REG -J$PROJ -Ljbr+w100k+f+o0.5c/0.9c+u --FONT_LABEL=8p,Helvetica,black

  # 12. vertical earthquake-depth colour bar on the LEFT, close to the map
  gmt colorbar -Ceqdepth.cpt -Dx-1.30c/0.6c+w9.6c/0.34c+e+ml \
      -Bxa15f1+l"Earthquake focal depth (km)"

  # 13. horizontal topography colour bar, full map width
  gmt colorbar -Ctopo.cpt -Dx0c/-1.7c+w17c/0.34c+h+e+ml \
      -Bxa1000f250+l"Topography / bathymetry (m, SRTM/GEBCO 15@+\042@+)" -By

  # 14. legend (top-left), anchored by lower-left corner at 35 deg 06 min E
  gmt legend -R$REG -J$PROJ -Dg35:06E/38:24N+jBL+w3.5c \
      -F+gwhite@12+p0.6p,gray40+r2p --FONT=7p,Helvetica,black << 'EOF'
H 8 Helvetica-Bold Legend
D 0.04c 0.4p
S 0.24c - 0.45c - 2.0p,purple 0.62c Eurasian plate boundary
S 0.24c - 0.45c - 2.0p,darkorange 0.62c Arabian plate boundary
S 0.24c - 0.45c - 0.9p,black 0.62c Active faults
S 0.24c a 0.26c yellow 0.7p,black 0.62c Major earthquake
S 0.24c c 0.30c magenta 0.5p,black 0.62c Gaziantep
S 0.24c c 0.13c magenta 0.3p,black 0.62c City
D 0.04c 0.4p
S 0.24c c 0.10c gray70 0.2p 0.62c M 3
S 0.24c c 0.22c gray70 0.2p 0.62c M 5
S 0.24c c 0.34c gray70 0.2p 0.62c M 6.5
S 0.24c c 0.46c gray70 0.2p 0.62c M 7.8
EOF

  # 15. inset globe (top-right), plate boundaries + study-area box
  gmt inset begin -Dg41:06E/39:00N+jBL+w2.6c -F+gwhite@10+p0.8p,gray30+r2p
     gmt coast -Rg -JG38.75/37.75/2.45c -Gtan -Slightsteelblue1 -A8000 \
         -Bg30 -W0.12p,gray40 --MAP_GRID_PEN_PRIMARY=0.2p,gray80
     gmt plot TP_Eurasian.txt -Rg -JG38.75/37.75/2.45c -W1.0p,purple
     gmt plot TP_Arabian.txt  -Rg -JG38.75/37.75/2.45c -W1.0p,darkorange
     gmt plot studybox.txt    -Rg -JG38.75/37.75/2.45c -W1.1p,red -Gred@65 -L
  gmt inset end

  # 16. projection note + credit below the map
  gmt text -R0/17/0/10 -Jx1c -N -F+f7.5p,Helvetica-Oblique,gray25+jBC << 'EOF'
8.5 -2.62 Projection: Transverse Mercator (central meridian 38.75@.E). SRTM/GEBCO relief; faults: active-fault traces; seismicity: EarthScope IEB (M@-w@- 2.5-7.8, 1970-2026).
EOF
gmt end
