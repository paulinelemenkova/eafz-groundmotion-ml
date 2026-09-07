#!/usr/bin/env bash
# =====================================================================
# fig14_interpretation -- 3-panel interpretive hazard comparison (EAFZ)
#   (a) ML seismic-hazard estimate   (b) independent reference hazard
#   (c) difference (ML - reference)
# Rendered with GMT 6 (modern mode).  Run:  bash fig14_interpretation.sh
#
# >>> IMPORTANT <<<
# The three grids below (ml.nc, ref.nc, diff.nc) are produced here from an
# ILLUSTRATIVE synthetic field so the figure renders standalone. Replace the
# "DATA" block with YOUR real grids before publication, e.g.:
#     gmt xyz2grd ml_hazard.xyz  -Gml.nc  -R$R -I0.02
#     gmt xyz2grd ref_hazard.xyz -Gref.nc -R$R -I0.02
#     gmt grdmath ml.nc ref.nc SUB = diff.nc
# Keep the same -R/-I so the overlays line up. faults.txt should hold your
# digitised active-fault traces (multisegment, ">"-separated lon lat).
# =====================================================================
set -e
export GMT_DATA_SERVER=oceania
R=35.3/39.7/35.8/39.2

# ---------------------- DATA (replace with real grids) ----------------
python3 - <<'PY'
import numpy as np
W,E,S,N,dx=35.3,39.7,35.8,39.2,0.02
lons=np.linspace(W,E,round((E-W)/dx)+1); lats=np.linspace(S,N,round((N-S)/dx)+1)
LON,LAT=np.meshgrid(lons,lats); cl=np.cos(np.radians(37.5))
main=[(39.55,39.05),(39.0,38.7),(38.5,38.4),(38.1,38.2),(37.8,38.0),(37.55,37.75),
      (37.35,37.45),(37.25,37.15),(37.0,36.9),(36.7,36.65),(36.4,36.45),(36.2,36.25),(36.05,36.0)]
branch=[(38.1,38.2),(37.7,38.12),(37.3,38.08),(36.9,38.0),(36.6,37.95)]
def dens(p,st=0.01):
    o=[]
    for (x0,y0),(x1,y1) in zip(p[:-1],p[1:]):
        n=max(2,int(np.hypot(x1-x0,y1-y0)/st)); o+=list(zip(np.linspace(x0,x1,n),np.linspace(y0,y1,n)))
    return np.array(o)
F=np.vstack([dens(main),dens(branch)]); d=np.full(LON.shape,1e9)
for fx,fy in F: np.minimum(d,np.hypot((LON-fx)*cl,(LAT-fy)),out=d)
dk=d*111.0
cities=[("Kahramanmaras",37.58,37.57,.35,.32,.18),("Gaziantep",37.38,37.07,.30,.30,.16),("Antakya",36.16,36.20,.33,.34,.16)]
def bump(k,ex):
    b=np.zeros(LON.shape)
    for _,cx,cy,aml,aref,s in cities:
        a=aml if k=="ml" else aref; b+=a*np.exp(-(((LON-cx)*cl)**2+(LAT-cy)**2)/(2*(s+ex)**2))
    return b
g=0.06*((E-LON)/(E-W))+0.04*((N-LAT)/(N-S))
ml=np.clip(0.15+0.80*np.exp(-(dk/12.)**2)+bump("ml",0)+g,0.05,1.2)
ref=np.clip(0.15+0.70*np.exp(-(dk/18.)**2)+bump("ref",0.06)+g,0.05,1.2)
def sv(n,Z): np.savetxt(n,np.column_stack([LON.ravel(),LAT.ravel(),Z.ravel()]),fmt="%.4f %.4f %.5f")
sv("ml.xyz",ml); sv("ref.xyz",ref); sv("diff.xyz",ml-ref)
open("faults.txt","w").write(">\n"+"\n".join("%.3f %.3f"%p for p in main)+"\n>\n"+"\n".join("%.3f %.3f"%p for p in branch)+"\n")
open("cities.txt","w").write("".join("%.3f %.3f %s\n"%(cx,cy,nm) for nm,cx,cy,*_ in cities))
open("epi.txt","w").write("37.29 37.17 Pazarcik\n37.22 38.02 Elbistan\n")
PY
for v in ml ref diff; do gmt xyz2grd $v.xyz -G$v.nc -R$R -I0.02; done
# ---------------------- RENDER ---------------------------------------
gmt begin fig14_interpretation png,pdf
  gmt set FONT_ANNOT_PRIMARY 7p FONT_LABEL 8p MAP_FRAME_TYPE plain \
          FORMAT_GEO_MAP dddF MAP_FRAME_PEN 0.8p MAP_TICK_LENGTH_PRIMARY 0.07c MAP_ANNOT_OFFSET 0.08c
  gmt makecpt -Cinferno -T0.1/1.0/0.05 -H > haz.cpt
  gmt makecpt -Cvik -T-0.3/0.3/0.025 -H > diff.cpt
  gmt subplot begin 1x3 -Fs6.0c/0 -M0.4c/1.5c -R$R -JM6.0c \
      -Bxa1f0.5 -Bya1f0.5 -BWSne -SRl -A"(a)+jTL+gwhite@20+p0.4p,black+o0.1c"
    LAYER='gmt coast -N1/0.9p,magenta3 -W0.3p,gray30 -Di
           gmt plot faults.txt -W2.2p,white
           gmt plot faults.txt -W0.9p,black
           gmt plot epi.txt -i0,1 -Sa0.34c -Gyellow -W0.4p,black
           gmt plot cities.txt -i0,1 -Ss0.15c -Gwhite -W0.5p,black
           gmt text cities.txt -F+f6.2p,Helvetica-Bold,white+jBL -Gblack@35 -D0.09c/0.07c'
    gmt subplot set 0
      gmt grdimage ml.nc -Chaz.cpt; eval "$LAYER"
      echo "ML hazard estimate" | gmt text -F+cTC+jTC+f8p,Helvetica-Bold -Gwhite@20 -W0.3p,black -D0/-0.16c
      gmt colorbar -Chaz.cpt -DJBC+w4.8c/0.22c+h -Bxa0.2f0.1+l"PGA (g)"
    gmt subplot set 1
      gmt grdimage ref.nc -Chaz.cpt; eval "$LAYER"
      echo "Reference hazard" | gmt text -F+cTC+jTC+f8p,Helvetica-Bold -Gwhite@20 -W0.3p,black -D0/-0.16c
      gmt colorbar -Chaz.cpt -DJBC+w4.8c/0.22c+h -Bxa0.2f0.1+l"PGA (g)"
    gmt subplot set 2
      gmt grdimage diff.nc -Cdiff.cpt; eval "$LAYER"
      echo "Difference (ML @%2%-@%% reference)" | gmt text -F+cTC+jTC+f8p,Helvetica-Bold -Gwhite@20 -W0.3p,black -D0/-0.16c
      gmt colorbar -Cdiff.cpt -DJBC+w4.8c/0.22c+h -Bxa0.15f0.05+l"@~D@~PGA (g)"
  gmt subplot end
gmt end
