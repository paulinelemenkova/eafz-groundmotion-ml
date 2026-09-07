#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fig11_classified.py -- ML seismic hazard as classified surface-PGA zones.

The surface-PGA hazard field (rock-level ML PGA x site amplification) is binned
into six named classes by PGA(g) thresholds 0.10/0.20/0.40/0.70/1.00:
    Low (<0.1) | Moderate (0.1-0.2) | High (0.2-0.4) |
    V.High (0.4-0.7) | Severe (0.7-1.0) | Extreme (>=1.0)
Rendered as discrete zones with a categorical colour bar; threshold isolines
are drawn and labelled on the map.

DATA HOOKS:  PRED_FILE = lon lat log10(PGA_rock[g])   (None -> representative)
             VS30_FILE = lon lat Vs30[m/s]            (None -> representative)
REQUIREMENTS: GMT>=6.5, pygmt, numpy      USAGE: python fig11_classified.py
"""
import os, tempfile
import numpy as np
import pygmt

# ---------------------------- configuration ----------------------------
PRED_FILE  = None
VS30_FILE  = None
REGION     = [36.0, 41.0, 36.0, 40.0]
PROJECTION = "T38.75/15c"
SPACING    = 0.02
VREF, KAMP = 760.0, 0.45
THRESH     = [0.10, 0.20, 0.40, 0.70, 1.00]     # PGA(g) class boundaries
OUTSTEM    = "fig11_classified"
DPI        = 400

# ---------------------------- overlays ----------------------------
N_STRAND=[(40.58,39.30),(40.40,39.10),(40.15,38.92),(39.92,38.78),(39.60,38.62),(39.35,38.45)]
S_STRAND=[(39.35,38.45),(39.00,38.25),(38.62,38.06),(38.20,37.88),(37.78,37.78),
          (37.55,37.72),(37.30,37.52),(37.05,37.42),(36.80,37.20),(36.55,36.92),(36.35,36.62),(36.22,36.40)]
CARDAK=[(38.00,38.02),(37.70,38.08),(37.40,38.12),(37.10,38.16),(36.85,38.20)]
DST=[(36.55,36.92),(36.42,36.55),(36.30,36.22),(36.20,36.02)]
FAULTS=[N_STRAND+S_STRAND[1:], CARDAK, DST]
EQ_BIG  =[(37.22,37.17),(37.20,38.02)]
EQ_SMALL=[(39.07,38.36),(40.42,39.04),(40.05,38.80),(36.20,36.20)]
EQ_LABELS=[(37.22,37.17,"Pazarcik M7.8","RM","-0.24c/0.05c","8p"),
           (37.20,38.02,"Elbistan M7.5","RB","-0.22c/0.10c","8p"),
           (39.07,38.36,"Elazig M6.7 (2020)","TL","0.14c/-0.12c","7.5p"),
           (40.42,39.04,"Bingol M6.4 (2003)","BC","0c/0.16c","7.5p")]
CITIES=[("Gaziantep",37.38,37.07,"LB","0.14c/0.10c"),("Kahramanmaras",36.92,37.58,"LB","0.14c/0.10c"),
        ("Malatya",38.31,38.36,"LB","0.14c/0.10c"),("Adiyaman",38.28,37.76,"LB","0.14c/0.10c"),
        ("Elazig",39.22,38.68,"LB","0.14c/0.10c"),("Sanliurfa",38.79,37.17,"LB","0.14c/0.10c"),
        ("Antakya",36.16,36.20,"LB","0.14c/0.10c"),("Diyarbakir",40.23,37.91,"RB","-0.14c/0.10c"),
        ("Bingol",40.50,38.88,"RB","-0.14c/0.10c"),("Tunceli",39.55,39.10,"RB","-0.14c/0.10c"),
        ("Mardin",40.74,37.31,"RB","-0.14c/0.10c"),("Siverek",39.32,37.75,"LB","0.14c/0.10c"),
        ("Osmaniye",36.25,37.07,"TC","0c/-0.12c"),("Kilis",37.12,36.72,"RB","-0.14c/0.10c"),
        ("Golbasi",37.64,37.78,"LB","0.14c/0.10c"),("Celikhan",38.49,38.03,"BC","0c/0.12c"),
        ("Palu",39.93,38.69,"TC","0c/-0.12c"),("Turkoglu",36.85,37.39,"TC","0c/-0.12c")]
COUNTRIES=[("T\u00fcrkiye",37.95,39.62),("Syria",37.85,36.36)]

# ---------------------------- field builders ----------------------------
LAT0=38.0
_lon=np.round(np.arange(REGION[0],REGION[1]+1e-9,SPACING),3)
_lat=np.round(np.arange(REGION[2],REGION[3]+1e-9,SPACING),3)
LON,LAT=np.meshgrid(_lon,_lat)
def _to_km(lo,la):
    return ((np.asarray(lo)-36.)*111.320*np.cos(np.radians(LAT0)),
            (np.asarray(la)-36.)*110.574)
def rock_pga():
    def densify(poly,step=1.0):
        pts=[]
        for (l0,a0),(l1,a1) in zip(poly[:-1],poly[1:]):
            x0,y0=_to_km(l0,a0); x1,y1=_to_km(l1,a1)
            n=max(2,int(np.hypot(x1-x0,y1-y0)/step)+1)
            for t in np.linspace(0,1,n): pts.append((x0+(x1-x0)*t,y0+(y1-y0)*t))
        return np.array(pts)
    gx,gy=_to_km(LON,LAT)
    sysdef={"main":(S_STRAND+CARDAK,0.838),"north":(N_STRAND,0.428),"dst":(DST,0.578)}
    b,h=1.05,6.0; f=np.full(LON.shape,-9.9)
    for poly,a in sysdef.values():
        r2=np.full(LON.shape,np.inf)
        for px,py in densify(poly): np.minimum(r2,(gx-px)**2+(gy-py)**2,out=r2)
        np.maximum(f,a-b*np.log10(np.sqrt(r2)+h),out=f)
    rng=np.random.default_rng(7); co=rng.normal(size=(11,14))
    cy,cx=np.linspace(36,40,11),np.linspace(36,41,14)
    tm=np.vstack([np.interp(_lon,cx,co[i]) for i in range(11)])
    pe=np.vstack([np.interp(_lat,cy,tm[:,j]) for j in range(len(_lon))]).T
    return np.clip(f+0.05*pe,-1.55,0.10)
def vs30_model():
    logV=np.full(LON.shape,np.log10(480.0))
    rng=np.random.default_rng(11); co=rng.normal(size=(14,18))
    cy,cx=np.linspace(36,40,14),np.linspace(36,41,18)
    tm=np.vstack([np.interp(_lon,cx,co[i]) for i in range(14)])
    logV+=0.10*np.vstack([np.interp(_lat,cy,tm[:,j]) for j in range(len(_lon))]).T
    G=lambda lo,la,s: np.exp(-(((LON-lo)**2+(LAT-la)**2)/(2*s**2)))
    for lo,la,s,d in [(36.25,36.35,0.30,0.42),(36.95,37.55,0.26,0.30),(37.40,37.10,0.28,0.22),
        (38.35,38.35,0.26,0.28),(39.25,38.62,0.24,0.25),(40.20,37.90,0.34,0.22),(38.30,37.75,0.22,0.20),
        (38.00,37.40,0.24,0.24),(38.55,37.65,0.22,0.20),(38.85,37.95,0.22,0.18),(37.30,37.50,0.20,0.18)]:
        logV-=d*G(lo,la,s)
    for lo,la,s,d in [(36.30,36.90,0.34,0.18),(36.50,38.20,0.40,0.20),(40.50,39.20,0.45,0.22),
        (39.60,39.30,0.40,0.20),(37.80,38.65,0.34,0.15),(40.55,37.55,0.40,0.12)]:
        logV+=d*G(lo,la,s)
    return np.clip(10**logV,150,900)

if PRED_FILE and os.path.exists(PRED_FILE):
    d=np.loadtxt(PRED_FILE)
    pga=pygmt.surface(x=d[:,0],y=d[:,1],z=d[:,2],region=REGION,spacing=SPACING,tension=0.35).values
else:
    pga=rock_pga()
if VS30_FILE and os.path.exists(VS30_FILE):
    v=np.loadtxt(VS30_FILE)
    vg=pygmt.surface(x=v[:,0],y=v[:,1],z=v[:,2],region=REGION,spacing=SPACING,tension=0.25).values
else:
    vg=vs30_model()
loga=np.clip(KAMP*(np.log10(VREF)-np.log10(vg)),-0.32,0.32)
H=np.clip(pga+loga,-1.5,0.40)
PGAg=10**H                                              # surface PGA (g)
CLS=(np.digitize(PGAg,THRESH)+1).astype(float)          # class index 1..6

pgag_grid=pygmt.xyz2grd(data=np.column_stack([LON.ravel(),LAT.ravel(),PGAg.ravel()]),
                        region=REGION,spacing=SPACING)
cls_grid =pygmt.xyz2grd(data=np.column_stack([LON.ravel(),LAT.ravel(),CLS.ravel()]),
                        region=REGION,spacing=SPACING)

# ---- categorical CPT, class-name colourbar annotations, threshold contours ----
_tmp=tempfile.mkdtemp()
CPT=os.path.join(_tmp,"cls.cpt"); ANNOT=os.path.join(_tmp,"cls_annot.txt"); BND=os.path.join(_tmp,"bounds.txt")
open(CPT,"w").write("0.5 199/233/192 1.5 199/233/192\n1.5 255/255/153 2.5 255/255/153\n"
    "2.5 253/174/97 3.5 253/174/97\n3.5 240/59/32 4.5 240/59/32\n"
    "4.5 178/24/43 5.5 178/24/43\n5.5 103/0/31 6.5 103/0/31\nN gray\n")
open(ANNOT,"w").write("0.5 f\n1 a Low\n1.5 f\n2 a Moderate\n2.5 f\n3 a High\n3.5 f\n"
    "4 a V.High\n4.5 f\n5 a Severe\n5.5 f\n6 a Extreme\n6.5 f\n")
open(BND,"w").write("".join(f"{t:.2f}\n" for t in THRESH))

# ---- render ----
pygmt.config(FONT_ANNOT_PRIMARY="10p,Helvetica",FONT_LABEL="11p,Helvetica",
             MAP_FRAME_TYPE="plain",MAP_FRAME_PEN="0.9p,black",
             MAP_GRID_PEN_PRIMARY="0.25p,gray75",MAP_TICK_LENGTH_PRIMARY="4p",
             PS_CHAR_ENCODING="ISOLatin1+")
fig=pygmt.Figure()
fig.grdimage(grid=cls_grid,cmap=CPT,projection=PROJECTION,region=REGION,frame=["a1f0.1g1","WESN"])
fig.grdcontour(grid=pgag_grid,levels=BND,annotation="+f7p,Helvetica-Bold,gray10+gwhite@30+u g",pen="0.5p,gray20")
fig.coast(water="173/205/226",shorelines="0.4p,gray40",borders="1/1.2p,magenta",area_thresh=30,resolution="i")
for seg in FAULTS:
    x,y=zip(*seg); fig.plot(x=x,y=y,pen="2.6p,white@15")
for seg in FAULTS:
    x,y=zip(*seg); fig.plot(x=x,y=y,pen="1.1p,gray10")
fig.plot(x=[e[0] for e in EQ_BIG],y=[e[1] for e in EQ_BIG],style="a0.60c",fill="yellow",pen="0.8p,black")
fig.plot(x=[e[0] for e in EQ_SMALL],y=[e[1] for e in EQ_SMALL],style="a0.46c",fill="yellow",pen="0.7p,black")
LBL=dict(fill="white@30",pen="0.25p,gray40",clearance="2p/2p+tO")
for name,lo,la,just,off in CITIES:
    fig.text(x=lo,y=la,text=name,font="8p,Helvetica,black",justify=just,offset=off,**LBL)
for lo,la,txt,just,off,fsz in EQ_LABELS:
    fig.text(x=lo,y=la,text=txt,font=f"{fsz},Helvetica-Bold,black",justify=just,offset=off,
             fill="white@20",pen="0.3p,gray40",clearance="2p/2p+tO")
fig.text(x=39.80,y=38.95,text="EAFZ",font="8p,Helvetica-BoldOblique,gray10",justify="LM",**LBL)
for txt,lo,la in COUNTRIES:
    fig.text(x=lo,y=la,text=txt,font="12p,Helvetica-Bold,gray20",justify="CM",fill="white@40",clearance="2p/2p+tO")
fig.colorbar(cmap=CPT,position="JBC+w11c/0.40c+o0/1.05c+h",
             frame=[f"xc{ANNOT}+lSeismic-hazard class  (surface PGA thresholds: 0.1 / 0.2 / 0.4 / 0.7 / 1.0 g)"])
with pygmt.config(FONT_ANNOT_PRIMARY="9p,Helvetica"):
    fig.basemap(map_scale="g39.5/36.25+w100k+f+u+c38")
fig.basemap(rose="g36.25/39.65+w0.9c+l,,,N")
fig.savefig(f"{OUTSTEM}.png",dpi=DPI); fig.savefig(f"{OUTSTEM}.pdf")
print(f"wrote {OUTSTEM}.png and {OUTSTEM}.pdf")
