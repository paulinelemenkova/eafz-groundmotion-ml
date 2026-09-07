# eafz-gmm-ml

Machine-learning prediction and cartographic mapping of ground-motion intensity
measures along the **East Anatolian Fault Zone (EAFZ)**.

This repository contains the figure-generation code for a study that predicts
ground-motion intensity measures (PGA, PGV, spectral acceleration, Arias intensity
and significant duration) from open-access strong-motion data with gradient-boosted
regression trees, and renders the resulting seismic-hazard fields as continuous maps
with GMT / PyGMT.

## Repository layout

```
eafz-gmm-ml/
├── scripts/                 figure-generation code
│   ├── fig01_studyarea.gmt.sh      study-area / seismotectonic map (GMT)
│   ├── fig02_seismicity.gmt.sh     regional seismicity + focal mechanisms (GMT)
│   ├── fig03_stations.gmt.sh       strong-motion station map (GMT)
│   ├── fig04_datadist.py           magnitude–distance–depth data distribution
│   ├── fig05_workflow.py           methodological workflow diagram (schematic)
│   ├── fig06_mlarch.py             model-architecture diagram (schematic)
│   ├── fig08_pga.py                predicted PGA map (PyGMT)
│   ├── fig09_sa.py                 predicted PGV / Sa panels (PyGMT)
│   ├── fig10_vs30.py               Vs30 site model + amplification (PyGMT)
│   ├── fig11_hazard.py             integrated hazard map (PyGMT)
│   ├── fig12_mmi.py                macroseismic intensity (MMI) map (PyGMT)
│   └── fig14_interpretation.sh     hazard comparison / difference map (GMT)
├── data/
│   ├── fig03_stations.csv          station coordinates and record counts
│   └── meca.txt                    focal-mechanism (GCMT-style) parameters
├── requirements.txt
├── LICENSE
└── README.md
```

## Requirements

- Python 3.12 with the packages in `requirements.txt`
  (`numpy`, `pandas`, `scikit-learn`, `xgboost`, `matplotlib`, `pygmt`)
- **GMT 6.6.0** (for the `*.gmt.sh` and `*.sh` scripts and for PyGMT)

Install the Python side with:

```bash
pip install -r requirements.txt
```

## Usage

Each script is self-contained and writes a `.png` and `.pdf`:

```bash
python scripts/fig08_pga.py           # Python / PyGMT figures
bash   scripts/fig01_studyarea.gmt.sh # GMT shell figures
```

### Using your own model outputs

The map scripts read the predicted fields from plain three-column tables
(`lon lat value`). Point the file variable at the top of each script to your
output and re-run; if it is left as `None`, the script falls back to a
**representative** field so the figure renders standalone. For example, in
`fig08_pga.py`:

```python
PRED_FILE = "predictions_pga.xyz"   # lon  lat  log10(PGA[g])
```

and in `fig10_vs30.py`:

```python
VS30_FILE = "vs30_model.xyz"        # lon  lat  Vs30[m/s]
```

`fig14_interpretation.sh` documents, in its header, the `xyz2grd` commands for
substituting your real hazard grids. Keep the same region (`-R`) and spacing
(`-I0.02`) so the overlays line up.

## Data sources

All input data are openly available: pan-European Engineering Strong-Motion
(ESM) / TADAS strong-motion records, the EarthScope IEB instrumental earthquake
catalogue, Global CMT moment tensors, and the SRTM and GEBCO terrain grids.

## License

Released under the MIT License (see `LICENSE`).

## Citation

If you use this code, please cite the accompanying article
(Lemenkova & Zülfikar, in review).
