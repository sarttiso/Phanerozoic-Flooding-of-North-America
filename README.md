# Phanerozoic Flooding of North America

[![DOI](https://zenodo.org/badge/666557676.svg)](https://zenodo.org/badge/latestdoi/666557676)

This repository contains the data and Jupyter notebooks for performing the analysis and
creating the figures for the article "Phanerozoic Flooding of North America and the
Great Unconformity".

The `environment.yml` will produce a `conda` environment with almost everything
necessary to run these notebooks. A working installation of
[PyGPlates](https://www.gplates.org/docs/pygplates/) is necessary, and an example
directory for a Windows installation is included with this repository.

The include notebooks are listed below:
- **phanerozoic_coverage.ipynb** : This notebook contains most of the analysis presented
  in the main manuscript, including the figures of fractional flooding area,
  extrapolation of flooding to eustasy, and the statistical analysis of the spatial
  footprint of flooding.
- data_download_load.ipynb : This notebook downloads Macrostrat data via their API and
  saves it as `.gpkg` files.
- modify_macrostrat.ipynb : This notebook performs the modifications described in the
  manuscript (reincorporation of eroded units, updating of Sauk ages, etc.).
- rotate_Greenland.ipynb
- sediment_thickness_estimation.ipynb : This notebook explains the methodology for
  estimating min/max unit thicknesses within time intervals for specific Macrostrat
  columns
- macrostrat_visualizer.ipynb
- cordilleran_validation.ipynb : This notebook performs the comparison between
  Macrostrat and the USGS isopach dataset for the Western Interior Basin.
- gondwana.ipynb : Visualizations of Gondwana.
