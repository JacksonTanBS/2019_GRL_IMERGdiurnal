# 2019_GRL_IMERGdiurnal
Codes in support of "Diurnal Cycle of IMERG V06 Precipitation" (10.1029/2019GL085395).

This repository provides the codes used in the paper *Diurnal Cycle of IMERG V06 Precipitation* by Jackson Tan, George J. Huffman, David T. Bolvin, and Eric J. Nelkin, published in the Geophysical Research Letters. These codes allow you to process the data and reproduce the figures in the paper. You will need to retrieve the data separately and update the paths in the codes yourself.

Here are a quick description of the purpose of each script:
* `process_data.py`: produces an intermediate file for a month from the half-hour IMERG HDF5 files.
* `compile_data.py`: reads the monthly intermediate file and computes the diurnal cycle for a specified period and season.
* `plot.py`: plots the figures used in the paper (specific plots will use different outputs from `compile_data.py`).
* `func.py`: a file containing supporting functions for `plot.py`.

These codes have been tested with Python 3.6.8 on Ubuntu 16.04.4.

Jackson Tan  
10 Dec 2019
