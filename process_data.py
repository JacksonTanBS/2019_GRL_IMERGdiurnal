#!/usr/bin/python

import numpy as np
import h5py
import os
import sys
from glob import glob
from calendar import monthrange

# Preliminaries.

year, month = [int(ii) for ii in sys.argv[1 : 3]]
try:
    ver = sys.argv[3]
except IndexError:
    ver = 'V06B'

# Set up the directories.

imergpath = '/path/to/imerg/'
datapath = '/path/to/processed/data/'

os.makedirs(datapath, exist_ok = True)

# Set up the grid.

files = sorted(glob('%s2014/06/01/*.HDF5' % imergpath))

with h5py.File(files[0], 'r') as f:
    lats = f['Grid/lat'][:]
    lons = f['Grid/lon'][:]
    fv = f['Grid/precipitationCal'].attrs['_FillValue']

nlon, nlat = len(lons), len(lats)

# Process the month of data.

filename = '%sIMERG.%s.%4d%02d.npz' % (datapath, ver, year, month)

p = np.zeros([48, nlon, nlat], 'f4')
n = np.zeros([48, nlon, nlat], 'i4')

for day in range(1, monthrange(year, month)[1] + 1):

    fname = '{0:s}{1:4d}/{2:02d}/{3:02d}/*.HDF5'
    files = sorted(glob(fname.format(imergpath, year, month, day)))

    for hhr in range(48):

        with h5py.File(files[hhr], 'r') as f:
            ptmp = f['Grid/precipitationCal'][0, :]
            ptmp[ptmp == fv] = 0
            p[hhr] += ptmp
            n[hhr] += (f['Grid/precipitationCal'][0, :] != fv)

np.savez_compressed(filename, p = p, n = n)
