#!/usr/bin/python

import numpy as np
import h5py
import sys
from glob import glob

# Preliminaries.

year0, month0, year1, month1 = [int(ii) for ii in sys.argv[1 : 5]]
season = sys.argv[5]
try:
    ver = sys.argv[6]
except IndexError:
    ver = 'V06B'

# Set up the directories.

imergpath = '/path/to/imerg/'
datapath = '/path/to/processed/data/'

# Set up the grid.

files = sorted(glob('%s2014/06/01/*.HDF5' % imergpath))

with h5py.File(files[0], 'r') as f:
    lats = f['Grid/lat'][:]
    lons = f['Grid/lon'][:]

nlon, nlat = len(lons), len(lats)

# Define the required functions.

def getMonths(year, year0, month0, year1, month1):
    if   year0 == year1: return range(month0, month1 + 1)
    elif year  == year0: return range(month0, 13)
    elif year  == year1: return range(1, month1 + 1)
    else               : return range(1, 13)

def monthInSeason(month, season):
    if ((season == 'DJF' and month in (12, 1, 2)) or
        (season == 'MAM' and month in (3, 4, 5)) or
        (season == 'JJA' and month in (6, 7, 8)) or
        (season == 'SON' and month in (9, 10, 11))):
        return True
    else:
        return False

# Compile the data.

compiled_file = ('%sIMERG/compiled.%s.%4d%02d-%4d%02d.%s.npy' % 
                 (datapath, ver, year0, month0, year1, month1, season))

P = np.zeros([48, nlon, nlat], 'f4')
N = np.zeros([48, nlon, nlat], 'i4')

for year in range(year0, year1 + 1):
    for month in getMonths(year, year0, month0, year1, month1):

        if monthInSeason(month, season):

            filename = ('%sIMERG/IMERG.%s.%4d%02d.npz' % 
                        (datapath, ver, year, month))

            data = np.load(filename)
            P += data['p'][:]
            N += data['n'][:]

meanPrecip = np.divide(P, N, where = N != 0, out = np.full(P.shape, np.nan))

np.save(compiled_file, meanPrecip)
