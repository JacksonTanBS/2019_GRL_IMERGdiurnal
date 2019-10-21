#!/usr/bin/python

'''This script defines the function for the plotting script.'''

import numpy as np
import h5py
import os
from glob import glob

datapath = '/path/to/processed/data/'

def getGrid():

    imergpath = '/path/to/imerg/'
    files = sorted(glob('%s2014/06/01/*.HDF5' % imergpath))

    with h5py.File(files[0], 'r') as f:
        lats = f['Grid/lat'][:]
        lons = f['Grid/lon'][:]
        fv = f['Grid/precipitationCal'].attrs['_FillValue']

    nlon, nlat = len(lons), len(lats)

    lonedges = np.linspace(-180, 180, nlon + 1)
    latedges = np.linspace(-90, 90, nlat + 1)

    return nlon, nlat, lons, lats, lonedges, latedges

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

def getData(year0, month0, year1, month1, lon0, lon1, lat0, lat1, season, ver):

    compiled_file = ('%sIMERG/compiled.%s.%4d%02d-%4d%02d.%s.npy' % 
                     (datapath, ver, year0, month0, year1, month1, season))

    if os.path.exists(compiled_file):

        meanPrecip = np.load(compiled_file)

    else:

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

        meanPrecip = np.divide(P, N, where = N != 0, 
                               out = np.full(P.shape, np.nan))

        np.save(compiled_file, meanPrecip)

    return meanPrecip[:, lon0 : lon1, lat0 : lat1]

def mapToLST(P0, lons, lon0, lon1):

    utc2lst = np.int32(np.round(lons[lon0 : lon1] / 7.5))    # 7.5° / hhr
    P1 = np.zeros(P0.shape, 'f4')

    for lon in range(lon1 - lon0):
        for hhr in range(48):
            P1[(hhr + utc2lst[lon]) % 48, lon, :] = P0[hhr, lon, :]

    return P1

def getMRMSData(year0, month0, year1, month1, lonl, lonr, latb, latt, 
                season, ver):

    from calendar import monthrange

    mrmspath = '/path/to/MRMS/'

    files = sorted(glob('%s2015/01/01/*' % mrmspath))

    with h5py.File(files[0], 'r') as f:
        lats_mrms = f['lat'][:]
        lons_mrms = f['lon'][:]
        fv = f['MRMS/meanPrecip4'].attrs['_FillValue']

    lat0 = np.where(lats_mrms == latb)[0][0]
    lat1 = np.where(lats_mrms == latt)[0][0] + 1
    lon0 = np.where(lons_mrms == lonl)[0][0]
    lon1 = np.where(lons_mrms == lonr)[0][0] + 1

    P = np.zeros([48, lon1 - lon0, lat1 - lat0], 'f4')
    N = np.zeros([48, lon1 - lon0, lat1 - lat0], 'i4')

    for year in range(year0, year1 + 1):
        for month in getMonths(year, year0, month0, year1, month1):

            if monthInSeason(month, season):

                for day in range(1, monthrange(year, month)[1] + 1):

                    files = sorted(glob('%s%4d/%02d/%02d/*' % 
                                        (mrmspath, year, month, day)))
    
                    meanPrecip = np.zeros([48, lon1 - lon0, lat1 - lat0], 
                                          dtype = 'f4')

                    for ts in range(48):
                        with h5py.File(files[ts], 'r') as f:
                            meanPrecip[ts] = f['MRMS/meanPrecip4'][lon0 : lon1, 
                                                                   lat0 : lat1]

                    N += (meanPrecip != fv)
                    meanPrecip[meanPrecip == fv] = 0
                    P += meanPrecip

    Pm = np.divide(P, N, where = N != 0, out = np.full(P.shape, np.nan), 
                   dtype = 'f4')

    return Pm

def fft_fit(y, cutoff0 = 0, cutoff1 = 5):

    from scipy import fftpack

    w = fftpack.rfft(y)
    w[:cutoff0] = 0
    w[cutoff1:] = 0

    return fftpack.irfft(w)

def getIRData(year0, month0, year1, month1, lon0, lon1, lat0, lat1, season):

    from calendar import monthrange

    imergpath = '/path/to/imerg/'

    files = sorted(glob('%s2015/01/01/*' % imergpath))

    with h5py.File(files[0], 'r') as f:
        fv = f['Grid/IRprecipitation'].attrs['_FillValue']

    P = np.zeros([48, lon1 - lon0, lat1 - lat0], 'f4')
    N = np.zeros([48, lon1 - lon0, lat1 - lat0], 'i4')

    for year in range(year0, year1 + 1):
        for month in getMonths(year, year0, month0, year1, month1):

            if monthInSeason(month, season):

                for day in range(1, monthrange(year, month)[1] + 1):

                    files = sorted(glob('%s%4d/%02d/%02d/*' % 
                                        (imergpath, year, month, day)))
    
                    meanPrecip = np.zeros([48, lon1 - lon0, lat1 - lat0], 
                                          dtype = 'f4')

                    for ts in range(48):
                        with h5py.File(files[ts], 'r') as f:
                            meanPrecip[ts] = \
                        f['Grid/IRprecipitation'][0, lon0 : lon1, lat0 : lat1]

                    N += (meanPrecip != fv)
                    meanPrecip[meanPrecip == fv] = 0
                    P += meanPrecip

    Pi = np.divide(P, N, where = N != 0, out = np.full(P.shape, np.nan), 
                   dtype = 'f4')

    return Pi

def plotMapCoords(ax, lonl, lonr, latb, latt):

    import cartopy.crs as ccrs
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    import matplotlib.ticker as mticker

    gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels = True,
                      linewidth = 0)
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xlocator = mticker.FixedLocator([lonl - 0.05, lonr + 0.05])
    gl.ylocator = mticker.FixedLocator([latb - 0.05, latt + 0.05])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    return None

def cmap_add_white(cmap, N = None):
    from matplotlib.colors import LinearSegmentedColormap
    if N == None: N = cmap.N
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmaplist.insert(0, (1.0, 1.0, 1.0, 1.0))
    return LinearSegmentedColormap.from_list('mcm', cmaplist, N)
