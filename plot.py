#!/usr/bin/python

'''This script plots the figures for the paper.'''

import numpy as np
import matplotlib.pyplot as plt
import sys
import h5py
import gzip
import cartopy.crs as ccrs
import cartopy
import shapefile
from glob import glob
from calendar import monthrange
from func import *

options = sys.argv[1:]
figtype = 'pdf'

# plot configurations
scol = 3.503    # single column (89 mm)
dcol = 7.204    # double column (183 mm)
flpg = 9.724    # full page length (247 mm)
plt.rcParams['figure.figsize'] = (scol, 0.75 * scol)
plt.rcParams['font.size'] = 9
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['axes.titlesize'] = 'medium'
plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['font.sans-serif'] = ['TeX Gyre Heros', 'Helvetica',
                                   'DejaVu Sans']
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['image.cmap'] = 'viridis'
subplttxt = ('(a)', '(b)', '(c)', '(d)', '(e)', '(f)')

mrmspath = '/path/to/mrms/'
datapath = '/path/to/processed/data/'


if 'V05vsV06' in options:

    year0, year1 = 2014, 2015
    month0, month1 = 12, 2
    lonl, lonr, latb, latt = -179.95, 179.95, -59.95, -40.05
    imergpath = '/path/to/imerg/'

    nlon, nlat, lons, lats, _, _ = getGrid()

    lon0 = np.where(lons == lonl)[0][0]
    lon1 = np.where(lons == lonr)[0][0] + 1
    lat0 = np.where(lats == latb)[0][0]
    lat1 = np.where(lats == latt)[0][0] + 1

    Pi5 = mapToLST(getData(year0, month0, year1, month1, 
                           lon0, lon1, lat0, lat1, 'DJF', 'V05B'),
                   lons, lon0, lon1)
    Pi6 = mapToLST(getData(year0, month0, year1, month1, 
                           lon0, lon1, lat0, lat1, 'DJF', 'V06B'),
                   lons, lon0, lon1)

    # count the `HQprecipSource` for each half-hour.

    utc2lst = np.int32(np.round(lons[lon0 : lon1] / 7.5))
    HQcount = np.zeros([48, 12], 'i8')
    for year in range(year0, year1 + 1):
        for month in getMonths(year, year0, month0, year1, month1):
            if monthInSeason(month, 'DJF'):
                for day in range(1, monthrange(year, month)[1] + 1):

                    files = sorted(glob('%s%d/%02d/%02d/*.HDF5' % 
                                   (imergpath, year, month, day)))
                    for hhr in range(48):
                        with h5py.File(files[hhr], 'r') as f:
                            HQ = f['Grid/HQprecipSource'][0, lon0 : lon1, 
                                                             lat0 : lat1]
                        for lon in range(lon1 - lon0):
                            HQcount[(hhr + utc2lst[lon]) % 48] += \
                                np.bincount(HQ[lon], minlength = 12)

    sensors = np.where(np.sum(HQcount, 0))[0][1:]

    # plot the diurnal cycle

    cols = {3: 'C0', 5: 'C1', 7: 'C2', 11: 'C3'}
    sensorlabel = {3: 'AMSR2', 5: 'SSMIS', 7: 'MHS', 11: 'ATMS'}

    T = np.arange(0.25, 24, 0.5)

    plt.figure(figsize = (scol, scol))
    ax1 = plt.subplot2grid((3, 1), (0, 0))
    for sensor in sensors:
        if sensor == 9: continue    # skip GMI
        ax1.plot(T, HQcount[:, sensor], c = cols[sensor], lw = 1, 
                 label = sensorlabel[sensor])
    ax1.set_xticks(np.linspace(0, 24, 5))
    ax1.set_xticklabels(())
    ax1.set_ylim([-2e6, 4e7 + 2e6])
    ax1.set_yticks((0, 2e7, 4e7))
    ax1.set_yticklabels(('0', '2e7', '4e7'))
    ax1.set_ylabel('no. PMW obs.')
    ax1.legend(bbox_to_anchor = (0., 1.02, 1., .102), loc = 'lower left',
               ncol = 4, mode = 'expand', borderaxespad = 0.)
    ax1.grid()

    ax2 = plt.subplot2grid((3, 1), (1, 0), rowspan = 2)
    ax2.plot(T, np.nanmean(Pi5, (1, 2)), c = 'k', ls = ':', label = 'V05B')
    ax2.plot(T, np.nanmean(Pi6, (1, 2)), c = 'k', ls = '-', label = 'V06B')
    ax2.set_xlabel('LST hour')
    ax2.set_xticks(np.linspace(0, 24, 5))
    ax2.set_ylim([0.096, 0.184])
    ax2.set_ylabel('mean precipitation rate (mm / h)')
    ax2.legend()
    ax2.grid()

    plt.savefig('fig.V05vsV06.%s' % figtype)
    plt.close()


if 'CONUS' in options:

    year0, year1 = 2014, 2018
    month0, month1 = 1, 12
    lonl, lonr, latb, latt = -124.95, -65.05, 25.05, 49.95

    states = ('Louisiana', 'Mississippi', 'Alabama', 'Georgia', 'Florida',
              'South Carolina', 'North Carolina')
    cols = ['C%d' % ii for ii in range(len(states))]

    nlon, nlat, lons, lats, lonedges, latedges = getGrid()

    lon0 = np.where(lons == lonl)[0][0]
    lon1 = np.where(lons == lonr)[0][0] + 1
    lat0 = np.where(lats == latb)[0][0]
    lat1 = np.where(lats == latt)[0][0] + 1

    Pi = mapToLST(getData(year0, month0, year1, month1, 
                          lon0, lon1, lat0, lat1, 'JJA', 'V06B'),
                  lons, lon0, lon1)

    Pm = mapToLST(getMRMSData(year0, month0, year1, month1, 
                              lonl, lonr, latb, latt, 'JJA', 'V06B'),
                  lons, lon0, lon1)

    # find the diurnal peak

    peaki = np.full([lon1 - lon0, lat1 - lat0], np.nan, 'f4')
    peakm = np.full([lon1 - lon0, lat1 - lat0], np.nan, 'f4')
    for lon in range(lon1 - lon0):
        for lat in range(lat1 - lat0):
            peaki[lon, lat] = np.nanargmax(fft_fit(Pi[:, lon, lat]), 0) / 2
            if not np.any(np.isnan(Pm[:, lon, lat]), 0):
                peakm[lon, lat] = np.nanargmax(fft_fit(Pm[:, lon, lat]), 0) / 2

    # pick out desired region

    with gzip.open('%sGrids/gridUSState.txt.gz' % datapath, 'rb') as f:
        gridShape = np.loadtxt(f, dtype = 'i4')
    shp = shapefile.Reader('%sShapefiles/cb_2017_us_state_500k' % datapath)
    shpRecords = shp.records()
    stateIndices = [[rr + 1 for rr, rec in enumerate(shpRecords) 
                     if rec[5] == state][0] for state in states]

    inRegion = np.zeros([lon1 - lon0, lat1 - lat0], np.bool)
    for ii, index in enumerate(stateIndices):
        inRegion[gridShape[lon0 : lon1, lat0 : lat1] == index] = True

    # compute the peak differences in desired region

    diffpeak = np.full([lon1 - lon0, lat1 - lat0], np.nan, 'f4')
    for lon, lat in zip(*np.where(inRegion)):
        if not np.isnan(peaki[lon, lat]) and not np.isnan(peakm[lon, lat]):
            diff = peaki[lon, lat] - peakm[lon, lat] 
        if   diff > 12:
            diffpeak[lon, lat] = 24 - diff
        elif diff < -12:
            diffpeak[lon, lat] = 24 + diff
        else:
            diffpeak[lon, lat] = diff

    # plot the diurnal peaks

    def set_map():
        ax.coastlines(resolution = '50m', color = 'C2', lw = 1.25)
        ax.add_feature(cartopy.feature.BORDERS, edgecolor = 'C2', lw = 1.0)
        ax.add_feature(cartopy.feature.STATES, edgecolor = 'C2', lw = 0.5);
        ax.contour(lons[lon0 : lon1], lats[lat0 : lat1], inRegion.T, (0.5, 1.5), 
                   colors = 'k', linewidths = 1.0)
        ax.set_extent([lonl - 0.1, lonr + 0.1, latb - 0.1, latt + 0.1], 
                      ccrs.PlateCarree())
        plotMapCoords(ax, lonl, lonr, latb, latt)
        return None

    T = np.arange(0, 24, 0.5)

    plt.figure(figsize = (dcol, 0.7 * dcol))
    plt.subplots_adjust(hspace = 0.3, wspace = 0.25)

    ax = plt.subplot(221, projection = ccrs.PlateCarree())
    ax.set_title(subplttxt[0], loc = 'left')
    mp = plt.pcolormesh(lonedges[lon0 : lon1 + 1], latedges[lat0 : lat1 + 1], 
                        peaki.T, vmin = 0, vmax = 24, rasterized = True,
                        cmap = plt.cm.twilight_shifted, 
                        transform = ccrs.PlateCarree())
    set_map()

    ax = plt.subplot(222, projection = ccrs.PlateCarree())
    ax.set_title(subplttxt[1], loc = 'left')
    mp = plt.pcolormesh(lonedges[lon0 : lon1 + 1], latedges[lat0 : lat1 + 1], 
                        peakm.T, vmin = 0, vmax = 24, rasterized = True,
                        cmap = plt.cm.twilight_shifted, 
                        transform = ccrs.PlateCarree())
    set_map()

    cb = plt.colorbar(cax = plt.axes([0.125, 0.55, 0.775, 0.02]), mappable = mp, 
                      orientation = 'horizontal', ticks = np.arange(0, 24.1, 6))
    cb.set_label('LST hour')

    plt.subplot(223)
    plt.title(subplttxt[2], loc = 'left')
    plt.hist(diffpeak.flatten(), np.arange(-12.25, 12.26, 0.5), color = 'k')
    plt.xticks(np.arange(-12, 12.01, 6))
    plt.xlabel('difference in diurnal peak (h)')
    plt.ylabel('no. of pixels')
    plt.grid()

    plt.subplot(224)
    plt.title(subplttxt[3], loc = 'left')
    plt.plot(T, np.nanmean(Pi[:, inRegion], 1), 'k-')
    plt.plot(T, np.nanmean(Pm[:, inRegion], 1), 'k--')
    plt.plot(T, fft_fit(np.nanmean(Pi[:, inRegion], 1)), c = '0.5', ls = '-')
    plt.plot(T, fft_fit(np.nanmean(Pm[:, inRegion], 1)), c = '0.5', ls = '--')
    plt.xticks(np.arange(0, 24.1, 6))
    plt.xlabel('LST hour')
    plt.ylabel('mean precipitation rate (mm / h)')
    plt.grid()

    plt.savefig('fig.CONUS.%s' % figtype)
    plt.close()

    with open('data.CONUS.txt', 'w') as f:
        f.write('mean = %f\n' % np.mean(diffpeak[~np.isnan(diffpeak)]))
        f.write('median = %f\n' % np.median(diffpeak[~np.isnan(diffpeak)]))
        f.write('p25 = %f\n' % np.percentile(diffpeak[~np.isnan(diffpeak)], 25))
        f.write('p75 = %f\n' % np.percentile(diffpeak[~np.isnan(diffpeak)], 75))


if 'MaritimeContinent' in options:

    year0, year1 = 2001, 2018
    month0, month1 = 1, 12
    lonl, lonr, latb, latt = 90.05, 159.95, -14.95, 24.95

    nlon, nlat, lons, lats, lonedges, latedges = getGrid()

    lon0 = np.where(lons == lonl)[0][0]
    lon1 = np.where(lons == lonr)[0][0] + 1
    lat0 = np.where(lats == latb)[0][0]
    lat1 = np.where(lats == latt)[0][0] + 1

    Pi = mapToLST(getData(year0, month0, year1, month1, 
                          lon0, lon1, lat0, lat1, 'MAM', 'V06B'),
                  lons, lon0, lon1)

    # find the diurnal peak

    peak = np.full([lon1 - lon0, lat1 - lat0], np.nan, 'f4')
    for lon in range(lon1 - lon0):
        for lat in range(lat1 - lat0):
            peak[lon, lat] = np.nanargmax(fft_fit(Pi[:, lon, lat]), 0) / 2

    # plot the diurnal peaks

    plt.figure(figsize = (dcol, 0.6 * dcol))

    ax = plt.axes(projection = ccrs.PlateCarree())
    mp = plt.pcolormesh(lonedges[lon0 : lon1 + 1], latedges[lat0 : lat1 + 1], 
                        peak.T, vmin = 0, vmax = 24, rasterized = True,
                        cmap = plt.cm.twilight_shifted, 
                        transform = ccrs.PlateCarree())
    ax.coastlines(resolution = '50m', color = 'C2', lw = 0.8)
    ax.add_feature(cartopy.feature.BORDERS, edgecolor = 'C2', lw = 0.5)
    ax.set_extent([lonl - 0.06, lonr + 0.06, latb - 0.06, latt + 0.06], 
                  ccrs.PlateCarree())
    plotMapCoords(ax, lonl, lonr, latb, latt)

    cb = plt.colorbar(cax = plt.axes([0.125, 0.02, 0.775, 0.02]), mappable = mp, 
                      orientation = 'horizontal', ticks = np.arange(0, 24.1, 6))
    cb.set_label('LST hour')

    plt.savefig('fig.MaritimeContinent.%s' % figtype)
    plt.close()


if 'DiurnalCycles' in options:

    year0, year1 = 2001, 2018
    month0, month1 = 1, 12

    seasons = ('DJF', 'MAM', 'JJA', 'SON')
    T = np.arange(0, 24, 0.5)
    cols = {'DJF': 'C0', 'MAM': 'C2', 'JJA': 'C3', 'SON': 'C1'}

    _, _, lons, _, _, _ = getGrid()

    # read the shape map

    with gzip.open('%sGrids/gridCountry.txt.gz' % datapath, 'rb') as f:
        gridShape1 = np.loadtxt(f, dtype = 'i4')
    with gzip.open('%sGrids/gridLake.txt.gz' % datapath, 'rb') as f:
        gridShape2 = np.loadtxt(f, dtype = 'i4')

    # plot the diurnal cycle

    titles = {0: 'Singapore', 1: 'Bangladesh', 2: 'Lake Victoria'}

    plt.figure(figsize = (dcol, 0.5 * scol))
    plt.subplots_adjust(wspace = 0.25)

    sbplts = {}
    sbplts[0] = plt.subplot(131)
    sbplts[1] = plt.subplot(132, sharex = sbplts[0], sharey = sbplts[0])
    sbplts[2] = plt.subplot(133, sharex = sbplts[0], sharey = sbplts[0])

    for season in seasons:
        P = mapToLST(getData(year0, month0, year1, month1, 
                             0, 3600, 0, 1800, season, 'V06B'),
                     lons, 0, 3600)
        sbplts[0].plot(T, np.nanmean(P[:, gridShape1 == 204], 1), 
                       cols[season], label = season)
        sbplts[1].plot(T, np.nanmean(P[:, gridShape1 == 22], 1), cols[season])
        sbplts[2].plot(T, np.nanmean(P[:, gridShape2 == 479], 1), cols[season])

    sbplts[0].set_ylim([-0.05, 1.05])
    sbplts[0].set_yticks(np.arange(0, 1.01, 0.25))
    sbplts[0].set_xticks(np.arange(0, 24.1, 6))
    sbplts[0].set_ylabel('mean precip. rate (mm / h)')
    sbplts[0].legend()
    for ss, sbplt in sbplts.items():
        sbplt.set_title('%s %s' % (subplttxt[ss], titles[ss]), loc = 'left')
        sbplt.set_xlabel('LST hour')
        sbplt.grid()

    plt.savefig('fig.DiurnalCycles.%s' % figtype)
    plt.close()


#--- Supplementary ---#


if 'IRprecip' in options:

    year0, year1 = 2014, 2018
    month0, month1 = 6, 8
    lonl, lonr, latb, latt = -124.95, -65.05, 25.05, 49.95

    states = ('Louisiana', 'Mississippi', 'Alabama', 'Georgia', 'Florida',
              'South Carolina', 'North Carolina')
    cols = ['C%d' % ii for ii in range(len(states))]

    nlon, nlat, lons, lats, lonedges, latedges = getGrid()

    lon0 = np.where(lons == lonl)[0][0]
    lon1 = np.where(lons == lonr)[0][0] + 1
    lat0 = np.where(lats == latb)[0][0]
    lat1 = np.where(lats == latt)[0][0] + 1

    Pi = mapToLST(getIRData(year0, month0, year1, month1, 
                            lon0, lon1, lat0, lat1, 'JJA'),
                  lons, lon0, lon1)

    Pm = mapToLST(getMRMSData(year0, month0, year1, month1, 
                              lonl, lonr, latb, latt, 'JJA', 'V06B'),
                  lons, lon0, lon1)

    # find the diurnal peak

    peaki = np.full([lon1 - lon0, lat1 - lat0], np.nan, 'f4')
    peakm = np.full([lon1 - lon0, lat1 - lat0], np.nan, 'f4')
    for lon in range(lon1 - lon0):
        for lat in range(lat1 - lat0):
            peaki[lon, lat] = np.nanargmax(fft_fit(Pi[:, lon, lat]), 0) / 2
            if not np.any(np.isnan(Pm[:, lon, lat]), 0):
                peakm[lon, lat] = np.nanargmax(fft_fit(Pm[:, lon, lat]), 0) / 2

    # pick out desired region

    with gzip.open('%sGrids/gridUSState.txt.gz' % datapath, 'rb') as f:
        gridShape = np.loadtxt(f, dtype = 'i4')
    shp = shapefile.Reader('%sShapefiles/cb_2017_us_state_500k' % datapath)
    shpRecords = shp.records()
    stateIndices = [[rr + 1 for rr, rec in enumerate(shpRecords) 
                     if rec[5] == state][0] for state in states]

    inRegion = np.zeros([lon1 - lon0, lat1 - lat0], np.bool)
    for ii, index in enumerate(stateIndices):
        inRegion[gridShape[lon0 : lon1, lat0 : lat1] == index] = True

    # compute the peak differences in desired region

    diffpeak = np.full([lon1 - lon0, lat1 - lat0], np.nan, 'f4')
    for lon, lat in zip(*np.where(inRegion)):
        if not np.isnan(peaki[lon, lat]) and not np.isnan(peakm[lon, lat]):
            diff = peaki[lon, lat] - peakm[lon, lat] 
        if   diff > 12:
            diffpeak[lon, lat] = 24 - diff
        elif diff < -12:
            diffpeak[lon, lat] = 24 + diff
        else:
            diffpeak[lon, lat] = diff

    # plot the diurnal peaks

    def set_map():
        ax.coastlines(resolution = '50m', color = 'C2', lw = 1.25)
        ax.add_feature(cartopy.feature.BORDERS, edgecolor = 'C2', lw = 1.0)
        ax.add_feature(cartopy.feature.STATES, edgecolor = 'C2', lw = 0.5);
        ax.contour(lons[lon0 : lon1], lats[lat0 : lat1], inRegion.T, (0.5, 1.5), 
                   colors = 'k', linewidths = 1.0)
        ax.set_extent([lonl - 0.1, lonr + 0.1, latb - 0.1, latt + 0.1], 
                      ccrs.PlateCarree())
        plotMapCoords(ax, lonl, lonr, latb, latt)
        return None

    T = np.arange(0, 24, 0.5)

    plt.figure(figsize = (dcol, 0.7 * dcol))
    plt.subplots_adjust(hspace = 0.3, wspace = 0.25)

    ax = plt.subplot(221, projection = ccrs.PlateCarree())
    ax.set_title(subplttxt[0], loc = 'left')
    mp = plt.pcolormesh(lonedges[lon0 : lon1 + 1], latedges[lat0 : lat1 + 1], 
                        peaki.T, vmin = 0, vmax = 24, rasterized = True,
                        cmap = plt.cm.twilight_shifted, 
                        transform = ccrs.PlateCarree())
    set_map()


    ax = plt.subplot(222, projection = ccrs.PlateCarree())
    ax.set_title(subplttxt[1], loc = 'left')
    mp = plt.pcolormesh(lonedges[lon0 : lon1 + 1], latedges[lat0 : lat1 + 1], 
                        peakm.T, vmin = 0, vmax = 24, rasterized = True,
                        cmap = plt.cm.twilight_shifted, 
                        transform = ccrs.PlateCarree())
    set_map()

    cb = plt.colorbar(cax = plt.axes([0.125, 0.55, 0.775, 0.02]), mappable = mp, 
                      orientation = 'horizontal', ticks = np.arange(0, 24.1, 6))
    cb.set_label('LST hour')

    plt.subplot(223)
    plt.title(subplttxt[2], loc = 'left')
    plt.hist(diffpeak.flatten(), np.arange(-12.25, 12.26, 0.5), color = 'k')
    plt.xticks(np.arange(-12, 12.01, 6))
    plt.xlabel('difference in diurnal peak (h)')
    plt.ylabel('no. of pixels')
    plt.grid()

    plt.subplot(224)
    plt.title(subplttxt[3], loc = 'left')
    plt.plot(T, np.nanmean(Pi[:, inRegion], 1), 'k-')
    plt.plot(T, np.nanmean(Pm[:, inRegion], 1), 'k--')
    plt.plot(T, fft_fit(np.nanmean(Pi[:, inRegion], 1)), c = '0.5', ls = '-')
    plt.plot(T, fft_fit(np.nanmean(Pm[:, inRegion], 1)), c = '0.5', ls = '--')
    plt.xticks(np.arange(0, 24.1, 6))
    plt.xlabel('LST hour')
    plt.ylabel('mean precipitation rate (mm / h)')
    plt.grid()

    plt.savefig('sfig.IRprecip.%s' % figtype)
    plt.close()

    with open('sdata.IRprecip.txt', 'w') as f:
        f.write('mean = %f\n' % np.mean(diffpeak[~np.isnan(diffpeak)]))
        f.write('median = %f\n' % np.median(diffpeak[~np.isnan(diffpeak)]))
        f.write('p25 = %f\n' % np.percentile(diffpeak[~np.isnan(diffpeak)], 25))
        f.write('p75 = %f\n' % np.percentile(diffpeak[~np.isnan(diffpeak)], 75))


if 'Maps' in options:

    import itertools
    from cartopy.io import shapereader
    from shapely.geometry import asShape
    from shapely.ops import unary_union

    _, _, lons, lats, lonedges, latedges = getGrid()

    # read the shape map

    with gzip.open('%sGrids/gridCountry.txt.gz' % datapath, 'rb') as f:
        gridShape1 = np.loadtxt(f, dtype = 'i4')
    with gzip.open('%sGrids/gridLake.txt.gz' % datapath, 'rb') as f:
        gridShape2 = np.loadtxt(f, dtype = 'i4')

    # define the setings for each region

    titles = {0: 'Singapore', 1: 'Bangladesh', 2: 'Lake Victoria'}
    files = {0: 'MP14_PLNG_AREA_NO_SEA_PL', 1: 'Countries_WGS84',
             2: 'ne_10m_lakes'}
    projs = {0: ccrs.TransverseMercator(central_longitude = 103.8333333333333, 
                                        central_latitude = 1.366666666666667, 
                                        false_easting = 28001.642, 
                                        false_northing = 38744.572),
             1: ccrs.PlateCarree(), 2: ccrs.PlateCarree()}   # proj. of shpfile
    indices = {0: 204, 1: 22, 2: 479}
    exts = {0: 3, 1: 20, 2: 20}
    grids = {0: gridShape1, 1: gridShape1, 2: gridShape2}
    cres = {0: '10m', 1: '50m', 2: '50m'}

    # plot each region

    for ii in range(3):

        if ii:    # standard treatment
            shp = shapereader.Reader('%sShapefiles/%s' % (datapath, files[ii]))
            shape = next(itertools.islice(shp.records(), 
                                          indices[ii] - 1, None)).geometry
        else:    # special handling for Singapore due to shapefile specs.
            shp = shapefile.Reader('%sShapefiles/%s' % (datapath, files[ii]))
            shape = unary_union([asShape(shape) for shape in shp.shapes()])

        # plot the diurnal cycle

        alllon, alllat = np.where(grids[ii] == indices[ii])
        lon0, lon1 = alllon.min() - exts[ii], alllon.max() + exts[ii]
        lat0, lat1 = alllat.min() - exts[ii], alllat.max() + exts[ii]
        lonl, lonr, latb, latt = lons[lon0], lons[lon1], lats[lat0], lats[lat1]

        plt.figure(figsize = (scol, 0.75 * scol))
        ax = plt.axes(projection = ccrs.PlateCarree())
        ax.set_title('%s (n = %d)' % (titles[ii], len(alllon)), loc = 'left')
        ax.set_extent([lonl - 0.1, lonr + 0.1, latb - 0.1, latt + 0.1], 
                      ccrs.PlateCarree())
        ax.add_feature(cartopy.feature.BORDERS.with_scale('50m'), 
                       edgecolor = '0.75', lw = 0.5, zorder = 1)
        ax.coastlines(resolution = cres[ii], color = '0.75', 
                      lw = 0.5, zorder = 1)
        ax.add_geometries(shape, projs[ii], facecolor = 'w', 
                          edgecolor = 'k', lw = 1.0, zorder = 2)
        ax.scatter([lons[lon] for lon in alllon], [lats[lat] for lat in alllat],
                   s = (75 / (lon1 - lon0)) ** 2, facecolor = 'k', 
                   edgecolor = 'none', zorder = 3, 
                   transform = ccrs.PlateCarree())

        plotMapCoords(ax, lonl, lonr, latb, latt)

        plt.savefig('sfig.Map.%s.%s' % (titles[ii].replace(' ', ''), figtype))
        plt.close()


if 'V05vsV06NEPac' in options:

    year0, year1 = 2014, 2014
    month0, month1 = 1, 12
    lonl, lonr, latb, latt = -179.95, -130.05, 30.05, 49.95
    imergpath = '/path/to/imerg/'

    nlon, nlat, lons, lats, _, _ = getGrid()

    lon0 = np.where(lons == lonl)[0][0]
    lon1 = np.where(lons == lonr)[0][0] + 1
    lat0 = np.where(lats == latb)[0][0]
    lat1 = np.where(lats == latt)[0][0] + 1

    Pi5 = mapToLST(getData(year0, month0, year1, month1, 
                           lon0, lon1, lat0, lat1, 'JJA', 'V05B'),
                   lons, lon0, lon1)
    Pi6 = mapToLST(getData(year0, month0, year1, month1, 
                           lon0, lon1, lat0, lat1, 'JJA', 'V06B'),
                   lons, lon0, lon1)

    # count the `HQprecipSource` for each half-hour.

    utc2lst = np.int32(np.round(lons[lon0 : lon1] / 7.5))
    HQcount = np.zeros([48, 12], 'i8')
    for year in range(year0, year1 + 1):
        for month in getMonths(year, year0, month0, year1, month1):
            if monthInSeason(month, 'JJA'):
                for day in range(1, monthrange(year, month)[1] + 1):

                    files = sorted(glob('%s%d/%02d/%02d/*.HDF5' % 
                                   (imergpath, year, month, day)))
                    for hhr in range(48):
                        with h5py.File(files[hhr], 'r') as f:
                            HQ = f['Grid/HQprecipSource'][0, lon0 : lon1, 
                                                             lat0 : lat1]
                        for lon in range(lon1 - lon0):
                            HQcount[(hhr + utc2lst[lon]) % 48] += \
                                np.bincount(HQ[lon], minlength = 12)

    sensors = np.where(np.sum(HQcount, 0))[0][1:]

    # plot the diurnal cycle

    cols = {sensor: 'C%d' % ss for ss, sensor in enumerate(sensors[1:])}
    sensorlabel = {3: 'AMSR2', 5: 'SSMIS', 7: 'MHS', 11: 'ATMS'}

    T = np.arange(0.25, 24, 0.5)

    plt.figure(figsize = (scol, scol))
    ax1 = plt.subplot2grid((3, 1), (0, 0))
    for sensor in sensors:
        if sensor == 1 or sensor == 9: continue    # skip TMI & GMI
        ax1.plot(T, HQcount[:, sensor], c = cols[sensor], lw = 1, 
                 label = sensorlabel[sensor])
    ax1.set_xticks(np.linspace(0, 24, 5))
    ax1.set_xticklabels(())
    ax1.set_ylim([-3e5, 6e6 + 3e5])
    ax1.set_yticks((0, 3e6, 6e6))
    ax1.set_yticklabels(('0', '3e7', '6e7'))
    ax1.set_ylabel('no. PMW obs.')
    ax1.legend(bbox_to_anchor = (0., 1.02, 1., .102), loc = 'lower left',
               ncol = 4, mode = 'expand', borderaxespad = 0.)
    ax1.grid()

    ax2 = plt.subplot2grid((3, 1), (1, 0), rowspan = 2)
    ax2.plot(T, np.nanmean(Pi5, (1, 2)), c = 'k', ls = ':', label = 'V05B')
    ax2.plot(T, np.nanmean(Pi6, (1, 2)), c = 'k', ls = '-', label = 'V06B')
    ax2.set_xlabel('LST hour')
    ax2.set_xticks(np.linspace(0, 24, 5))
    ax2.set_ylabel('mean precipitation rate (mm / h)')
    ax2.legend()
    ax2.grid()

    plt.savefig('sfig.V05vsV06NEPac.%s' % figtype)
    plt.close()


if 'CONUSGPstates' in options:

    year0, year1 = 2014, 2018
    month0, month1 = 1, 12
    lonl, lonr, latb, latt = -124.95, -65.05, 25.05, 49.95

    states = ('Oklahoma', 'Kansas', 'Nebraska', 'Colorado')
    cols = ['C%d' % ii for ii in range(len(states))]

    nlon, nlat, lons, lats, lonedges, latedges = getGrid()

    lon0 = np.where(lons == lonl)[0][0]
    lon1 = np.where(lons == lonr)[0][0] + 1
    lat0 = np.where(lats == latb)[0][0]
    lat1 = np.where(lats == latt)[0][0] + 1

    Pi = mapToLST(getData(year0, month0, year1, month1, 
                          lon0, lon1, lat0, lat1, 'JJA', 'V06B'),
                  lons, lon0, lon1)

    Pm = mapToLST(getMRMSData(year0, month0, year1, month1, 
                              lonl, lonr, latb, latt, 'JJA', 'V06B'),
                  lons, lon0, lon1)

    # find the diurnal peak

    peaki = np.full([lon1 - lon0, lat1 - lat0], np.nan, 'f4')
    peakm = np.full([lon1 - lon0, lat1 - lat0], np.nan, 'f4')
    for lon in range(lon1 - lon0):
        for lat in range(lat1 - lat0):
            peaki[lon, lat] = np.nanargmax(fft_fit(Pi[:, lon, lat]), 0) / 2
            if not np.any(np.isnan(Pm[:, lon, lat]), 0):
                peakm[lon, lat] = np.nanargmax(fft_fit(Pm[:, lon, lat]), 0) / 2

    # pick out desired region

    with gzip.open('%sGrids/gridUSState.txt.gz' % datapath, 'rb') as f:
        gridShape = np.loadtxt(f, dtype = 'i4')
    shp = shapefile.Reader('%sShapefiles/cb_2017_us_state_500k' % datapath)
    shpRecords = shp.records()
    stateIndices = [[rr + 1 for rr, rec in enumerate(shpRecords) 
                     if rec[5] == state][0] for state in states]

    inRegion = np.zeros([lon1 - lon0, lat1 - lat0], np.bool)
    for ii, index in enumerate(stateIndices):
        inRegion[gridShape[lon0 : lon1, lat0 : lat1] == index] = True

    # compute the peak differences in desired region

    diffpeak = np.full([lon1 - lon0, lat1 - lat0], np.nan, 'f4')
    for lon, lat in zip(*np.where(inRegion)):
        if not np.isnan(peaki[lon, lat]) and not np.isnan(peakm[lon, lat]):
            diff = peaki[lon, lat] - peakm[lon, lat] 
        if   diff > 12:
            diffpeak[lon, lat] = 24 - diff
        elif diff < -12:
            diffpeak[lon, lat] = 24 + diff
        else:
            diffpeak[lon, lat] = diff

    # plot the diurnal peaks

    def set_map():
        ax.coastlines(resolution = '50m', color = 'C2', lw = 1.25)
        ax.add_feature(cartopy.feature.BORDERS, edgecolor = 'C2', lw = 1.0)
        ax.add_feature(cartopy.feature.STATES, edgecolor = 'C2', lw = 0.5);
        ax.contour(lons[lon0 : lon1], lats[lat0 : lat1], inRegion.T, (0.5, 1.5), 
                   colors = 'k', linewidths = 1.0)
        ax.set_extent([lonl - 0.1, lonr + 0.1, latb - 0.1, latt + 0.1], 
                      ccrs.PlateCarree())
        plotMapCoords(ax, lonl, lonr, latb, latt)
        return None

    T = np.arange(0, 24, 0.5)

    plt.figure(figsize = (dcol, 0.7 * dcol))
    plt.subplots_adjust(hspace = 0.3, wspace = 0.25)

    ax = plt.subplot(221, projection = ccrs.PlateCarree())
    ax.set_title(subplttxt[0], loc = 'left')
    mp = plt.pcolormesh(lonedges[lon0 : lon1 + 1], latedges[lat0 : lat1 + 1], 
                        peaki.T, vmin = 0, vmax = 24, rasterized = True,
                        cmap = plt.cm.twilight_shifted, 
                        transform = ccrs.PlateCarree())
    set_map()

    ax = plt.subplot(222, projection = ccrs.PlateCarree())
    ax.set_title(subplttxt[1], loc = 'left')
    mp = plt.pcolormesh(lonedges[lon0 : lon1 + 1], latedges[lat0 : lat1 + 1], 
                        peakm.T, vmin = 0, vmax = 24, rasterized = True,
                        cmap = plt.cm.twilight_shifted, 
                        transform = ccrs.PlateCarree())
    set_map()

    cb = plt.colorbar(cax = plt.axes([0.125, 0.55, 0.775, 0.02]), mappable = mp, 
                      orientation = 'horizontal', ticks = np.arange(0, 24.1, 6))
    cb.set_label('LST hour')

    plt.subplot(223)
    plt.title(subplttxt[2], loc = 'left')
    plt.hist(diffpeak.flatten(), np.arange(-12.25, 12.26, 0.5), color = 'k')
    plt.xticks(np.arange(-12, 12.01, 6))
    plt.xlabel('difference in diurnal peak (h)')
    plt.ylabel('no. of pixels')
    plt.grid()

    plt.subplot(224)
    plt.title(subplttxt[3], loc = 'left')
    plt.plot(T, np.nanmean(Pi[:, inRegion], 1), 'k-')
    plt.plot(T, np.nanmean(Pm[:, inRegion], 1), 'k--')
    plt.plot(T, fft_fit(np.nanmean(Pi[:, inRegion], 1)), c = '0.5', ls = '-')
    plt.plot(T, fft_fit(np.nanmean(Pm[:, inRegion], 1)), c = '0.5', ls = '--')
    plt.xticks(np.arange(0, 24.1, 6))
    plt.xlabel('LST hour')
    plt.ylabel('mean precipitation rate (mm / h)')
    plt.grid()

    plt.savefig('sfig.CONUSGPstates.%s' % figtype)
    plt.close()

    with open('sdata.CONUSGPstates.txt', 'w') as f:
        f.write('mean = %f\n' % np.mean(diffpeak[~np.isnan(diffpeak)]))
        f.write('median = %f\n' % np.median(diffpeak[~np.isnan(diffpeak)]))
        f.write('p25 = %f\n' % np.percentile(diffpeak[~np.isnan(diffpeak)], 25))
        f.write('p75 = %f\n' % np.percentile(diffpeak[~np.isnan(diffpeak)], 75))


if 'CONUSamplitude' in options:

    year0, year1 = 2014, 2018
    month0, month1 = 1, 12
    lonl, lonr, latb, latt = -124.95, -65.05, 25.05, 49.95

    states = ('Louisiana', 'Mississippi', 'Alabama', 'Georgia', 'Florida',
              'South Carolina', 'North Carolina')
    cols = ['C%d' % ii for ii in range(len(states))]

    nlon, nlat, lons, lats, lonedges, latedges = getGrid()

    lon0 = np.where(lons == lonl)[0][0]
    lon1 = np.where(lons == lonr)[0][0] + 1
    lat0 = np.where(lats == latb)[0][0]
    lat1 = np.where(lats == latt)[0][0] + 1

    Pi = mapToLST(getData(year0, month0, year1, month1, 
                          lon0, lon1, lat0, lat1, 'JJA', 'V06B'),
                  lons, lon0, lon1)

    Pm = mapToLST(getMRMSData(year0, month0, year1, month1, 
                              lonl, lonr, latb, latt, 'JJA', 'V06B'),
                  lons, lon0, lon1)

    # calculate the amplitudes

    Ai = np.full([lon1 - lon0, lat1 - lat0], np.nan, 'f4')
    Am = np.full([lon1 - lon0, lat1 - lat0], np.nan, 'f4')
    for lon in range(lon1 - lon0):
        for lat in range(lat1 - lat0):
            p = fft_fit(Pi[:, lon, lat])
            Ai[lon, lat] = (np.max(p) - np.min(p)) / np.mean(p)
            if not np.any(np.isnan(Pm[:, lon, lat]), 0):
                p = fft_fit(Pm[:, lon, lat])
                Am[lon, lat] = (np.max(p) - np.min(p)) / np.mean(p)

    # plot the amplitudes

    def set_map():
        ax.coastlines(resolution = '50m', color = 'C2', lw = 1.25)
        ax.add_feature(cartopy.feature.BORDERS, edgecolor = 'C2', lw = 1.0)
        ax.add_feature(cartopy.feature.STATES, edgecolor = 'C2', lw = 0.5);
        ax.set_extent([lonl - 0.1, lonr + 0.1, latb - 0.1, latt + 0.1], 
                      ccrs.PlateCarree())
        plotMapCoords(ax, lonl, lonr, latb, latt)
        return None

    plt.figure(figsize = (dcol, 0.35 * dcol))
    plt.subplots_adjust(hspace = 0.3, wspace = 0.25)

    ax = plt.subplot(121, projection = ccrs.PlateCarree())
    ax.set_title(subplttxt[0], loc = 'left')
    mp = plt.pcolormesh(lonedges[lon0 : lon1 + 1], latedges[lat0 : lat1 + 1], 
                        Ai.T, vmin = 0, vmax = 4, rasterized = True,
                        cmap = plt.cm.inferno_r, 
                        transform = ccrs.PlateCarree())
    set_map()

    ax = plt.subplot(122, projection = ccrs.PlateCarree())
    ax.set_title(subplttxt[1], loc = 'left')
    mp = plt.pcolormesh(lonedges[lon0 : lon1 + 1], latedges[lat0 : lat1 + 1], 
                        Am.T, vmin = 0, vmax = 4, rasterized = True,
                        cmap = plt.cm.inferno_r, 
                        transform = ccrs.PlateCarree())
    set_map()

    cb = plt.colorbar(cax = plt.axes([0.125, 0.1, 0.775, 0.04]), mappable = mp, 
                      orientation = 'horizontal', extend = 'max')
    cb.set_label('(max − min) / mean')

    plt.savefig('sfig.CONUSamplitude.%s' % figtype)
    plt.close()


if 'BiennialChange' in options:

    year0, year1 = 2001, 2018
    month0, month1 = 1, 12
    lonl, lonr, latb, latt = -124.95, -65.05, 25.05, 49.95

    nlon, nlat, lons, lats, lonedges, latedges = getGrid()

    lon0 = np.where(lons == lonl)[0][0]
    lon1 = np.where(lons == lonr)[0][0] + 1
    lat0 = np.where(lats == latb)[0][0]
    lat1 = np.where(lats == latt)[0][0] + 1

    # pick out desired region

    states = ('Louisiana', 'Mississippi', 'Alabama', 'Georgia', 'Florida',
              'South Carolina', 'North Carolina')

    with gzip.open('%sGrids/gridUSState.txt.gz' % datapath, 'rb') as f:
        gridShape = np.loadtxt(f, dtype = 'i4')
    shp = shapefile.Reader('%sShapefiles/cb_2017_us_state_500k' % datapath)
    shpRecords = shp.records()
    stateIndices = [[rr + 1 for rr, rec in enumerate(shpRecords) 
                     if rec[5] == state][0] for state in states]

    inRegion = np.zeros([lon1 - lon0, lat1 - lat0], np.bool)
    for ii, index in enumerate(stateIndices):
        inRegion[gridShape[lon0 : lon1, lat0 : lat1] == index] = True

    # get the diurnal cycle of the biennial period

    P = {}
    for year in range(year0, year1 + 1, 2):
        Pi = mapToLST(getData(year, month0, year + 1, month1, 
                              lon0, lon1, lat0, lat1, 'JJA', 'V06B'),
                      lons, lon0, lon1)
        P[year] = fft_fit(np.nanmean(Pi[:, inRegion], 1))


    T = np.arange(0, 24, 0.5)

    plt.figure(figsize = (scol, 0.75 * scol))
    for yy, year in enumerate(range(year0, year1 + 1, 2)):
        plt.plot(T, P[year], c = 'C%d' % yy, ls = '-', lw = 1.0, 
                 label = '%s–%s' % (year, year + 1))
    plt.xticks(np.arange(0, 24.1, 6))
    plt.xlabel('LST hour')
    plt.ylabel('mean precipitation rate (mm / h)')
    plt.legend(ncol = 1)
    plt.grid()

    plt.savefig('sfig.BiennialChange.%s' % figtype)
    plt.close()


if 'LakeVictoria' in options:

    import itertools
    from cartopy.io import shapereader

    year0, year1 = 2001, 2018
    month0, month1 = 1, 12
    lonl, lonr, latb, latt = 29.65, 36.75, -4.95, 2.45

    nlon, nlat, lons, lats, lonedges, latedges = getGrid()

    lon0 = np.where(lons == lonl)[0][0]
    lon1 = np.where(lons == lonr)[0][0] + 1
    lat0 = np.where(lats == latb)[0][0]
    lat1 = np.where(lats == latt)[0][0] + 1

    # find the diurnal peak

    Pi = mapToLST(getData(year0, month0, year1, month1, 
                          lon0, lon1, lat0, lat1, 'MAM', 'V06B'),
                  lons, lon0, lon1)

    peak = np.full([lon1 - lon0, lat1 - lat0], np.nan, 'f4')
    for lon in range(lon1 - lon0):
        for lat in range(lat1 - lat0):
            peak[lon, lat] = np.nanargmax(fft_fit(Pi[:, lon, lat]), 0) / 2

    # get the shape of Lake Victoria

    shp = shapereader.Reader('%sShapefiles/ne_10m_lakes' % (datapath,))
    shape = next(itertools.islice(shp.records(), 479 - 1, None)).geometry

    # plot the diurnal peaks

    plt.figure(figsize = (scol, 1.1 * scol))

    ax = plt.axes(projection = ccrs.PlateCarree())
    mp = plt.pcolormesh(lonedges[lon0 : lon1 + 1], latedges[lat0 : lat1 + 1], 
                        peak.T, vmin = 0, vmax = 24, rasterized = True,
                        cmap = plt.cm.twilight_shifted, 
                        transform = ccrs.PlateCarree())
    ax.coastlines(resolution = '10m', color = 'C2', lw = 0.8)
    ax.add_feature(cartopy.feature.BORDERS.with_scale('50m'), 
                   edgecolor = 'C2', lw = 0.5)
    ax.add_geometries(shape, ccrs.PlateCarree(), facecolor = 'none', 
                      edgecolor = 'g', lw = 1.0, zorder = 2)
    ax.set_extent([lonl - 0.06, lonr + 0.06, latb - 0.06, latt + 0.06], 
                  ccrs.PlateCarree())
    plotMapCoords(ax, lonl, lonr, latb, latt)

    cb = plt.colorbar(cax = plt.axes([0.125, 0.02, 0.775, 0.02]), mappable = mp, 
                      orientation = 'horizontal', ticks = np.arange(0, 24.1, 6))
    cb.set_label('LST hour')

    plt.savefig('sfig.LakeVictoria.%s' % figtype)
    plt.close()


#--- Not used ---#


if 'GlobalPeak' in options:

    year0, year1 = 2001, 2018
    month0, month1 = 1, 12
    lonl, lonr, latb, latt = -179.95, 179.95, -89.95, 89.95

    nlon, nlat, lons, lats, lonedges, latedges = getGrid()

    lon0 = np.where(lons == lonl)[0][0]
    lon1 = np.where(lons == lonr)[0][0] + 1
    lat0 = np.where(lats == latb)[0][0]
    lat1 = np.where(lats == latt)[0][0] + 1

    for season in ('JJA', 'DJF', 'MAM', 'SON'):

        Pi = mapToLST(getData(year0, month0, year1, month1, 
                              lon0, lon1, lat0, lat1, season, 'V06B'),
                      lons, lon0, lon1)

        # find the diurnal peak

        peak = np.full([lon1 - lon0, lat1 - lat0], np.nan, 'f4')
        for lon in range(lon1 - lon0):
            for lat in range(lat1 - lat0):
                if np.all(~np.isnan(Pi[:, lon, lat])):
                    peak[lon, lat] = (np.nanargmax(fft_fit(Pi[:, lon, lat]), 0) 
                                      / 2)

        # plot the diurnal peaks

        plt.figure(figsize = (30, 17))

        ax = plt.axes([0.1, 0.14, 0.8, 0.77], projection = ccrs.Mollweide())
        ax.set_title(season, loc = 'left', size = 25)
        mp = plt.pcolormesh(lonedges[lon0 : lon1 + 1], 
                            latedges[lat0 : lat1 + 1], 
                            peak.T, vmin = 0, vmax = 24, rasterized = True,
                            cmap = plt.cm.twilight_shifted, 
                            transform = ccrs.PlateCarree())
        ax.coastlines(resolution = '50m', color = 'C2', lw = 0.6)
        ax.add_feature(cartopy.feature.BORDERS, edgecolor = 'C2', lw = 0.3)
        ax.set_global()

        cb = plt.colorbar(cax = plt.axes([0.1, 0.12, 0.8, 0.02]), 
                          mappable = mp,  orientation = 'horizontal', 
                          ticks = np.arange(0, 24.1, 3))
        cb.ax.tick_params(labelsize = 15)
        cb.set_label('LST hour', size = 15)

        plt.savefig('sfig.GlobalPeak.%s.%s' % (season, figtype), 
                    bbox_inches = '')
        plt.close()


if 'SEstates' in options:

    year0, year1 = 2014, 2018
    month0, month1 = 1, 12
    lonl, lonr, latb, latt = -124.95, -65.05, 25.05, 49.95

    states = ('Alabama', 'Georgia', 'Florida', 'South Carolina')
    cols = ['C%d' % ii for ii in range(len(states))]

    nlon, nlat, lons, lats, lonedges, latedges = getGrid()

    lon0 = np.where(lons == lonl)[0][0]
    lon1 = np.where(lons == lonr)[0][0] + 1
    lat0 = np.where(lats == latb)[0][0]
    lat1 = np.where(lats == latt)[0][0] + 1

    Pi = mapToLST(getData(year0, month0, year1, month1, 
                          lon0, lon1, lat0, lat1, 'JJA', 'V06B'),
                  lons, lon0, lon1)

    Pm = mapToLST(getMRMSData(year0, month0, year1, month1, 
                              lonl, lonr, latb, latt, 'JJA', 'V06B'),
                  lons, lon0, lon1)

    # pick out desired region

    with gzip.open('%sGrids/gridUSState.txt.gz' % datapath, 'rb') as f:
        gridShape = np.loadtxt(f, dtype = 'i4')
    shp = shapefile.Reader('%sShapefiles/cb_2017_us_state_500k' % datapath)
    shpRecords = shp.records()
    stateIndices = [[rr + 1 for rr, rec in enumerate(shpRecords) 
                     if rec[5] == state][0] for state in states]

    # plot the diurnal cycle

    T = np.arange(0, 24, 0.5)

    plt.figure(figsize = (scol, 0.75 * scol))

    for ii, index in enumerate(stateIndices):
        inRegion = gridShape[lon0 : lon1, lat0 : lat1] == index
        plt.plot(T, np.nanmean(Pi[:, inRegion], 1), 
                 c = cols[ii], label = shpRecords[index - 1][5])
        plt.plot(T, np.nanmean(Pm[:, inRegion], 1), 
                 c = cols[ii], ls = '--')

    plt.xticks(np.arange(0, 24.1, 6))
    plt.xlabel('LST hour')
    plt.ylabel('mean precipitation rate (mm / h)')
    plt.legend(loc = 2)
    plt.grid()
    plt.savefig('fig.SEstates.%s' % figtype)
    plt.close()
