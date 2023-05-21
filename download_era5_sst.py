#!/usr/bin/env python3

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "26/04/2019"
__description__ = "This script downloads the ERA5  dataset in surface"


import os
import cdsapi
import calendar

cal = calendar.Calendar( )
apiclient = cdsapi.Client( )

ys = 1979
ye = 2018

vname = { 'sea_surface_temperature' : 'sst' }

for year in range(ys,ye+1):
    yy = '%04d' % year
    try:
        os.mkdir('SST')
    except OSError:
        pass
    for month in range(1,13):
        mm = '%02d' % month
        dlist = [ ]
        for d in cal.itermonthdays(year, month):
            if d > 0:
                dlist.append('%02d' % d)
        var = 'sea_surface_temperature'
        netcdf = os.path.join('SST',(vname[var]+"_"+yy+'_'+mm+'.nc'))
        if not os.path.isfile(netcdf):
            apiclient.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type':'reanalysis',
        'format':'netcdf',
        'variable':[ var ],
        'year': yy,
        'month': mm,
        'day': dlist,
        'time':[ '00:00','06:00','12:00','18:00' ]
    },
    netcdf)
        else:
            print('File '+netcdf+' already on disk.')
