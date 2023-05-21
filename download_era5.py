#!/usr/bin/env python3

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "26/04/2019"
__description__ = "This script downloads the ERA5  dataset in pressure levels"


import os
import cdsapi
import calendar

cal = calendar.Calendar( )
apiclient = cdsapi.Client( )

ys = 1979
ye = 2018

vname = { 'relative_humidity' : 'rhum',
          'specific_humidity' : 'qhum',
          'u_component_of_wind' : 'uwnd',
          'v_component_of_wind' : 'vwnd' }

for year in range(ys,ye+1):
    yy = '%04d' % year
    try:
        os.mkdir(yy)
    except OSError:
        pass
    for month in range(4,6):
        mm = '%02d' % month
        dlist = [ ]
        for d in cal.itermonthdays(year, month):
            if d > 0:
                dlist.append('%02d' % d)
        for var in ( 'relative_humidity','specific_humidity','u_component_of_wind',
            'v_component_of_wind' ):
            netcdf = os.path.join(str(year),(vname[var]+"_"+yy+'_'+mm+'.nc'))
            if not os.path.isfile(netcdf):
                apiclient.retrieve(
    'reanalysis-era5-pressure-levels',
    {
        'product_type':'reanalysis',
        'format':'netcdf',
        'variable':[ var ],
        'pressure_level':[
            '1','2','3',
            '5','7','10',
            '20','30','50',
            '70','100','125',
            '150','175','200',
            '225','250','300',
            '350','400','450',
            '500','550','600',
            '650','700','750',
            '775','800','825',
            '850','875','900',
            '925','950','975',
            '1000'
        ],
        'year': yy,
        'month': mm,
        'day': dlist,
        'time':[ '00:00','06:00','12:00','18:00' ]
    },
    netcdf)
            else:
                print('File '+netcdf+' already on disk.')
