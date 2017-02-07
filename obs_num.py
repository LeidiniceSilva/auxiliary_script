# -*- coding: utf-8 -*-
"""Comparison between Thiessen computed by different databases."""

__author__ = "Leidinice Silva"
__copyright__ = "Copyright 2016, Funceme Hydropy Project"
__credits__ = ["Jarbas Camurca", "Diogenes Fontenele"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Jarbas Camurca"
__email__ = "leidinice.silva@funceme.br"
__date__ = 9 / 18 / 2016

# Import datas
import netCDF4
import numpy as np
import pandas as pd
import os
import calendar
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
from pandas import Series, DataFrame
from hidropy.preprocessing.utils import basin_dict
from matplotlib import pyplot as plt

scale = 'daily'
param = 'pr'
period = 'calibration'
start_date = '20090101'
end_date = '20141231'
hidropy_path = "/home/leidinice/Documentos/musf"

home = os.path.expanduser("~")

folders = os.listdir("{0}/hidropy/hidropy/shapes/basins/".format(hidropy_path))
basins = sorted(basin_dict(micro=True, basin_name='tiete'))  # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< change here!!

for basin in basins:
    basin_fullname = basin_dict(basin)[2]
    macro_name = basin_dict(basin)[1]

    st1 = []
    st2 = []
    st3 = []
    st4 = []

    # open netcdf
    link1 = home+"/io/chirps/calibration/{0}/{1}_thiessen/{2}".format(scale, param, macro_name)
    arq1 = "{0}/{1}_{2}_chirps_obs_19810101_20141231_thiessen_{3}.nc".format(link1, param, scale, macro_name)
    data_chirps = netCDF4.Dataset(arq1)
    variable_chirps = data_chirps.variables[param][:].T
    time_chirps = data_chirps.variables['time']
    st1 = variable_chirps[10227:12418]

    link2 = home+"/io/inmet_ana/calibration/{0}/{1}_thiessen/{2}".format(scale, param, macro_name)
    arq2 = "{0}/{1}_{2}_inmet_ana_obs_19610101_20141231_thiessen_{3}.nc".format(link2, param, scale, macro_name)
    data_inmet_ana = netCDF4.Dataset(arq2)
    variable_inmet_ana = data_inmet_ana.variables[param][:].T
    time_inmet_ana = data_inmet_ana.variables['time']
    st2 = variable_inmet_ana[17533:19716]

    link3 = home+"/io/merge/calibration/{0}/{1}_thiessen/{2}".format(scale, param, macro_name)
    arq3 = "{0}/{1}_{2}_merge_obs_19980101_20141231_thiessen_{3}.nc".format(link3, param, scale, macro_name)
    data_merge = netCDF4.Dataset(arq3)
    variable_merge = data_merge.variables[param][:].T
    time_merge = data_merge.variables['time']
    st3 = variable_merge[4018:6209]

    stg = []

    # open fcst
    for year in range(2009, 2015):
        for month in range(1, 13):
            link4 = home+"/io/gfs05/{0:04d}/{0:04d}{1:02d}/pr_thiessen/{2}".format(year, month, macro_name)

            dias_mes = calendar.monthrange(year, month)[1]
            for daily in range(1, dias_mes + 1):
                if year == 2014 and month == 12 and daily == 24:
                    break
                init = np.datetime64('{0:04d}-{1:02d}-{2:02d}'.format(year, month, daily))
                start = init + np.timedelta64(1, 'D')
                end = init + np.timedelta64(7, 'D')

                init_y = str(init)[0:4] + str(init)[5:7] + str(init)[8:10]
                start_y = str(start)[0:4] + str(start)[5:7] + str(start)[8:10]
                end_y = str(end)[0:4] + str(end)[5:7] + str(end)[8:10]
                # print init_y, start_y, end_y

                arq4 = "{0}/pr_daily_gfs05_fcst_{1}00_{2}_{3}_thiessen_{4}.nc".format(link4, init_y, start_y, end_y,
                                                                                      macro_name)
                data_gfs05 = netCDF4.Dataset(arq4)
                variable_gfs05 = data_gfs05.variables[param][:]
                time_gfs05 = data_gfs05.variables['time']

                st4 = variable_gfs05[0:7]

                stg.append(np.sum(st4))

    # Creating weekly mobile media daily data
    stc = []
    sti = []
    stm = []

    for i in xrange(len(st1)):
        stc.append(np.sum(st1[i:i + 7]))
        if i == len(st1):
            print "acabou"
            exit()

    for j in xrange(len(st2)):
        sti.append(np.sum(st2[j:j + 7]))

        if j == len(st2):
            print "terminou"
            exit()

    for k in xrange(len(st3)):
        stm.append(np.sum(st3[k:k + 7]))
        if k == len(st3):
            print "fim"
            exit()

    dates = pd.date_range('2009-01-02', '2014-12-24', freq='D')
    GFS05 = Series(np.float64(stg), index=dates)
    INMET_ANA = Series(sti[:], index=dates)
    MERGE = Series(stm[:], index=dates)

    # gfs05chirps = DataFrame({'GFS05': GFS05, 'CHIRPS': CHIRPS})
    gfs05inmet_ana = DataFrame({'GFS05': GFS05, 'INMET_ANA': INMET_ANA})
    # gfs05merge = DataFrame({'GFS05': GFS05', 'MERGE': MERGE'})

    # Generating comparative graphs of pr thiessen
    gfs05inmet_ana.plot(kind='scatter', x='GFS05', y='INMET_ANA')
    plt.xlim(0, 350)
    plt.ylim(0, 350)
    plt.title(u'Dispersão do Thiessen de precipitação para a bacia {0}\n Peíodo: 01/01/2009 - 31/12/2014'.format(macro_name))
    plt.savefig('scatter_pr_thiessen_gfs05_inmet_ana_{0}.png'.format(macro_name))
    plt.show()

    gfs05inmet_ana['2009':'2014'].plot()
    plt.title(u'Thiessen de precipitação acumulada de 7 dias para a bacia {0}\n Peíodo: 02/01/2009 - 24/12/2014'.format(macro_name))
    plt.xlabel(u'Ano')
    plt.ylabel(u'mm')
    legenda = ('GFS05', 'INMET_ANA')
    plt.legend(legenda, frameon=False)
    plt.ylim(0, 350)
    plt.savefig('pr_thiessen_gfs05_inmet_ana_{0}.png'.format(macro_name))
    plt.show()
    exit()


