# -*- coding: utf-8 -*-
"""Comparison between Thiessen computed by different databases."""

__author__ = "Leidinice Silva"
__copyright__ = "Copyright 2016, Funceme Hydropy Project"
__credits__ = ["Jarbinhas Camurca", "Marcelo Rodrigues", "Diogenes Fontenele", "Enzo Pinheiro"]
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
import matplotlib.gridspec as gridspec
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


def txtbox(text, xpos, ypos, fontsize, (col, lin, pos)):
    txtsp = plt.subplot(col, lin, pos)
    txt = text
    props = dict(boxstyle='round', facecolor='wheat', alpha=0)
    txtsp.text(xpos, ypos, txt, transform=txtsp.transAxes, fontsize=fontsize, verticalalignment='top',
               bbox=props)


folders = os.listdir("{0}/hidropy/hidropy/shapes/basins/".format(hidropy_path))
basins = sorted(basin_dict(micro=True, basin_name='tiete'))  # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< change here!!

for basin in basins:
    basin_fullname = basin_dict(basin)[2]
    macro_name = basin_dict(basin)[1]

    st1 = []
    st2 = []

    # open netcdf
    link1 = home+"/io/inmet_ana_chirps_merge/calibration/{0}/{1}_thiessen/{2}".format(scale, param, macro_name)
    arq1 = "{0}/{1}_{2}_inmet_ana_chirps_merge_obs_19610101_20141231_thiessen_{3}.nc".format(link1, param, scale, basin_fullname)
    data1 = netCDF4.Dataset(arq1)
    variable1 = data1.variables[param][:].T
    time_chirps = data1.variables['time']
    st1 = variable1[17533:19716]

    stc = []

    for i in xrange(len(st1)):
        stc.append(np.sum(st1[i:i + 7]))
        if i == len(st1):
            print "acabou"
            exit()

    stg = []

    # open fcst
    for year in range(2009, 2014+1):
        for month in range(1, 12+1):
            link2 = home+"/io/gfs05/{0:04d}/{0:04d}{1:02d}/pr_thiessen/{2}".format(year, month, macro_name)

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

                arq2 = "{0}/pr_daily_gfs05_fcst_{1}00_{2}_{3}_thiessen_{4}.nc".format(link2, init_y, start_y, end_y,
                                                                                      basin_fullname)

                data_gfs05 = netCDF4.Dataset(arq2)
                variable_gfs05 = data_gfs05.variables[param][:]
                time_gfs05 = data_gfs05.variables['time']

                st2 = variable_gfs05[0:7]

                stg.append(np.sum(st2))

    # Creating weekly mobile media daily data
    dates = pd.date_range('2009-01-02', '2014-12-24', freq='D')
    print dates
    exit()
    OBS = pd.Series(stc[:], index=dates)
    GFS05 = pd.Series(np.float64(stg), index=dates)
    obsgfs05 = pd.DataFrame({'OBS': OBS, 'GFS05': GFS05})

    # Generating comparative graphs of pr thiessen
    obs_nan = OBS.isnull().sum()
    gfs_nan = GFS05.isnull().sum()
    obsgfs05['2009':'2014'].plot()
    txtbox(u'N° de falhas\nGFS05 = {0}\nOBS = {1}'.format(gfs_nan, obs_nan), 0.79, 0.87, 9, (1, 1, 1))
    plt.title(u'Thiessen de precipitação semanal da bacia {0}\n Período: 02/01/2009 - 24/12/2014'.format(basin_fullname))
    plt.ylabel(u'mm')
    legenda = ('GFS05', 'OBS')
    plt.legend(legenda, frameon=False)
    plt.ylim(0, 500)
    plt.savefig('pr_thiessen_obs_gfs05_{0}.png'.format(basin_fullname))

    corr = OBS.corr(GFS05)
    obsgfs05.plot(kind='scatter', x='GFS05', y='OBS')
    txtbox(u'r = {0}'.format(round(corr, 3)), 0.1, 0.87, 9, (1, 1, 1))
    plt.xlim(0, 500)
    plt.ylim(0, 500)
    plt.title(u'Dispersão do Thiessen de precipitação \n Período: 02/01/2009 - 24/12/2014')
    plt.savefig('scatter_pr_thiessen_obs_gfs05_{0}.png'.format(basin_fullname))
exit()


