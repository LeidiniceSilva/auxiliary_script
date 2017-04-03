# -*- coding: utf-8 -*-

__author__ = "Leidinice Silva"
__copyright__ = "Copyright 2016, Funceme Hydropy Project"
__credits__ = ["Jarbas Camurca", "Diogenes Fontenele", "Enzo Pinheiro"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Jarbas Camurca"
__email__ = "leidinice.silva@funceme.br"
__date__ = 9 / 18 / 2016

# Import Datas
import netCDF4
import numpy as np
import calendar
import os
from datetime import datetime
from hidropy.utils.write_thiessen import write_thiessen
from hidropy.utils.hidropy_utils import basin_dict, create_path

# Define parameters to calculate thiessen
scale = 'daily'
param = 'pet'
start_date = '19610101'
end_date = '20141231'
hidropy_path = "/home/leidinice/Documentos/musf"
home = os.path.expanduser("~")


init_date = datetime.strptime(start_date, '%Y%m%d')
final_date = datetime.strptime(start_date, '%Y%m%d')
ano = range(1961, 2015)


def climate(variable, time): # Setting function to calculate the climatology
    dummy = np.empty((54, 365))
    clim = np.empty((1, 365))

    for j, i in enumerate(ano):

        if calendar.isleap(i):

            index_01jan = netCDF4.date2index(datetime(i, 1, 1, 12), time)
            index_27fev = netCDF4.date2index(datetime(i, 2, 27, 12), time)
            index_28fev = netCDF4.date2index(datetime(i, 2, 28, 12), time)
            index_29fev = netCDF4.date2index(datetime(i, 2, 29, 12), time)
            index_01mar = netCDF4.date2index(datetime(i, 3, 1, 12), time)
            index_31dez = netCDF4.date2index(datetime(i, 12, 31, 12), time)

            mean_28 = np.nanmean(variable[index_28fev:index_29fev + 1])

            dummy[j, 0:58] = variable[index_01jan:index_27fev + 1]
            dummy[j, 58] = mean_28
            dummy[j, 59:] = variable[index_01mar:index_31dez + 1]

        else:
            index_01_jan = netCDF4.date2index(datetime(i, 1, 1, 12), time)
            index_31_dez = netCDF4.date2index(datetime(i, 12, 31, 12), time)
            dummy[j, 0:365] = variable[index_01_jan:index_31_dez +1]

    clim = np.nanmean(dummy[0:54, 0:365], axis=0)

    return clim

folders = os.listdir("{0}/hidropy/hidropy/shapes/basins/".format(hidropy_path))
basins = sorted(basin_dict(micro=True, basin_name='grande'))  # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< change here!!

for basin in basins:
    basin_fullname = basin_dict(basin)[2]
    macro_name = basin_dict(basin)[1]
    print "Running", basin

    link = "http://opendap4.funceme.br:8001/io/inmet/calibration/{0}/{1}_thiessen/{2}".format(scale, param, macro_name)

    arq = "{0}/{1}_{2}_inmet_obs_19610101_20141231_thiessen_{3}.nc".format(link, param, scale, basin_fullname)

    # open netcdf
    dataetp = netCDF4.Dataset(arq)
    variable_etp = dataetp.variables[param][:].T
    variable_etp[np.where(np.ma.getmask(variable_etp))] = np.NaN
    time_etp = dataetp.variables['time']

    # print variable.shape
    clim = climate(variable_etp, time_etp)
    print clim

    ic = 0

    for year in range(1961, 2014 + 1):
        icc = 0
        for month in range(0, 11 + 1):
            last_daymon = calendar.monthrange(year, month + 1)[1]
            for id in range(1, last_daymon + 1):
                cook = variable_etp[ic]

                if (month == 1 and id == 29) and np.isnan(cook):
                    variable_etp[ic] = clim[icc - 1]
                elif np.isnan(cook):
                    variable_etp[ic] = clim[icc]

                if month == 1 and id == 29:
                    icc = icc - 1

                ic += 1
                icc += 1

    variable_etp = variable_etp.reshape((-1, 1))
    np.savetxt('climatology_{0}.asc'.format(basin_fullname), clim)

    local_dir = '{0}/io/inmet/calibration/{0}/{1}_thiessen/{2}/'.format(home, scale, param, basin_dict(basin)[1])
    create_path(local_dir)
    write_thiessen(variable_etp, start_date, end_date, scale, param, 'inmet', 'obs', '{0}'.format(basin_fullname),
                   output_path=local_dir)
