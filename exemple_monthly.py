# -*- coding: utf-8 -*-

__author__ = "Leidinice Silva"
__copyright__ = "Copyright 2016, Funceme Hydropy Project"
__credits__ = ["Francisco Vasconcelos Junior", "Marcelo Rodrigues", "Enzo Pinheiro"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marcelo Rodrigues"
__email__ = "leidinice.silvae@funceme.br"
__date__ = 20/9/2016

# skill echam
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

scale = 'monthly'
param = 'pr'
period = 'calibration'
start_date = '1810101'
end_date = '20101231'
hidropy_path = "/home/leidinice/Documentos/musf"
home = os.path.expanduser("~")


def function_next_month(date):
    return date + timedelta(days=calendar.monthrange(date.year, date.month)[1])


def function_last_daily(end_date):
    return datetime(end_date.year, end_date.month, calendar.monthrange(end_date.year, end_date.month)[1])


folders = os.listdir("{0}/hidropy/hidropy/shapes/basins/".format(hidropy_path))
basins = sorted(basin_dict(micro=True, basin_name='atlantico_sul'))  # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< change here!!

for basin in basins:
    basin_fullname = basin_dict(basin)[2]
    macro_name = basin_dict(basin)[1]


ste = []
ste_jan = []

mon = ['jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec']

for monthly in range(1, 12 + 1):

    link2 = home + "/io/echam46/hind8110/{0}/monthly/pr_thiessen/{1}".format(mon[monthly -1], macro_name)

    for year in range(1981, 2010 + 1):

        init_date = datetime(year, monthly, 1)
        next_date = init_date + relativedelta(months=+1)
        end_date = init_date + relativedelta(months=+3)
        end_date = function_last_daily(end_date)

        # print '{0}{1:02d}{2:02d}'.format(init_date.year,init_date.month, init_date.day)
        init_y = str(init_date)[0:4] + str(init_date)[5:7] + str(init_date)[8:10]
        start_y = str(next_date)[0:4] + str(next_date)[5:7] + str(next_date)[8:10]
        end_y = str(end_date)[0:4] + str(end_date)[5:7] + str(end_date)[8:10]
        # print init_y, start_y, end_y

        arq2 = "{0}/pr_monthly_echam46_hind8110_fcst_{1}_{2}_{3}_thiessen_{4}.nc".format(link2, init_y,
                                                                                           start_y, end_y,
                                                                                          basin_fullname)
        data_echam46 = netCDF4.Dataset(arq2)
        variable_echam46 = data_echam46.variables[param][:]
        time_echam46 = data_echam46.variables['time']
        st2 = variable_echam46[0:3]
        ste.append(np.sum(st2))
        ste_jan.append(st2[0::12])
        print ste_jan

