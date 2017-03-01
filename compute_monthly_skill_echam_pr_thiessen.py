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
import os
import numpy
import numpy as np
import calendar
from datetime import date, timedelta
from dateutil.relativedelta import relativedelta
from hidropy.utils.hidropy_utils import basin_dict
from matplotlib import pyplot as plt

scale = 'monthly'
param = 'pr'
period = 'calibration'
start_date = '1810101'
end_date = '20101231'
hidropy_path = "/home/leidinice/Documentos/musf"
home = os.path.expanduser("~")

folders = os.listdir("{0}/hidropy/hidropy/shapes/basins/".format(hidropy_path))
basins = sorted(basin_dict(micro=True, basin_name='uruguai'))  # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< change here!!

for basin in basins:
    basin_fullname = basin_dict(basin)[2]
    macro_name = basin_dict(basin)[1]

    st1 = []
    stc1 = []
    stc2 = []
    stc3 = []

    # open netcdf
    link1 = home+"/io/inmet_ana_chirps/calibration/{0}/{1}_thiessen/{2}".format(scale, param, macro_name)
    arq1 = "{0}/{1}_{2}_inmet_ana_chirps_obs_19610101_20141231_thiessen_{3}.nc".format(link1, param, scale,
                                                                                       basin_fullname)
    data1 = netCDF4.Dataset(arq1)
    variable1 = data1.variables[param][:].T
    time_obs = data1.variables['time']
    st1 = variable1[240:600]
    stc1.append(st1[0::12])
    stc2.append(st1[1::12])
    stc3.append(st1[2::12])

    # open fcst
    st2 = []
    ste1 = []
    ste2 = []
    ste3 = []

    link2 = home + "/io/echam46/hind8110/jan/monthly/{0}_thiessen/{1}".format(param, macro_name)
    for year in range(1981, 2010 + 1):

        arq2 = "{0}/{1}_{2}_echam46_hind8110_fcst_{3}0101_{3}0201_{3}0430_thiessen_{4}.nc".format(link2, param, scale,
                                                                                          year, basin_fullname)

        data_echam46 = netCDF4.Dataset(arq2)
        variable_echam46 = data_echam46.variables[param][:]
        time_echam46 = data_echam46.variables['time']
        st2 = variable_echam46
        ste1.append(st2[0::3])
        ste2.append(st2[1::3])
        ste3.append(st2[2::3])

    # Calculate vies e corr pr_thiessen
    print basin_fullname
    print ste1
    print ste2
    print ste3
    # vies1 = np.nanmean(np.squeeze(ste1) - np.squeeze(stc1))
    # vies2 = np.nanmean(np.squeeze(ste2) - np.squeeze(stc2))
    # vies3 = np.nanmean(np.squeeze(ste3) - np.squeeze(stc3))
    #
    #
    # corr1 = numpy.corrcoef(np.squeeze(stc1), np.squeeze(ste1))
    # corr2 = numpy.corrcoef(np.squeeze(stc2), np.squeeze(ste2))
    # corr3 = numpy.corrcoef(np.squeeze(stc3), np.squeeze(ste3))
    #
    # print vies1, round(corr1[0][1], 3)
    # print vies2, round(corr2[0][1], 3)
    # print vies3, round(corr3[0][1], 3)

exit()









