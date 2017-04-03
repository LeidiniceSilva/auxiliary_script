# -*- coding: utf-8 -*-

""" Verification of the precipitation ability of ECHAM46 model. """

__author__ = "Leidinice Silva"
__email__ = "leidinice.silvae@funceme.br"
__date__ = "20/09/2016"
__description__ = "Verification of the precipitation ability of ECHAM46 models"

import netCDF4
import os
import numpy
import matplotlib as mpl ; mpl.use('Agg')
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

    st = []
    st1 = []
    stc = []

    # open netcdf
    link1 = home+"/io/inmet_ana_chirps/calibration/{0}/{1}_thiessen/{2}".format(scale, param, macro_name)
    arq1 = "{0}/{1}_{2}_inmet_ana_chirps_obs_19610101_20141231_thiessen_{3}.nc".format(link1, param, scale,
                                                                                       basin_fullname)
    data1 = netCDF4.Dataset(arq1)
    variable1 = data1.variables[param][:].T
    time_obs = data1.variables['time']
    st = variable1[252:603]
    st1 = variable1[240:600]

    stc.append(st1[11::12][:] + st[0::12][:] + st[1::12][:])

    # open fcst
    st2 = []
    ste = []

    link2 = home + "/io/echam46/hind8110/nov/monthly/{0}_thiessen/{1}".format(param, macro_name)
    for year in range(1981, 2010 + 1):

        last_day = calendar.monthrange(year, 1)[1]
        start_date = date(year, 10, last_day)
        new_year = start_date + relativedelta(months=4)
        new_y = str(new_year)[0:4] + str(new_year)[5:7] + str(new_year)[8:10]

        arq2 = "{0}/{1}_{2}_echam46_hind8110_fcst_{3}1101_{3}1201_{4}_thiessen_{5}.nc".format(link2, param, scale,
                                                                                              year, new_y,
                                                                                              basin_fullname)

        data_echam46 = netCDF4.Dataset(arq2)
        variable_echam46 = data_echam46.variables[param][:]
        time_echam46 = data_echam46.variables['time']
        st2 = variable_echam46
        ste.append(np.sum(st2))

    # Generating comparative graphs of pr thiessen
    vies = np.nanmean(ste - stc[0])
    dates = np.arange(1981, 2011)
    plt.plot(dates, stc[0], 'b', dates, ste, 'r')
    txtbox(u'vies = {0}'.format(round(vies, 3)), 0.1, 0.87, 9, (1, 1, 1))
    plt.title(u'Thiessen de precipitação - Nov - Trim: DJF\n bacia {0}'.format(basin_fullname))
    plt.ylim(0, 1500)
    plt.xlabel(u'anos')
    plt.ylabel(u'mm')
    legenda = ('OBS', 'ECHAM46')
    plt.legend(legenda, frameon=False)
    plt.savefig('pet_thiessen_obs_echam46_{0}.png'.format(basin_fullname))
    plt.close('all')
    plt.cla()

    corr = numpy.corrcoef(stc[0], ste)
    plt.scatter(stc[0], ste)
    txtbox(u'r = {0}'.format(round(corr[0][1], 3)), 0.1, 0.87, 9, (1, 1, 1))
    plt.xlim(0, 1500)
    plt.ylim(0, 1500)
    plt.xlabel(u'OBS')
    plt.ylabel(u'ECHAM46')
    plt.title(u'Dispersão do Thiessen de precipitação\n Jan - Trim: FMA')
    plt.savefig('scatter_pet_thiessen_obs_echam46_{0}.png'.format(basin_fullname))
    plt.close('all')
    plt.cla()

    print vies, corr[0][1]
exit()








