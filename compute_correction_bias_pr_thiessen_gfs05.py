# -*- coding: utf-8 -*-

""" Bias correction of pr_thiessen gfs05 by gamma """

import os
import requests
import calendar
import argparse

import numpy as np
import pandas as pd
import scipy.stats as ss
import netCDF4
import math
from matplotlib import pyplot as plt

from hidropy.utils.hidropy_utils import basin_dict
from hidropy.utils.write_thiessen import write_thiessen

__author__ = "Leidinice Silva"
__email__ = "leidinice.silvae@funceme.br"
__date__ = "19/12/2016"
__description__ = " Bias correction of pr_thiessen gfs05 "

scale = 'daily'
param = 'pr'
period = 'calibration'
home = os.path.expanduser("~")
hidropy_path = "/home/leidinice/Documentos/musf"


def gamma_correction(hind, clim_obs, model):

    indx = np.where(hind >= 0.0)
    alpha_mod, loc_mod, beta_mod = ss.gamma.fit(hind[indx], loc=0)
    indy = np.where(clim_obs >= 0.0)
    alpha_obs, loc_obs, beta_obs = ss.gamma.fit(clim_obs[indy], loc=0)

    corrected_hind = []
    for i in hind[indx]:
        prob = ss.gamma.cdf(i, alpha_mod, scale=beta_mod)
        corrected_hind.append(ss.gamma.ppf(prob, alpha_obs, scale=beta_obs))

    return corrected_hind


def arguments():
    global args

    parser = argparse.ArgumentParser(description=__description__)
    args = parser.parse_args()

if __name__ == '__main__':
    arguments()

    folders = os.listdir("{0}/hidropy/hidropy/shapes/basins/".format(hidropy_path))
    basins = sorted(basin_dict(micro=True))

    bas_new = []
    for bas in basins:
        if bas != "amazonas_cachoeira_caldeirao":
            bas_new.append(bas)

    for basin in bas_new:
        basin_fullname = basin_dict(basin)[2]
        macro_name = basin_dict(basin)[1]
        print basin_fullname

        # open netcdf obs
        # print "read obs"
        st = []
        st1 = []
        link1 = home+"/io/inmet_ana_chirps_merge/calibration/{0}/{1}_thiessen/{2}".format(scale, param, macro_name)
        arq1 = "{0}/{1}_{2}_inmet_ana_chirps_merge_obs_19610101_20141231_thiessen_{3}.nc".format(link1, param, scale, basin_fullname)
        data1 = netCDF4.Dataset(arq1)
        variable1 = data1.variables[param][:].T
        time_chirps = data1.variables['time']
        st = variable1[17533:19723]

        for i in xrange(len(st)):
            st1.append(np.sum(st[i:i + 7]))
            if i == len(st):
                print "acabou"
                exit()

        stc = np.full((2184), np.nan)
        stc = st1[0:2184]

        # open netcdf mod
        # print "read sim"
        stg = []
        for year in range(2009, 2014+1):
            for month in range(1, 12+1):
                link2 = home+"/io/gfs05/{0:04d}/{0:04d}{1:02d}/pr_thiessen/{2}".format(year, month, macro_name)
                day_month = calendar.monthrange(year, month)[1]
                for daily in range(1, day_month + 1):
                    if year == 2014 and month == 12 and daily == 25:
                        break
                    init = np.datetime64('{0:04d}-{1:02d}-{2:02d}'.format(year, month, daily))
                    start = init + np.timedelta64(1, 'D')
                    end = init + np.timedelta64(7, 'D')
                    init_y = str(init)[0:4] + str(init)[5:7] + str(init)[8:10]
                    start_y = str(start)[0:4] + str(start)[5:7] + str(start)[8:10]
                    end_y = str(end)[0:4] + str(end)[5:7] + str(end)[8:10]

                    arq2 = "{0}/pr_daily_gfs05_fcst_{1}00_{2}_{3}_thiessen_{4}.nc".format(link2, init_y, start_y, end_y,
                                                                                          basin_fullname)
                    data_gfs05 = netCDF4.Dataset(arq2)

                    variable_gfs05 = data_gfs05.variables[param][:]
                    time_gfs05 = data_gfs05.variables['time']
                    stg.append(np.sum(variable_gfs05[0:7]))

        # Calculate pr_corrected of the gfs_bru
        obser = np.full((2184), np.nan)
        obser[0:2184] = stc

        gfs_bru = np.full((2184), np.nan)
        gfs_bru[0:2184] = stg
        gfs_bru = [0.0 if math.isnan(x) else x for x in gfs_bru]
        pr_corrected = np.full((2184), np.nan)
        pr_corrected = gamma_correction((np.array(gfs_bru)), obser, (np.array(gfs_bru)))

        # Ploting graphs obser x echam
        dates = pd.date_range('2009-01-02', '2014-12-25', freq='D')
        OBS = pd.Series(obser, index=dates)
        BRU = pd.Series(gfs_bru, index=dates)
        COR = pd.Series(pr_corrected, index=dates)
        obsbrucor = pd.DataFrame({'OBS': OBS, 'BRU': BRU, 'COR': COR})
        obsbrucor['2009':'2014'].plot(figsize=(18, 6), legend=True, fontsize=12)

        plt.title(u'Comparação Pr_Thiessen - OBS x BRUTO x CORRIGIDO - 02/01/2009 - 24/12/2014\n'
                  u' bacia {0}'.format(basin_fullname))
        plt.ylim(0, 700)
        plt.ylabel(u'mm')
        plt.xlabel(u'anos')
        legenda = ('OBS', 'BRUTO', 'CORRIGIDO')
        plt.legend(legenda, frameon=False)
        path_out1 = ('{0}/results/results_gfs05_basins/pr_thiessen_corrected_calib/gamma_norm/figuras/'
                     '{1}/'.format(hidropy_path, basin_dict(basin)[1]))
        plt.savefig(os.path.join(path_out1, 'pr_thiessen_corrigido_{0}.png'.format(basin_fullname)))
        plt.close('all')
        plt.cla()

        # Write output thiessen in netCDF4 file
        m = 0
        for yea in range(2009, 2014+1):
            for mon in range(1, 12+1):
                last_day_mon = calendar.monthrange(yea, mon)[1]
                for day in range(1, last_day_mon + 1):
                    if yea == 2014 and mon == 12 and day == 25:
                        break

                    path_out2 = ('{0}/results/results_gfs05_basins/pr_thiessen_corrected_calib/gamma_norm/'
                                 '{1}/{2}/{3}/'.format(hidropy_path, yea, mon, basin_dict(basin)[1]))
                    aux = np.full((7), np.nan)
                    aux[0:1] = pr_corrected[m]
                    m += 1

                    initt = np.datetime64('{0:04d}-{1:02d}-{2:02d}'.format(yea, mon, day))
                    startt = initt + np.timedelta64(1, 'D')
                    endd = initt + np.timedelta64(7, 'D')

                    initt_y = str(initt)[0:4] + str(initt)[5:7] + str(initt)[8:10] + '00'
                    startt_y = str(startt)[0:4] + str(startt)[5:7] + str(startt)[8:10]
                    endd_y = str(endd)[0:4] + str(endd)[5:7] + str(endd)[8:10]

                    name_nc = write_thiessen(aux, startt_y, endd_y, 'daily', 'pr', 'gfs05_acc', 'fcst',
                                             'correc_{0}'.format(basin_fullname), init_date=initt_y,
                                             output_path=path_out2)
