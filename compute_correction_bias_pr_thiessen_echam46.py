# -*- coding: utf-8 -*-

""" Bias correction of pr_thiessen echam46 """

import os
import requests
import calendar
import argparse
from datetime import date, timedelta
from dateutil.relativedelta import relativedelta
from datetime import datetime

import numpy as np
import scipy.stats as ss
import plotly.plotly as py
import plotly.tools as tls
import netCDF4
import math
from netCDF4 import Dataset
from matplotlib import pyplot as plt

from hidropy.utils.hidropy_utils import basin_dict
from hidropy.utils.write_thiessen import write_thiessen

__author__ = "Leidinice Silva"
__email__ = "leidinice.silvae@funceme.br"
__date__ = "19/12/2016"
__description__ = " Bias correction of pr_thiessen echam46 "

scale = 'monthly'
param = 'pr'
period = 'calibration'
home = os.path.expanduser("~")
hidropy_path = "/home/leidinice/Documentos/musf"


def gamma_correction(hind, clim_obs, model):
    mod = np.sort(hind)
    alpha_mod, loc_mod, beta_mod = ss.gamma.fit(hind, loc=0)
    obs = np.sort(clim_obs)
    alpha_obs, loc_obs, beta_obs = ss.gamma.fit(obs, loc=0)

    corrected_hind = []
    for i in hind:
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
    basins = sorted(basin_dict(micro=True))  # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< change here!!

    for basin in basins:
        basin_fullname = basin_dict(basin)[2]
        macro_name = basin_dict(basin)[1]

        # open netcdf obs
        st = []
        st1 = []
        stc1 = []
        stc2 = []
        stc3 = []

        link1 = home+"/io/inmet_ana_chirps/calibration/{0}/{1}_thiessen/{2}".format(scale, param, macro_name)
        arq1 = "{0}/{1}_{2}_inmet_ana_chirps_obs_19610101_20141231_thiessen_{3}.nc".format(link1, param, scale, basin_fullname)
        data1 = netCDF4.Dataset(arq1)
        variable1 = data1.variables[param][:].T
        time_obs = data1.variables['time']
        st = variable1[252:603]
        st1 = variable1[240:600]
        stc1.append(st[0::12])
        stc2.append(st[1::12])
        stc3.append(st[2::12])

        # open netcdf mod
        st2 = []
        ste1 = []
        ste2 = []
        ste3 = []

        link2 = home + "/io/echam46/hind8110/dec/monthly/{0}_thiessen/{1}".format(param, macro_name)
        for year in range(1981, 2011):

            last_day = calendar.monthrange(year, 3)[1]

            start_date = date(year, 12, 1)
            new_year = start_date + relativedelta(months=1)
            new_start_date = date(year, 12, last_day)
            end_year = new_start_date + relativedelta(months=3)

            new_year_y = str(new_year)[0:4] + str(new_year)[5:7] + str(new_year)[8:10]
            end_year_y = str(end_year)[0:4] + str(end_year)[5:7] + str(end_year)[8:10]

            arq2 = "{0}/{1}_{2}_echam46_hind8110_fcst_{3}1201_{4}_{5}_thiessen_{6}.nc".format(link2, param, scale,
                                                                                                  year, new_year_y,
                                                                                                  end_year_y,
                                                                                                  basin_fullname)
            data_echam46 = netCDF4.Dataset(arq2)
            variable_echam46 = data_echam46.variables[param][:]
            time_echam46 = data_echam46.variables['time']
            ste1.append(variable_echam46[0::3])
            ste2.append(variable_echam46[1::3])
            ste3.append(variable_echam46[2::3])

        # Calculate vies and pr_correction echam46
        observado = np.full((90), np.nan)
        observado[0:30] = stc1[0]
        observado[30:60] = stc2[0]
        observado[60:90] = stc3[0]

        echam_bru = np.full((90), np.nan)
        echam_bru[0:30] = ste1
        echam_bru[30:60] = ste2
        echam_bru[60:90] = ste3

        pr_corrected1 = gamma_correction(ste1, stc1[0], ste1)
        pr_corrected2 = gamma_correction(ste1, stc2[0], ste2)
        pr_corrected3 = gamma_correction(ste1, stc3[0], ste2)

        echam_corri = np.full((90), np.nan)
        echam_corri[0:30] = pr_corrected1
        echam_corri[30:60] = pr_corrected2
        echam_corri[60:90] = pr_corrected3

        obser = []
        for m in range(0, 30):
            obser = np.concatenate((obser, observado[m::30]), axis=0)

        echam_b = []
        for n in range(0, 30):
            echam_b = np.concatenate((echam_b, echam_bru[n::30]), axis=0)

        echam_c = []
        for o in range(0, 30):
            echam_c = np.concatenate((echam_c, echam_corri[o::30]), axis=0)

        # Ploting graphs obser x echam
        print basin_fullname

        data = []
        for ano in range(1981, 2011):
            for mes in range(1, 4):
                data.append(datetime(ano, mes, 1))

        fig = plt.figure(figsize=(18, 6))
        plt.plot(np.array(data), obser, 'b', np.array(data), echam_b, '--k', np.array(data), echam_c, 'r')
        plt.title(u'Pr_Thiessen - viés corrigido - Dez (JFM)\n bacia {0}'.format(basin_fullname))
        plt.ylim(0, 700)
        plt.ylabel(u'mm')
        plt.xlabel(u'anos')
        legenda = ('OBS', 'ECHAM_bru', 'ECHAM46_corri')
        plt.legend(legenda, frameon=False)
        path_out1 = ('{0}/check_echam46_obs_basins/pr_thiessen/pr_thiessen_monthly_echam46_corrected/figures/dec/{1}/'.format(hidropy_path, basin_dict(basin)[1]))
        path_out2 = ('{0}/check_echam46_obs_basins/pr_thiessen/pr_thiessen_monthly_echam46_corrected/dec/{1}/'.format(hidropy_path, basin_dict(basin)[1]))
        plt.savefig(os.path.join(path_out1, 'pr_thiessen_obs_echam46_corrigido_{0}.png'.format(basin_fullname)))
        plt.close('all')
        plt.cla()

        # Write output thiessen in netCDF4 file
        for k, yea in enumerate(range(1981, 2011)):
            aux = echam_corri[k::30]

            dat1 = date(yea, 12, 1)
            last_day_mon = calendar.monthrange(yea, 1)[1]
            dat2 = date(yea, 12, last_day_mon)
            new_start = dat1 + relativedelta(months=1)
            new_endd = dat2 + relativedelta(months=3)

            start_y = str(dat1)[0:4] + str(dat1)[5:7] + str(dat1)[8:10]
            new_y = str(new_start)[0:4] + str(new_start)[5:7] + str(new_start)[8:10]
            end_y = str(new_endd)[0:4] + str(new_endd)[5:7] + str(new_endd)[8:10]

            name_nc = write_thiessen(aux, new_y, end_y, 'monthly', 'pr', 'echam46_hind8110', 'fcst', 'corrigido_{0}'.format(basin_fullname),
                                     init_date=start_y, output_path=path_out2)
