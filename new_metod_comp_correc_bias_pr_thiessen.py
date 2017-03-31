# -*- coding: utf-8 -*-

""" Bias correction of pr_thiessen echam46 by desagre gamma """

import os
import requests
import calendar
import argparse
from datetime import date
from dateutil.relativedelta import relativedelta
from datetime import datetime

import numpy as np
import scipy.stats as ss
import netCDF4
from matplotlib import pyplot as plt

from hidropy.utils.hidropy_utils import basin_dict, create_path
from hidropy.utils.write_thiessen import write_thiessen

__author__ = "Leidinice Silva"
__email__ = "leidinice.silvae@funceme.br"
__date__ = "19/12/2016"
__description__ = " Bias correction of pr_thiessen echam46 "

scale = 'monthly'
param = 'pr'
period = 'calibration'
home = os.path.expanduser("~")
hidropy_path = "/home/leidinice/documentos/projetos_git_funceme"


def gamma_correction(hind, clim_obs, fcst):

    mod = np.sort(hind)
    alpha_mod, loc_mod, beta_mod = ss.gamma.fit(hind, loc=0)
    obs = np.sort(clim_obs)
    alpha_obs, loc_obs, beta_obs = ss.gamma.fit(obs, loc=0)

    corrected_fcst = []

    for i in fcst:

        prob = ss.gamma.cdf(i, alpha_mod, scale=beta_mod)
        corrected_fcst.append(ss.gamma.ppf(prob, alpha_obs, scale=beta_obs))

    return corrected_fcst


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
        st1 = []
        st2 = []
        st3 = []
        st4 = []

        stc_obs = []

        link1 = home+"/io/inmet_ana_chirps/calibration/{0}/{1}_thiessen/{2}".format(scale, param, macro_name)
        arq1 = "{0}/{1}_{2}_inmet_ana_chirps_obs_19610101_20141231_thiessen_{3}.nc".format(link1, param, scale,
                                                                                           basin_fullname)
        data_obs = netCDF4.Dataset(arq1)
        variable_obs = data_obs.variables[param][:].T
        time_obs = data_obs.variables['time']
        st1 = variable_obs[240:600]
        st2 = variable_obs[252:603]
        st3 = variable_obs[600:648]
        stc_obs.append(st1[8::12] + st1[9::12] + st1[10::12])

        sto = []
        sto1 = []
        sto2 = []
        sto3 = []
        link = home + "/io/inmet_ana_chirps/operation/{0}/{1}_thiessen/{2}".format(scale, param, macro_name)
        arq = "{0}/{1}_{2}_inmet_ana_chirps_obs_20150101_20161130_thiessen_{3}.nc".format(link, param, scale,
                                                                                          basin_fullname)
        data_obs2 = netCDF4.Dataset(arq)
        variable_obs2 = data_obs2.variables[param][:].T
        time_obs2 = data_obs2.variables['time']
        st4 = variable_obs2[0:23]
        observ = np.full(71, np.nan)
        observ[0:48] = st3
        observ[48:71] = st4
        sto1.append(observ[8::12])
        sto2.append(observ[9::12])
        sto3.append(observ[10::12])
        sto.append(observ[8::12] + observ[9::12] + observ[10::12])

        # open netcdf hind
        ste_hind = []
        link2 = home + "/io/echam46/hind8110/aug/monthly/{0}_thiessen/{1}".format(param, macro_name)
        for year2 in range(1981, 2011):
            arq2 = "{0}/{1}_{2}_echam46_hind8110_fcst_{3}0801_{3}0901_{3}1130_thiessen_{4}.nc".format(link2, param,
                                                                                                      scale, year2,
                                                                                                      basin_fullname)
            data_hind = netCDF4.Dataset(arq2)
            variable_hind = data_hind.variables[param][:]
            time_hind = data_hind.variables['time']
            ste_hind.append(np.sum(variable_hind))

        # open netcdf fcst
        ste_fcst = []
        ste_fcst1 = []
        ste_fcst2 = []
        ste_fcst3 = []
        link3 = home + "/io/echam46/hind8110/aug/monthly/{0}_thiessen/{1}".format(param, macro_name)
        for year3 in range(2011, 2017):
            arq3 = "{0}/{1}_{2}_echam46_hind8110_fcst_{3}0801_{3}0901_{3}1130_thiessen_{4}.nc".format(link3, param,
                                                                                                      scale, year3,
                                                                                                      basin_fullname)
            data_fcst = netCDF4.Dataset(arq3)
            variable_fcst = data_fcst.variables[param][:]
            time_fcst = data_fcst.variables['time']
            ste_fcst1.append(variable_fcst[0::3])
            ste_fcst2.append(variable_fcst[1::3])
            ste_fcst3.append(variable_fcst[2::3])
            ste_fcst.append(np.sum(variable_fcst))

        # Calculate vies and pr_correction echam46
        mon1_corr = []
        mon2_corr = []
        mon3_corr = []

        pr_corrected = gamma_correction(np.squeeze(ste_hind), np.squeeze(stc_obs), np.squeeze(ste_fcst))

        mon1_corr.append((np.squeeze(ste_fcst1) / np.squeeze(ste_fcst)) * pr_corrected)
        mon2_corr.append((np.squeeze(ste_fcst2) / np.squeeze(ste_fcst)) * pr_corrected)
        mon3_corr.append((np.squeeze(ste_fcst3) / np.squeeze(ste_fcst)) * pr_corrected)

        obser = np.full((18), np.nan)
        obser[0:6] = np.squeeze(sto1)
        obser[6:12] = np.squeeze(sto2)
        obser[12:18] = np.squeeze(sto3)

        echam_fcst = np.full((18), np.nan)
        echam_fcst[0:6] = np.squeeze(ste_fcst1)
        echam_fcst[6:12] = np.squeeze(ste_fcst2)
        echam_fcst[12:18] = np.squeeze(ste_fcst3)

        echam_corri = np.full((18), np.nan)
        echam_corri[0:6] = np.squeeze(mon1_corr)
        echam_corri[6:12] = np.squeeze(mon2_corr)
        echam_corri[12:18] = np.squeeze(mon3_corr)

        observado = []
        for m in range(0, 6):
            observado = np.concatenate((observado, obser[m::6]), axis=0)

        echam_b = []
        for n in range(0, 6):
            echam_b = np.concatenate((echam_b, echam_fcst[n::6]), axis=0)

        echam_c = []
        for o in range(0, 6):
            echam_c = np.concatenate((echam_c, echam_corri[o::6]), axis=0)

        # Ploting graphs obser x echam
        print basin_fullname
        data = []
        for ano in range(2011, 2017):
            for mes in range(1, 4):
                data.append(datetime(ano, mes, 1))

        fig = plt.figure(figsize=(14, 8))
        plt.plot(np.array(data), observado, 'b', np.array(data), echam_b, 'k', np.array(data), echam_c, 'r')
        plt.title(u'Comparação da Pr_Thiessen - OBS x BRUTO x CORRIGIDO - Ago (S-O-N)'
                  u'\n bacia {0}'.format(basin_fullname))
        plt.ylim(0, 700)
        plt.ylabel(u'mm')
        plt.xlabel(u'anos')
        legenda = ('OBS', 'BRUTO', 'CORRIGIDO')
        plt.legend(legenda, frameon=False)
        path_out1 = ('{0}/results/results_echam46_basins/pr_thiessen/monthly_corrected_operation/acc_gamma/figures/'
                     'ago/{1}/'.format(hidropy_path, basin_dict(basin)[1]))
        path_out2 = ('{0}/results/results_echam46_basins/pr_thiessen/monthly_corrected_operation/acc_gamma/'
                     'ago/{1}/'.format(hidropy_path, basin_dict(basin)[1]))

        if not os.path.exists(path_out1):
            create_path(path_out1)

        plt.savefig(os.path.join(path_out1, 'pr_thiessen_corrigido_{0}.png'.format(basin_fullname)))
        plt.close('all')
        plt.cla()

        # Write output thiessen in netCDF4 file
        for k, yea in enumerate(range(2011, 2017)):
            aux = echam_corri[k::6]

            dat1 = date(yea, 8, 1)
            last_day_mon = calendar.monthrange(yea, 8)[1]
            dat2 = date(yea, 8, last_day_mon)
            new_start = dat1 + relativedelta(months=1)
            new_endd = dat2 + relativedelta(months=3)

            start_y = str(dat1)[0:4] + str(dat1)[5:7] + str(dat1)[8:10]
            new_y = str(new_start)[0:4] + str(new_start)[5:7] + str(new_start)[8:10]
            end_y = str(new_endd)[0:4] + str(new_endd)[5:7] + str(new_endd)[8:10]
            # print start_y, new_y, end_y
            # exit()

            if not os.path.exists(path_out2):
                create_path(path_out2)

            name_nc = write_thiessen(aux, new_y, end_y, 'monthly', 'pr', 'echam46_hind8110', 'fcst',
                                     'correc_{0}'.format(basin_fullname), init_date=start_y, output_path=path_out2)















