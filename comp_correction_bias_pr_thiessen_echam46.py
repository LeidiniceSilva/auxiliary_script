# -*- coding: utf-8 -*-

""" Bias correction of pr_thiessen echam46 by gumbel / Gamma. """

import os
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


def gamma_correction(hind, clim_obs, fcst):  # Função Gumbel para correção de viés

    sh_mod = np.nanstd(hind) * np.pi / np.sqrt(6)
    sh_obs = np.nanmean(clim_obs) - 0.57721 * sh_mod

    corrected_fcst = []

    for i in fcst:

        prob = ss.genextreme.cdf(i, 1.5, loc=sh_mod, scale=sh_mod)
        corrected_fcst.append(ss.genextreme.ppf(prob, 1.5, loc=sh_obs, scale=sh_obs))

    return corrected_fcst


def gamma_correction(hind, clim_obs, fcst):  # Função Gamma para correção de viés

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

        stc_obs1 = []
        stc_obs2 = []
        stc_obs3 = []

        link1 = home+"/io/inmet_ana_chirps/calibration/{0}/{1}_thiessen/{2}".format(scale, param, macro_name)
        arq1 = "{0}/{1}_{2}_inmet_ana_chirps_obs_19610101_20141231_thiessen_{3}.nc".format(link1, param, scale,
                                                                                           basin_fullname)
        data_obs = netCDF4.Dataset(arq1)
        variable_obs = data_obs.variables[param][:].T
        time_obs = data_obs.variables['time']
        st1 = variable_obs[240:600]
        st2 = variable_obs[252:603]
        st3 = variable_obs[600:648]
        stc_obs1.append(st1[1::12])
        stc_obs2.append(st1[2::12])
        stc_obs3.append(st1[3::12])

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
        sto1.append(observ[1::12])
        sto2.append(observ[2::12])
        sto3.append(observ[3::12])

        # open netcdf hind
        ste_hind1 = []
        ste_hind2 = []
        ste_hind3 = []

        link2 = home + "/io/echam46/hind8110/jan/monthly/{0}_thiessen/{1}".format(param, macro_name)
        for year2 in range(1981, 2011):
            arq2 = "{0}/{1}_{2}_echam46_hind8110_fcst_{3}0101_{3}0201_{3}0430_thiessen_{4}.nc".format(link2, param,
                                                                                                      scale, year2,
                                                                                                      basin_fullname)
            data_hind = netCDF4.Dataset(arq2)
            variable_hind = data_hind.variables[param][:]
            time_hind = data_hind.variables['time']
            ste_hind1.append(variable_hind[0::3])
            ste_hind2.append(variable_hind[1::3])
            ste_hind3.append(variable_hind[2::3])

        # open netcdf fcst
        ste_fcst1 = []
        ste_fcst2 = []
        ste_fcst3 = []

        link3 = home + "/io/echam46/hind8110/jan/monthly/{0}_thiessen/{1}".format(param, macro_name)
        for year3 in range(2011, 2017):
            arq3 = "{0}/{1}_{2}_echam46_hind8110_fcst_{3}0101_{3}0201_{3}0430_thiessen_{4}.nc".format(link3, param,
                                                                                                      scale, year3,
                                                                                                      basin_fullname)
            data_fcst = netCDF4.Dataset(arq3)
            variable_fcst = data_fcst.variables[param][:]
            time_fcst = data_fcst.variables['time']
            ste_fcst1.append(variable_fcst[0::3])
            ste_fcst2.append(variable_fcst[1::3])
            ste_fcst3.append(variable_fcst[2::3])

        # Calculate vies and pr_correction echam46
        pr_corrected1 = gamma_correction(ste_hind1, stc_obs1[0], np.squeeze(ste_fcst1))
        pr_corrected2 = gamma_correction(ste_hind2, stc_obs2[0], np.squeeze(ste_fcst2))
        pr_corrected3 = gamma_correction(ste_hind3, stc_obs3[0], np.squeeze(ste_fcst3))

        obser = np.full((18), np.nan)
        obser[0:6] = sto1[0]
        obser[6:12] = sto2[0]
        obser[12:18] = sto3[0]

        echam_fcst = np.full((18), np.nan)
        echam_fcst[0:6] = ste_fcst1
        echam_fcst[6:12] = ste_fcst2
        echam_fcst[12:18] = ste_fcst3

        echam_corri = np.full((18), np.nan)
        echam_corri[0:6] = pr_corrected1
        echam_corri[6:12] = pr_corrected2
        echam_corri[12:18] = pr_corrected3

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

        obs_ini = np.full((90), np.nan)
        obs_ini[0:30] = stc_obs1[0]
        obs_ini[30:60] = stc_obs2[0]
        obs_ini[60:90] = stc_obs3[0]

        echam_bru = np.full((90), np.nan)
        echam_bru[0:30] = ste_hind1
        echam_bru[30:60] = ste_hind2
        echam_bru[60:90] = ste_hind3

        fig = plt.figure(figsize=(14, 8))
        plt.plot(np.array(data), observado, 'b', np.array(data), echam_b, '--k', np.array(data), echam_c, 'r')
        plt.title(u'Comparação Pr_Thiessen - OBS x BRUTO x CORRIGIDO - Jan (F-M-A)\n bacia {0}'.format(basin_fullname))
        plt.ylim(0, 700)
        plt.ylabel(u'mm')
        plt.xlabel(u'anos')
        legenda = ('OBS', 'BRUTO', 'CORRIGIDO')
        plt.legend(legenda, frameon=False)
        path_out1 = ('{0}/check_echam46_obs_basins/pr_thiessen/pr_thiessen_monthly_echam46_corrected_operation/gumbel/'
                     'figures/jan/{1}/'.format(hidropy_path, basin_dict(basin)[1]))
        path_out2 = ('{0}/check_echam46_obs_basins/pr_thiessen/pr_thiessen_monthly_echam46_corrected_operation/gumbel/'
                     'jan/{1}/'.format(hidropy_path, basin_dict(basin)[1]))

        if not os.path.exists(path_out1):
            create_path(path_out1)

        plt.savefig(os.path.join(path_out1, 'pr_thiessen_corrigido_{0}.png'.format(basin_fullname)))
        plt.close('all')
        plt.cla()

        # Write output thiessen in netCDF4 file
        for k, yea in enumerate(range(2011, 2017)):
            aux = echam_corri[k::6]

            dat1 = date(yea, 1, 1)
            last_day_mon = calendar.monthrange(yea, 1)[1]
            dat2 = date(yea, 1, last_day_mon)
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
                                     'corrigido_{0}'.format(basin_fullname), init_date=start_y, output_path=path_out2)















