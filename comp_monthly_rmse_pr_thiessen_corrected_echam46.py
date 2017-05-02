# -*- coding: utf-8 -*-

""" RMSE of bias correction by many different types for pr_thiessen echam46. """

import os
import argparse
import numpy as np
import scipy.stats as ss
import netCDF4
from datetime import date
from dateutil.relativedelta import relativedelta

from matplotlib import pyplot as plt
from matplotlib.font_manager import FontProperties
from PyFuncemeClimateTools import ClimateStats as cs
from hidropy.utils.hidropy_utils import basin_dict, create_path


__author__ = "Leidinice Silva"
__email__ = "leidinice.silvae@funceme.br"
__date__ = "30/03/2017"
__description__ = " RMSE of bias correction by many different types for pr_thiessen echam46 "

scale = 'monthly'
month = 'jan'
mon = 'jan'
model = 'echam46'
param = 'pr'
period = 'hind8110'
home = os.path.expanduser("~")
hidropy_path = "/home/leidinice/documentos/projetos_git_funceme"


def gamma_correction(hind, clim_obs, fcst):  # Função Gumbel para correção de viés

    sh_mod = np.nanstd(hind) * np.pi / np.sqrt(6)
    sh_obs = np.nanmean(clim_obs) - 0.57721 * sh_mod

    corrected_fcst_gamma = []

    for i in fcst:

        prob = ss.genextreme.cdf(i, 1.5, loc=sh_mod, scale=sh_mod)
        corrected_fcst_gamma.append(ss.genextreme.ppf(prob, 1.5, loc=sh_obs, scale=sh_obs))

    return corrected_fcst_gamma


def gumbel_correction(hind, clim_obs, fcst):  # Função Gumbel para correção de viés

    mod = np.sort(hind)
    alpha_mod, loc_mod, beta_mod = ss.gamma.fit(hind, loc=0)
    obs = np.sort(clim_obs)
    alpha_obs, loc_obs, beta_obs = ss.gamma.fit(obs, loc=0)

    corrected_fcst_gumbel = []

    for i in fcst:

        prob = ss.gamma.cdf(i, alpha_mod, scale=beta_mod)
        corrected_fcst_gumbel.append(ss.gamma.ppf(prob, alpha_obs, scale=beta_obs))

    return corrected_fcst_gumbel


def gamma_desag_correction(hind, clim_obs, fcst):

    mod = np.sort(hind)
    alpha_mod, loc_mod, beta_mod = ss.gamma.fit(hind, loc=0)
    obs = np.sort(clim_obs)
    alpha_obs, loc_obs, beta_obs = ss.gamma.fit(obs, loc=0)

    corrected_fcst_desag = []

    for i in fcst:
        prob = ss.gamma.cdf(i, alpha_mod, scale=beta_mod)
        corrected_fcst_desag.append(ss.gamma.ppf(prob, alpha_obs, scale=beta_obs))

    return corrected_fcst_desag


def arguments():
    global args

    parser = argparse.ArgumentParser(description=__description__)
    args = parser.parse_args()

if __name__ == '__main__':
    arguments()

    print month
    print "inc"

    macros_total = ['amazonas', 'atlantico_leste', 'atlantico_sudeste', 'atlantico_sul', 'doce', 'grande', 'iguacu',
                    'jacui', 'paraguai', 'paraiba_do_sul', 'parana', 'paranaiba', 'paranapanema', 'parnaiba',
                    'sao_francisco', 'tiete', 'tocantins', 'uruguai']

    macros_inc = ['amazonas', 'atlantico_leste', 'doce', 'grande', 'iguacu', 'jacui', 'paraiba_do_sul', 'parana',
                  'paranaiba', 'paranapanema', 'sao_francisco', 'tiete', 'tocantins', 'uruguai']

    for macro in macros_total:
        print macro

        folders = os.listdir("{0}/hidropy/hidropy/shapes/basins/".format(hidropy_path))
        basins = sorted(basin_dict(micro=True, basin_name=macro))
        basin_name = macro

        # total
        bas_new = []
        for bas in basins:
            if '_inc' not in (bas):
                bas_new.append(bas)

        # inc
        # bas_new = []
        # for bas in basins:
        #     if '_inc' in (bas):
        #         bas_new.append(bas)

        len_bas = len(bas_new)

        GAMMA = []
        GUMBEL = []
        GAMMA_DESAG = []
        plot_lines = []

        for basin in bas_new:
            basin_fullname = basin_dict(basin)[2]
            macro_name = basin_dict(basin)[1]
            print basin

            # open netcdf obs
            stc_obs1 = []
            stc_obs2 = []
            stc_obs3 = []
            stc_obs_acc = []

            link1 = home+"/io/inmet_ana_chirps/calibration/{0}/{1}_thiessen/{2}".format(scale, param, macro_name)
            arq1 = "{0}/{1}_{2}_inmet_ana_chirps_obs_19610101_20141231_thiessen_{3}.nc".format(link1, param, scale,
                                                                                               basin_fullname)

            data1 = netCDF4.Dataset(arq1)
            variable1 = data1.variables[param][:].T
            time1 = data1.variables['time']

            st1 = variable1[240:600]
            st2 = variable1[252:603]
            st3 = variable1[600:648]
            stc_obs1.append(st1[1::12])
            stc_obs2.append(st1[2::12])
            stc_obs3.append(st2[3::12])
            stc_obs_acc.append(st1[1::12] + st1[2::12] + st2[3::12])

            sto1 = []
            sto2 = []
            sto3 = []

            link2 = home + "/io/inmet_ana_chirps/operation/{0}/{1}_thiessen/{2}".format(scale, param, macro_name)
            arq2 = "{0}/{1}_{2}_inmet_ana_chirps_obs_20150101_20170131_thiessen_{3}.nc".format(link2, param, scale,
                                                                                               basin_fullname)
            data2 = netCDF4.Dataset(arq2)
            variable2 = data2.variables[param][:].T
            time2 = data2.variables['time']

            st4 = variable2[0:25]
            observ = np.full(73, np.nan)
            observ[0:48] = st3
            observ[48:73] = st4

            sto1.append(observ[1::12])
            sto2.append(observ[2::12])
            sto3.append(observ[3::12])

            # open netcdf hind and fcst
            ste_hind1 = []
            ste_hind2 = []
            ste_hind3 = []
            ste_hind_acc = []

            link3 = home + "/io/{0}/{1}/{2}/{3}/{4}_thiessen/{5}".format(model, period, mon, scale, param, macro_name)
            for year1 in range(1981, 2011):

                start_date = date(year1, 12, 31)
                new_year = start_date + relativedelta(months=1)
                new_year_y = str(new_year)[0:4] + str(new_year)[5:7] + str(new_year)[8:10]

                arq3 = "{0}/{1}_{2}_{3}_{4}_fcst_{5}0101_{5}0201_{5}0430_thiessen_{6}.nc".format(link3, param,  scale,
                                                                                                 model, period, year1,
                                                                                                 basin_fullname)
                data3 = netCDF4.Dataset(arq3)
                variable3 = data3.variables[param][:]
                time3 = data3.variables['time']

                ste_hind1.append(variable3[0::3])
                ste_hind2.append(variable3[1::3])
                ste_hind3.append(variable3[2::3])
                ste_hind_acc.append(np.sum(variable3))

            # open netcdf fcst
            ste_fcst1 = []
            ste_fcst2 = []
            ste_fcst3 = []
            ste_fcst_acc = []

            pr_corrected_gamma_desag1 = []
            pr_corrected_gamma_desag2 = []
            pr_corrected_gamma_desag3 = []

            link4 = home + "/io/{0}/{1}/{2}/{3}/{4}_thiessen/{5}".format(model, period, mon, scale, param, macro_name)
            for year2 in range(2011, 2017):

                start_date1 = date(year2, 12, 31)
                new_year1 = start_date1 + relativedelta(months=1)
                new_year_y1 = str(new_year1)[0:4] + str(new_year1)[5:7] + str(new_year1)[8:10]

                arq4 = "{0}/{1}_{2}_{3}_{4}_fcst_{5}0101_{5}0201_{5}0430_thiessen_{5}.nc".format(link4, param,  scale,
                                                                                             model, period, year2,
                                                                                             basin_fullname)
                data4 = netCDF4.Dataset(arq4)
                variable4 = data4.variables[param][:]
                time4 = data4.variables['time']

                ste_fcst1.append(variable4[0::3])
                ste_fcst2.append(variable4[1::3])
                ste_fcst3.append(variable4[2::3])
                ste_fcst_acc.append(np.sum(variable4))

            # Calculate vies by gammma, gumbe and desag. gamma
            pr_corrected_gamma1 = gamma_correction(np.squeeze(ste_hind1), stc_obs1[0], np.squeeze(ste_fcst1))
            pr_corrected_gamma2 = gamma_correction(np.squeeze(ste_hind2), stc_obs2[0], np.squeeze(ste_fcst2))
            pr_corrected_gamma3 = gamma_correction(np.squeeze(ste_hind3), stc_obs3[0], np.squeeze(ste_fcst3))

            pr_corrected_gumbel1 = gumbel_correction(np.squeeze(ste_hind1), stc_obs1[0], np.squeeze(ste_fcst1))
            pr_corrected_gumbel2 = gumbel_correction(np.squeeze(ste_hind2), stc_obs2[0], np.squeeze(ste_fcst2))
            pr_corrected_gumbel3 = gumbel_correction(np.squeeze(ste_hind3), stc_obs3[0], np.squeeze(ste_fcst3))

            pr_corrected_gamma_desag = gamma_desag_correction(np.squeeze(ste_hind_acc), np.squeeze(stc_obs_acc[0]), np.squeeze(ste_fcst_acc))
            pr_corrected_gamma_desag1.append((np.squeeze(ste_fcst1) / np.squeeze(ste_fcst_acc)) * pr_corrected_gamma_desag)
            pr_corrected_gamma_desag2.append((np.squeeze(ste_fcst2) / np.squeeze(ste_fcst_acc)) * pr_corrected_gamma_desag)
            pr_corrected_gamma_desag3.append((np.squeeze(ste_fcst3) / np.squeeze(ste_fcst_acc)) * pr_corrected_gamma_desag)

            # Calculate rmse obs x corrected (gamma, gumbel and desag. gamma)
            rmse_gamma1 = cs.compute_rmse(np.squeeze(ste_fcst1), np.squeeze(pr_corrected_gamma1))
            rmse_gamma2 = cs.compute_rmse(np.squeeze(ste_fcst2), np.squeeze(pr_corrected_gamma2))
            rmse_gamma3 = cs.compute_rmse(np.squeeze(ste_fcst3), np.squeeze(pr_corrected_gamma3))
            pr_corrected_gamma = np.array([rmse_gamma1, rmse_gamma2, rmse_gamma3])
            GAMMA.append(pr_corrected_gamma)

            rmse_gumbel1 = cs.compute_rmse(np.squeeze(ste_fcst1), np.squeeze(pr_corrected_gumbel1))
            rmse_gumbel2 = cs.compute_rmse(np.squeeze(ste_fcst2), np.squeeze(pr_corrected_gumbel2))
            rmse_gumbel3 = cs.compute_rmse(np.squeeze(ste_fcst3), np.squeeze(pr_corrected_gumbel3))
            pr_corrected_gumbel = np.array([rmse_gumbel1, rmse_gumbel2, rmse_gumbel3])
            GUMBEL.append(pr_corrected_gumbel)

            rmse_gamma_desag1 = cs.compute_rmse(np.squeeze(ste_fcst1), np.squeeze(pr_corrected_gamma_desag1[0]))
            rmse_gamma_desag2 = cs.compute_rmse(np.squeeze(ste_fcst2), np.squeeze(pr_corrected_gamma_desag2[0]))
            rmse_gamma_desag3 = cs.compute_rmse(np.squeeze(ste_fcst3), np.squeeze(pr_corrected_gamma_desag3[0]))
            pr_corrected_gamma_desag2 = np.array([rmse_gamma_desag1, rmse_gamma_desag2, rmse_gamma_desag3])
            GAMMA_DESAG.append(pr_corrected_gamma_desag2)

            # print GAMMA
            # print GUMBEL
            # print GAMMA_DESAG
            # exit()

        # Generating comparative bar of pr thiessen
        cmap = plt.cm.gray
        fig, ax = plt.subplots(3, figsize=(210, 180))
        w = 0.14

        GAMMA = np.array(GAMMA).T
        GUMBEL = np.array(GUMBEL).T
        GAMMA_DESAG = np.array(GAMMA_DESAG).T

        bar1 = len(GAMMA[:, 0])
        bar1w = w / bar1

        bar2 = len(GUMBEL[:, 0])
        bar2w = w / bar2

        bar3 = len(GAMMA_DESAG[:, 0])
        bar3w = w / bar3

        x = np.arange(len(GAMMA))
        y = np.arange(len(GUMBEL))
        z = np.arange(len(GAMMA_DESAG))

        line_colors = cmap(np.linspace(0, 3))
        labels = np.arange(len(bas_new))

        aux = []
        list_basins = []
        objects = ['FEV', 'MAR', 'ABR']

        for i in range(len(bas_new)):

            list_basins.append(bas_new[i].split('{0}_'.format(basin_name))[1].capitalize())

            fig.suptitle(u'Correções: GAMMA, GUMBEL e GAMMA DESAG. \n Modelo: {0} / Previsão: {1} - {2} \n Bacia:'
                         u' {3} - Usinas: Totais'.format(model.upper(), month.capitalize(), period.upper(),
                                                         macro_name.capitalize()), fontsize=250, fontweight='bold')

            count = (1 - len_bas * bar1w) / 2

            a = GAMMA[:, i]
            b = ax[0].bar(x + count + i * bar1w, a, bar1w, color=line_colors[i], label=i)
            ax[0].set_ylim([0, 400])
            ax[0].set_xlim([0, 3])
            ax[0].text(-0.20, 0, u'A)', fontsize=200, fontweight='bold')
            ax[0].set_ylabel(u'RMSE', fontsize=200, fontweight='bold')
            ax[0].set_xticks([0.5, 1.5, 2.5])
            ax[0].set_xticklabels(objects, fontsize=200, fontweight='bold')

            c = GUMBEL[:, i]
            d = ax[1].bar(y + count + i * bar2w, c, bar2w, color=line_colors[i], label=i)
            ax[1].set_ylim([0, 400])
            ax[1].set_xlim([0, 3])
            ax[1].text(-0.20, 0, u'B)', fontsize=200, fontweight='bold')
            ax[1].set_ylabel(u'RMSE', fontsize=200, fontweight='bold')
            ax[1].set_xticks([0.5, 1.5, 2.5])
            ax[1].set_xticklabels(objects, fontsize=200, fontweight='bold')

            e = GAMMA_DESAG[:, i]
            f = ax[2].bar(z + count + i * bar3w, e, bar3w, color=line_colors[i], label=i)
            ax[2].set_ylim([0, 400])
            ax[2].set_xlim([0, 3])
            ax[2].text(-0.20, 0, u'C)', fontsize=200, fontweight='bold')
            ax[2].set_ylabel(u'RMSE', fontsize=200, fontweight='bold')
            ax[2].set_xticks([0.5, 1.5, 2.5])
            ax[2].set_xticklabels(objects, fontsize=200, fontweight='bold')
            aux.append(f)

            ax[0].axvline(0, color='k')
            ax[0].axvline(1, color='#808080')
            ax[0].axvline(2, color='#808080')
            ax[0].axvline(3, color='k')
            ax[0].axhline(0, color='k')
            ax[0].axhline(400, color='k')

            ax[1].axvline(0, color='k')
            ax[1].axvline(1, color='#808080')
            ax[1].axvline(2, color='#808080')
            ax[1].axvline(3, color='k')
            ax[1].axhline(0, color='k')
            ax[1].axhline(400, color='k')

            ax[2].axvline(0, color='k')
            ax[2].axvline(1, color='#808080')
            ax[2].axvline(2, color='#808080')
            ax[2].axvline(3, color='k')
            ax[2].axhline(0, color='k')
            ax[2].axhline(400, color='k')

            ax[0].tick_params(axis='both', which='major', labelsize=160, length=50, width=10, pad=50, labelcolor='k')
            ax[0].grid(True, which='major', linestyle='-.', linewidth='5')

            ax[1].tick_params(axis='both', which='major', labelsize=160, length=50, width=10, pad=50, labelcolor='k')
            ax[1].grid(True, which='major', linestyle='-.', linewidth='5')

            ax[2].tick_params(axis='both', which='major', labelsize=160, length=50, width=10, pad=50, labelcolor='k')
            ax[2].grid(True, which='major', linestyle='-.', linewidth='5')

        font = FontProperties(weight='bold', size=155)
        plt.figlegend(aux, list_basins, loc=8, ncol=6, prop=font)

        path_out = ('/home/leidinice/documentos/results/{0}/correc_rmse/{1}/{2}/{3}/{4}_thiessen/{5}/'
                    'totais/'.format(model, period, mon, scale, param, macro_name))

        if not os.path.exists(path_out):
            create_path(path_out)

        plt.savefig(os.path.join(path_out, '{0}_totais_indices_RMSE_{1}_thiessen_{2}_{3}_'
                                           '{4}.png'.format(macro_name, param, model, month, period)), dpi=25)

        plt.close('all')
        plt.cla()