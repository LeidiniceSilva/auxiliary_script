# -*- coding: utf-8 -*-

""" Verification of the precipitation ability of flow. """

import matplotlib as mpl; mpl.use('Agg')
import scipy.stats as st
import netCDF4
import os
from datetime import date
from pylab import *
from matplotlib import pyplot as plt
from PyFuncemeClimateTools import ClimateStats as cs
from hidropy.utils.hidropy_utils import basin_dict, create_path
from matplotlib.font_manager import FontProperties
from basins_names_code import basins_flow

__author__ = "Leidinice Silva"
__email__ = "leidinice.silvae@funceme.br"
__date__ = "07/04/2017"
__description__ = "Verification of the precipitation ability of flow"

scale = 'monthly'
month = 'abr'
mon = 'apr'
param = 'flow'
period = 'hind8210'
home = os.path.expanduser("~")
hidropy_path = "/home/leidinice/documentos/projetos_git_funceme"

print month


def all_indices(basin_code, list_idx):
    indices = []
    idx = -1
    while True:
        try:
            idx = list_idx.index(basin_code, idx+1)
            indices.append(idx)
        except ValueError:
            break
    return indices


def load_flow_obs(basin_name, mon):

    path_in = '{0}/io/flow/smap_monthly/'.format(home)
    flow_obs = 'flow_obs_fluviometric_stations.asc'

    imp_data = np.loadtxt(path_in+flow_obs)
    list_idx = list(imp_data[:, 0])

    basin_code = basins_flow[basin_name]
    idx_list = all_indices(basin_code, list_idx)

    data_aux = imp_data[idx_list[0]:idx_list[-1] + 1, :]
    list_years = list(imp_data[idx_list[0]:idx_list[-1] + 1, 1])
    idx_flow = []

    for y in range(1982, 2010 + 1):
        idx_flow.append((all_indices(y, list_years))[0])

    var = data_aux[idx_flow[0]:idx_flow[-1]+1, 2:][:, mon-1]

    return var


def load_flow_obs_(basin_name, mon):

    path_in = '{0}/io/flow/smap_monthly/'.format(home)
    flow_obs = 'flow_obs_fluviometric_stations.asc'

    imp_data = np.loadtxt(path_in+flow_obs)
    list_idx = list(imp_data[:, 0])

    basin_code = basins_flow[basin_name]
    idx_list = all_indices(basin_code, list_idx)

    data_aux = imp_data[idx_list[0]:idx_list[-1] + 1, :]
    list_years = list(imp_data[idx_list[0]:idx_list[-1] + 1, 1])
    idx_flow = []

    for yy in range(1983, 2011 + 1):
        idx_flow.append((all_indices(yy, list_years))[0])

    var_ = data_aux[idx_flow[0]:idx_flow[-1]+1, 2:][:, mon-1]

    return var_

macros_total = ['amazonas', 'atlantico_leste', 'atlantico_sudeste', 'atlantico_sul', 'grande', 'iguacu',
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
    bas_new1 = []
    for bas1 in basins:
        if '_guapore' not in(bas1):
            bas_new1.append(bas1)

    bas_new2 = []
    for bas2 in bas_new1:
        if '_henry' not in(bas2):
            bas_new2.append(bas2)

    bas_new22 = []
    for bas22 in bas_new2:
        if '_jauru' not in (bas22):
            bas_new22.append(bas22)

    bas_new222 = []
    for bas222 in bas_new22:
        if '_xingo' not in (bas222):
            bas_new222.append(bas222)

    bas_new3 = []
    for bas3 in bas_new222:
        if '_inc' not in(bas3):
            bas_new3.append(bas3)

    len_bas = len(bas_new3)

    models = ['cmc1-cancm3', 'cmc2-cancm4', 'cola-rsmas-ccsm3', 'cola-rsmas-ccsm4', 'gfdl-cm2p5-flor-b01',
              'nasa-gmao-062012']

    for model in models:
        print model

        CORREL = []
        MBE = []
        RMSE = []
        plot_lines = []

        for basin in bas_new3:
            basin_fullname = basin_dict(basin)[2]
            macro_name = basin_dict(basin)[1]
            print basin

            stc1 = load_flow_obs(basin_fullname, 5)
            stc2 = load_flow_obs(basin_fullname, 6)
            stc3 = load_flow_obs(basin_fullname, 7)

            # open fcst
            ste1 = []
            ste2 = []
            ste3 = []

            link1 = home + "/io/{0}/smap_{1}/nmme/{2}/{3}/{4}".format(param, scale,  model, mon, macro_name)
            for year in range(1982, 2011):

                start_date = date(year, 12, 15)
                new_year = start_date + relativedelta(months=1)
                new_year_ = start_date + relativedelta(months=3)
                new_year_y = str(new_year)[0:4] + str(new_year)[5:7] + str(new_year)[8:10]
                new_year_y_ = str(new_year_)[0:4] + str(new_year_)[5:7] + str(new_year_)[8:10]

                arq1 = "{0}/{1}_{2}_{3}_{4}_fcst_{5}0401_{5}0515_{5}0715_smap_{6}".format(link1, param, scale, model,
                                                                                     period, year, basin_fullname)
                arq_new = arq1+"_petclim.nc"
                data1 = netCDF4.Dataset(arq_new)
                variable1 = data1.variables[param][:]
                time1 = data1.variables['time']
                st1 = variable1
                ste1.append(st1[0::3])
                ste2.append(st1[1::3])
                ste3.append(st1[2::3])

                print arq_new

            # Calculate correlacao, bias and rmse of pr_thiessen
            corr1 = np.corrcoef(np.squeeze(stc1), np.squeeze(ste1))
            corr2 = np.corrcoef(np.squeeze(stc2), np.squeeze(ste2))
            corr3 = np.corrcoef(np.squeeze(stc3), np.squeeze(ste3))
            CORR = np.array([round(corr1[0][1], 3), round(corr2[0][1], 3), round(corr3[0][1], 3)])
            CORREL.append(CORR)

            vies1 = np.nanmean(np.squeeze(ste1) - np.squeeze(stc1))
            vies2 = np.nanmean(np.squeeze(ste2) - np.squeeze(stc2))
            vies3 = np.nanmean(np.squeeze(ste3) - np.squeeze(stc3))
            BIAS = np.array([vies1, vies2, vies3])
            MBE.append(BIAS)

            desv = (np.squeeze(ste1) - np.squeeze(stc1)) ** 2
            myrmse = np.sqrt(np.nanmean(desv))

            rmse1 = cs.compute_rmse(np.squeeze(ste1), np.squeeze(stc1))
            rmse2 = cs.compute_rmse(np.squeeze(ste2), np.squeeze(stc2))
            rmse3 = cs.compute_rmse(np.squeeze(ste3), np.squeeze(stc3))
            ERRO = np.array([rmse1, rmse2, rmse3])
            RMSE.append(ERRO)

        # Generating comparative bar of pr thiessen
        cmap = plt.cm.gray
        fig, ax = plt.subplots(3, figsize=(210, 180))
        w = 0.14

        CORREL = np.array(CORREL).T
        MBE = np.array(MBE).T
        RMSE = np.array(RMSE).T

        bar1 = len(CORREL[:, 0])
        bar1w = w / bar1

        bar2 = len(MBE[:, 0])
        bar2w = w / bar2

        bar3 = len(RMSE[:, 0])
        bar3w = w / bar3

        x = np.arange(len(CORREL))
        y = np.arange(len(MBE))
        z = np.arange(len(RMSE))

        line_colors = cmap(np.linspace(0, 3))
        labels = np.arange(len(bas_new3))

        teste = []
        list_basins = []
        objects = ['MAI', 'JUN', 'JUL']

        for i in range(len(bas_new3)):

            list_basins.append(bas_new3[i].split('{0}_'.format(basin_name))[1].capitalize())

            fig.suptitle(u'Índices: CORRELAÇÃO, MBE e RMSE \n Modelo: {0} / Previsão: {1} - {2} \n Bacia: {3} - Usinas:'
                         u' Petclim'.format(model.upper(), month.capitalize(), period.upper(), macro_name.capitalize()),
                         fontsize=250, fontweight='bold')

            count = (1 - len_bas * bar1w) / 2
            a = CORREL[:, i]
            b = ax[0].bar(x + count + i * bar1w, a, bar1w, color=line_colors[i], label=i)
            ax[0].set_ylim([-1, 1])
            ax[0].set_xlim([0, 3])
            ax[0].axhline(y=0, linewidth=1, color='k')
            ax[0].text(-0.20, -1.0, u'A)', fontsize=200, fontweight='bold')
            ax[0].set_ylabel(u'Correlação',  fontsize=200, fontweight='bold')
            ax[0].set_xticks([0.5, 1.5, 2.5])
            ax[0].set_xticklabels(objects, fontsize=200, fontweight='bold')

            c = MBE[:, i]
            d = ax[1].bar(y + count + i * bar2w, c, bar2w, color=line_colors[i], label=i)
            ax[1].set_xlim([0, 3])
            ax[1].set_ylim([-1500, 1500])
            ax[1].axhline(y=0, linewidth=1, color='k')
            ax[1].text(-0.20, -1000, u'B)', fontsize=200, fontweight='bold')
            ax[1].set_ylabel(u'MBE', fontsize=200, fontweight='bold')
            ax[1].set_xticks([0.5, 1.5, 2.5])
            ax[1].set_xticklabels(objects, fontsize=200, fontweight='bold')

            e = RMSE[:, i]
            f = ax[2].bar(z + count + i * bar3w, e, bar3w, color=line_colors[i], label=i)
            ax[2].set_ylim([0, 1500])
            ax[2].set_xlim([0, 3])
            ax[2].text(-0.20, 0, u'C)', fontsize=200, fontweight='bold')
            ax[2].set_ylabel(u'RMSE',  fontsize=200, fontweight='bold')
            ax[2].set_xticks([0.5, 1.5, 2.5])
            ax[2].set_xticklabels(objects, fontsize=200, fontweight='bold')
            teste.append(f)

            ax[0].axvline(0, color='k')
            ax[0].axvline(1, color='#808080')
            ax[0].axvline(2, color='#808080')
            ax[0].axvline(3, color='k')
            ax[0].axhline(-1, color='k')
            ax[0].axhline(1, color='k')

            ax[1].axvline(0, color='k')
            ax[1].axvline(1, color='#808080')
            ax[1].axvline(2, color='#808080')
            ax[1].axvline(3, color='k')
            ax[1].axhline(-1500, color='k')
            ax[1].axhline(1500, color='k')

            ax[2].axvline(0, color='k')
            ax[2].axvline(1, color='#808080')
            ax[2].axvline(2, color='#808080')
            ax[2].axvline(3, color='k')
            ax[2].axhline(0, color='k')
            ax[2].axhline(1500, color='k')

            ax[0].tick_params(axis='both', which='major', labelsize=160, length=50, width=10, pad=50, labelcolor='k')
            ax[0].grid(True, which='major', linestyle='-.', linewidth='5')

            ax[1].tick_params(axis='both', which='major', labelsize=160, length=50, width=10, pad=50, labelcolor='k')
            ax[1].grid(True, which='major', linestyle='-.', linewidth='5')

            ax[2].tick_params(axis='both', which='major', labelsize=160, length=50, width=10, pad=50, labelcolor='k')
            ax[2].grid(True, which='major', linestyle='-.', linewidth='5')

        font = FontProperties(weight='bold', size=155)
        plt.figlegend(teste, list_basins, loc=8, ncol=6, prop=font)

        path_out = ('/home/leidinice/documentos/results/nmme/skill/{0}/{1}/{2}/{3}/{4}_smap/{5}/'
                    'petclim/'.format(model, period, mon, scale, param, macro_name))

        if not os.path.exists(path_out):
            create_path(path_out)

        plt.savefig(os.path.join(path_out, '{0}_petclim_indices_CORREL_MBE_RMSE_{1}_smap_{2}_{3}_{4}'
                                           '.png'.format(macro_name, param, model, month, period)), dpi=25)

        plt.close('all')
        plt.cla()
