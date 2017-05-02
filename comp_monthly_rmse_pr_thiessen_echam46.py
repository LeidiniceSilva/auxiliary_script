# -*- coding: utf-8 -*-

""" Verification of the precipitation ability of ECHAM46. """

import matplotlib as mpl; mpl.use('Agg')
import netCDF4
import calendar
import os
from datetime import date
from pylab import *
from matplotlib import pyplot as plt
from PyFuncemeClimateTools import ClimateStats as cs
from hidropy.utils.hidropy_utils import basin_dict, create_path
from matplotlib.font_manager import FontProperties

__author__ = "Leidinice Silva"
__email__ = "leidinice.silvae@funceme.br"
__date__ = "04/04/2017"
__description__ = "Verification of the precipitation ability of ECHAM46"

scale = 'monthly'
month = 'jan'
mon = 'jan'
param = 'pr'
period = 'hind8110'
model = "echam46"
home = os.path.expanduser("~")
hidropy_path = "/home/leidinice/documentos/projetos_git_funceme"

print month
print "int"

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
        if '_inc' not in(bas):
            bas_new.append(bas)

    # inc
    # bas_new = []
    # for bas in basins:
    #     if '_inc' in (bas):
    #         bas_new.append(bas)

    len_bas = len(bas_new)

    CORREL = []
    MBE = []
    RMSE = []
    plot_lines = []

    for basin in bas_new:
        basin_fullname = basin_dict(basin)[2]
        macro_name = basin_dict(basin)[1]
        print basin

       # open netcdf
        link1 = home + "/io/inmet_ana_chirps/calibration/{0}/{1}_thiessen/{2}".format(scale, param, macro_name)
        arq1 = "{0}/{1}_{2}_inmet_ana_chirps_obs_19610101_20141231_thiessen_{3}.nc".format(link1, param, scale,
                                                                                           basin_fullname)
        data1 = netCDF4.Dataset(arq1)
        variable1 = data1.variables[param][:].T
        time1 = data1.variables['time']
        st1 = variable1[600:648]
        st11 = variable1[612:648]

        sto_obs1 = []
        sto_obs2 = []
        sto_obs3 = []

        link2 = home + "/io/inmet_ana_chirps/operation/{0}/{1}_thiessen/{2}".format(scale, param, macro_name)
        arq2 = "{0}/{1}_{2}_inmet_ana_chirps_obs_20150101_20170228_thiessen_{3}.nc".format(link2, param, scale,
                                                                                           basin_fullname)
        data2 = netCDF4.Dataset(arq2)
        variable2 = data2.variables[param][:].T
        time2 = data2.variables['time']
        st2 = variable2[0:24]
        st22 = variable2[0:26]

        observ = np.full(72, np.nan)
        observ[0:48] = st1
        observ[48:72] = st2

        observa = np.full(62, np.nan)
        observa[0:36] = st11
        observa[36:62] = st22

        sto_obs1.append(observ[1::12])
        sto_obs2.append(observ[2::12])
        sto_obs3.append(observ[3::12])

        # open fcst
        sto_sim1 = []
        sto_sim2 = []
        sto_sim3 = []

        link3 = home + "/io/{0}/{1}/{2}/{3}/{4}_thiessen/{5}".format(model, period, mon, scale, param, macro_name)
        for year in range(2011, 2017):

            start_date = date(year, 12, 31)
            new_year = start_date + relativedelta(months=1)
            new_year_y = str(new_year)[0:4] + str(new_year)[5:7] + str(new_year)[8:10]

            arq3 = "{0}/{1}_{2}_{3}_{4}_fcst_{5}0101_{5}0201_{5}0430_thiessen_{6}.nc".format(link3, param, scale, model,
                                                                                             period, year,
                                                                                             basin_fullname)
            data3 = netCDF4.Dataset(arq3)
            variable3 = data3.variables[param][:]
            time3 = data3.variables['time']
            st3 = variable3
            sto_sim1.append(st3[0::3])
            sto_sim2.append(st3[1::3])
            sto_sim3.append(st3[2::3])

        # Calculate correlacao, bias and rmse of pr_thiessen
        corr1 = np.corrcoef(np.squeeze(sto_obs1), np.squeeze(sto_sim1))
        corr2 = np.corrcoef(np.squeeze(sto_obs2), np.squeeze(sto_sim2))
        corr3 = np.corrcoef(np.squeeze(sto_obs3), np.squeeze(sto_sim3))
        CORR = np.array([round(corr1[0][1], 3), round(corr2[0][1], 3), round(corr3[0][1], 3)])
        CORREL.append(CORR)

        vies1 = np.nanmean(np.squeeze(sto_sim1) - np.squeeze(sto_obs1))
        vies2 = np.nanmean(np.squeeze(sto_sim2) - np.squeeze(sto_obs2))
        vies3 = np.nanmean(np.squeeze(sto_sim3) - np.squeeze(sto_obs3))
        BIAS = np.array([vies1, vies2, vies3])
        MBE.append(BIAS)

        rmse1 = cs.compute_rmse(np.squeeze(sto_sim1), np.squeeze(sto_obs1))
        rmse2 = cs.compute_rmse(np.squeeze(sto_sim2), np.squeeze(sto_obs2))
        rmse3 = cs.compute_rmse(np.squeeze(sto_sim3), np.squeeze(sto_obs3))
        ERRO = np.array([rmse1, rmse2, rmse3])
        RMSE.append(ERRO)

        # print CORREL
        # print MBE
        # print RMSE
        # exit()

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
    labels = np.arange(len(bas_new))

    teste = []
    list_basins = []
    objects = ['FEV', 'MAR', 'ABR']

    for i in range(len(bas_new)):

        list_basins.append(bas_new[i].split('{0}_'.format(basin_name))[1].capitalize())

        fig.suptitle(u'Índices: CORRELAÇÃO, MBE e RMSE \n Modelo: {0} / Previsão: {1} - 2011_2016 \n Bacia: {2} - Usinas:'
                     u' Totais'.format(model.upper(), month.capitalize(), macro_name.capitalize()),
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
        ax[1].set_ylim([-400, 400])
        ax[1].axhline(y=0, linewidth=1, color='k')
        ax[1].text(-0.20, -400, u'B)', fontsize=200, fontweight='bold')
        ax[1].set_ylabel(u'MBE', fontsize=200, fontweight='bold')
        ax[1].set_xticks([0.5, 1.5, 2.5])
        ax[1].set_xticklabels(objects, fontsize=200, fontweight='bold')

        e = RMSE[:, i]
        f = ax[2].bar(z + count + i * bar3w, e, bar3w, color=line_colors[i], label=i)
        ax[2].set_ylim([0, 400])
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
        ax[1].axhline(-400, color='k')
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
    plt.figlegend(teste, list_basins, loc=8, ncol=6, prop=font)

    path_out = ('/home/leidinice/documentos/results/{0}/correct_rmse/{1}/{2}/{3}_thiessen/{4}/'
                'totais/'.format(model, mon, scale, param, macro_name))

    if not os.path.exists(path_out):
        create_path(path_out)

    plt.savefig(os.path.join(path_out, '{0}_totais_indices_RMSE_{1}_thiessen_{2}_{3}.png'.format(macro_name, param,
                                                                                                 model, month)), dpi=25)

    plt.close('all')
    plt.cla()
