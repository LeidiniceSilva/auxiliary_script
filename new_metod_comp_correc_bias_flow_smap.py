# -*- coding: utf-8 -*-

""" Bias correction of flow py desagre gamma """

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
from hidropy.utils.write_flow import write_flow


__author__ = "Leidinice Silva"
__email__ = "leidinice.silvae@funceme.br"
__date__ = "19/12/2016"
__description__ = " Bias correction of flow echam46 "

scale = 'monthly'
param = 'flow'
period = 'calibration'
home = os.path.expanduser("~")
hidropy_path = "/home/leidinice/Documentos/musf"


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
    basins = sorted(basin_dict(micro=True))  # <<<<<<<<<<<<<<<<<<<<<<<<<<<< change here!!

    bas_new = []
    for bas in basins:
        if not os.path.exists(bas):
            bas_new.append(bas)

    for basin in bas_new:
        basin_fullname = basin_dict(basin)[2]
        macro_name = basin_dict(basin)[1]

        # open netcdf obs
        stc = []
        sto1 = []
        sto2 = []
        sto3 = []

        link1 = home+"/io/flow/smap_monthly/obs/{0}".format(macro_name)
        arq1 = "{0}/{1}_{2}_inmet_ana_chirps_obs_19770215_20161115_smap_{3}.nc".format(link1, param, scale, basin_fullname)
        data1 = netCDF4.Dataset(arq1)
        variable1 = data1.variables[param][:].T
        time1 = data1.variables['time']
        st_obs1 = variable1[47:407]
        st_obs2 = variable1[407:479]

        stc.append(st_obs1[1::12] + st_obs1[2::12] + st_obs1[3::12])
        sto1.append(st_obs2[1::12])
        sto2.append(st_obs2[2::12])
        sto3.append(st_obs2[3::12])

        # open netcdf mod
        ste = []

        link2 = home + "/io/{0}/smap_monthly/echam46/jan/{1}".format(param, macro_name)
        for year1 in range(1981, 2010+1):
            arq2 = "{0}/{1}_{2}_echam46_hind8110_fcst_{3}0101_{3}0215_{3}0415_smap_{4}.nc".format(link2, param, scale,
                                                                                              year1, basin_fullname)
            data2 = netCDF4.Dataset(arq2)
            variable2 = data2.variables[param][:]
            time2 = data2.variables['time']
            ste.append(np.sum(variable2))

        ste_fcst = []
        ste1 = []
        ste2 = []
        ste3 = []

        link3 = home + "/io/{0}/smap_monthly/echam46/jan/{1}".format(param, macro_name)
        for year2 in range(2011, 2017):
            arq3 = "{0}/{1}_{2}_echam46_hind8110_fcst_{3}0101_{3}0215_{3}0415_smap_{4}.nc".format(link3, param,
                                                                                                      scale, year2,
                                                                                                      basin_fullname)
            data3 = netCDF4.Dataset(arq3)
            variable3 = data3.variables[param][:]
            time3 = data3.variables['time']
            ste_fcst.append(np.sum(variable3))

            ste1.append(variable3[0::3])
            ste2.append(variable3[1::3])
            ste3.append(variable3[2::3])

        # Calculate vies and pr_correction echam46
        mon1_corr = []
        mon2_corr = []
        mon3_corr = []

        pr_corrected = gamma_correction(np.squeeze(ste), np.squeeze(stc), np.squeeze(ste_fcst))

        mon1_corr.append((np.squeeze(ste1) / np.squeeze(ste_fcst)) * pr_corrected)
        mon2_corr.append((np.squeeze(ste2) / np.squeeze(ste_fcst)) * pr_corrected)
        mon3_corr.append((np.squeeze(ste3) / np.squeeze(ste_fcst)) * pr_corrected)

        obser = np.full((18), np.nan)
        obser[0:6] = np.squeeze(sto1)
        obser[6:12] = np.squeeze(sto2)
        obser[12:18] = np.squeeze(sto3)

        echam_fcst = np.full((18), np.nan)
        echam_fcst[0:6] = np.squeeze(ste1)
        echam_fcst[6:12] = np.squeeze(ste2)
        echam_fcst[12:18] = np.squeeze(ste3)

        echam_corri = np.full((18), np.nan)
        echam_corri[0:6] = np.squeeze(mon1_corr)
        echam_corri[6:12] = np.squeeze(mon2_corr)
        echam_corri[12:18] = np.squeeze(mon3_corr)

        # print basin_fullname
        # print len(obser)
        # print type(obser)
        # print obser
        #
        # print len(echam_fcst)
        # print type(echam_fcst)
        # print echam_fcst
        #
        # print len(echam_corri)
        # print type(echam_corri)
        # print echam_corri

        observado = []
        for m in range(0, 6):
            observado = np.concatenate((observado, obser[m::6]), axis=0)

        vazao = []
        for n in range(0, 6):
            vazao = np.concatenate((vazao, echam_fcst[n::6]), axis=0)

        vazao_corr = []
        for o in range(0, 6):
            vazao_corr = np.concatenate((vazao_corr, echam_corri[o::6]), axis=0)

        # Ploting graphs obser x echam
        data = []
        for ano in range(2011, 2017):
            for mes in range(1, 4):
                data.append(datetime(ano, mes, 1))

        fig = plt.figure(figsize=(18, 6))
        plt.plot(np.array(data), observado, 'b', np.array(data), vazao, 'k', np.array(data), vazao_corr, 'r')
        plt.title(u'Comparação da Vazão - OBS x BRUTO x CORRIGIDO - Jan (F-M-A)\n bacia {0}'.format(basin_fullname))
        plt.ylim(0, 5000)
        plt.ylabel(u'mm')
        plt.xlabel(u'anos')
        legenda = ('OBS', 'BRUTO', 'CORRIGIDO')
        plt.legend(legenda, frameon=False)
        path_out1 = ('{0}/results/results_flow_basins/monthly_corrected_gamma_desag/figures/jan/'
                     '{1}/'.format(hidropy_path, basin_dict(basin)[1]))
        path_out2 = ('{0}/results/results_flow_basins/monthly_corrected_gamma_desag/jan/'
                     '{1}/'.format(hidropy_path, basin_dict(basin)[1]))

        if not os.path.exists(path_out1):
            create_path(path_out1)

        plt.savefig(os.path.join(path_out1, 'vazao_smap_corrigido_{0}.png'.format(basin_fullname)))
        plt.close('all')
        plt.cla()

        # Write output thiessen in netCDF4 file
        for k, yea in enumerate(range(2011, 2017)):
            aux = echam_corri[k::6]

            inid = date(yea, 1, 1)
            dat1 = date(yea, 1, 15)
            dat2 = date(yea, 1, 15)
            new_start = dat1 + relativedelta(months=1)
            new_endd = dat2 + relativedelta(months=3)

            start_y = str(inid)[0:4] + str(inid)[5:7] + str(inid)[8:10]
            new_y = str(new_start)[0:4] + str(new_start)[5:7] + str(new_start)[8:10]
            end_y = str(new_endd)[0:4] + str(new_endd)[5:7] + str(new_endd)[8:10]
            # print start_y, new_y, end_y
            # exit()

            if not os.path.exists(path_out2):
                create_path(path_out2)

            name_nc = write_flow(aux, new_y, end_y, 'monthly', 'flow', 'echam46_hind8110', 'fcst',
                                 'correc_{0}'.format(basin_fullname), 'smap', init_date=start_y, output_path=path_out2)
