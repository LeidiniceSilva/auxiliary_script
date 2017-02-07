# -*- coding: utf-8 -*-

""" Bias correction of flow echam46 """

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
from netCDF4 import Dataset
from matplotlib import pyplot as plt

from hidropy.utils.hidropy_utils import basin_dict
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
    basins = sorted(basin_dict(micro=True, basin_name='iguacu'))  # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< change here!!

    for basin in basins:
        basin_fullname = basin_dict(basin)[2]
        macro_name = basin_dict(basin)[1]

        # open netcdf obs
        # st = []
        st1 = []
        stc1 = []
        stc2 = []
        stc3 = []

        link1 = home+"/io/flow/smap_monthly/obs/{0}".format(macro_name)
        arq1 = "{0}/{1}_{2}_inmet_ana_chirps_obs_19770215_20161115_smap_{3}.nc".format(link1, param, scale, basin_fullname)
        data1 = netCDF4.Dataset(arq1)
        variable1 = data1.variables[param][:].T
        time_obs = data1.variables['time']

        # st = variable1[252:603]
        st1 = variable1[47:407]
        stc1.append(st1[1::12])
        stc2.append(st1[2::12])
        stc3.append(st1[3::12])

        # open netcdf mod
        st2 = []
        ste1 = []
        ste2 = []
        ste3 = []

        link2 = home + "/io/{0}/smap_monthly/echam46/jan/{1}".format(param, macro_name)
        for year in range(1981, 2010+1):
            arq2 = "{0}/{1}_{2}_echam46_hind8110_fcst_{3}0101_{3}0215_{3}0415_smap_{4}.nc".format(link2, param, scale,
                                                                                              year, basin_fullname)

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

        print basin_fullname
        print len(observado), len(echam_bru), len(echam_corri)
        print type(observado), type(echam_bru), type(echam_corri)
        print np.min(observado), np.max(observado)
        print np.min(echam_bru), np.max(echam_bru)
        print np.min(echam_corri), np.max(echam_corri)
        print observado, echam_bru, echam_corri

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
        plt.title(u'Vazão - viés corrigido - Jan (FMA)\n bacia {0}'.format(basin_fullname))
        plt.ylim(0, 5000)
        plt.ylabel(u'mm')
        plt.xlabel(u'anos')
        legenda = ('OBS', 'VAZAO_bru', 'VAZAO_corri')
        plt.legend(legenda, frameon=False)
        path_out1 = ('{0}/check_echam46_obs_basins/flow_corrected/figures/jan/{1}/'.format(hidropy_path, basin_dict(basin)[1]))
        path_out2 = ('{0}/check_echam46_obs_basins/flow_corrected/jan/{1}/'.format(hidropy_path, basin_dict(basin)[1]))
        plt.savefig(os.path.join(path_out1, 'vazao_echam46_corrigido_{0}.png'.format(basin_fullname)))
        plt.close('all')
        plt.cla()

        # Write output thiessen in netCDF4 file
        for k, yea in enumerate(range(1981, 2010+1)):
            aux = echam_corri[k::30]

            inid = date(yea, 1, 1)
            dat1 = date(yea, 1, 15)
            dat2 = date(yea, 1, 15)
            new_start = dat1 + relativedelta(months=1)
            new_endd = dat2 + relativedelta(months=3)

            start_y = str(inid)[0:4] + str(inid)[5:7] + str(inid)[8:10]
            new_y = str(new_start)[0:4] + str(new_start)[5:7] + str(new_start)[8:10]
            end_y = str(new_endd)[0:4] + str(new_endd)[5:7] + str(new_endd)[8:10]

            name_nc = write_flow(aux, new_y, end_y, 'monthly', 'flow', 'echam46_hind8110', 'fcst', 'corrected_{0}'.format(basin_fullname),
                                 'smap', init_date=start_y, output_path=path_out2)

