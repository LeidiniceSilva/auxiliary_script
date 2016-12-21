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
from netCDF4 import Dataset
from matplotlib import pyplot as plt

from hidropy.preprocessing.utils import basin_dict
from hidropy.preprocessing.write_thiessen import write_thiessen

__author__ = "Leidinice Silva"
__email__ = "leidinice.silvae@funceme.br"
__date__ = "19/12/2016"
__description__ = " Bias correction of pr_thiessen echam46 "

scale = 'monthly'
param = 'pr'
period = 'calibration'
start_date = '19810101'
end_date = '20101231'
home = os.path.expanduser("~")
hidropy_path = "/home/leidinice/Documentos/musf"

init_date = datetime.strptime(start_date, '%Y%m%d')
final_date = datetime.strptime(end_date, '%Y%m%d')


def gamma_correction(hind, clim_obs, model):
    mod = np.sort(hind)
    alpha_mod, loc_mod, beta_mod = ss.gamma.fit(model, loc=0)
    obs = np.sort(clim_obs)
    alpha_obs, loc_obs, beta_obs = ss.gamma.fit(obs, loc=0)

    corrected_model = []

    for i in model:
        prob = ss.gamma.cdf(i, alpha_mod, scale=beta_mod)
        corrected_model.append(ss.gamma.ppf(prob, alpha_obs, scale=beta_obs))

    return corrected_model


def arguments():
    global args

    parser = argparse.ArgumentParser(description=__description__)

    args = parser.parse_args()

if __name__ == '__main__':
    arguments()

    folders = os.listdir("{0}/hidropy/hidropy/shapes/basins/".format(hidropy_path))
    basins = sorted(basin_dict(micro=True, basin_name='doce'))  # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< change here!!

    for basin in basins:
        basin_fullname = basin_dict(basin)[2]
        macro_name = basin_dict(basin)[1]

        # open netcdf
        # st = []
        st1 = []
        stc1 = []
        stc2 = []
        stc3 = []

        link1 = home+"/io/inmet_ana_chirps/calibration/{0}/{1}_thiessen/{2}".format(scale, param, macro_name)
        arq1 = "{0}/{1}_{2}_inmet_ana_chirps_obs_19610101_20141231_thiessen_{3}.nc".format(link1, param, scale, basin_fullname)
        data1 = netCDF4.Dataset(arq1)
        variable1 = data1.variables[param][:].T
        time_obs = data1.variables['time']
        # st = variable1[252:603]
        st1 = variable1[240:600]
        stc1.append(st1[1::12])
        stc2.append(st1[2::12])
        stc3.append(st1[3::12])

        # open fcst
        st2 = []
        ste1 = []
        ste2 = []
        ste3 = []

        link2 = home + "/io/echam46/hind8110/jan/monthly/{0}_thiessen/{1}".format(param, macro_name)
        for year in range(1981, 2011):

            arq2 = "{0}/{1}_{2}_echam46_hind8110_fcst_{3}0101_{3}0201_{3}0430_thiessen_{4}.nc".format(link2, param, scale, year, basin_fullname)
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

        dates = np.arange(1981, 2011)
        pr_corrected1 = gamma_correction(dates, stc1[0], ste1)
        pr_corrected2 = gamma_correction(dates, stc2[0], ste2)
        pr_corrected3 = gamma_correction(dates, stc3[0], ste2)

        echam_corri = np.full((90), np.nan)
        echam_corri[0:30] = pr_corrected1
        echam_corri[30:60] = pr_corrected2
        echam_corri[60:90] = pr_corrected3

        obser = []
        for m in range(0, 30):
            obser = np.concatenate((obser, observado[m::30]), axis=0)

        echam = []
        for j in range(0, 30):
            echam = np.concatenate((echam, echam_corri[j::30]), axis=0)

        # Ploting graphs obser x echam
        print basin_fullname
        data = []
        for ano in range(1981, 2011):
            for mes in range(2, 5):
                data.append(datetime(ano, mes, 1))

        fig = plt.figure(figsize=(18, 6))
        plt.plot(np.array(data), obser, 'b', np.array(data), echam, 'r')
        plt.title(u'Pr_Thiessen - viés corrigido - Jan (FMA)\n bacia {0}'.format(basin_fullname))
        plt.ylim(0, 700)
        plt.ylabel(u'mm')
        plt.xlabel(u'anos')
        legenda = ('OBS', 'ECHAM46_corrigido')
        plt.legend(legenda, frameon=False)
        path_out = ('{0}/verificar_echam46_obs_usinas/pr_thiessen/pr_thiessen_mensal_echam46_corrigido/jan/{1}/'.format(hidropy_path, basin_dict(basin)[1]))
        plt.savefig(os.path.join(path_out, 'pr_thiessen_obs_echam46_corrigido_{0}.png'.format(basin_fullname)))
        plt.close('all')
        plt.cla()

        # Write output thiessen in netCDF4 file
        for k, yea in enumerate(range(1981, 2011)):
            aux = echam_corri[k::30]
            # print aux

            inid = date(yea, 1, 1)
            stad = date(yea, 2, 1)
            last_day_month = calendar.monthrange(yea, 4)[1]
            endd = date(yea, 4, last_day_month)

            start_y = str(inid)[0:4] + str(inid)[5:7] + str(inid)[8:10]
            new_y = str(stad)[0:4] + str(stad)[5:7] + str(stad)[8:10]
            end_y = str(endd)[0:4] + str(endd)[5:7] + str(endd)[8:10]
            # print start_y, new_y, end_y

            write_thiessen(aux, new_y, end_y, 'monthly', 'pr', 'echam46_hind8110', 'fcst', 'corrigido_{0}'.format(basin_fullname),
                       init_date=start_y, output_path=path_out)