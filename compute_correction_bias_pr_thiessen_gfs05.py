# -*- coding: utf-8 -*-

__author__ = "Leidinice Silva"
__copyright__ = "Copyright 2016, Funceme Hydropy Project"
__credits__ = ["Francisco Vasconcelos Junior", "Marcelo Rodrigues", "Diogenes Fontenele"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marcelo Rodrigues"
__email__ = "leidinice.silvae@funceme.br"
__date__ = 19/12/2016

# Bias correction of pr_thiessen gfs05

import os
import requests
import netCDF4
import calendar
import numpy as np
import scipy.stats as ss
import plotly.plotly as py
import plotly.tools as tls
from netCDF4 import Dataset
from datetime import date, timedelta
from dateutil.relativedelta import relativedelta
from datetime import datetime
from matplotlib import pyplot as plt
from hidropy.preprocessing.utils import basin_dict
from hidropy.preprocessing.write_thiessen import write_thiessen

scale = 'daily'
param = 'pr'
period = 'calibration'
start_date = '20090101'
end_date = '20141231'
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


folders = os.listdir("{0}/hidropy/hidropy/shapes/basins/".format(hidropy_path))
basins = sorted(basin_dict(micro=True, basin_name='doce'))  # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< change here!!

for basin in basins:
    basin_fullname = basin_dict(basin)[2]
    macro_name = basin_dict(basin)[1]

    # open netcdf
    st1 = []
    st2 = []

    link1 = home+"/io/inmet_ana_chirps_merge/calibration/{0}/{1}_thiessen/{2}".format(scale, param, macro_name)
    arq1 = "{0}/{1}_{2}_inmet_ana_chirps_merge_obs_19610101_20141231_thiessen_{3}.nc".format(link1, param, scale, basin_fullname)
    data1 = netCDF4.Dataset(arq1)
    variable1 = data1.variables[param][:].T
    time_chirps = data1.variables['time']
    st1 = variable1[17533:19716]

    stc = []

    for i in xrange(len(st1)):
        stc.append(np.sum(st1[i:i + 7]))
        if i == len(st1):
            print "acabou"
            exit()

    stg = []

    for year in range(2009, 2014+1):
        for month in range(1, 12+1):
            link2 = home+"/io/gfs05/{0:04d}/{0:04d}{1:02d}/pr_thiessen/{2}".format(year, month, macro_name)

            dias_mes = calendar.monthrange(year, month)[1]
            for daily in range(1, dias_mes + 1):
                if year == 2014 and month == 12 and daily == 24:
                    break
                init = np.datetime64('{0:04d}-{1:02d}-{2:02d}'.format(year, month, daily))
                start = init + np.timedelta64(1, 'D')
                end = init + np.timedelta64(7, 'D')

                init_y = str(init)[0:4] + str(init)[5:7] + str(init)[8:10]
                start_y = str(start)[0:4] + str(start)[5:7] + str(start)[8:10]
                end_y = str(end)[0:4] + str(end)[5:7] + str(end)[8:10]
                # print init_y, start_y, end_y

                arq2 = "{0}/pr_daily_gfs05_fcst_{1}00_{2}_{3}_thiessen_{4}.nc".format(link2, init_y, start_y, end_y,
                                                                                      basin_fullname)

                data_gfs05 = netCDF4.Dataset(arq2)
                variable_gfs05 = data_gfs05.variables[param][:]
                time_gfs05 = data_gfs05.variables['time']

                st2 = variable_gfs05[0:7]

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
    for i in range(0, 30):
        obser = np.concatenate((obser, observado[i::30]), axis=0)

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
    legenda = ('OBS', 'ECHAM46')
    plt.legend(legenda, frameon=False)
    path_out = ('{0}/verificar_echam46_obs_usinas/pr_thiessen/pr_thiessen_mensal_echam46_corrigido/jan/{1}/'.format(hidropy_path, basin_dict(basin)[1]))
    plt.savefig(os.path.join(path_out, 'pr_thiessen_obs_echam46_corrigido_{0}.png'.format(basin_fullname)))
    plt.close('all')
    plt.cla()

    # Write output thiessen in netCDF4 file
    for k, yea in enumerate(range(1981, 2011)):
        aux = echam_corri[k::30]
        print aux

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

exit()


















