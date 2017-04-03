# -*- coding: utf-8 -*-

""" Bias correction of pet_thiessen echam46 by gamma"""

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
__description__ = "Bias correction of pet thiessen echam46"

scale = 'monthly'
param = 'pet'
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
        stc1 = []
        stc2 = []
        stc3 = []

        link1 = home+"/io/inmet_filled/calibration/{0}/{1}_thiessen/{2}".format(scale, param, macro_name)
        arq1 = "{0}/{1}_{2}_inmet_filled_obs_19610101_20141231_thiessen_{3}.nc".format(link1, param, scale,
                                                                                       basin_fullname)
        data1 = netCDF4.Dataset(arq1)
        variable1 = data1.variables[param][:].T
        time1 = data1.variables['time']
        st1 = variable1[252:603]
        st2 = variable1[240:600]
        stc1.append(st2[0::12])
        stc2.append(st2[1::12])
        stc3.append(st2[2::12])

        sto1 = []
        sto2 = []
        sto3 = []

        link2 = home + "/io/inmet_filled/operation/{0}/{1}_thiessen/{2}".format(scale, param, macro_name)
        arq2 = "{0}/{1}_{2}_inmet_filled_obs_20150101_20161231_thiessen_{3}.nc".format(link2, param, scale,
                                                                                       basin_fullname)
        data2 = netCDF4.Dataset(arq2)
        variable2 = data2.variables[param][:].T
        time2 = data2.variables['time']
        st3 = variable2[0:23]

        observ = np.full(71, np.nan)
        observ[0:48] = st2
        observ[48:71] = st3
        sto1.append(observ[8::12])
        sto2.append(observ[9::12])
        sto3.append(observ[10::12])

        # open netcdf mod
        ste1 = []
        ste2 = []
        ste3 = []

        link3 = home + "/io/echam46/hind8110/dec/monthly/{0}_thiessen/{1}".format(param, macro_name)
        for year1 in range(1981, 2011):

            last_day = calendar.monthrange(year1, 3)[1]
            start_date = date(year1, 12, 1)
            new_year = start_date + relativedelta(months=1)
            new_start_date = date(year1, 12, last_day)
            end_year = new_start_date + relativedelta(months=3)

            new_year_y = str(new_year)[0:4] + str(new_year)[5:7] + str(new_year)[8:10]
            end_year_y = str(end_year)[0:4] + str(end_year)[5:7] + str(end_year)[8:10]

            arq3 = "{0}/{1}_{2}_echam46_hind8110_fcst_{3}1201_{4}_{5}_thiessen_{6}.nc".format(link2, param, scale,
                                                                                              year1, new_year_y,
                                                                                              end_year_y,
                                                                                              basin_fullname)
            data3 = netCDF4.Dataset(arq3)
            variable3 = data3.variables[param][:]
            time3 = data3.variables['time']
            ste1.append(variable3[0::3])
            ste2.append(variable3[1::3])
            ste3.append(variable3[2::3])

        ste_fcst1 = []
        ste_fcst2 = []
        ste_fcst3 = []

        link4 = home + "/io/echam46/hind8110/dec/monthly/{0}_thiessen/{1}".format(param, macro_name)
        for year2 in range(2011, 2017):

            last_day1 = calendar.monthrange(year2, 3)[1]
            start_date1 = date(year2, 12, 1)
            new_year1 = start_date1 + relativedelta(months=1)
            new_start_date1 = date(year2, 12, last_day1)
            end_year1 = new_start_date1 + relativedelta(months=3)

            new_year_1 = str(new_year1)[0:4] + str(new_year1)[5:7] + str(new_year1)[8:10]
            end_year_1 = str(end_year1)[0:4] + str(end_year1)[5:7] + str(end_year1)[8:10]

            arq4 = "{0}/{1}_{2}_echam46_hind8110_fcst_{3}1201_{4}_{5}_thiessen_{6}.nc".format(link3, param, scale,
                                                                                              year2, new_year_1,
                                                                                              end_year_1,
                                                                                              basin_fullname)
            data4 = netCDF4.Dataset(arq4)
            variable4 = data4.variables[param][:]
            time4 = data4.variables['time']
            ste_fcst1.append(variable4[0::3])
            ste_fcst2.append(variable4[1::3])
            ste_fcst3.append(variable4[2::3])

        # Calculate vies and pr_correction echam46
        pr_corrected1 = gamma_correction(ste1, stc1[0], np.squeeze(ste_fcst1))
        pr_corrected2 = gamma_correction(ste2, stc2[0], np.squeeze(ste_fcst2))
        pr_corrected3 = gamma_correction(ste3, stc3[0], np.squeeze(ste_fcst3))

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

        fig = plt.figure(figsize=(18, 6))
        plt.plot(np.array(data), obser, 'b', np.array(data), echam_b, '--k', np.array(data), echam_c, 'r')
        plt.title(u'Comparação Pet_Thiessen - OBS x BRUTO x CORRIGIDO - Dez (J-F-M)\n bacia {0}'.format(basin_fullname))
        plt.ylim(0, 5000)
        plt.ylabel(u'mm')
        plt.xlabel(u'anos')
        legenda = ('OBS', 'BRUTO', 'CORRIGIDO')
        plt.legend(legenda, frameon=False)
        path_out1 = ('{0}/results/results_echam46_basins/pet_thiessen/monthly_corrected_gamma_norm/figures/'
                     'dec/{1}/'.format(hidropy_path, basin_dict(basin)[1]))
        path_out2 = ('{0}/results/results_echam46_basins/pet_thiessen/monthly_corrected_gamma_norm/'
                     'dec/{1}/'.format(hidropy_path, basin_dict(basin)[1]))

        if not os.path.exists(path_out1):
            create_path(path_out1)

        plt.savefig(os.path.join(path_out1, 'pet_thiessen_corrigido_{0}.png'.format(basin_fullname)))
        plt.close('all')
        plt.cla()

        # Write output thiessen in netCDF4 file
        for k, yea in enumerate(range(1981, 2011)):
            aux = echam_corri[k::6]

            dat1 = date(yea, 12, 1)
            last_day_mon = calendar.monthrange(yea, 12)[1]
            dat2 = date(yea, 12, last_day_mon)
            new_start = dat1 + relativedelta(months=1)
            new_endd = dat2 + relativedelta(months=3)

            start_y = str(dat1)[0:4] + str(dat1)[5:7] + str(dat1)[8:10]
            new_y = str(new_start)[0:4] + str(new_start)[5:7] + str(new_start)[8:10]
            end_y = str(new_endd)[0:4] + str(new_endd)[5:7] + str(new_endd)[8:10]
            # print start_y, new_y, end_y
            # exit()

            if not os.path.exists(path_out2):
                create_path(path_out2)

            name_nc = write_thiessen(aux, new_y, end_y, 'monthly', 'pet', 'echam46_hind8110', 'fcst', 'correc_{0}'
                                     .format(basin_fullname), init_date=start_y, output_path=path_out2)


