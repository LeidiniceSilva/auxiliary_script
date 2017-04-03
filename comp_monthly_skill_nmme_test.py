# -*- coding: utf-8 -*-

""" Verification of the precipitation ability of NMME models """

import os
import requests
import calendar
import argparse
from datetime import date
from dateutil.relativedelta import relativedelta
from datetime import datetime

import numpy as np
import netCDF4
import itertools
import warnings
from matplotlib import pyplot as plt
from PyFuncemeClimateTools import ClimateStats as cs

from hidropy.utils.hidropy_utils import basin_dict, create_path


__author__ = "Leidinice Silva"
__email__ = "leidinice.silvae@funceme.br"
__date__ = "19/12/2016"
__description__ = " Bias correction of flow echam46 "

scale = 'monthly'
param = 'pr'
period = 'hind8210'
home = os.path.expanduser("~")
hidropy_path = "/home/leidinice/documentos/projetos_git_funceme"


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

    for basin in basins:
        basin_fullname = basin_dict(basin)[2]
        macro_name = basin_dict(basin)[1]
        print basin

        models = {'cmc1-cancm3': [0], 'cmc2-cancm4': [1], 'cola-rsmas-ccsm3': [2], 'cola-rsmas-ccsm4': [3],
                  'gfdl-cm2p5-flor-b01': [4], 'nasa-gmao-062012': [5], 'ncep-cfsv2': [6]}

        for model in ['cmc1-cancm3', 'cmc2-cancm4', 'cola-rsmas-ccsm3', 'cola-rsmas-ccsm4', 'gfdl-cm2p5-flor-b01',
                      'nasa-gmao-062012', 'ncep-cfsv2']:
            print model

            mes_prev = {'jan': [1, 2, 3], 'feb': [2, 3, 4], 'mar': [3, 4, 5], 'apr': [4, 5, 6], 'may': [5, 6, 7],
                        'jun': [6, 7, 8], 'jul': [7, 8, 9], 'aug': [8, 9, 10], 'sep': [9, 10, 11], 'oct': [10, 11, 0],
                        'nov': [11, 0, 1], 'dec': [0, 1, 2]}

            for monthly in ['jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec']:

                mon1 = mes_prev[monthly][0]
                mon2 = mes_prev[monthly][1]
                mon3 = mes_prev[monthly][2]

                # open netcdf
                stc1 = []
                stc2 = []
                stc3 = []

                link1 = home+"/io/inmet_ana_chirps/calibration/{0}/{1}_thiessen/{2}".format(scale, param, macro_name)
                arq1 = "{0}/{1}_{2}_inmet_ana_chirps_obs_19610101_20141231_thiessen_{3}.nc".format(link1, param, scale,
                                                                                                   basin_fullname)
                print arq1

                data1 = netCDF4.Dataset(arq1)
                variable1 = data1.variables[param][:].T
                time_obs = data1.variables['time']
                st1 = variable1[252:600]
                stc1.append(st1[mon1::12])
                stc2.append(st1[mon2::12])
                stc3.append(st1[mon3::12])
                print stc1
                exit()

                # open fcst
                ste1 = []
                ste2 = []
                ste3 = []

                link2 = home + "/io/NMME/{0}/hind8210/{1}/{2}/{3}/{4}".format(model, monthly, period, param, macro_name)
                for year in range(1981, 2010 + 1):
                    arq2 = "{0}/{1}_{2}_{3}_hind8210_fcst_{4}0201_{4}0301_{4}0531_thiessen_{5}.nc".format(link2, param,
                                                                                                          scale, model,
                                                                                                          year,
                                                                                                          basin_fullname)
                    data_echam46 = netCDF4.Dataset(arq2)
                    variable_echam46 = data_echam46.variables[param][:]
                    time_echam46 = data_echam46.variables['time']
                    st2 = variable_echam46
                    ste1.append(st2[0::3])
                    ste2.append(st2[1::3])
                    ste3.append(st2[2::3])

                # Calculate vies e corr pr_thiessen
                print basin_fullname

                corr1 = np.corrcoef(np.squeeze(stc1), np.squeeze(ste1))
                corr2 = np.corrcoef(np.squeeze(stc2), np.squeeze(ste2))
                corr3 = np.corrcoef(np.squeeze(stc3), np.squeeze(ste3))

                vies1 = np.nanmean(np.squeeze(ste1) - np.squeeze(stc1))
                vies2 = np.nanmean(np.squeeze(ste2) - np.squeeze(stc2))
                vies3 = np.nanmean(np.squeeze(ste3) - np.squeeze(stc3))

                rmse1 = cs.compute_rmse(np.squeeze(ste1), np.squeeze(stc1))
                rmse2 = cs.compute_rmse(np.squeeze(ste1), np.squeeze(stc1))
                rmse3 = cs.compute_rmse(np.squeeze(ste1), np.squeeze(stc1))

                print vies1, rmse1, round(corr1[0][1], 3)
                print vies2, rmse2, round(corr2[0][1], 3)
                print vies3, rmse3, round(corr3[0][1], 3)
                exit()

                # Generating comparative graphs of pr thiessen
                fig, axes = plt.subplots(figsize=(10, 5))

                ax1 = plt.subplot(311)
                plt.bar(pos, df['pre_score'], width, alpha=0.5, color='#EE3224', label=df['first_name'][0]))
                ax1.set_ylabel('Score')
                ax1.set_title('Test Subject Scores')
                ax1.set_xticks([p + 1.5 * width for p in pos])
                ax1.set_xticklabels(df['first_name'])
                plt.xlim(min(pos) - width, max(pos) + width * 4)
                plt.ylim([0, max(df['pre_score'] + df['mid_score'] + df['post_score'])])
                plt.legend(['Pre Score', 'Mid Score', 'Post Score'], loc='upper left')

                # share x only
                ax2 = plt.subplot(312, sharex=axes[0])
                plt.bar([p + width for p in pos], df['mid_score'], width, alpha=0.5, color='#F78F1E', label=df['first_name'][1])
                ax2.set_ylabel('Score')
                ax2.set_title('Test Subject Scores')
                ax2.set_xticks([p + 1.5 * width for p in pos])
                ax2.set_xticklabels(df['first_name'])
                plt.xlim(min(pos) - width, max(pos) + width * 4)
                plt.ylim([0, max(df['pre_score'] + df['mid_score'] + df['post_score'])])
                plt.legend(['Pre Score', 'Mid Score', 'Post Score'], loc='upper left')

                # share x and y
                ax3 = plt.subplot(313, sharex=ax1, sharey=ax1)
                plt.bar([p + width for p in pos], df['mid_score'], width, alpha=0.5, color='#FFC222'', label=df['first_name'][2])
                ax3.set_ylabel('Score')
                ax3.set_title('Test Subject Scores')
                ax3.set_xticks([p + 1.5 * width for p in pos])
                ax3.set_xticklabels(df['first_name'])
                plt.xlim(min(pos) - width, max(pos) + width * 4)
                plt.ylim([0, max(df['pre_score'] + df['mid_score'] + df['post_score'])])
                plt.legend(['Pre Score', 'Mid Score', 'Post Score'], loc='upper left')

                plt.show()
                exit()
