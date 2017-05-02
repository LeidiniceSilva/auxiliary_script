# -*- coding: utf-8 -*-

"""
This script remove monthly mean flow bias from NMME, ECHAM4.6 and
RSM2008 models. Ths scritp focus on operational period and creates files in
the same flow_smap pattern.
"""

import argparse
import calendar
from dateutil.relativedelta import *
from datetime import datetime, date
import glob
from hidropy.utils.hidropy_utils import date2index, basin_dict
from hidropy.utils.write_thiessen import write_thiessen
from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
from PyFuncemeClimateTools.Thiessen import thiessen
import os
from os.path import expanduser
import scipy.stats as ss
home = expanduser("~")
import shutil

__author__ = 'Leidinice Silva'
__email__ = 'leidinicesilva@gmail.com'
__date__ = '26/04/2017'
__description__ = 'This script remove monthly mean flow bias from NMME,' \
                ' ECHAM4.6 and RSM2008 models. Ths scritp focus on operational' \
                ' period and creates files in the same flow_smap pattern.'


def arguments():
    global args

# Usage:
# script.py --month_target=jan --year_target=1981 --localdir=/home/musf/io --model_name=echam46 --basin=amazonas --micro

    parser = argparse.ArgumentParser(description=__description__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--month_target', default='', help='Insert initial forecast month.')
    parser.add_argument('--year_target', default='', help='Insert initial forecast year.')
    parser.add_argument('--local_dir', default='/io', help='Path to data directory.')
    parser.add_argument('--model_name', choices=['echam46', 'rsm2008', 'NMME', 'all'], type=str, nargs='?',
                        help='Model options.')
    parser.add_argument('--basin', help='Name of input basin.', default='')
    parser.add_argument('--macro', action='store_true', help='Compute flow for all macro-basins.')
    parser.add_argument('--micro', action='store_true', help='Compute flow for all micro-basins of input basin.')
    parser.add_argument('--all_basins', action='store_true', help='True to compute flow for all macro-basins '
                                                                  'and micro-basins.')
    args = parser.parse_args()
    return args


def define_initial_parameters(modelname):
    """ Define target data and models list.

    :param: modelname: str
    :return: target_date: datatime
    :return: model name list: list
    """

    # List Models
    if modelname == 'all':
        list_model = ['echam46', 'cmc1-cancm3', 'cmc2-cancm4',
                      'cola-rsmas-ccsm3','cola-rsmas-ccsm4', 'gfdl-cm2p5-flor-b01', 
                      'ncep-cfsv2', 'nasa-gmao-062012', 'rsm2008']
    elif modelname == 'NMME':
        list_model = ['cmc1-cancm3', 'cmc2-cancm4',
                      'cola-rsmas-ccsm3','cola-rsmas-ccsm4', 'gfdl-cm2p5-flor-b01', 
                      'ncep-cfsv2', 'nasa-gmao-062012']

    else:
        list_model = [modelname]

    if not month:
        target_date = datetime(date.today().year, date.today().month, 15) 
    else:
        target_date = datetime(int(year), int(month), 15)

    return target_date, list_model


def define_dates(target_date):
    """ Define date start run, and dates that the forecast starts and ends.

    :param: target_date: datatime
    :return: str_mon: str
    :return: start_rundate: str
    :return: start_fcstdate: str
    :return: end_fcstdate: str
    """

    target_ifcstdate = target_date + relativedelta(months=1)
    target_efcstdate = target_date + relativedelta(months=3)

    str_mon = target_date.strftime("%b").lower()

    start_rundate = '{0}{1:02d}01'.format(target_date.year, target_date.month)
    start_fcstdate = '{0}{1:02d}15'.format(target_ifcstdate.year, target_ifcstdate.month)
    end_fcstdate = '{0}{1:02d}15'.format(target_efcstdate.year, target_efcstdate.month)

    return str_mon, start_rundate, start_fcstdate, end_fcstdate

def import_fcstmodel_data(model_name, target_date, basin):

    """ Import fcst model data.

    :param: model_name: str ('ECHAM4.6', 'RSM2008', 'NMME')
    :param: target_date: datetime
    :param: basin: str
    :return: flow Smap time serie
    :rtype:list
    """

    str_mon, start_rundate, start_fcstdate, end_fcstdate = define_dates(target_date)
    basin_full_name = basin_dict(basin)[2]
    basin_name = basin_dict(basin)[1]

    if model_name=='echam46' or model_name=='rsm2008':
        hind_period = 'hind8110'
        dir_modelname = '{0}{1}/{2}/{3}/{4}/monthly/{5}_smap/{6}/'.format(home, localdir, model_name, hind_period,
                                                                          str_mon, var_name, basin_name)

    else:
        hind_period = 'hind8210'
        dir_modelname = '{0}{1}/nmme/{2}/{3}/{4}/monthly/{5}_smap/{6}/'.format(home, localdir, model_name, hind_period,
                                                                               str_mon, var_name, basin_name)

    file_modelname = "{0}_monthly_{1}_{2}_fcst_{3}_{4}_{5}_smap_{6}.nc". \
            format(var_name, model_name, hind_period, start_rundate,
                   start_fcstdate, end_fcstdate, basin_full_name)

    try:
        input_data = Dataset(os.path.join(dir_modelname, file_modelname))
        var = input_data.variables[var_name][:]
        input_data.close()
        flag = True
    except:
        print 'File: {0} is not available'.\
            format(os.path.join(dir_modelname, file_modelname))
        var = []
        flag = False
    
    return var, flag

def import_hindmodel_data(model_name, target_date, basin, idx):
    """ Import hind model data.

        :param: model_name: str ('ECHAM4.6', 'RSM2008', 'NMME')
        :param: target_date: datetime
        :param: basin: str
        :return: precipitation Thiessen time serie
        :rtype: list
        """

    if model_name == 'echam46' or model_name == 'rsm2008':
        ihind = 1981
        hind_list = np.full((30), np.nan)
    else:
        hind_list = np.full((29), np.nan)
        ihind = 1982

    year_list = np.arange(ihind, 2010 + 1)

    for ii, y in enumerate(year_list):

        target_rundate = datetime(y, target_date.month, 15)
        str_mon, start_rundate, start_fcstdate, end_fcstdate = define_dates(target_rundate)

        basin_full_name = basin_dict(basin)[2]
        basin_name = basin_dict(basin)[1]

        if model_name=='echam46' or model_name=='rsm2008':
            hind_period = 'hind8110'
            dir_modelname = '{0}{1}/{2}/{3}/{4}/monthly/{5}_smap/{6}/'.format(home, localdir, model_name, hind_period,
                                                                              str_mon, var_name,
                       basin_name)
        else:
            hind_period = 'hind8210'
            dir_modelname = '{0}{1}/nmme/{2}/{3}/{4}/monthly/{5}_smap/{6}/'.format(home, localdir, model_name,
                                                                                   hind_period, str_mon, var_name,
                                                                                   basin_name)

        file_modelname = "{0}_monthly_{1}_{2}_fcst_{3}_{4}_{5}_smap_" \
                         "{6}.nc".format(var_name, model_name, hind_period, start_rundate, start_fcstdate, end_fcstdate,
                                         basin_full_name)

        try:
            input_data = Dataset(os.path.join(dir_modelname, file_modelname))
            hind_list[ii] = input_data.variables[var_name][idx]
            input_data.close()
        except:
            print 'Hindcast File: {0} is not available'.\
                format(os.path.join(dir_modelname, file_modelname))
            hind_list[ii] = np.nan
    
    return hind_list

def import_hindobs_data(model_name, target_date, basin, idx):
    """ Import obs model data.

    :param: model_name: str ('ECHAM4.6', 'RSM2008', 'NMME')
    :param: target_date: datetime
    :param: basin: str
    :param: idx: int
    :return: precipitation Thiessen time serie
    """
    
    basin_full_name = basin_dict(basin)[2]
    basin_name = basin_dict(basin)[1]
    fcstdate = target_date + relativedelta(months=idx + 1)

    if model_name == 'echam46' or model_name == 'rsm2008':
        ihind = 1981
    else:
        ihind = 1982

    fhind = 2010

    if fcstdate.year != target_date.year:
        ihind = ihind + 1
        fhind = fhind + 1

    file_name = '{0}/io/flow/calibration/smap_{1}/obs/{2}/{3}_{1}_{4}_obs_{5}_smap_' \
                '{6}.nc'.format(home, time_freq, basin_name, var_name, data_base, hind_obs_period, basin_full_name)

    try:
        inc = Dataset(file_name)
        time_var = inc.variables['time']
        iidx = date2index(datetime(ihind, 01, 15), time_var)
        fidx = date2index(datetime(fhind, 12, 15), time_var)
        var = inc.variables[var_name][iidx:fidx + 1]
        return var[fcstdate.month-1:len(var):12]
    except:
        print 'Observational Calibration File: {0} is not available'.format(os.path.join(file_name))
        exit()

def gamma_correction(hind, clim_obs, fcst):
    """ Remove precipitaiton bias via Gama distribution.

    :param: hind: list
    :param: clim_obs: list
    :param: fcst: list
    :return: precipitation
    """

    # Accumulating the Quarter
    hind = np.sum(hind, axis=0)
    obs = np.sum(clim_obs, axis=0)

    # Deleting zeros and np.nan of hindcast serie!
    hind = hind[np.nonzero(hind)]
    hind = hind[~np.isnan(hind)]

    zeros_obs = obs[np.where(obs == 0.)]

    # Probability to find zeros in observational serie!
    q = len(zeros_obs) / np.float(len(obs))
    obs_nozero = obs[np.nonzero(obs)]

    if (len(obs_nozero) < 20) or (len(hind) < 20):
        print 'Error in the hindcast length'
        exit()

    mod = np.sort(hind)
    alpha_mod, loc_mod, beta_mod = ss.gamma.fit(mod, floc=0)
    obs_nozero = np.sort(obs_nozero)
    alpha_obs, loc_obs, beta_obs = ss.gamma.fit(obs_nozero, floc=0)

    if fcst == 0.0:
       corrected_fcst = np.copy(fcst)
    else:
       prob = q + (1-q) * (ss.gamma.cdf(fcst, alpha_mod, scale=beta_mod))
       corrected_fcst = ss.gamma.ppf(prob, alpha_obs, scale=beta_obs)

    return np.squeeze(corrected_fcst)

if __name__ == "__main__":

    arguments()
    global month, year, localdir, var_name, time_freq, data_base, hind_obs_period

    # Target basin
    bname = args.basin
    macros = args.macro
    micros = args.micro
    allbasins = args.all_basins

    # Local Directory,variable name and Model
    localdir = args.local_dir
    modelname = args.model_name

    bname = 'amazonas'
    modelname = 'echam46'
    month = 'jan'
    year = 1981
    localdir = '/home/leidinice/io'
    var_name = 'flow'
    time_freq = 'monthly'
    data_base = 'inmet_ana_chirps'
    hind_obs_period = '19770215_20161215'

    # Target Period
    month = args.month_target
    year = args.year_target

    # Target Basin
    if macros or micros or allbasins:
        basin_list = basin_dict(bname, macro=macros, micro=micros, all_basins=allbasins)
    else:
        basin_list = [bname]

    target_date, list_model = define_initial_parameters(modelname)

    for modname in list_model:
        
        print "Run Model: {0}".format(modname)

        for basin in basin_list:
            
            print 'Processing Basin: {0}'.format(basin)

            basin_full_name = basin_dict(basin)[2]
            basin_name = basin_dict(basin)[1]

            str_mon, start_rundate, start_fcstdate, end_fcstdate = define_dates(target_date)
            fcst_list, flag = import_fcstmodel_data(modname, target_date, basin)

            if flag:

                if modname == 'echam46' or modname == 'rsm2008':
                    hind = np.full((3,30), np.nan)
                    obs = np.full((3,30), np.nan)
                    hindname = 'hind8110'
                    path_out = '{0}{1}/{2}/{3}/{4}/monthly/{5}_smap_cor/{6}/'.format(home, localdir, modname, hindname,
                                                                                     str_mon, var_name, basin_name)
                    print path_out
                    exit()
                else:
                    hind = np.full((3,29), np.nan)
                    obs = np.full((3,29), np.nan)
                    hindname = 'hind8210'
                    path_out = '{0}{1}/nmme/{2}/{3}/{4}/monthly/{5}_smap_cor/{6}/'.format(home, localdir, modname,
                                                                                          hindname, str_mon, var_name,
                                                                                          basin_name)
                
                if not os.path.exists(path_out):
                    os.makedirs(path_out)
                    
                fcst = np.nansum(fcst_list)
                for idx in range(0, 3):
                    hind[idx, :] = import_hindmodel_data(modname, target_date, basin, idx)
                    obs[idx, :] = import_hindobs_data(modname, target_date, basin, idx)
                
                # Removing bias - Gamma Disaggregation Method
                fcst_cor = gamma_correction(hind, obs, fcst)

                # Disaggregating preciptation
                corr = []
                for i in range(0, 3):
                    corr.append((fcst_list[i] / float(fcst)) * fcst_cor)

                corr = np.array(corr)
                corr[np.isnan(corr)] = -999.
                
                # Saving data
                write_thiessen(corr, start_fcstdate, end_fcstdate, time_freq, var_name,
                               '{0}_{1}'.format(modname, hindname), 'fcst', '{0}_cor'.format(basin_full_name),
                               init_date=start_rundate, output_path=path_out)
            else:
                continue

