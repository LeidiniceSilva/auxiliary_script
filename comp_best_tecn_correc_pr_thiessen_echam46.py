# -*- coding: utf-8 -*-

""" Bias correction of pr_thiessen echam46 """

import os
import warnings
import argparse

import numpy as np
import pandas as pd
import scipy.stats as st
import netCDF4
import matplotlib
from matplotlib import pyplot as plt
from datetime import datetime

from hidropy.utils.hidropy_utils import basin_dict, create_path

__author__ = "Leidinice Silva"
__email__ = "leidinice.silvae@funceme.br"
__date__ = "19/12/2016"
__description__ = " Bias correction of pr_thiessen echam46 "

freq = 'monthly'
param = 'pr'
period = 'calibration'
home = os.path.expanduser("~")
hidropy_path = "/home/leidinice/documentos/projetos_git_funceme"

matplotlib.rcParams['figure.figsize'] = (16.0, 12.0)
matplotlib.style.use('ggplot')


def best_fit_distribution(data, bins='', ax=None):
    """Model data by finding best fit distribution to data"""
    # Get histogram of original data
    y, x = np.histogram(data, bins=bins, normed=True)
    x = (x + np.roll(x, -1))[:-1] / 2.0


    # Distributions to check
    DISTRIBUTIONS = [st.alpha, st.anglit, st.arcsine, st.beta, st.betaprime, st.bradford, st.burr, st.cauchy, st.chi,
                     st.chi2, st.cosine, st.dgamma, st.dweibull, st.erlang, st.expon, st.exponnorm, st.exponweib,
                     st.exponpow, st.f, st.fatiguelife, st.fisk, st.foldcauchy, st.foldnorm, st.frechet_r, st.frechet_l,
                     st.genlogistic, st.genpareto, st.gennorm, st.genexpon, st.genextreme, st.gausshyper, st.gamma,
                     st.gengamma, st.genhalflogistic, st.gilbrat, st.gompertz, st.gumbel_r, st.gumbel_l, st.halfcauchy,
                     st.halflogistic, st.halfnorm, st.halfgennorm, st.hypsecant, st.invgamma, st.invgauss,
                     st.invweibull, st.johnsonsb, st.johnsonsu, st.ksone, st.kstwobign, st.laplace, st.levy, st.levy_l,
                     st.levy_stable, st.logistic, st.loggamma, st.loglaplace, st.lognorm, st.lomax, st.maxwell,
                     st.mielke, st.nakagami, st.ncx2, st.ncf, st.nct, st.norm, st.pareto, st.pearson3, st.powerlaw,
                     st.powerlognorm, st.powernorm, st.rdist, st.reciprocal, st.rayleigh, st.rice, st.recipinvgauss,
                     st.semicircular, st.t, st.triang, st.truncexpon, st.truncnorm, st.tukeylambda, st.uniform,
                     st.vonmises, st.vonmises_line, st.wald, st.weibull_min, st.weibull_max, st.wrapcauchy]

    # Best holders
    best_distribution = st.norm
    best_params = (0.0, 1.0)
    best_sse = np.inf

    # Estimate distribution parameters from data
    for distribution in DISTRIBUTIONS:

        # Try to fit the distribution
        try:
            # Ignore warnings from data that can't be fit
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore')
                # fit dist to data
                params = distribution.fit(data)

                # Separate parts of parameters
                arg = params[:-2]
                loc = params[-2]
                scale = params[-1]

                # Calculate fitted PDF and error with fit in distribution
                pdf = distribution.pdf(x, loc=loc, scale=scale, *arg)
                sse = np.sum(np.power(y - pdf, 2.0))

                # if axis pass in add to plot
                try:
                    if ax:
                        pd.Series(pdf, x).plot(ax=ax)

                except Exception:
                    pass

                # identify if this distribution is better
                if best_sse > sse > 0:
                    best_distribution = distribution
                    best_params = params
                    best_sse = sse
        except Exception:
            pass

    return (best_distribution, best_distribution.name, best_params)


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
        st1 = []
        st2 = []
        st3 = []
        st4 = []

        stc_obs1 = []
        stc_obs2 = []
        stc_obs3 = []

        link1 = home+"/io/inmet_ana_chirps/calibration/{0}/{1}_thiessen/{2}".format(freq, param, macro_name)
        arq1 = "{0}/{1}_{2}_inmet_ana_chirps_obs_19610101_20141231_thiessen_{3}.nc".format(link1, param, freq, basin_fullname)

        data_obs = netCDF4.Dataset(arq1)
        variable_obs = data_obs.variables[param][:].T
        time_obs = data_obs.variables['time']
        st1 = variable_obs[240:600]
        st2 = variable_obs[252:603]
        st3 = variable_obs[600:648]
        stc_obs1.append(st1[1::12])
        stc_obs2.append(st1[2::12])
        stc_obs3.append(st1[3::12])

        sto1 = []
        sto2 = []
        sto3 = []

        link = home + "/io/inmet_ana_chirps/operation/{0}/{1}_thiessen/{2}".format(freq, param, macro_name)
        arq = "{0}/{1}_{2}_inmet_ana_chirps_obs_20150101_20161130_thiessen_{3}.nc".format(link, param, freq,
                                                                                          basin_fullname)
        data_obs2 = netCDF4.Dataset(arq)
        variable_obs2 = data_obs2.variables[param][:].T
        time_obs2 = data_obs2.variables['time']
        st4 = variable_obs2[0:23]
        observ = np.full(71, np.nan)
        observ[0:48] = st3
        observ[48:71] = st4
        sto1.append(observ[1::12])
        sto2.append(observ[2::12])
        sto3.append(observ[3::12])

        # open netcdf hind
        ste_hind1 = []
        ste_hind2 = []
        ste_hind3 = []

        link2 = home + "/io/echam46/hind8110/jan/monthly/{0}_thiessen/{1}".format(param, macro_name)
        for year2 in range(1981, 2011):
            arq2 = "{0}/{1}_{2}_echam46_hind8110_fcst_{3}0101_{3}0201_{3}0430_thiessen_{4}.nc".format(link2, param,
                                                                                                      freq, year2,
                                                                                                      basin_fullname)
            data_hind = netCDF4.Dataset(arq2)
            variable_hind = data_hind.variables[param][:]
            time_hind = data_hind.variables['time']
            ste_hind1.append(variable_hind[0::3])
            ste_hind2.append(variable_hind[1::3])
            ste_hind3.append(variable_hind[2::3])

        # open netcdf fcst
        ste_fcst1 = []
        ste_fcst2 = []
        ste_fcst3 = []

        link3 = home + "/io/echam46/hind8110/jan/monthly/{0}_thiessen/{1}".format(param, macro_name)
        for year3 in range(2011, 2017):
            arq3 = "{0}/{1}_{2}_echam46_hind8110_fcst_{3}0101_{3}0201_{3}0430_thiessen_{4}.nc".format(link3, param,
                                                                                                      freq, year3,
                                                                                                      basin_fullname)
            data_fcst = netCDF4.Dataset(arq3)
            variable_fcst = data_fcst.variables[param][:]
            time_fcst = data_fcst.variables['time']
            ste_fcst1.append(variable_fcst[0::3])
            ste_fcst2.append(variable_fcst[1::3])
            ste_fcst3.append(variable_fcst[2::3])

        q_gamma1 = []
        q_gamma2 = []
        q_gamma3 = []

        dic = home + "/documentos/musf/check_echam46_obs_basins/pr_thiessen/pr_thiessen_monthly_echam46_corrected_operation/gamma/jan/{0}".format(macro_name)

        for yea in range(2011, 2017):
            arq_nc = "{0}/{1}_{2}_echam46_hind8110_fcst_{3}0101_{3}0201_{3}0430_thiessen_corrigido_{4}.nc".format(dic, param, freq, yea, basin_fullname)
            data_fcst_gamma = netCDF4.Dataset(arq_nc)
            variable_fcst_gamma = data_fcst_gamma.variables[param][:]
            time_fcst_gamma = data_fcst_gamma.variables['time']
            q_gamma1.append(variable_fcst_gamma[0::3])
            q_gamma2.append(variable_fcst_gamma[1::3])
            q_gamma3.append(variable_fcst_gamma[2::3])

        # Calculate vies and pr_correction echam46
        print basin_fullname

        data_obs1 = pd.Series(np.squeeze(stc_obs1))
        data_obs2 = pd.Series(np.squeeze(stc_obs2))
        data_obs3 = pd.Series(np.squeeze(stc_obs3))

        obs_1 = pd.Series(np.squeeze(sto1))
        obs_2 = pd.Series(np.squeeze(sto1))
        obs_3 = pd.Series(np.squeeze(sto1))

        data_mod1 = pd.Series(np.squeeze(ste_hind1))
        data_mod2 = pd.Series(np.squeeze(ste_hind2))
        data_mod3 = pd.Series(np.squeeze(ste_hind3))

        data_fcst1 = pd.Series(np.squeeze(ste_fcst1))
        data_fcst2 = pd.Series(np.squeeze(ste_fcst2))
        data_fcst3 = pd.Series(np.squeeze(ste_fcst3))

        gamma1 = pd.Series(np.squeeze(q_gamma1))
        gamma2 = pd.Series(np.squeeze(q_gamma2))
        gamma3 = pd.Series(np.squeeze(q_gamma3))

        b = 10  # Parametro fixo
        best_distr, best_fit_name, param_obs = best_fit_distribution(data_obs3, b)

        param_mod = best_distr.fit(data_mod3)
        arg = param_mod[:-2]
        loc = param_mod[-2]
        scale = param_mod[-1]

        prob = best_distr.cdf(data_fcst3, loc=loc, scale=scale, *arg)
        q_corr = best_distr.ppf(prob, loc=param_obs[-2], scale=param_obs[-1], *param_obs[:-2])

        plt.plot(obs_3, 'r', q_corr, 'b', gamma2, 'k')
        plt.title(u'Pr_Thiessen - viés corrigido - Jan (Fev-Mar-Abr)\n bacia {0} \n {1}'.format(basin_fullname, best_fit_name))
        plt.ylim(0, 700)
        plt.ylabel(u'mm')
        plt.xlabel(u'anos')
        legenda = ('obs_3', 'q_corr', 'q_gamma')
        plt.legend(legenda, frameon=False)

        path_out = ('{0}/check_echam46_obs_basins/pr_thiessen/best_dist_corrected_operation_/figures/'
                    'jan/{1}/'.format(hidropy_path, basin_dict(basin)[1]))
        if not os.path.exists(path_out):
            create_path(path_out)

        plt.savefig(os.path.join(path_out, 'data_fcst3_q_corr_gamma3_pr_thiessen_{0}.png'.format(basin_fullname)))
        plt.close('all')
        plt.cla()