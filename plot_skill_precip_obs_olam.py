# -*- coding: utf-8 -*-

__author__ = "Leidinice Silva"
__email__ = "leidinice.silvae@funceme.br"
__date__ = "08/06/2017"
__description__ = " Compute real precipitation of olamv.3.3 model "

import os
import netCDF4
import datetime
import numpy as np
import pandas as pd
import scipy.stats as st
import matplotlib.pyplot as plt

from sklearn import metrics
from matplotlib import pyplot
from matplotlib.font_manager import FontProperties


def filter_nan(s, o):
    data = np.array([s.flatten(), o.flatten()])
    data = np.transpose(data)
    data = data[~np.isnan(data).any(1)]
    return data[:, 0], data[:, 1]


def pc_bias(s, o):
    s, o = filter_nan(s, o)
    return 100.0 * sum(s - o) / sum(o)


def mbe(s, o):
    s, o = filter_nan(s, o)
    return np.mean(s - o)


def mae(s, o):
    s, o = filter_nan(s, o)
    return np.mean(abs(s - o))


def rmse(s, o):
    s, o = filter_nan(s, o)
    return np.sqrt(np.mean((s - o) ** 2))
    
    
def r(s, o):
    s, o = filter_nan(s, o)
    if s.size == 0:
        corr = np.NaN
    else:
        corr = st.pearsonr(s, o)[1]
    return corr


def nse(s, o):
    s, o = filter_nan(s, o)
    return 1 - sum((s - o) ** 2) / sum((o - np.mean(o)) ** 2)
    

def import_sim(path, exp):

	arq  = '{0}/precip_controle_1982_2012_{1}_g2_neb_new_REAL_ok_full_negcor_monsum_noocean.nc'.format(path, exp)
	data = netCDF4.Dataset(arq)
	var  = data.variables['precip'][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	
	value = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)
	clim_exp = []
	for mon in range(1, 12 + 1):
		exp = np.nanmean(value[mon::12], axis=0)
		clim_exp.append(exp)

	seasonal = np.nanmean(np.nanmean(var[:][2:372:3,:,:], axis=1), axis=1)
	mam_exp = seasonal[0:120:4]
	jja_exp = seasonal[1:120:4]
	son_exp = seasonal[2:120:4]
	djf_exp = seasonal[3:120:4]

	um_sim = (value[0:12])
	dois_sim = (value[12:24])
	tres_sim = (value[24:36])
	quatro_sim = (value[36:48])
	cinco_sim = (value[48:60])
	seis_sim = (value[60:72])
	sete_sim = (value[72:84])
	oito_sim = (value[84:96])
	nove_sim = (value[96:108])
	dez_sim = (value[108:120])
	onze_sim = (value[120:132])
	doze_sim = (value[132:144])
	treze_sim = (value[144:156])
	quatorze_sim = (value[156:168])
	quize_sim = (value[168:180])
	dezeseis_sim = (value[180:192])
	dezesete_sim = (value[192:204])
	dezoito_sim = (value[204:216])
	dezenove_sim = (value[216:228])
	vinte_sim = (value[228:240])
	vinte_um_sim = (value[240:252])
	vinte_dois_sim = (value[252:264])
	vinte_tres_sim = (value[264:276])
	vinte_quatro_sim = (value[276:288])
	vinte_cinco_sim = (value[288:300])
	vinte_seis_sim = (value[300:312])
	vinte_sete_sim = (value[312:324])
	vinte_oito_sim = (value[324:336])
	vinte_nove_sim = (value[336:348])
	trinta_sim = (value[348:360])
	trinta_um_sim = (value[360:372])
	
	annual_exp = [um_sim, dois_sim, tres_sim, quatro_sim, cinco_sim, seis_sim, sete_sim, oito_sim, nove_sim, 
	dez_sim, onze_sim, doze_sim, treze_sim, quatorze_sim, quize_sim, dezeseis_sim, dezesete_sim, dezoito_sim,
	dezenove_sim, vinte_sim, vinte_um_sim, vinte_dois_sim, vinte_tres_sim, vinte_quatro_sim,vinte_cinco_sim,
	vinte_seis_sim, vinte_sete_sim, vinte_oito_sim, vinte_nove_sim, trinta_sim,
	trinta_um_sim]
	
	return clim_exp, djf_exp, mam_exp, jja_exp, son_exp, annual_exp


def import_obs(path):

	arq  = '{0}/pr_Amon_CRU-TS3.22_observation_198201-201212_new_mmm_neb.nc'.format(path)
	data = netCDF4.Dataset(arq)
	var  = data.variables['pr'][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)

	clim_obs = []
	for mon in range(1, 12 + 1):
		obs = np.nanmean(value[mon::12], axis=0)
		clim_obs.append(obs)

	seasonal = add(add(var[:][2:372:3,:,:], axis=1), axis=1)
	mam_obs = seasonal[0:120:4]
	jja_obs = seasonal[1:120:4]
	son_obs = seasonal[2:120:4]
	djf_obs = seasonal[3:120:4]
	print(mam_obs)
	print(jja_obs)
	exit()
	
	um_obs = np.nanmean(value[0:12])
	dois_obs = np.nanmean(value[12:24])
	tres_obs = np.nanmean(value[24:36])
	quatro_obs = np.nanmean(value[36:48])
	cinco_obs = np.nanmean(value[48:60])
	seis_obs = np.nanmean(value[60:72])
	sete_obs = np.nanmean(value[72:84])
	oito_obs = np.nanmean(value[84:96])
	nove_obs = np.nanmean(value[96:108])
	dez_obs = np.nanmean(value[108:120])
	onze_obs = np.nanmean(value[120:132])
	doze_obs = np.nanmean(value[132:144])
	treze_obs = np.nanmean(value[144:156])
	quatorze_obs = np.nanmean(value[156:168])
	quize_obs = np.nanmean(value[168:180])
	dezeseis_obs = np.nanmean(value[180:192])
	dezesete_obs = np.nanmean(value[192:204])
	dezoito_obs = np.nanmean(value[204:216])
	dezenove_obs = np.nanmean(value[216:228])
	vinte_obs = np.nanmean(value[228:240])
	vinte_um_obs = np.nanmean(value[240:252])
	vinte_dois_obs = np.nanmean(value[252:264])
	vinte_tres_obs = np.nanmean(value[264:276])
	vinte_quatro_obs = np.nanmean(value[276:288])
	vinte_cinco_obs = np.nanmean(value[288:300])
	vinte_seis_obs = np.nanmean(value[300:312])
	vinte_sete_obs = np.nanmean(value[312:324])
	vinte_oito_obs = np.nanmean(value[324:336])
	vinte_nove_obs = np.nanmean(value[336:348])
	trinta_obs = np.nanmean(value[348:360])
	trinta_um_obs = np.nanmean(value[360:372])
	
	annual_obs = [um_obs, dois_obs, tres_obs, quatro_obs, cinco_obs, seis_obs, sete_obs, oito_obs, nove_obs,
	dez_obs, onze_obs, doze_obs, treze_obs, quatorze_obs, quize_obs, dezeseis_obs, dezesete_obs, dezoito_obs,
	dezenove_obs, vinte_obs, vinte_um_obs, vinte_dois_obs, vinte_tres_obs, vinte_quatro_obs, vinte_cinco_obs,
	vinte_seis_obs, vinte_sete_obs, vinte_oito_obs, vinte_nove_obs, trinta_obs, trinta_um_obs]
	
	return clim_obs, djf_obs, mam_obs, jja_obs, son_obs, annual_obs
	

# Import exp model end obs database 
home = os.path.expanduser("~")
path = home + "/Documents/ufrn/papers/olam/datas"

exp1  = u'chen'
clim_exp1, djf_exp1, mam_exp1, jja_exp1, son_exp1, annual_exp1 = import_sim(path, exp1)
		
exp2  = u'harr'
clim_exp2, djf_exp2, mam_exp2, jja_exp2, son_exp2, annual_exp2 = import_sim(path, exp2)

clim_obs1, djf_obs1, mam_obs1, jja_obs1, son_obs1, annual_obs1 = import_obs(path)

# Calculate statistic index - Chen
pc_bias_djf1 = pc_bias(djf_exp1, djf_obs1)
pc_bias_mam1 = pc_bias(mam_exp1, mam_obs1)
pc_bias_jja1 = pc_bias(jja_exp1, jja_obs1)
pc_bias_son1 = pc_bias(son_exp1, son_obs1)

mbe_djf1 = mbe(djf_exp1, djf_obs1)
mbe_mam1 = mbe(mam_exp1, mam_obs1)
mbe_jja1 = mbe(jja_exp1, jja_obs1)
mbe_son1 = mbe(son_exp1, son_obs1)

mae_djf1 = mae(djf_exp1, djf_obs1)
mae_mam1 = mae(mam_exp1, mam_obs1)
mae_jja1 = mae(jja_exp1, jja_obs1)
mae_son1 = mae(son_exp1, son_obs1)

rmse_djf1 = rmse(djf_exp1, djf_obs1)
rmse_mam1 = rmse(mam_exp1, mam_obs1)
rmse_jja1 = rmse(jja_exp1, jja_obs1)
rmse_son1 = rmse(son_exp1, son_obs1)

r_djf1 = r(djf_exp1, djf_obs1)
r_mam1 = r(mam_exp1, mam_obs1)
r_jja1 = r(jja_exp1, jja_obs1)
r_son1 = r(son_exp1, son_obs1)

nse_djf1 = nse(djf_exp1, djf_obs1)
nse_mam1 = nse(mam_exp1, mam_obs1)
nse_jja1 = nse(jja_exp1, jja_obs1)
nse_son1 = nse(son_exp1, son_obs1)

pc_bias1 = np.array([pc_bias_djf1, pc_bias_mam1, pc_bias_jja1, pc_bias_son1])
mbe1 = np.array([mbe_djf1, mbe_mam1, mbe_jja1, mbe_son1])
mae1 = np.array([mae_djf1, mae_mam1, mae_jja1, mae_son1])
rmse1 = np.array([rmse_djf1, rmse_mam1, rmse_jja1, rmse_son1])
r1 = np.array([r_djf1, r_mam1, r_jja1, r_son1])
nse1 = np.array([nse_djf1, nse_mam1, nse_jja1, nse_son1])

# Calculate statistic index - Chen
pc_bias_djf2 = pc_bias(djf_exp2, djf_obs1)
pc_bias_mam2 = pc_bias(mam_exp2, mam_obs1)
pc_bias_jja2 = pc_bias(jja_exp2, jja_obs1)
pc_bias_son2 = pc_bias(son_exp2, son_obs1)

mbe_djf2 = mbe(djf_exp2, djf_obs1)
mbe_mam2 = mbe(mam_exp2, mam_obs1)
mbe_jja2 = mbe(jja_exp2, jja_obs1)
mbe_son2 = mbe(son_exp2, son_obs1)

mae_djf2 = mae(djf_exp2, djf_obs1)
mae_mam2 = mae(mam_exp2, mam_obs1)
mae_jja2 = mae(jja_exp2, jja_obs1)
mae_son2 = mae(son_exp2, son_obs1)

rmse_djf2 = rmse(djf_exp2, djf_obs1)
rmse_mam2 = rmse(mam_exp2, mam_obs1)
rmse_jja2 = rmse(jja_exp2, jja_obs1)
rmse_son2 = rmse(son_exp2, son_obs1)

r_djf2 = r(djf_exp2, djf_obs1)
r_mam2 = r(mam_exp2, mam_obs1)
r_jja2 = r(jja_exp2, jja_obs1)
r_son2 = r(son_exp2, son_obs1)

nse_djf2 = nse(djf_exp2, djf_obs1)
nse_mam2 = nse(mam_exp2, mam_obs1)
nse_jja2 = nse(jja_exp2, jja_obs1)
nse_son2 = nse(son_exp2, son_obs1)

pc_bias2 = np.array([pc_bias_djf2, pc_bias_mam2, pc_bias_jja2, pc_bias_son2])
mbe2 = np.array([mbe_djf2, mbe_mam2, mbe_jja2, mbe_son2])
mae2 = np.array([mae_djf2, mae_mam2, mae_jja2, mae_son2])
rmse2 = np.array([rmse_djf2, rmse_mam2, rmse_jja2, rmse_son2])
r2 = np.array([r_djf2, r_mam2, r_jja2, r_son2])
nse2 = np.array([nse_djf2, nse_mam2, nse_jja2, nse_son2])

# Print statistic index (Chen and Harr)
print(pc_bias1)
print(mbe1)
print(mae1)
print(rmse1)
print(r1)
print(nse1)
print()

print(pc_bias2)
print(mbe2)
print(mae2)
print(rmse2)
print(r2)
print(nse2)
print()

#~ # Plot climatology obs x model
#~ fig = plt.figure(1)
#~ plt.subplot(211)

#~ time = np.arange(0.5, 12 + 0.5)
#~ a = plt.plot(time, clim_exp1, time, clim_obs1)

#~ plt.fill_between(time, clim_exp1, clim_obs1, facecolor='slategray', alpha=0.2, interpolate=True)
#~ l1, l2 = a
#~ plt.setp(l1, linewidth=2, markeredgewidth=2, marker='+', color='blue')
#~ plt.setp(l2, linewidth=2, markeredgewidth=2, marker='+', color='red')
#~ plt.title(u'Climatologia de Precipitação 1982-2012', fontweight='bold')
#~ plt.ylabel(u'Precipitação (mm)', fontweight='bold')
#~ plt.xticks(time, [u'Jan', u'Fev', u'Mar', u'Abr', u'Mai', u'Jun', u'Jul', u'Ago', u'Set', u'Out', u'Nov', u'Dez'])
#~ plt.yticks(np.arange(0, 220, 20))
#~ plt.tick_params(axis='both', which='major', labelsize=8, length=5, width=1.5, pad=5, labelcolor='k')
#~ plt.legend([u'OLAMv.3.3_Chen', u'CRU'], loc='best', ncol=2)
#~ plt.grid()

#~ plt.subplot(212)
#~ time = np.arange(0.5, 12 + 0.5)
#~ a = plt.plot(time, clim_exp2, time, clim_obs1)

#~ plt.fill_between(time, clim_exp2, clim_obs1, facecolor='slategray', alpha=0.2, interpolate=True)
#~ l1, l2 = a
#~ plt.setp(l1, linewidth=2, markeredgewidth=2, marker='+', color='green')
#~ plt.setp(l2, linewidth=2, markeredgewidth=2, marker='+', color='red')
#~ plt.xlabel(u'Meses', fontweight='bold')
#~ plt.ylabel(u'Precipitação (mm)', fontweight='bold')
#~ plt.xticks(time, [u'Jan', u'Fev', u'Mar', u'Abr', u'Mai', u'Jun', u'Jul', u'Ago', u'Set', u'Out', u'Nov', u'Dez'])
#~ plt.yticks(np.arange(0, 220, 20))
#~ plt.tick_params(axis='both', which='major', labelsize=8, length=5, width=1.5, pad=5, labelcolor='k')
#~ plt.legend([u'OLAMv.3.3_Harr', u'CRU'], loc='best', ncol=2)
#~ plt.grid()

#~ path_out = home + "/Documents/ufrn/papers/olam/results/"
#~ if not os.path.exists(path_out):
	#~ create_path(path_out)
#~ plt.savefig(os.path.join(path_out, 'clim_chen_harr_cru.png'), bbox_inches='tight', dpi=400)
#~ exit()

# Boxplot anual cicle obs x model
fig = plt.figure()

plt.subplot(211)
time = np.arange(1, 32)
a = plt.plot(time, annual_obs1)
plt.setp(a, linewidth=2, markeredgewidth=2, marker='+', color='red')

plt_bp = plt.boxplot(annual_exp1, patch_artist=True, bootstrap=10000, vert=1)
# Change outline and fill color
for box in plt_bp['boxes']:
    box.set( color='black', linewidth=2)
    box.set( facecolor = 'gray' )

# Change color and linewidth of the whiskers
for whisker in plt_bp['whiskers']:
    whisker.set(color='blue', linewidth=2)

# Change color and linewidth of the caps
for cap in plt_bp['caps']:
    cap.set(color='black', linewidth=2)

# Change color and linewidth of the medians
for median in plt_bp['medians']:
    median.set(color='black', linewidth=2)

# Change the style of fliers and their fill
for flier in plt_bp['fliers']:
    flier.set(marker='+', color='black', alpha=1)

plt.title(u'Boxplot de precipitação anual 1982-2012', fontweight='bold')
plt.ylabel(u'Precipitação (mm)', fontweight='bold')
objects = [u'1982', u'1984', u'1986', u'1988', u'1990', u'1992', 
           u'1994', u'1996', u'1998', u'2000', u'2002', u'2004',
           u'2006', u'2008', u'2010', u'2012']
plt.xticks(np.arange(1, 32, 2), objects)
plt.yticks(np.arange(0, 220, 20))
plt.tick_params(axis='both', which='major', labelsize=8, length=5, width=1.5, pad=5, labelcolor='k')
plt.legend([u'CRU', u'OLAMv.3.3_Chen'], loc='best', ncol=2)
plt.grid()

plt.subplot(212)
time = np.arange(1, 32)
a = plt.plot(time, annual_obs1)
plt.setp(a, linewidth=2, markeredgewidth=2, marker='+', color='red')
b = plt.xlim([1, 32])

plt_bp = plt.boxplot(annual_exp2, patch_artist=True, bootstrap=10000, vert=1)
# Change outline and fill color
for box in plt_bp['boxes']:
    box.set( color='black', linewidth=2)
    box.set( facecolor = 'gray' )

# Change color and linewidth of the whiskers
for whisker in plt_bp['whiskers']:
    whisker.set(color='green', linewidth=2)

# Change color and linewidth of the caps
for cap in plt_bp['caps']:
    cap.set(color='green', linewidth=2)

# Change color and linewidth of the medians
for median in plt_bp['medians']:
    median.set(color='black', linewidth=2)

# Change the style of fliers and their fill
for flier in plt_bp['fliers']:
    flier.set(marker='o', color='green', alpha=1)

plt.xlabel(u'Anos', fontweight='bold')    
plt.ylabel(u'Precipitação (mm)', fontweight='bold')
objects = [u'1982', u'1984', u'1986', u'1988', u'1990', u'1992', 
           u'1994', u'1996', u'1998', u'2000', u'2002', u'2004',
           u'2006', u'2008', u'2010', u'2012']
plt.xticks(np.arange(1, 32, 2), objects)
plt.yticks(np.arange(0, 220, 20))
plt.tick_params(axis='both', which='major', labelsize=8, length=5, width=1.5, pad=5, labelcolor='k')
plt.legend([u'CRU', u'OLAMv.3.3_Harr'], loc='best', ncol=2)
plt.grid()

path_out = home + "/Documents/ufrn/papers/olam/results/"
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, 'boxplot_anual_chen_harr_cru.png'), bbox_inches='tight', dpi=400)
exit()






