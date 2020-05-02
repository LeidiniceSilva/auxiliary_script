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
import matplotlib.gridspec as gridspec

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

	arq  = '{0}/rol_controle_{1}_g2_mon_1982-2012_r720x360_neb.nc'.format(path, exp)
	data = netCDF4.Dataset(arq)
	var  = data.variables['olr_top'][:]
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
	
	return value, clim_exp, djf_exp, mam_exp, jja_exp, son_exp


def import_obs(path):

	arq  = '{0}/rol_noaa_mon_1982-2012_r720x360_neb.nc'.format(path)
	data = netCDF4.Dataset(arq)
	var  = data.variables['olr'][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)

	clim_obs = []
	for mon in range(1, 12 + 1):
		obs = np.nanmean(value[mon::12], axis=0)
		clim_obs.append(obs)

	seasonal = np.nanmean(np.nanmean(var[:][2:372:3,:,:], axis=1), axis=1)
	mam_obs = seasonal[0:120:4]
	jja_obs = seasonal[1:120:4]
	son_obs = seasonal[2:120:4]
	djf_obs = seasonal[3:120:4]
	
	return value, clim_obs, djf_obs, mam_obs, jja_obs, son_obs
	

# Import exp model end obs database 
home = os.path.expanduser("~")
path = home + "/Documents/ufrn/papers/olam/datas"

x, clim_obs1, djf_obs1, mam_obs1, jja_obs1, son_obs1 = import_obs(path)

exp1  = u'chen'
y, clim_exp1, djf_exp1, mam_exp1, jja_exp1, son_exp1 = import_sim(path, exp1)
		
exp2  = u'harr'
z, clim_exp2, djf_exp2, mam_exp2, jja_exp2, son_exp2 = import_sim(path, exp2)

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

exit()

#~ # Plot histogram and scatter obs x model
#~ fig = plt.figure()

#~ # First subplot
#~ ax1 = fig.add_subplot(311)
#~ ax1.hist(x, color='gray', edgecolor='black')
#~ ax1.set_title(u'Histograma de Precipitação Mensal (1982-2012)', fontweight='bold')
#~ ax1.yaxis.tick_left()
#~ ax1.yaxis.set_label_position('left')
#~ ax1.set_ylim(0,150)

#~ mu = x.mean()
#~ median = np.median(x)
#~ sigma = x.std()
#~ textstr = '\n'.join((r'$\mu=%.2f$' % (mu, ), r'$\mathrm{median}=%.2f$' % (median, ), r'$\sigma=%.2f$' % (sigma, )))
#~ ax1.text(280.,120., 'A) CRU', fontweight='bold', size=8)
#~ ax1.text(280.,60., textstr, size=8)

#~ ax2 = fig.add_subplot(311, sharex=ax1, frameon=False)
#~ sortedtime = np.sort(x)
#~ p = 1. * np.arange(len(x))/(len(x) - 1)
#~ ax2.plot(sortedtime, p, color='red')
#~ ax2.tick_params(axis="y", labelcolor="r")
#~ ax2.yaxis.tick_right()
#~ ax2.yaxis.set_label_position('right')
        
#~ # Second subplot
#~ ax1 = fig.add_subplot(312)
#~ ax1.hist(y, color='gray', edgecolor='black')
#~ ax1.yaxis.tick_left()
#~ ax1.yaxis.set_label_position('left')
#~ ax1.set_ylabel('Frequência', fontweight='bold')
#~ ax1.set_ylim(0,150)

#~ mu = y.mean()
#~ median = np.median(y)
#~ sigma = y.std()
#~ textstr = '\n'.join((r'$\mu=%.2f$' % (mu, ), r'$\mathrm{median}=%.2f$' % (median, ), r'$\sigma=%.2f$' % (sigma, )))
#~ ax1.text(10.,130., 'B) OLAMv.3.3_Chen', fontweight='bold', size=8)
#~ ax1.text(10.,70., textstr, size=8)

#~ ax2 = fig.add_subplot(312, sharex=ax1, frameon=False)
#~ sortedtime = np.sort(y)
#~ p = 1. * np.arange(len(y))/(len(y) - 1)
#~ ax2.plot(sortedtime, p, color='red')
#~ ax2.tick_params(axis="y", labelcolor="r")
#~ ax2.yaxis.tick_right()
#~ ax2.yaxis.set_label_position('right')
#~ ax2.set_ylabel('Função de Distribuição Cumulativa', color='red', fontweight='bold')

#~ # Thirth subplot
#~ ax1 = fig.add_subplot(313)
#~ ax1.hist(z, color='gray', edgecolor='black')
#~ ax1.yaxis.tick_left()
#~ ax1.yaxis.set_label_position('left')
#~ ax1.set_ylim(0,150)

#~ mu = z.mean()
#~ median = np.median(z)
#~ sigma = z.std()
#~ textstr = '\n'.join((r'$\mu=%.2f$' % (mu, ), r'$\mathrm{median}=%.2f$' % (median, ), r'$\sigma=%.2f$' % (sigma, )))
#~ ax1.text(120.,120., 'C) OLAMv.3.3_Harr', fontweight='bold', size=8)
#~ ax1.text(120.,60., textstr, size=8)

#~ ax2 = fig.add_subplot(313, sharex=ax1, frameon=False)
#~ sortedtime = np.sort(z)
#~ p = 1. * np.arange(len(z))/(len(z) - 1)
#~ ax2.plot(sortedtime, p, color='red')
#~ ax2.tick_params(axis="y", labelcolor="r")
#~ ax2.yaxis.tick_right()
#~ ax2.yaxis.set_label_position('right')
#~ ax2.set_xlabel('Precipitação (mm)', fontweight='bold')

#~ path_out = home + "/Documentos/ufrn/papers/olam/results/"
#~ if not os.path.exists(path_out):
	#~ create_path(path_out)
#~ plt.savefig(os.path.join(path_out, 'hist_chen_harr_cru.png'), bbox_inches='tight', dpi=400)
#~ plt.show()
#~ exit()


#~ # Plot climatology obs x model
#~ fig = plt.figure()
#~ plt.subplot(211)

#~ time = np.arange(0.5, 12 + 0.5)
#~ a = plt.plot(time, clim_exp1, time, clim_obs1)

#~ plt.fill_between(time, clim_exp1, clim_obs1, facecolor='slategray', alpha=0.2, interpolate=True)
#~ l1, l2 = a
#~ plt.setp(l1, linewidth=2, markeredgewidth=2, marker='+', color='blue')
#~ plt.setp(l2, linewidth=2, markeredgewidth=2, marker='+', color='red')
#~ plt.title(u'Climatologia de Precipitação (1982-2012)', fontweight='bold')
#~ plt.ylabel(u'Precipitação (mm)', fontweight='bold')
#~ plt.xticks(time, [u'Jan', u'Fev', u'Mar', u'Abr', u'Mai', u'Jun', u'Jul', u'Ago', u'Set', u'Out', u'Nov', u'Dez'])
#~ plt.yticks(np.arange(0, 220, 20))
#~ plt.tick_params(axis='both', which='major', labelsize=8, length=5, width=1.5, pad=5, labelcolor='k')
#~ plt.legend([u'OLAMv.3.3_Chen', u'CRU'], loc='best', ncol=1)
#~ plt.text(0.5, 180., 'A)', fontweight='bold', size=8)
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
#~ plt.legend([u'OLAMv.3.3_Harr', u'CRU'], loc='best', ncol=1)
#~ plt.text(0.5, 180., 'B)', fontweight='bold', size=8)
#~ plt.grid()

#~ path_out = home + "/Documentos/ufrn/papers/olam/results/"
#~ if not os.path.exists(path_out):
	#~ create_path(path_out)
#~ plt.savefig(os.path.join(path_out, 'clim_chen_harr_cru.png'), bbox_inches='tight', dpi=400)
#~ plt.show()
#~ exit()

#~ # Plot time siries whit bar erro
#~ City_A=x
#~ City_B=y

#~ Mean_City_A=np.mean(x)
#~ Mean_City_B=np.mean(y)

#~ STDV_City_A=np.std(x)
#~ STDV_City_B=np.std(y)

#~ # Create a figure with customized size
#~ fig = plt.figure(figsize=(12, 6))
#~ ax = fig.add_subplot(111)

#~ # Set the axis lables
#~ ax.set_xlabel('Anos')
#~ ax.set_ylabel('Precipitação (mm)')

#~ # X axis is day numbers from 1 to 15
#~ xaxis = np.array(range(1,32))

#~ # Line color for error bar
#~ color_City_A = 'red'
#~ color_City_B = 'darkblue'

#~ # Line style for each dataset
#~ lineStyle_City_A={"linestyle":"--", "linewidth":2, "markeredgewidth":2, "elinewidth":2, "capsize":3}
#~ lineStyle_City_B={"linestyle":"-", "linewidth":2, "markeredgewidth":2, "elinewidth":2, "capsize":3}

#~ # Create an error bar for each dataset
#~ line_City_A=ax.errorbar(xaxis, City_A, yerr=STDV_City_A, **lineStyle_City_A, color=color_City_A, label='CRU')
#~ line_City_B=ax.errorbar(xaxis, City_B, yerr=STDV_City_B, **lineStyle_City_A, color=color_City_B, label='OLAM_CHEN')

#~ # Label each dataset on the graph, xytext is the label's position 
#~ for i, txt in enumerate(City_A):
        #~ ax.annotate(txt, xy=(xaxis[i], City_A[i]), xytext=(xaxis[i]+0.03, City_A[i]+0.3),color=color_City_A)

#~ for i, txt in enumerate(City_B):
        #~ ax.annotate(txt, xy=(xaxis[i], City_B[i]), xytext=(xaxis[i]+0.03, City_B[i]+0.3),color=color_City_B)
        

#~ # Draw a legend bar
#~ plt.legend(handles=[line_City_A, line_City_B], loc='upper right')

#~ # Customize the tickes on the graph
#~ plt.xticks(xaxis)               
#~ plt.yticks(np.arange(200, 2200, 100))

#~ # Customize the legend font and handle length
#~ params = {'legend.fontsize': 12,
          #~ 'legend.handlelength': 2}
#~ plt.rcParams.update(params)

#~ # Customize the font
#~ font = {'family' : 'Arial',
        #~ 'weight':'bold',
        #~ 'size'   : 12}
#~ plt.rc('font', **font)

#~ # Draw a grid for the graph
#~ ax.grid(color='lightgrey', linestyle='-')
#~ ax.set_facecolor('w')
#~ plt.show()
#~ exit()

#~ # Boxplot anual cicle obs x model
#~ fig = plt.figure()
#~ time = np.arange(1, 4)
#~ data = [x, y, z]


#~ plt_box = plt.boxplot(data, patch_artist=True, bootstrap=10000, vert=1)

#~ # Change outline and fill color
#~ for box in plt_box['boxes']:
    #~ box.set( color='black', linewidth=2)
    #~ box.set( facecolor = 'gray' )

#~ # Change color and linewidth of the whiskers
#~ for whisker in plt_box['whiskers']:
    #~ whisker.set(color='black', linewidth=2)

#~ # Change color and linewidth of the caps
#~ for cap in plt_box['caps']:
    #~ cap.set(color='black', linewidth=2)

#~ # Change color and linewidth of the medians
#~ for median in plt_box['medians']:
    #~ median.set(color='red', linewidth=2)

#~ # Change the style of fliers and their fill
#~ for flier in plt_box['fliers']:
    #~ flier.set(marker='+', color='black', alpha=4)

#~ mux = x.mean()
#~ medianx = np.median(x)
#~ sigmax = x.std()
#~ textstrx = '\n'.join((r'$\mu=%.2f$' % (mux, ), r'$\mathrm{median}=%.2f$' % (medianx, ), r'$\sigma=%.2f$' % (sigmax, )))
#~ plt.text(.75,1600., 'A) CRU', fontweight='bold')
#~ plt.text(.75,1400., textstrx)

#~ muy = y.mean()
#~ mediany = np.median(y)
#~ sigmay = y.std()
#~ textstry = '\n'.join((r'$\mu=%.2f$' % (muy, ), r'$\mathrm{median}=%.2f$' % (mediany, ), r'$\sigma=%.2f$' % (sigmay, )))
#~ plt.text(1.75,1300., 'B) OLAMv.3.3_Chen', fontweight='bold')
#~ plt.text(1.75,1100., textstry)

#~ muz = z.mean()
#~ medianz = np.median(z)
#~ sigmaz = z.std()
#~ textstrz = '\n'.join((r'$\mu=%.2f$' % (muz, ), r'$\mathrm{median}=%.2f$' % (medianz, ), r'$\sigma=%.2f$' % (sigmaz, )))
#~ plt.text(2.55,1000., 'C) OLAMv.3.3_Harr', fontweight='bold')
#~ plt.text(2.55,800., textstrz)
	    
#~ plt.title(u'Boxplot de Precipitação Anual (1982-2012)', fontweight='bold')
#~ plt.xlabel(u'Observação e Experimentos', fontweight='bold')	
#~ plt.ylabel(u'Precipitação (mm)', fontweight='bold')
#~ plt.yticks(np.arange(350, 1850, 100))
#~ plt.xticks(time, [u'CRU', u'OLAMv.3.3_Chen', u'OLAMv.3.3_Harr'])
#~ plt.tick_params(axis='both', which='major', length=5, width=1.5, pad=5, labelcolor='black')
#~ plt.grid()

#~ path_out = home + "/Documentos/ufrn/papers/olam/results/"
#~ if not os.path.exists(path_out):
	#~ create_path(path_out)
#~ plt.savefig(os.path.join(path_out, 'boxplot_anual_chen_harr_cru.png'), bbox_inches='tight', dpi=400)
#~ plt.show()
#~ exit()    








