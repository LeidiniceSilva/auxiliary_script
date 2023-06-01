# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 01, 2023"
__description__ = "This script plot rank of cmip6 models"

import os
import netCDF4
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from dict_cmip6_models_name import cmip6
from comp_stats_metrics import compute_mbe
from comp_stats_metrics import compute_rmse
from comp_stats_metrics import compute_tss
from comp_stats_metrics import compute_pcc
from comp_stats_metrics import compute_ivs


def import_obs_latlon(param, area, period, date):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/obs'
	arq   = '{0}/{1}_{2}_BR-DWGD_UFES_UTEXAS_v_3.0_{3}_{4}_lonlat.nc'.format(path, param, area, period, date)	
		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	fld_mean = np.nanmean(value, axis=0)
	
	latlon = []
	for i in range(0, fld_mean.shape[0]):
		for ii in fld_mean[i]:
			latlon.append(ii)
	ts_latlon = np.array(latlon)
	
	return ts_latlon
	
	
def import_obs_mon(param, area, period, date):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/obs'
	arq   = '{0}/{1}_{2}_BR-DWGD_UFES_UTEXAS_v_3.0_{3}_{4}_lonlat.nc'.format(path, param, area, period, date)	
		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	mon_mean = np.nanmean(np.nanmean(value, axis=1), axis=1)
	
	ts_mon = []
	for i in range(0, 11 + 1):
		clim = np.nanmean(mon_mean[i::12], axis=0)
		ts_mon.append(clim)
	
	return ts_mon


def import_obs_ann(param, area, period, date):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/obs'
	arq   = '{0}/{1}_{2}_BR-DWGD_UFES_UTEXAS_v_3.0_{3}_{4}_lonlat.nc'.format(path, param, area, period, date)	
		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	ts_ann = np.nanmean(np.nanmean(value, axis=1), axis=1)
	
	return ts_ann
	
	
def import_cmip_latlon(param, area, model, exp, period, date):

	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/cmip6'
	arq   = '{0}/{1}_{2}_{3}_historical_{4}_{5}_{6}_lonlat.nc'.format(path, param, area, model, exp, period, date)	
				
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	fld_mean = np.nanmean(value, axis=0)
	
	latlon = []
	for i in range(0, fld_mean.shape[0]):
		for ii in fld_mean[i]:
			latlon.append(ii)
	ts_latlon = np.array(latlon)
	
	return ts_latlon
	              

def import_cmip_mon(param, area, model, exp, period, date):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/cmip6'
	arq   = '{0}/{1}_{2}_{3}_historical_{4}_{5}_{6}_lonlat.nc'.format(path, param, area, model, exp, period, date)	
				
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	mon_mean = np.nanmean(np.nanmean(value, axis=1), axis=1)
	
	ts_mon = []
	for i in range(0, 11 + 1):
		clim = np.nanmean(mon_mean[i::12], axis=0)
		ts_mon.append(clim)
	
	return ts_mon
	

def import_cmip_ann(param, area, model, exp, period, date):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/cmip6'
	arq   = '{0}/{1}_{2}_{3}_historical_{4}_{5}_{6}_lonlat.nc'.format(path, param, area, model, exp, period, date)	
				
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	ts_ann = np.nanmean(np.nanmean(value, axis=1), axis=1)
	
	return ts_ann
	

def sort_list(data_list):
	
	li = []
	for i in range(len(data_list)):
		li.append([data_list[i], i])
	  
	li.sort()
	sort_index = []
	for x in li:
		sort_index.append(x[1])
	
	return sort_index
	
	
# Import cmip models and obs database 
idx = 'ivs'
var_obs = 'Tmin'
var_cmip6 = 'tasmin'
dt = '1986-2005'

namz_obs_latlon = import_obs_latlon(var_obs, 'NAMZ', 'ANN', dt)
samz_obs_latlon = import_obs_latlon(var_obs, 'SAMZ', 'ANN', dt)
neb_obs_latlon = import_obs_latlon(var_obs, 'NEB', 'ANN', dt)
sam_obs_latlon = import_obs_latlon(var_obs, 'SAM', 'ANN', dt)
lpb_obs_latlon = import_obs_latlon(var_obs, 'LPB', 'ANN', dt)
br_obs_latlon = import_obs_latlon(var_obs, 'BR', 'ANN', dt)

namz_obs_mon_ts = import_obs_mon(var_obs, 'NAMZ', 'MON', dt)
samz_obs_mon_ts = import_obs_mon(var_obs, 'SAMZ', 'MON', dt)
neb_obs_mon_ts  = import_obs_mon(var_obs, 'NEB', 'MON', dt)
sam_obs_mon_ts = import_obs_mon(var_obs, 'SAM', 'MON', dt)
lpb_obs_mon_ts = import_obs_mon(var_obs, 'LPB', 'MON', dt)
br_obs_mon_ts = import_obs_mon(var_obs, 'BR', 'MON', dt)

namz_obs_ann_ts = import_obs_ann(var_obs, 'NAMZ', 'ANN', dt)
samz_obs_ann_ts = import_obs_ann(var_obs, 'SAMZ', 'ANN', dt)
neb_obs_ann_ts  = import_obs_ann(var_obs, 'NEB', 'ANN', dt)
sam_obs_ann_ts = import_obs_ann(var_obs, 'SAM', 'ANN', dt)
lpb_obs_ann_ts = import_obs_ann(var_obs, 'LPB', 'ANN', dt)
br_obs_ann_ts = import_obs_ann(var_obs, 'BR', 'ANN', dt)

mbe_namz_cmip6 = []
rmse_namz_cmip6 = []
tss_namz_cmip6 = []
pcc_namz_cmip6 = []
ivs_namz_cmip6 = []

mbe_samz_cmip6 = []
rmse_samz_cmip6 = []
tss_samz_cmip6 = []
pcc_samz_cmip6 = []
ivs_samz_cmip6 = []

mbe_neb_cmip6 = []
rmse_neb_cmip6 = []
tss_neb_cmip6 = []
pcc_neb_cmip6 = []
ivs_neb_cmip6 = []

mbe_sam_cmip6 = []
rmse_sam_cmip6 = []
tss_sam_cmip6 = []
pcc_sam_cmip6 = []
ivs_sam_cmip6 = []

mbe_lpb_cmip6 = []
rmse_lpb_cmip6 = []
tss_lpb_cmip6 = []
pcc_lpb_cmip6 = []
ivs_lpb_cmip6 = []

mbe_br_cmip6 = []
rmse_br_cmip6 = []
tss_br_cmip6 = []
pcc_br_cmip6 = []
ivs_br_cmip6 = []

legend = []

for i in range(1, 18):

	namz_cmip_latlon = import_cmip_latlon(var_cmip6, 'NAMZ', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	samz_cmip_latlon = import_cmip_latlon(var_cmip6, 'SAMZ', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	neb_cmip_latlon = import_cmip_latlon(var_cmip6, 'NEB', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	sam_cmip_latlon = import_cmip_latlon(var_cmip6, 'SAM', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	lpb_cmip_latlon = import_cmip_latlon(var_cmip6, 'LPB', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	br_cmip_latlon = import_cmip_latlon(var_cmip6, 'BR', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	
	namz_cmip_mon_ts = import_cmip_mon(var_cmip6, 'NAMZ', cmip6[i][0], cmip6[i][1], 'MON', dt)
	samz_cmip_mon_ts = import_cmip_mon(var_cmip6, 'SAMZ', cmip6[i][0], cmip6[i][1], 'MON', dt)
	neb_cmip_mon_ts = import_cmip_mon(var_cmip6, 'NEB', cmip6[i][0], cmip6[i][1], 'MON', dt)
	sam_cmip_mon_ts = import_cmip_mon(var_cmip6, 'SAM', cmip6[i][0], cmip6[i][1], 'MON', dt)
	lpb_cmip_mon_ts = import_cmip_mon(var_cmip6, 'LPB', cmip6[i][0], cmip6[i][1], 'MON', dt)
	br_cmip_mon_ts = import_cmip_mon(var_cmip6, 'BR', cmip6[i][0], cmip6[i][1], 'MON', dt)
	
	namz_cmip_ann_ts = import_cmip_ann(var_cmip6, 'NAMZ', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	samz_cmip_ann_ts = import_cmip_ann(var_cmip6, 'SAMZ', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	neb_cmip_ann_ts = import_cmip_ann(var_cmip6, 'NEB', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	sam_cmip_ann_ts = import_cmip_ann(var_cmip6, 'SAM', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	lpb_cmip_ann_ts = import_cmip_ann(var_cmip6, 'LPB', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	br_cmip_ann_ts = import_cmip_ann(var_cmip6, 'BR', cmip6[i][0], cmip6[i][1], 'ANN', dt)
			
	mbe_namz_cmip6.append(compute_mbe(namz_cmip_latlon, namz_obs_latlon))
	rmse_namz_cmip6.append(compute_rmse(namz_cmip_latlon, namz_obs_latlon))	
	tss_namz_cmip6.append(compute_tss(namz_obs_latlon, namz_cmip_latlon))
	pcc_namz_cmip6.append(compute_pcc(namz_obs_mon_ts, namz_cmip_mon_ts))
	ivs_namz_cmip6.append(compute_ivs(namz_obs_ann_ts, namz_cmip_ann_ts))
	
	mbe_samz_cmip6.append(compute_mbe(samz_cmip_latlon, samz_obs_latlon))
	rmse_samz_cmip6.append(compute_rmse(samz_cmip_latlon, samz_obs_latlon))
	tss_samz_cmip6.append(compute_tss(samz_obs_latlon, samz_cmip_latlon))
	pcc_samz_cmip6.append(compute_pcc(samz_obs_mon_ts, samz_cmip_mon_ts))
	ivs_samz_cmip6.append(compute_ivs(samz_obs_ann_ts, samz_cmip_ann_ts))
		
	mbe_neb_cmip6.append(compute_mbe(neb_cmip_latlon, neb_obs_latlon))
	rmse_neb_cmip6.append(compute_rmse(neb_cmip_latlon, neb_obs_latlon))
	tss_neb_cmip6.append(compute_tss(neb_obs_latlon, neb_cmip_latlon))
	pcc_neb_cmip6.append(compute_pcc(neb_obs_mon_ts, neb_cmip_mon_ts))
	ivs_neb_cmip6.append(compute_ivs(neb_obs_ann_ts, neb_cmip_ann_ts))
		
	mbe_sam_cmip6.append(compute_mbe(sam_cmip_latlon, sam_obs_latlon))
	rmse_sam_cmip6.append(compute_rmse(sam_cmip_latlon, sam_obs_latlon))
	tss_sam_cmip6.append(compute_tss(sam_obs_latlon, sam_cmip_latlon))
	pcc_sam_cmip6.append(compute_pcc(sam_obs_mon_ts, sam_cmip_mon_ts))
	ivs_sam_cmip6.append(compute_ivs(sam_obs_ann_ts, sam_cmip_ann_ts))
		
	mbe_lpb_cmip6.append(compute_mbe(lpb_cmip_latlon, lpb_obs_latlon))
	rmse_lpb_cmip6.append(compute_rmse(lpb_cmip_latlon, lpb_obs_latlon))
	tss_lpb_cmip6.append(compute_tss(lpb_obs_latlon, lpb_cmip_latlon))
	pcc_lpb_cmip6.append(compute_pcc(lpb_obs_mon_ts, lpb_cmip_mon_ts))
	ivs_lpb_cmip6.append(compute_ivs(lpb_obs_ann_ts, lpb_cmip_ann_ts))
	
	mbe_br_cmip6.append(compute_mbe(br_cmip_latlon, br_obs_latlon))
	rmse_br_cmip6.append(compute_rmse(br_cmip_latlon, br_obs_latlon))
	tss_br_cmip6.append(compute_tss(br_obs_latlon, br_cmip_latlon))
	pcc_br_cmip6.append(compute_pcc(br_obs_mon_ts, br_cmip_mon_ts))
	ivs_br_cmip6.append(compute_ivs(br_obs_ann_ts, br_cmip_ann_ts))

	legend.append(cmip6[i][0])

if idx == 'mbe':
	namz_cmip6 = mbe_namz_cmip6
	samz_cmip6 = mbe_samz_cmip6
	neb_cmip6 = mbe_neb_cmip6
	sam_cmip6 = mbe_sam_cmip6
	lpb_cmip6 = mbe_lpb_cmip6
	br_cmip6 = mbe_br_cmip6
	idx_label = 'MBE'
elif idx == 'rmse':
	namz_cmip6 = rmse_namz_cmip6
	samz_cmip6 = rmse_samz_cmip6
	neb_cmip6 = rmse_neb_cmip6
	sam_cmip6 = rmse_sam_cmip6
	lpb_cmip6 = rmse_lpb_cmip6
	br_cmip6 = rmse_br_cmip6
	idx_label = 'RMSE'
elif idx == 'tss':
	namz_cmip6 = tss_namz_cmip6
	samz_cmip6 = tss_samz_cmip6
	neb_cmip6 = tss_neb_cmip6
	sam_cmip6 = tss_sam_cmip6
	lpb_cmip6 = tss_lpb_cmip6
	br_cmip6 = tss_br_cmip6
	idx_label = 'TSS'
elif idx == 'pcc':
	namz_cmip6 = pcc_namz_cmip6
	samz_cmip6 = pcc_samz_cmip6
	neb_cmip6 = pcc_neb_cmip6
	sam_cmip6 = pcc_sam_cmip6
	lpb_cmip6 = pcc_lpb_cmip6
	br_cmip6 = pcc_br_cmip6
	idx_label = 'PCC'
else:
	namz_cmip6 = ivs_namz_cmip6
	samz_cmip6 = ivs_samz_cmip6
	neb_cmip6 = ivs_neb_cmip6
	sam_cmip6 = ivs_sam_cmip6
	lpb_cmip6 = ivs_lpb_cmip6
	br_cmip6 = ivs_br_cmip6
	idx_label = 'IVS'
		
sort_list_namz = sort_list(namz_cmip6)
model_list_namz = []
value_list_namz = []
for i in sort_list_namz:
	model_list_namz.append(cmip6[i+1][0])
	value_list_namz.append(namz_cmip6[i])

sort_list_samz = sort_list(samz_cmip6)
model_list_samz = []
value_list_samz = []
for ii in sort_list_samz:
	model_list_samz.append(cmip6[ii+1][0])
	value_list_samz.append(samz_cmip6[ii])

sort_list_neb = sort_list(neb_cmip6)
model_list_neb = []
value_list_neb = []
for iii in sort_list_neb:
	model_list_neb.append(cmip6[iii+1][0])
	value_list_neb.append(neb_cmip6[iii])

sort_list_sam = sort_list(sam_cmip6)
model_list_sam = []
value_list_sam = []
for iv in sort_list_sam:
	model_list_sam.append(cmip6[iv+1][0])
	value_list_sam.append(sam_cmip6[iv])

sort_list_lpb = sort_list(lpb_cmip6)
model_list_lpb = []
value_list_lpb = []
for v in sort_list_lpb:
	model_list_lpb.append(cmip6[v+1][0])
	value_list_lpb.append(lpb_cmip6[v])

sort_list_br = sort_list(br_cmip6)
model_list_br = []
value_list_br = []
for vi in sort_list_br:
	model_list_br.append(cmip6[vi+1][0])
	value_list_br.append(br_cmip6[vi])

# Plot cmip models and obs database 
fig = plt.figure(figsize=(10, 8))

if var_cmip6 == 'pr':
	color = 'blue'
else:
	color = 'red'
	
ax = fig.add_subplot(3, 2, 1)  
ax.barh(model_list_namz, value_list_namz, color=color, edgecolor='white')
plt.title(u'(a) NAMZ', loc='left', fontsize=8, fontweight='bold')
plt.yticks(fontsize=7)
plt.xticks(fontsize=7)
ax.xaxis.set_tick_params(pad=-5)
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('none')
ax.grid(b=True, color ='gray', linestyle='--', linewidth=0.5, alpha = 0.2)
for s in ['top', 'bottom', 'left', 'right']:
    ax.spines[s].set_visible(False)
for i in ax.patches:
    plt.text(i.get_width()+0.02, i.get_y()+0.5, str(round((i.get_width()), 2)), fontsize=7, fontweight='bold', color='gray')
if idx == 'rmse':
	ax.invert_yaxis()
elif idx == 'ivs':
	ax.invert_yaxis()
	             
ax = fig.add_subplot(3, 2, 2)  
ax.barh(model_list_samz, value_list_samz, color=color, edgecolor='white')
plt.title(u'(b) SAMZ', loc='left', fontsize=8, fontweight='bold')
plt.yticks(fontsize=7)
plt.xticks(fontsize=7)
ax.xaxis.set_tick_params(pad=-5)
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('none')
ax.grid(b=True, color ='gray', linestyle='--', linewidth=0.5, alpha = 0.2)
for s in ['top', 'bottom', 'left', 'right']:
    ax.spines[s].set_visible(False)
for i in ax.patches:
    plt.text(i.get_width()+0.02, i.get_y()+0.5, str(round((i.get_width()), 2)), fontsize=7, fontweight='bold', color='gray')
if idx == 'rmse':
	ax.invert_yaxis()
elif idx == 'ivs':
	ax.invert_yaxis()
	        
ax = fig.add_subplot(3, 2, 3)  
ax.barh(model_list_neb, value_list_neb, color=color, edgecolor='white')
plt.title(u'(c) NEB', loc='left', fontsize=8, fontweight='bold')
plt.yticks(fontsize=7)
plt.xticks(fontsize=7)
ax.xaxis.set_tick_params(pad=-5)
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('none')
ax.grid(b=True, color ='gray', linestyle='--', linewidth=0.5, alpha = 0.2)
for s in ['top', 'bottom', 'left', 'right']:
    ax.spines[s].set_visible(False)
for i in ax.patches:
    plt.text(i.get_width()+0.02, i.get_y()+0.5, str(round((i.get_width()), 2)), fontsize=7, fontweight='bold', color='gray')
if idx == 'rmse':
	ax.invert_yaxis()
elif idx == 'ivs':
	ax.invert_yaxis()
	        
ax = fig.add_subplot(3, 2, 4)  
ax.barh(model_list_sam, value_list_sam, color=color, edgecolor='white')
plt.title(u'(d) SAM', loc='left', fontsize=8, fontweight='bold')
plt.yticks(fontsize=7)
plt.xticks(fontsize=7)
ax.xaxis.set_tick_params(pad=-5)
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('none')
ax.grid(b=True, color ='gray', linestyle='--', linewidth=0.5, alpha = 0.2)
for s in ['top', 'bottom', 'left', 'right']:
    ax.spines[s].set_visible(False)
for i in ax.patches:
    plt.text(i.get_width()+0.02, i.get_y()+0.5, str(round((i.get_width()), 2)), fontsize=7, fontweight='bold', color='gray')
if idx == 'rmse':
	ax.invert_yaxis()
elif idx == 'ivs':
	ax.invert_yaxis()
	        
ax = fig.add_subplot(3, 2, 5)  
ax.barh(model_list_lpb, value_list_lpb, color=color, edgecolor='white')
plt.title(u'(e) LPB', loc='left', fontsize=8, fontweight='bold')
plt.xlabel('{0}'.format(idx_label), fontsize=8, fontweight='bold')
plt.yticks(fontsize=7)
plt.xticks(fontsize=7)
ax.xaxis.set_tick_params(pad=-5)
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('none')
ax.grid(b=True, color ='gray', linestyle='--', linewidth=0.5, alpha = 0.2)
for s in ['top', 'bottom', 'left', 'right']:
    ax.spines[s].set_visible(False)
for i in ax.patches:
    plt.text(i.get_width()+0.02, i.get_y()+0.5, str(round((i.get_width()), 2)), fontsize=7, fontweight='bold', color='gray')
if idx == 'rmse':
	ax.invert_yaxis()
elif idx == 'ivs':
	ax.invert_yaxis()
	
ax = fig.add_subplot(3, 2, 6)  
ax.barh(model_list_br, value_list_br, color=color, edgecolor='white')
plt.title(u'(f) BR', loc='left', fontsize=8, fontweight='bold')
plt.xlabel('{0}'.format(idx_label), fontsize=8, fontweight='bold')
plt.yticks(fontsize=7)
plt.xticks(fontsize=7)
ax.xaxis.set_tick_params(pad=-5)
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('none')
ax.grid(b=True, color ='gray', linestyle='--', linewidth=0.5, alpha = 0.2)
for s in ['top', 'bottom', 'left', 'right']:
    ax.spines[s].set_visible(False)
for i in ax.patches:
    plt.text(i.get_width()+0.02, i.get_y()+0.5, str(round((i.get_width()), 2)), fontsize=7, fontweight='bold', color='gray')
if idx == 'rmse':
	ax.invert_yaxis()
elif idx == 'ivs':
	ax.invert_yaxis()

plt.subplots_adjust(wspace=0.35)
                    	                
# Path out to save figure
path_out = '/home/nice/Documentos/AdaptaBrasil_MCTI/figs/figs_report-II'
name_out = 'pyplt_rank_{0}_cmip6_{1}_{2}.png'.format(idx, var_cmip6, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()

