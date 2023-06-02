# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jun 01, 2023"
__description__ = "This script plot statistical rank of cmip6 models"

import os
import netCDF4
import numpy as np
import numpy.ma as ma
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors


def import_obs(param, area):
	
	path  = '/home/nice/Documentos/paper_mari/database/obs'
	arq   = '{0}/{1}_19860101_20050131_BR-DWGD_UFES_UTEXAS_v_3.0_0.5.nc_{2}.nc'.format(path, param, area)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	value = var[:][:,:,:]

	fld_mean = np.nanmean(value, axis=0)	
	latlon = []
	for i in range(0, fld_mean.shape[0]):
		for ii in fld_mean[i]:
			latlon.append(ii)
			
	ts_latlon = np.array(latlon)
	ts_time = np.nanmean(np.nanmean(value, axis=1), axis=1)
	
	return ts_latlon, ts_time

	
def import_cmip(param, model, area):
	
	if model == 'CMCC-ESM2':
		dict_var = {u'cddETCCDI': u'cdd',
		u'r95pETCCDI': u'r95p',
		u'rx5dayETCCDI': u'rx5day',
		u'tx90pETCCDI': u'tx90p',
		u'wsdiETCCDI': u'wsdi'}
	else:
		dict_var = {u'cddETCCDI': u'cddETCCDI',
		u'r95pETCCDI': u'r95pETCCDI',
		u'rx5dayETCCDI': u'rx5dayETCCDI',
		u'tx90pETCCDI': u'tx90pETCCDI',
		u'wsdiETCCDI': u'wsdiETCCDI'}
	
	path  = '/home/nice/Documentos/paper_mari/database/cmip6'
	arq   = '{0}/{1}_yr_{2}_historical_1986-2005.nc_{3}.nc'.format(path, param, model, area)		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[dict_var[param]][:] 
	value = var[:][:,:,:]
	
	fld_mean = np.nanmean(value, axis=0)
	latlon = []
	for i in range(0, fld_mean.shape[0]):
		for ii in fld_mean[i]:
			latlon.append(ii)
			
	ts_latlon = np.array(latlon)
	ts_time = np.nanmean(np.nanmean(value, axis=1), axis=1)
	
	return ts_latlon, ts_time


def compute_ree(model, obs):

    p1 = np.nanmean(np.array(model) - np.array(obs))
    p2 = np.nanmean(np.array(model))
    p3 = p1/p2
    ree = p3 * 100.0
    
    return ree
    
    
def compute_pcc(obs, model):
   
    pcc = ma.corrcoef(ma.masked_invalid(obs), ma.masked_invalid(model))[0][1]
    
    return pcc
           
               
def compute_ivs(obs, model):

    p1 = np.nanstd(obs, ddof=0)
    p2 = np.nanstd(model, ddof=0)
    p3 = p2 / p1
    p4 = p1 / p2
    ivs = (p3 - p4)**2  
    
    return ivs   
    

def sort_list_x(data_list):
	
	li = []
	for i in range(len(data_list)):
		
		if data_list[i] < 0:
			data_list[i] = data_list[i]*(-1)
		
		li.append([data_list[i], i])
	
	li.sort()
	sort_index = []
	for x in li:
		sort_index.append(x[1])

	model_list = []
	value_list = []
	for ii in sort_index:
		model_list.append(cmip6_i[ii+1][0])
		value_list.append(data_list[ii])
	
	data_argsort = np.argsort(model_list)

	data_argsort_i = []
	for idx in data_argsort:
		data_argsort_i.append(idx+1)

	return data_argsort_i


def sort_list_y(data_list):
	
	li = []
	for i in range(len(data_list)):
		li.append([data_list[i], i])
	
	li.sort(reverse=True)
	sort_index = []
	for x in li:
		sort_index.append(x[1])

	model_list = []
	value_list = []
	for ii in sort_index:
		model_list.append(cmip6_i[ii+1][0])
		value_list.append(data_list[ii])

	data_argsort = np.argsort(model_list)
	
	data_argsort_i = []
	for idx in data_argsort:
		data_argsort_i.append(idx+1)
	
	return data_argsort_i
		

def sort_list_z(data_list):
	
	li = []
	for i in range(len(data_list)):		
		li.append([data_list[i], i])
	
	li.sort()
	sort_index = []
	for x in li:
		sort_index.append(x[1])

	model_list = []
	value_list = []
	for ii in sort_index:
		model_list.append(cmip6_i[ii+1][0])
		value_list.append(data_list[ii])
	
	data_argsort = np.argsort(model_list)

	data_argsort_i = []
	for idx in data_argsort:
		data_argsort_i.append(idx+1)

	return data_argsort_i


cmip6 = {1	:['ACCESS-CM2'],
		 2	:['BCC-CSM2-MR'],
		 3	:['CanESM5'],
		 4	:['CMCC-ESM2'],
		 5	:['CNRM-CM6-1'],
		 6	:['CNRM-ESM2-1'],
		 7	:['GFDL-ESM4'],
		 8	:['INM-CM4-8'],
		 9	:['INM-CM5-0'],
		 10	:['KIOST-ESM'],
		 11	:['MIROC6'],
		 12 :['MIROC-ES2L'],
		 13	:['MPI-ESM1-2-HR'],
		 14	:['MPI-ESM1-2-LR'],
		 15	:['MRI-ESM2-0'],
		 16	:['NESM3'],
		 17	:['NorESM2-MM']}

cmip6_i = {1 :['ACCESS-CM2'],
		 2	:['BCC-CSM2-MR'],
		 3	:['CANESM5'],
		 4	:['CMCC-ESM2'],
		 5	:['CNRM-CM6-1'],
		 6	:['CNRM-ESM2-1'],
		 7	:['GFDL-ESM4'],
		 8	:['INM-CM4-8'],
		 9	:['INM-CM5-0'],
		 10	:['KIOST-ESM'],
		 11	:['MIROC6'],
		 12 :['MIROCES2L'],
		 13	:['MPI-ESM1-2-HR'],
		 14	:['MPI-ESM1-2-LR'],
		 15	:['MRI-ESM2-0'],
		 16	:['NESM3'],
		 17	:['NORESM2-MM']}
		 	
# Import cmip models and obs database 
cdd_naz_obs_x, cdd_naz_obs_y = import_obs('cdd', 'NAZ')
cdd_saz_obs_x, cdd_saz_obs_y = import_obs('cdd', 'SAZ')
cdd_neb_obs_x, cdd_neb_obs_y = import_obs('cdd', 'NEB')
cdd_sam_obs_x, cdd_sam_obs_y = import_obs('cdd', 'SAM')
cdd_lpb_obs_x, cdd_lpb_obs_y = import_obs('cdd', 'LPB')

r95p_naz_obs_x, r95p_naz_obs_y = import_obs('r95p', 'NAZ')
r95p_saz_obs_x, r95p_saz_obs_y = import_obs('r95p', 'SAZ')
r95p_neb_obs_x, r95p_neb_obs_y = import_obs('r95p', 'NEB')
r95p_sam_obs_x, r95p_sam_obs_y = import_obs('r95p', 'SAM')
r95p_lpb_obs_x, r95p_lpb_obs_y = import_obs('r95p', 'LPB')

rx5day_naz_obs_x, rx5day_naz_obs_y = import_obs('rx5day', 'NAZ')
rx5day_saz_obs_x, rx5day_saz_obs_y = import_obs('rx5day', 'SAZ')
rx5day_neb_obs_x, rx5day_neb_obs_y = import_obs('rx5day', 'NEB')
rx5day_sam_obs_x, rx5day_sam_obs_y = import_obs('rx5day', 'SAM')
rx5day_lpb_obs_x, rx5day_lpb_obs_y = import_obs('rx5day', 'LPB')

tx90p_naz_obs_x, tx90p_naz_obs_y = import_obs('tx90p', 'NAZ')
tx90p_saz_obs_x, tx90p_saz_obs_y = import_obs('tx90p', 'SAZ')
tx90p_neb_obs_x, tx90p_neb_obs_y = import_obs('tx90p', 'NEB')
tx90p_sam_obs_x, tx90p_sam_obs_y = import_obs('tx90p', 'SAM')
tx90p_lpb_obs_x, tx90p_lpb_obs_y = import_obs('tx90p', 'LPB')

wsdi_naz_obs_x, wsdi_naz_obs_y = import_obs('wsdi', 'NAZ')
wsdi_saz_obs_x, wsdi_saz_obs_y = import_obs('wsdi', 'SAZ')
wsdi_neb_obs_x, wsdi_neb_obs_y = import_obs('wsdi', 'NEB')
wsdi_sam_obs_x, wsdi_sam_obs_y = import_obs('wsdi', 'SAM')
wsdi_lpb_obs_x, wsdi_lpb_obs_y = import_obs('wsdi', 'LPB')

ree_cdd_naz, ree_cdd_saz, ree_cdd_neb, ree_cdd_sam, ree_cdd_lpb = [], [], [], [], []
pcc_cdd_naz, pcc_cdd_saz, pcc_cdd_neb, pcc_cdd_sam, pcc_cdd_lpb = [], [], [], [], []
ivs_cdd_naz, ivs_cdd_saz, ivs_cdd_neb, ivs_cdd_sam, ivs_cdd_lpb = [], [], [], [], []

ree_r95p_naz, ree_r95p_saz, ree_r95p_neb, ree_r95p_sam, ree_r95p_lpb = [], [], [], [], []
pcc_r95p_naz, pcc_r95p_saz, pcc_r95p_neb, pcc_r95p_sam, pcc_r95p_lpb = [], [], [], [], []
ivs_r95p_naz, ivs_r95p_saz, ivs_r95p_neb, ivs_r95p_sam, ivs_r95p_lpb = [], [], [], [], []

ree_rx5day_naz, ree_rx5day_saz, ree_rx5day_neb, ree_rx5day_sam, ree_rx5day_lpb = [], [], [], [], []
pcc_rx5day_naz, pcc_rx5day_saz, pcc_rx5day_neb, pcc_rx5day_sam, pcc_rx5day_lpb = [], [], [], [], []
ivs_rx5day_naz, ivs_rx5day_saz, ivs_rx5day_neb, ivs_rx5day_sam, ivs_rx5day_lpb = [], [], [], [], []

ree_tx90p_naz, ree_tx90p_saz, ree_tx90p_neb, ree_tx90p_sam, ree_tx90p_lpb = [], [], [], [], []
pcc_tx90p_naz, pcc_tx90p_saz, pcc_tx90p_neb, pcc_tx90p_sam, pcc_tx90p_lpb = [], [], [], [], []
ivs_tx90p_naz, ivs_tx90p_saz, ivs_tx90p_neb, ivs_tx90p_sam, ivs_tx90p_lpb = [], [], [], [], []

ree_wsdi_naz, ree_wsdi_saz, ree_wsdi_neb, ree_wsdi_sam, ree_wsdi_lpb = [], [], [], [], []
pcc_wsdi_naz, pcc_wsdi_saz, pcc_wsdi_neb, pcc_wsdi_sam, pcc_wsdi_lpb = [], [], [], [], []
ivs_wsdi_naz, ivs_wsdi_saz, ivs_wsdi_neb, ivs_wsdi_sam, ivs_wsdi_lpb = [], [], [], [], []

legend = []

for i in range(1, 18):
	
	print(cmip6[i])
	
	cdd_naz_cmip_x, cdd_naz_cmip_y = import_cmip('cddETCCDI', cmip6[i][0], 'NAZ')
	cdd_saz_cmip_x, cdd_saz_cmip_y = import_cmip('cddETCCDI', cmip6[i][0], 'SAZ')
	cdd_neb_cmip_x, cdd_neb_cmip_y = import_cmip('cddETCCDI', cmip6[i][0], 'NEB')
	cdd_sam_cmip_x, cdd_sam_cmip_y = import_cmip('cddETCCDI', cmip6[i][0], 'SAM')
	cdd_lpb_cmip_x, cdd_lpb_cmip_y = import_cmip('cddETCCDI', cmip6[i][0], 'LPB')

	r95p_naz_cmip_x, r95p_naz_cmip_y = import_cmip('r95pETCCDI', cmip6[i][0], 'NAZ')
	r95p_saz_cmip_x, r95p_saz_cmip_y = import_cmip('r95pETCCDI', cmip6[i][0], 'SAZ')
	r95p_neb_cmip_x, r95p_neb_cmip_y = import_cmip('r95pETCCDI', cmip6[i][0], 'NEB')
	r95p_sam_cmip_x, r95p_sam_cmip_y = import_cmip('r95pETCCDI', cmip6[i][0], 'SAM')
	r95p_lpb_cmip_x, r95p_lpb_cmip_y = import_cmip('r95pETCCDI', cmip6[i][0], 'LPB')

	rx5day_naz_cmip_x, rx5day_naz_cmip_y = import_cmip('rx5dayETCCDI', cmip6[i][0], 'NAZ')
	rx5day_saz_cmip_x, rx5day_saz_cmip_y = import_cmip('rx5dayETCCDI', cmip6[i][0], 'SAZ')
	rx5day_neb_cmip_x, rx5day_neb_cmip_y = import_cmip('rx5dayETCCDI', cmip6[i][0], 'NEB')
	rx5day_sam_cmip_x, rx5day_sam_cmip_y = import_cmip('rx5dayETCCDI', cmip6[i][0], 'SAM')
	rx5day_lpb_cmip_x, rx5day_lpb_cmip_y = import_cmip('rx5dayETCCDI', cmip6[i][0], 'LPB')

	tx90p_naz_cmip_x, tx90p_naz_cmip_y = import_cmip('tx90pETCCDI', cmip6[i][0], 'NAZ')
	tx90p_saz_cmip_x, tx90p_saz_cmip_y = import_cmip('tx90pETCCDI', cmip6[i][0], 'SAZ')
	tx90p_neb_cmip_x, tx90p_neb_cmip_y = import_cmip('tx90pETCCDI', cmip6[i][0], 'NEB')
	tx90p_sam_cmip_x, tx90p_sam_cmip_y = import_cmip('tx90pETCCDI', cmip6[i][0], 'SAM')
	tx90p_lpb_cmip_x, tx90p_lpb_cmip_y = import_cmip('tx90pETCCDI', cmip6[i][0], 'LPB')

	wsdi_naz_cmip_x, wsdi_naz_cmip_y = import_cmip('wsdiETCCDI', cmip6[i][0], 'NAZ')
	wsdi_saz_cmip_x, wsdi_saz_cmip_y = import_cmip('wsdiETCCDI', cmip6[i][0], 'SAZ')
	wsdi_neb_cmip_x, wsdi_neb_cmip_y = import_cmip('wsdiETCCDI', cmip6[i][0], 'NEB')
	wsdi_sam_cmip_x, wsdi_sam_cmip_y = import_cmip('wsdiETCCDI', cmip6[i][0], 'SAM')
	wsdi_lpb_cmip_x, wsdi_lpb_cmip_y = import_cmip('wsdiETCCDI', cmip6[i][0], 'LPB')
	
	# NAZ
	ree_cdd_naz.append(compute_ree(cdd_naz_cmip_x, cdd_naz_obs_x))
	pcc_cdd_naz.append(compute_pcc(cdd_naz_obs_x, cdd_naz_cmip_x))
	ivs_cdd_naz.append(compute_ivs(cdd_naz_obs_y, cdd_naz_cmip_y))
	
	ree_r95p_naz.append(compute_ree(r95p_naz_cmip_x, r95p_naz_obs_x))
	pcc_r95p_naz.append(compute_pcc(r95p_naz_obs_x, r95p_naz_cmip_x))
	ivs_r95p_naz.append(compute_ivs(r95p_naz_obs_y, r95p_naz_cmip_y))

	ree_rx5day_naz.append(compute_ree(rx5day_naz_cmip_x, rx5day_naz_obs_x))
	pcc_rx5day_naz.append(compute_pcc(rx5day_naz_obs_x, rx5day_naz_cmip_x))
	ivs_rx5day_naz.append(compute_ivs(rx5day_naz_obs_y, rx5day_naz_cmip_y))

	ree_tx90p_naz.append(compute_ree(tx90p_naz_cmip_x, tx90p_naz_obs_x))
	pcc_tx90p_naz.append(compute_pcc(tx90p_naz_obs_x, tx90p_naz_cmip_x))
	ivs_tx90p_naz.append(compute_ivs(tx90p_saz_obs_y, tx90p_naz_cmip_y))

	ree_wsdi_naz.append(compute_ree(wsdi_naz_cmip_x, wsdi_naz_obs_x))
	pcc_wsdi_naz.append(compute_pcc(wsdi_naz_obs_x, wsdi_naz_cmip_x))
	ivs_wsdi_naz.append(compute_ivs(wsdi_saz_obs_y, wsdi_naz_cmip_y))

	# SAZ
	ree_cdd_saz.append(compute_ree(cdd_saz_cmip_x, cdd_saz_obs_x))
	pcc_cdd_saz.append(compute_pcc(cdd_saz_obs_x, cdd_saz_cmip_x))
	ivs_cdd_saz.append(compute_ivs(cdd_saz_obs_y, cdd_saz_cmip_y))
	
	ree_r95p_saz.append(compute_ree(r95p_saz_cmip_x, r95p_saz_obs_x))
	pcc_r95p_saz.append(compute_pcc(r95p_saz_obs_x, r95p_saz_cmip_x))
	ivs_r95p_saz.append(compute_ivs(r95p_saz_obs_y, r95p_saz_cmip_y))

	ree_rx5day_saz.append(compute_ree(rx5day_saz_cmip_x, rx5day_saz_obs_x))
	pcc_rx5day_saz.append(compute_pcc(rx5day_saz_obs_x, rx5day_saz_cmip_x))
	ivs_rx5day_saz.append(compute_ivs(rx5day_saz_obs_y, rx5day_saz_cmip_y))

	ree_tx90p_saz.append(compute_ree(tx90p_saz_cmip_x, tx90p_saz_obs_x))
	pcc_tx90p_saz.append(compute_pcc(tx90p_saz_obs_x, tx90p_saz_cmip_x))
	ivs_tx90p_saz.append(compute_ivs(tx90p_saz_obs_y, tx90p_saz_cmip_y))

	ree_wsdi_saz.append(compute_ree(wsdi_saz_cmip_x, wsdi_saz_obs_x))
	pcc_wsdi_saz.append(compute_pcc(wsdi_saz_obs_x, wsdi_saz_cmip_x))
	ivs_wsdi_saz.append(compute_ivs(wsdi_saz_obs_y, wsdi_saz_cmip_y))
	
	# NEB
	ree_cdd_neb.append(compute_ree(cdd_neb_cmip_x, cdd_neb_obs_x))
	pcc_cdd_neb.append(compute_pcc(cdd_neb_obs_x, cdd_neb_cmip_x))
	ivs_cdd_neb.append(compute_ivs(cdd_neb_obs_y, cdd_neb_cmip_y))
	
	ree_r95p_neb.append(compute_ree(r95p_neb_cmip_x, r95p_neb_obs_x))
	pcc_r95p_neb.append(compute_pcc(r95p_neb_obs_x, r95p_neb_cmip_x))
	ivs_r95p_neb.append(compute_ivs(r95p_neb_obs_y, r95p_neb_cmip_y))

	ree_rx5day_neb.append(compute_ree(rx5day_neb_cmip_x, rx5day_neb_obs_x))
	pcc_rx5day_neb.append(compute_pcc(rx5day_neb_obs_x, rx5day_neb_cmip_x))
	ivs_rx5day_neb.append(compute_ivs(rx5day_neb_obs_y, rx5day_neb_cmip_y))
	
	ree_tx90p_neb.append(compute_ree(tx90p_neb_cmip_x, tx90p_neb_obs_x))
	pcc_tx90p_neb.append(compute_pcc(tx90p_neb_obs_x, tx90p_neb_cmip_x))
	ivs_tx90p_neb.append(compute_ivs(tx90p_neb_obs_y, tx90p_neb_cmip_y))

	ree_wsdi_neb.append(compute_ree(wsdi_neb_cmip_x, wsdi_neb_obs_x))
	pcc_wsdi_neb.append(compute_pcc(wsdi_neb_obs_x, wsdi_neb_cmip_x))
	ivs_wsdi_neb.append(compute_ivs(wsdi_neb_obs_y, wsdi_neb_cmip_y))

	# SAM
	ree_cdd_sam.append(compute_ree(cdd_sam_cmip_x, cdd_sam_obs_x))
	pcc_cdd_sam.append(compute_pcc(cdd_sam_obs_x, cdd_sam_cmip_x))
	ivs_cdd_sam.append(compute_ivs(cdd_sam_obs_y, cdd_sam_cmip_y))
	
	ree_r95p_sam.append(compute_ree(r95p_sam_cmip_x, r95p_sam_obs_x))
	pcc_r95p_sam.append(compute_pcc(r95p_sam_obs_x, r95p_sam_cmip_x))
	ivs_r95p_sam.append(compute_ivs(r95p_sam_obs_y, r95p_sam_cmip_y))

	ree_rx5day_sam.append(compute_ree(rx5day_sam_cmip_x, rx5day_sam_obs_x))
	pcc_rx5day_sam.append(compute_pcc(rx5day_sam_obs_x, rx5day_sam_cmip_x))
	ivs_rx5day_sam.append(compute_ivs(rx5day_sam_obs_y, rx5day_sam_cmip_y))

	ree_tx90p_sam.append(compute_ree(tx90p_sam_cmip_x, tx90p_sam_obs_x))
	pcc_tx90p_sam.append(compute_pcc(tx90p_sam_obs_x, tx90p_sam_cmip_x))
	ivs_tx90p_sam.append(compute_ivs(tx90p_sam_obs_y, tx90p_sam_cmip_y))
	
	ree_wsdi_sam.append(compute_ree(wsdi_sam_cmip_x, wsdi_sam_obs_x))
	pcc_wsdi_sam.append(compute_pcc(wsdi_sam_obs_x, wsdi_sam_cmip_x))
	ivs_wsdi_sam.append(compute_ivs(wsdi_sam_obs_y, wsdi_sam_cmip_y))

	# LPB
	ree_cdd_lpb.append(compute_ree(cdd_lpb_cmip_x, cdd_lpb_obs_x))
	pcc_cdd_lpb.append(compute_pcc(cdd_lpb_obs_x, cdd_lpb_cmip_x))
	ivs_cdd_lpb.append(compute_ivs(cdd_lpb_obs_y, cdd_lpb_cmip_y))
	
	ree_r95p_lpb.append(compute_ree(r95p_lpb_cmip_x, r95p_lpb_obs_x))
	pcc_r95p_lpb.append(compute_pcc(r95p_lpb_obs_x, r95p_lpb_cmip_x))
	ivs_r95p_lpb.append(compute_ivs(r95p_lpb_obs_y, r95p_lpb_cmip_y))

	ree_rx5day_lpb.append(compute_ree(rx5day_lpb_cmip_x, rx5day_lpb_obs_x))
	pcc_rx5day_lpb.append(compute_pcc(rx5day_lpb_obs_x, rx5day_lpb_cmip_x))
	ivs_rx5day_lpb.append(compute_ivs(rx5day_lpb_obs_y, rx5day_lpb_cmip_y))

	ree_tx90p_lpb.append(compute_ree(tx90p_lpb_cmip_x, tx90p_lpb_obs_x))
	pcc_tx90p_lpb.append(compute_pcc(tx90p_lpb_obs_x, tx90p_lpb_cmip_x))
	ivs_tx90p_lpb.append(compute_ivs(tx90p_lpb_obs_y, tx90p_lpb_cmip_y))
	
	ree_wsdi_lpb.append(compute_ree(wsdi_lpb_cmip_x, wsdi_lpb_obs_x))
	pcc_wsdi_lpb.append(compute_pcc(wsdi_lpb_obs_x, wsdi_lpb_cmip_x))
	ivs_wsdi_lpb.append(compute_ivs(wsdi_lpb_obs_y, wsdi_lpb_cmip_y))
		
	legend.append(cmip6[i][0])

# CDD
sort_ree_cdd_naz = sort_list_x(ree_cdd_naz)
sort_pcc_cdd_naz = sort_list_y(pcc_cdd_naz)
sort_ivs_cdd_naz = sort_list_z(ivs_cdd_naz)
	
sort_ree_cdd_saz = sort_list_x(ree_cdd_saz)
sort_pcc_cdd_saz = sort_list_y(pcc_cdd_saz)
sort_ivs_cdd_saz = sort_list_z(ivs_cdd_saz)

sort_ree_cdd_neb = sort_list_x(ree_cdd_neb)
sort_pcc_cdd_neb = sort_list_y(pcc_cdd_neb)
sort_ivs_cdd_neb = sort_list_z(ivs_cdd_neb)

sort_ree_cdd_sam = sort_list_x(ree_cdd_sam)
sort_pcc_cdd_sam = sort_list_y(pcc_cdd_sam)
sort_ivs_cdd_sam = sort_list_z(ivs_cdd_sam)

sort_ree_cdd_lpb = sort_list_x(ree_cdd_lpb)
sort_pcc_cdd_lpb = sort_list_y(pcc_cdd_lpb)
sort_ivs_cdd_lpb = sort_list_z(ivs_cdd_lpb)

# R95P
sort_ree_r95p_naz = sort_list_x(ree_r95p_naz)
sort_pcc_r95p_naz = sort_list_y(pcc_r95p_naz)
sort_ivs_r95p_naz = sort_list_z(ivs_r95p_naz)

sort_ree_r95p_saz = sort_list_x(ree_r95p_saz)
sort_pcc_r95p_saz = sort_list_y(pcc_r95p_saz)
sort_ivs_r95p_saz = sort_list_z(ivs_r95p_saz)

sort_ree_r95p_neb = sort_list_x(ree_r95p_neb)
sort_pcc_r95p_neb = sort_list_y(pcc_r95p_neb)
sort_ivs_r95p_neb = sort_list_z(ivs_r95p_neb)

sort_ree_r95p_sam = sort_list_x(ree_r95p_sam)
sort_pcc_r95p_sam = sort_list_y(pcc_r95p_sam)
sort_ivs_r95p_sam = sort_list_z(ivs_r95p_sam)

sort_ree_r95p_lpb = sort_list_x(ree_r95p_lpb)
sort_pcc_r95p_lpb = sort_list_y(pcc_r95p_lpb)
sort_ivs_r95p_lpb = sort_list_z(ivs_r95p_lpb)

# RX5DAY
sort_ree_rx5day_naz = sort_list_x(ree_rx5day_naz)
sort_pcc_rx5day_naz = sort_list_y(pcc_rx5day_naz)
sort_ivs_rx5day_naz = sort_list_z(ivs_rx5day_naz)

sort_ree_rx5day_saz = sort_list_x(ree_rx5day_saz)
sort_pcc_rx5day_saz = sort_list_y(pcc_rx5day_saz)
sort_ivs_rx5day_saz = sort_list_z(ivs_rx5day_saz)

sort_ree_rx5day_neb = sort_list_x(ree_rx5day_neb)
sort_pcc_rx5day_neb = sort_list_y(pcc_rx5day_neb)
sort_ivs_rx5day_neb = sort_list_z(ivs_rx5day_neb)

sort_ree_rx5day_sam = sort_list_x(ree_rx5day_sam)
sort_pcc_rx5day_sam = sort_list_y(pcc_rx5day_sam)
sort_ivs_rx5day_sam = sort_list_z(ivs_rx5day_sam)

sort_ree_rx5day_lpb = sort_list_x(ree_rx5day_lpb)
sort_pcc_rx5day_lpb = sort_list_y(pcc_rx5day_lpb)
sort_ivs_rx5day_lpb = sort_list_z(ivs_rx5day_lpb)

# TX90P
sort_ree_tx90p_naz = sort_list_x(ree_tx90p_naz)
sort_pcc_tx90p_naz = sort_list_y(pcc_tx90p_naz)
sort_ivs_tx90p_naz = sort_list_z(ivs_tx90p_naz)

sort_ree_tx90p_saz = sort_list_x(ree_tx90p_saz)
sort_pcc_tx90p_saz = sort_list_y(pcc_tx90p_saz)
sort_ivs_tx90p_saz = sort_list_z(ivs_tx90p_saz)

sort_ree_tx90p_neb = sort_list_x(ree_tx90p_neb)
sort_pcc_tx90p_neb = sort_list_y(pcc_tx90p_neb)
sort_ivs_tx90p_neb = sort_list_z(ivs_tx90p_neb)

sort_ree_tx90p_sam = sort_list_x(ree_tx90p_sam)
sort_pcc_tx90p_sam = sort_list_y(pcc_tx90p_sam)
sort_ivs_tx90p_sam = sort_list_z(ivs_tx90p_sam)

sort_ree_tx90p_lpb = sort_list_x(ree_tx90p_lpb)
sort_pcc_tx90p_lpb = sort_list_y(pcc_tx90p_lpb)
sort_ivs_tx90p_lpb = sort_list_z(ivs_tx90p_lpb)

# WSDI
sort_ree_wsdi_naz = sort_list_x(ree_wsdi_naz)
sort_pcc_wsdi_naz = sort_list_y(pcc_wsdi_naz)
sort_ivs_wsdi_naz = sort_list_z(ivs_wsdi_naz)

sort_ree_wsdi_saz = sort_list_x(ree_wsdi_saz)
sort_pcc_wsdi_saz = sort_list_y(pcc_wsdi_saz)
sort_ivs_wsdi_saz = sort_list_z(ivs_wsdi_saz)

sort_ree_wsdi_neb = sort_list_x(ree_wsdi_neb)
sort_pcc_wsdi_neb = sort_list_y(pcc_wsdi_neb)
sort_ivs_wsdi_neb = sort_list_z(ivs_wsdi_neb)

sort_ree_wsdi_sam = sort_list_x(ree_wsdi_sam)
sort_pcc_wsdi_sam = sort_list_y(pcc_wsdi_sam)
sort_ivs_wsdi_sam = sort_list_z(ivs_wsdi_sam)

sort_ree_wsdi_lpb = sort_list_x(ree_wsdi_lpb)
sort_pcc_wsdi_lpb = sort_list_y(pcc_wsdi_lpb)
sort_ivs_wsdi_lpb = sort_list_z(ivs_wsdi_lpb)

rank_cdd_naz = np.array([sort_ivs_cdd_naz, sort_pcc_cdd_naz, sort_ree_cdd_naz])
rank_cdd_saz = np.array([sort_ivs_cdd_saz, sort_pcc_cdd_saz, sort_ree_cdd_saz])
rank_cdd_neb = np.array([sort_ivs_cdd_neb, sort_pcc_cdd_neb, sort_ree_cdd_neb])
rank_cdd_sam = np.array([sort_ivs_cdd_sam, sort_pcc_cdd_sam, sort_ree_cdd_sam])
rank_cdd_lpb = np.array([sort_ivs_cdd_lpb, sort_pcc_cdd_lpb, sort_ree_cdd_lpb])

rank_r95p_naz = np.array([sort_ivs_r95p_naz, sort_pcc_r95p_naz, sort_ree_r95p_naz])
rank_r95p_saz = np.array([sort_ivs_r95p_saz, sort_pcc_r95p_saz, sort_ree_r95p_saz])
rank_r95p_neb = np.array([sort_ivs_r95p_neb, sort_pcc_r95p_neb, sort_ree_r95p_neb])
rank_r95p_sam = np.array([sort_ivs_r95p_sam, sort_pcc_r95p_sam, sort_ree_r95p_sam])
rank_r95p_lpb = np.array([sort_ivs_r95p_lpb, sort_pcc_r95p_lpb, sort_ree_r95p_lpb])

rank_rx5day_naz = np.array([sort_ivs_rx5day_naz, sort_pcc_rx5day_naz, sort_ree_rx5day_naz])
rank_rx5day_saz = np.array([sort_ivs_rx5day_saz, sort_pcc_rx5day_saz, sort_ree_rx5day_saz])
rank_rx5day_neb = np.array([sort_ivs_rx5day_neb, sort_pcc_rx5day_neb, sort_ree_rx5day_neb])
rank_rx5day_sam = np.array([sort_ivs_rx5day_sam, sort_pcc_rx5day_sam, sort_ree_rx5day_sam])
rank_rx5day_lpb = np.array([sort_ivs_rx5day_lpb, sort_pcc_rx5day_lpb, sort_ree_rx5day_lpb])

rank_tx90p_naz = np.array([sort_ivs_tx90p_naz, sort_pcc_tx90p_naz, sort_ree_tx90p_naz])
rank_tx90p_saz = np.array([sort_ivs_tx90p_saz, sort_pcc_tx90p_saz, sort_ree_tx90p_saz])
rank_tx90p_neb = np.array([sort_ivs_tx90p_neb, sort_pcc_tx90p_neb, sort_ree_tx90p_neb])
rank_tx90p_sam = np.array([sort_ivs_tx90p_sam, sort_pcc_tx90p_sam, sort_ree_tx90p_sam])
rank_tx90p_lpb = np.array([sort_ivs_tx90p_lpb, sort_pcc_tx90p_lpb, sort_ree_tx90p_lpb])

rank_wsdi_naz = np.array([sort_ivs_wsdi_naz, sort_pcc_wsdi_naz, sort_ree_wsdi_naz])
rank_wsdi_saz = np.array([sort_ivs_wsdi_saz, sort_pcc_wsdi_saz, sort_ree_wsdi_saz])
rank_wsdi_neb = np.array([sort_ivs_wsdi_neb, sort_pcc_wsdi_neb, sort_ree_wsdi_neb])
rank_wsdi_sam = np.array([sort_ivs_wsdi_sam, sort_pcc_wsdi_sam, sort_ree_wsdi_sam])
rank_wsdi_lpb = np.array([sort_ivs_wsdi_lpb, sort_pcc_wsdi_lpb, sort_ree_wsdi_lpb])

# Plot cmip models and obs database 
fig = plt.figure(figsize=(10, 8))

xlabels = legend
ylabels = ['IVS', 'PCC', 'RE']

norm = colors.BoundaryNorm(boundaries=np.arange(1, 18, 1), ncolors=256)
color = cm.rainbow

ax = fig.add_subplot(5, 5, 1)  
pcm = ax.pcolormesh(rank_cdd_naz, edgecolors='white', linewidths=0., norm=norm, cmap=color)
ax.set_title(u'(a) CDD', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(rank_cdd_naz.shape[1]) + 0.5)
ax.set_yticks(np.arange(rank_cdd_naz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8)
ax.set_yticklabels(ylabels, fontsize=8)
ax.set_ylabel('NAZ', fontweight='bold', fontsize=8, rotation=-90, labelpad=455)
ax.yaxis.set_label_position("right")
plt.setp(ax.get_xticklabels(), visible=False)

ax = fig.add_subplot(5, 5, 2)  
pcm = ax.pcolormesh(rank_r95p_naz, edgecolors='white', linewidths=0., norm=norm, cmap=color)
ax.set_title(u'(b) R95p', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(rank_r95p_naz.shape[1]) + 0.5)
ax.set_yticks(np.arange(rank_r95p_naz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8)
ax.set_yticklabels(ylabels, fontsize=8)
ax.get_yaxis().set_visible(False)
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)

ax = fig.add_subplot(5, 5, 3)  
pcm = ax.pcolormesh(rank_rx5day_naz, edgecolors='white', linewidths=0., norm=norm, cmap=color)
ax.set_title(u'(c) Rx5day', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(rank_rx5day_naz.shape[1]) + 0.5)
ax.set_yticks(np.arange(rank_rx5day_naz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8)
ax.set_yticklabels(ylabels, fontsize=8)
ax.get_yaxis().set_visible(False)
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)

ax = fig.add_subplot(5, 5, 4)  
pcm = ax.pcolormesh(rank_tx90p_naz, edgecolors='white', linewidths=0., norm=norm, cmap=color)
ax.set_title(u'(d) Tx90p', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(rank_tx90p_naz.shape[1]) + 0.5)
ax.set_yticks(np.arange(rank_tx90p_naz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8)
ax.set_yticklabels(ylabels, fontsize=8)
ax.get_yaxis().set_visible(False)
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)

ax = fig.add_subplot(5, 5, 5)  
pcm = ax.pcolormesh(rank_wsdi_naz, edgecolors='white', linewidths=0., norm=norm, cmap=color)
ax.set_title(u'(e) WSDI', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(rank_wsdi_naz.shape[1]) + 0.5)
ax.set_yticks(np.arange(rank_wsdi_naz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8)
ax.set_yticklabels(ylabels, fontsize=8)
ax.get_yaxis().set_visible(False)
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)
cbar = plt.colorbar(pcm, cax=fig.add_axes([0.92, 0.28, 0.02, 0.43]), pad=0.01)
cbar.ax.tick_params(labelsize=8)

ax = fig.add_subplot(5, 5, 6)  
pcm = ax.pcolormesh(rank_cdd_saz, edgecolors='white', linewidths=0., norm=norm, cmap=color)
ax.set_title(u'(f)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(rank_cdd_saz.shape[1]) + 0.5)
ax.set_yticks(np.arange(rank_cdd_saz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8)
ax.set_yticklabels(ylabels, fontsize=8)
ax.set_ylabel('SAZ', fontweight='bold', fontsize=8, rotation=-90, labelpad=455)
ax.yaxis.set_label_position("right")
plt.setp(ax.get_xticklabels(), visible=False)

ax = fig.add_subplot(5, 5, 7)  
pcm = ax.pcolormesh(rank_r95p_saz, edgecolors='white', linewidths=0., norm=norm, cmap=color)
ax.set_title(u'(g)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(rank_r95p_saz.shape[1]) + 0.5)
ax.set_yticks(np.arange(rank_r95p_saz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8)
ax.set_yticklabels(ylabels, fontsize=8)
ax.get_yaxis().set_visible(False)
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)

ax = fig.add_subplot(5, 5, 8)  
pcm = ax.pcolormesh(rank_rx5day_saz, edgecolors='white', linewidths=0., norm=norm, cmap=color)
ax.set_title(u'(h)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(rank_rx5day_saz.shape[1]) + 0.5)
ax.set_yticks(np.arange(rank_rx5day_saz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8)
ax.set_yticklabels(ylabels, fontsize=8)
ax.get_yaxis().set_visible(False)
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)

ax = fig.add_subplot(5, 5, 9)  
pcm = ax.pcolormesh(rank_tx90p_saz, edgecolors='white', linewidths=0., norm=norm, cmap=color)
ax.set_title(u'(i)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(rank_tx90p_saz.shape[1]) + 0.5)
ax.set_yticks(np.arange(rank_tx90p_saz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8)
ax.set_yticklabels(ylabels, fontsize=8)
ax.get_yaxis().set_visible(False)
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)

ax = fig.add_subplot(5, 5, 10)  
pcm = ax.pcolormesh(rank_wsdi_saz, edgecolors='white', linewidths=0., norm=norm, cmap=color)
ax.set_title(u'(j)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(rank_wsdi_saz.shape[1]) + 0.5)
ax.set_yticks(np.arange(rank_wsdi_saz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8)
ax.set_yticklabels(ylabels, fontsize=8)
ax.get_yaxis().set_visible(False)
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)

ax = fig.add_subplot(5, 5, 11)  
pcm = ax.pcolormesh(rank_cdd_neb, edgecolors='white', linewidths=0., norm=norm, cmap=color)
ax.set_title(u'(k)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(rank_cdd_neb.shape[1]) + 0.5)
ax.set_yticks(np.arange(rank_cdd_neb.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8)
ax.set_yticklabels(ylabels, fontsize=8)
ax.set_ylabel('NEB', fontweight='bold', fontsize=8, rotation=-90, labelpad=455)
ax.yaxis.set_label_position("right")
plt.setp(ax.get_xticklabels(), visible=False)

ax = fig.add_subplot(5, 5, 12)  
pcm = ax.pcolormesh(rank_r95p_neb, edgecolors='white', linewidths=0., norm=norm, cmap=color)
ax.set_title(u'(l)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(rank_r95p_neb.shape[1]) + 0.5)
ax.set_yticks(np.arange(rank_r95p_neb.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8)
ax.set_yticklabels(ylabels, fontsize=8)
ax.get_yaxis().set_visible(False)
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)

ax = fig.add_subplot(5, 5, 13)  
pcm = ax.pcolormesh(rank_rx5day_neb, edgecolors='white', linewidths=0., norm=norm, cmap=color)
ax.set_title(u'(m)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(rank_rx5day_neb.shape[1]) + 0.5)
ax.set_yticks(np.arange(rank_rx5day_neb.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8)
ax.set_yticklabels(ylabels, fontsize=8)
ax.get_yaxis().set_visible(False)
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)

ax = fig.add_subplot(5, 5, 14)  
pcm = ax.pcolormesh(rank_tx90p_neb, edgecolors='white', linewidths=0., norm=norm, cmap=color)
ax.set_title(u'(n)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(rank_tx90p_neb.shape[1]) + 0.5)
ax.set_yticks(np.arange(rank_tx90p_neb.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8)
ax.set_yticklabels(ylabels, fontsize=8)
ax.get_yaxis().set_visible(False)
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)

ax = fig.add_subplot(5, 5, 15)  
pcm = ax.pcolormesh(rank_wsdi_neb, edgecolors='white', linewidths=0., norm=norm, cmap=color)
ax.set_title(u'(o)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(rank_wsdi_neb.shape[1]) + 0.5)
ax.set_yticks(np.arange(rank_wsdi_neb.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8)
ax.set_yticklabels(ylabels, fontsize=8)
ax.get_yaxis().set_visible(False)
plt.setp(ax.get_xticklabels(), visible=False)

ax = fig.add_subplot(5, 5, 16)  
pcm = ax.pcolormesh(rank_cdd_sam, edgecolors='white', linewidths=0., norm=norm, cmap=color)
ax.set_title(u'(p)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(rank_cdd_sam.shape[1]) + 0.5)
ax.set_yticks(np.arange(rank_cdd_sam.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8)
ax.set_yticklabels(ylabels, fontsize=8)
ax.set_ylabel('SAM', fontweight='bold', fontsize=8, rotation=-90, labelpad=455)
ax.yaxis.set_label_position("right")
plt.setp(ax.get_xticklabels(), visible=False)

ax = fig.add_subplot(5, 5, 17)  
pcm = ax.pcolormesh(rank_r95p_sam, edgecolors='white', linewidths=0., norm=norm, cmap=color)
ax.set_title(u'(q)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(rank_r95p_sam.shape[1]) + 0.5)
ax.set_yticks(np.arange(rank_r95p_sam.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8)
ax.set_yticklabels(ylabels, fontsize=8)
ax.get_yaxis().set_visible(False)
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)

ax = fig.add_subplot(5, 5, 18)  
pcm = ax.pcolormesh(rank_rx5day_sam, edgecolors='white', linewidths=0., norm=norm, cmap=color)
ax.set_title(u'(r)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(rank_rx5day_saz.shape[1]) + 0.5)
ax.set_yticks(np.arange(rank_rx5day_saz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8)
ax.set_yticklabels(ylabels, fontsize=8)
ax.get_yaxis().set_visible(False)
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)

ax = fig.add_subplot(5, 5, 19)  
pcm = ax.pcolormesh(rank_tx90p_sam, edgecolors='white', linewidths=0., norm=norm, cmap=color)
ax.set_title(u'(s)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(rank_tx90p_saz.shape[1]) + 0.5)
ax.set_yticks(np.arange(rank_tx90p_saz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8)
ax.set_yticklabels(ylabels, fontsize=8)
ax.get_yaxis().set_visible(False)
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)

ax = fig.add_subplot(5, 5, 20)  
pcm = ax.pcolormesh(rank_wsdi_sam, edgecolors='white', linewidths=0., norm=norm, cmap=color)
ax.set_title(u'(t)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(rank_wsdi_saz.shape[1]) + 0.5)
ax.set_yticks(np.arange(rank_wsdi_saz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8)
ax.set_yticklabels(ylabels, fontsize=8)
ax.get_yaxis().set_visible(False)
ax.yaxis.set_label_position("right")
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)

ax = fig.add_subplot(5, 5, 21)  
pcm = ax.pcolormesh(rank_cdd_lpb, edgecolors='white', linewidths=0., norm=norm, cmap=color)
ax.set_title(u'(u)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(rank_cdd_saz.shape[1]) + 0.5)
ax.set_yticks(np.arange(rank_cdd_saz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=7, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
ax.set_ylabel('LPB', fontweight='bold', fontsize=8, rotation=-90, labelpad=455)
ax.yaxis.set_label_position("right")

ax = fig.add_subplot(5, 5, 22)  
pcm = ax.pcolormesh(rank_r95p_lpb, edgecolors='white', linewidths=0., norm=norm, cmap=color)
ax.set_title(u'(v)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(rank_r95p_saz.shape[1]) + 0.5)
ax.set_yticks(np.arange(rank_r95p_saz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=7, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
ax.get_yaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
plt.setp(ax.get_yticklabels(), visible=False)

ax = fig.add_subplot(5, 5, 23)  
pcm = ax.pcolormesh(rank_rx5day_lpb, edgecolors='white', linewidths=0., norm=norm, cmap=color)
ax.set_title(u'(w)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(rank_rx5day_saz.shape[1]) + 0.5)
ax.set_yticks(np.arange(rank_rx5day_saz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=7, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
ax.get_yaxis().set_visible(False)
plt.setp(ax.get_yticklabels(), visible=False)

ax = fig.add_subplot(5, 5, 24)  
pcm = ax.pcolormesh(rank_tx90p_lpb, edgecolors='white', linewidths=0., norm=norm, cmap=color)
ax.set_title(u'(x)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(rank_tx90p_saz.shape[1]) + 0.5)
ax.set_yticks(np.arange(rank_tx90p_saz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=7, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
ax.get_yaxis().set_visible(False)
plt.setp(ax.get_yticklabels(), visible=False)

ax = fig.add_subplot(5, 5, 25)  
pcm = ax.pcolormesh(rank_wsdi_lpb, edgecolors='white', linewidths=0., norm=norm, cmap=color)
ax.set_title(u'(y)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(rank_wsdi_saz.shape[1]) + 0.5)
ax.set_yticks(np.arange(rank_wsdi_saz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=7, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
ax.get_yaxis().set_visible(False)
ax.set_ylabel('LPB', fontweight='bold', fontsize=8, rotation=-90, labelpad=10)
ax.yaxis.set_label_position("right")
plt.setp(ax.get_yticklabels(), visible=False)

plt.subplots_adjust(wspace=0.)
plt.subplots_adjust(hspace=0.25)

# Path out to save figure
path_out = '/home/nice/Documentos/paper_mari/figs'
name_out = 'pyplt_portrait_rank_etccdi_indices_cmip6.png'
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()











	


