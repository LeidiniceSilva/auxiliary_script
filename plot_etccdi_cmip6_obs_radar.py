# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jun 01, 2023"
__description__ = "This script plot cri rank of cmip6 models"

import os
import netCDF4
import numpy as np
import numpy.ma as ma
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors


def import_obs(param, area):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/paper_mari/database/obs'
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
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/paper_mari/database/cmip6'
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


def compute_cri(rank1,rank2,rank3):
	
	p1 = (rank1+rank2+rank3)
	p2 = p1/51
	cri = 1 - p2
	
	return cri
	

	
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
	
	print(i, ':  ', cmip6[i])
	
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

cri_cdd_naz, cri_cdd_saz, cri_cdd_neb, cri_cdd_sam, cri_cdd_lpb = [], [], [], [], []
cri_r95p_naz, cri_r95p_saz, cri_r95p_neb,cri_r95p_sam, cri_r95p_lpb = [], [], [], [], []
cri_rx5day_naz, cri_rx5day_saz, cri_rx5day_neb, cri_rx5day_sam, cri_rx5day_lpb = [], [], [], [], []
cri_tx90p_naz, cri_tx90p_saz, cri_tx90p_neb, cri_tx90p_sam, cri_tx90p_lpb = [], [], [], [], []
cri_wsdi_naz, cri_wsdi_saz, cri_wsdi_neb, cri_wsdi_sam, cri_wsdi_lpb = [], [], [], [], []

for i in range(0, 17):

	cri_cdd_naz.append(compute_cri(sort_ree_cdd_naz[i],sort_pcc_cdd_naz[i],sort_ivs_cdd_naz[i]))
	cri_cdd_saz.append(compute_cri(sort_ree_cdd_saz[i],sort_pcc_cdd_saz[i],sort_ivs_cdd_saz[i]))
	cri_cdd_neb.append(compute_cri(sort_ree_cdd_neb[i],sort_pcc_cdd_neb[i],sort_ivs_cdd_neb[i]))
	cri_cdd_sam.append(compute_cri(sort_ree_cdd_sam[i],sort_pcc_cdd_sam[i],sort_ivs_cdd_sam[i]))
	cri_cdd_lpb.append(compute_cri(sort_ree_cdd_lpb[i],sort_pcc_cdd_lpb[i],sort_ivs_cdd_lpb[i]))

	cri_r95p_naz.append(compute_cri(sort_ree_r95p_naz[i],sort_pcc_r95p_naz[i],sort_ivs_r95p_naz[i]))
	cri_r95p_saz.append(compute_cri(sort_ree_r95p_saz[i],sort_pcc_r95p_saz[i],sort_ivs_r95p_saz[i]))
	cri_r95p_neb.append(compute_cri(sort_ree_r95p_neb[i],sort_pcc_r95p_neb[i],sort_ivs_r95p_neb[i]))
	cri_r95p_sam.append(compute_cri(sort_ree_r95p_sam[i],sort_pcc_r95p_sam[i],sort_ivs_r95p_sam[i]))
	cri_r95p_lpb.append(compute_cri(sort_ree_r95p_lpb[i],sort_pcc_r95p_lpb[i],sort_ivs_r95p_lpb[i]))

	cri_rx5day_naz.append(compute_cri(sort_ree_rx5day_naz[i],sort_pcc_rx5day_naz[i],sort_ivs_rx5day_naz[i]))
	cri_rx5day_saz.append(compute_cri(sort_ree_rx5day_saz[i],sort_pcc_rx5day_saz[i],sort_ivs_rx5day_saz[i]))
	cri_rx5day_neb.append(compute_cri(sort_ree_rx5day_neb[i],sort_pcc_rx5day_neb[i],sort_ivs_rx5day_neb[i]))
	cri_rx5day_sam.append(compute_cri(sort_ree_rx5day_sam[i],sort_pcc_rx5day_sam[i],sort_ivs_rx5day_sam[i]))
	cri_rx5day_lpb.append(compute_cri(sort_ree_rx5day_lpb[i],sort_pcc_rx5day_lpb[i],sort_ivs_rx5day_lpb[i]))

	cri_tx90p_naz.append(compute_cri(sort_ree_tx90p_naz[i],sort_pcc_tx90p_naz[i],sort_ivs_tx90p_naz[i]))
	cri_tx90p_saz.append(compute_cri(sort_ree_tx90p_saz[i],sort_pcc_tx90p_saz[i],sort_ivs_tx90p_saz[i]))
	cri_tx90p_neb.append(compute_cri(sort_ree_tx90p_neb[i],sort_pcc_tx90p_neb[i],sort_ivs_tx90p_neb[i]))
	cri_tx90p_sam.append(compute_cri(sort_ree_tx90p_sam[i],sort_pcc_tx90p_sam[i],sort_ivs_tx90p_sam[i]))
	cri_tx90p_lpb.append(compute_cri(sort_ree_tx90p_lpb[i],sort_pcc_tx90p_lpb[i],sort_ivs_tx90p_lpb[i]))

	cri_wsdi_naz.append(compute_cri(sort_ree_wsdi_naz[i],sort_pcc_wsdi_naz[i],sort_ivs_wsdi_naz[i]))
	cri_wsdi_saz.append(compute_cri(sort_ree_wsdi_saz[i],sort_pcc_wsdi_saz[i],sort_ivs_wsdi_saz[i]))
	cri_wsdi_neb.append(compute_cri(sort_ree_wsdi_neb[i],sort_pcc_wsdi_neb[i],sort_ivs_wsdi_neb[i]))
	cri_wsdi_sam.append(compute_cri(sort_ree_wsdi_sam[i],sort_pcc_wsdi_sam[i],sort_ivs_wsdi_sam[i]))
	cri_wsdi_lpb.append(compute_cri(sort_ree_wsdi_lpb[i],sort_pcc_wsdi_lpb[i],sort_ivs_wsdi_lpb[i]))

sort_cri_cdd_naz = sort_list_y(cri_cdd_naz)
sort_cri_cdd_saz = sort_list_y(cri_cdd_saz)
sort_cri_cdd_neb = sort_list_y(cri_cdd_neb)
sort_cri_cdd_sam = sort_list_y(cri_cdd_sam)
sort_cri_cdd_lpb = sort_list_y(cri_cdd_lpb)

sort_cri_r95p_naz = sort_list_y(cri_r95p_naz)
sort_cri_r95p_saz = sort_list_y(cri_r95p_saz)
sort_cri_r95p_neb = sort_list_y(cri_r95p_neb)
sort_cri_r95p_sam = sort_list_y(cri_r95p_sam)
sort_cri_r95p_lpb = sort_list_y(cri_r95p_lpb)

sort_cri_rx5day_naz = sort_list_y(cri_rx5day_naz)
sort_cri_rx5day_saz = sort_list_y(cri_rx5day_saz)
sort_cri_rx5day_neb = sort_list_y(cri_rx5day_neb)
sort_cri_rx5day_sam = sort_list_y(cri_rx5day_sam)
sort_cri_rx5day_lpb = sort_list_y(cri_rx5day_lpb)

sort_cri_tx90p_naz = sort_list_y(cri_tx90p_naz)
sort_cri_tx90p_saz = sort_list_y(cri_tx90p_saz)
sort_cri_tx90p_neb = sort_list_y(cri_tx90p_neb)
sort_cri_tx90p_sam = sort_list_y(cri_tx90p_sam)
sort_cri_tx90p_lpb = sort_list_y(cri_tx90p_lpb)

sort_cri_wsdi_naz = sort_list_y(cri_wsdi_naz)
sort_cri_wsdi_saz = sort_list_y(cri_wsdi_saz)
sort_cri_wsdi_neb = sort_list_y(cri_wsdi_neb)
sort_cri_wsdi_sam = sort_list_y(cri_wsdi_sam)
sort_cri_wsdi_lpb = sort_list_y(cri_wsdi_lpb)

print('CRI wsdi NAZ: ', sort_cri_wsdi_naz)
print('CRI wsdi SAZ: ', sort_cri_wsdi_saz)
print('CRI wsdi NEB: ', sort_cri_wsdi_neb)
print('CRI wsdi SAM: ', sort_cri_wsdi_sam)
print('CRI wsdi LPB: ', sort_cri_wsdi_lpb)
exit()

# Plot cmip models and obs database 
fig = plt.figure(figsize=(9, 7))

labels = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17']
labels = np.array(labels)
angles = np.arange(0, 2*np.pi, 2*np.pi/len(labels)) 
angles = np.concatenate((angles,[angles[0]]))

ax = fig.add_subplot(2, 3, 1, polar=True)
stats1 = np.concatenate((sort_cri_cdd_naz,[sort_cri_cdd_naz[0]]))
stats2 = np.concatenate((sort_cri_r95p_naz,[sort_cri_r95p_naz[0]]))
stats3 = np.concatenate((sort_cri_rx5day_naz,[sort_cri_rx5day_naz[0]]))
stats4 = np.concatenate((sort_cri_tx90p_naz,[sort_cri_tx90p_naz[0]]))
stats5 = np.concatenate((sort_cri_wsdi_naz,[sort_cri_wsdi_naz[0]]))
ax.plot(angles, stats1, color='green', linewidth=1, label='CDD')
ax.plot(angles, stats2, color='blue', linewidth=1, label='R95p')
ax.plot(angles, stats3, color='purple', linewidth=1, label='Rx5day')
ax.plot(angles, stats4, color='brown', linewidth=1, label='Tx90p')
ax.plot(angles, stats5, color='orange', linewidth=1, label='WSDI')
# ~ ax.plot(angles, stats1, color='green', linewidth=1)
# ~ ax.plot(angles, stats2, color='blue', linewidth=1)
# ~ ax.plot(angles, stats3, color='purple', linewidth=1)
# ~ ax.plot(angles, stats4, color='brown', linewidth=1)
# ~ ax.plot(angles, stats5, color='orange', linewidth=1)
# ~ ax.plot(angles, stats1, color='green', linewidth=0.1)
# ~ ax.fill(angles, stats1, color='green', alpha=0.25, label='CDD')
# ~ ax.plot(angles, stats2, color='blue', linewidth=0.1)
# ~ ax.fill(angles, stats2, color='blue', alpha=0.25, label='R95p')
# ~ ax.plot(angles, stats3, color='purple', linewidth=0.1)
# ~ ax.fill(angles, stats3, color='purple', alpha=0.25, label='Rx5day')
# ~ ax.plot(angles, stats4, color='brown', linewidth=0.1)
# ~ ax.fill(angles, stats4, color='brown', alpha=0.25, label='Tx90p')
# ~ ax.plot(angles, stats5, color='orange', linewidth=0.1)
# ~ ax.fill(angles, stats5, color='orange', alpha=0.25, label='WSDI')
ax.set_title(u'(a) NAZ', loc='left', fontweight='bold', fontsize=8)
ax.set_thetagrids(angles[:-1] * 180/np.pi, labels, fontsize=8)
ax.grid(True)

ax = fig.add_subplot(2, 3, 2, polar=True)
stats1 = np.concatenate((sort_cri_cdd_saz,[sort_cri_cdd_saz[0]]))
stats2 = np.concatenate((sort_cri_r95p_saz,[sort_cri_r95p_saz[0]]))
stats3 = np.concatenate((sort_cri_rx5day_saz,[sort_cri_rx5day_saz[0]]))
stats4 = np.concatenate((sort_cri_tx90p_saz,[sort_cri_tx90p_saz[0]]))
stats5 = np.concatenate((sort_cri_wsdi_saz,[sort_cri_wsdi_saz[0]]))
ax.plot(angles, stats1, color='green', linewidth=1, label='CDD')
ax.plot(angles, stats2, color='blue', linewidth=1, label='R95p')
ax.plot(angles, stats3, color='purple', linewidth=1, label='Rx5day')
ax.plot(angles, stats4, color='brown', linewidth=1, label='Tx90p')
ax.plot(angles, stats5, color='orange', linewidth=1, label='WSDI')
# ~ ax.plot(angles, stats1, color='green', linewidth=0.1)
# ~ ax.fill(angles, stats1, color='green', alpha=0.25, label='CDD')
# ~ ax.plot(angles, stats2, color='blue', linewidth=0.1)
# ~ ax.fill(angles, stats2, color='blue', alpha=0.25, label='R95p')
# ~ ax.plot(angles, stats3, color='purple', linewidth=0.1)
# ~ ax.fill(angles, stats3, color='purple', alpha=0.25, label='Rx5day')
# ~ ax.plot(angles, stats4, color='brown', linewidth=0.1)
# ~ ax.fill(angles, stats4, color='brown', alpha=0.25, label='Tx90p')
# ~ ax.plot(angles, stats5, color='orange', linewidth=0.1)
# ~ ax.fill(angles, stats5, color='orange', alpha=0.25, label='WSDI')
ax.set_title(u'(b) SAZ', loc='left', fontweight='bold', fontsize=8)
ax.set_thetagrids(angles[:-1] * 180/np.pi, labels, fontsize=8)
ax.grid(True)

ax = fig.add_subplot(2, 3, 3, polar=True)
stats1 = np.concatenate((sort_cri_cdd_neb,[sort_cri_cdd_neb[0]]))
stats2 = np.concatenate((sort_cri_r95p_neb,[sort_cri_r95p_neb[0]]))
stats3 = np.concatenate((sort_cri_rx5day_neb,[sort_cri_rx5day_neb[0]]))
stats4 = np.concatenate((sort_cri_tx90p_neb,[sort_cri_tx90p_neb[0]]))
stats5 = np.concatenate((sort_cri_wsdi_neb,[sort_cri_wsdi_neb[0]]))
ax.plot(angles, stats1, color='green', linewidth=1, label='CDD')
ax.plot(angles, stats2, color='blue', linewidth=1, label='R95p')
ax.plot(angles, stats3, color='purple', linewidth=1, label='Rx5day')
ax.plot(angles, stats4, color='brown', linewidth=1, label='Tx90p')
ax.plot(angles, stats5, color='orange', linewidth=1, label='WSDI')
# ~ ax.plot(angles, stats1, color='green', linewidth=0.1)
# ~ ax.fill(angles, stats1, color='green', alpha=0.25, label='CDD')
# ~ ax.plot(angles, stats2, color='blue', linewidth=0.1)
# ~ ax.fill(angles, stats2, color='blue', alpha=0.25, label='R95p')
# ~ ax.plot(angles, stats3, color='purple', linewidth=0.1)
# ~ ax.fill(angles, stats3, color='purple', alpha=0.25, label='Rx5day')
# ~ ax.plot(angles, stats4, color='brown', linewidth=0.1)
# ~ ax.fill(angles, stats4, color='brown', alpha=0.25, label='Tx90p')
# ~ ax.plot(angles, stats5, color='orange', linewidth=0.1)
# ~ ax.fill(angles, stats5, color='orange', alpha=0.25, label='WSDI')
ax.set_title(u'(c) NEB', loc='left', fontweight='bold', fontsize=8)
ax.set_thetagrids(angles[:-1] * 180/np.pi, labels, fontsize=8)
ax.grid(True)

ax = fig.add_subplot(2, 3, 4, polar=True)
stats1 = np.concatenate((sort_cri_cdd_sam,[sort_cri_cdd_sam[0]]))
stats2 = np.concatenate((sort_cri_r95p_sam,[sort_cri_r95p_sam[0]]))
stats3 = np.concatenate((sort_cri_rx5day_sam,[sort_cri_rx5day_sam[0]]))
stats4 = np.concatenate((sort_cri_tx90p_sam,[sort_cri_tx90p_sam[0]]))
stats5 = np.concatenate((sort_cri_wsdi_sam,[sort_cri_wsdi_sam[0]]))
ax.plot(angles, stats1, color='green', linewidth=1, label='CDD')
ax.plot(angles, stats2, color='blue', linewidth=1, label='R95p')
ax.plot(angles, stats3, color='purple', linewidth=1, label='Rx5day')
ax.plot(angles, stats4, color='brown', linewidth=1, label='Tx90p')
ax.plot(angles, stats5, color='orange', linewidth=1, label='WSDI')
# ~ ax.plot(angles, stats1, color='green', linewidth=0.1)
# ~ ax.fill(angles, stats1, color='green', alpha=0.25, label='CDD')
# ~ ax.plot(angles, stats2, color='blue', linewidth=0.1)
# ~ ax.fill(angles, stats2, color='blue', alpha=0.25, label='R95p')
# ~ ax.plot(angles, stats3, color='purple', linewidth=0.1)
# ~ ax.fill(angles, stats3, color='purple', alpha=0.25, label='Rx5day')
# ~ ax.plot(angles, stats4, color='brown', linewidth=0.1)
# ~ ax.fill(angles, stats4, color='brown', alpha=0.25, label='Tx90p')
# ~ ax.plot(angles, stats5, color='orange', linewidth=0.1)
# ~ ax.fill(angles, stats5, color='orange', alpha=0.25, label='WSDI')
ax.set_title(u'(d) SAM', loc='left', fontweight='bold', fontsize=8)
ax.set_thetagrids(angles[:-1] * 180/np.pi, labels, fontsize=8)
ax.grid(True)

ax = fig.add_subplot(2, 3, 5, polar=True)
stats1 = np.concatenate((sort_cri_cdd_lpb,[sort_cri_cdd_lpb[0]]))
stats2 = np.concatenate((sort_cri_r95p_lpb,[sort_cri_r95p_lpb[0]]))
stats3 = np.concatenate((sort_cri_rx5day_lpb,[sort_cri_rx5day_lpb[0]]))
stats4 = np.concatenate((sort_cri_tx90p_lpb,[sort_cri_tx90p_lpb[0]]))
stats5 = np.concatenate((sort_cri_wsdi_lpb,[sort_cri_wsdi_lpb[0]]))
ax.plot(angles, stats1, color='green', linewidth=1, label='CDD')
ax.plot(angles, stats2, color='blue', linewidth=1, label='R95p')
ax.plot(angles, stats3, color='purple', linewidth=1, label='Rx5day')
ax.plot(angles, stats4, color='brown', linewidth=1, label='Tx90p')
ax.plot(angles, stats5, color='orange', linewidth=1, label='WSDI')
# ~ ax.plot(angles, stats1, color='green', linewidth=0.1)
# ~ ax.fill(angles, stats1, color='green', alpha=0.25, label='CDD')
# ~ ax.plot(angles, stats2, color='blue', linewidth=0.1)
# ~ ax.fill(angles, stats2, color='blue', alpha=0.25, label='R95p')
# ~ ax.plot(angles, stats3, color='purple', linewidth=0.1)
# ~ ax.fill(angles, stats3, color='purple', alpha=0.25, label='Rx5day')
# ~ ax.plot(angles, stats4, color='brown', linewidth=0.1)
# ~ ax.fill(angles, stats4, color='brown', alpha=0.25, label='Tx90p')
# ~ ax.plot(angles, stats5, color='orange', linewidth=0.1)
# ~ ax.fill(angles, stats5, color='orange', alpha=0.25, label='WSDI')
ax.set_title(u'(e) LPB', loc='left', fontweight='bold', fontsize=8)
ax.set_thetagrids(angles[:-1] * 180/np.pi, labels, fontsize=8)
ax.grid(True)
plt.legend(ncol=6, loc=(-0.85, -0.35), fontsize=10)

ax = fig.add_subplot(2, 3, 6, polar=True)
legend = ['ACCESS-CM2', 'BCC-CSM2-MR', 'CanESM5', 'CMCC-ESM2', 'CNRM-CM6-1', 
'CNRM-ESM2-1', 'GFDL-ESM4', 'INM-CM4-8', 'INM-CM5-0', 'KIOST-ESM', 'MIROC6', 'MIROC-ES2L',
'MPI-ESM1-2-HR', 'MPI-ESM1-2-LR', 'MRI-ESM2-0', 'NESM3', 'NorESM2-MM']
stats1 = np.concatenate((sort_cri_cdd_lpb,[sort_cri_cdd_lpb[0]]))
ax.plot(angles, stats1, color='white', linewidth=0.1)
ax.set_thetagrids(angles[:-1] * 180/np.pi, legend, fontsize=6)
ax.grid(True)

plt.subplots_adjust(hspace=0.02)
plt.subplots_adjust(wspace=0.35)

# Path out to save figure
path_out = '/home/nice/Documentos/paper_mari/figs'
name_out = 'pyplt_rank_etccdi_indices_cmip6.png'
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()
















	


