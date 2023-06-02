# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jun 01, 2023"
__description__ = "This script plot all statistical of cmip6 models"

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

ree_etccdi_naz = np.array([ree_wsdi_naz, ree_tx90p_naz, ree_rx5day_naz, ree_r95p_naz, ree_cdd_naz])
ree_etccdi_saz = np.array([ree_wsdi_saz, ree_tx90p_saz, ree_rx5day_saz, ree_r95p_saz, ree_cdd_saz])
ree_etccdi_neb = np.array([ree_wsdi_neb, ree_tx90p_neb, ree_rx5day_neb, ree_r95p_neb, ree_cdd_neb])
ree_etccdi_sam = np.array([ree_wsdi_sam, ree_tx90p_sam, ree_rx5day_sam, ree_r95p_sam, ree_cdd_sam])
ree_etccdi_lpb = np.array([ree_wsdi_lpb, ree_tx90p_lpb, ree_rx5day_lpb, ree_r95p_lpb, ree_cdd_lpb])

pcc_etccdi_naz = np.array([pcc_wsdi_naz, pcc_tx90p_naz, pcc_rx5day_naz, pcc_r95p_naz, pcc_cdd_naz])
pcc_etccdi_saz = np.array([pcc_wsdi_saz, pcc_tx90p_saz, pcc_rx5day_saz, pcc_r95p_saz, pcc_cdd_saz])
pcc_etccdi_neb = np.array([pcc_wsdi_neb, pcc_tx90p_neb, pcc_rx5day_neb, pcc_r95p_neb, pcc_cdd_neb])
pcc_etccdi_sam = np.array([pcc_wsdi_sam, pcc_tx90p_sam, pcc_rx5day_sam, pcc_r95p_sam, pcc_cdd_sam])
pcc_etccdi_lpb = np.array([pcc_wsdi_lpb, pcc_tx90p_lpb, pcc_rx5day_lpb, pcc_r95p_lpb, pcc_cdd_lpb])

ivs_etccdi_naz = np.array([ivs_wsdi_naz, ivs_tx90p_naz, ivs_rx5day_naz, ivs_r95p_naz, ivs_cdd_naz])
ivs_etccdi_saz = np.array([ivs_wsdi_saz, ivs_tx90p_saz, ivs_rx5day_saz, ivs_r95p_saz, ivs_cdd_saz])
ivs_etccdi_neb = np.array([ivs_wsdi_neb, ivs_tx90p_neb, ivs_rx5day_neb, ivs_r95p_neb, ivs_cdd_neb])
ivs_etccdi_sam = np.array([ivs_wsdi_sam, ivs_tx90p_sam, ivs_rx5day_sam, ivs_r95p_sam, ivs_cdd_sam])
ivs_etccdi_lpb = np.array([ivs_wsdi_lpb, ivs_tx90p_lpb, ivs_rx5day_lpb, ivs_r95p_lpb, ivs_cdd_lpb])

# Plot cmip models and obs database 
fig = plt.figure(figsize=(8, 8))

xlabels = legend
ylabels = ['WSDI', 'Tx90p',  'Rx5day', 'R95p', 'CDD']

norm1 = colors.BoundaryNorm(boundaries=np.arange(-100, 110, 10), ncolors=256)
norm2 = colors.BoundaryNorm(boundaries=np.arange(-1, 1.1, 0.1), ncolors=256)
norm3 = colors.BoundaryNorm(boundaries=np.arange(0, 11, 1), ncolors=256)

color1 = cm.PuOr
color2 = cm.PiYG
color3 = cm.Blues

ax = fig.add_subplot(5, 3, 1)  
pcm = ax.pcolormesh(ree_etccdi_naz, edgecolors='white', linewidths=1., norm=norm1, cmap=color1)
ax.set_title(u'(a) RE', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(ree_etccdi_naz.shape[1]) + 0.5)
ax.set_yticks(np.arange(ree_etccdi_naz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
clb = fig.colorbar(pcm, ax=ax, extend='both', pad=0.03)
clb.ax.yaxis.set_label_position('right')
clb.ax.tick_params(labelsize=8)

ax = fig.add_subplot(5, 3, 2)  
pcm = ax.pcolormesh(pcc_etccdi_naz, edgecolors='white', linewidths=1., norm=norm2, cmap=color2)
ax.set_title(u'(b) PCC', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(pcc_etccdi_naz.shape[1]) + 0.5)
ax.set_yticks(np.arange(pcc_etccdi_naz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)
clb = fig.colorbar(pcm, ax=ax, extend='both', pad=0.03)
clb.ax.yaxis.set_label_position('right')
clb.ax.tick_params(labelsize=8)

ax = fig.add_subplot(5, 3, 3)  
pcm = ax.pcolormesh(ivs_etccdi_naz, edgecolors='white', linewidths=1., norm=norm3, cmap=color3)
ax.set_title(u'(c) IVS', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(ivs_etccdi_naz.shape[1]) + 0.5)
ax.set_yticks(np.arange(ivs_etccdi_naz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)
ax.set_ylabel('NAZ', fontweight='bold', fontsize=8, rotation=-90, labelpad=30)
ax.yaxis.set_label_position("right")
clb = fig.colorbar(pcm, ax=ax, extend='max', pad=0.03)
clb.ax.yaxis.set_label_position('right')
clb.ax.tick_params(labelsize=8)

ax = fig.add_subplot(5, 3, 4)  
pcm = ax.pcolormesh(ree_etccdi_saz, edgecolors='white', linewidths=1., norm=norm1, cmap=color1)
ax.set_title(u'(d)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(ree_etccdi_saz.shape[1]) + 0.5)
ax.set_yticks(np.arange(ree_etccdi_saz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
clb = fig.colorbar(pcm, ax=ax, extend='both', pad=0.03)
clb.ax.yaxis.set_label_position('right')
clb.ax.tick_params(labelsize=8)

ax = fig.add_subplot(5, 3, 5)  
pcm = ax.pcolormesh(pcc_etccdi_saz, edgecolors='white', linewidths=1., norm=norm2, cmap=color2)
ax.set_title(u'(e)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(pcc_etccdi_saz.shape[1]) + 0.5)
ax.set_yticks(np.arange(pcc_etccdi_saz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)
clb = fig.colorbar(pcm, ax=ax, extend='both', pad=0.03)
clb.ax.yaxis.set_label_position('right')
clb.ax.tick_params(labelsize=8)

ax = fig.add_subplot(5, 3, 6)  
pcm = ax.pcolormesh(ivs_etccdi_saz, edgecolors='white', linewidths=1., norm=norm3, cmap=color3)
ax.set_title(u'(f)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(ivs_etccdi_saz.shape[1]) + 0.5)
ax.set_yticks(np.arange(ivs_etccdi_saz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)
ax.set_ylabel('SAZ', fontweight='bold', fontsize=8, rotation=-90, labelpad=30)
ax.yaxis.set_label_position("right")
clb = fig.colorbar(pcm, ax=ax, extend='max', pad=0.03)
clb.ax.yaxis.set_label_position('right')
clb.ax.tick_params(labelsize=8)

ax = fig.add_subplot(5, 3, 7)  
pcm = ax.pcolormesh(ree_etccdi_neb, edgecolors='white', linewidths=1., norm=norm1, cmap=color1)
ax.set_title(u'(g)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(ree_etccdi_neb.shape[1]) + 0.5)
ax.set_yticks(np.arange(ree_etccdi_neb.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
clb = fig.colorbar(pcm, ax=ax, extend='both', pad=0.03)
clb.ax.yaxis.set_label_position('right')
clb.ax.tick_params(labelsize=8)

ax = fig.add_subplot(5, 3, 8)  
pcm = ax.pcolormesh(pcc_etccdi_neb, edgecolors='white', linewidths=1., norm=norm2, cmap=color2)
ax.set_title(u'(h)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(pcc_etccdi_neb.shape[1]) + 0.5)
ax.set_yticks(np.arange(pcc_etccdi_neb.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)
clb = fig.colorbar(pcm, ax=ax, extend='both', pad=0.03)
clb.ax.yaxis.set_label_position('right')
clb.ax.tick_params(labelsize=8)

ax = fig.add_subplot(5, 3, 9)  
pcm = ax.pcolormesh(ivs_etccdi_neb, edgecolors='white', linewidths=1., norm=norm3, cmap=color3)
ax.set_title(u'(i)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(ivs_etccdi_neb.shape[1]) + 0.5)
ax.set_yticks(np.arange(ivs_etccdi_neb.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)
ax.set_ylabel('NEB', fontweight='bold', fontsize=8, rotation=-90, labelpad=30)
ax.yaxis.set_label_position("right")
clb = fig.colorbar(pcm, ax=ax, extend='max', pad=0.03)
clb.ax.yaxis.set_label_position('right')
clb.ax.tick_params(labelsize=8)

ax = fig.add_subplot(5, 3, 10)  
pcm = ax.pcolormesh(ree_etccdi_sam, edgecolors='white', linewidths=1., norm=norm1, cmap=color1)
ax.set_title(u'(j)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(ree_etccdi_sam.shape[1]) + 0.5)
ax.set_yticks(np.arange(ree_etccdi_sam.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
clb = fig.colorbar(pcm, ax=ax, extend='both', pad=0.03)
clb.ax.yaxis.set_label_position('right')
clb.ax.tick_params(labelsize=8)

ax = fig.add_subplot(5, 3, 11)  
pcm = ax.pcolormesh(pcc_etccdi_sam, edgecolors='white', linewidths=1., norm=norm2, cmap=color2)
ax.set_title(u'(k)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(pcc_etccdi_sam.shape[1]) + 0.5)
ax.set_yticks(np.arange(pcc_etccdi_sam.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)
clb = fig.colorbar(pcm, ax=ax, extend='both', pad=0.03)
clb.ax.yaxis.set_label_position('right')
clb.ax.tick_params(labelsize=8)

ax = fig.add_subplot(5, 3, 12)  
pcm = ax.pcolormesh(ivs_etccdi_sam, edgecolors='white', linewidths=1., norm=norm3, cmap=color3)
ax.set_title(u'(l)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(ivs_etccdi_sam.shape[1]) + 0.5)
ax.set_yticks(np.arange(ivs_etccdi_sam.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)
ax.set_ylabel('SAM', fontweight='bold', fontsize=8, rotation=-90, labelpad=30)
ax.yaxis.set_label_position("right")
clb = fig.colorbar(pcm, ax=ax, extend='max', pad=0.03)
clb.ax.yaxis.set_label_position('right')
clb.ax.tick_params(labelsize=8)

ax = fig.add_subplot(5, 3, 13)  
pcm = ax.pcolormesh(ree_etccdi_lpb, edgecolors='white', linewidths=1., norm=norm1, cmap=color1)
ax.set_title(u'(m)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(ree_etccdi_lpb.shape[1]) + 0.5)
ax.set_yticks(np.arange(ree_etccdi_lpb.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
clb = fig.colorbar(pcm, ax=ax, extend='both', pad=0.03)
clb.ax.yaxis.set_label_position('right')
clb.ax.tick_params(labelsize=8)

ax = fig.add_subplot(5, 3, 14)  
pcm = ax.pcolormesh(pcc_etccdi_lpb, edgecolors='white', linewidths=1., norm=norm2, cmap=color2)
ax.set_title(u'(n)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(pcc_etccdi_lpb.shape[1]) + 0.5)
ax.set_yticks(np.arange(pcc_etccdi_lpb.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_yticklabels(), visible=False)
clb = fig.colorbar(pcm, ax=ax, extend='both', pad=0.03)
clb.ax.yaxis.set_label_position('right')
clb.ax.tick_params(labelsize=8)

ax = fig.add_subplot(5, 3, 15)  
pcm = ax.pcolormesh(ivs_etccdi_lpb, edgecolors='white', linewidths=1., norm=norm3, cmap=color3)
ax.set_title(u'(o)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(ivs_etccdi_lpb.shape[1]) + 0.5)
ax.set_yticks(np.arange(ivs_etccdi_lpb.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_yticklabels(), visible=False)
ax.set_ylabel('LPB', fontweight='bold', fontsize=8, rotation=-90, labelpad=30)
ax.yaxis.set_label_position("right")
clb = fig.colorbar(pcm, ax=ax, extend='max', pad=0.03)
clb.ax.yaxis.set_label_position('right')
clb.ax.tick_params(labelsize=8)

plt.subplots_adjust(hspace=0.25)
plt.subplots_adjust(wspace=0.15)

# Path out to save figure
path_out = '/home/nice/Documentos/paper_mari/figs'
name_out = 'pyplt_portrait_etccdi_indices_cmip6.png'
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()

