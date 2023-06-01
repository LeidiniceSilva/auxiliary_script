# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 01, 2023"
__description__ = "This script plot cri of cmip6 models"

import os
import netCDF4
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from dict_cmip6_models_name import cmip6, cmip6_i


def sort_list_reverse(data_list):
	
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
		

def compute_cri(rank1,rank2,rank3,rank4,rank5):
	
	p1 = (rank1+rank2+rank3+rank4+rank5)
	p2 = p1/85
	cri = 1 - p2
	
	return cri
	
	
# Import cmip models and obs database 
var_cmip6 = 'tasmin'
dt = '1986-2005'

mbe_namz_pr = [9, 16, 15, 7, 11, 10, 8, 2, 6, 17, 4, 3, 12, 14, 1, 13, 5]
rmse_namz_pr = [9, 16, 15, 7, 11, 10, 8, 4, 3, 17, 6, 1, 13, 14, 2, 12, 5]
tss_namz_pr = [3, 10, 1, 6, 17, 16, 12, 14, 15, 13, 7, 5, 11, 8, 2, 4, 9]
pcc_namz_pr = [16, 13, 12, 8, 10, 11, 4, 9, 5, 15, 2, 6, 7, 14, 3, 17, 1]
ivs_namz_pr = [2, 11, 8, 14, 13, 10, 1, 16, 15, 17, 3, 12, 6, 7, 4, 5, 9]

mbe_samz_pr = [9, 7, 15, 12, 13, 14, 6, 5, 10, 17, 11, 4, 3, 2, 16, 8, 1]
rmse_samz_pr = [9, 10, 16, 15, 11, 14, 6, 2, 8, 17, 3, 4, 1, 5, 12, 13, 7]
tss_samz_pr = [14, 9, 12, 17, 11, 13, 10, 5, 6, 7, 2, 15, 3, 4, 1, 16, 8]
pcc_samz_pr = [4, 8, 13, 7, 9, 11, 12, 14, 10, 16, 5, 15, 1, 6, 2, 17, 3]
ivs_samz_pr = [15, 12, 8, 5, 4, 2, 6, 3, 1, 17, 11, 10, 13, 9, 14, 16, 7]

mbe_neb_pr = [14, 15, 7, 16, 9, 8, 11, 13, 12, 1, 6, 10, 4, 5, 2, 17, 3]
rmse_neb_pr = [12, 15, 7, 16, 10, 8, 11, 13, 14, 1, 6, 9, 3, 5, 2, 17, 4]
tss_neb_pr = [9, 6, 16, 7, 13, 10, 3, 14, 11, 5, 15, 12, 4, 8, 1, 17, 2]
pcc_neb_pr = [7, 6, 13, 3, 5, 4, 8, 1, 2, 15, 12, 16, 9, 11, 14, 17, 10]
ivs_neb_pr = [16, 4, 5, 15, 11, 2, 7, 14, 1, 17, 10, 8, 12, 6, 9, 3, 13]

mbe_sam_pr = [16, 12, 10, 6, 8, 11, 7, 3, 1, 15, 13, 14, 4, 9, 2, 17, 5]
rmse_sam_pr = [16, 11, 12, 6, 9, 10, 7, 3, 1, 14, 13, 15, 4, 8, 5, 17, 2]
tss_sam_pr = [10, 7, 16, 11, 13, 14, 8, 6, 3, 9, 15, 17, 2, 5, 1, 12, 4]
pcc_sam_pr = [7, 1, 14, 9, 15, 17, 13, 16, 12, 11, 2, 4, 8, 6, 10, 5, 3]
ivs_sam_pr = [12, 4, 7, 2, 1, 11, 13, 9, 3, 17, 16, 15, 8, 6, 14, 5, 10]

mbe_lpb_pr = [16, 13, 3, 9, 15, 14, 11, 12, 10, 17, 2, 7, 4, 1, 6, 5, 8]
rmse_lpb_pr = [17, 11, 2, 3, 14, 15, 8, 12, 10, 16, 6, 13, 1, 4, 5, 7, 9]
tss_lpb_pr = [10, 11, 4, 3, 6, 7, 2, 15, 14, 5, 9, 17, 1, 12, 8, 16, 13]
pcc_lpb_pr = [15, 10, 7, 8, 13, 6, 1, 2, 3, 14, 5, 11, 16, 9, 17, 12, 4]
ivs_lpb_pr = [17, 4, 7, 9, 8, 13, 5, 6, 1, 15, 12, 14, 10, 3, 2, 11, 16]

mbe_br_pr = [14, 5, 16, 6, 13, 15, 1, 11, 12, 17, 9, 7, 2, 4, 10, 8, 3]
rmse_br_pr = [13, 16, 14, 11, 10, 12, 6, 5, 4, 17, 3, 8, 7, 9, 1, 15, 2]
tss_br_pr = [14, 15, 16, 11, 12, 10, 8, 5, 3, 13, 4, 7, 6, 9, 1, 17, 2]
pcc_br_pr = [6, 12, 8, 2, 3, 1, 7, 9, 4, 16, 13, 15, 10, 14, 11, 17, 5]
ivs_br_pr = [11, 1, 7, 14, 13, 12, 3, 15, 8, 17, 9, 5, 4, 2, 10, 16, 6]

mbe_namz_tasmax = [4, 11, 12, 3, 10, 7, 5, 16, 15, 17, 2, 13, 6, 8, 9, 14, 1]
rmse_namz_tasmax = [2, 11, 12, 5, 7, 3, 4, 16, 15, 17, 10, 13, 9, 8, 6, 14, 1]
tss_namz_tasmax = [7, 9, 12, 13, 1, 2, 11, 5, 6, 17, 15, 3, 16, 14, 4, 8, 10]
pcc_namz_tasmax = [16, 11, 15, 6, 9, 7, 2, 8, 4, 14, 3, 12, 5, 13, 10, 17, 1]
ivs_namz_tasmax = [8, 4, 3, 1, 2, 9, 7, 17, 13, 15, 16, 12, 11, 6, 10, 5, 14]

mbe_samz_tasmax = [2, 7, 4, 1, 9, 5, 8, 13, 14, 16, 6, 11, 12, 15, 10, 17, 3]
rmse_samz_tasmax = [1, 8, 5, 4, 7, 3, 6, 13, 14, 16, 9, 11, 12, 15, 10, 17, 2]
tss_samz_tasmax = [11, 13, 15, 16, 2, 1, 7, 5, 6, 17, 10, 8, 4, 9, 3, 12, 14]
pcc_samz_tasmax = [16, 5, 14, 2, 1, 3, 10, 7, 9, 4, 8, 6, 12, 15, 13, 17, 11]
ivs_samz_tasmax = [12, 3, 13, 9, 5, 6, 10, 11, 4, 14, 16, 17, 7, 2, 8, 1, 15]

mbe_neb_tasmax = [7, 2, 3, 6, 13, 9, 12, 15, 14, 11, 5, 4, 10, 16, 1, 17, 8]
rmse_neb_tasmax = [8, 2, 3, 5, 12, 6, 11, 15, 13, 14, 9, 4, 10, 16, 1, 17, 7]
tss_neb_tasmax = [8, 7, 10, 9, 1, 3, 4, 15, 16, 17, 14, 12, 5, 11, 2, 13, 6]
pcc_neb_tasmax = [9, 14, 10, 7, 13, 11, 8, 4, 6, 16, 17, 15, 12, 1, 5, 2, 3]
ivs_neb_tasmax = [7, 9, 5, 2, 12, 6, 4, 16, 10, 14, 17, 15, 1, 3, 8, 11, 13]

mbe_sam_tasmax = [8, 2, 3, 4, 7, 5, 13, 12, 11, 14, 1, 10, 15, 16, 6, 17, 9]
rmse_sam_tasmax = [8, 1, 2, 4, 7, 3, 12, 11, 10, 15, 5, 13, 14, 16, 6, 17, 9]
tss_sam_tasmax = [8, 11, 12, 10, 5, 7, 2, 14, 13, 15, 16, 17, 1, 3, 9, 4, 6]
pcc_sam_tasmax = [3, 13, 10, 5, 8, 9, 4, 7, 14, 11, 2, 1, 6, 15, 12, 17, 16]
ivs_sam_tasmax = [6, 3, 7, 11, 10, 12, 14, 4, 9, 17, 15, 16, 2, 1, 8, 5, 13]

mbe_lpb_tasmax = [12, 9, 7, 5, 3, 1, 10, 2, 6, 14, 11, 13, 15, 16, 8, 17, 4]
rmse_lpb_tasmax = [11, 9, 4, 2, 1, 3, 10, 7, 8, 14, 12, 13, 15, 16, 6, 17, 5]
tss_lpb_tasmax = [7, 10, 6, 4, 3, 5, 2, 15, 14, 9, 11, 17, 8, 13, 1, 16, 12]
pcc_lpb_tasmax = [2, 7, 14, 11, 12, 13, 8, 16, 15, 17, 10, 5, 3, 1, 6, 9, 4]
ivs_lpb_tasmax = [1, 8, 12, 6, 15, 16, 5, 7, 13, 17, 11, 10, 14, 9, 3, 2, 4]

mbe_br_tasmax = [5, 7, 3, 2, 10, 6, 11, 14, 13, 16, 1, 9, 12, 15, 8, 17, 4]
rmse_br_tasmax = [5, 7, 6, 1, 8, 2, 9, 14, 13, 16, 10, 11, 12, 15, 4, 17, 3]
tss_br_tasmax = [7, 10, 9, 6, 4, 3, 2, 15, 14, 13, 16, 17, 8, 11, 1, 12, 5]
pcc_br_tasmax = [7, 8, 3, 12, 13, 11, 6, 14, 16, 17, 9, 15, 2, 1, 4, 10, 5]
ivs_br_tasmax = [12, 1, 8, 11, 4, 7, 9, 15, 10, 13, 17, 16, 6, 2, 3, 5, 14]

mbe_namz_tasmin = [7, 5, 17, 9, 12, 6, 1, 3, 2, 11, 15, 14, 13, 10, 8, 16, 4]
rmse_namz_tasmin = [6, 7, 17, 9, 11, 4, 5, 3, 1, 10, 15, 13, 14, 12, 8, 16, 2]
tss_namz_tasmin = [3, 6, 13, 8, 2, 1, 9, 12, 10, 5, 16, 15, 17, 14, 4, 11, 7]
pcc_namz_tasmin = [8, 7, 9, 6, 3, 4, 1, 10, 12, 16, 15, 14, 11, 13, 5, 17, 2]
ivs_namz_tasmin = [3, 9, 2, 4, 7, 12, 11, 16, 17, 5, 15, 14, 8, 6, 1, 10, 13]

mbe_samz_tasmin = [8, 12, 17, 13, 7, 1, 2, 3, 4, 10, 15, 14, 11, 5, 9, 16, 6]
rmse_samz_tasmin = [8, 12, 17, 15, 6, 1, 2, 4, 5, 11, 13, 14, 9, 3, 10, 16, 7]
tss_samz_tasmin = [10, 6, 11, 17, 2, 5, 4, 13, 14, 9, 3, 16, 1, 7, 12, 8, 15]
pcc_samz_tasmin = [2, 13, 14, 7, 8, 11, 10, 4, 5, 17, 6, 16, 12, 9, 3, 15, 1]
ivs_samz_tasmin = [4, 9, 3, 10, 11, 12, 6, 13, 16, 7, 14, 17, 2, 8, 5, 1, 15]

mbe_neb_tasmin = [9, 7, 15, 17, 10, 5, 6, 3, 1, 4, 16, 12, 8, 2, 13, 14, 11]
rmse_neb_tasmin = [10, 8, 15, 17, 9, 1, 2, 4, 5, 6, 16, 12, 7, 3, 13, 14, 11]
tss_neb_tasmin = [7, 9, 10, 17, 1, 2, 3, 15, 16, 13, 6, 14, 4, 11, 12, 8, 5]
pcc_neb_tasmin = [7, 10, 4, 3, 2, 6, 1, 15, 13, 12, 17, 8, 14, 11, 5, 16, 9]
ivs_neb_tasmin = [10, 4, 1, 13, 8, 5, 3, 9, 14, 12, 16, 17, 7, 11, 2, 6, 15]

mbe_sam_tasmin = [6, 14, 15, 17, 8, 3, 4, 1, 2, 9, 12, 10, 7, 5, 13, 16, 11]
rmse_sam_tasmin = [7, 14, 15, 17, 8, 2, 1, 4, 3, 9, 13, 11, 6, 5, 12, 16, 10]
tss_sam_tasmin = [8, 15, 13, 10, 4, 3, 2, 16, 14, 11, 9, 17, 1, 6, 5, 7, 12]
pcc_sam_tasmin = [4, 11, 13, 14, 15, 16, 12, 2, 3, 17, 8, 5, 9, 1, 6, 10, 7]
ivs_sam_tasmin = [11, 3, 1, 14, 10, 6, 12, 2, 4, 13, 15, 17, 8, 7, 5, 9, 16]

mbe_lpb_tasmin = [4, 13, 14, 10, 8, 7, 6, 2, 3, 12, 15, 17, 5, 1, 9, 16, 11]
rmse_lpb_tasmin = [4, 13, 15, 12, 5, 3, 2, 8, 9, 11, 14, 17, 1, 6, 7, 16, 10]
tss_lpb_tasmin = [11, 7, 14, 10, 2, 3, 4, 15, 16, 9, 8, 17, 1, 13, 5, 12, 6]
pcc_lpb_tasmin = [1, 4, 12, 6, 14, 16, 10, 7, 8, 17, 15, 13, 3, 2, 9, 11, 5]
ivs_lpb_tasmin = [1, 12, 9, 8, 16, 17, 13, 6, 3, 10, 7, 5, 4, 2, 15, 14, 11]

mbe_br_tasmin = [6, 10, 16, 13, 8, 5, 3, 1, 2, 12, 15, 14, 9, 4, 11, 17, 7]
rmse_br_tasmin = [6, 11, 16, 14, 7, 1, 2, 3, 4, 12, 13, 15, 8, 5, 10, 17, 9]
tss_br_tasmin = [6, 13, 15, 16, 1, 2, 3, 11, 10, 7, 8, 17, 4, 9, 5, 12, 14]
pcc_br_tasmin = [4, 9, 12, 6, 11, 13, 7, 3, 5, 17, 16, 15, 8, 10, 1, 14, 2]
ivs_br_tasmin = [6, 7, 1, 12, 5, 9, 10, 13, 15, 8, 14, 17, 4, 3, 2, 11, 16]

cri_namz_pr = []
cri_samz_pr = []
cri_neb_pr = []
cri_sam_pr = []
cri_lpb_pr = []
cri_br_pr = []

cri_namz_tasmax = []
cri_samz_tasmax = []
cri_neb_tasmax = []
cri_sam_tasmax = []
cri_lpb_tasmax = []
cri_br_tasmax = []

cri_namz_tasmin = []
cri_samz_tasmin = []
cri_neb_tasmin = []
cri_sam_tasmin = []
cri_lpb_tasmin = []
cri_br_tasmin = []

for i in range(0, 17):

	cri_namz_pr.append(compute_cri(mbe_namz_pr[i],rmse_namz_pr[i],tss_namz_pr[i],pcc_namz_pr[i],ivs_namz_pr[i]))
	cri_samz_pr.append(compute_cri(mbe_samz_pr[i],rmse_samz_pr[i],tss_samz_pr[i],pcc_samz_pr[i],ivs_samz_pr[i]))
	cri_neb_pr.append(compute_cri(mbe_neb_pr[i],rmse_neb_pr[i],tss_neb_pr[i],pcc_neb_pr[i],ivs_neb_pr[i]))
	cri_sam_pr.append(compute_cri(mbe_sam_pr[i],rmse_sam_pr[i],tss_sam_pr[i],pcc_sam_pr[i],ivs_sam_pr[i]))
	cri_lpb_pr.append(compute_cri(mbe_lpb_pr[i],rmse_lpb_pr[i],tss_lpb_pr[i],pcc_lpb_pr[i],ivs_lpb_pr[i]))
	cri_br_pr.append(compute_cri(mbe_br_pr[i],rmse_br_pr[i],tss_br_pr[i],pcc_br_pr[i],ivs_br_pr[i]))

	cri_namz_tasmax.append(compute_cri(mbe_namz_tasmax[i],rmse_namz_tasmax[i],tss_namz_tasmax[i],pcc_namz_tasmax[i],ivs_namz_tasmax[i]))
	cri_samz_tasmax.append(compute_cri(mbe_samz_tasmax[i],rmse_samz_tasmax[i],tss_samz_tasmax[i],pcc_samz_tasmax[i],ivs_samz_tasmax[i]))
	cri_neb_tasmax.append(compute_cri(mbe_neb_tasmax[i],rmse_neb_tasmax[i],tss_neb_tasmax[i],pcc_neb_tasmax[i],ivs_neb_tasmax[i]))
	cri_sam_tasmax.append(compute_cri(mbe_sam_tasmax[i],rmse_sam_tasmax[i],tss_sam_tasmax[i],pcc_sam_tasmax[i],ivs_sam_tasmax[i]))
	cri_lpb_tasmax.append(compute_cri(mbe_lpb_tasmax[i],rmse_lpb_tasmax[i],tss_lpb_tasmax[i],pcc_lpb_tasmax[i],ivs_lpb_tasmax[i]))
	cri_br_tasmax.append(compute_cri(mbe_br_tasmax[i],rmse_br_tasmax[i],tss_br_tasmax[i],pcc_br_tasmax[i],ivs_br_tasmax[i]))

	cri_namz_tasmin.append(compute_cri(mbe_namz_tasmin[i],rmse_namz_tasmin[i],tss_namz_tasmin[i],pcc_namz_tasmin[i],ivs_namz_tasmin[i]))
	cri_samz_tasmin.append(compute_cri(mbe_samz_tasmin[i],rmse_samz_tasmin[i],tss_samz_tasmin[i],pcc_samz_tasmin[i],ivs_samz_tasmin[i]))
	cri_neb_tasmin.append(compute_cri(mbe_neb_tasmin[i],rmse_neb_tasmin[i],tss_neb_tasmin[i],pcc_neb_tasmin[i],ivs_neb_tasmin[i]))
	cri_sam_tasmin.append(compute_cri(mbe_sam_tasmin[i],rmse_sam_tasmin[i],tss_sam_tasmin[i],pcc_sam_tasmin[i],ivs_sam_tasmin[i]))
	cri_lpb_tasmin.append(compute_cri(mbe_lpb_tasmin[i],rmse_lpb_tasmin[i],tss_lpb_tasmin[i],pcc_lpb_tasmin[i],ivs_lpb_tasmin[i]))
	cri_br_tasmin.append(compute_cri(mbe_br_tasmin[i],rmse_br_tasmin[i],tss_br_tasmin[i],pcc_br_tasmin[i],ivs_br_tasmin[i]))

if var_cmip6 == 'pr':
	sort_cri_namz = sort_list_reverse(cri_namz_pr)
	sort_cri_samz = sort_list_reverse(cri_samz_pr)
	sort_cri_neb = sort_list_reverse(cri_neb_pr)
	sort_cri_sam = sort_list_reverse(cri_sam_pr)
	sort_cri_lpb = sort_list_reverse(cri_lpb_pr)
	sort_cri_br = sort_list_reverse(cri_br_pr)
elif var_cmip6 == 'tasmax':
	sort_cri_namz = sort_list_reverse(cri_namz_tasmax)
	sort_cri_samz = sort_list_reverse(cri_samz_tasmax)
	sort_cri_neb = sort_list_reverse(cri_neb_tasmax)
	sort_cri_sam = sort_list_reverse(cri_sam_tasmax)
	sort_cri_lpb = sort_list_reverse(cri_lpb_tasmax)
	sort_cri_br = sort_list_reverse(cri_br_tasmax)
else:
	sort_cri_namz = sort_list_reverse(cri_namz_tasmin)
	sort_cri_samz = sort_list_reverse(cri_samz_tasmin)
	sort_cri_neb = sort_list_reverse(cri_neb_tasmin)
	sort_cri_sam = sort_list_reverse(cri_sam_tasmin)
	sort_cri_lpb = sort_list_reverse(cri_lpb_tasmin)
	sort_cri_br = sort_list_reverse(cri_br_tasmin)

# Plot cmip models and obs database 
fig = plt.figure(figsize=(9, 7))

if var_cmip6 == 'pr':
	cm = plt.get_cmap('jet')
	colors = cm(np.linspace(0.1, 0.9, 17))
elif var_cmip6 == 'tasmax':
	cm = plt.get_cmap('jet')
	colors = cm(np.linspace(0.1, 0.9, 17))
else:
	cm = plt.get_cmap('jet')
	colors = cm(np.linspace(0.1, 0.9, 17))

ax = fig.add_subplot(231, polar=True)
theta = np.arange(0, 2*np.pi, 2*np.pi/len(sort_cri_namz)) 
bars = ax.bar(theta, sort_cri_namz, color=colors, width=0.4)
ax.set_title(u'(a) NAMZ', loc='left', fontsize=8, fontweight='bold')
ax.set_xticks(theta)
ax.set_xticklabels(range(1, len(theta)+1), fontsize=8)
ax.yaxis.grid(True)

ax = fig.add_subplot(232, polar=True)
theta = np.arange(0, 2*np.pi, 2*np.pi/len(sort_cri_samz)) 
bars = bars = ax.bar(theta, sort_cri_samz, color=colors, width=0.4)
ax.set_title(u'(b) SAMZ', loc='left', fontsize=8, fontweight='bold')
ax.set_xticks(theta)
ax.set_xticklabels(range(1, len(theta)+1), fontsize=8)
ax.yaxis.grid(True)

ax = fig.add_subplot(233, polar=True)
theta = np.arange(0, 2*np.pi, 2*np.pi/len(sort_cri_neb)) 
bars = ax.bar(theta, sort_cri_neb, color=colors, width=0.4)
ax.set_title(u'(c) NEB', loc='left', fontsize=8, fontweight='bold')
ax.set_xticks(theta)
ax.set_xticklabels(range(1, len(theta)+1), fontsize=8)
ax.yaxis.grid(True)

ax = fig.add_subplot(234, polar=True)
theta = np.arange(0, 2*np.pi, 2*np.pi/len(sort_cri_sam)) 
bars = ax.bar(theta, sort_cri_sam, color=colors, width=0.4)
ax.set_title(u'(d) SAM', loc='left', fontsize=8, fontweight='bold')
ax.set_xticks(theta)
ax.set_xticklabels(range(1, len(theta)+1), fontsize=8)
ax.yaxis.grid(True)

ax = fig.add_subplot(235, polar=True)
theta = np.arange(0, 2*np.pi, 2*np.pi/len(sort_cri_lpb)) 
bars = ax.bar(theta, sort_cri_lpb, color=colors, width=0.4)
ax.set_title(u'(e) LPB', loc='left', fontsize=8, fontweight='bold')
ax.set_xticks(theta)
ax.set_xticklabels(range(1, len(theta)+1), fontsize=8)
ax.yaxis.grid(True)

ax = fig.add_subplot(236, polar=True)
theta = np.arange(0, 2*np.pi, 2*np.pi/len(sort_cri_br)) 
bars = ax.bar(theta, sort_cri_br, color=colors, width=0.4)
ax.set_title(u'(f) BR', loc='left', fontsize=8, fontweight='bold')
ax.set_xticks(theta)
ax.set_xticklabels(range(1, len(theta)+1), fontsize=8)
ax.yaxis.grid(True)

legend = ['ACCESS-CM2', 'BCC-CSM2-MR', 'CanESM5', 'CMCC-ESM2', 'CNRM-CM6-1', 
'CNRM-ESM2-1', 'GFDL-ESM4', 'INM-CM4-8', 'INM-CM5-0', 'KIOST-ESM', 'MIROC6', 'MIROC-ES2L',
'MPI-ESM1-2-HR', 'MPI-ESM1-2-LR', 'MRI-ESM2-0', 'NESM3', 'NorESM2-MM']
plt.legend(bars, legend, ncol=6, loc=(-2.53, -0.43), fontsize=8)

# Path out to save figure
path_out = '/home/nice/Documentos/AdaptaBrasil_MCTI/figs/figs_report-II'
name_out = 'pyplt_radar_rank_cmip6_{0}_{1}.png'.format(var_cmip6, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()















	


