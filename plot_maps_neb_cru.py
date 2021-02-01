# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "06/03/2019"
__description__ = "This script plot Rclimdex based in HadGEM2-ES model CMIP5"

import os
import netCDF4
import numpy as np
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from netCDF4 import Dataset



def import_accum_database():
	
	path = '/home/nice/Documents/daniele'
	arq  = '{0}/pre_cru_ts4.04_neb_yearsum_1901-2019.nc'.format(path)

	data = netCDF4.Dataset(arq)
	var  = data.variables['pre'][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	year_accum = np.nanmean(np.nanmean(var[:,:,:], axis=1), axis=1) 
	
	return year_accum
	

def import_anom_database():
	
	path = '/home/nice/Documents/daniele'
	arq1  = '{0}/pre_cru_ts4.04_neb_clim_1981-2010.nc'.format(path)

	data1 = netCDF4.Dataset(arq1)
	var1  = data1.variables['pre'][:] 
	lat  = data1.variables['lat'][:]
	lon  = data1.variables['lon'][:]
	clim = np.nanmean(np.nanmean(var1[:,:,:], axis=1), axis=1) 

	arq2  = '{0}/pre_cru_ts4.04_neb_yearavg_1901-2019.nc'.format(path)
	data2 = netCDF4.Dataset(arq2)
	var2  = data2.variables['pre'][:] 
	lat  = data2.variables['lat'][:]
	lon  = data2.variables['lon'][:]
	year_mean = np.nanmean(np.nanmean(var2[:,:,:], axis=1), axis=1)
	
	year_anom = []
	for i in range(0, 119):
		anomaly = year_mean[i] - clim
		year_anom.append(np.squeeze(anomaly))
		
	return year_anom
	
		
# Import obs database 	
pre_accum = import_accum_database()
pre_mean = np.nanmean(pre_accum)

pre_anom = import_anom_database()

# Plot obs database 	
fig = plt.figure(figsize=(14, 7))
time = np.arange(1., 120.)
bar_width = .50

# Subplot one
plt.subplot(211)
plt.bar(time, pre_accum, alpha=1., color='black', width=0.50, edgecolor='darkgray')
plt.ylabel('Precipitação (mm/ano)', fontweight='bold')
plt.axhline(pre_mean, linewidth=1.5, linestyle='dashed', color='black', alpha=1.)
plt.title(u'A) Acumulado anual de precipitação', fontweight='bold', loc='left')
plt.text(116., 1160., u'1122.6', fontweight='bold')
objects = [u'1901', u'1910', u'1919', u'1928', u'1937', u'1946', u'1955', u'1964', u'1973', u'1982', u'1991', u'2000', u'2009', u'2018']
plt.xticks(np.arange(1., 120., 9), objects)

# Subplot two
plt.subplot(212)

norm = matplotlib.colors.Normalize(vmin=-20.1, vmax=20.1)
m = cm.ScalarMappable(norm=norm, cmap=cm.RdBu)
colors = m.to_rgba(pre_anom)

plt.bar(time, pre_anom, alpha=1., width=0.50, color=colors, edgecolor='darkgray')
plt.xlabel('Anos', fontweight='bold')
plt.axvline(29.5, linewidth=1.5, linestyle='dashed', color='black', alpha=1.)
plt.axvline(32.5, linewidth=1.5, linestyle='dashed', color='black', alpha=1.)
plt.axvline(49.5, linewidth=1.5, linestyle='dashed', color='black', alpha=1.)
plt.axvline(55.5, linewidth=1.5, linestyle='dashed', color='black', alpha=1.)
plt.axvline(80.5, linewidth=1.5, linestyle='dashed', color='black', alpha=1.)
plt.axvline(83.5, linewidth=1.5, linestyle='dashed', color='black', alpha=1.)
plt.axvline(89.5, linewidth=1.5, linestyle='dashed', color='black', alpha=1.)
plt.axvline(93.5, linewidth=1.5, linestyle='dashed', color='black', alpha=1.)
plt.axvline(111.5, linewidth=1.5, linestyle='dashed', color='black', alpha=1.)
plt.axvline(116.5, linewidth=1.5, linestyle='dashed', color='black', alpha=1.)
plt.axhline(0, linewidth=1., color='black', alpha=1.)
plt.title(u'B) Anomalia anual de precipitação', fontweight='bold', loc='left')
objects = [u'1901', u'1910', u'1919', u'1928', u'1937', u'1946', u'1955', u'1964', u'1973', u'1982', u'1991', u'2000', u'2009', u'2018']
plt.xticks(np.arange(1., 120., 9), objects)

path_out = '/home/nice/Documents/daniele'
name_out = 'graph_clim_anom_cru_neb.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')

plt.show()
exit()

