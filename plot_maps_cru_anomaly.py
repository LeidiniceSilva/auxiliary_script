# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "06/03/2019"
__description__ = "This script plot climatology maps from obs database"

import os
import conda
import netCDF4
import numpy as np
import matplotlib as mpl 
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# mpl.use('Agg')

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from os.path import expanduser
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import BoundaryNorm
from sklearn.linear_model import LinearRegression
from comp_statist_indices import compute_added_value


def import_obs(var):
	
	path = '/home/nice/Downloads'
	arq  = '{0}/cru_ts4.04.1901.2019.{1}.dat_lonlat.nc'.format(path, var)	
			
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	
	clim_8110  = np.nanmean(var[:][960:1319,:,:], axis=0)
	clim_0119  = np.nanmean(var[:][:,:,:], axis=0)
	clim_anomaly = clim_8110 - clim_0119
	
	value_8110 = np.nanmean(np.nanmean(var[960:1319,:,:], axis=1), axis=1)
	value_0119 = np.nanmean(np.nanmean(var[:,:,:], axis=1), axis=1) 

	year_anomaly = []
	for i in range(0, 1428):
		anomaly = value_0119[i] - value_8110
		year_anomaly.append(anomaly)

	return lat, lon, clim_anomaly, year_anomaly


def import_function_trend(data):
	
	data_x = [i for i in range(0, len(data))]
	data_x = np.reshape(data_x, (len(data_x), 1))
	data_y = data
	model = LinearRegression()
	model.fit(data_x, data_y)
	data_trend = model.predict(data_x)
	data_z = np.polyfit(data_y.flatten(), data_trend.flatten(), 1)
	data_r2 = r2_score(data_y, data_trend)
	data_median = statistics.median(data_y)

	return data_trend, data_median
		

def basemap(lat, lon):
	
	aux_lon1 = []
	aux_lon2 = []
	for l in lon:
		if l <= 180:
			aux_lon1.append(l)
		else:
			aux_lon2.append(l-360)
		
	lon = np.array(aux_lon1[::-1] + aux_lon2[::-1])
	new_lat = lat
	new_lon = lon[::-1]
	
	map = Basemap(projection='cyl', llcrnrlon=-85., llcrnrlat=-20., urcrnrlon=-15.,urcrnrlat=10., resolution='c')
	map.drawmeridians(np.arange(-85.,-5.,20.), size=6, labels=[0,0,0,1], linewidth=0.4)
	map.drawparallels(np.arange(-20.,15.,10.), size=6, labels=[1,0,0,0], linewidth=0.4)
	
	xin = np.linspace(map.xmin,map.xmax,20) 
	yin = np.linspace(map.ymin,map.ymax,20) 
	lons = np.arange(-85.,-5.,0.25) 
	lats = np.arange(-20.,15.,-0.25) 
	lons, lats = np.meshgrid(new_lon, new_lat)

	xx, yy = map(lons,lats)
	
	path = '/home/nice/Documents/github_projects/shp'
	map.readshapefile('{0}/shp_world/world'.format(path), 'world', drawbounds=True, color='gray', linewidth=.5)
	map.readshapefile('{0}/lim_unid_fed/lim_unid_fed'.format(path), 'lim_unid_fed', drawbounds=True, color='black', linewidth=.5)
	
	return map, xx, yy
	
	
def plot_maps_mean(clim_anomaly_pre, anomaly_pre, anomaly_pre_trend, anomaly_pre_median, clim_anomaly_tas, anomaly_tas):
	
	fig = plt.figure(figsize=(8, 4))

	levs1 = [-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10]
	levs2 = [-1, -0.60, -0.40, -0.20, -0.10, 0, 0.10, 0.20, 0.40, 0.60, 1]
	time = np.arange(1., 1429.)
	bar_width = .50
	
	ax = fig.add_subplot(2, 2, 1)
	plt.title(u'A) CRU 1981-2010 (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.ylabel(u'Latitude', fontsize=6, labelpad=15, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, clim_anomaly_pre, levels=levs1, latlon=True, cmap=cm.BrBG)
	cbar = map.colorbar(ticks=levs1,  drawedges=True, ax=ax, extend='both', shrink=0.8)
	cbar.ax.tick_params(labelsize=6)

	ax = fig.add_subplot(2, 2, 2)
	norm = mpl.colors.Normalize(vmin=-75, vmax=75)
	m = cm.ScalarMappable(norm=norm, cmap=cm.BrBG)
	colors = m.to_rgba(anomaly_pre)
	plt.bar(time, anomaly_pre, alpha=1., width=0.50, color=colors)
	plt.axhline(0, linewidth=1., color='black', alpha=1.)

	ax = fig.add_subplot(2, 2, 3)
	plt.title(u'B) CRU 1981-2010 (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, clim_anomaly_tas, levels=levs2, latlon=True, cmap=cm.seismic)
	cbar = map.colorbar(ticks=levs2,  drawedges=True, ax=ax, extend='both', shrink=0.8)
	cbar.ax.tick_params(labelsize=6)
	
	ax = fig.add_subplot(2, 2, 4)
	norm = mpl.colors.Normalize(vmin=-1.75, vmax=1.75)
	m = cm.ScalarMappable(norm=norm, cmap=cm.seismic)
	colors = m.to_rgba(anomaly_tas)
	plt.bar(time, anomaly_tas, alpha=1., width=0.50, color=colors)
	plt.axhline(0, linewidth=1., color='black', alpha=1.)

	return plot_maps_mean


# Import regcm exp and cru databases 	   
lat, lon, clim_anomaly_pre, year_anomaly_pre = import_obs('pre')
lat, lon, clim_anomaly_tas, year_anomaly_tas = import_obs('tmp')
anomaly_pre = np.nanmean(year_anomaly_pre, axis=1)
anomaly_tas = np.nanmean(year_anomaly_tas, axis=1)

anomaly_pre_trend, anomaly_pre_median = import_function_trend(year_anomaly_pre)

plt_map = plot_maps_mean(clim_anomaly_pre, anomaly_pre, anomaly_pre_trend, anomaly_pre_median, clim_anomaly_tas, anomaly_tas)
plt.subplots_adjust(left=0.15, bottom=0.15, right=0.93, top=0.93, wspace=0.10, hspace=0.10)

# Path out to save bias figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_anomaly_cru.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')
plt.show()
exit()



