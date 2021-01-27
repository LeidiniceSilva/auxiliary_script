# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/29/2020"
__description__ = "This script plot climatology maps from Reg and Had models end obs database"

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


#~ def import_obs(var, area):
	
	#~ path = '/home/nice/Documents/ufrn/papers/eneb_nikolai'	
	#~ arq  = '{0}/{1}_chirps-v2.0-p25_{2}_mm_day_201705_lonlat.nc'.format(path, var, area)	
						
	#~ data = netCDF4.Dataset(arq)
	#~ var  = data.variables[var][:]
	#~ lat  = data.variables['lat'][:]
	#~ lon  = data.variables['lon'][:]
	#~ obs = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)
	
	#~ return obs


#~ def import_rcm(var, area):
	
	#~ path = '/home/nice/Documents/ufrn/papers/eneb_nikolai'	
	#~ arq  = '{0}/{1}_RegCM4.7_{2}_mm_day_201705_lonlat.nc'.format(path, var, area)	
	
	#~ data = netCDF4.Dataset(arq)
	#~ var  = data.variables[var][:]
	#~ lat  = data.variables['lat'][:]
	#~ lon  = data.variables['lon'][:]
	#~ rcm = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)

	#~ return rcm


#~ # Import regcm exp and chirps databases 	   
#~ obs_pre = import_obs('precip', 'eneb')
#~ rcm_pre = import_rcm('pr', 'eneb')

#~ # Plot maps 
#~ fig = plt.figure()
#~ time = np.arange(1, 31 + 1)

#~ ts = plt.plot(time, obs_pre, time, rcm_pre)
#~ l1, l2 = ts
#~ plt.setp(l1, linewidth=1, color='black', marker='.', linestyle='dashed')
#~ plt.setp(l2, linewidth=1, color='blue', marker='+', linestyle='dashed')
#~ plt.title(u'A) 201705 - ENEB', fontweight='bold')
#~ plt.xlabel('Days', fontweight='bold')
#~ plt.ylabel('Rainfall (mm d⁻¹)', fontweight='bold')
#~ plt.grid(True, which='major', linestyle='--')
#~ plt.legend(ts, ['Chirps','Reg'], loc='best', shadow=True, ncol=1)

#~ # Path out to save bias figure
#~ path_out = '/home/nice/Documents/ufrn/papers/eneb_nikolai'	
#~ name_out = 'pyplt_maps_eneb_201705.png'
#~ if not os.path.exists(path_out):
	#~ create_path(path_out)
#~ plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

#~ plt.show()
#~ exit()


def import_obs(var, area):
	
	path = '/home/nice/Documents/ufrn/papers/eneb_nikolai'	
	arq  = '{0}/{1}_chirps-v2.0-p25_{2}_mm_day_201705_lonlat.nc'.format(path, var, area)	
		
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	obs  = var[:][:,:,:]
	
	return lat, lon, obs


def import_rcm(var, area):
	
	path = '/home/nice/Documents/ufrn/papers/eneb_nikolai'	
	arq  = '{0}/{1}_RegCM4.7_{2}_mm_day_201705_lonlat.nc'.format(path, var, area)	
		
	data = netCDF4.Dataset(arq)
	var  = data.variables[var][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	rcm  = var[:][:,:,:]

	return lat, lon, rcm
	

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
	
	map = Basemap(projection='cyl', llcrnrlon=-49.5, llcrnrlat=-19., urcrnrlon=-20.5,urcrnrlat=3.5, resolution='c')
	map.drawmeridians(np.arange(-49.5,-20.5,10.), size=6, labels=[0,0,0,1], linewidth=0.4)
	map.drawparallels(np.arange(-20.5,3.5,5.), size=6, labels=[1,0,0,0], linewidth=0.4)
	
	#~ lons = np.arange(-49.5,-5.,0.25) 
	#~ lats = np.arange(-20.,15.,-0.25) 
	lons, lats = np.meshgrid(new_lon, new_lat)

	xx, yy = map(lons,lats)
	
	path = '/home/nice/Documents/github_projects/shp'
	map.readshapefile('{0}/shp_world/world'.format(path), 'world', drawbounds=True, color='gray', linewidth=.5)
	map.readshapefile('{0}/lim_unid_fed/lim_unid_fed'.format(path), 'lim_unid_fed', drawbounds=True, color='black', linewidth=.5)
	
	return map, xx, yy
	
	
def plot_maps_mean(obs_pre, rcm_pre):
		
	fig = plt.figure()

	levs1 = [0.5, 1, 2, 4, 6, 8, 10, 12, 16, 20]

	ax = fig.add_subplot(2, 4, 1)
	plt.title(u'A) 26/05 Chirps (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')	
	plt.ylabel(u'Latitude', fontsize=6, labelpad=20, fontweight='bold')
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, obs_pre[25,:,:], levels=levs1, latlon=True, cmap=cm.Blues)
	
	ax = fig.add_subplot(2, 4, 2)
	plt.title(u'B) 27/05 Chirps (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')	
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)	
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, obs_pre[26,:,:], levels=levs1, latlon=True, cmap=cm.Blues)

	ax = fig.add_subplot(2, 4, 3)
	plt.title(u'C) 28/05 Chirps (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')	
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, obs_pre[27,:,:], levels=levs1, latlon=True, cmap=cm.Blues) 
	
	ax = fig.add_subplot(2, 4, 4)
	plt.title(u'D) 29/05 Chirps (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')	
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, obs_pre[28,:,:], levels=levs1, latlon=True, cmap=cm.Blues)
	cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	
	ax = fig.add_subplot(2, 4, 5)
	plt.title(u'E) 26/05 Reg (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')	
	plt.ylabel(u'Latitude', fontsize=6, labelpad=20, fontweight='bold')	
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, rcm_pre[25,:,:], levels=levs1, latlon=True, cmap=cm.Blues) 
	
	ax = fig.add_subplot(2, 4, 6)
	plt.title(u'G) 27/05 Reg (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')	
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, rcm_pre[26,:,:], levels=levs1, latlon=True, cmap=cm.Blues)
	
	ax = fig.add_subplot(2, 4, 7)
	plt.title(u'H) 28/05 Reg (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')	
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, rcm_pre[27,:,:], levels=levs1, latlon=True, cmap=cm.Blues) 

	ax = fig.add_subplot(2, 4, 8)
	plt.title(u'H) 29/05 Reg (mm d⁻¹)', fontsize=6, fontweight='bold')
	plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')	
	plt.text(-23, -17, u'\u25B2 \nN ', fontsize=6)
	map, xx, yy = basemap(lat, lon)
	plot_maps_mean = map.contourf(xx, yy, rcm_pre[28,:,:], levels=levs1, latlon=True, cmap=cm.Blues) 
	cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	
	fig.tight_layout()
	
	return plot_maps_mean


# Import regcm exp and chirps databases 	   
lat, lon, obs_pre = import_obs('precip', 'neb')
lat, lon, rcm_pre = import_rcm('pr', 'neb')

# Plot maps with the function
plt_map = plot_maps_mean(obs_pre, rcm_pre)
plt.subplots_adjust(left=0.15, bottom=0.15, right=0.93, top=0.93, wspace=0.30, hspace=0.10)

# Path out to save bias figure
path_out = '/home/nice/Documents/ufrn/papers/eneb_nikolai'	
name_out = 'pyplt_maps_neb_201705.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()


