# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 04, 2024"
__description__ = "This script plot bias maps"

import os
import netCDF4
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from mpl_toolkits.basemap import Basemap

font_size = 8
path='/marconi/home/userexternal/mdasilva'


def import_mean(param, hw):

	arq   = '/home/esp-shared-a/Distribution/Users/epichell/mdasilva/ERA5/{}_era5_{}_mean.nc'.format(param, hw)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['latitude'][:]	
	lon   = data.variables['longitude'][:]
	mean = var[:][:,:,:]
	
	return lat, lon, mean
	

def import_anomaly(param, hw):

	arq   = '/home/esp-shared-a/Distribution/Users/epichell/mdasilva/ERA5/{}_era5_{}_anom.nc'.format(param, hw)	
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['latitude'][:]	
	lon   = data.variables['longitude'][:]
	anomaly = var[:][:,:,:]
	
	return lat, lon, anomaly


def basemap(lat, lon):
	
	
	map = Basemap(projection='cyl', llcrnrlon=180, llcrnrlat=-70., urcrnrlon=360,urcrnrlat=10, resolution='c')
	map.drawmeridians(np.arange(180, 360, 20), size=font_size, labels=[0,0,0,1], linewidth=0.5, color='black')
	map.drawparallels(np.arange(-70, 10, 10), size=font_size, labels=[1,0,0,0], linewidth=0.5, color='black')
	map.drawcoastlines(linewidth=0.5, color='black')
	
	lons, lats = np.meshgrid(lon, lat)
	xx, yy = map(lons,lats)
	
	return map, xx, yy
	
	
# Import model and obs dataset 
period = 'Nov'

if period == 'Sep':
	hw = 'hw_1_17-27Sep2023'
	legend = '17-27 September 2023'
else:
	hw = 'hw_2_11-18Nov2023'
	legend = '11-18 November 2023'

lat, lon, msl_anom = import_anomaly('msl', hw)
lat, lon, t2m_anom = import_anomaly('t2m', hw)
lat, lon, q_anom = import_anomaly('q', hw)
lat, lon, u_anom = import_anomaly('u', hw)
lat, lon, v_anom = import_anomaly('v', hw)
lat, lon, z_anom = import_anomaly('z', hw)

lat, lon, q_mean = import_mean('q', hw)
lat, lon, u_mean = import_mean('u', hw)
lat, lon, v_mean = import_mean('v', hw)
lat, lon, z_mean = import_mean('z', hw)

# Plot figure 
fig = plt.figure(figsize=(10, 6))

# Subplot 1
ax = fig.add_subplot(2, 2, 1)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contour(xx, yy, msl_anom[0,:,:]/100, colors='black', levels=[-15,-10,-5,-1,1,5,10,15], linewidths=0.5)
plt.clabel(plt_map, inline=1, fontsize=font_size)
plt_map = map.contourf(xx, yy, t2m_anom[0,:,:], cmap=cm.bwr, levels=[-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10], extend='both') 
plt.title(u'(a)', loc='left', fontsize=font_size, fontweight='bold')
# Set colobar
cbar = plt.colorbar(plt_map, orientation='horizontal', pad=0.12)
cbar.ax.tick_params(labelsize=font_size)

# Subplot 2
ax = fig.add_subplot(2, 2, 2) 
map, xx, yy = basemap(lat, lon)
plt_map = map.contour(xx, yy, z_mean[0,1,:,:]/10, 6, colors='black', linewidths=0.5)
plt.clabel(plt_map, inline=1, fontsize=font_size)
plt_map = map.contourf(xx, yy, z_anom[0,1,:,:]/10, cmap=cm.PiYG, levels=[-200,-180,-160,-140,-120,-100,-80,-60,-40,-20,20,40,60,80,100,120,140,160,180,200], extend='both') 
plt.title(u'(b)', loc='left', fontsize=font_size, fontweight='bold')
# Set colobar
cbar = plt.colorbar(plt_map, orientation='horizontal', pad=0.12)
cbar.ax.tick_params(labelsize=font_size)

# Subplot 3
ax = fig.add_subplot(2, 2, 3) 
map, xx, yy = basemap(lat, lon)
plt_map = map.contour(xx, yy, u_mean[0,0,:,:], colors='black', linewidths=0.5)
plt.clabel(plt_map, inline=1, fontsize=font_size)
plt_map = map.contourf(xx, yy, u_anom[0,0,:,:], cmap=cm.RdBu_r, levels=[-30,-27,-24,-21,-18,-15,-12,-9,-6,-3,3,6,9,12,15,18,21,24,27,30], extend='both') 
plt.title(u'(c)', loc='left', fontsize=font_size, fontweight='bold')
# Set colobar
cbar = plt.colorbar(plt_map, orientation='horizontal', pad=0.12)
cbar.ax.tick_params(labelsize=font_size)
	
# Subplot 4
ax = fig.add_subplot(2, 2, 4) 
map, xx, yy = basemap(lat, lon)
#plt_map = map.contour(xx, yy, q_mean[0,3,:,:]*1000, colors='black', linewidths=0.5)
#plt.clabel(plt_map, inline=1, fontsize=font_size)
plt_map = map.contourf(xx, yy, q_anom[0,3,:,:]*1000, cmap=cm.BrBG, levels=[-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5], extend='both') 
plt.title(u'(d)', loc='left', fontsize=font_size, fontweight='bold')
# Set colobar
cbar = plt.colorbar(plt_map, orientation='horizontal', pad=0.12)
cbar.ax.tick_params(labelsize=font_size)

plt_map = map.quiver(xx[::25, ::25], yy[::25, ::25], u_mean[0,3,:,:][::25, ::25], v_mean[0,3,:,:][::25, ::25])
plt.quiverkey(plt_map, X=0.9, Y=1.05, U=np.max(u_mean[0,3,:,:]), label='5m/s', labelpos='E', fontproperties={'size': 8})

# Path out to save figure
path_out = '/home/esp-shared-a/Distribution/Users/epichell/mdasilva/ERA5/'
name_out = 'pyplt_maps_anom_era5_{0}.png'.format(hw)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
