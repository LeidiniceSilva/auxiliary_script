# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "05/30/2021"
__description__ = "This script plot compoosites maps from era5"

import os
import conda
import netCDF4
import numpy as np
import numpy.ma as ma
import matplotlib.cm as cm
import scipy.stats as stats
import matplotlib.pyplot as plt
import warnings ; warnings.filterwarnings("ignore")

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from matplotlib.path import Path
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import PathPatch
from scipy.stats import t


def import_era5(var):

	dict_var = {u'msl': u'msl'}

	print('Open data')
	
	path = '/home/nice/Documents/era5'
	arq  = '{0}/{1}_era5_br_day_2011-2020.nc'.format(path, var)	
	data = netCDF4.Dataset(arq)		
	var  = data.variables[dict_var[var]][:]
	lat  = data.variables['latitude'][:]
	lon  = data.variables['longitude'][:]
	dataset = var[:][:,:,:]	
	
	D_ii  = np.mean(var[:][1766:1886,:,:], axis=0) - np.mean(var[:][:,:,:], axis=0)
	D_i   = np.mean(var[:][1781:1901,:,:], axis=0) - np.mean(var[:][:,:,:], axis=0)
	D     = np.mean(var[:][1795:1916,:,:], axis=0) - np.mean(var[:][:,:,:], axis=0)
	Di    = np.mean(var[:][1810:1930,:,:], axis=0) - np.mean(var[:][:,:,:], axis=0)
	Dii   = np.mean(var[:][1825:1945,:,:], axis=0) - np.mean(var[:][:,:,:], axis=0)

	return lat, lon, dataset, D_ii, D_i, D, Di, Dii

	
def ttest(variable, composite):

	print('Compute ttest')

	# Standard error	
	std = np.std(variable, axis=0)
	
	# Compute t statistics
	t = (composite * np.sqrt(121))/std

	# Compute p value
	p_value = 1 - stats.t.cdf(t, df=121)
	
	return p_value
	
		
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
	
	map = Basemap(projection='cyl', llcrnrlat=-40, urcrnrlat=10, llcrnrlon=-80, urcrnrlon=-30, resolution=None, suppress_ticks=True, lon_0=0, celestial=False)
	lons, lats = np.meshgrid(new_lon, new_lat)
	xx, yy = map(lons,lats)
	xin = np.linspace(map.xmin,map.xmax,10) 
	yin = np.linspace(map.ymin,map.ymax,10) 
	lons = np.arange(-80.,-30.,0.25) 
	lats = np.arange(-40.,10.,-0.25) 

	path = '/home/nice/Documents/github_projects/shp'
	map.readshapefile('{0}/shp_world/world'.format(path), 'world', drawbounds=True, color='gray', linewidth=.5)
	map.readshapefile('{0}/lim_unid_fed/lim_unid_fed'.format(path), 'lim_unid_fed', drawbounds=True, color='black', linewidth=.5)
	
	return map, xx, yy
	
	
# Import era5 database
print('Import database')
lat, lon, var, D_ii, D_i, D, Di, Dii = import_era5('msl')   

# Plot maps with the function
print('Plot maps with the function')
fig = plt.figure(figsize=(8, 2))

levs1 = [-4, -3, -2, -1, 1, 2, 3, 4]
#~ levs1 = [-20, -15, -10, -5, 5, 10, 15, 20]
#~ levs1 = [-40, -30, -20, -10, 10, 20, 30, 40]

ax = fig.add_subplot(1, 5, 1)
map, xx, yy = basemap(lat, lon)
plt.title(u'K) Evento de OC (D-2) \n PNMM (hPa)', fontsize=6, fontweight='bold')
map.contourf(xx, yy, D_ii, levels=levs1, latlon=True, cmap=cm.PuOr)	
map.drawmeridians(np.arange(-80.,-30.,10.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-40.,10.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value = ttest(var, D_ii)
p_value = ma.masked_where(p_value >= 0.01, p_value) 
plt.rcParams['hatch.linewidth'] = 0.1
map.contourf(xx, yy, p_value, colors='none', hatches=['....'])

ax = fig.add_subplot(1, 5, 2)
map, xx, yy = basemap(lat, lon)
plt.title(u'L) Evento de OC (D-1) \n PNMM (hPa)', fontsize=6, fontweight='bold')
map.contourf(xx, yy, D_i, levels=levs1, latlon=True, cmap=cm.PuOr)
map.drawmeridians(np.arange(-80.,-30.,10.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-40.,10.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value = ttest(var, D_i)
p_value = ma.masked_where(p_value >= 0.01, p_value) 
map.contourf(xx, yy, p_value, colors='none', hatches=['....'])

ax = fig.add_subplot(1, 5, 3)
map, xx, yy = basemap(lat, lon)
plt.title(u'M) Evento de OC (D0) \n PNMM (hPa)', fontsize=6, fontweight='bold')
map.contourf(xx, yy, D, levels=levs1, latlon=True, cmap=cm.PuOr)
map.drawmeridians(np.arange(-80.,-30.,10.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-40.,10.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value = ttest(var, D)
p_value = ma.masked_where(p_value >= 0.01, p_value) 
map.contourf(xx, yy, p_value, colors='none', hatches=['....'])

ax = fig.add_subplot(1, 5, 4)
map, xx, yy = basemap(lat, lon)
plt.title(u'N) Evento de OC (D+1) \n PNMM (hPa)', fontsize=6, fontweight='bold')
map.contourf(xx, yy, Di, levels=levs1, latlon=True, cmap=cm.PuOr)
map.drawmeridians(np.arange(-80.,-30.,10.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-40.,10.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
p_value = ttest(var, Di)
p_value = ma.masked_where(p_value >= 0.01, p_value) 
map.contourf(xx, yy, p_value, colors='none', hatches=['....'])

ax = fig.add_subplot(1, 5, 5)
map, xx, yy = basemap(lat, lon)
plt.title(u'O) Evento de OC (D+2) \n PNMM (hPa)', fontsize=6, fontweight='bold')
map.contourf(xx, yy, Dii, levels=levs1, latlon=True, cmap=cm.PuOr)  
map.drawmeridians(np.arange(-80.,-30.,10.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-40.,10.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6)
p_value = ttest(var, Dii)
p_value = ma.masked_where(p_value >= 0.01, p_value) 
map.contourf(xx, yy, p_value, colors='none', hatches=['....'])

# Path out to save bias figure
print('Path out to save bias figure')
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_composites_pnmm_oc_era5.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')

plt.show()
exit()	
	
	
	
	
