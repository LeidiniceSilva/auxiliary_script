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
import matplotlib as mpl
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

	dict_var = {u'prec': u'prec',
	u'mx2t': u'mx2t',
	u'mtnlwrf': u'mtnlwrf',
	u'msl': u'msl',
	u'r': u'ur10',
	u'z': u'z',
	u'uv10': u'u10'}

	print('Open data')
	
	path = '/home/nice/Documents/janio/dados'
	arq  = '{0}/{1}_xavier_br_day_2011-2016_lonlat.nc'.format(path, var)	
	data = netCDF4.Dataset(arq)		
	var  = data.variables[dict_var[var]][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	dataset = var[:][:,:,:]	

	std = np.std(dataset, axis=0)	
	mean = np.nanmean(dataset, axis=0)
	
	# ENOC
	# sumer
	D_ii  = np.nanmean(var[:][1401:1520,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	D_i   = np.nanmean(var[:][1416:1535,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	D     = np.nanmean(var[:][1431:1551,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	Di    = np.nanmean(var[:][1446:1566,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	Dii   = np.nanmean(var[:][1461:1580,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	
	std_D_ii  = np.std(var[:][1401:1520,:,:], axis=0) 
	std_D_i   = np.std(var[:][1416:1535,:,:], axis=0)
	std_D     = np.std(var[:][1431:1551,:,:], axis=0)
	std_Di    = np.std(var[:][1446:1566,:,:], axis=0) 
	std_Dii   = np.std(var[:][1461:1580,:,:], axis=0) 

	mean_D_ii  = np.nanmean(var[:][1401:1520,:,:], axis=0) 
	mean_D_i   = np.nanmean(var[:][1416:1535,:,:], axis=0)
	mean_D     = np.nanmean(var[:][1431:1551,:,:], axis=0)
	mean_Di    = np.nanmean(var[:][1446:1566,:,:], axis=0) 
	mean_Dii   = np.nanmean(var[:][1461:1580,:,:], axis=0) 

	std_clim  = np.std(var[:][:,:,:], axis=0)
	mean_clim  = np.nanmean(var[:][:,:,:], axis=0)

	# autumn
	#~ D_ii  = np.nanmean(var[:][1766:1886,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ D_i   = np.nanmean(var[:][1781:1901,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ D     = np.nanmean(var[:][1795:1926,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ Di    = np.nanmean(var[:][1810:1930,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ Dii   = np.nanmean(var[:][1825:1945,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)

	# EEOC
	# sumer
	#~ D_ii  = np.nanmean(var[:][1401:1520,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ D_i   = np.nanmean(var[:][1416:1535,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ D     = np.nanmean(var[:][1431:1551,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ Di    = np.nanmean(var[:][1446:1566,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ Dii   = np.nanmean(var[:][1461:1580,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)

	# autumn
	#~ D_ii  = np.nanmean(var[:][1766:1886,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ D_i   = np.nanmean(var[:][1781:1901,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ D     = np.nanmean(var[:][1795:1926,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ Di    = np.nanmean(var[:][1810:1930,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ Dii   = np.nanmean(var[:][1825:1945,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	
	return lat, lon, dataset, std_D_ii, std_clim, mean_D_ii, mean_clim, D_ii, D_i, D, Di, Dii

	
def ttest(mean_sample1, mean_sample2, std_sample1, std_sample2):

	# Calculate t statistics	
	p1 = mean_sample1 - mean_sample2 
	p2 = (std_sample1 - std_sample2) / np.sqrt(240)

	ttest = p1 / p2

	# Calculate p value
	p_value = 1 - stats.t.cdf(ttest, df=240)

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
	map.drawmeridians(np.arange(-80.,-30.,10.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
	map.drawparallels(np.arange(-40.,10.,10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')

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
	
	
# Plot maps with the function
print('Plot maps with the function')
fig = plt.figure(figsize=(8, 2))
	
# Import era5 database
variable = u'prec'
event    = u'ENOC'
season   = u'DJF'

dict_unit = {u'prec': u'PRE (mm d⁻¹)',
u'mx2t': u'Tmax (°C)',
u'mtnlwrf': u'ROL (W m⁻²)',
u'msl': u'PNMM (hPa)',
u'r': u'UR (%)',
u'z': u'GEO 500hPa (m)',
u'uv10': u'UV10m (m s⁻¹)'}

lat, lon, var, std_D_ii, std_clim, mean_D_ii, mean_clim, D_ii, D_i, D, Di, Dii = import_era5(variable) 

# PRE
cor_map = mpl.cm.BrBG
levs1  = [-4, -3, -2, -1, 1, 2, 3, 4]
#~ # TEMP
#~ cor_map = mpl.cm.brw
#~ levs1   = [-4, -3, -2, -1, 1, 2, 3, 4]
#~ # ROL
#~ cor_map = mpl.cm.bwr
#~ levs1   = [-40, -30, -20, -10, 10, 20, 30, 40]
#~ # PNMM
#~ cmap    = mpl.cm.PiYG
#~ cor_map = cmap.reversed()
#~ levs1 = [-4, -3, -2, -1, 1, 2, 3, 4]
#~ # UR
#~ cmap = mpl.cm.bwr
#~ cmap_r = cmap.reversed()
#~ levs1 = [-8, -6, -4, -2, 2, 4, 6, 8]
#~ # GEO
#~ cor_map = mpl.cm.PiYG
#~ levs1 = [-40, -30, -20, -10, 10, 20, 30, 40]
#~ # Rol
#~ cor_map = mpl.cm.RdGy
#~ levs1 = [-2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2]

ax = fig.add_subplot(1, 5, 1)
map, xx, yy = basemap(lat, lon)
plt.title(u'A) Evento de {0} (D-2) \n {1} {2}'.format(event, season, dict_unit[variable]), fontsize=6, fontweight='bold')
map.contourf(xx, yy, D_ii, levels=levs1, latlon=True, cmap=cor_map)	
p_value = ttest(std_D_ii, std_clim, mean_D_ii, mean_clim)
p_value = ma.masked_where(p_value >= 0.05, p_value) 
map.contourf(xx, yy, p_value, colors='none', hatches=['....'])

ax = fig.add_subplot(1, 5, 2)
map, xx, yy = basemap(lat, lon)
plt.title(u'B) Evento de ENOC (D-1) \n {1} {2}'.format(event, season, dict_unit[variable]), fontsize=6, fontweight='bold')
map.contourf(xx, yy, D_i, levels=levs1, latlon=True, cmap=cor_map)
p_value = ttest(std_D_ii, std_clim, mean_D_ii, mean_clim)
p_value = ma.masked_where(p_value >= 0.05, p_value) 
map.contourf(xx, yy, p_value, colors='none', hatches=['....'])

ax = fig.add_subplot(1, 5, 3)
map, xx, yy = basemap(lat, lon)
plt.title(u'C) Evento de ENOC (D0) \n {1} {2}'.format(event, season, dict_unit[variable]), fontsize=6, fontweight='bold')
map.contourf(xx, yy, D, levels=levs1, latlon=True, cmap=cor_map)
p_value = ttest(std_D_ii, std_clim, mean_D_ii, mean_clim)
p_value = ma.masked_where(p_value >= 0.05, p_value) 
map.contourf(xx, yy, p_value, colors='none', hatches=['....'])

ax = fig.add_subplot(1, 5, 4)
map, xx, yy = basemap(lat, lon)
plt.title(u'D) Evento de ENOC (D+1) \n {1} {2}'.format(event, season, dict_unit[variable]), fontsize=6, fontweight='bold')
map.contourf(xx, yy, Di, levels=levs1, latlon=True, cmap=cor_map)
p_value = ttest(std_D_ii, std_clim, mean_D_ii, mean_clim)
p_value = ma.masked_where(p_value >= 0.05, p_value) 
map.contourf(xx, yy, p_value, colors='none', hatches=['....'])

ax = fig.add_subplot(1, 5, 5)
map, xx, yy = basemap(lat, lon)
plt.title(u'E) Evento de ENOC (D+2) \n {1} {2}'.format(event, season, dict_unit[variable]), fontsize=6, fontweight='bold')
map.contourf(xx, yy, Dii, levels=levs1, latlon=True, cmap=cor_map)  
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6)
p_value = ttest(std_D_ii, std_clim, mean_D_ii, mean_clim)
p_value = ma.masked_where(p_value >= 0.05, p_value) 
map.contourf(xx, yy, p_value, colors='none', hatches=['....'])

# Path out to save figure
print('Path out to save figure')
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_composites_{0}_{1}_{2}_era5.png'.format(event, season, dict_unit[variable])
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')
plt.show()
exit()	
