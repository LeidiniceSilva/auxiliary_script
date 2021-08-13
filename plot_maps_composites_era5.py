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
	u'q': u'q',
	u'z': u'z',
	u'uv10': u'u10'}

	path = '/home/nice/Documents/janio/dados'
	arq  = '{0}/{1}_xavier_br_day_2011-2015_lonlat.nc'.format(path, var)	
	data = netCDF4.Dataset(arq)		
	var  = data.variables[dict_var[var]][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]

	std_clim  = np.std(var[:][:,:,:], axis=0)
	mean_clim  = np.nanmean(var[:][:,:,:], axis=0)
	
	# summer
	# ENOC
	D_ii  = np.nanmean(var[:][1046:1135,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	D_i   = np.nanmean(var[:][1056:1145,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	D     = np.nanmean(var[:][1066:1155,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	Di    = np.nanmean(var[:][1076:1165,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	Dii   = np.nanmean(var[:][1086:1175,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)

	std_D_ii  = np.std(var[:][1046:1135,:,:], axis=0) 
	std_D_i   = np.std(var[:][1056:1145,:,:], axis=0)
	std_D     = np.std(var[:][1066:1155,:,:], axis=0)
	std_Di    = np.std(var[:][1076:1165,:,:], axis=0) 
	std_Dii   = np.std(var[:][1086:1175,:,:], axis=0) 
	
	mean_D_ii  = np.nanmean(var[:][1046:1135,:,:], axis=0) 
	mean_D_i   = np.nanmean(var[:][1056:1145,:,:], axis=0)
	mean_D     = np.nanmean(var[:][1066:1155,:,:], axis=0)
	mean_Di    = np.nanmean(var[:][1076:1165,:,:], axis=0) 
	mean_Dii   = np.nanmean(var[:][1086:1175,:,:], axis=0) 

	# EEOC
	#~ D_ii  = np.nanmean(var[:][1411:1500,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ D_i   = np.nanmean(var[:][1426:1510,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ D     = np.nanmean(var[:][1431:1520,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ Di    = np.nanmean(var[:][1446:1530,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ Dii   = np.nanmean(var[:][1451:1540,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)

	#~ std_D_ii  = np.std(var[:][1411:1500,:,:], axis=0) 
	#~ std_D_i   = np.std(var[:][1426:1510,:,:], axis=0)
	#~ std_D     = np.std(var[:][1431:1520,:,:], axis=0)
	#~ std_Di    = np.std(var[:][1446:1530,:,:], axis=0) 
	#~ std_Dii   = np.std(var[:][1451:1540,:,:], axis=0) 

	#~ mean_D_ii  = np.nanmean(var[:][1411:1500,:,:], axis=0) 
	#~ mean_D_i   = np.nanmean(var[:][1426:1510,:,:], axis=0)
	#~ mean_D     = np.nanmean(var[:][1431:1520,:,:], axis=0)
	#~ mean_Di    = np.nanmean(var[:][1446:1530,:,:], axis=0) 
	#~ mean_Dii   = np.nanmean(var[:][1451:1540,:,:], axis=0) 

	# autumn
	# ENOC
	#~ D_ii  = np.nanmean(var[:][1136:1227,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ D_i   = np.nanmean(var[:][1146:1237,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ D     = np.nanmean(var[:][1156:1247,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ Di    = np.nanmean(var[:][1166:1257,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ Dii   = np.nanmean(var[:][1176:1267,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)

	#~ std_D_ii  = np.std(var[:][1136:1227,:,:], axis=0) 
	#~ std_D_i   = np.std(var[:][1146:1237,:,:], axis=0)
	#~ std_D     = np.std(var[:][1156:1247,:,:], axis=0)
	#~ std_Di    = np.std(var[:][1166:1257,:,:], axis=0) 
	#~ std_Dii   = np.std(var[:][1086:1175,:,:], axis=0) 

	#~ mean_D_ii  = np.nanmean(var[:][1136:1227,:,:], axis=0) 
	#~ mean_D_i   = np.nanmean(var[:][1146:1237,:,:], axis=0)
	#~ mean_D     = np.nanmean(var[:][1156:1247,:,:], axis=0)
	#~ mean_Di    = np.nanmean(var[:][1166:1257,:,:], axis=0) 
	#~ mean_Dii   = np.nanmean(var[:][1176:1267,:,:], axis=0) 

	# EEOC
	#~ D_ii  = np.nanmean(var[:][1500:1602,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ D_i   = np.nanmean(var[:][1510:1602,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ D     = np.nanmean(var[:][1520:1612,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ Di    = np.nanmean(var[:][1530:1622,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ Dii   = np.nanmean(var[:][1540:1632,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)

	#~ std_D_ii  = np.std(var[:][1500:1602,:,:], axis=0) 
	#~ std_D_i   = np.std(var[:][1510:1602,:,:], axis=0)
	#~ std_D     = np.std(var[:][1520:1612,:,:], axis=0)
	#~ std_Di    = np.std(var[:][1530:1622,:,:], axis=0) 
	#~ std_Dii   = np.std(var[:][1540:1632,:,:], axis=0) 

	#~ mean_D_ii  = np.nanmean(var[:][1500:1602,:,:], axis=0) 
	#~ mean_D_i   = np.nanmean(var[:][1510:1602,:,:], axis=0)
	#~ mean_D     = np.nanmean(var[:][1520:1612,:,:], axis=0)
	#~ mean_Di    = np.nanmean(var[:][1530:1622,:,:], axis=0) 
	#~ mean_Dii   = np.nanmean(var[:][1540:1632,:,:], axis=0) 

	# spring
	# ENOC
	#~ D_ii  = np.nanmean(var[:][1330:1420,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ D_i   = np.nanmean(var[:][1335:1425,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ D     = np.nanmean(var[:][1340:1430,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ Di    = np.nanmean(var[:][1345:1435,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ Dii   = np.nanmean(var[:][1350:1440,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)

	#~ std_D_ii  = np.std(var[:][1330:1420,:,:], axis=0) 
	#~ std_D_i   = np.std(var[:][1335:1425,:,:], axis=0)
	#~ std_D     = np.std(var[:][1340:1430,:,:], axis=0)
	#~ std_Di    = np.std(var[:][1345:1435,:,:], axis=0) 
	#~ std_Dii   = np.std(var[:][1350:1440,:,:], axis=0) 

	#~ mean_D_ii  = np.nanmean(var[:][1330:1420,:,:], axis=0) 
	#~ mean_D_i   = np.nanmean(var[:][1335:1425,:,:], axis=0)
	#~ mean_D     = np.nanmean(var[:][1340:1430,:,:], axis=0)
	#~ mean_Di    = np.nanmean(var[:][1345:1435,:,:], axis=0) 
	#~ mean_Dii   = np.nanmean(var[:][1350:1440,:,:], axis=0) 

	# EEOC
	#~ D_ii  = np.nanmean(var[:][1701:1785,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ D_i   = np.nanmean(var[:][1703:1790,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ D     = np.nanmean(var[:][1705:1795,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ Di    = np.nanmean(var[:][1707:1800,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ Dii   = np.nanmean(var[:][1709:1805,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)

	#~ std_D_ii  = np.std(var[:][1701:1785,:,:], axis=0) 
	#~ std_D_i   = np.std(var[:][1703:1790,:,:], axis=0)
	#~ std_D     = np.std(var[:][1705:1795,:,:], axis=0)
	#~ std_Di    = np.std(var[:][1707:1800,:,:], axis=0) 
	#~ std_Dii   = np.std(var[:][1709:1805,:,:], axis=0) 

	#~ mean_D_ii  = np.nanmean(var[:][1701:1785,:,:], axis=0) 
	#~ mean_D_i   = np.nanmean(var[:][1703:1790,:,:], axis=0)
	#~ mean_D     = np.nanmean(var[:][1705:1795,:,:], axis=0)
	#~ mean_Di    = np.nanmean(var[:][1707:1800,:,:], axis=0) 
	#~ mean_Dii   = np.nanmean(var[:][1709:1805,:,:], axis=0) 
		
	return lat, lon, std_D_ii, std_D_i, std_D, std_Di, std_Dii, std_clim, mean_D_ii, mean_D_i, mean_D, mean_Di, mean_Dii, mean_clim, D_ii, D_i, D, Di, Dii

	
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
	
	
# Import era5 database
print('Import era5 database')
event    = u'ENOC'
season   = u'DJF'
variable = u'prec'

dict_unit = {u'prec': u'PRE (mm d⁻¹)',
u'mx2t': u'Tmax (°C)',
u'mtnlwrf': u'ROL (W m⁻²)',
u'msl': u'PNMM (hPa)',
u'q': u'Q (g kg⁻¹)',
u'z': u'GEO 500hPa (m)',
u'uv10': u'UV10m (m s⁻¹)'}

lat, lon, std_D_ii, std_D_i, std_D, std_Di, std_Dii, std_clim, mean_D_ii, mean_D_i, mean_D, mean_Di, mean_Dii, mean_clim, D_ii, D_i, D, Di, Dii = import_era5(variable) 

# Plot maps with the function
print('Plot maps with the function')
fig = plt.figure(figsize=(8, 2))

# PRE
cor_map = mpl.cm.BrBG
levs1  = [-8, -6, -4, -2, 2, 4, 6, 8]
# TEMP
#~ cor_map = mpl.cm.bwr
#~ levs1   = [-4, -3, -2, -1, 1, 2, 3, 4]
# ROL
#~ cor_map = mpl.cm.bwr
#~ levs1   = [-40, -30, -20, -10, 10, 20, 30, 40]
# PNMM
#~ cmap    = mpl.cm.PiYG
#~ cor_map = cmap.reversed()
#~ levs1 = [-4, -3, -2, -1, 1, 2, 3, 4]
# Q
#~ cmap = mpl.cm.bwr
#~ cor_map = cmap.reversed()
#~ levs1 = [-4, -3, -2, -1, 1, 2, 3, 4]
# GEO
#~ cor_map = mpl.cm.PiYG
#~ levs1 = [-40, -30, -20, -10, 10, 20, 30, 40]
# UV10m
#~ cor_map = mpl.cm.RdGy
#~ levs1 = [-2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2]

ax = fig.add_subplot(1, 5, 1)
map, xx, yy = basemap(lat, lon)
plt.title(u'U) {0} (D-2) \n {1} {2}'.format(event, season, dict_unit[variable]), fontsize=6, fontweight='bold')
map.contourf(xx, yy, D_ii, levels=levs1, latlon=True, cmap=cor_map)	
p_value = ttest(std_D_ii, std_clim, mean_D_ii, mean_clim)
p_value = ma.masked_where(p_value >= 0.05, p_value) 
map.contourf(xx, yy, p_value, colors='none', hatches=['....'])

ax = fig.add_subplot(1, 5, 2)
map, xx, yy = basemap(lat, lon)
plt.title(u'V) {0} (D-1) \n {1} {2}'.format(event, season, dict_unit[variable]), fontsize=6, fontweight='bold')
map.contourf(xx, yy, D_i, levels=levs1, latlon=True, cmap=cor_map)
p_value = ttest(std_D_i, std_clim, mean_D_i, mean_clim)
p_value = ma.masked_where(p_value >= 0.05, p_value) 
map.contourf(xx, yy, p_value, colors='none', hatches=['....'])

ax = fig.add_subplot(1, 5, 3)
map, xx, yy = basemap(lat, lon)
plt.title(u'W) {0} (D0) \n {1} {2}'.format(event, season, dict_unit[variable]), fontsize=6, fontweight='bold')
map.contourf(xx, yy, D, levels=levs1, latlon=True, cmap=cor_map)
p_value = ttest(std_D, std_clim, mean_D, mean_clim)
p_value = ma.masked_where(p_value >= 0.05, p_value) 
map.contourf(xx, yy, p_value, colors='none', hatches=['....'])

ax = fig.add_subplot(1, 5, 4)
map, xx, yy = basemap(lat, lon)
plt.title(u'X) {0} (D+1) \n {1} {2}'.format(event, season, dict_unit[variable]), fontsize=6, fontweight='bold')
map.contourf(xx, yy, Di, levels=levs1, latlon=True, cmap=cor_map)
p_value = ttest(std_Di, std_clim, mean_Di, mean_clim)
p_value = ma.masked_where(p_value >= 0.05, p_value) 
map.contourf(xx, yy, p_value, colors='none', hatches=['....'])

ax = fig.add_subplot(1, 5, 5)
map, xx, yy = basemap(lat, lon)
plt.title(u'Y) {0} (D+2) \n {1} {2}'.format(event, season, dict_unit[variable]), fontsize=6, fontweight='bold')
map.contourf(xx, yy, Dii, levels=levs1, latlon=True, cmap=cor_map)  
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6)
p_value = ttest(std_Dii, std_clim, mean_Dii, mean_clim)
p_value = ma.masked_where(p_value >= 0.05, p_value) 
map.contourf(xx, yy, p_value, colors='none', hatches=['....'])

# Path out to save figure
print('Path out to save figure')
path_out = '/home/nice/Documents/janio/fig'
name_out = 'pyplt_maps_composites_{0}_{1}_{2}_era5.png'.format(event, season, variable)
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')
plt.show()
exit()	
