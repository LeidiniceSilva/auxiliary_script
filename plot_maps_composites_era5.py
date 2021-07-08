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
	u'uv10': u'u10',
	u'u10': u'u10',
	u'v10': u'v10'}

	print('Open data')
	
	path = '/home/nice/Documents/era5'
	arq  = '{0}/{1}_xavier_br_day_2011-2016_lonlat.nc'.format(path, var)	
	data = netCDF4.Dataset(arq)		
	var  = data.variables[dict_var[var]][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	dataset = var[:][:,:,:]	

	#~ # OC
	#~ D_ii  = np.nanmean(var[:][1401:1520,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ D_i   = np.nanmean(var[:][1416:1535,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ D     = np.nanmean(var[:][1431:1551,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ Di    = np.nanmean(var[:][1446:1566,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ Dii   = np.nanmean(var[:][1461:1580,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	
	# EOC
	D_ii  = np.nanmean(var[:][1766:1886,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	D_i   = np.nanmean(var[:][1781:1901,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	D     = np.nanmean(var[:][1795:1926,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	Di    = np.nanmean(var[:][1810:1930,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	Dii   = np.nanmean(var[:][1825:1945,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)

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
print('Import database: PRE')
lat, lon, var, D_ii, D_i, D, Di, Dii = import_era5('prec')   
levs1 = [-4, -3, -2, -1, 1, 2, 3, 4]

ax = fig.add_subplot(1, 5, 1)
map, xx, yy = basemap(lat, lon)
plt.title(u'A) Evento de EOC (D-2) \n PRE (mm d⁻¹)', fontsize=6, fontweight='bold')
map.contourf(xx, yy, D_ii, levels=levs1, latlon=True, cmap=cm.BrBG)	
p_value = ttest(var, D_ii)
p_value = ma.masked_where(p_value >= 0.05, p_value) 
map.contourf(xx, yy, p_value, colors='none', hatches=['....'])

ax = fig.add_subplot(1, 5, 2)
map, xx, yy = basemap(lat, lon)
plt.title(u'B) Evento de EOC (D-1) \n PRE (mm d⁻¹)', fontsize=6, fontweight='bold')
map.contourf(xx, yy, D_i, levels=levs1, latlon=True, cmap=cm.BrBG)
p_value = ttest(var, D_i)
p_value = ma.masked_where(p_value >= 0.05, p_value) 
map.contourf(xx, yy, p_value, colors='none', hatches=['....'])

ax = fig.add_subplot(1, 5, 3)
map, xx, yy = basemap(lat, lon)
plt.title(u'C) Evento de EOC (D0) \n PRE (mm d⁻¹)', fontsize=6, fontweight='bold')
map.contourf(xx, yy, D, levels=levs1, latlon=True, cmap=cm.BrBG)
p_value = ttest(var, D)
p_value = ma.masked_where(p_value >= 0.05, p_value) 
map.contourf(xx, yy, p_value, colors='none', hatches=['....'])

ax = fig.add_subplot(1, 5, 4)
map, xx, yy = basemap(lat, lon)
plt.title(u'D) Evento de EOC (D+1) \n PRE (mm d⁻¹)', fontsize=6, fontweight='bold')
map.contourf(xx, yy, Di, levels=levs1, latlon=True, cmap=cm.BrBG)
p_value = ttest(var, Di)
p_value = ma.masked_where(p_value >= 0.05, p_value) 
map.contourf(xx, yy, p_value, colors='none', hatches=['....'])

ax = fig.add_subplot(1, 5, 5)
map, xx, yy = basemap(lat, lon)
plt.title(u'E) Evento de EOC (D+2) \n PRE (mm d⁻¹)', fontsize=6, fontweight='bold')
map.contourf(xx, yy, Dii, levels=levs1, latlon=True, cmap=cm.BrBG)  
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6)
p_value = ttest(var, Dii)
p_value = ma.masked_where(p_value >= 0.05, p_value) 
map.contourf(xx, yy, p_value, colors='none', hatches=['....'])

# Path out to save figure
print('Path out to save figure')
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_composites_pre_eoc_era5.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')
plt.show()
exit()	

#~ # Import era5 database
#~ print('Import database: Tmax')
#~ lat, lon, var, D_ii, D_i, D, Di, Dii = import_era5('mx2t')   
#~ levs1 = [-4, -3, -2, -1, 1, 2, 3, 4]

#~ ax = fig.add_subplot(1, 5, 1)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'F) Evento de EOC (D-2) \n Tmax (°C)', fontsize=6, fontweight='bold')
#~ map.contourf(xx, yy, D_ii, levels=levs1, latlon=True, cmap=cm.bwr)	
#~ p_value = ttest(var, D_ii)
#~ p_value = ma.masked_where(p_value >= 0.05, p_value) 
#~ map.contourf(xx, yy, p_value, colors='none', hatches=['.....'])

#~ ax = fig.add_subplot(1, 5, 2)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'G) Evento de EOC (D-1) \n Tmax (°C)', fontsize=6, fontweight='bold')
#~ map.contourf(xx, yy, D_i, levels=levs1, latlon=True, cmap=cm.bwr)
#~ p_value = ttest(var, D_i)
#~ p_value = ma.masked_where(p_value >= 0.05, p_value) 
#~ map.contourf(xx, yy, p_value, colors='none', hatches=['.....'])

#~ ax = fig.add_subplot(1, 5, 3)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'H) Evento de EOC (D0) \n Tmax (°C)', fontsize=6, fontweight='bold')
#~ map.contourf(xx, yy, D, levels=levs1, latlon=True, cmap=cm.bwr)
#~ p_value = ttest(var, D)
#~ p_value = ma.masked_where(p_value >= 0.05, p_value) 
#~ map.contourf(xx, yy, p_value, colors='none', hatches=['.....'])

#~ ax = fig.add_subplot(1, 5, 4)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'I) Evento de EOC (D+1) \n Tmax (°C)', fontsize=6, fontweight='bold')
#~ map.contourf(xx, yy, Di, levels=levs1, latlon=True, cmap=cm.bwr)
#~ p_value = ttest(var, Di)
#~ p_value = ma.masked_where(p_value >= 0.05, p_value) 
#~ map.contourf(xx, yy, p_value, colors='none', hatches=['.....'])

#~ ax = fig.add_subplot(1, 5, 5)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'J) Evento de EOC (D+2) \n Tmax (°C)', fontsize=6, fontweight='bold')
#~ map.contourf(xx, yy, Dii, levels=levs1, latlon=True, cmap=cm.bwr)  
#~ cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
#~ cbar.ax.tick_params(labelsize=6)
#~ p_value = ttest(var, Dii)
#~ p_value = ma.masked_where(p_value >= 0.05, p_value) 
#~ map.contourf(xx, yy, p_value, colors='none', hatches=['.....'])

#~ # Path out to save figure
#~ print('Path out to save figure')
#~ path_out = '/home/nice/Downloads'
#~ name_out = 'pyplt_maps_composites_tmax_eoc_era5.png'
#~ if not os.path.exists(path_out):
	#~ create_path(path_out)
#~ plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')
#~ plt.show()
#~ exit()	

#~ # Import era5 database
#~ print('Import database: ROL')
#~ lat, lon, var, D_ii, D_i, D, Di, Dii = import_era5('mtnlwrf')   
#~ levs1 = [-40, -30, -20, -10, 10, 20, 30, 40]

#~ ax = fig.add_subplot(1, 5, 1)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'K) Evento de EOC (D-2) \n ROL (W m⁻²)', fontsize=6, fontweight='bold')
#~ map.contourf(xx, yy, D_ii, levels=levs1, latlon=True, cmap=cm.bwr)	
#~ p_value = ttest(var, D_ii)
#~ p_value = ma.masked_where(p_value >= 0.05, p_value) 
#~ map.contourf(xx, yy, p_value, colors='none', hatches=['.....'])

#~ ax = fig.add_subplot(1, 5, 2)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'L) Evento de EOC (D-1) \n ROL (W m⁻²)', fontsize=6, fontweight='bold')
#~ map.contourf(xx, yy, D_i, levels=levs1, latlon=True, cmap=cm.bwr)
#~ p_value = ttest(var, D_i)
#~ p_value = ma.masked_where(p_value >= 0.05, p_value) 
#~ map.contourf(xx, yy, p_value, colors='none', hatches=['.....'])

#~ ax = fig.add_subplot(1, 5, 3)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'M) Evento de EOC (D0) \n ROL (W m⁻²)', fontsize=6, fontweight='bold')
#~ map.contourf(xx, yy, D, levels=levs1, latlon=True, cmap=cm.bwr)
#~ p_value = ttest(var, D)
#~ p_value = ma.masked_where(p_value >= 0.05, p_value) 
#~ map.contourf(xx, yy, p_value, colors='none', hatches=['.....'])

#~ ax = fig.add_subplot(1, 5, 4)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'N) Evento de EOC (D+1) \n ROL (W m⁻²)', fontsize=6, fontweight='bold')
#~ map.contourf(xx, yy, Di, levels=levs1, latlon=True, cmap=cm.bwr)
#~ p_value = ttest(var, Di)
#~ p_value = ma.masked_where(p_value >= 0.05, p_value) 
#~ map.contourf(xx, yy, p_value, colors='none', hatches=['.....'])

#~ ax = fig.add_subplot(1, 5, 5)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'O) Evento de EOC (D+2) \n ROL (W m⁻²)', fontsize=6, fontweight='bold')
#~ map.contourf(xx, yy, Dii, levels=levs1, latlon=True, cmap=cm.bwr)  
#~ cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
#~ cbar.ax.tick_params(labelsize=6)
#~ p_value = ttest(var, Dii)
#~ p_value = ma.masked_where(p_value >= 0.05, p_value) 
#~ map.contourf(xx, yy, p_value, colors='none', hatches=['.....'])

#~ # Path out to save figure
#~ print('Path out to save figure')
#~ path_out = '/home/nice/Downloads'
#~ name_out = 'pyplt_maps_composites_rol_eoc_era5.png'
#~ if not os.path.exists(path_out):
	#~ create_path(path_out)
#~ plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')
#~ plt.show()
#~ exit()	

#~ # Import era5 database
#~ print('Import database: PNMM')
#~ lat, lon, var, D_ii, D_i, D, Di, Dii = import_era5('msl')   
#~ levs1 = [-4, -3, -2, -1, 1, 2, 3, 4]
#~ cmap = mpl.cm.PiYG
#~ cmap_r = cmap.reversed()

#~ ax = fig.add_subplot(1, 5, 1)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'P) Evento de EOC (D-2) \n PNMM (hPa)', fontsize=6, fontweight='bold')
#~ map.contourf(xx, yy, D_ii, levels=levs1, latlon=True, cmap=cmap_r)	
#~ p_value = ttest(var, D_ii)
#~ p_value = ma.masked_where(p_value >= 0.05, p_value) 
#~ map.contourf(xx, yy, p_value, colors='none', hatches=['.....'])

#~ ax = fig.add_subplot(1, 5, 2)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'Q) Evento de EOC (D-1) \n PNMM (hPa)', fontsize=6, fontweight='bold')
#~ map.contourf(xx, yy, D_i, levels=levs1, latlon=True, cmap=cmap_r)
#~ p_value = ttest(var, D_i)
#~ p_value = ma.masked_where(p_value >= 0.05, p_value) 
#~ map.contourf(xx, yy, p_value, colors='none', hatches=['.....'])

#~ ax = fig.add_subplot(1, 5, 3)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'R) Evento de EOC (D0) \n PNMM (hPa)', fontsize=6, fontweight='bold')
#~ map.contourf(xx, yy, D, levels=levs1, latlon=True, cmap=cmap_r)
#~ p_value = ttest(var, D)
#~ p_value = ma.masked_where(p_value >= 0.05, p_value) 
#~ map.contourf(xx, yy, p_value, colors='none', hatches=['.....'])

#~ ax = fig.add_subplot(1, 5, 4)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'S) Evento de EOC (D+1) \n PNMM (hPa)', fontsize=6, fontweight='bold')
#~ map.contourf(xx, yy, Di, levels=levs1, latlon=True, cmap=cmap_r)
#~ p_value = ttest(var, Di)
#~ p_value = ma.masked_where(p_value >= 0.05, p_value) 
#~ map.contourf(xx, yy, p_value, colors='none', hatches=['.....'])

#~ ax = fig.add_subplot(1, 5, 5)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'T) Evento de EOC (D+2) \n PNMM (hPa)', fontsize=6, fontweight='bold')
#~ map.contourf(xx, yy, Dii, levels=levs1, latlon=True, cmap=cmap_r)  
#~ cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
#~ cbar.ax.tick_params(labelsize=6)
#~ p_value = ttest(var, Dii)
#~ p_value = ma.masked_where(p_value >= 0.05, p_value) 
#~ map.contourf(xx, yy, p_value, colors='none', hatches=['.....'])

#~ # Path out to save figure
#~ print('Path out to save figure')
#~ path_out = '/home/nice/Downloads'
#~ name_out = 'pyplt_maps_composites_pnmm_eoc_era5.png'
#~ if not os.path.exists(path_out):
	#~ create_path(path_out)
#~ plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')
#~ plt.show()
#~ exit()	

#~ # Import era5 database
#~ print('Import database: UR')
#~ lat, lon, var, D_ii, D_i, D, Di, Dii = import_era5('r')   
#~ levs1 = [-8, -6, -4, -2, 2, 4, 6, 8]
#~ cmap = mpl.cm.bwr
#~ cmap_r = cmap.reversed()

#~ ax = fig.add_subplot(1, 5, 1)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'U) Evento de EOC (D-2) \n UR (%)', fontsize=6, fontweight='bold')
#~ map.contourf(xx, yy, D_ii, levels=levs1, latlon=True, cmap=cmap_r)	
#~ p_value = ttest(var, D_ii)
#~ p_value = ma.masked_where(p_value >= 0.05, p_value) 
#~ map.contourf(xx, yy, p_value, colors='none', hatches=['.....'])

#~ ax = fig.add_subplot(1, 5, 2)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'V) Evento de EOC (D-1) \n UR (%)', fontsize=6, fontweight='bold')
#~ map.contourf(xx, yy, D_i, levels=levs1, latlon=True, cmap=cmap_r)
#~ p_value = ttest(var, D_i)
#~ p_value = ma.masked_where(p_value >= 0.05, p_value) 
#~ map.contourf(xx, yy, p_value, colors='none', hatches=['.....'])

#~ ax = fig.add_subplot(1, 5, 3)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'W) Evento de EOC (D0) \n UR (%)', fontsize=6, fontweight='bold')
#~ map.contourf(xx, yy, D, levels=levs1, latlon=True, cmap=cmap_r)
#~ p_value = ttest(var, D)
#~ p_value = ma.masked_where(p_value >= 0.05, p_value) 
#~ map.contourf(xx, yy, p_value, colors='none', hatches=['.....'])

#~ ax = fig.add_subplot(1, 5, 4)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'X) Evento de EOC (D+1) \n UR (%)', fontsize=6, fontweight='bold')
#~ map.contourf(xx, yy, Di, levels=levs1, latlon=True, cmap=cmap_r)
#~ p_value = ttest(var, Di)
#~ p_value = ma.masked_where(p_value >= 0.05, p_value) 
#~ map.contourf(xx, yy, p_value, colors='none', hatches=['.....'])

#~ ax = fig.add_subplot(1, 5, 5)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'Y) Evento de EOC (D+2) \n UR (%)', fontsize=6, fontweight='bold')
#~ map.contourf(xx, yy, Dii, levels=levs1, latlon=True, cmap=cmap_r)  
#~ cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
#~ cbar.ax.tick_params(labelsize=6)
#~ p_value = ttest(var, Dii)
#~ p_value = ma.masked_where(p_value >= 0.05, p_value) 
#~ map.contourf(xx, yy, p_value, colors='none', hatches=['.....'])

#~ # Path out to save figure
#~ print('Path out to save figure')
#~ path_out = '/home/nice/Downloads'
#~ name_out = 'pyplt_maps_composites_ur_eoc_era5.png'
#~ if not os.path.exists(path_out):
	#~ create_path(path_out)
#~ plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')
#~ plt.show()
#~ exit()	

#~ # Import era5 database
#~ print('Import database: GEO')
#~ lat, lon, var, D_ii, D_i, D, Di, Dii = import_era5('z')   
#~ levs1 = [-40, -30, -20, -10, 10, 20, 30, 40]

#~ ax = fig.add_subplot(1, 5, 1)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'Z) Evento de EOC (D-2) \n GEO 500hPa (m)', fontsize=6, fontweight='bold')
#~ map.contourf(xx, yy, D_ii, levels=levs1, latlon=True, cmap=cm.PiYG)	
#~ p_value = ttest(var, D_ii)
#~ p_value = ma.masked_where(p_value >= 0.05, p_value) 
#~ map.contourf(xx, yy, p_value, colors='none', hatches=['.....'])

#~ ax = fig.add_subplot(1, 5, 2)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'A.1) Evento de EOC (D-1) \n GEO 500hPa (m)', fontsize=6, fontweight='bold')
#~ map.contourf(xx, yy, D_i, levels=levs1, latlon=True, cmap=cm.PiYG)
#~ p_value = ttest(var, D_i)
#~ p_value = ma.masked_where(p_value >= 0.05, p_value) 
#~ map.contourf(xx, yy, p_value, colors='none', hatches=['.....'])

#~ ax = fig.add_subplot(1, 5, 3)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'B.1) Evento de EOC (D0) \n GEO 500hPa (m)', fontsize=6, fontweight='bold')
#~ map.contourf(xx, yy, D, levels=levs1, latlon=True, cmap=cm.PiYG)
#~ p_value = ttest(var, D)
#~ p_value = ma.masked_where(p_value >= 0.05, p_value) 
#~ map.contourf(xx, yy, p_value, colors='none', hatches=['.....'])

#~ ax = fig.add_subplot(1, 5, 4)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'C.1) Evento de EOC (D+1) \n GEO 500hPa (m)', fontsize=6, fontweight='bold')
#~ map.contourf(xx, yy, Di, levels=levs1, latlon=True, cmap=cm.PiYG)
#~ p_value = ttest(var, Di)
#~ p_value = ma.masked_where(p_value >= 0.05, p_value) 
#~ map.contourf(xx, yy, p_value, colors='none', hatches=['.....'])

#~ ax = fig.add_subplot(1, 5, 5)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'D.1) Evento de EOC (D+2) \n GEO 500 hPa (m)', fontsize=6, fontweight='bold')
#~ map.contourf(xx, yy, Dii, levels=levs1, latlon=True, cmap=cm.PiYG)  
#~ cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
#~ cbar.ax.tick_params(labelsize=6)
#~ p_value = ttest(var, Dii)
#~ p_value = ma.masked_where(p_value >= 0.05, p_value) 
#~ map.contourf(xx, yy, p_value, colors='none', hatches=['.....'])

#~ # Path out to save figure
#~ print('Path out to save figure')
#~ path_out = '/home/nice/Downloads'
#~ name_out = 'pyplt_maps_composites_z_eoc_era5.png'
#~ if not os.path.exists(path_out):
	#~ create_path(path_out)
#~ plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')
#~ plt.show()
#~ exit()	

#~ # Import era5 database
#~ print('Import database: UV10m')
#~ lat, lon, var, D_ii, D_i, D, Di, Dii = import_era5('uv10')    
#~ levs1 = [-2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2]

#~ ax = fig.add_subplot(1, 5, 1)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'E.1) Evento de EOC (D-2) \n UV10m (m s⁻¹)', fontsize=6, fontweight='bold')
#~ map.contourf(xx, yy, D_ii, levels=levs1, latlon=True, cmap=cm.RdGy)	
#~ p_value = ttest(var, D_ii)
#~ p_value = ma.masked_where(p_value >= 0.05, p_value) 
#~ map.contourf(xx, yy, p_value, colors='none', hatches=['.....'])

#~ ax = fig.add_subplot(1, 5, 2)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'F.1) Evento de EOC (D-1) \n UV10m (m s⁻¹)', fontsize=6, fontweight='bold')
#~ map.contourf(xx, yy, D_i, levels=levs1, latlon=True, cmap=cm.RdGy)
#~ p_value = ttest(var, D_i)
#~ p_value = ma.masked_where(p_value >= 0.05, p_value) 
#~ map.contourf(xx, yy, p_value, colors='none', hatches=['.....'])

#~ ax = fig.add_subplot(1, 5, 3)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'G.1) Evento de EOC (D0) \n UV10m (m s⁻¹)', fontsize=6, fontweight='bold')
#~ map.contourf(xx, yy, D, levels=levs1, latlon=True, cmap=cm.RdGy)
#~ p_value = ttest(var, D)
#~ p_value = ma.masked_where(p_value >= 0.05, p_value) 
#~ map.contourf(xx, yy, p_value, colors='none', hatches=['.....'])

#~ ax = fig.add_subplot(1, 5, 4)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'H.1) Evento de EOC (D+1) \n UV10m (m s⁻¹)', fontsize=6, fontweight='bold')
#~ map.contourf(xx, yy, Di, levels=levs1, latlon=True, cmap=cm.RdGy)
#~ p_value = ttest(var, Di)
#~ p_value = ma.masked_where(p_value >= 0.05, p_value) 
#~ map.contourf(xx, yy, p_value, colors='none', hatches=['.....'])

#~ ax = fig.add_subplot(1, 5, 5)
#~ map, xx, yy = basemap(lat, lon)
#~ plt.title(u'I.1) Evento de EOC (D+2) \n UV10m (m s⁻¹)', fontsize=6, fontweight='bold')
#~ map.contourf(xx, yy, Dii, levels=levs1, latlon=True, cmap=cm.RdGy)  
#~ cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
#~ cbar.ax.tick_params(labelsize=6)
#~ p_value = ttest(var, Dii)
#~ p_value = ma.masked_where(p_value >= 0.05, p_value) 
#~ map.contourf(xx, yy, p_value, colors='none', hatches=['.....'])

#~ # Path out to save figure
#~ print('Path out to save figure')
#~ path_out = '/home/nice/Downloads'
#~ name_out = 'pyplt_maps_composites_uv10m_eoc_era5.png'
#~ if not os.path.exists(path_out):
	#~ create_path(path_out)
#~ plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')
#~ plt.show()
#~ exit()	

