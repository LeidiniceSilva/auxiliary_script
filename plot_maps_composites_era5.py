# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "05/30/2021"
__description__ = "This script plot compoosites maps from era5"

import os
import conda
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import warnings ; warnings.filterwarnings("ignore")
import matplotlib.cm as cm

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from matplotlib.path import Path
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import PathPatch


def import_era5(var):

	dict_var = {u'tmax': u'mx2t', 
	u'tmin': u'mn2t',
	u'pre': u'tp', 
	u'rad': u'mtnlwrf'}
	
	path = '/home/nice/Downloads'
	arq  = '{0}/{1}_era5_br_day_2011-2020.nc'.format(path, var)	
	data = netCDF4.Dataset(arq)
	var  = data.variables[dict_var[var]][:]
	lat  = data.variables['latitude'][:]
	lon  = data.variables['longitude'][:]
	
	D  = np.nanmean(var[:][2495:2616,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	DO = np.nanmean(var[:][2525:2646,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	DI = np.nanmean(var[:][2565:2676,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
		
	return lat, lon, D, DO, DI
	
	
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
	map.drawmapboundary(color='white')
 
	lons, lats = np.meshgrid(new_lon, new_lat)

	xx, yy = map(lons,lats)

	# Import shapefile from word and matopiba 
	map.readshapefile('/home/nice/Documents/github_projects/shp/shp_world/world', 'world', drawbounds=True, color='white')
	map.readshapefile('/home/nice/Documents/github_projects/shp/lim_unid_fed/lim_unid_fed', 'lim_unid_fed', drawbounds=True, color='black')
	x0, x1 = plt.xlim()
	y0, y1 = plt.ylim()
	map_edges = np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]])
	polys = [map_edges]
	map.readshapefile('/home/nice/Documents/github_projects/shp/lim_unid_fed/lim_unid_fed', 'lim_unid_fed2', drawbounds=False)
	polys = polys + getattr(map, 'lim_unid_fed2')
	codes = [[Path.MOVETO] + [Path.LINETO for p in p[1:]] for p in polys] # creating a PathPatch
	polys_lin = [v for p in polys for v in p]
	codes_lin = [cdg for cdgs in codes for cdg in cdgs]
	path  = Path(polys_lin, codes_lin)
	patch = PathPatch(path, facecolor='white', lw=0)
	plt.gca().add_patch(patch)
	
	return map, xx, yy
	
	
# Import regcm exp and cru databases 	
lat, lon, D_mx2t, DO_mx2t, DI_mx2t = import_era5('tmax')   
lat, lon, D_mn2t, DO_mn2t, DI_mn2t = import_era5('tmin')   
lat, lon, D_pre, DO_pre, DI_pre = import_era5('pre')   
lat, lon, D_rad, DO_rad, DI_rad = import_era5('rad')   

# Plot maps with the function
fig = plt.figure(figsize=(5, 7))
levs1 = [-3, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 3]
levs2 = [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5]
levs3 = [-40, -30, -20, -10, -5, 0, 5, 10, 20, 30, 40]

ax = fig.add_subplot(4, 3, 1)
map, xx, yy = basemap(lat, lon)
plt.title(u'A) Evento de OC (D-1) \n Tmax (°C)', fontsize=6, fontweight='bold')
map.contourf(xx, yy, D_mx2t, levels=levs1, latlon=True, cmap=cm.bwr)
	
ax = fig.add_subplot(4, 3, 2)
map, xx, yy = basemap(lat, lon)
plt.title(u'B) Evento de OC (D0) \n Tmax (°C)', fontsize=6, fontweight='bold')
map.contourf(xx, yy, DO_mx2t, levels=levs1, latlon=True, cmap=cm.bwr)

ax = fig.add_subplot(4, 3, 3)
map, xx, yy = basemap(lat, lon)
plt.title(u'C) Evento de OC (D+1) \n Tmax (°C)', fontsize=6, fontweight='bold')
map.contourf(xx, yy, DI_mx2t, levels=levs1, latlon=True, cmap=cm.bwr) 
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(4, 3, 4)
map, xx, yy = basemap(lat, lon)
plt.title(u'D) Evento de OC (D-1) \n Tmin (°C)', fontsize=6, fontweight='bold')
map.contourf(xx, yy, D_mn2t, levels=levs1, latlon=True, cmap=cm.bwr)

ax = fig.add_subplot(4, 3, 5)
map, xx, yy = basemap(lat, lon)
plt.title(u'E) Evento de OC (D0) \n Tmin (°C)', fontsize=6, fontweight='bold')
map.contourf(xx, yy, DO_mn2t, levels=levs1, latlon=True, cmap=cm.bwr) 

ax = fig.add_subplot(4, 3, 6)
map, xx, yy = basemap(lat, lon)
plt.title(u'F) Evento de OC (D+1) \n Tmin (°C)', fontsize=6, fontweight='bold')
map.contourf(xx, yy, DI_mn2t, levels=levs1, latlon=True, cmap=cm.bwr)
cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(4, 3, 7)
map, xx, yy = basemap(lat, lon)
plt.title(u'G) Evento de OC (D-1) \n Pre (mm)', fontsize=6, fontweight='bold')
map.contourf(xx, yy, D_pre, levels=levs2, latlon=True, cmap=cm.BrBG) 

ax = fig.add_subplot(4, 3, 8)
map, xx, yy = basemap(lat, lon)
plt.title(u'H) Evento de OC (D0) \n Pre (mm)', fontsize=6, fontweight='bold')
map.contourf(xx, yy, DO_pre, levels=levs2, latlon=True, cmap=cm.BrBG)

ax = fig.add_subplot(4, 3, 9)
map, xx, yy = basemap(lat, lon)
plt.title(u'I) Evento de OC (D+1) \n Pre (mm)', fontsize=6, fontweight='bold')
map.contourf(xx, yy, DI_pre, levels=levs2, latlon=True, cmap=cm.BrBG)
cbar = map.colorbar(ticks=levs2, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

ax = fig.add_subplot(4, 3, 10) 
map, xx, yy = basemap(lat, lon)
plt.title(u'J) Evento de OC (D-1) \n Rad (W m⁻²)', fontsize=6, fontweight='bold')
map.contourf(xx, yy, D_rad, levels=levs3, latlon=True, cmap=cm.bwr)

ax = fig.add_subplot(4, 3, 11)
map, xx, yy = basemap(lat, lon)
plt.title(u'K) Evento de OC (D0) \n Rad (W m⁻²)', fontsize=6, fontweight='bold')
map.contourf(xx, yy, DO_rad, levels=levs3, latlon=True, cmap=cm.bwr) 

ax = fig.add_subplot(4, 3, 12)
map, xx, yy = basemap(lat, lon)
plt.title(u'L) Evento de OC (D+1) \n Rad (W m⁻²)', fontsize=6, fontweight='bold')
map.contourf(xx, yy, DI_rad, levels=levs3, latlon=True, cmap=cm.bwr) 
cbar = map.colorbar(ticks=levs3, drawedges=True, ax=ax)
cbar.ax.tick_params(labelsize=6) 

# Path out to save bias figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_composites_oc_era5.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=100, bbox_inches='tight')

plt.show()
exit()	
	
	
	
	
	
