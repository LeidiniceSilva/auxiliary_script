# -*- coding: utf-8 -*-

__author__ = "Leidinice Silva"
__email__ = "leidinice.silvae@funceme.br"
__date__ = "08/06/2017"
__description__ = " Plot map precipitation of olamv.3.3 model and obs database"

import os
import conda
import netCDF4
import numpy as np
import pandas as pd
import matplotlib.cm as cm
import matplotlib.pyplot as plt

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
    

def import_sim(path, exp):

	arq  = '{0}/precip_controle_1982_2012_{1}_g2_neb_new_REAL_ok_full_negcor_yearsum_noocean.nc'.format(path, exp)
	data = netCDF4.Dataset(arq)
	var  = data.variables['precip'][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	ymean_exp = var[:][:,:,:]
	
	return lat, lon, ymean_exp


def import_obs(path):

	arq  = '{0}/pr_Amon_CRU-TS3.22_observation_198201-201212_new_mma_neb.nc'.format(path)
	data = netCDF4.Dataset(arq)
	var  = data.variables['pr'][:]
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	ymean_obs = var[:][:,:,:]
	
	return lat, lon, ymean_obs
	

# Import exp model end obs database 
home = os.path.expanduser("~")
path = home + "/Documents/ufrn/papers/olam/datas"

exp1  = u'chen'
lat, lon, ymean_exp1 = import_sim(path, exp1)

exp2  = u'harr'
lat, lon, ymean_exp2 = import_sim(path, exp2)

lat, lon, ymean_obs1 = import_obs(path)


# Plot maps olam model cru database
fig = plt.figure(1, figsize=(34,24))

plt.subplot(131)
plt.title(u'OLAMv.3.3_g2_Chen', fontsize=20)
plt.ylabel(u'Latitude', fontsize=20, labelpad=40)
map = Basemap(projection='cyl', llcrnrlon=-43.25, llcrnrlat=-9.75, urcrnrlon=-33.25, urcrnrlat=-2.25, resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawstates()
xx, yy = map(lon,lat)
levs   = [50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
plot_maps = plt.contourf(xx, yy, ymean_exp1[0,:,:], levels=levs, cmap=cm.jet, extend='both')

plt.subplot(132)
plt.title(u'OLAMv.3.3_g2_Harr', fontsize=20)
plt.xlabel(u'Longitude', fontsize=20, labelpad=30)
map = Basemap(projection='cyl', llcrnrlon=-43.25, llcrnrlat=-9.75, urcrnrlon=-33.25, urcrnrlat=-2.25, resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawstates()
xx, yy = map(lon,lat)
levs   = [50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
plot_maps = plt.contourf(xx, yy, ymean_exp2[0,:,:], levels=levs, cmap=cm.jet, extend='both')

plt.subplot(133)
plt.title(u'CRU', fontsize=20)
map = Basemap(projection='cyl', llcrnrlon=-43.25, llcrnrlat=-9.75, urcrnrlon=-33.25, urcrnrlat=-2.25, resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawstates()
xx, yy = map(lon,lat)
levs   = [50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
plot_maps = plt.contourf(xx, yy, ymean_obs1[0,:,:], levels=levs, cmap=cm.jet, extend='both')
bar = fig.colorbar(plot_maps, spacing='uniform', ticks=levs, drawedges=True)

path_out = home + "/Documents/ufrn/papers/olam/results/"
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, 'mapa_anual_chen_harr_cru.png'), bbox_inches='tight', dpi=300)
exit()






