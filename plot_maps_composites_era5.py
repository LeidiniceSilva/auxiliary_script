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

	dict_var = {u'mtpr': u'mtpr',
	u'mx2t': u'mx2t',
	u'mn2t': u'mn2t',
	u'mtnlwrf': u'mtnlwrf',
	u'q': u'q'}

	path = '/home/nice/Downloads/janio/dados'
	arq  = '{0}/{1}_era5_br_day_2011-2020.nc'.format(path, var)	
	data = netCDF4.Dataset(arq)		
	var  = data.variables[dict_var[var]][:]
	lat  = data.variables['latitude'][:]
	lon  = data.variables['longitude'][:]

	std_clim  = np.std(var[:][:,:,:], axis=0)
	mean_clim  = np.nanmean(var[:][:,:,:], axis=0)
	
	# OCD
	# SUMMER DJF
	#~ D_ii  = np.nanmean(var[:][2872:2961,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ D_i   = np.nanmean(var[:][2882:2971,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ D     = np.nanmean(var[:][2892:2981,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ Di    = np.nanmean(var[:][2902:2991,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ Dii   = np.nanmean(var[:][2912:3001,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)

	#~ std_D_ii  = np.std(var[:][2872:2961,:,:], axis=0) 
	#~ std_D_i   = np.std(var[:][2882:2971,:,:], axis=0)
	#~ std_D     = np.std(var[:][2892:2981,:,:], axis=0)
	#~ std_Di    = np.std(var[:][2902:2991,:,:], axis=0) 
	#~ std_Dii   = np.std(var[:][2912:3001,:,:], axis=0) 
	
	#~ mean_D_ii  = np.nanmean(var[:][2872:2961,:,:], axis=0) 
	#~ mean_D_i   = np.nanmean(var[:][2882:2971,:,:], axis=0)
	#~ mean_D     = np.nanmean(var[:][2892:2985,:,:], axis=0)
	#~ mean_Di    = np.nanmean(var[:][2902:2991,:,:], axis=0) 
	#~ mean_Dii   = np.nanmean(var[:][2912:3001,:,:], axis=0) 

	# AUTUMN MAM
	#~ D_ii  = np.nanmean(var[:][2962:3052,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ D_i   = np.nanmean(var[:][2972:3062,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ D     = np.nanmean(var[:][2982:3072,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ Di    = np.nanmean(var[:][2992:3082,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ Dii   = np.nanmean(var[:][3002:3092,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)

	#~ std_D_ii  = np.std(var[:][2962:3052,:,:], axis=0) 
	#~ std_D_i   = np.std(var[:][2972:3062,:,:], axis=0)
	#~ std_D     = np.std(var[:][2982:3072,:,:], axis=0)
	#~ std_Di    = np.std(var[:][2992:3082,:,:], axis=0) 
	#~ std_Dii   = np.std(var[:][3002:3092,:,:], axis=0) 

	#~ mean_D_ii  = np.nanmean(var[:][2962:3052,:,:], axis=0) 
	#~ mean_D_i   = np.nanmean(var[:][2972:3062,:,:], axis=0)
	#~ mean_D     = np.nanmean(var[:][2982:3072,:,:], axis=0)
	#~ mean_Di    = np.nanmean(var[:][2992:3082,:,:], axis=0) 
	#~ mean_Dii   = np.nanmean(var[:][3002:3092,:,:], axis=0) 

	# WINTER JJA
	#~ D_ii  = np.nanmean(var[:][326:330,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ D_i   = np.nanmean(var[:][328:332,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ D     = np.nanmean(var[:][330:334,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ Di    = np.nanmean(var[:][332:336,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ Dii   = np.nanmean(var[:][334:338,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)

	#~ std_D_ii  = np.std(var[:][326:330,:,:], axis=0) 
	#~ std_D_i   = np.std(var[:][328:332,:,:], axis=0)
	#~ std_D     = np.std(var[:][330:334,:,:], axis=0)
	#~ std_Di    = np.std(var[:][332:336,:,:], axis=0) 
	#~ std_Dii   = np.std(var[:][334:338,:,:], axis=0) 

	#~ mean_D_ii  = np.nanmean(var[:][326:330,:,:], axis=0) 
	#~ mean_D_i   = np.nanmean(var[:][328:332,:,:], axis=0)
	#~ mean_D     = np.nanmean(var[:][330:334,:,:], axis=0)
	#~ mean_Di    = np.nanmean(var[:][332:336,:,:], axis=0) 
	#~ mean_Dii   = np.nanmean(var[:][334:338,:,:], axis=0) 

	# SPRING SON
	#~ D_ii  = np.nanmean(var[:][2781:2871,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ D_i   = np.nanmean(var[:][2791:2881,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ D     = np.nanmean(var[:][2801:2891,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ Di    = np.nanmean(var[:][2811:2901,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ Dii   = np.nanmean(var[:][2821:2911,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)

	#~ std_D_ii  = np.std(var[:][2781:2871,:,:], axis=0) 
	#~ std_D_i   = np.std(var[:][2791:2881,:,:], axis=0)
	#~ std_D     = np.std(var[:][2801:2891,:,:], axis=0)
	#~ std_Di    = np.std(var[:][2811:2901,:,:], axis=0) 
	#~ std_Dii   = np.std(var[:][2821:2911,:,:], axis=0) 

	#~ mean_D_ii  = np.nanmean(var[:][2781:2871,:,:], axis=0) 
	#~ mean_D_i   = np.nanmean(var[:][2791:2881,:,:], axis=0)
	#~ mean_D     = np.nanmean(var[:][2801:2891,:,:], axis=0)
	#~ mean_Di    = np.nanmean(var[:][2811:2901,:,:], axis=0) 
	#~ mean_Dii   = np.nanmean(var[:][2821:2911,:,:], axis=0) 
	
	# OCN
	# SUMMER DJF
	#~ D_ii  = np.nanmean(var[:][680:770,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ D_i   = np.nanmean(var[:][690:780,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ D     = np.nanmean(var[:][700:790,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ Di    = np.nanmean(var[:][710:800,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ Dii   = np.nanmean(var[:][720:810,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)

	#~ std_D_ii  = np.std(var[:][680:770,:,:], axis=0) 
	#~ std_D_i   = np.std(var[:][690:780,:,:], axis=0)
	#~ std_D     = np.std(var[:][700:790,:,:], axis=0)
	#~ std_Di    = np.std(var[:][710:800,:,:], axis=0) 
	#~ std_Dii   = np.std(var[:][720:810,:,:], axis=0) 
	
	#~ mean_D_ii  = np.nanmean(var[:][680:770,:,:], axis=0) 
	#~ mean_D_i   = np.nanmean(var[:][690:780,:,:], axis=0)
	#~ mean_D     = np.nanmean(var[:][700:790,:,:], axis=0)
	#~ mean_Di    = np.nanmean(var[:][710:800,:,:], axis=0) 
	#~ mean_Dii   = np.nanmean(var[:][720:810,:,:], axis=0) 

	# AUTUMN MAM
	#~ D_ii  = np.nanmean(var[:][770:860,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ D_i   = np.nanmean(var[:][780:870,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ D     = np.nanmean(var[:][790:880,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ Di    = np.nanmean(var[:][800:890,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ Dii   = np.nanmean(var[:][810:900,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)

	#~ std_D_ii  = np.std(var[:][770:860,:,:], axis=0) 
	#~ std_D_i   = np.std(var[:][780:870,:,:], axis=0)
	#~ std_D     = np.std(var[:][790:880,:,:], axis=0)
	#~ std_Di    = np.std(var[:][800:890,:,:], axis=0) 
	#~ std_Dii   = np.std(var[:][810:900,:,:], axis=0) 

	#~ mean_D_ii  = np.nanmean(var[:][770:860,:,:], axis=0) 
	#~ mean_D_i   = np.nanmean(var[:][780:870,:,:], axis=0)
	#~ mean_D     = np.nanmean(var[:][790:880,:,:], axis=0)
	#~ mean_Di    = np.nanmean(var[:][800:890,:,:], axis=0) 
	#~ mean_Dii   = np.nanmean(var[:][810:900,:,:], axis=0) 

	# WINTER JJA
	#~ D_ii  = np.nanmean(var[:][770:860,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ D_i   = np.nanmean(var[:][780:870,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ D     = np.nanmean(var[:][790:880,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ Di    = np.nanmean(var[:][800:890,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	#~ Dii   = np.nanmean(var[:][810:900,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)

	#~ std_D_ii  = np.std(var[:][770:860,:,:], axis=0) 
	#~ std_D_i   = np.std(var[:][780:870,:,:], axis=0)
	#~ std_D     = np.std(var[:][790:880,:,:], axis=0)
	#~ std_Di    = np.std(var[:][800:890,:,:], axis=0) 
	#~ std_Dii   = np.std(var[:][810:900,:,:], axis=0) 

	#~ mean_D_ii  = np.nanmean(var[:][770:860,:,:], axis=0) 
	#~ mean_D_i   = np.nanmean(var[:][780:870,:,:], axis=0)
	#~ mean_D     = np.nanmean(var[:][790:880,:,:], axis=0)
	#~ mean_Di    = np.nanmean(var[:][800:890,:,:], axis=0) 
	#~ mean_Dii   = np.nanmean(var[:][810:900,:,:], axis=0) 

	# SPRING SON
	D_ii  = np.nanmean(var[:][955:1005,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	D_i   = np.nanmean(var[:][965:1015,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	D     = np.nanmean(var[:][975:1025,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	Di    = np.nanmean(var[:][985:1035,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)
	Dii   = np.nanmean(var[:][995:1035,:,:], axis=0) - np.nanmean(var[:][:,:,:], axis=0)

	std_D_ii  = np.std(var[:][955:1005,:,:], axis=0) 
	std_D_i   = np.std(var[:][965:1015,:,:], axis=0)
	std_D     = np.std(var[:][975:1025,:,:], axis=0)
	std_Di    = np.std(var[:][985:1035,:,:], axis=0) 
	std_Dii   = np.std(var[:][995:1035,:,:], axis=0) 

	mean_D_ii  = np.nanmean(var[:][955:1005,:,:], axis=0) 
	mean_D_i   = np.nanmean(var[:][965:1015,:,:], axis=0)
	mean_D     = np.nanmean(var[:][975:1025,:,:], axis=0)
	mean_Di    = np.nanmean(var[:][985:1035,:,:], axis=0) 
	mean_Dii   = np.nanmean(var[:][995:1035,:,:], axis=0) 
		
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
	map.drawmeridians(np.arange(-80.,-30.,10.), size=5.5, labels=[0,0,0,1], linewidth=0.4, color='black')
	map.drawparallels(np.arange(-40.,10.,10.), size=5.5, labels=[1,0,0,0], linewidth=0.4, color='black') 

	lons, lats = np.meshgrid(new_lon, new_lat)
	xx, yy = map(lons,lats)

	path = '/home/nice/Documents/github_projects/shp'
	map.readshapefile('{0}/shp_world/world'.format(path), 'world', drawbounds=True, color='gray', linewidth=.5)
	map.readshapefile('{0}/lim_unid_fed/lim_unid_fed'.format(path), 'lim_unid_fed', drawbounds=True, color='black', linewidth=.5)
	
	return map, xx, yy
	
	
# Define event, season and variable
event    = u'OCN'
season   = u'SON'

variable_list = ['mtpr', 'mx2t', 'mn2t', 'mtnlwrf', 'q']
for variable in variable_list:
	
	dict_unit = {u'mtpr': u'PRE (mm d⁻¹)',
	u'mx2t': u'Tmax (°C)',
	u'mn2t': u'Tmin (°C)',
	u'mtnlwrf': u'ROL (W m⁻²)',
	u'q': u'Q (g kg⁻¹)'}

	# Import era5 database
	print(variable)
	print('Import era5 database')
	lat, lon, std_D_ii, std_D_i, std_D, std_Di, std_Dii, std_clim, mean_D_ii, mean_D_i, mean_D, mean_Di, mean_Dii, mean_clim, D_ii, D_i, D, Di, Dii = import_era5(variable) 

	# Plot maps with the function
	print('Plot maps with the function')
	fig = plt.figure(figsize=(8, 2))

	if variable == 'mtpr':
		cor_map = mpl.cm.BrBG
		levs1  = [-8, -6, -4, -2, 2, 4, 6, 8]

	elif variable == 'mx2t':
		cor_map = mpl.cm.bwr
		levs1   = [-4, -3, -2, -1, 1, 2, 3, 4]
		
	elif variable == 'mn2t':
		cor_map = mpl.cm.bwr
		levs1   = [-4, -3, -2, -1, 1, 2, 3, 4]

	elif variable == 'mtnlwrf':
		cor_map = mpl.cm.RdBu
		levs1   = [-40, -30, -20, -10, 10, 20, 30, 40]
			
	else:
		cmap = mpl.cm.PiYG
		levs1 = [-4, -3, -2, -1, 1, 2, 3, 4]
		
	ax = fig.add_subplot(1, 5, 1)
	map, xx, yy = basemap(lat, lon)
	plt.title(u'{0} {1} \n {2} (D-2)'.format(dict_unit[variable], season, event), fontsize=6, fontweight='bold')
	map.contourf(xx, yy, D_ii, levels=levs1, latlon=True, cmap=cor_map)	
	p_value = ttest(std_D_ii, std_clim, mean_D_ii, mean_clim)
	p_value = ma.masked_where(p_value >= 0.05, p_value) 
	map.contourf(xx, yy, p_value, colors='none', hatches=['....'])

	ax = fig.add_subplot(1, 5, 2)
	map, xx, yy = basemap(lat, lon)
	plt.title(u'{0} {1} \n {2} (D-1)'.format(dict_unit[variable], season, event), fontsize=6, fontweight='bold')
	map.contourf(xx, yy, D_i, levels=levs1, latlon=True, cmap=cor_map)
	p_value = ttest(std_D_i, std_clim, mean_D_i, mean_clim)
	p_value = ma.masked_where(p_value >= 0.05, p_value) 
	map.contourf(xx, yy, p_value, colors='none', hatches=['....'])

	ax = fig.add_subplot(1, 5, 3)
	map, xx, yy = basemap(lat, lon)
	plt.title(u'{0} {1} \n {2} (D0)'.format(dict_unit[variable], season, event), fontsize=6, fontweight='bold')
	map.contourf(xx, yy, D, levels=levs1, latlon=True, cmap=cor_map)
	p_value = ttest(std_D, std_clim, mean_D, mean_clim)
	p_value = ma.masked_where(p_value >= 0.05, p_value) 
	map.contourf(xx, yy, p_value, colors='none', hatches=['....'])

	ax = fig.add_subplot(1, 5, 4)
	map, xx, yy = basemap(lat, lon)
	plt.title(u'{0} {1} \n {2} (D+1)'.format(dict_unit[variable], season, event), fontsize=6, fontweight='bold')
	map.contourf(xx, yy, Di, levels=levs1, latlon=True, cmap=cor_map)
	p_value = ttest(std_Di, std_clim, mean_Di, mean_clim)
	p_value = ma.masked_where(p_value >= 0.05, p_value) 
	map.contourf(xx, yy, p_value, colors='none', hatches=['....'])

	ax = fig.add_subplot(1, 5, 5)
	map, xx, yy = basemap(lat, lon)
	plt.title(u'{0} {1} \n {2} (D+2)'.format(dict_unit[variable], season, event), fontsize=6, fontweight='bold')
	map.contourf(xx, yy, Dii, levels=levs1, latlon=True, cmap=cor_map, extend='both') 
	cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
	cbar.ax.tick_params(labelsize=6) 
	p_value = ttest(std_Dii, std_clim, mean_Dii, mean_clim)
	p_value = ma.masked_where(p_value >= 0.05, p_value) 
	map.contourf(xx, yy, p_value, colors='none', hatches=['....'])

	# Path out to save figure
	print('Path out to save figure')
	path_out = '/home/nice/Downloads/janio/figs'
	name_out = 'pyplt_maps_composites_{0}_{1}_{2}_era5.png'.format(event, season, variable)
	if not os.path.exists(path_out):
		create_path(path_out)
	plt.savefig(os.path.join(path_out, name_out), dpi=200, bbox_inches='tight')
exit()	
