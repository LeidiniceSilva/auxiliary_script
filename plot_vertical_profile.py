# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "26/04/2019"
__description__ = "Plot potential temperature ERA5 and RegCM dataset"

import os
import netCDF4
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.font_manager import FontProperties

def import_era5(area, season):
	
	mod_path = '/home/nice/Downloads'
	arq1     = '{0}/tatm_{1}_era5_{2}_2001-2010.nc'.format(mod_path, area, season)
	data1    = netCDF4.Dataset(arq1)
	var1     = data1.variables['t']
	lev      = data1.variables['level']
	lat      = data1.variables['lat']
	lon      = data1.variables['lon']
	
	level_list = [36, 35, 34, 33, 32, 30, 28, 26, 24, 23, 21, 19, 17, 16, 14, 11, 9, 6]

	clim_season = []
	for level in level_list:
		era5 = np.nanmean(np.nanmean((var1[:,level,:,:]), axis=1), axis=1)
		clim_season.append(era5)

	return np.squeeze(clim_season)


clim_namz_djf = import_era5('nam', 'djf')
clim_samz_djf = import_era5('sam', 'djf')
clim_neb_djf = import_era5('neb', 'djf')

clim_namz_jja = import_era5('nam', 'jja')
clim_samz_jja = import_era5('sam', 'jja')
clim_neb_jja = import_era5('neb', 'jja')

print(clim_namz_djf, clim_samz_djf, clim_neb_djf, clim_namz_jja, clim_samz_jja, clim_neb_jja)
exit()

# Plot model end obs data climatology
fig = plt.figure(figsize=(8,8))
y =  [20, 70, 125, 200, 250, 350, 450, 500, 600, 650, 750, 800, 850, 900, 925, 950, 975, 1000]
line_labels = ['Reg_Exp1', 'Reg_Exp2', 'ERA5']

# Subplot one
plt.subplot(231)
l1 = plt.plot(clim_namz_djf, y, marker='o', color="blue")[0]
l2 = plt.plot(clim_namz_djf, y, marker='o', color="red")[0]
l3 = plt.plot(clim_namz_djf, y, marker='o', color="black")[0]
plt.ylabel('Pressão (mbar)', fontsize=8, fontweight='bold')
plt.yticks(y, ('1000', ' ', '950', ' ', '900', ' 850', '750', ' ', '600', ' ', '450', ' ', '250', ' ', '125', ' ', ' ', '10'))
plt.tick_params(axis='both', which='major', labelsize=8, length=4., width=1., pad=1., labelcolor='black')
plt.text(240., 975., u'A) NAMZ (DJF)', fontsize=8, fontweight='bold')
plt.legend([l1, l2, l3], labels=line_labels,loc='lower left', shadow=True, ncol=1, prop=FontProperties(size=8))

# Subplot two
plt.subplot(232)
l4 = plt.plot(clim_samz_djf, y, marker='o', color="blue")[0]
l5 = plt.plot(clim_samz_djf, y, marker='o', color="red")[0]
l6 = plt.plot(clim_samz_djf, y, marker='o', color="black")[0]
plt.title('Perfil Vertical de Temperatura Potencial', fontsize=8, fontweight='bold')
plt.yticks(y, ('1000', ' ', '950', ' ', '900', ' 850', '750', ' ', '600', ' ', '450', ' ', '250', ' ', '125', ' ', ' ', '10'))
plt.tick_params(axis='both', which='major', labelsize=8, length=4., width=1., pad=1., labelcolor='black')
plt.text(240., 975., u'B) SAMZ (DJF)', fontsize=8, fontweight='bold')

# Subplot three
plt.subplot(233)
l7 = plt.plot(clim_neb_djf, y, marker='o', color="blue")[0]
l8 = plt.plot(clim_neb_djf, y, marker='o', color="red")[0]
l9 = plt.plot(clim_neb_djf, y, marker='o', color="black")[0]
plt.yticks(y, ('1000', ' ', '950', ' ', '900', ' 850', '750', ' ', '600', ' ', '450', ' ', '250', ' ', '125', ' ', ' ', '10'))
plt.tick_params(axis='both', which='major', labelsize=8, length=4., width=1., pad=1., labelcolor='black')
plt.text(240., 975., u'C) NEB (DJF)', fontsize=8, fontweight='bold')

# Subplot four
plt.subplot(234)
l1 = plt.plot(clim_namz_jja, y, marker='o', color="blue")[0]
l2 = plt.plot(clim_namz_jja, y, marker='o', color="red")[0]
l3 = plt.plot(clim_namz_jja, y, marker='o', color="black")[0]
plt.xlabel('Temperatura Potencial (K)', fontsize=8, fontweight='bold')
plt.ylabel('Pressão (mbar)', fontsize=8, fontweight='bold')
plt.yticks(y, ('1000', ' ', '950', ' ', '900', ' 850', '750', ' ', '600', ' ', '450', ' ', '250', ' ', '125', ' ', ' ', '10'))
plt.tick_params(axis='both', which='major', labelsize=8, length=4., width=1., pad=1., labelcolor='black')
plt.text(240., 975., u'C) NAMZ (JJA)', fontsize=8, fontweight='bold')

# Subplot five
plt.subplot(235)
l4 = plt.plot(clim_samz_jja, y, marker='o', color="blue")[0]
l5 = plt.plot(clim_samz_jja, y, marker='o', color="red")[0]
l6 = plt.plot(clim_samz_jja, y, marker='o', color="black")[0]
plt.xlabel('Temperatura Potencial (K)', fontsize=8, fontweight='bold')
plt.yticks(y, ('1000', ' ', '950', ' ', '900', ' 850', '750', ' ', '600', ' ', '450', ' ', '250', ' ', '125', ' ', ' ', '10'))
plt.tick_params(axis='both', which='major', labelsize=8, length=4., width=1., pad=1., labelcolor='black')
plt.text(240., 975., u'D) SAMZ (JJA)', fontsize=8, fontweight='bold')

# Subplot six
plt.subplot(236)
l7 = plt.plot(clim_neb_jja, y, marker='o', color="blue")[0]
l8 = plt.plot(clim_neb_jja, y, marker='o', color="red")[0]
l9 = plt.plot(clim_neb_jja, y, marker='o', color="black")[0]
plt.xlabel('Temperatura Potencial (K)', fontsize=8, fontweight='bold')
plt.yticks(y, ('1000', ' ', '950', ' ', '900', ' 850', '750', ' ', '600', ' ', '450', ' ', '250', ' ', '125', ' ', ' ', '10'))
plt.tick_params(axis='both', which='major', labelsize=8, length=4., width=1., pad=1., labelcolor='black')
plt.text(240., 975., u'E) NEB (JJA)', fontsize=8, fontweight='bold')

# Path out to save figure
path_out = '/home/nice'
name_out = 'pyplt_perf_vert_regcm_pbl_obs_2001-2010.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')

plt.show()
exit()
