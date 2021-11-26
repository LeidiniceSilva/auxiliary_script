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

level_list1 = [34, 33, 32, 30, 28, 26, 24, 23, 21, 19, 17, 16, 14, 11, 9, 6]
level_list2 = [33, 32, 31, 29, 28, 27, 26, 24, 23, 22, 20, 17, 16, 13, 10, 8]

clim1= level_list1
clim2 = level_list2

# Plot model end obs data climatology
fig = plt.figure()
y =  [100, 200, 250, 350, 450, 500, 600, 650, 750, 800, 850, 900, 925, 950, 975, 1000]
line_labels = ['RegCM5', 'OBS']

# Subplot one
plt.subplot(121)
l1 = plt.plot(clim1, y, marker='o', color="blue")
l2 = plt.plot(clim2, y, marker='o', color="red")
plt.title('A) Perfil Vertical - AMZ', fontsize=8, fontweight='bold')
plt.xlabel('Temperatura (°C)', fontsize=8, fontweight='bold')
plt.ylabel('Pressão (hPa)', fontsize=8, fontweight='bold')
plt.yticks(y, ('1000', ' ', '950', ' ', '900', ' 850', ' ', '750', ' ', '500', ' ', ' 350', ' ', ' ', ' ', '100'))
plt.tick_params(axis='both', which='major', labelsize=8, length=4., width=1., pad=1., labelcolor='black')
plt.legend([l1, l2], labels=line_labels,loc='lower left', shadow=True, ncol=1, prop=FontProperties(size=8))

# Subplot two
plt.subplot(122)
l1 = plt.plot(clim1, y, marker='o', color="blue")
l2 = plt.plot(clim2, y, marker='o', color="red")
plt.title('B) Perfil Vertical - LPB', fontsize=8, fontweight='bold')
plt.xlabel('Temperatura (°C)', fontsize=8, fontweight='bold')
plt.ylabel('Pressão (hPa)', fontsize=8, fontweight='bold')
plt.yticks(y, ('1000', ' ', '950', ' ', '900', ' 850', ' ', '750', ' ', '500', ' ', ' 350', ' ', ' ', ' ', '100'))
plt.tick_params(axis='both', which='major', labelsize=8, length=4., width=1., pad=1., labelcolor='black')

# Path out to save figure
path_out = '/home/nice'
name_out = 'pyplt_perf_vert_temp.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')

plt.show()
exit()
