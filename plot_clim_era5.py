# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "11/22/2019"
__description__ = "Plotting average monthly data on Prec, OLR and UV10m"

import os
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
	
# Export variable
path_var = '/home/nice/Downloads/era5/'

tp_var = xr.open_mfdataset(path_var + 'tp_era5_mon_1979_2019.nc', combine='by_coords')['tp']
mtnlwrf_var = xr.open_mfdataset(path_var + 'mtnlwrf_era5_mon_1979_2019.nc', combine='by_coords')['mtnlwrf']
si10_var = xr.open_mfdataset(path_var + 'si10_era5_mon_1979_2019.nc', combine='by_coords')['si10']

# Select station
Name = ['Fortaleza-CE']
lat_lon = [[-3.75, -38.55],]               
                                  
# Calculate month mean and select lat and lon
for n, name in enumerate(Name):

	tp_mon = [242,243,271,167,55,20,12,8,53,144,223,263]

	mtnlwrf_lat_lon = mtnlwrf_var.sel(latitude=lat_lon[n][0], longitude=lat_lon[n][1], method='nearest')
	mtnlwrf_mon = mtnlwrf_lat_lon.groupby('time.month').mean('time')

	si10_lat_lon = si10_var.sel(latitude=lat_lon[n][0], longitude=lat_lon[n][1], method='nearest')
	si10_mon = si10_lat_lon.groupby('time.month').mean('time')

# plot data ERA5
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax3 = ax1.twinx()

# Make space on the right for the additional axes
fig.subplots_adjust(right=0.8)

# Move additional axes to free space on the right
ax3.spines["right"].set_position(("axes", 1.175))

time = np.arange(0.5, 12 + 0.5)
line1, = ax1.plot(time, tp_mon, '--', color='Black', label=u'Prec')
line2, = ax2.plot(time, mtnlwrf_mon, '-.', color='red', label=u'OLR')
line3, = ax3.plot(time, si10_mon, ':', color='blue', label=u'Wind 10m')

# Set axes limits
ax1.set_xlim(0, 12)
ax1.set_ylim(0, 300)
ax2.set_ylim(-300, -200)
ax3.set_ylim(2, 12)

# Set axes colors
ax1.spines["left"].set_color(line1.get_color())
ax2.spines["right"].set_color(line2.get_color())
ax3.spines["right"].set_color(line3.get_color())

ax1.set_xlabel(u'Months')
ax1.set_ylabel(r'Precipitation (mm)')
ax2.set_ylabel(r'Outgoing Longwave Radiation (W m⁻²)', color='red')
ax3.set_ylabel(r'10 Metre Wind Speed (m s⁻¹)', color='blue')

ax1.set_xticks(np.arange(0.5, 12 + 0.5))
ax1.set_xticklabels(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])
plt.setp(ax1.get_xticklabels())
plt.title('A) Fortaleza-CE (Lat=-3.75, Lon=-38.55) - ERA5 (1979-2019)')

# Set up ticks and grid lines
ax1.minorticks_on()
ax2.minorticks_on()
ax3.minorticks_on()
ax1.tick_params(direction="in", which="both", colors=line1.get_color())
ax2.tick_params(direction="in", which="both", colors=line2.get_color())
ax3.tick_params(direction="in", which="both", colors=line3.get_color())

ax1.xaxis.grid(True, which='major', linestyle='--', linewidth='0.8', zorder=0.4)
ax1.yaxis.grid(True, which='major', linestyle='--', linewidth='0.8', zorder=0.4)
	
lines = line1, line2, line3
labels = [l.get_label() for l in lines]
plt.legend(lines, labels, loc='upper center', ncol=3)

# Path out to save figure
path_out = '/home/nice/Downloads/era5'
name_out = 'pyplt_clim_era5_fortaleza_1979-2019.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()






