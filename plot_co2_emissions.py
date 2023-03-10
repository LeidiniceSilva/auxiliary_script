# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "03/03/2023"
__description__ = "This script CO2 emissions pathway"

import os
import numpy as np
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt

print('Read dataset')
hist = [23.4, 23.3, 23.1, 23.6, 24.8, 25.2, 25.6, 26.2, 27.0, 27.10, 27.4, 
27.8, 27.2, 27.4, 28.0, 28.5, 29.2, 31.1, 29.4, 29.9, 30.3, 30.1, 31.1, 32.9, 
33.6, 33.9, 35.2, 35.4, 36.2, 36.1, 37.6, 38.8, 39.2, 39.3, 39.8, 40.2, 39.3, 
39.7, 40.2, 40.9, 38.5, 40.2, 40.5]

b1 = [29.4, 34.3, 39.3, 41.1, 43.36, 41.69, 36.0, 30.2, 24.7, 19.6, 15.6]
a1 = [27.09, 35.5, 45.5, 55.5, 64.1, 71.0, 72.6, 74.8, 77.62, 68.19, 59.0]

rcp26 = [28.89, 29.88, 30.77, 31.71, 32.64, 33.58, 34.11, 34.63, 35.15, 
35.67, 36.19, 36.33, 36.47, 36.61, 36.7, 36.89, 37.07, 37.17, 37.31, 37.45,
37.59, 29.11, 18.41, 12.41, 7.45, 2.40, 0.43, -0.98, -1.54]

rcp85 = [28.89, 29.88, 30.77, 31.71, 32.64, 33.58, 34.11, 34.63, 35.35, 
35.94, 36.53, 37.43, 38.34, 39.25, 40.15, 41.06, 41.97, 42.87,43.78, 44.69, 
45.59, 53.33, 63.87, 76.14, 88.29, 96.63, 101.55, 104.54, 105.59]

ssp126 = [39.15, 39.28, 39.41, 39.54, 39.67, 39.80, 34.73, 26.51, 17.96, 10.53, 4.48, -3.29, -8.39, -8.62]
ssp585 = [39.15, 40.06, 40.98, 41.89, 42.80, 43.71, 55.30, 68.78, 83.30, 100.34, 116.81, 129.65, 130.58, 126.29]

hist_dt = pd.date_range(start="19800101", end="20230101", freq="Y")
hist_xy = pd.Series(data=hist, index=hist_dt)

sres_dt = pd.date_range(start="20000101", end="21000101", freq="10YS")
a1_xy = pd.Series(data=a1, index=sres_dt)
b1_xy = pd.Series(data=b1, index=sres_dt)

rcp_dt_i = pd.date_range(start="20000101", end="20200101", freq="Y")
rcp26_xy_i = pd.Series(data=rcp26[:20], index=rcp_dt_i)
rcp85_xy_i = pd.Series(data=rcp85[:20], index=rcp_dt_i)

rcp_dt = pd.date_range(start="20200101", end="21000101", freq="10YS")
rcp26_xy = pd.Series(data=rcp26[20:], index=rcp_dt)
rcp85_xy = pd.Series(data=rcp85[20:], index=rcp_dt)

ssp_dt_i = pd.date_range(start="20150101", end="20200101", freq="Y")
ssp126_xy_i = pd.Series(data=ssp126[:5], index=ssp_dt_i)
ssp585_xy_i = pd.Series(data=ssp585[:5], index=ssp_dt_i)

ssp_dt = pd.date_range(start="20200101", end="21000101", freq="10YS")
ssp126_xy = pd.Series(data=ssp126[5:], index=ssp_dt)
ssp585_xy = pd.Series(data=ssp585[5:], index=ssp_dt)

print('Plot figure')
# Plot figure
fig = plt.figure(figsize=(10, 6))

ax = fig.add_subplot(1, 1, 1)
plt.plot(hist_xy, linewidth=1.5, linestyle='-', color='black', label='Historical')

plt.plot(b1_xy, linewidth=1.5, linestyle='-.', color='blue', label='B1')
plt.plot(a1_xy, linewidth=1.5, linestyle='-.', color='red', label='A1')

plt.plot(rcp26_xy_i, linewidth=1.5, linestyle='--', color='blue')
plt.plot(rcp85_xy_i, linewidth=1.5, linestyle='--', color='red')
plt.plot(rcp26_xy, linewidth=1.5, linestyle='--', color='blue', label='RCP2.6')
plt.plot(rcp85_xy, linewidth=1.5, linestyle='--', color='red', label='RCP8.5')

plt.plot(ssp126_xy_i, linewidth=1.5, linestyle=':', color='blue')
plt.plot(ssp585_xy_i, linewidth=1.5, linestyle=':', color='red')
plt.plot(ssp126_xy, linewidth=1.5, linestyle=':', color='blue', label='SSP1-2.6')
plt.plot(ssp585_xy, linewidth=1.5, linestyle=':', color='red', label='SSP5-8.5')

plt.axvline(19400, linewidth=1., linestyle='--',  color='black')
plt.axhline(linewidth=1., linestyle='--',  color='black')
plt.axhspan(0, -20, color='blue', alpha=0.2, lw=0)
plt.axhspan(0, 140, color='red', alpha=0.2, lw=0)
plt.ylim(-20, 140)
plt.yticks(np.arange(-20, 160, 20))
plt.xlabel('Year')
plt.ylabel('Global emissions (GtCO2/yr)')
plt.text(20000, 135, '2022')
plt.text(4000, 35, 'Historical emissions')
plt.text(20000, 15, 'Projections')
plt.text(2060, -10, 'Net-negative global emissions')
plt.grid(axis='y')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
# ~ ax.tick_params(axis='both',which='both',left=False)
plt.legend(bbox_to_anchor=(0.5, 1.075), loc=9, ncol=7, frameon=False)

print('Path out to save figure')
# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_co2_emissions.png'
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()
