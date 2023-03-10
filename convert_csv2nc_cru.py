# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "01/05/2023"
__description__ = "This script convert .csv to .nc"

import os
import numpy as np
import pandas as pd
import xarray as xr
import netCDF4 as nc
import matplotlib.pyplot as plt

from netCDF4 import Dataset
from dict_stations_inmet import inmet

# create date list
dt = pd.date_range('1978-01-01','2020-12-31', freq='M').strftime('%Y-%m-%d').tolist()
																																																											
# This load the file (date, variable)
data = np.loadtxt(os.path.join('/home/nice/Downloads', 'PRE_PET_PARANA_CIC_78_20_CRU405_csv.csv'), dtype='str', delimiter=';', unpack=True)
data = data[:,1:]
data_values = np.where(data[4,:] == 'null', -999., data[4,:])
data_values = np.array(data_values, dtype=float)

data_dates  = []
for i in range(len(dt)):
		
	data_dates.append('{0}'.format(dt[i]))

	print('Date organized:', data_dates[i], data_values[i])

nc_output = '/home/nice/Downloads/pre_parana_CRU405_mon_1978-2020.nc'

# create netcdf
ds = Dataset(nc_output, mode='w', format='NETCDF4_CLASSIC')

ds.Conventions 	= 'CF-1.6'
ds.title 		= 'CRU TS4.05 Precipitation'
ds.institution 	= 'Data held at British Atmospheric Data Centre, RAL, UK'
ds.source 		= 'Run ID = 2103051243. Data generated from:pre.2103041709.dtb'
ds.history 		= 'Rewrote via python script'
ds.references 	= 'Information on the data is available at http://badc.nerc.ac.uk/data/cru/'
ds.comment 		= 'This script convert .csv to .nc'
	
ds.createDimension('time', None)

time 				= ds.createVariable('time', float, ('time'))
time.axis 			= 'L'
time.calendar 		= '360_day'
time.units			= 'months since 1978-01-01'
time[:]				= range(len(data_dates))

var 				= ds.createVariable('pre',  float, ('time'))
var.units 			= 'mm'
var.long_name 		= 'Precipitation'
var.standard_name 	= 'Precipitation'
var.missing_value 	= -999
var[:] 				= data_values[:]
	
ds.close()

if os.path.exists(nc_output): 
	print('Done -->', nc_output)
exit()

# open netcdf
dn = xr.open_dataset('/home/nice/Downloads/' + 'pet_parana_CRU405_mon_1978-2020.nc')
dn = dn.pet.sel(time=slice('1978-01-01','2020-12-01'))
dn = dn.groupby('time.month').mean('time')
dn_pet = dn.values
print(dn_pet)

ds = xr.open_dataset('/home/nice/Downloads/' + 'pre_parana_CRU405_mon_1978-2020.nc')
ds = ds.pre.sel(time=slice('1978-01-01','2020-12-01'))
ds = ds.groupby('time.month').mean('time')
ds_pre = ds.values
print(ds_pre)

print('Plot figure')
# Plot figure
fig = plt.figure()
time = np.arange(0.5, 12 + 0.5)
plt.plot(time, dn_pet, linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='black', label='PET')
plt.plot(time, ds_pre, linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='blue', label='PRE')
plt.title('Parana - CRU TS4.05 (1978-2020)')
plt.yticks(np.arange(0, 220, 20))
plt.xticks(time, ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'))
plt.xlabel('Annual Cycle', fontsize=8, fontweight='bold')
plt.legend(fontsize=8)
plt.grid()

print('Path out to save figure')
# Path out to save figure
path_out = '/home/nice/Downloads/'
name_out = 'pyplt_annual_cycle.png'
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()












