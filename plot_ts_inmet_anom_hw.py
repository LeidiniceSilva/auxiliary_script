# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 04, 2024"
__description__ = "This script plot bias maps"

import os
import netCDF4
import numpy as np
import pandas as pd
import pandas as pd
import matplotlib.cm as cm
import matplotlib.pyplot as plt

var = 'Tmax'
period = '1991_2023'

path = '/afs/ictp.it/home/m/mda_silv/Documents/Heat_wave'

def comp_anom(region):

	# Open file 
	data = pd.read_csv(os.path.join('{0}/INMET'.format(path), '{0}_inmet_{1}_diario_{2}.csv'.format(var, period, region)))
	columns = ['Date'] + [f'Station_{i}' for i in range(1, data.shape[1])]
	data.columns = columns

	# Convert the Date column to datetime format
	data['Date'] = pd.to_datetime(data['Date'])

	# Extract year, month, and day from the Date column
	data['Year'] = data['Date'].dt.year
	data['Month'] = data['Date'].dt.month
	data['Day'] = data['Date'].dt.day

	# List of station columns
	station_columns = [col for col in data.columns if col.startswith('Station_')]

	# Compute daily mean temperatures across all stations
	data['Mean_Temperature'] = data[station_columns].mean(axis=1)

	# Filter data to exclude 2023 for baseline calculation
	baseline_data = data[(data['Year'] >= 1991) & (data['Year'] <= 2020)]

	# Calculate daily mean temperatures for the baseline period (1991-2020) across all stations
	daily_mean_temps = baseline_data.groupby(['Month', 'Day'])['Mean_Temperature'].mean().reset_index()

	# Filter data for SON (September, October, November) 2023
	son_2023 = data[(data['Year'] == 2023) & (data['Month'].isin([9, 10, 11]))]

	# Calculate daily average temperatures for SON 2023 across all stations
	daily_avg_temps_2023 = son_2023.groupby(['Month', 'Day'])['Mean_Temperature'].mean().reset_index()

	# Merge SON 2023 daily averages with baseline daily mean temperatures
	daily_avg_temps_2023 = daily_avg_temps_2023.merge(daily_mean_temps, on=['Month', 'Day'], suffixes=('_2023', '_mean'))

	# Calculate anomalies
	daily_avg_temps_2023['Anomaly'] = daily_avg_temps_2023['Mean_Temperature_2023'] - daily_avg_temps_2023['Mean_Temperature_mean']

	# Create a Date column for plotting
	daily_avg_temps_2023['Date'] = pd.to_datetime(daily_avg_temps_2023[['Month', 'Day']].assign(Year=2023))

	return daily_avg_temps_2023


# Import dataset
day_anom_son_2023_CO = comp_anom('CO')
day_anom_son_2023_SE = comp_anom('SE')

# Plot figure 
fig = plt.figure()

colors_CO = ['red' if value > 0 else 'blue' for value in day_anom_son_2023_CO['Anomaly']]
colors_SE = ['red' if value > 0 else 'blue' for value in day_anom_son_2023_SE['Anomaly']]

# Subplot 1
ax=fig.add_subplot(2, 1, 1)
plt.bar(day_anom_son_2023_CO['Date'], day_anom_son_2023_CO['Anomaly'], color=colors_CO, edgecolor='black', linewidth=1)
plt.title('(a)', loc='left', fontweight='bold')
plt.ylabel('Tmax anomaly (°C)', fontweight='bold')
plt.ylim(-10,10)
plt.grid(True, linestyle='--', alpha=0.8)
plt.axhline(0, color='black', linewidth=0.8)
plt.xticks(rotation=45)
plt.setp(ax.get_xticklabels(), visible=False)

# Subplot 2
ax=fig.add_subplot(2, 1, 2)
plt.bar(day_anom_son_2023_SE['Date'], day_anom_son_2023_SE['Anomaly'], color=colors_SE, edgecolor='black', linewidth=1)
plt.title('(b)', loc='left', fontweight='bold')
plt.ylabel('Tmax anomaly (°C)', fontweight='bold')
plt.ylim(-10,10)
plt.grid(True, linestyle='--', alpha=0.8)
plt.axhline(0, color='black', linewidth=0.8)
plt.xticks(rotation=45)
	
# Path out to save figure
path_out = '{0}/figs'.format(path)
name_out = 'pyplt_regions_anom_hw_SON-2023.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
