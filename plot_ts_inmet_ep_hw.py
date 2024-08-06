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

path = '/home/esp-shared-a/Distribution/Users/epichell/mdasilva/Study_case/Heat_wave'

def comp_days_hw(region):

	# Open file 
	df = pd.read_csv(os.path.join('{0}/INMET'.format(path), '{0}_inmet_{1}_diario_{2}.csv'.format(var, period, region)))
	columns = ['Date'] + [f'Station_{i}' for i in range(1, df.shape[1])]
	df.columns = columns

	# Convert the Date column to datetime format
	df['Date'] = pd.to_datetime(df['Date'])
	df.set_index('Date', inplace=True)

	# List of station columns
	stations = [col for col in df.columns if col.startswith('Station_')]

	# Separate the data into baseline and target years
	baseline_data = df[(df.index.year >= 1991) & (df.index.year <= 2020)]
	target_data = df[df.index.year >= 1991]

	# Compute the 90th percentile for each day of the year in the baseline period
	baseline_data['day_of_year'] = baseline_data.index.dayofyear
	percentile_90th = baseline_data.groupby('day_of_year')[stations].quantile(0.90)

	# Add a 'day_of_year' column to the target data
	target_data['day_of_year'] = target_data.index.dayofyear

	# Merge target data with the 90th percentile thresholds
	target_data = target_data.join(percentile_90th, on='day_of_year', rsuffix='_90th')

	# Initialize dictionary to count heat wave episodes per year
	heat_wave_counts = {year: 0 for year in range(1991, 2024)}

	# Identify heat wave episodes
	for station in stations:
    		target_data[f'{station}_5day_ma'] = target_data[station].rolling(window=5).mean()
    		target_data[f'{station}_heat_wave'] = target_data[f'{station}_5day_ma'] > target_data[f'{station}_90th']

    		for year in range(1991, 2024):
       			year_data = target_data[target_data.index.year == year]
        		in_heat_wave = False
        		for date, row in year_data.iterrows():
            			if row[f'{station}_heat_wave']:
                			if not in_heat_wave:
                    				heat_wave_counts[year] += 1
                    				in_heat_wave = True
            			else:
                			in_heat_wave = False

	# Compute average number of heat waves across all stations for each year
	average_heat_wave_counts = {year: count / len(stations) for year, count in heat_wave_counts.items()}
	
	return average_heat_wave_counts

# Import dataset
avg_hw_counts_CO = comp_days_hw('CO')
avg_hw_counts_N = comp_days_hw('N')
avg_hw_counts_NE = comp_days_hw('NE')
avg_hw_counts_SE = comp_days_hw('SE')

# Plot figure 
# Plot figure 
fig = plt.figure(figsize=(12, 6))

years_CO = list(avg_hw_counts_CO.keys())
years_N = list(avg_hw_counts_N.keys())
years_NE = list(avg_hw_counts_NE.keys())
years_SE = list(avg_hw_counts_SE.keys())

counts_CO = list(avg_hw_counts_CO.values())
counts_N = list(avg_hw_counts_N.values())
counts_NE = list(avg_hw_counts_NE.values())
counts_SE = list(avg_hw_counts_SE.values())

# Subplot 1
ax=fig.add_subplot(2, 2, 1)
plt.bar(years_CO, counts_CO, color='blue', edgecolor='black', linewidth=1)
plt.title('(a)', loc='left', fontweight='bold')
plt.ylabel('Number of Heat Waves', fontweight='bold')
plt.xticks(rotation=45)
plt.grid(True, linestyle='--', alpha=0.8)
plt.setp(ax.get_xticklabels(), visible=False)

# Subplot 2
ax=fig.add_subplot(2, 2, 2)
plt.bar(years_N, counts_N, color='blue', edgecolor='black', linewidth=1)
plt.title('(b)', loc='left', fontweight='bold')
plt.xticks(rotation=45)
plt.grid(True, linestyle='--', alpha=0.8)
plt.setp(ax.get_xticklabels(), visible=False)

# Subplot 3
ax=fig.add_subplot(2, 2, 3)
plt.bar(years_NE, counts_NE, color='blue', edgecolor='black', linewidth=1)
plt.title('(c)', loc='left', fontweight='bold')
plt.xlabel('Year', fontweight='bold')
plt.ylabel('Number of Heat Waves', fontweight='bold')
plt.xticks(rotation=45)
plt.grid(True, linestyle='--', alpha=0.8)

# Subplot 4
ax=fig.add_subplot(2, 2, 4)
plt.bar(years_SE, counts_SE, color='blue', edgecolor='black', linewidth=1)
plt.title('(d)', loc='left', fontweight='bold')
plt.xlabel('Year', fontweight='bold')
plt.xticks(rotation=45)
plt.grid(True, linestyle='--', alpha=0.8)

# Path out to save figure
path_out = '{0}/figs'.format(path)
name_out = 'pyplt_regions_ep_hw_SON-2023.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
