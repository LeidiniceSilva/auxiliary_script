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

from scipy.stats import linregress

var = 'Tmax'
region = 'SE'
period = '1991_2023'
path = '/afs/ictp.it/home/m/mda_silv/Documents/Heat_wave'


def import_weather_stations(region):

	weather_data = pd.read_csv('{0}/INMET/{1}_inmet_{2}_diario_{3}.csv'.format(path, var, period, region), parse_dates=['Data'])	
	columns = ['Date'] + [f'Station_{i}' for i in range(1, weather_data.shape[1])]
	weather_data.columns = columns

	# Extract year and day of year
	weather_data['Year'] = weather_data['Date'].dt.year
	weather_data['DayOfYear'] = weather_data['Date'].dt.dayofyear
	
	# Calculate the daily average maximum temperature across all stations
	stations = [col for col in weather_data.columns if col.startswith('Station_')]
	weather_data['AvgMaxTemp'] = weather_data[stations].mean(axis=1)

	# Define the baseline period
	baseline_start = '1991-01-01'
	baseline_end = '2020-12-31'

	# Filter data for the baseline period
	baseline_data = weather_data[(weather_data['Date'] >= baseline_start) & (weather_data['Date'] <= baseline_end)]

	# Calculate the 90th percentile centered on a 5-day window for the baseline period
	weather_percentiles = baseline_data.groupby('DayOfYear')['AvgMaxTemp'].apply(lambda x: x.rolling(window=5, center=True, min_periods=1).quantile(0.9)).groupby('DayOfYear').mean()
		
	return weather_data, weather_percentiles
	
	
def identify_heat_waves(data, percentiles, min_days=6):

	heat_wave_dates = []
	heat_wave_counts = {year: 0 for year in range(1991, 2024)}
	
	for year in range(1991, 2024):
		yearly_data = data[data['Year'] == year]
		above_percentile = yearly_data.apply(lambda row: row['AvgMaxTemp'] > percentiles[row['DayOfYear']], axis=1)

		# Calculate consecutive days above the 90th percentile
		consecutive_days = above_percentile * (above_percentile.groupby((above_percentile != np.roll(above_percentile, 1)).cumsum()).cumcount() + 1)
		
		# Identify heat wave periods
		heat_waves = consecutive_days >= min_days
	
		# Identify the start and end dates of each heat wave
		current_wave = []
		for i, is_heat_wave in enumerate(heat_waves):
			
			if is_heat_wave:
				current_wave.append(yearly_data.iloc[i]['Date'])
			elif current_wave:
				if len(current_wave) >= min_days:
					heat_wave_dates.append((current_wave[0], current_wave[-1]))
				current_wave = []
			
		# If a heat wave continues until the end of the year
		if current_wave and len(current_wave) >= min_days:
			heat_wave_dates.append((current_wave[0], current_wave[-1]))
	    
	return heat_wave_dates


# Import weather stations and 90th percentile
weather_data, weather_percentile = import_weather_stations(region)

# Import heat wave dates
hw_dates = identify_heat_waves(weather_data, weather_percentile)

# Save heat waves episodes date file 
with open('{0}/figs/heat_wave_episodes_{1}.txt'.format(path, region), 'w') as file:
	
	for start_date, end_date in hw_dates:
		
		episode = f"Heat wave from {start_date.date()} to {end_date.date()}"
		print(episode)
		
		file.write(episode + '\n')
		
print('Heat wave episodes have been saved')
exit()

