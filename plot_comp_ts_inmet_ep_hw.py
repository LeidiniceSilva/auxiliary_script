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
		consecutive_days = above_percentile * (above_percentile.groupby((above_percentile != above_percentile.shift()).cumsum()).cumcount() + 1)
		
		# Identify heat wave periods
		heat_waves = consecutive_days >= min_days
		
		# Count and identify the start and end dates of each heat wave
		heat_wave_count = 0
		current_wave = []
		
		for i, is_heat_wave in enumerate(heat_waves):
			
			if is_heat_wave:
				current_wave.append(yearly_data.iloc[i]['Date'])
			elif current_wave:
				if len(current_wave) >= min_days:
					heat_wave_dates.append((current_wave[0], current_wave[-1]))
					heat_wave_count += 1

				current_wave = []
		
		# If a heat wave continues until the end of the year
		if current_wave and len(current_wave) >= min_days:
			
			heat_wave_dates.append((current_wave[0], current_wave[-1]))
			heat_wave_count += 1
			
		heat_wave_counts[year] = heat_wave_count
	
	return heat_wave_dates, heat_wave_counts
    
    
# Import weather stations and 90th percentile
weather_data_CO, weather_percentile_CO = import_weather_stations('CO')
weather_data_SE, weather_percentile_SE = import_weather_stations('SE')

# Identify heat wave dates
hw_dates_CO, hw_counts_CO = identify_heat_waves(weather_data_CO, weather_percentile_CO)
hw_dates_SE, hw_counts_SE = identify_heat_waves(weather_data_SE, weather_percentile_SE)

# Convert to DataFrame 
hw_df_CO = pd.DataFrame(list(hw_counts_CO.items()), columns=['Year', 'Heat_Waves'])
hw_df_SE = pd.DataFrame(list(hw_counts_SE.items()), columns=['Year', 'Heat_Waves'])

# Calculate trend line
slope_CO, intercept_CO, r_value_CO, p_value_CO, std_err_CO = linregress(hw_df_CO['Year'], hw_df_CO['Heat_Waves'])
slope_SE, intercept_SE, r_value_SE, p_value_SE, std_err_SE = linregress(hw_df_SE['Year'], hw_df_SE['Heat_Waves'])

print(weather_data_CO)
print()
print(weather_percentile_CO)
print()
print(hw_dates_CO)
print()
print(hw_counts_CO)
print()

# Save heat waves episodes date file 
with open('{0}/figs/heat_wave_episodes_CO.txt'.format(path), 'w') as file:
	for start_date, end_date in hw_dates_CO:
		episode = f"Heat wave from {start_date.date()} to {end_date.date()}"
		print(episode)
		file.write(episode + '\n')	
print('Heat wave episodes have been saved')

# Save heat waves episodes date file 
with open('{0}/figs/heat_wave_episodes_SE.txt'.format(path), 'w') as file:
	for start_date, end_date in hw_dates_SE:
		episode = f"Heat wave from {start_date.date()} to {end_date.date()}"
		print(episode)
		file.write(episode + '\n')	
print('Heat wave episodes have been saved')

# Plot figure 
fig = plt.figure()

# Subplot 1
ax=fig.add_subplot(2, 1, 1)
plt.bar(hw_df_CO['Year'], hw_df_CO['Heat_Waves'], color='blue', edgecolor='black', linewidth=1)
plt.plot(hw_df_CO['Year'], intercept_CO + slope_CO * hw_df_CO['Year'], '--r', label=f'Trend line (slope={slope_CO:.2f} p_value < 0.05)')
plt.title('(a)', loc='left', fontweight='bold')
plt.ylabel('Number of Heat Waves', fontweight='bold')
plt.ylim(0,4)
plt.xticks(rotation=45)
plt.grid(True, linestyle='--', alpha=0.8)
plt.setp(ax.get_xticklabels(), visible=False)
plt.legend(loc=2)

# Subplot 2
ax=fig.add_subplot(2, 1, 2)
plt.bar(hw_df_SE['Year'], hw_df_SE['Heat_Waves'], color='blue', edgecolor='black', linewidth=1)
plt.plot(hw_df_SE['Year'], intercept_SE + slope_SE * hw_df_SE['Year'], '--r', label=f'Trend line (slope={slope_SE:.2f} p_value < 0.05)')
plt.title('(b)', loc='left', fontweight='bold')
plt.ylabel('Number of Heat Waves', fontweight='bold')
plt.ylim(0,4)
plt.xticks(rotation=45)
plt.grid(True, linestyle='--', alpha=0.8)
plt.legend(loc=2)

# Path out to save figure
path_out = '{0}/figs'.format(path)
name_out = 'pyplt_regions_ep_hw_SON-2023.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
