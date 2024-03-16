#!/bin/bash

#__author__      = 'Leidinice Silva'
#__email__       = 'leidinicesilva@gmail.com'
#__date__        = 'Mar 04, 2024'
#__description__ = 'Preprocessing the OBS datasets with CDO'

var_list='msl q t2m u v z'

for var in ${var_list[@]}; do

	# Temperature
	cdo -timmean -selseas,SON ${var}_era5_mon_1991-2021.nc ${var}_era5_son_1991-2021.nc

	cdo seldate,2023-09-17,2023-09-27 ${var}_era5_day_2023.nc ${var}_era5_hw_1_17-27Sep2023.nc
	cdo seldate,2023-11-11,2023-11-18 ${var}_era5_day_2023.nc ${var}_era5_hw_2_11-18Nov2023.nc
	
	cdo timmean ${var}_era5_hw_1_17-27Sep2023.nc ${var}_era5_hw_1_17-27Sep2023_mean.nc
	cdo timmean ${var}_era5_hw_2_11-18Nov2023.nc ${var}_era5_hw_2_11-18Nov2023_mean.nc

	cdo -b F32 sub ${var}_era5_hw_1_17-27Sep2023_mean.nc ${var}_era5_son_1991-2021.nc ${var}_era5_hw_1_17-27Sep2023_anom.nc
	cdo -b F32 sub ${var}_era5_hw_2_11-18Nov2023_mean.nc ${var}_era5_son_1991-2021.nc ${var}_era5_hw_2_11-18Nov2023_anom.nc

done
