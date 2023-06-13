#!/bin/bash

#__author__      = 'Leidinice Silva'
#__email__       = 'leidinicesilva@gmail.com'
#__date__        = 'Mar 01, 2023'
#__description__ = 'Posprocessing the CMIP6 models'

echo
echo "--------------- INIT POSPROCESSING CMIP6 MODELS ----------------"

# Models list
model_list=( 'ACCESS-CM2' 'AWI-CM-1-1-MR' 'BCC-CSM2-MR' 'CAMS-CSM1-0' 'CanESM5-CanOE' 'CESM2' 'CMCC-CM2-SR5' 
'CMCC-ESM2' 'CNRM-CM6-1' 'CNRM-CM6-1-HR' 'CNRM-ESM2-1' 'EC-Earth3-CC' 'EC-Earth3-Veg-LR' 'FGOALS-f3-L' 'FGOALS-g3' 
'FIO-ESM-2-0' 'GFDL-ESM4' 'HadGEM3-GC31-LL' 'IITM-ESM' 'INM-CM4-8' 'INM-CM5-0' 'IPSL-CM6A-LR' 'KACE-1-0-G' 
'MCM-UA-1-0' 'MIROC6' 'MIROC-ES2L' 'MPI-ESM1-2-LR' 'MRI-ESM2-0' 'NESM3' 'NorESM2-MM' 'TaiESM1' 'UKESM1-0-LL')

# Variables list
var_list=( 'pr' )     

for var in ${var_list[@]}; do

	echo
	echo ${var}

    for model in ${model_list[@]}; do

		path="/home/nice/Downloads/"
		cd ${path}
		
		echo
		echo ${model}
		
		# Experiment name
		exp='historical'

		# Member name
		if [ ${model} == 'CNRM-CM6-1' ]
		then
		member='r1i1p1f2_gr'
		elif [ ${model} == 'CNRM-ESM2-1' ]
		then
		member='r1i1p1f2_gr'
		elif [ ${model} == 'GFDL-ESM4' ]
		then
		member='r1i1p1f1_gr1'
		elif [ ${model} == 'INM-CM4-8' ]
		then
		member='r1i1p1f1_gr1'
		elif [ ${model} == 'INM-CM5-0' ]
		then
		member='r1i1p1f1_gr1'
		elif [ ${model} == 'MIROC-ES2L' ]
		then
		member='r1i1p1f2_gn'
		else
		member='r1i1p1f1_gn'
		fi
		
		# Datetime
		dt0='185001-201412'
		dt1='198501-201412'
		
		echo
		echo "1. Set domain"
		cdo sellonlatbox,260,340,-50,10 ${var}_Amon_${model}_${exp}_${member}_${dt0}.nc ${var}_SA_Amon_${model}_${exp}_${member}_${dt0}.nc

		echo
		echo "2. Select date"
		cdo seldate,1985-01-01T00:00:00,2014-12-31T00:00:00 ${var}_SA_Amon_${model}_${exp}_${member}_${dt0}.nc ${var}_SA_Amon_${model}_${exp}_${member}_${dt1}.nc

		echo 
		echo "3. Deleting file"
		rm ${var}_Amon_${model}_${exp}_${member}_${dt0}.nc
		rm ${var}_SA_Amon_${model}_${exp}_${member}_${dt0}.nc
	
	done
done

echo
echo "--------------- INIT POSPROCESSING CMIP6 MODELS ----------------"
