#!/bin/bash

#__author__      = 'Leidinice Silva'
#__email__       = 'leidinicesilva@gmail.com'
#__date__        = 'Mar 01, 2023'
#__description__ = 'Posprocessing the CMIP6 models'

echo
echo "--------------- INIT POSPROCESSING CMIP6 MODELS ----------------"

# Models list
model_list=( 'ACCESS-CM2' 'BCC-CSM2-MR' 'CanESM5' 'CMCC-ESM2' 'CNRM-CM6-1' 'CNRM-ESM2-1' 'GFDL-ESM4' 'INM-CM4-8' 'INM-CM5-0' 'KIOST-ESM' 'MIROC6' 'MIROC-ES2L' 'MPI-ESM1-2-HR' 'MPI-ESM1-2-LR' 'MRI-ESM2-0' 'NESM3' 'NorESM2-MM' )

# Variables list
var_list=('pr' 'tas')     

for var in ${var_list[@]}; do

	echo
	echo ${var}

    for model in ${model_list[@]}; do

		path="/home/nice/Documentos/paper_nete/cmip6/"
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
		elif [ ${model} == 'KIOST-ESM' ]
		then
		member='r1i1p1f1_gr1'
		elif [ ${model} == 'MIROC-ES2L' ]
		then
		member='r1i1p1f2_gn'
		else
		member='r1i1p1f1_gn'
		fi
		
		# Datetime
		dt='185001-201412'
		
		echo
		echo "1. Set domain"
		cdo sellonlatbox,260,340,-50,10 ${var}_Amon_${model}_${exp}_${member}_${dt}.nc ${var}_SA_Amon_${model}_${exp}_${member}_${dt}.nc

		echo 
		echo "2. Deleting file"
		rm ${var}_Amon_${model}_${exp}_${member}_${dt}.nc
	
	done
done

echo
echo "--------------- INIT POSPROCESSING CMIP6 MODELS ----------------"
