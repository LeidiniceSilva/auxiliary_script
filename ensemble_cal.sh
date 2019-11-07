#!/bin/bash

#####calculating ensemble we should use the output from timeslices

dir=/home/pesquisa/Documentos/Analises_teste/timeslices/

model1=MOHC-HadGEM2-ES
model2=MPI-M-MPI-ESM-MR
model3=NCC-NorESM1-M

var=pr_

#scenario=rcp26
scenario=rcp85

#timeslice=1995-2014
#timeslice=2041-2060
timeslice=2080-2099

#season=annual
#season=djf
season=jja

domain=SAM-22


cdo add ${dir}${var}${domain}_${model1}_${scenario}_r1i1p1_ICTP-RegCM4-7_v0_${timeslice}_${season}.nc ${dir}${var}${domain}_${model2}_${scenario}_r1i1p1_ICTP-RegCM4-7_v0_${timeslice}_${season}.nc ensemble1.nc

cdo add ${dir}${var}${domain}_${model3}_${scenario}_r1i1p1_ICTP-RegCM4-7_v0_${timeslice}_${season}.nc ensemble1.nc ensemble_${scenario}_${season}${timeslice}.nc

rm -rf ensemble1.nc

cdo divc,3 ensemble_${scenario}_${season}${timeslice}.nc ensemble_${scenario}_${season}${timeslice}.ncn

mv ensemble_${scenario}_${season}${timeslice}.ncn ensemble_${scenario}_${season}${timeslice}.nc

cdo mulc,86400 ensemble_${scenario}_${season}${timeslice}.nc ensemble_${scenario}_${season}${timeslice}.ncn

mv ensemble_${scenario}_${season}${timeslice}.ncn ensemble_${scenario}_${season}${timeslice}.nc






