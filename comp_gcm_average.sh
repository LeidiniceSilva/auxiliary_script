#!/bin/bash

## you will need the output data from data.sh

dir=/home/pesquisa/Documentos/Analises_teste/data/
var=pr_
#model=MOHC-HadGEM2-ES
#model=MPI-M-MPI-ESM-MR
model=NCC-NorESM1-M

period1=19700101-20060101
period2=20060101-21000101

scenario1=historical
scenario2=rcp26
scenario3=rcp85
domain=SAM-22

allperiod=1970-2100
season1=annual
season2=djf
season3=jja

timeslice1=1995-2014
timeslice2=2041-2060
timeslice3=2080-2099

######### add historical plus scenario - it is better to make the averages - annual, djf and jja ########

#RCP2.6
cdo mergetime ${dir}${var}${domain}_${model}_${scenario1}_r1i1p1_ICTP-RegCM4-7_v0_mon_${period1}.nc ${dir}${var}${domain}_${model}_${scenario2}_r1i1p1_ICTP-RegCM4-7_v0_mon_${period2}.nc  ${var}${domain}_${model}_${scenario2}_r1i1p1_ICTP-RegCM4-7_v0_${allperiod}.nc 

#RCP8.5
cdo mergetime ${dir}${var}${domain}_${model}_${scenario1}_r1i1p1_ICTP-RegCM4-7_v0_mon_${period1}.nc  ${dir}${var}${domain}_${model}_${scenario3}_r1i1p1_ICTP-RegCM4-7_v0_mon_${period2}.nc  ${var}${domain}_${model}_${scenario3}_r1i1p1_ICTP-RegCM4-7_v0_${allperiod}.nc


#################### annual average #################################

##RCP2.6
cdo -yearavg ${var}${domain}_${model}_${scenario2}_r1i1p1_ICTP-RegCM4-7_v0_${allperiod}.nc ${var}${domain}_${model}_${scenario2}_r1i1p1_ICTP-RegCM4-7_v0_${allperiod}_${season1}.nc

###RCP8.5
cdo -yearavg ${var}${domain}_${model}_${scenario3}_r1i1p1_ICTP-RegCM4-7_v0_${allperiod}.nc ${var}${domain}_${model}_${scenario3}_r1i1p1_ICTP-RegCM4-7_v0_${allperiod}_${season1}.nc


################# djf ##########################################

#RCP2.6
cdo -r -timselavg,3 -selmon,1,2,12 ${var}${domain}_${model}_${scenario2}_r1i1p1_ICTP-RegCM4-7_v0_${allperiod}.nc ${var}${domain}_${model}_${scenario2}_r1i1p1_ICTP-RegCM4-7_v0_${allperiod}_${season2}.nc

###RCP8.5
cdo -r -timselavg,3 -selmon,1,2,12 ${var}${domain}_${model}_${scenario3}_r1i1p1_ICTP-RegCM4-7_v0_${allperiod}.nc ${var}${domain}_${model}_${scenario3}_r1i1p1_ICTP-RegCM4-7_v0_${allperiod}_${season2}.nc

################# jja ##########################################

#RCP2.6
cdo -r -timselavg,3 -selmon,6,7,8 ${var}${domain}_${model}_${scenario2}_r1i1p1_ICTP-RegCM4-7_v0_${allperiod}.nc ${var}${domain}_${model}_${scenario2}_r1i1p1_ICTP-RegCM4-7_v0_${allperiod}_${season3}.nc

###RCP8.5
cdo -r -timselavg,3 -selmon,6,7,8 ${var}${domain}_${model}_${scenario3}_r1i1p1_ICTP-RegCM4-7_v0_${allperiod}.nc ${var}${domain}_${model}_${scenario3}_r1i1p1_ICTP-RegCM4-7_v0_${allperiod}_${season3}.nc



