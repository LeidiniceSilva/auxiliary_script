#!/bin/bash

## you will need the output data from averages.sh

dir=/home/pesquisa/Documentos/Analises_teste/averages/
var=pr_
#model=MOHC-HadGEM2-ES
#model=MPI-M-MPI-ESM-MR
model=NCC-NorESM1-M


scenario=rcp26
#scenario=rcp85
domain=SAM-22

allperiod=1970-2100
season1=annual
season2=djf
season3=jja

timeslice1=1995-2014
timeslice2=2041-2060
timeslice3=2080-2099

######################### separating for periods - historical (1995-2014), near future (2041-2060) and far future (2080-2099) ###################

### RCP2.6

#annual

cdo -seldate,1995-01-01,2014-12-31 ${dir}${var}${domain}_${model}_${scenario}_r1i1p1_ICTP-RegCM4-7_v0_${allperiod}_${season1}.nc ${var}${domain}_${model}_${scenario}_r1i1p1_ICTP-RegCM4-7_v0_${timeslice1}_${season1}.nc 

cdo -seldate,2041-01-01,2060-12-31 ${dir}${var}${domain}_${model}_${scenario}_r1i1p1_ICTP-RegCM4-7_v0_${allperiod}_${season1}.nc ${var}${domain}_${model}_${scenario}_r1i1p1_ICTP-RegCM4-7_v0_${timeslice2}_${season1}.nc

cdo -seldate,2080-01-01,2099-12-31 ${dir}${var}${domain}_${model}_${scenario}_r1i1p1_ICTP-RegCM4-7_v0_${allperiod}_${season1}.nc ${var}${domain}_${model}_${scenario}_r1i1p1_ICTP-RegCM4-7_v0_${timeslice3}_${season1}.nc

#djf

cdo -seldate,1995-01-01,2014-12-31 ${dir}${var}${domain}_${model}_${scenario}_r1i1p1_ICTP-RegCM4-7_v0_${allperiod}_${season2}.nc ${var}${domain}_${model}_${scenario}_r1i1p1_ICTP-RegCM4-7_v0_${timeslice1}_${season2}.nc 

cdo -seldate,2041-01-01,2060-12-31 ${dir}${var}${domain}_${model}_${scenario}_r1i1p1_ICTP-RegCM4-7_v0_${allperiod}_${season2}.nc ${var}${domain}_${model}_${scenario}_r1i1p1_ICTP-RegCM4-7_v0_${timeslice2}_${season2}.nc


cdo -seldate,2080-01-01,2099-12-31 ${dir}${var}${domain}_${model}_${scenario}_r1i1p1_ICTP-RegCM4-7_v0_${allperiod}_${season2}.nc ${var}${domain}_${model}_${scenario}_r1i1p1_ICTP-RegCM4-7_v0_${timeslice3}_${season2}.nc


#jja

cdo -seldate,1995-01-01,2014-12-31 ${dir}${var}${domain}_${model}_${scenario}_r1i1p1_ICTP-RegCM4-7_v0_${allperiod}_${season3}.nc ${var}${domain}_${model}_${scenario}_r1i1p1_ICTP-RegCM4-7_v0_${timeslice1}_${season3}.nc 

cdo -seldate,2041-01-01,2060-12-31 ${dir}${var}${domain}_${model}_${scenario}_r1i1p1_ICTP-RegCM4-7_v0_${allperiod}_${season3}.nc ${var}${domain}_${model}_${scenario}_r1i1p1_ICTP-RegCM4-7_v0_${timeslice2}_${season3}.nc


cdo -seldate,2080-01-01,2099-12-31 ${dir}${var}${domain}_${model}_${scenario}_r1i1p1_ICTP-RegCM4-7_v0_${allperiod}_${season3}.nc ${var}${domain}_${model}_${scenario}_r1i1p1_ICTP-RegCM4-7_v0_${timeslice3}_${season3}.nc



