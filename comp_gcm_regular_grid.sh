
#this script is useful to:
#regular grid using bilinear method
#change the calendar - from 360 days to 365 days otherwise it will not open on grands

#!/bin/bash

############################################ Historial ########################################

#historical MOHC-HadGEM2-ES
datadir=/home/pesquisa/Documentos/Analises_teste/Tratados/South_America_HadGEM_MENSAIS/Historical/

#historical MPI-M-MPI-ESM-MR
#datadir=/home/pesquisa/Documentos/Analises_teste/Tratados/South_America_MPI_MENSAIS/Historical/

#historical NCC-NorESM1-M
#datadir=/home/pesquisa/Documentos/Analises_teste/Tratados/South_America_NorESM_MENSAIS/Historical/


################################################### RCP2.6 ##############################################

#RCP2.6 MOHC-HadGEM2-ES
#datadir=/home/pesquisa/Documentos/Analises_teste/Tratados/South_America_HadGEM_MENSAIS/RCP2.6/

#RCP2.6 MPI-M-MPI-ESM-MR
#datadir=/home/pesquisa/Documentos/Analises_teste/Tratados/South_America_MPI_MENSAIS/RCP2.6/

#RCP2.6 NCC-NorESM1-M
#datadir=/home/pesquisa/Documentos/Analises_teste/Tratados/South_America_NorESM_MENSAIS/RCP2.6/

################################################ RCP8.5 ######################################

#RCP8.5
#datadir=/home/pesquisa/Documentos/Analises_teste/Tratados/South_America_HadGEM_MENSAIS/RCP8.5/

#RCP8.5 MPI-M-MPI-ESM-MR
#datadir=/home/pesquisa/Documentos/Analises_teste/Tratados/South_America_MPI_MENSAIS/RCP8.5/

#RCP8.5 NCC-NorESM1-MS
#datadir=/home/pesquisa/Documentos/Analises_teste/Tratados/South_America_NorESM_MENSAIS/RCP8.5/

var=pr_
model=MOHC-HadGEM2-ES
#model=MPI-M-MPI-ESM-MR
#model=NCC-NorESM1-M

period=19700101-20060101
#period=20060101-21000101

scenario=historical
#scenario=rcp26
#scenario=rcp85
domain=SAM-22



#################################
### setting calendar 

cdo setcalendar,standard  ${datadir}${var}${domain}_${model}_${scenario}_r1i1p1_ICTP-RegCM4-7_v0_mon_${period}.nc  ${var}${domain}_${model}_${scenario}_r1i1p1_ICTP-RegCM4-7_v0_mon_${period}.ncn 

mv ${var}${domain}_${model}_${scenario}_r1i1p1_ICTP-RegCM4-7_v0_mon_${period}.ncn  ${var}${domain}_${model}_${scenario}_r1i1p1_ICTP-RegCM4-7_v0_mon_${period}.nc

#### regular grid (22km). 

./regrid.sh  ${var}${domain}_${model}_${scenario}_r1i1p1_ICTP-RegCM4-7_v0_mon_${period}.nc -59.38107,18.97321,0.22 -111.9328,-6.067216,0.22  bil


