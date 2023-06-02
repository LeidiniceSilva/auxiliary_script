#!/bin/bash

#__author__      = 'Leidinice Silva'
#__email__       = 'leidinicesilva@gmail.com'
#__date__        = '05/31/2023'
#__description__ = 'Download TRMM satelite data'

for YEAR in 2005 2006; do
    for MON in 01 02 03 04 05 06 07 08 09 10 11 12; do
		for DAY in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31; do
			for HOUR in 00 03 06 09 12 15 18 21; do
			
				PATH="/home/nice/Downloads/trmm/"${YEAR}    	    
				cd ${PATH}
				
				echo
				echo "1. Download file: trmm_3B42_${YEAR}${MON}${DAY}.${HOUR}.7.nc"
				/usr/bin/wget -N http://clima-dods.ictp.it/CHyM_Data/museo/PREC/TRMM/${YEAR}/trmm_3B42_${YEAR}${MON}${DAY}.${HOUR}.7.nc

				echo
				echo "2. Select South America domain"
				/usr/bin/cdo -sellonlatbox,275,330,-40,15 trmm_3B42_${YEAR}${MON}${DAY}.${HOUR}.7.nc trmm_3B42_SA_${YEAR}${MON}${DAY}.${HOUR}.7.nc
				
				/usr/bin/rm trmm_3B42_${YEAR}${MON}${DAY}.${HOUR}.7.nc
				
			done
        done
    done
done
