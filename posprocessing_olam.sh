#!/bin/bash

echo
echo "--------------- INICIO ----------------"


rm precip_controle_1982_2012_chen_neb_new_-6hours.nc
rm precip_controle_1982_2012_chen_neb_new_1982-01-01-00_2012-12-30-12.nc
rm precip_controle_1982_2012_chen_neb_new_-6hours_1982-01-01-00_2012-12-30-12.nc
rm precip_controle_1982_2012_chen_neb_new_REAL.nc
rm precip_controle_1982_2012_chen_neb_new_REAL_ok.nc
rm precip_controle_1982_2012_chen_neb_new_1982-01-01-00.nc
rm precip_controle_1982_2012_chen_neb_new_REAL_ok_full.nc
rm precip_controle_1982_2012_chen_neb_new_REAL_ok_full_monsum.nc
rm mask_positivos.nc
rm precip_controle_1982_2012_chen_neb_new_REAL_ok_full_negcor.nc
rm precip_controle_1982_2012_chen_neb_new_REAL_ok_full_negcor_monsum.nc


echo 
echo "1 - shifttime (-6hours)"
cdo shifttime,-6hours \
precip_controle_1982_2012_chen_g1_neb.nc \
precip_controle_1982_2012_chen_neb_new_-6hours.nc


echo 
echo "2 - seldate"
cdo seldate,1982-01-01T00:00,2012-12-30T12:00 \
precip_controle_1982_2012_chen_g1_neb.nc \
precip_controle_1982_2012_chen_neb_new_1982-01-01-00_2012-12-30-12.nc


echo 
echo "3 - seldate"
cdo seldate,1982-01-01T00:00,2012-12-30T12:00 \
precip_controle_1982_2012_chen_neb_new_-6hours.nc \
precip_controle_1982_2012_chen_neb_new_-6hours_1982-01-01-00_2012-12-30-12.nc


echo 
echo "4 - sub"
cdo sub \
precip_controle_1982_2012_chen_neb_new_-6hours_1982-01-01-00_2012-12-30-12.nc \
precip_controle_1982_2012_chen_neb_new_1982-01-01-00_2012-12-30-12.nc \
precip_controle_1982_2012_chen_neb_new_REAL.nc


echo 
echo "5 - shifttime (+6hours)"
cdo shifttime,+6hours \
precip_controle_1982_2012_chen_neb_new_REAL.nc \
precip_controle_1982_2012_chen_neb_new_REAL_ok.nc


echo 
echo "6 - seldate"
cdo seldate,1982-01-01T00:00,1982-01-01T4:00 \
precip_controle_1982_2012_chen_g1_neb.nc \
precip_controle_1982_2012_chen_neb_new_1982-01-01-00.nc


echo 
echo "7 - cat"
cdo cat \
precip_controle_1982_2012_chen_neb_new_1982-01-01-00.nc \
precip_controle_1982_2012_chen_neb_new_REAL_ok.nc \
precip_controle_1982_2012_chen_neb_new_REAL_ok_full.nc


echo 
echo "8 - gec"
cdo gec,0 \
precip_controle_1982_2012_chen_neb_new_REAL_ok_full.nc \
mask_positivos.nc


echo 
echo "9 - gec x full"
cdo mul \
precip_controle_1982_2012_chen_neb_new_REAL_ok_full.nc \
mask_positivos.nc \
precip_controle_1982_2012_chen_neb_new_REAL_ok_full_negcor.nc


echo 
echo "10 - monsum"
cdo monsum \
precip_controle_1982_2012_chen_neb_new_REAL_ok_full.nc \
precip_controle_1982_2012_chen_neb_new_REAL_ok_full_monsum.nc


echo 
echo "11 - monsum (gec x full)"
cdo monsum \
precip_controle_1982_2012_chen_neb_new_REAL_ok_full_negcor.nc \
precip_controle_1982_2012_chen_neb_new_REAL_ok_full_negcor_monsum.nc


echo
echo "--------------- FIM ----------------"
