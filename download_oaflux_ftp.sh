#!/bin/bash

# Download via ftp --------------------- evapr
path="/home/leidinice/documentos/oaflux/dados/evapr/"
echo
cd ${path}
for year in `seq 1981 2016`; do
    /usr/bin/wget -N ftp://ftp.whoi.edu/pub/science/oaflux/data_v3/monthly/evaporation/evapr_oaflux_${year}.nc.gz
    gunzip evapr_oaflux_${year}.nc.gz
done
/usr/bin/wget -N ftp://ftp.whoi.edu/pub/science/oaflux/data_v3/monthly/evaporation/evapr_oaflux_2017rt.nc.gz
gunzip evapr_oaflux_2017rt.nc.gz



# Download via ftp ----------------------- lh
path="/home/leidinice/documentos/oaflux/dados/lh/"
echo
cd ${path}
for year in `seq 1981 2016`; do
    /usr/bin/wget -N ftp://ftp.whoi.edu/pub/science/oaflux/data_v3/monthly/turbulence/lh_oaflux_${year}.nc.gz
    gunzip lh_oaflux_${year}.nc.gz
done
/usr/bin/wget -N ftp://ftp.whoi.edu/pub/science/oaflux/data_v3/monthly/turbulence/lh_oaflux_2017rt.nc.gz
gunzip lh_oaflux_2017rt.nc.gz



# Download via ftp ----------------------- qa
path="/home/leidinice/documentos/oaflux/dados/qa/"
echo
cd ${path}
for year in `seq 1981 2016`; do
    /usr/bin/wget -N ftp://ftp.whoi.edu/pub/science/oaflux/data_v3/monthly/turbulence/qa_oaflux_${year}.nc.gz
    gunzip qa_oaflux_${year}.nc.gz
done
/usr/bin/wget -N ftp://ftp.whoi.edu/pub/science/oaflux/data_v3/monthly/turbulence/qa_oaflux_2017rt.nc.gz
gunzip qa_oaflux_2017rt.nc.gz



# Download via ftp ----------------------- sh
path="/home/leidinice/documentos/oaflux/dados/sh/"
echo
cd ${path}
for year in `seq 2010 2016`; do
    /usr/bin/wget -N ftp://ftp.whoi.edu/pub/science/oaflux/data_v3/monthly/turbulence/sh_oaflux_${year}.nc.gz
    gunzip sh_oaflux_${year}.nc.gz
done
/usr/bin/wget -N ftp://ftp.whoi.edu/pub/science/oaflux/data_v3/monthly/turbulence/sh_oaflux_2017rt.nc.gz
gunzip sh_oaflux_2017rt.nc.gz



# Download via ftp ----------------------- ts
path="/home/leidinice/documentos/oaflux/dados/ts/"
echo
cd ${path}
for year in `seq 1981 2016`; do
    /usr/bin/wget -N ftp://ftp.whoi.edu/pub/science/oaflux/data_v3/monthly/turbulence/ts_oaflux_${year}.nc.gz
    gunzip ts_oaflux_${year}.nc.gz
done
/usr/bin/wget -N ftp://ftp.whoi.edu/pub/science/oaflux/data_v3/monthly/turbulence/ts_oaflux_2017rt.nc.gz
gunzip ts_oaflux_2017rt.nc.gz



# Download via ftp ---------------------- ws
path="/home/leidinice/documentos/oaflux/dados/ws/"
echo
cd ${path}
for year in `seq 1981 2016`; do
    /usr/bin/wget -N ftp://ftp.whoi.edu/pub/science/oaflux/data_v3/monthly/turbulence/ws_oaflux_${year}.nc.gz
    gunzip ws_oaflux_${year}.nc.gz
done
/usr/bin/wget -N ftp://ftp.whoi.edu/pub/science/oaflux/data_v3/monthly/turbulence/ws_oaflux_2017rt.nc.gz
gunzip ws_oaflux_2017rt.nc.gz


