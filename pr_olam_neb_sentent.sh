#!/bin/bash
# gerar arquivo .nc apartir do .ctl 
# Leidinice Silva / FUNCEME

#lats4d.sh -q -format netcdf -mxtimes 1461 -lat -20 -2 -lon -50 -30 -i precip_controle1982harr_g2.ctl -o precip_controle_1982_harr_g2_neb

#cdo sellonlatbox,-44,-33,2,-10 precip_controle_1982_harr_g2_neb.nc precip_controle_1982_harr_g2_neb.nc

#cdo cat precip_controle_????_harr_g2_neb.nc precip_controle_1982_2012_harr_g2_neb.nc

#cdo remapbil,pr_Amon_CRU-TS3.22_observation_198201_201212_neb.nc precip_controle_1982_2012_harr_g2_neb.nc precip_controle_1982_2012_harr_g2_neb_grad_obs.nc

#cdo precip_controle_1982_2012_harr_g2_neb_grad_obs.nc obs_um.nc precip_controle_1982_2012_harr_g2_neb_noocean.nc

#cdo mulc,86400 pr_Amon_CRU-TS3.22_observation_198201_201212_neb.nc pr_Amon_CRU-TS3.22_observation_198201_201212_neb_mmd.nc

#cdo mulc,30.25 pr_Amon_CRU-TS3.22_observation_198201_201212_neb_mmd.nc pr_Amon_CRU-TS3.22_observation_198201_201212_neb_mm/m.nc


