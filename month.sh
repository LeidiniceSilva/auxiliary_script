#!/bin/bash

# gerar arquivo month .nc 
# Leidinice Silva / FUNCEME

for AAAA in `seq 2003 2015` ; do

	cd $AAAA

		for m in `seq -w 1 12` ; do
			
			echo
			echo $AAAA $m
			echo
			ls prec_${AAAA}${m}??.nc

			cdo cat prec_${AAAA}${m}??.nc prec_${AAAA}${m}.nc

			cdo timmean prec_${AAAA}${m}.nc prec_mon_${AAAA}${m}.nc

		done
	

	cd ..


done






