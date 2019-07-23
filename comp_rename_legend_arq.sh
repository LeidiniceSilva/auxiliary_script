#!/bin/bash

cd /home/leidinice/Documentos/esi/code/esi/12WK/SA/AMJ

for ARQ in `ls -1 *.png`; do # Pega todas as figuras da pasta

     echo ${ARQ}
     P1=$(echo ${ARQ} | awk -F '_' '{ print $6 }')

     echo ${P1}                
     mv -v $ARQ ESI_Regiao_SA_12WK_AMJ_${P1}.png

done

#cd /home/junior/funceme/tsub/data/godas/fig/atsub_mean

#for ARQ in `ls -1 *.png`; do # Pega todas as figuras da pasta
#     
#     P1=$(echo $ARQ | awk -F '_' '{ print $1 }')
#     P2=$(echo $ARQ | awk -F '_' '{ print $2 }')
#     
#     echo ${ARQ} "->" ${P1}_${P2}_PT_ATSUB_AVG4PENT_PCF.png
#     mv -v $ARQ ${P1}_${P2}_PT_ATSUB_AVG4PENT_PCF.png
#     
#  done
  

#cd /home/junior/funceme/tsub/data/godas/fig/hovmoller

#for ARQ in `ls -1 *.png`; do # Pega todas as figuras da pasta
#     
#     P1=$(echo $ARQ | awk -F '_' '{ print $1 }')
#     P2=$(echo $ARQ | awk -F '_' '{ print $2 }')
#     
#     echo ${ARQ} "->" ${P1}_${P2}_PT_ATSUB_EVOL_PCF.png
#     mv -v $ARQ ${P1}_${P2}_PT_ATSUB_EVOL_PCF.png
#     
#  done
#  


