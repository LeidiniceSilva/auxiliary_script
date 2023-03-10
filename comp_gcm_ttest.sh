#!/bin/bash

####you should use the output from ensemble - 20 time steps for each time slice

#Near future RCP8.5 (don't need to uncommented)
area=NearFuture85
modelRCP85=/media/christiana/Seagate/MLLOPART/pr/pr_AFR-22_rcp85_r1i1p1_ICTP-RegCM4-7_v0_mon_20410101-20601231_ensemble_HadGEM_NorESM.nc

#Far future RCP8.5 (don't need to uncommented)
#area=FarFuture85
#modelRCP85=/home/pesquisa/Documentos/Analises_teste/global_data/EUR_22/ensemble/pr_Amon_ensemble_1970-2100_r1i1p1_rcp85_lonlat_annual_2080-2099.nc

#historical RCP8.5 (don't need to uncommented)
modelpres=/media/christiana/Seagate/MLLOPART/pr/pr_AFR-22_rcp85_r1i1p1_ICTP-RegCM4-7_v0_mon_19950101-20141231_ensemble_HadGEM_NorESM.nc

###definindo o numero de tempos do arquivo
cdo mulc,0 ${modelRCP85} num.nc
cdo addc,1 num.nc num1.nc
cdo timsum num1.nc 00cont.nc
cdo chname,pr,n 00cont.nc ndata_GCM.nc

####calcula a media
cdo timmean ${modelRCP85} medfut.nc
cdo timmean ${modelpres}  medpres.nc

####calcula a diferença futuro menos o presente
cdo sub medfut.nc medpres.nc  0diff_GCM.nc

####calcula o desvio-padrao
cdo timstd ${modelRCP85} 0std_GCM_fut.nc
cdo timstd ${modelpres} 0std_GCM_pres.nc

#### eleva o desvio-padrao ao quadrado
cdo mul 0std_GCM_fut.nc 0std_GCM_fut.nc std2_GCM_fut.nc
cdo mul 0std_GCM_pres.nc 0std_GCM_pres.nc std2_GCM_pres.nc

#### altera o nome da variavel para os calculos futuros
cdo chname,pr,stdfut std2_GCM_fut.nc std_GCM_fut.nc
cdo chname,pr,stdpres std2_GCM_pres.nc std_GCM_pres.nc
cdo chname,pr,dif 0diff_GCM.nc diff_GCM.nc

#### junta os calculos em um unico nc (a diff da media, std^2, e o n)
cdo merge diff_GCM.nc std_GCM_fut.nc std_GCM_pres.nc ndata_GCM.nc 0ttest_GCM.nc
cdo expr,'t=abs(dif/sqrt((stdfut/n)+(stdpres/n)))' 0ttest_GCM.nc ttest_GCM.nc
###Esse eh a saida que nos interessa que tem a diferenca e o teste
cdo merge diff_GCM.nc ttest_GCM.nc sig_GCM.nc

echo GCM
cat <<EOF> Plot_sig.gs
'reinit'
say GCM
*ESCOLHA PORC OU DIFF
"sdfopen sig_GCM.nc"

*** UTILIZAR PARA PLOT PARA BR OU AS
*'enable print teste.gmf'
'set display color white'
'c'
'set mpdset mres'
'set ylint 10'
'set xlint 20'
'set grads off'
'set gxout shaded'

*'set xlopts 1 8 0.3'
*'set ylopts 1 8 0.3'
'set lon -41.54 65.54'
'set lat -49.92 47.27'
'run mycolors'
'set rgb 90 244 245 252'
'set rgb 91 221 221 253'
'set rgb 92 160 215 160'
'set rgb 93 119 245 113'
'set rgb 94 250 231 120'
'set rgb 95 253 94   0'
'set rgb 96 165 0    0'
'mycolors.gs'
'set clevs -2 -1.5 -1.0 -0.5  0  0.5 1.0 1.5 2 '
*'set clevs  -75 -50 -25 -5  0 5 25 50 75 '
'set ccols  40 41 42 43 44   50 51 52 53 54  '
*'set lon -41.54 65.54'
*'set lat -49.92 47.27'
'd dif'
*'cbarn 1 1  9.7 4.1'
'cbarn.gs'

*Southern African
*Southern African

'drawbox 18 30 -30 -10'

*West Africa subzone
'drawbox -15 4 6 16'

*Africa subzone
'drawbox -20 52 -35 38'

***** NESSA ETAPA E DEFINIDA O VALOR DE T SIGNIFICANTE (TABELA TESTT, CONSULTE O VALOR DE ACORDO COM SUA AMOSTRA)
***** O VALOR PADRAO ESTA 1.96
***** O t calculado é significativo onde os valores sao superiores a 1.96

*****95%
'hatch.gs t 1.96 1000'

*****************************************

'draw title RCMEns $area AFR-22'


#############################################################
'printim $area.png'

'gxprint $area.pdf'
'disable gxprint'
'!gxeps -c -i teste.gmf -o $area.eps'
'!/bin/rm teste.gmf'


EOF

grads -l -c Plot_sig.gs





