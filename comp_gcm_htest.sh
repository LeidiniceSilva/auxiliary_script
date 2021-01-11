#!/bin/bash

###Faz o teste de hipotese para medias considerando que ambos os conjuntos possuam o mesmo numero de dados
#you should use the output from ensemble


#Near future RCP2.6
#modelRCP85=/home/pesquisa/Documentos/Analises_teste/ensemble/ensemble_rcp26_annual2041-2060.nc
#Far future RCP2.6
#modelRCP85=/home/pesquisa/Documentos/Analises_teste/ensemble/ensemble_rcp26_annual2080-2099.nc


#Near future RCP8.5
#modelRCP85=/home/pesquisa/Documentos/Analises_teste/ensemble/ensemble_rcp85_annual2041-2060.nc
#Far future RCP8.5
modelRCP85=/home/pesquisa/Documentos/Analises_teste/ensemble/ensemble_rcp85_annual2080-2099.nc


#historical RCP2.6
#modelpres=/home/pesquisa/Documentos/Analises_teste/ensemble/ensemble_rcp26_annual1995-2014.nc
#historical RCP8.5
modelpres=/home/pesquisa/Documentos/Analises_teste/ensemble/ensemble_rcp85_annual1995-2014.nc

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
'enable print teste.gmf'
'set display color white'
'c'
'set mpdset mres'
'set ylint 10'
'set xlint 20'
'set grads off'
'set gxout shaded'

'set xlopts 1 8 0.3'
'set ylopts 1 8 0.3'  
*'set lon -88 -25'
*'set lat -58 10'
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
*'set lon -88 -33'
*'set lat -57 15'
'd dif'
'cbarn 1 1  9.7 4.1'
*'cbarn.gs'
*Caixa America Sul

'drawbox -88 -30 -56 13'

*caixa AMZ
'drawbox -65 -50 -15 0'
*Caixa LPB
'drawbox -63 -48.9 -32.5 -20'
*Caixa NDE
'drawbox -46 -35 -15 -3'
*Caixa AND
'drawbox -75 -68 -18 -13'
*Caixa PATAGONIA
'drawbox -71 -66 -50 -35'


***** NESSA ETAPA E DEFINIDA O VALOR DE T SIGNIFICANTE (TABELA TESTT, CONSULTE O VALOR DE ACORDO COM SUA AMOSTRA)
***** O VALOR PADRAO ESTA 1.96	
***** O t calculado é significativo onde os valores sao superiores a 1.96	

*****95%
'hatch.gs t 1.96 1000' 

*****************************************

'draw title EnsembleRCP8.5 (Far Future)'


#############################################################

'printim diff_rcp85_farfuture.png'

'print'
'disable print'
'!gxeps -c -i teste.gmf -o diff_rcp85_farfuture.eps'
'!/bin/rm teste.gmf'


EOF

grads -l -c Plot_sig.gs











