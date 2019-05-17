#!/bin/bash

#Nesse script feita a interpolacao da grade do dado1 para a grade do dado2
#o resultado eh alocado na variavel dado_int.nc

dado1='model.nc'
dado2='obs.nc'

cdo genbil,dado2 dado1 weight.nc # gerando os pesos
cdo remap,dados1,weight.nc dado2 dada_int.nc

##############################################################################

# Caso o dado a ser interpolado tenha dominio diferente eh interssante recortar em um regiao que os dois sejam proximos
# selecione a regiao no dado a ser interpolado

se ele compreender as lon2,lon1 e latitudes lat2,lat1

cdo sellonlatbox,lon2,lon1,lat2,lat1 dado2.nc dados2_cutted.nc
