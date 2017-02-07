# -*- coding: utf-8 -*-

# Libs ---------------------------------------------
import numpy as np
from netCDF4 import Dataset
# --------------------------------------------------


"""
Exemplo de como criar um netcdf usando a função Dataset do módulo netCDF4

As dimensões da matriz serão (time, lat, lon)

"""


""" Criando uma matriz qualquer para salvar em um netcdf """

dimensoes_da_matriz = (2, 5, 7)

matriz = np.full(dimensoes_da_matriz, np.nan)

matriz[0, 1, 1] = 157
matriz[0, 2, 2] = 99


""" Vamos definir latitudes e longitudes quaisquer para o exemplo """

latitudes  = np.array( [ -7,  -6,  -4,  -3,  -2] )
longitudes = np.array( [-43, -42, -41, -40, -39, -38, -37] )


""" Salvando o netcdf da matriz do exemplo """

# Nome do arquivo de saída
nc_output_name = 'netcdf_criado_como_exemplo.nc'

# Nome da variável a ser salva
variavel = 'var_exemplo'

# Criando o netcdf

foo = Dataset(nc_output_name, 'w', format='NETCDF3_CLASSIC')

foo.createDimension('time', None)
foo.createDimension('latitude', matriz.shape[1])
foo.createDimension('longitude', matriz.shape[2])

T                   = foo.createVariable('time', float, ('time'))
T.units             = 'months since 1982-01-15 00:00:00'
T.calendar          = 'standard'
T.standard_name     = 'time'
T.axis              = 'T'
T[:]                = range(matriz.shape[0])

Y                   = foo.createVariable('latitude', float, ('latitude'))
Y.axis              = 'Y'
Y.units             = 'degrees_north'
Y.long_name         = 'latitude'
Y[:]                = latitudes

X                   = foo.createVariable('longitude', float, ('longitude'))
X.axis              = 'X'
X.units             = 'degrees_east'
X.long_name         = 'longitude'
X[:]                = longitudes

VAR                 = foo.createVariable(variavel, float, ('time', 'latitude', 'longitude'))
VAR.units           = 'mm'
VAR.standard_name   = 'precipitation'
VAR.missing_value   = np.nan
VAR[:]              = matriz[:]

foo.close()

