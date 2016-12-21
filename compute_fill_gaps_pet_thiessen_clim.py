# -*- coding: utf-8 -*-

__author__ = "Leidinice Silva"
__copyright__ = "Copyright 2016, Funceme Hydropy Project"
__credits__ = ["Francisco Vasconcelos Junior", "Marcelo Rodrigues", "Enzo Pinheiro", "Diogenes Fontenele"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marcelo Rodrigues"
__email__ = "leidinice.silva@funceme.br"
__date__ = 8/17/2016

# Import data
import netCDF4
import requests
import calendar
import argparse
import datetime
import numpy as np
import numpy.ma as ma
import os
from os.path import expanduser
from netCDF4 import Dataset
from PyFuncemeClimateTools import Thiessen
from hidropy.preprocessing.write_thiessen import write_thiessen
from hidropy.preprocessing.utils import basin_dict

# Directories of input and output data
home = os.path.expanduser("~")
hidropy_path = "/home/leidinice/Documentos/musf"
path_out = "/home/leidinice/Documentos/musf/results/pet_daily_oper/"

folders = os.listdir("{0}/hidropy/hidropy/shapes/basins/".format(hidropy_path))
basins = sorted(basin_dict(micro=True, basin_name='grande'))  # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< change here!!

for basin in basins:
    basin_fullname = basin_dict(basin)[2]
    macro_name = basin_dict(basin)[1]

pet = np.full((1, 245, 313), np.nan) # Declaring variable  (time step, lat and lon)

ano = range(2015, 2016+1)
mes = range(0, 12)

if not home:
    pth = home+"/home/leidinice/io/inmet/operation/daily/pet/"

pth = "{0}{1}/prec_{1}{2:02d}{3:02d}".format(pth, ano, mes)

print "Acumulado de", mes, " ano:", ano
dias_mes = calendar.monthrange(ano, mes)[1]
petnw = []

for dia in xrange(1, dias_mes + 1):
    lead0_date = datetime.date(ano, mes, dia)
    start_date = lead0_date.strftime('%Y%m%d')
    end_date = lead0_date.strftime('%Y%m%d')
    arq = "{0}{1:02d}.ctl".format(pth, dia)
    r = requests.head(arq)

    # testa se arquivo existe no opendap
    if r.status_code == requests.codes.ok:
        df = Dataset(arq, 'r')
        pcp = df.variables['pet'][:]
        lat = df.variables['latitude'][:]
        lon = df.variables['longitude'][:]
        df.close()
    else:
        print "erro"

    # Calculate Thiessen
    arq_aux = Thiessen.thiessen(pet, lat, lon, basin_fullname, -1., sep=',', usenc=True)
    print arq

    # Save Thiessen in .asc and .nc through the declared variables above
    write_thiessen(arq_aux, start_date, end_date, 'daily', 'pet', 'inmet', 'obs', 'grande',
                   lead0_date.strftime('%Y%m%d') + '00', path_out)









