# -*- coding: utf-8 -*-

__author__ = "Leidinice Silva"
__copyright__ = "Copyright 2016, Funceme Hydropy Project"
__credits__ = ["Francisco Vasconcelos Junior", "Diogenes Fontenele"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marcelo Rodrigues"
__email__ = "leidinice.silva@funceme.br"
__date__ = 8/17/2016

# Import data
import requests
import calendar
import datetime
import numpy as np
from netCDF4 import Dataset
from PyFuncemeClimateTools import Thiessen
from hidropy.utils.write_thiessen import write_thiessen

# Directories of input and output data
pth = "http://opendap2.funceme.br:8001/data/dados-obs/merge/diario/"
pathsh = "/home/leidinice/Documentos/hidropy/hidropy/shapes/basins/amazonas/amazonas.asc" # Directory of the study area
path_out = "/home/leidinice/Documentos/musf/results/" #  Directory output of calculation pet
pcp = np.full((1, 245, 313), np.nan) # Declaring variable  (time step, lat and lon)

ano = 1998
mes = 2

if not pth:
    pth = "http://opendap2.funceme.br:8001/data/dados-obs/merge/diario/"

pth = "{0}{1}/prec_{1}{2:02d}".format(pth, ano, mes)

print "Acumulado de", mes, " ano:", ano
dias_mes = calendar.monthrange(ano, mes)[1]
pcpnw = []

for dia in xrange(1, dias_mes + 1):
    lead0_date = datetime.date(ano, mes, dia)
    start_date = lead0_date.strftime('%Y%m%d')
    end_date = lead0_date.strftime('%Y%m%d')
    arq = "{0}{1:02d}.ctl".format(pth, dia)
    r = requests.head(arq)

    # testa se arquivo existe no opendap
    if r.status_code == requests.codes.ok:
        df = Dataset(arq, 'r')
        pcp = df.variables['prec'][:]
        lat = df.variables['latitude'][:]
        lon = df.variables['longitude'][:]
        df.close()
    else:
        print "erro"

    # Calculate Thiessen
    arq_aux = Thiessen.thiessen(pcp, lat, lon, pathsh, -1., sep=',', usenc=True)
    print arq

    # Save Thiessen in .asc and .nc through the declared variables above
    write_thiessen(arq_aux, start_date, end_date, 'daily', 'prec', 'merge', 'fcst', 'amazonas',
                   lead0_date.strftime('%Y%m%d') + '00', path_out)









