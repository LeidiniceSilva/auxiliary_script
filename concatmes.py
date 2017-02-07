# -*- coding: utf-8 -*-

__author__ = ["Paulo Jarbas Camurca"]
__license__ = "GPL"
__version__ = "1.0"
__email__ = "pjarbas312@gmail.com"

import os
import requests
from netCDF4 import Dataset
import numpy as np
from datetime import datetime, timedelta
import calendar
from dateutil.relativedelta import relativedelta


# Script para concatenar os dados mensais do MERGE

def concatmes(dataini, datafim):
    """
    Concatena os acumulados mensais.

    :param dataini: string com a data inicial. (ex: aaaa/mm)
    :param datafim: string com a data final. (ex: aaaa/mm)
    :return pcpnw: Numpy array com os meses concatenados.
    :return lon: Numpy array com as longitutes
    :return lat: Numpy array com as latitudes

    """

    d1 = datetime.strptime(dataini, "%Y%m")
    d2 = datetime.strptime(datafim, "%Y%m")

    nmes = 12*(d2.year - d1.year) + (d2.month - d1.month) + 1

    link = "http://opendap2.funceme.br:8001/data/dados-obs/merge/mensal"

    pcpnw = []

    for m in range(nmes):

        dini = d1 + relativedelta(months=+m)

        arq = "{0}/merge-pcp-monthly-{1}-{2:02d}.nc".format(link, dini.year,
                                                            dini.month)

        df = Dataset(arq, 'r')
        pcp = df.variables['prec'][:]
        lon = df.variables['longitude'][:]
        lat = df.variables['latitude'][:]
        df.close()

        # Concatenar os valores de passo de tempo
        pcpnw.append(pcp[0, :, :])

    pcpnw = np.array(pcpnw)

    return pcpnw, lat, lon
