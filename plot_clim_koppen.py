# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "01/15/2020"
__description__ = "This script plot annual climatology with koppen classification"

import comp_koppen
import xarray as xr
import numpy as np
import scipy.interpolate
import cartopy.crs as ccrs
import matplotlib.pylab as plt
import matplotlib.patches as mpatches

from cartopy.feature import NaturalEarthFeature

# Testando para uma localidade, neste exemplo, cidade de Verdelandia na Bahia.

# precipitacao media mensal
prec = [140.1, 87.3, 115.9, 35.7, 5.6,  1.7, 0.6, 1.8, 10.8, 53.7, 142.9, 189.7]

# temperatura media mensal
avgtemp = [25.5, 25.9, 25.76, 25.1, 23.7, 22.2, 21.9, 23.0, 25.0, 26.3, 25.6, 25.3]

# latitude. So´ serve para a definicao dos hemisferios
lat = -15.5

# Definindo o clima
clima_localidade = comp_koppen.koppen_classification(prec, avgtemp, lat)
print(clima_localidade)  # Cfa

# plotando os dados de precipitacao e temperatura
fig, ax1 = plt.subplots()
ax1.plot(prec)
ax1.set_ylabel('Precipitação mensal [mm]', color='b')
ax1.tick_params('y', colors='b')
ax2 = ax1.twinx()
ax2.plot(avgtemp, 'r')
ax2.set_ylabel('Temperatura média', color='r')
ax2.tick_params('y', colors='r')
fig.suptitle(clima_localidade)

plt.show()
exit()
