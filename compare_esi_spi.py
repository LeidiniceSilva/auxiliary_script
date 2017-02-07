# -*- coding: utf-8 -*-

__author__ = "Leidinice Silva"
__copyright__ = "Copyright 2016, Funceme Hydropy Project"
__credits__ = ["Francisco Vasconcelos Junior", "Marcelo Rodrigues"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marcelo Rodrigues"
__email__ = "leidinice.silvae@funceme.br"
__date__ = 19/8/2016

# Plot graphs ESI x SPI

# Import Datas
import netCDF4
import numpy as np
import pandas as pd
import numpy.ma as ma
import scipy.stats as st
import matplotlib as mpl ; mpl.use('Agg')  # Descomente para não mostrar a janela em cada plot
import re
import wget
import os
import matplotlib.pyplot as plt
import datetime as dt
from datetime import date
from PyFuncemeClimateTools import DefineGrid as Dg
from PyFuncemeClimateTools import PlotMaps as pm
from pandas import Series, DataFrame


def caso_shape(data, lat, lon, shp):

    shapeurl = "http://opendap2.funceme.br:8001/data/utils/shapes/{0}".format(shp)
    line = re.sub('[/]', ' ', shapeurl); line = line.split(); Ptshape= '/tmp/'+line[-1]; shpn = line[-2].title()
    if not os.path.exists(Ptshape): Ptshape = wget.download(shapeurl, '/tmp', bar=None)

    xy = np.loadtxt(Ptshape)
    xx = xy[:, 0]
    yy = xy[:, 1]

    plt.plot(xx, yy, color='k',  linewidth=1.0)
    plt.show()

    plt.savefig('Contorno_{0}.png'.format(shpn))
    plt.close()

    shpn = line[-2].title()

    # Quando lon de 0:360 mudando para -180:180
    if not np.any(lon<0): lon = np.where(lon>180, lon-360, lon)

    # PONTOS DO POLIGONO QUE SERA MASCARADO
    Ptsgrid, lonlatgrid, Ptmask = Dg.pointinside(lat, lon, shapefile=Ptshape)

    # APLICANDO MASCARA DO POLIGONO NO DADO DE ENTRADA
    VarMasked_data = np.ma.array(data[:, :, :], mask=np.tile(Ptmask, (data.shape[0], 1))) # Array masked

    return VarMasked_data


def mapa_show(data3d, lat, lon):

    y1, y2, x1, x2 = -25, 0, -52, -34
    cor1 = ('#750000', '#ff0000', '#ff8000', '#fcd17d', '#ffff00',
            '#ffffff', '#00ffff', '#7dd1fa', '#0080ff', '#0000ff', '#000075')  # Colors
    lev1 = (-2., -1.6, -1.3, -0.8, -0.5, 0.5, 0.8, 1.3, 1.6, 2.)

    pm.plotmap(
        data3d[0, :, :], lat, lon,
        latsouthpoint=y1, latnorthpoint=y2, lonwestpoint=x1, loneastpoint=x2, ocean_mask=0,
        fig_name="figou1_teste.png.", barcolor=cor1, barlevs=lev1, barinf='both', barloc='right')


def txtbox(text, xpos, ypos, fontsize, (col, lin, pos)):
    txtsp = plt.subplot(col, lin, pos)
    txt = text
    props = dict(boxstyle='round', facecolor='wheat', alpha=0)
    txtsp.text(xpos, ypos, txt, transform=txtsp.transAxes, fontsize=fontsize, verticalalignment='top',
               bbox=props)


# Directories of input and output data
path_in1 = "/home/leidinice/Documentos/esi/results/data/esi/clim/dado_4WK/"
path_in2 = "/home/leidinice/Documentos/esi/results/data/spi/"

mes_esi = ['Jan', 'Fev', 'Mar', 'Abr', 'Mai', 'Jun', 'Jul', 'Ago', 'Set', 'Out', 'Nov', 'Dez']
dic_esi = {'Jan': 3, 'Fev': 7, 'Mar': 12, 'Abr': 16, 'Mai': 20, 'Jun': 25,
           'Jul': 29, 'Ago': 34, 'Set': 38, 'Out': 42, 'Nov': 47, 'Dez': 51}

esi_st = []
spi_st = []

for year in range(2003, 2015 + 1):

    for month in mes_esi:

        print year, month

        # Declaring variables
        name1 = 'esi_4WK_SMN_{0:04d}_SA.nc'.format(year)
        data1 = netCDF4.Dataset(str(path_in1) + name1)
        lat = data1.variables['latitude'][:]  # Declaring latitude
        lon = data1.variables['longitude'][:]  # Declaring longitude

        min_lat, min_lon, min_lat_index, min_lon_index = Dg.gridpoint(lat, lon, -22., -48.0)
        max_lat, max_lon, max_lat_index, max_lon_index = Dg.gridpoint(lat, lon, -6., -36.)
        lats = lat[min_lat_index:max_lat_index]
        lons = lon[min_lon_index:max_lon_index]

        aux_in1 = data1.variables['esi'][dic_esi[month], min_lat_index:max_lat_index, min_lon_index:max_lon_index]  # Declaring variable under study to calculate the thiessen aux_in
        time_esi = data1.variables['time']
        aux_in1 = np.expand_dims(aux_in1, axis=0)
        aux_in1 = ma.masked_where(aux_in1 <= -6., aux_in1)
        aux_in1 = ma.masked_where(aux_in1 >= 6., aux_in1)

        print "esi_sao_francisco"
        Dmasked1 = caso_shape(aux_in1[:, :, :], lats, lons, 'bacias/sao_francisco/sao_francisco.txt')
        print np.max(Dmasked1)
        fldmean1 = np.nanmean(Dmasked1)

        esi_st.append(fldmean1)

for ano in range(2003, 2015 + 1):

    for mes in range(0, 12):

        print ano, mes

        name2 = 'spi_1_{0:04d}.nc'.format(ano)
        data2 = netCDF4.Dataset(str(path_in2) + name2)

        lat = data2.variables['lat'][:]  # Declaring latitude
        lon = data2.variables['lon'][:]  # Declaring longitude

        min_lat, min_lon, min_lat_index, min_lon_index = Dg.gridpoint(lat, lon, -22., -48.0)
        max_lat, max_lon, max_lat_index, max_lon_index = Dg.gridpoint(lat, lon, -6., -36.)
        lats = lat[min_lat_index:max_lat_index]
        lons = lon[min_lon_index:max_lon_index]

        aux_in2 = data2.variables['spi'][mes, min_lat_index:max_lat_index, min_lon_index:max_lon_index]  # Declaring variable under study to calculate the thiessen aux_in
        time_spi = data2.variables['time']
        aux_in2 = np.expand_dims(aux_in2, axis=0)
        aux_in2 = ma.masked_where(aux_in2 <= -6., aux_in2)
        aux_in2 = ma.masked_where(aux_in2 >= 6., aux_in2)

        print "spi_sao_francisco"
        Dmasked2 = caso_shape(aux_in2[:, :, :], lats, lons, 'bacias/sao_francisco/sao_francisco.txt')
        print np.max(Dmasked2)
        fldmean2 = np.nanmean(Dmasked2)

        spi_st.append(fldmean2)

print "plot figure"
# Creating weekly mobile media daily data
dates = pd.date_range('2003-01', '2016-01', freq='M')
ESI = Series(esi_st[:], index=dates)
SPI = Series(spi_st[:], index=dates)
esispi = DataFrame({'ESI' : ESI, 'SPI' : SPI})

# Generating comparative graphs of pr thiessen
corr = ESI.corr(SPI)
esispi['2003':'2015'].plot()
txtbox(u'r = {0}'.format(round(corr, 3)), 0.1, 0.88, 9, (1, 1, 1))
plt.title(u'a. ESI-4WK x SPI-1 para a bacia do rio São Francisco\n Período: 2003 - 2015')
plt.xlabel(u'Anos')
plt.ylabel(u'Anomalia Padronizada')
legenda = ('ESI', 'SPI')
plt.legend(legenda, frameon=False)
plt.ylim(-2, 2)
plt.savefig('comparacao_esi_4wk_spi_1_2003_2015.png')
plt.close('all')
plt.cla()
exit()
