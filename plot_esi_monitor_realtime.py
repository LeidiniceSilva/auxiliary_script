# -*- coding: utf-8 -*-

__author__ = "Leidinice Silva"
__copyright__ = "Copyright 2016, Funceme Hydropy Project"
__credits__ = ["Francisco Vasconcelos Junior", "Marcelo Rodrigues"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marcelo Rodrigues"
__email__ = "leidinice.silvae@funceme.br"
__date__ = 07/25/2016

# Monitoring Evaporative Stress Index

# Import Datas
import netCDF4
import matplotlib as mpl; mpl.use('Agg')  # Descomente para não mostrar a janela em cada plot
import re
import wget
import numpy as np
import os
import matplotlib.pyplot as plt
import numpy.ma as ma
import scipy.stats as st
from PyFuncemeClimateTools import DefineGrid as Dg
from PyFuncemeClimateTools import PlotMaps as pm


def caso_shape(data, lat, lon, shp):

    shapeurl = "http://opendap2.funceme.br:8001/data/utils/shapes/regioes/{0}".format(shp)
    line = re.sub('[/]', ' ', shapeurl); line = line.split(); Ptshape = '/tmp/' + line[-1]; shpn = line[-2].title()
    if not os.path.exists(Ptshape): Ptshape = wget.download(shapeurl, '/tmp', bar=None)

    xy = np.loadtxt(Ptshape)
    xx = xy[:, 0]
    yy = xy[:, 1]

    plt.plot(xx, yy, color='k', linewidth=0.5)
    plt.show()
    plt.savefig('Contorno_{0}.png'.format(shpn))
    # plt.close()

    shpn  = line[-2].title()

    # Quando lon de 0:360 mudando para -180:180
    if not np.any(lon<0): lon=np.where(lon>180, lon-360, lon)

    # PONTOS DO POLIGONO QUE SERA MASCARADO
    Ptsgrid, lonlatgrid, Ptmask = Dg.pointinside(lat, lon, shapefile=Ptshape)

    # APLICANDO MASCARA DO POLIGONO NO DADO DE ENTRADA
    VarMasked_data = np.ma.array(data[:, :, :], mask=np.tile(Ptmask, (data.shape[0], 1)))

    return VarMasked_data


def caso_shape_so_plot(shp):

    shapeurl = "http://opendap2.funceme.br:8001/data/utils/shapes/{0}".format(shp)
    line = re.sub('[/]', ' ', shapeurl); line = line.split(); Ptshape = '/tmp/' + line[-1]; shpn = line[-2].title()
    if not os.path.exists(Ptshape): Ptshape = wget.download(shapeurl, '/tmp', bar=None)
    xy = np.loadtxt(Ptshape)
    xx = xy[:, 0]
    yy = xy[:, 1]

    plt.plot(xx, yy, color='k', linewidth=0.5)
    plt.show()


def txtbox(text, xpos, ypos, fontsize, (col, lin, pos)):
    txtsp = plt.subplot(col, lin, pos)
    txt = text
    props = dict(boxstyle='round', facecolor='wheat', alpha=0)
    txtsp.text(xpos, ypos, txt, transform=txtsp.transAxes, fontsize=fontsize, verticalalignment='top',
               bbox=props)


# Directories of input and output data
path_in = "/home/leidinice/Documentos/esi/results/data/esi/realtime/2016/"
name = 'esi_4WK_2016274_SA.nc'
data = netCDF4.Dataset(str(path_in) + name)

lat = data.variables['latitude'][:]  # Declaring latitude
lon = data.variables['longitude'][:]  # Declaring longitude
min_lat, min_lon, min_lat_index, min_lon_index = Dg.gridpoint(lat, lon, -20., -50.0)
max_lat, max_lon, max_lat_index, max_lon_index = Dg.gridpoint(lat, lon, 0., -34.)
lats = lat[min_lat_index:max_lat_index]
lons = lon[min_lon_index:max_lon_index]

aux_in = data.variables['esi'][:, min_lat_index:max_lat_index, min_lon_index:max_lon_index]

aux_in = ma.masked_where(aux_in <= -6., aux_in)
aux_in = ma.masked_where(aux_in >= 6., aux_in)

txtbox(u'Fonte: NOAA/ESSIC and USDA-ARS\nElaboração: FUNCEME', 0.46, 0.06, 8, (1, 1, 1))

Dmasked = caso_shape(aux_in[:], lats, lons, '/nordeste_do_brasil/nordeste_do_brasil.txt')
caso_shape_so_plot('/estados/alagoas/alagoas.txt')
caso_shape_so_plot('/estados/bahia/bahia.txt')
caso_shape_so_plot('/estados/ceara/ceara.txt')
caso_shape_so_plot('/estados/maranhao/maranhao.txt')
caso_shape_so_plot('/estados/paraiba/paraiba.txt')
caso_shape_so_plot('/estados/pernambuco/pernambuco.txt')
caso_shape_so_plot('/estados/piaui/piaui.txt')
caso_shape_so_plot('/estados/rio_grande_do_norte/rio_grande_do_norte.txt')
caso_shape_so_plot('/estados/sergipe/sergipe.txt')

y1, y2, x1, x2 = -20., 0., -50., -34.
cor1 = ('#750000', '#ff0000', '#ff8000', '#fcd17d', '#ffff00', '#ffffff', '#00ffff', '#7dd1fa',
        '#0080ff', '#0000ff', '#000075')  # Paleta
lev1 = (-2., -1.6, -1.3, -0.8, -0.5, 0.5, 0.8, 1.3, 1.6, 2.)

figou1 = 'ESI_Regiao_NEB_4WK_Setembro.png'
title1 = u'Índice de Estresse Evaporativo \nSetembro'

pm.plotmap(
    Dmasked[0, :, :], lats, lons,
    latsouthpoint=y1, latnorthpoint=y2, lonwestpoint=x1, loneastpoint=x2, ocean_mask=1, shapefile=None,
    fig_name=figou1, fig_title=title1, barcolor=cor1, barlevs=lev1, barinf='both', barloc='right')

plt.close('all')
plt.cla()









