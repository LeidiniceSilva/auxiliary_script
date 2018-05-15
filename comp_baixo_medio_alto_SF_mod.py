# -*- coding: utf-8 -*-

__author__ = "Leidinice Silva"
__copyright__ = "Copyright 2016, Funceme Hydropy Project"
__credits__ = ["Francisco Vasconcelos Junior", "Marcelo Rodrigues"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marcelo Rodrigues"
__email__ = "leidinice.silvae@funceme.br"
__date__ = 07/11/2016


# Plot maps

# Import Datas
import netCDF4
import matplotlib as mpl ; mpl.use('Agg')  # Descomente para não mostrar a janela em cada plot
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

    shapeurl = "http://opendap2.funceme.br:8001/data/utils/shapes/{0}".format(shp)

    line = re.sub('[/]', ' ', shapeurl) ; line = line.split() ; Ptshape='/tmp/'+line[-1] ; shpn = line[-2].title()

    if not os.path.exists(Ptshape): Ptshape = wget.download(shapeurl,'/tmp',bar=None)

    xy = np.loadtxt(Ptshape)
    xx = xy[:,0]
    yy = xy[:,1]

    # plt.plot(xx, yy, color='k',  linewidth=1.0)
    # plt.show()

    #plt.savefig('Contorno_{0}.png'.format(shpn))
    #plt.close()

    shpn  = line[-2].title()

    # Quando lon de 0:360 mudando para -180:180
    if not np.any(lon<0): lon=np.where(lon>180, lon-360, lon)

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


# Directories of input and output data

path_in = "/home/leidinice/Documentos/esi/results/data/esi/clim/dado_12WK/" # Directory entry of  data

trim = ['NDJ', 'DJF', 'JFM', 'FMA', 'MAM', 'AMJ', 'MJJ', 'JJA', 'JAS', 'ASO', 'SON', 'OND']
dic_esi = {'NDJ': 3, 'DJF': 7, 'JFM': 12, 'FMA': 16, 'MAM': 20, 'AMJ': 25, 'MJJ': 29, 'JJA': 34, 'JAS': 38, 'ASO': 42, 'SON': 47, 'OND': 51}

for season in trim:
    sf_baixo = []
    sf_medio = []
    sf_alto = []

    for ANO in range(2003, 2015 + 1):

        print "ANO", ANO, season

        #Declaring variables

        name = 'esi_12WK_SMN_{0:04d}_SA.nc'.format(ANO)
        data = netCDF4.Dataset(str(path_in) + name)


        lat = data.variables['latitude'][:]  # Declaring latitude
        lon = data.variables['longitude'][:]  # Declaring longitude

        min_lat, min_lon, min_lat_index, min_lon_index = Dg.gridpoint(lat, lon, -22., -48.0)
        max_lat, max_lon, max_lat_index, max_lon_index = Dg.gridpoint(lat, lon, -6., -36.)

        # min_lat, min_lon, min_lat_index, min_lon_index = Dg.gridpoint(lat, lon, -20., -50.0)
        # max_lat, max_lon, max_lat_index, max_lon_index = Dg.gridpoint(lat, lon, -0., -34.)

        lats = lat[min_lat_index:max_lat_index]
        lons = lon[min_lon_index:max_lon_index]

        aux_in = data.variables['esi'][dic_esi[season], min_lat_index:max_lat_index, min_lon_index:max_lon_index] # Declaring variable under study to calculate the thiessen aux_in
        aux_in = np.expand_dims(aux_in, axis=0)

        aux_in = ma.masked_where(aux_in <= -6., aux_in)
        aux_in = ma.masked_where(aux_in >=  6., aux_in)

        # print aux_in.shape
        print "baixo_sao_francisco"
        Dmasked1 = caso_shape(aux_in[:, :, :], lats, lons, 'bacias/sao_francisco_baixo/baixo_sao_francisco.txt')
        # mapa_show(Dmasked1, lat, lon)
        fldmean1 = np.nanmean(Dmasked1)
        # print fldmean1

        print "medio_sao_francisco"
        Dmasked2 = caso_shape(aux_in[:, :, :], lats, lons, 'bacias/sao_francisco_medio/medio_sao_francisco.txt')
        # mapa_show(Dmasked2, lat, lon)
        fldmean2 = np.nanmean(Dmasked2)
        # print fldmean2

        print "alto_sao_francisco"
        Dmasked3 = caso_shape(aux_in[:, :, :], lats, lons, 'bacias/sao_francisco_alto/alto_sao_francisco.txt')
        # mapa_show(Dmasked3, lat, lon)
        fldmean3 = np.nanmean(Dmasked3)
        # print fldmean3

        sf_baixo.append(fldmean1)
        sf_medio.append(fldmean2)
        sf_alto.append(fldmean3)

    print "Fazendo Plot"
    t = np.arange(2003, 2015+1)
    graf1 = plt.plot(t, sf_baixo, 'r--', t, sf_medio, 'b--', t, sf_alto, 'k--')
    plt.title(u'a. Série Temporal de ESI \n Baixo, Médio e Alto São Francisco - {0}'.format(season))
    plt.xlabel(u'Período'.format(season))
    plt.ylabel(u'Anomalia Padronizada'.format(season))
    legenda = ('sf_baixo', 'sf_medio', 'sf_alto')
    plt.legend(graf1, legenda, frameon=False)

    ex_x = np.arange(2003, 2015 + 1)
    eixo_x = ['2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014', '2015']

    plt.xticks(ex_x, eixo_x, fontsize=10)

    # plt.xlim(2003, 2015)
    plt.ylim(0, 300)
    print "salvando"
    plt.savefig(u'Evolucao_de_ESI_no_Baixo_Medio_Alto_Sao_Francisco_{0}.png'.format(season))
    plt.close('all')
    plt.cla()






















