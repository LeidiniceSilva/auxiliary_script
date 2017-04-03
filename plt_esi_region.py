# -*- coding: utf-8 -*-

"""
"""

__author__ = "Leidinice Silva"
__copyright__ = "Copyright 2016, Funceme Hydropy Project"
__credits__ = ["Leidinice Silva", "Francisco Vasconcelos Junior", "Marcelo Rodrigues", "Enzo Pinheiro"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Leidinice Silva"
__email__ = "leidinice.silvae@funceme.br"
__date__ = 07/11/2016


# DESCRIÇÃO
# Plot dos mapas

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
from PyFuncemeClimateTools import DefineGrid as dg
from PyFuncemeClimateTools import PlotMaps as pm


def caso_shape(data, lat, lon, shp):

    shapeurl = "http://opendap2.funceme.br:8001/data/utils/shapes/{0}".format(shp)

    line = re.sub('[/]', ' ', shapeurl) ; line = line.split() ; Ptshape='/tmp/'+line[-1] ; shpn = line[-2].title()
    if not os.path.exists(Ptshape): Ptshape = wget.download(shapeurl,'/tmp',bar=None)

    xx = [x.split(' ')[0] for x in open(Ptshape).readlines()]
    yy = [x.split(' ')[1] for x in open(Ptshape).readlines()]

    plt.plot(xx, yy, color='k',  linewidth=1.0)
    plt.show()
    plt.savefig('Contorno_{0}.png'.format(shpn))
    # plt.close()

    shpn  = line[-2].title()

    # Quando lon de 0:360 mudando para -180:180
    if not np.any(lon<0): lon=np.where(lon>180,lon-360,lon)

    # PONTOS DO POLIGONO QUE SERA MASCARADO
    Ptsgrid, lonlatgrid, Ptmask = dg.pointinside(lat, lon, shapefile=Ptshape)

    # APLICANDO MASCARA DO POLIGONO NO DADO DE ENTRADA
    VarMasked_data = np.ma.array(data[:,:,:], mask=np.tile(Ptmask, (data.shape[0], 1))) # Array mascarada!!!

    return VarMasked_data

# Directories of input and output data
y1, y2, x1, x2 = -25, 0, -52, -34
cor1 = ('#750000', '#ff0000', '#ff8000', '#fcd17d', '#ffff00', '#ffffff', '#00ffff', '#7dd1fa', '#0080ff', '#0000ff', '#000075') #Paleta
lev1 = (-2., -1.6, -1.3, -0.8, -0.5, 0.5, 0.8, 1.3, 1.6, 2.)

path_in = "/home/leidinice/Documentos/dados/esi/clim/12WK/" # Directory entry of  data
path_out = "/home/leidinice/Documentos/dados/esi/clim/12WK/fig/" #  Directory output of calculation pet

trim = ['OND', 'NDJ', 'DJF', 'JFM', 'FMA', 'MAM']
dic_esi = {'OND': 0, 'NDJ': 4, 'DJF': 8, 'JFM': 12, 'FMA': 17, 'MAM': 21}

for season in (trim):

    print season

    k = 0
    for ANO in range (2003, 2015 + 1):
        print k, ANO, dic_esi[season]

        #Declaring variables

        name = 'esi_12WK_{0:4d}_SA.nc'.format(ANO)
        data = netCDF4.Dataset(str(path_in) + name)
        aux_in = data.variables['esi'][dic_esi[season], :, :] # Declaring variable under study to calculate the thiessenaux_in
        aux_in = np.expand_dims(aux_in, axis=0)

        lat = data.variables['latitude'][:] # Declaring latitude
        lon = data.variables['longitude'][:] # Declaring longitude

        aux_in = ma.masked_where(aux_in <= -6., aux_in)
        aux_in = ma.masked_where(aux_in >= 6., aux_in)

        Dmasked = caso_shape(aux_in[:], lat, lon, 'bacias/sao_francisco/sao_francisco.txt')

        figou1 = 'plot_ESI_{0}_{1}.png'.format(season, ANO)
        figou2 = 'plot_ESI_NEB_{0}_{1}.png'.format(season, ANO)

        title1 = u'ESI - Bacia do Rio São Francisco - {0} - {1}'.format(season, ANO)
        title2 = u'ESI - NEB - {0} - {1}'.format(season, ANO)


        print ' ... Ploting {0} - {1} ...'.format(season, ANO)

        pm.plotmap(
        Dmasked[0, :, :], lat, lon,
        latsouthpoint=y1, latnorthpoint=y2, lonwestpoint=x1, loneastpoint=x2, ocean_mask=0,
        fig_name=figou1, fig_title=title1, barcolor=cor1, barlevs=lev1, barinf='both', barloc='right')

        pm.plotmap(aux_in[0, :, :], lat, lon, latsouthpoint=y1, latnorthpoint=y2, lonwestpoint=x1,
                   loneastpoint=x2, ocean_mask=1, fig_name=figou2, fig_title=title2, barcolor=cor1, barlevs=lev1, barinf='both', barloc='right')

        k += 1
        del aux_in





















