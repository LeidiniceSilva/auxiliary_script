# -*- coding: utf-8 -*-

__author__ = "Leidinice Silva"
__copyright__ = "Copyright 2016, Funceme Hydropy Project"
__credits__ = ["Francisco Vasconcelos Junior", "Marcelo Rodrigues"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marcelo Rodrigues"
__email__ = "leidinice.silvae@funceme.br"
__date__ = 07/25/2016

# Description
# Verificação por macrobacia do São Francisco

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


def txtbox(text, xpos, ypos, fontsize, (col, lin, pos)):
    txtsp = plt.subplot(col, lin, pos)
    txt = text
    props = dict(boxstyle='round', facecolor='wheat', alpha=0)
    txtsp.text(xpos, ypos, txt, transform=txtsp.transAxes, fontsize=fontsize, verticalalignment='top', bbox=props)

# Directories of input and output data

path_in = "/home/leidinice/Documentos/esi/dados/esi/clim/4WK/"
pathout_fig = "/home/leidinice/Documentos/esi/dados/esi/clim/4WK/"

mes = ['Jan', 'Fev', 'Mar', 'Abr']
dic_esi = {'Jan': 3, 'Fev': 7, 'Mar': 12, 'Abr': 16}

# mes = ['Jan', 'Fev', 'Mar', 'Abr', 'Mai', 'Jun', 'Jul', 'Ago', 'Set', 'Out', 'Nov', 'Dez']
# dic_esi = {'Jan': 3, 'Fev': 7, 'Mar': 12, 'Abr': 16, 'Mai': 20, 'Jun': 25, 'Jul': 29, 'Ago': 34, 'Set': 38, 'Out': 42, 'Nov': 47, 'Dez': 51}

for i, season in enumerate(mes):

    print season

    for ANO in range (2014, 2016 + 1):

        print ANO, dic_esi[season]

        name = 'esi_4WK_SMN_{0:4d}_SA.nc'.format(ANO)
        data = netCDF4.Dataset(str(path_in) + name)
        aux_in = data.variables['esi'][dic_esi[season], :, :] # Declaring variable under study to calculate the thiessenaux_in
        aux_in = np.expand_dims(aux_in, axis=0)
        lat = data.variables['latitude'][:] # Declaring latitude
        lon = data.variables['longitude'][:] # Declaring longitude
        aux_in = ma.masked_where(aux_in <= -6., aux_in)
        aux_in = ma.masked_where(aux_in >= 6., aux_in)

        y1, y2, x1, x2 = -60., 15., -90., -30.
        cor1 = ('#750000', '#ff0000', '#ff8000', '#fcd17d', '#ffff00', '#ffffff', '#00ffff', '#7dd1fa', '#0080ff', '#0000ff', '#000075')  # Paleta
        lev1 = (-2., -1.6, -1.3, -0.8, -0.5, 0.5, 0.8, 1.3, 1.6, 2.)

        fig = plt.figure()

        figou1 = 'ESI_Regiao_SA_4WK_{0}_{1}.png'.format(mes[i], ANO)
        title1 = u'Índice de Estresse Evaporativo \n{0} - {1}'.format(mes[i], ANO)

        txtbox(u'Fonte: NOAA/ESSIC and USDA-ARS\nElaboração: FUNCEME', 0.46, 0.06, 8, (1, 1, 1))

        pm.plotmap(
        aux_in[0, :, :], lat, lon,
        latsouthpoint=y1, latnorthpoint=y2, lonwestpoint=x1, loneastpoint=x2, ocean_mask=1,
        fig_name=figou1, fig_title=title1, barcolor=cor1, barlevs=lev1, barinf='both', barloc='right')

        plt.close('all')
        plt.cla()









