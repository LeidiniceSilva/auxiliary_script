# -*- coding: utf-8 -*-

""" Plot of climatology of etp. """

__author__ = "Leidinice Silva"
__email__ = "leidinice.silvae@funceme.br"
__date__ = "17/09/2016"
__description__ = "Plot of climatology of etp"


# Import Datas
import netCDF4
import matplotlib as mpl ; mpl.use('Agg')  # Descomente para não mostrar a janela em cada plot
from PyFuncemeClimateTools import PlotMaps as pm

y1, y2, x1, x2 = -60, 15, -90, -33

mes = ['Jan', 'Fev', 'Mar', 'Abr', 'Mai', 'Jun', 'Jul', 'Ago', 'Set', 'Out', 'Nov', 'Dez']
dic = {'Jan': 31., 'Fev': 28.5, 'Mar': 31., 'Abr': 30., 'Mai': 31., 'Jun': 30., 'Jul': 31., 'Ago': 31., 'Set': 30.,
       'Out': 31., 'Nov': 30., 'Dez': 31.}

cor1 = ('#750000', '#ff0000', '#ff8000', '#fcd17d', '#ffff00', '#ffffff', '#00ffff', '#7dd1fa', '#0080ff', '#0000ff',
        '#000075') #Paleta
lev1 = (10, 20, 30, 40, 50, 60., 70., 80., 90., 100.)

path_in = '/Users/diogenesfontenele/Documents/job/funceme/musf/eixo_tempo/dados/'
name = 'pet_CRU_SA_1981_2010_clim.nc'
data = netCDF4.Dataset(path_in + name)
var = data.variables['eto'][:, :, :] # Declaring variable under study to calculate the
lats = data.variables['lat'][:] # Declaring latitude
lons = data.variables['lon'][:] # Declaring longitude

for m in range(0, var.shape[0]):
    
    print m, mes[m], dic[mes[m]]

    title1 = u'Evapotranspiração Potencial\n Climatologia - {0}'.format(mes[m])
    figou1 = 'plot_ETO_CRU_{0}.png'.format(mes[m])
    pm.plotmap(var[m, :, :] * dic[mes[m]], lats, lons, latsouthpoint=y1, latnorthpoint=y2, lonwestpoint=x1,
               loneastpoint=x2, ocean_mask=1, fig_name=figou1, fig_title=title1, barcolor=cor1, barlevs=lev1,
               barinf='max', barloc='right')