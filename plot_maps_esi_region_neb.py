# -*- coding: utf-8 -*-

__author__ = "Leidinice Silva"
__copyright__ = "Copyright 2016, Funceme Hydropy Project"
__credits__ = ["Francisco Vasconcelos Junior", "Marcelo Rodrigues"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marcelo Rodrigues"
__email__ = "leidinice.silvae@funceme.br"
__date__ = 07/25/2016

import os
import re
import wget
import netCDF4
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib as mpl ; mpl.use('Agg')

from datetime import date
from gridfill import fill
from matplotlib import colors as c
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import maskoceans
from matplotlib.colors import BoundaryNorm
from PyFuncemeClimateTools import DefineGrid as Dg


def arguments():

    global args

    parser = argparse.ArgumentParser(description=__description__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--target_month', type=str, default='',
                        help='Insert target month')
    parser.add_argument('--target_year', type=int, default=date.today().year,
                        help='Insert target year')

    args = parser.parse_args()

    return args

def plotmap(data, lat, lon, **kwargs):

    deltalat = np.mean(np.diff(lat))/2.
    deltalon = np.mean(np.diff(lon))/2.

    lat = lat - deltalat
    lon = lon - deltalon

    # TODO: http://stackoverflow.com/a/7172970

    maptype       = kwargs.pop('maptype', 'shade')
    fig_name      = kwargs.pop('fig_name', 'my_picture.png')
    fig_title     = kwargs.pop('fig_title', 'My Title')
    barloc        = kwargs.pop('barloc', 'right')
    barinf        = kwargs.pop('barinf', 'neither')
    ocean_mask    = kwargs.pop('ocean_mask', 0)
    resol         = kwargs.pop('resol', 'c')
    parallels     = kwargs.pop('parallels', np.arange(-90., 91., 10.))
    meridians     = kwargs.pop('meridians', np.arange(-160., 161., 10.))
    dirshapefile  = kwargs.pop('dirshapefile', '/shp/brazil')  #TODO: Rever essa opção
    shapefile     = kwargs.pop('shapefile', 'brazil')
    barlevs       = kwargs.pop('barlevs', np.linspace(np.min(data),
                                np.max(data), 10))
    barcolor      = kwargs.pop('barcolor', ('#000092', '#0000ED',
                               '#0031FF', '#0059FF', '#0081FF',
                               '#00A5FF', '#00CDFF', '#2CFFCA',
                               '#6DFF8A', '#8AFF6D', '#CAFF2C',
                               '#EBFF0C', '#FFB900', '#FF7300',
                               '#FF2900', '#BF0000', '#920000'))
    latsouthpoint = kwargs.pop('latsouthpoint', np.nanmin(lat)+1)
    latnorthpoint = kwargs.pop('latnorthpoint', np.nanmax(lat))
    lonwestpoint  = kwargs.pop('lonwestpoint', np.min(lon)+1)
    loneastpoint  = kwargs.pop('loneastpoint', np.max(lon))
    # showplot      = kwargs.pop('showplot', 0)

    # verifica se o dado é 2d
    if data.ndim != 2:
        print "\n ... Erro ao plotar mapa! Verifique as dimensões ...\n"
        exit()

    mymap = Basemap(projection='cyl', llcrnrlat=latsouthpoint,
                    urcrnrlat=latnorthpoint, llcrnrlon=lonwestpoint,
                    urcrnrlon=loneastpoint, resolution=resol,
                    suppress_ticks=True,)

    mymap.drawmeridians(meridians, labels=[0, 0, 0, 1],
                        linewidth=0., fontsize=11)

    mymap.drawparallels(parallels, labels=[1, 0, 0, 0],
                        linewidth=0., fontsize=11)

    lons, lats = np.meshgrid(lon, lat)

    x, y = mymap(lons, lats)

    # Máscara do oceano
    if ocean_mask == 1:

        data = maskoceans(lons, lats, data, inlands=False)

    if barinf == 'both':

        cpalunder = barcolor[0]
        cpalover = barcolor[-1]
        barcolor = barcolor[:-1]
        my_cmap = c.ListedColormap(barcolor)
        my_cmap.set_under(cpalunder)
        my_cmap.set_over(cpalover)

    elif barinf == 'max':

        cpalover = barcolor[-1]
        barcolor = barcolor[:-1]
        my_cmap = c.ListedColormap(barcolor)
        my_cmap.set_over(cpalover)
        my_cmap.set_bad('w', 1.)

    elif barinf == 'min':

        cpalunder = barcolor[0]
        barcolor = barcolor[1:]
        my_cmap = c.ListedColormap(barcolor)
        my_cmap.set_under(cpalunder)
        # my_cmap.set_bad('#D3D3D3', 1.) # Se o valor for no valid data
        my_cmap.set_bad('w', 1.)

    elif barinf == 'neither':

        my_cmap = c.ListedColormap(barcolor)

    else:

        print 'Use neither, both, min, max para barinfo!'
        exit(1)

    norm = BoundaryNorm(barlevs, ncolors=my_cmap.N, clip=True)

    if maptype == 'shade':
        cs = plt.pcolormesh(x, y, data, cmap=my_cmap, norm=norm)

    elif maptype == 'fill':

        cs = plt.contourf(x, y, data, levels=barlevs, extend=barinf,
                          latlon=True, norm=norm, cmap=my_cmap)

    elif maptype == 'fillc':

        mpl.rcParams['contour.negative_linestyle'] = 'solid'

        cs = plt.contourf(x, y, data, levels=barlevs, extend=barinf,
                          latlon=True, norm=norm, cmap=my_cmap)

        cs = plt.contour(x, y, data, levels=barlevs, colors='k',
                         linewidths=0.5)

        plt.clabel(cs, inline=True, inline_spacing=2, fontsize=7.,
                   fmt='%1.1f')

    elif maptype == 'contour':

        mpl.rcParams['contour.negative_linestyle'] = 'solid'

        cs = plt.contour(x, y, data, levels=barlevs, colors='k', linewidths=0.5)

        plt.clabel(cs, inline=True, inline_spacing=2, fontsize=7., fmt='%1.1f')

    else:

        print "ERRO!: Use shade, fill or fillc para o tipo de mapa (maptype)!"
        exit(1)

    if maptype == 'shade' or maptype == 'fill':

        bar = mymap.colorbar(cs, location=barloc, spacing='uniform', ticks=barlevs, extendfrac='auto', extend=barinf, pad="5.2%")
        bar.ax.tick_params(labelsize=10)

    plt.title(fig_title, fontsize=16.9)

    local_dir = os.path.dirname(__file__)

    if shapefile != None:
        mymap.readshapefile(local_dir + dirshapefile, shapefile, drawbounds=True, linewidth=0.5, color='k')

    # plt.get_current_fig_manager().resize(6, 8)
    plt.savefig(fig_name, bbox_inches='tight')

    plt.show()

    plt.close()

def caso_shape(data, lat, lon, shp):

    shapeurl = "http://opendap2.funceme.br:8001/data/utils/shapes/regioes/{0}".format(shp)
    line = re.sub('[/]', ' ', shapeurl); line = line.split(); Ptshape = '/tmp/' + line[-1]; shpn = line[-2].title()
    if not os.path.exists(Ptshape): Ptshape = wget.download(shapeurl, '/tmp', bar=None)

    xy = np.loadtxt(Ptshape)
    xx = xy[:, 0]
    yy = xy[:, 1]

    plt.plot(xx, yy, color='k', linewidth=1.)

    shpn  = line[-2].title()

    # Quando lon de 0:360 mudando para -180:180
    if not np.any(lon<0): lon=np.where(lon>180, lon-360, lon)

    # PONTOS DO POLIGONO QUE SERA MASCARADO
    Ptsgrid, lonlatgrid, Ptmask = Dg.pointinside(lat, lon, shapefile=Ptshape)

    # APLICANDO MASCARA DO POLIGONO NO DADO DE ENTRADA
    VarMasked_data = np.ma.array(data[:, :, :], mask=np.tile(Ptmask, (data.shape[0], 1)))

    return VarMasked_data, Ptmask

def caso_shape_so_plot(shp):

    shapeurl = "http://opendap2.funceme.br:8001/data/utils/shapes/{0}".format(shp)
    line = re.sub('[/]', ' ', shapeurl); line = line.split(); Ptshape = '/tmp/' + line[-1]; shpn = line[-2].title()
    if not os.path.exists(Ptshape): Ptshape = wget.download(shapeurl, '/tmp', bar=None)

    xy = np.loadtxt(Ptshape)
    xx = xy[:, 0]
    yy = xy[:, 1]

    plt.plot(xx, yy, color='k', linewidth=1.)

def txtbox(text, xpos, ypos, fontsize, (col, lin, pos)):
    txtsp = plt.subplot(col, lin, pos)
    txt = text
    props = dict(boxstyle='round', facecolor='wheat', alpha=0)
    txtsp.text(xpos, ypos, txt, transform=txtsp.transAxes, fontsize=fontsize, verticalalignment='top',
               bbox=props)


print ""
print "Iniciando plote de ESI para o NEB!!!"
print ""

dic_mon = {'':'', 1:'Jan', 2:'Fev', 3:'Mar', 4:'Abr', 5:'Mai', 6:'Jun', 7:'Jul',  8:'Ago', 9:'Set', 10:'Out',11:'Nov', 12:'Dez'}
dic_esi = {'Jan': 3, 'Fev': 7, 'Mar': 12, 'Abr': 16, 'Mai': 20, 'Jun': 25,
           'Jul': 29, 'Ago': 34, 'Set': 38, 'Out': 42, 'Nov': 47, 'Dez': 51}
idx_mon = dict(map(reversed, dic_esi.iteritems()))

# Parametros
arguments()
month = args.target_month
newmonth = ''
if month:
    newmonth = dic_mon[int(month)]

path_in = '../data/mon/nc/'
pathout_fig = "../data/mon/fig/neb/"

mask_name = 'mask_esi_SA.nc'
maskin = netCDF4.Dataset(mask_name)
mask_data = maskin.variables['esi'][0, :, :]  # Declaring variable under study to calculate the thiessenaux_in
mask_sa = ma.getmask(mask_data)

curr_year = args.target_year

name = 'esi_4WK_{0:4d}_SA.nc'.format(curr_year)
data = netCDF4.Dataset(str(path_in) + name)
aux = data.variables['esi'][:, :, :]
lat = data.variables['latitude'][:] # Declaring latitude
lon = data.variables['longitude'][:] # Declaring longitude
min_lat, min_lon, min_lat_index, min_lon_index = Dg.gridpoint(lat, lon, -20., -50.0)
max_lat, max_lon, max_lat_index, max_lon_index = Dg.gridpoint(lat, lon, 0., -34.)
lats = lat[min_lat_index:max_lat_index]
lons = lon[min_lon_index:max_lon_index]

if not newmonth:
    new_idx = aux.shape[0] - 1
    flag = False
    while flag == False:
        try:
            newmonth = idx_mon[new_idx]
            flag = True
        except:
            new_idx = new_idx - 1

aux_in = data.variables['esi'][dic_esi[newmonth], :, :]
aux_in = ma.masked_where(aux_in <= -6., aux_in)
aux_in = ma.masked_where(aux_in >=  6., aux_in)

#Filling Gap
kw = dict(eps=1e-4, relax=0.6, itermax=1e4, initzonal=False,
          cyclic=False, verbose=True)
aux_in, converged = fill(aux_in, 1, 0, **kw)

aux_in = np.expand_dims(aux_in, axis=0)
aux_in = np.ma.array(aux_in, mask=mask_sa)
aux_in = aux_in[:, min_lat_index:max_lat_index, min_lon_index:max_lon_index]
Dmasked, mask = caso_shape(aux_in[:], lats, lons, '/nordeste_do_brasil/nordeste_do_brasil.txt')

# Ploting
caso_shape_so_plot('/estados/alagoas/alagoas.txt')
caso_shape_so_plot('/estados/bahia/bahia.txt')
caso_shape_so_plot('/estados/ceara/ceara.txt')
caso_shape_so_plot('/estados/maranhao/maranhao.txt')
caso_shape_so_plot('/estados/paraiba/paraiba.txt')
caso_shape_so_plot('/estados/pernambuco/pernambuco.txt')
caso_shape_so_plot('/estados/piaui/piaui.txt')
caso_shape_so_plot('/estados/rio_grande_do_norte/rio_grande_do_norte.txt')
caso_shape_so_plot('/estados/sergipe/sergipe.txt')

y1, y2, x1, x2 = -20, 0, -50, -34

cor1 = ('#750000', '#ff0000', '#ff8000', '#fcd17d', '#ffff00', '#ffffff')  # Paleta
lev1 = (-2., -1.6, -1.3, -0.8, -0.5, 0.5)

txtbox(u'Fonte: NOAA/USDA-ARS/NASA\nElaboração: FUNCEME', 0.46, 0.06, 8, (1, 1, 1))

figou1 = '{2}esi_neb_4WK_{0}_{1}.png'.format(newmonth.lower(), curr_year, pathout_fig)
title1 = u'Índice de Estresse Evaporativo \n{0} - {1}'.format(newmonth, curr_year)

pf = 3.

plotmap(Dmasked[0, :, :], lats, lons,
        latsouthpoint=y1, latnorthpoint=y2, lonwestpoint=x1, loneastpoint=x2, ocean_mask=1, shapefile=None,
        fig_name=figou1, fig_title=title1, barcolor=cor1, barlevs=lev1, barinf='min', barloc='right',
        parallels=np.arange(-90., 91., pf), meridians=np.arange(-160., 161., pf))

print ""
