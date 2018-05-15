# coding: utf-8

from mpl_toolkits.basemap import Basemap
from matplotlib import colors as c
from matplotlib.colors import BoundaryNorm
import numpy as np
import matplotlib.pyplot as plt

__author__ = 'Marcelo Rodrigues'
__email__ = 'marcelorodriguesss@gmail.com'
__date__ = '07 Nov 2017'


def plotmap(data, lat, lon, **kwargs):

    """
    => Parâmetros obrigatórios:
    :param data:          = matriz a ser plotada (np.array 2d)
    :param lat:           = pontos da grade na latitude (np.array) (-90:90)
    :param lon:           = pontos da grade na longitude (np.array) (-180:180)

    => Parâmetros opcionais:
    resol: (string)
        c = resolução muito baixa
        l = resolução baixa
        i = resolução intermediária (padrão)
        h = resolução alta
        f = resolução muito alta
    maptype       = Tipo do mapa: fillc, fill e shaded (padrão) (string)
    fig_name      = Nome da figura: 'my_picture.png' (padrão) (string)
    fig_title     = Título da figura: 'My Title' (padrão) (string)
    barloc        = Posição da paleta: 'bellow' ou 'right' (padrão) (string)
    barinf        = Se a paleta tende a infinito: 'both', 'max', 'min' ou 'neither' (padrão) (string)
    ocean_mask    = Mascara sobre o oceano: 1 ou 0 (padrão) (inteiro)
    parallels     = Range e intervalo das paralelas do globo:
                    Ex.: np.arange(-90., 91., 10.) (padrão) (np.array)
    meridians     = Range e intervalo dos meridianos do globo:
                    Ex.: np.arange(-160., 161., 10.) (padrão) (np.array)
    dirshapefile  = Diretório onde encontra-se o shape a ser plotado. (string)
    shapefile     = Arquivo do tipo shape a ser plotado no mapa. (string)
                    O padrão é o shape 'world' disponível nesse pacote.
    latsouthpoint = ponto mais ao sul   (float)
    latnorthpoint = ponto mais ao norte (float)
    lonwestpoint  = ponto mais ao oeste (float)
    loneastpoint  = ponto mais ao leste (float)
    barlevs       = Níveis da paleta de cores (list)
    barcolor      = Cores em hexadecimal da paleta de cores (list)
    """

    # colormap.set_bad('w', 1.0)

    # TODO: http://stackoverflow.com/a/7172970

    maptype = kwargs.get('maptype', 'shade')
    fig_name = kwargs.get('fig_name', 'my_picture.png')
    fig_title = kwargs.get('fig_title', 'My Title')
    barloc = kwargs.get('barloc', 'bellow')
    barinf = kwargs.get('barinf', 'neither')
    resol = kwargs.get('resol', 'c')
    parallels = kwargs.get('parallels', np.arange(-90., 91., 15.))
    meridians = kwargs.get('meridians', np.arange(0., 360., 15.))
    barlevs = kwargs.get('barlevs', np.linspace(np.min(data),
                                                np.max(data), 10))
    barcolor = kwargs.get('barcolor', ('#000092', '#0000ED',
                               '#0031FF', '#0059FF', '#0081FF',
                               '#00A5FF', '#00CDFF', '#2CFFCA',
                               '#6DFF8A', '#8AFF6D', '#CAFF2C',
                               '#EBFF0C', '#FFB900', '#FF7300',
                               '#FF2900', '#BF0000', '#920000'))
    latsouthpoint = kwargs.get('latsouthpoint', np.nanmin(lat)+1)
    latnorthpoint = kwargs.get('latnorthpoint', np.nanmax(lat))
    lonwestpoint = kwargs.get('lonwestpoint', np.min(lon)+1)
    loneastpoint = kwargs.get('loneastpoint', np.max(lon))
    dpi = kwargs.get('dpi', 200)

    if maptype == 'shade':
        # ajuste nos pontos do plot
        deltalat = np.mean(np.diff(lat))/2.
        deltalon = np.mean(np.diff(lon))/2.
        lat = lat - deltalat
        lon = lon - deltalon

    # verifica se o dado é 2d
    if data.ndim != 2:
        print('\n => Erro ao plotar mapa! Verifique dims do dado!')
        print('\n => O dado deve ser um Numpy Array 2d!')
        print('\n => Dims atual: {dims}\n'.format(dims=data.ndim))
        raise SystemExit  # changed here

    fig = plt.figure(figsize=(12, 6))

    mymap = Basemap(projection='cyl', llcrnrlat=latsouthpoint,
                    urcrnrlat=latnorthpoint, llcrnrlon=lonwestpoint,
                    urcrnrlon=loneastpoint, resolution=resol,
                    suppress_ticks=True)

    mymap.drawmeridians(meridians, labels=[0, 0, 0, 1],
                        linewidth=0.001, fontsize=6)

    mymap.drawparallels(parallels, labels=[1, 0, 0, 0],
                        linewidth=0.001, fontsize=6)

    lons, lats = np.meshgrid(lon, lat)

    x, y = mymap(lons, lats)

    if barinf == 'both':

        cpalunder = barcolor[0]
        cpalover = barcolor[-1]
        barcolor = barcolor[1:-1]
        my_cmap = c.ListedColormap(barcolor)
        my_cmap.set_under(cpalunder)
        my_cmap.set_over(cpalover)

    elif barinf == 'max':

        cpalover = barcolor[-1]
        barcolor = barcolor[:-1]
        my_cmap = c.ListedColormap(barcolor)
        my_cmap.set_over(cpalover)

    elif barinf == 'min':

        cpalunder = barcolor[0]
        barcolor = barcolor[1:]
        my_cmap = c.ListedColormap(barcolor)
        my_cmap.set_under(cpalunder)

    elif barinf == 'neither':

        my_cmap = c.ListedColormap(barcolor)

    else:

        print('Use neither, both, min, max para barinfo!')
        raise SystemExit

    norm = BoundaryNorm(barlevs, ncolors=my_cmap.N, clip=False)

    if maptype == 'shade':

        cs = plt.pcolormesh(x, y, data, cmap=my_cmap, norm=norm)

    elif maptype == 'fill':

        cs = plt.contourf(x, y, data, levels=barlevs, extend=barinf,
                          latlon=True, norm=norm, cmap=my_cmap)

    elif maptype == 'fillc':

        plt.rcParams['contour.negative_linestyle'] = 'solid'

        cs = plt.contourf(x, y, data, levels=barlevs, extend=barinf,
                          latlon=True, norm=norm, cmap=my_cmap)

        cs = plt.contour(x, y, data, levels=barlevs, colors='k',
                         linewidths=0.5)

        plt.clabel(cs, inline=True, inline_spacing=2, fontsize=7.,
                   fmt='%1.1f')

    elif maptype == 'contour':

        plt.rcParams['contour.negative_linestyle'] = 'solid'

        cs = plt.contour(x, y, data, levels=barlevs, colors='k', linewidths=0.5)

        plt.clabel(cs, inline=True, inline_spacing=2, fontsize=7., fmt='%1.1f')

    else:

        print('ERRO!: Use shade, fill or fillc para o tipo de mapa (maptype)!')
        raise SystemExit

    if maptype == 'shade' or maptype == 'fill':

        bar = mymap.colorbar(cs, location=barloc, spacing='uniform',
                             ticks=barlevs, extendfrac='auto',
                             extend=barinf, pad='8%')

        bar.ax.tick_params(labelsize=6)

    # mymap.drawcoastlines()

    mymap.fillcontinents(color='gray', lake_color='gray')

    # mymap.drawrivers(linewidth=0.2, color='gray')

    plt.title(fig_title, fontsize=12)

    plt.savefig(fig_name, bbox_inches='tight', dpi=dpi)
   
    plt.close()

