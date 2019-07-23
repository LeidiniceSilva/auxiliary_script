# -*- coding: utf-8 -*-

from netCDF4 import Dataset, num2date
from matplotlib import pyplot as plt
import matplotlib.colors as clr
import numpy as np
from PyFuncemeClimateTools import DefineGrid as Dg
import argparse
from mpl_toolkits.basemap import shiftgrid


def arguments():
    global args

    parser = argparse.ArgumentParser(prog='Merge daily files into weekly file')

    parser.add_argument('glorys_file', help='Caminho do arquivo GLORYS')
    parser.add_argument('glorys_var', help='Variável do GLORYS')
    parser.add_argument('comparison_file', help='Caminho do arquivo para '
                                                'comparação')
    parser.add_argument('comparison_var', help='Variável do arquivo para '
                                               'comparação')
    parser.add_argument('plot_title', help='Titulo da plotagem')
    parser.add_argument('comparison_dbase_name', help='Nome da base de dados '
                                                      'para comparação')
    parser.add_argument('paleta', help='Paleta deve ser obs ou anom')
    parser.add_argument('-levs', '--levs', help='Níveis para a paleta de cores',
                        nargs='+', type=float)
    args = parser.parse_args()
    return args


def hovmoller(array, lat_array, lon_array, latmin, latmax, lonmin,
              lonmax, mean):

    """Calcula a media meridional ou zonal
    :param array:
    :param lat_array:
    :param lon_array:
    :param latmin:
    :param latmax:
    :param lonmin:
    :param lonmax:
    :param mean:
    """

    latmin, lonmin, ilatmin, ilonmin = Dg.gridpoint(lat_array, lon_array,
                                                    latmin, lonmin)
    latmax, lonmax, ilatmax, ilonmax = Dg.gridpoint(lat_array, lon_array,
                                                    latmax, lonmax)

    lat = lat_array[ilatmin:ilatmax + 1]
    lon = lon_array[ilonmin:ilonmax + 1]

    if mean == 'lat':
        hov = np.nanmean(array[:, ilatmin:ilatmax + 1, ilonmin:ilonmax + 1],
                         axis=1)
    elif mean == 'lon':
        hov = np.nanmean(array[:, ilatmin:ilatmax + 1, ilonmin:ilonmax + 1],
                         axis=2)

    return hov, lat, lon


def plot_hovmoller(x, y, hov, levels='default', cls='default'):

    if levels == 'default':
        levels = np.linspace(np.min(hov), np.max(hov), 13)
        # levels = [-2, -1.6, -1.2, -0.8, -0.4, -0.2, 0, 0.2, 0.4, 0.8, 1.2,
        #               1.6, 2]
    else:
        levels = levels

    if cls == 'default':
        colors = ('#0033FF', '#007FFF', '#0099FF', '#00B2FF', '#00CCFF',
                  '#C8FFFF', '#FFFFAA', '#FFCC00', '#FF9900', '#FF7F00',
                  '#FF3300', '#A50000',)
        cmap = clr.ListedColormap(colors)
        cmap.set_over('#B48C82')
        cmap.set_under('#000044')
    else:
        colors = cls[1:-1]
        cmap = clr.ListedColormap(colors)
        cmap.set_over(cls[-1])
        cmap.set_under(cls[0])

    h = plt.contourf(x, y, hov, levels=levels, cmap=cmap, extend='both')
    return h


def int2str_lats(lat):
    str_lats = []
    for atllat in lat[:]:
        if atllat < 0:
            str_lats.append(str(int(-(atllat - 1))) + 'S')
        elif atllat > 1:
            str_lats.append(str(int(atllat)) + 'N')
        elif atllat > 0:
            if atllat < 1:
                str_lats.append(str(int(atllat)))
    return str_lats


arguments()
glorys_dbase = args.glorys_file
glorys_var_name = args.glorys_var
glorys_nc = Dataset(glorys_dbase)
gvar = glorys_nc.variables[glorys_var_name]
units = gvar.units
if units == 'degree_Celsius':
    units = '°C'
glat = glorys_nc.variables['lat']
glon = glorys_nc.variables['lon']
gtime = glorys_nc.variables['time']
ntime = gtime.shape[0]
if len(gvar.shape) == 4:
    gvar = gvar[:, 9, :, :]

comp_dbase = args.comparison_file
comp_var_name = args.comparison_var
comp_nc = Dataset(comp_dbase)
cvar = comp_nc.variables[comp_var_name]
clat = comp_nc.variables['lat']
clon = comp_nc.variables['lon']

if len(cvar.shape) == 4:
    cvar = cvar[:, 0, :, :]
if clat[0] > clat[-1]:
    clat = clat[::-1]
    cvar = cvar[:, ::-1, :]
if np.max(clon) >= 180.:
    cvar, clon = shiftgrid(180., cvar[:], clon[:], False)

title = args.plot_title
comp_dbase_name = args.comparison_dbase_name
paleta = args.paleta
lvls = args.levs

time_axis = np.arange(0, ntime)
if ntime == 12:
    taxis = np.arange(0, 12)
    xaxis = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep',
             'Oct', 'Nov', 'Dec']
    var_type = 'ltm'
else:
    taxis = time_axis[::12]
    xaxis = np.arange(1993, 2015)
    if paleta == 'anom':
        var_type = 'anom'
    elif paleta == 'obs':
        var_type = 'obs'

if paleta == 'anom':
    my_colors = ('#000044', '#0033FF', '#007FFF', '#0099FF', '#00B2FF',
                 '#00CCFF', '#C8FFFF', '#FFFFAA', '#FFCC00', '#FF9900',
                 '#FF7F00', '#FF3300', '#A50000', '#B48C82')
elif paleta == 'obs':
    my_colors = ('#4B004B', '#960096', '#B400B4', '#C800C8', '#6400FF',
                 '#0032E1', '#0096E1', '#00FFFF', '#FFFFFF', '#32E100',
                 '#96E100', '#E1E100', '#FFFF00', '#FFAA00', '#FF8700',
                 '#FF6400', '#C80000')


# 2.5S - 30N, 38W
# Glorys database
fig = plt.figure(figsize=(16, 12))
fig.add_subplot(1, 2, 1)
ghov, ghlat, ghlon = hovmoller(gvar[:], glat[:], glon[:], -3., 30., -38., -38.,
                               'lon')

ghov = np.transpose(ghov)
plot_hovmoller(time_axis, ghlat, ghov, lvls, my_colors)
plt.xticks(taxis, xaxis)
plt.title('Glorys')
plt.xlabel('Months')
plt.ylabel('Latitudes')

# Comparison database
fig.add_subplot(1, 2, 2)
chov, chlat, chlon = hovmoller(cvar[:], clat[:], clon[:], -3., 30., -38., -38.,
                               'lon')
chov = np.transpose(chov)
h = plot_hovmoller(time_axis, chlat, chov, lvls, my_colors)
plt.xticks(taxis, xaxis)
plt.title(comp_dbase_name)
plt.xlabel('Months')

cbar_ax = fig.add_axes([0.925, 0.1, 0.025, 0.8])
cbar = fig.colorbar(h, cax=cbar_ax, ticks=lvls)
cbar.set_label(units, fontsize=12, rotation=270)
fig.suptitle(title + u'\n Meridional section at 38°W', fontsize=16)

plt.savefig('../imagens/hov_{0}38w_mo_ltm.png'.format(glorys_var_name),
            bbox_inches='tight')

# 30S - 30N, 23W
# Glorys database
fig = plt.figure(figsize=(16, 12))
fig.add_subplot(1, 2, 1)
ghov, ghlat, ghlon = hovmoller(gvar[:], glat[:], glon[:], -29.9, 30., -23.,
                               -23., 'lon')

ghov = np.transpose(ghov)
plot_hovmoller(time_axis, ghlat, ghov, lvls, my_colors)
plt.xticks(taxis, xaxis)
plt.title('Glorys')
plt.xlabel('Months')
plt.ylabel('Latitudes')

# Comparison database
fig.add_subplot(1, 2, 2)
chov, chlat, chlon = hovmoller(cvar[:], clat[:], clon[:], -29.9, 30., -23.,
                               -23., 'lon')
chov = np.transpose(chov)
h = plot_hovmoller(time_axis, chlat, chov, lvls, my_colors)
plt.xticks(taxis, xaxis)
plt.title(comp_dbase_name)
plt.xlabel('Months')

cbar_ax = fig.add_axes([0.925, 0.1, 0.025, 0.8])
cbar = fig.colorbar(h, cax=cbar_ax, ticks=lvls)
cbar.set_label(units, fontsize=12, rotation=270)
fig.suptitle(title + u'\n Meridional section at 23°W', fontsize=16)

plt.savefig('../imagens/hov_{0}23w_mo_ltm.png'.format(glorys_var_name),
            bbox_inches='tight')

# 30S - 5N, 10W
# Glorys database
fig = plt.figure(figsize=(16, 12))
fig.add_subplot(1, 2, 1)
ghov, ghlat, ghlon = hovmoller(gvar[:], glat[:], glon[:], -29.9, 5., -10., -10.,
                               'lon')
ghov = np.transpose(ghov)
plot_hovmoller(time_axis, ghlat, ghov, lvls, my_colors)
plt.xticks(taxis, xaxis)
plt.title('Glorys')
plt.xlabel('Months')
plt.ylabel('Latitudes')

# Comparison database
fig.add_subplot(1, 2, 2)
chov, chlat, chlon = hovmoller(cvar[:], clat[:], clon[:], -29.9, 5., -10., -10.,
                               'lon')
chov = np.transpose(chov)
h = plot_hovmoller(time_axis, chlat, chov, lvls, my_colors)
plt.xticks(taxis, xaxis)
plt.title(comp_dbase_name)
plt.xlabel('Months')

cbar_ax = fig.add_axes([0.925, 0.1, 0.025, 0.8])
cbar = fig.colorbar(h, cax=cbar_ax, ticks=lvls)
cbar.set_label(units, fontsize=12, rotation=270)
fig.suptitle(title + u'\n Meridional section at 10°W', fontsize=16)

plt.savefig('../imagens/hov_{0}10w_mo_ltm.png'.format(glorys_var_name),
            bbox_inches='tight')
