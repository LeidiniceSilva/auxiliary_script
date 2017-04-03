# -*- coding: utf-8 -*-

import argparse
import matplotlib as mpl ; mpl.use('Agg')
from numpy import *
import matplotlib.pyplot as plt
from datetime import datetime, date
from concatmes import concatmes
from PyFuncemeClimateTools import DefineGrid as Dg
import os
import numpy as np
import scipy.interpolate
from dateutil.relativedelta import *

__author__ = 'Diogenes Fontenele'
__email__ = 'diogenesfontenele13@gmail.com'
__date__ = '13/3/2017'
__description__='Esse script acessa o dado do MERGE no opendap, calcula a' \
                'precipitação média por regiao de interesse e gera série ' \
                'dos últimos doze meses, as quais são disponibilizads via' \
                'arquivo em formato .asc e gráfico de barras.'

def caso_shape(data, lat, lon, shp):

    Ptshape = '../shapes/'

    xy = np.loadtxt(Ptshape + shp)
    xx = xy[:, 0]
    yy = xy[:, 1]

    plt.plot(xx, yy, color='k', linewidth=0.5)
    plt.show()

    if not np.any(lon<0): lon=np.where(lon>180,lon-360,lon)

    # PONTOS DO POLIGONO QUE SERA MASCARADO
    Ptsgrid, lonlatgrid, Ptmask = Dg.pointinside(lat, lon, shapefile=(Ptshape + shp))

    # APLICANDO MASCARA DO POLIGONO NO DADO DE ENTRADA
    VarMasked_data = np.ma.array(data[:,:,:], mask=np.tile(Ptmask, (data.shape[0], 1))) # Array mascarada!!!

    return VarMasked_data, xx, yy

def compute_axis(mon):

    dict_mon = {
     1: ['Fev', 'Mar', 'Abr', 'Mai', 'Jun', 'Jul', 'Ago', 'Set', 'Out', 'Nov', 'Dez', 'Jan'],
     2: ['Mar', 'Abr', 'Mai', 'Jun', 'Jul', 'Ago', 'Set', 'Out', 'Nov', 'Dez', 'Jan', 'Fev'],
     3: ['Abr', 'Mai', 'Jun', 'Jul', 'Ago', 'Set', 'Out', 'Nov', 'Dez', 'Jan', 'Fev', 'Mar'],
     4: ['Mai', 'Jun', 'Jul', 'Ago', 'Set', 'Out', 'Nov', 'Dez', 'Jan', 'Fev', 'Mar', 'Abr'],
     5: ['Jun', 'Jul', 'Ago', 'Set', 'Out', 'Nov', 'Dez', 'Jan', 'Fev', 'Mar', 'Abr', 'Mai'],
     6: ['Jul', 'Ago', 'Set', 'Out', 'Nov', 'Dez', 'Jan', 'Fev', 'Mar', 'Abr', 'Mai', 'Jun'],
     7: ['Ago', 'Set', 'Out', 'Nov', 'Dez', 'Jan', 'Fev', 'Mar', 'Abr', 'Mai', 'Jun', 'Jul'],
     8: ['Set', 'Out', 'Nov', 'Dez', 'Jan', 'Fev', 'Mar', 'Abr', 'Mai', 'Jun', 'Jul', 'Ago'],
     9: ['Out', 'Nov', 'Dez', 'Jan', 'Fev', 'Mar', 'Abr', 'Mai', 'Jun', 'Jul', 'Ago', 'Set'],
    10: ['Nov', 'Dez', 'Jan', 'Fev', 'Mar', 'Abr', 'Mai', 'Jun', 'Jul', 'Ago', 'Set', 'Out'],
    11: ['Dez', 'Jan', 'Fev', 'Mar', 'Abr', 'Mai', 'Jun', 'Jul', 'Ago', 'Set', 'Out', 'Nov'],
    12: ['Jan', 'Fev', 'Mar', 'Abr', 'Mai', 'Jun', 'Jul', 'Ago', 'Set', 'Out', 'Nov', 'Dez']
    }

    return dict_mon[mon]

def check_dir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)

def arguments():
    global args

    parser = argparse.ArgumentParser(description=__description__,
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--target_date', type=str,
                        help='Target date.')
    args = parser.parse_args()

    return args

if __name__ == "__main__":

    arguments()

    if args.target_date:
        curr_date = datetime.strptime(args.target_date, '%Y%m%d')
    else:
        curr_date = date.today() - relativedelta(months=1)
    init_date = curr_date - relativedelta(months=11)

    # Currrent Date

    month = curr_date.month
    year = curr_date.year

    data_ini = "{0}{1:02d}".format(init_date.year, init_date.month)
    data_final = "{0}{1:02d}".format(curr_date.year, curr_date.month)

    print "Data em Processamento: {1}".format(data_final)

    # Directories

    ptshape = '../shapes/'
    pathbar = '../../data/pr_merge/bar_pr/{0}/{1:02d}/'.format(year, month)
    pathasc = '../../data/pr_merge/asc_serie_pr/{0}/{1:02d}/'.format(year, month)

    check_dir(pathbar)
    check_dir(pathasc)

    # Processing Merge Data

    [pcp, oldlat, oldlon] = concatmes(data_ini, data_final)
    lat = np.arange(np.min(oldlat), np.max(oldlat) + 0.01, 0.01, dtype=np.float32)
    lon = np.arange(np.min(oldlon), np.max(oldlon) + 0.01, 0.01, dtype=np.float32)

    # Computing Mean Precipitation

    filessh = os.listdir(ptshape)
    for b in range(len(filessh)):
        bac = filessh[b]
        reg = bac.split(".asc")[0]

        print "Bacia: {0}, Month: {1:02d}, Year: {2}".format(reg, month, year)

        var = []
        for t in range(pcp.shape[0]):

            aux_in = pcp[t, :, :]
            aux_in = np.expand_dims(aux_in, axis=0)
            newgrid = np.full((aux_in.shape[0], lat.shape[0], lon.shape[0]), np.nan, dtype=np.float32)
            X, Y = np.meshgrid(oldlon, oldlat)
            XI, YI = np.meshgrid(lon,lat)
            newgrid[0,:,:]=scipy.interpolate.griddata((X.flatten(),Y.flatten()), aux_in[0, :, :].flatten(), (XI,YI), method='linear')
            newgrid = np.ma.masked_invalid(newgrid,copy=True)

            min_lat, min_lon, min_lat_index, min_lon_index = Dg.gridpoint(lat, lon, -8., -42.0)
            max_lat, max_lon, max_lat_index, max_lon_index = Dg.gridpoint(lat, lon, 0., -36.)

            lats = lat[min_lat_index:max_lat_index]
            lons = lon[min_lon_index:max_lon_index]
            var_aux = newgrid[:,min_lat_index:max_lat_index, min_lon_index:max_lon_index]

            Dmasked, xx, yy = caso_shape(var_aux[:], lats, lons, bac)

            var.append(np.nanmean(Dmasked[0,:,:]))

        var = np.transpose(np.array(var))

        # Plot Bar Chart
        fig_name = pathbar + reg + '.png'

        # Save Mean Precipitation Serie
        np.savetxt(pathasc + bac, var, delimiter=',', fmt='%1.1f')

        x = arange(var.shape[0])
        plt.figure(figsize=(120, 56))
        bar_width = 0.4
        plt.bar(x, var)

        plt.ylabel('P(mm)', fontsize=200, fontweight='bold', rotation=90)
        plt.xlabel('Meses', fontsize=200, weight='bold')

        xaxis = compute_axis(month)
        plt.yticks(fontweight='bold', fontsize=194)
        plt.xticks(x + bar_width, xaxis,fontweight='bold', fontsize=194)

        plt.tick_params(axis='y', which='major', length=30, width=5)
        plt.savefig(fig_name, bbox_inches='tight', dpi=30)
