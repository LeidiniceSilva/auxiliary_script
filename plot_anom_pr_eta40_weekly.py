# -*- coding: utf-8 -*-

"""
This script creates images of accumulated (7 and 10 days) 
    precipitation simulated by ETA40 model.
"""


import matplotlib as mpl
# mpl.use('Agg')

from matplotlib.colors import BoundaryNorm
from mpl_toolkits.basemap import Basemap
from datetime import datetime, timedelta
from matplotlib import colors as c
from netCDF4 import Dataset
import numpy as np
import argparse
# import calendar
import matplotlib.pyplot as plt
import os
from mpl_toolkits.basemap import interp, shiftgrid
from os.path import expanduser
from PyFuncemeClimateTools import DefineGrid as Dg
from PyFuncemeClimateTools.DefineDates import index_between_dates

home = expanduser("~")

__author__ = 'Leidinice Silva'
__email__ = 'leidinice.silva@funceme.br'
__date__ = '08/01/2018'
__description__='This script creates images of accumulated (0-7 and 7-15 days)' \
                    ' precipitation simulated by GFS025 model.'


def arguments():
    """
    pass arguments to __main__

    :return: args
    """

    parser = argparse.ArgumentParser(description=__description__,
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-d', '--date', type=str, nargs='?',
                        default=datetime.now().strftime('%Y-%m-%d'),
                        help='Initial date to run (YYYY-mm-dd)')

    argss = parser.parse_args()

    return argss


if __name__ == '__main__':

    args = arguments()


    # Define dates 
    now = datetime.strptime(args.date, '%Y-%m-%d')
    
    aux_dates = map(str, np.arange(args.date, 11, dtype='datetime64[D]'))
    dates = [datetime.strptime(d, '%Y-%m-%d') for d in aux_dates]
    
    i1, i2 = index_between_dates('19790101', '20170731', 
        '19810101', '20101231', time_frequency='months')

    # Open cmap file 
    cmap_data = Dataset('{0}/io/cmap/precip.mon.mean.nc'.format(home))
    cmap = cmap_data.variables['precip'][i1:i2+1, ::-1, :]
    clim_lat = cmap_data.variables['lat'][::-1]
    clim_lon = cmap_data.variables['lon'][:]
    cmap_data.close()

    cmap, clim_lon = shiftgrid(180., cmap, clim_lon, start=False)

    clim = []
    for i in range (0,12):
        month = cmap[i:360:12,:,:].mean(axis=0) 	
        clim.extend([month])
    clim = np.array(clim)

    # Open eta40 file 
    path_data = '{0}/io/eta40/{1}/{2}/pr/'.format(home, now.year, now.strftime('%Y%m'))
    
    pr = []
    nc = []
    
    # Pega os primeiros 7 arquivos correspondente a 10 dias de previsao
    for i in dates[1:]:

        name = 'pr_daily_eta40_fcst_{0}00_{1}.nc'.format(now.strftime('%Y%m%d'), i.strftime('%Y%m%d'))
        nc = path_data + name           
        aux = Dataset(nc)
        pr.append(aux.variables['pr'][0][:])
        aux.close()
	print nc
    exit()

    week_fcst = np.array(pr[:7])
    total_fcst = np.array(pr[:])
    del pr

    aux = Dataset(nc)
    lat = aux.variables['latitude'][:] 
    lon = aux.variables['longitude'][:]  
    aux.close()

    week_fcst  = week_fcst.sum(axis=0)
    total_fcst = total_fcst.sum(axis=0)
    
    # Week fcst
    z = [0, 0]
    aux_month = dates[0].month

    for d in dates[0:7]:
        if d.month == aux_month:
            z[0] += 1
        else:
            z[1] += 1
    
    cmap_clim = clim[dates[0].month-1,:,:] * z[0]
    
    if z[1] != 0:
        aux = [cmap_clim, clim[dates[:-3].month - 1, :, :] * z[1]]
        cmap_clim = np.array(aux).sum(axis=0)
    else:
        cmap_clim = np.array(cmap_clim)

    x, y = np.meshgrid(lon, lat)

    cmap_clim1 = np.array(interp(cmap_clim, clim_lon, clim_lat, x, y, order=1))
    del cmap_clim, aux_month, z

    # Total week 
    z = [0, 0]
    aux_month = dates[0].month
   
    for d in dates:
        if d.month == aux_month:
            z[0] += 1
        else:
            z[1] += 1

    cmap_clim = clim[dates[0].month - 1, :, :] * z[0]

    if z[1] != 0:
        aux = [cmap_clim, clim[dates[-1].month - 1, :, :] * z[1]]
        cmap_clim = np.array(aux).sum(axis=0)
    else:
        cmap_clim = np.array(cmap_clim)

    cmap_clim2 = np.array(interp(cmap_clim, clim_lon, clim_lat, x, y, order=1))
    del cmap_clim, aux_month, z

    # Calculate anomaly 
    week_anom_fcst = np.array(week_fcst - cmap_clim1)
    total_anom_fcst = np.array(total_fcst - cmap_clim2)
   
    print np.min(week_fcst), "week_fcst"
    print np.max(week_fcst), "week_fcst"
    print np.min(total_fcst), "total_fcst"
    print np.max(total_fcst), "total_fcst"
    print 
    print np.min(cmap_clim1), "cmap_clim1"
    print np.max(cmap_clim1), "cmap_clim1"
    print np.min(cmap_clim2), "cmap_clim2"
    print np.max(cmap_clim2), "cmap_clim2"
    exit()

    # Anomaly levs
    levsn  = (-170, -140, -110, -80, -50, -20, 20, 50, 80, 110, 140, 170)

    my_colorsn = ('#7F0000', '#EA2A2A', '#FF4500', '#FF8000', '#FCD17D', '#FFFF00', '#FFFFFF',
		  '#CCFFFF', '#6EFFFF', '#00FFFF', '#1199FF', '#2A2AEA', '#00007F')		  
		  
    print " ------- Calculate week pr anomaly from CE --------"

    title = u'Anomalia de Precipitação Acumulada (mm) \n{0} 00Z -- {1} 00Z'.format(now.strftime('%d/%m/%Y'),
                                    (now + timedelta(7)).strftime('%d/%m/%Y')) 

    s1 = u'Modelo: ETA40'
    s2 = u'Anomalia de precipitação ' \
         u'derivada climatologia mensal do CMAP (1981-2010)'
    s3 = u'Inicialização da previsão: {0} 00Z'.format(now.strftime('%d/%m/%Y'))

    # Plot figures 
    fig = plt.figure(figsize=(10,10))    
    plt.title(title, size=14)
    
    plt.text(-43., -9.6, s1, fontsize=9)
    plt.text(-43., -9.8, s2, fontsize=7)
    plt.text(-43., -10., s3, fontsize=7)
    
    maps = Basemap(projection='cyl', llcrnrlat=-9., urcrnrlat=1., llcrnrlon=-43., urcrnrlon=-36.)    

    maps.drawmeridians(np.arange(maps.lonmin,maps.lonmax+3,2),labels=[0,0,0,1], linewidth=0.0)  
    maps.drawparallels(np.arange(maps.latmin,maps.latmax+3,2),labels=[1,0,0,0], linewidth=0.0)
    
    x, y = maps(lon, lat)
    cpalunder = my_colorsn[0]
    cpalover = my_colorsn[-1]
    barcolor = my_colorsn[1:-1]
    my_cmap = c.ListedColormap(barcolor)
    my_cmap.set_under(cpalunder)
    my_cmap.set_over(cpalover)
    norml = BoundaryNorm(levsn, ncolors=my_cmap.N, clip=True)
    
    ff = plt.contourf(x, y, week_anom_fcst, cmap=my_cmap, norm=norml,levels=levsn, extend='both')

    shapedir = '{0}/apps/PyFuncemeClimateTools-0.1.12/' \
               'PyFuncemeClimateTools/shp/brazil'.format(home)
    
    maps.readshapefile(shapedir, 'brazil', drawbounds=True, linewidth=.5, color='k')
			    
    bar = fig.colorbar(ff, spacing='uniform', ticks=levsn, extendfrac='auto', extend='both', drawedges=True)
    bar.set_ticklabels(levsn)
    
    # Save figure 
    fig_name = 'ce_pr_anom_{0}_{0}_{1}_eta40.png'.format(now.strftime('%Y%m%d'),
							(now + timedelta(6)).strftime('%Y%m%d'))
    
    path_fig = './figs/{0}/{1}/'.format(now.year, now.strftime('%Y%m%d'))
    fig_save_name = path_fig + fig_name

    if not os.path.exists(path_fig):
        os.makedirs(path_fig)

    plt.savefig(fig_save_name, bbox_inches='tight')
    plt.close
    
    print " ------- Calculate total week pr anomaly from CE --------"
    
    title =  u'Anomalia de Precipitação Acumulada (mm) \n{0} 00Z -- {1} 00Z'.format(dates[0].strftime('%d/%m/%Y'), 
				    dates[-1].strftime('%d/%m/%Y'))    
    
    s1 = u'Modelo: ETA40'
    s2 = u'Anomalia de precipitação ' \
         u'derivada climatologia mensal do CMAP (1981-2010)'
    s3 = u'Inicialização da previsão: {0} 00Z'.format(now.strftime('%d/%m/%Y'))

    fig = plt.figure(figsize=(10,10))    
    plt.title(title, size=14)
    
    plt.text(-43., -9.6, s1, fontsize=9)
    plt.text(-43., -9.8, s2, fontsize=7)
    plt.text(-43., -10., s3, fontsize=7)        
    
    maps = Basemap(projection='cyl', llcrnrlat=-9., urcrnrlat=1., llcrnrlon=-43., urcrnrlon=-36.)    
			
    maps.drawmeridians(np.arange(maps.lonmin,maps.lonmax+3,2),labels=[0,0,0,1], linewidth=0.0)  
    maps.drawparallels(np.arange(maps.latmin,maps.latmax+3,2),labels=[1,0,0,0], linewidth=0.0)
    
    x, y = maps(lon, lat)
    cpalunder = my_colorsn[0]
    cpalover = my_colorsn[-1]
    barcolor = my_colorsn[1:-1]
    my_cmap = c.ListedColormap(barcolor)
    my_cmap.set_under(cpalunder)
    my_cmap.set_over(cpalover)
    norml = BoundaryNorm(levsn, ncolors=my_cmap.N, clip=True)
        		        
    ff = plt.contourf(x, y, total_anom_fcst, cmap=my_cmap, norm=norml,levels=levsn, extend='both',)

    shapedir = '{0}/apps/PyFuncemeClimateTools-0.1.12/' \
               'PyFuncemeClimateTools/shp/brazil'.format(home)
    
    maps.readshapefile(shapedir, 'brazil', drawbounds=True, linewidth=.5, color='k')
			    
    bar = fig.colorbar(ff, spacing='uniform', ticks=levsn, extendfrac='auto', extend='both', drawedges=True)
    bar.set_ticklabels(levsn)
    
    # Save figure 
    fig_name = 'ce_pr_anom_{0}_{1}_{2}_eta40.png'.format(now.strftime('%Y%'
			    'm%d'), now.strftime('%Y%m%d'),
			     dates[-2].strftime('%Y%m%d'))
    fig_save_name = path_fig + fig_name
    plt.savefig(fig_save_name, bbox_inches='tight')
    plt.close

    print " ------- Calculate week pr anomaly from NEB --------"
    
    title = u'Anomalia de Precipitação Acumulada (mm) \n{0} 00Z -- {1} 00Z'.format(now.strftime('%d/%m/%Y'),
                                    (now + timedelta(7)).strftime('%d/%m/%Y'))     
    
    s1 = u'Modelo: ETA40'
    s2 = u'Anomalia de precipitação ' \
         u'derivada climatologia mensal do CMAP (1981-2010)'
    s3 = u'Inicialização da previsão: {0} 00Z'.format(now.strftime('%d/%m/%Y'))

    fig = plt.figure(figsize=(10,10))    
    plt.title(title, size=14)
    
    plt.text(-52.1, -22.5, s1, fontsize=9)
    plt.text(-52.1, -23., s2, fontsize=7)
    plt.text(-52.1, -23.5, s3, fontsize=7)    
    
    maps = Basemap(projection='cyl', llcrnrlat=-21., urcrnrlat=4., llcrnrlon=-52., urcrnrlon=-34.)   

    maps.drawmeridians(np.arange(maps.lonmin,maps.lonmax+3,2),labels=[0,0,0,1], linewidth=0.0)  
    maps.drawparallels(np.arange(maps.latmin,maps.latmax+3,2),labels=[1,0,0,0], linewidth=0.0)
    
    x, y = maps(lon, lat)
    cpalunder = my_colorsn[0]
    cpalover = my_colorsn[-1]
    barcolor = my_colorsn[1:-1]
    my_cmap = c.ListedColormap(barcolor)
    my_cmap.set_under(cpalunder)
    my_cmap.set_over(cpalover)
    norml = BoundaryNorm(levsn, ncolors=my_cmap.N, clip=True)
        		        
    ff = plt.contourf(x, y, week_anom_fcst, cmap=my_cmap, norm=norml,levels=levsn, extend='both',)

    shapedir = '{0}/apps/PyFuncemeClimateTools-0.1.12/' \
               'PyFuncemeClimateTools/shp/brazil'.format(home)
    
    maps.readshapefile(shapedir, 'brazil', drawbounds=True, linewidth=.5, color='k')
			    
    bar = fig.colorbar(ff, spacing='uniform', ticks=levsn, extendfrac='auto', extend='both', drawedges=True)
    bar.set_ticklabels(levsn)
    
    # Save figure 
    fig_name = 'neb_pr_anom_{0}_{0}_{1}_eta40.png'.format(now.strftime('%Y%'
			    'm%d'), (now + timedelta(6)).strftime('%Y%m%d'))
    fig_save_name = path_fig + fig_name
    plt.savefig(fig_save_name, bbox_inches='tight')
    plt.close 
     
    print " ------- Calculate total week pr anomaly from NEB --------"
    
    title =  u'Anomalia de Precipitação Acumulada (mm) \n{0} 00Z -- {1} 00Z'.format(dates[0].strftime('%d/%m/%Y'), 
				    dates[-1].strftime('%d/%m/%Y'))  
    
    s1 = u'Modelo: ETA40'
    s2 = u'Anomalia de precipitação ' \
         u'derivada climatologia mensal do CMAP (1981-2010)'
    s3 = u'Inicialização da previsão: {0} 00Z'.format(now.strftime('%d/%m/%Y'))

    fig = plt.figure(figsize=(10,10))    
    plt.title(title, size=14)
    
    plt.text(-52.1, -22.5, s1, fontsize=9)
    plt.text(-52.1,  -23., s2, fontsize=7)
    plt.text(-52.1, -23.5, s3, fontsize=7)     
    
    maps = Basemap(projection='cyl', llcrnrlat=-21., urcrnrlat=4., llcrnrlon=-52., urcrnrlon=-34.)   
			
    maps.drawmeridians(np.arange(maps.lonmin,maps.lonmax+3,2),labels=[0,0,0,1], linewidth=0.0)  
    maps.drawparallels(np.arange(maps.latmin,maps.latmax+3,2),labels=[1,0,0,0], linewidth=0.0)
    
    x, y = maps(lon, lat)
    cpalunder = my_colorsn[0]
    cpalover = my_colorsn[-1]
    barcolor = my_colorsn[1:-1]
    my_cmap = c.ListedColormap(barcolor)
    my_cmap.set_under(cpalunder)
    my_cmap.set_over(cpalover)
    norml = BoundaryNorm(levsn, ncolors=my_cmap.N, clip=True)
        		        
    ff = plt.contourf(x, y, total_anom_fcst, cmap=my_cmap, norm=norml,levels=levsn, extend='both',)

    shapedir = '{0}/apps/PyFuncemeClimateTools-0.1.12/' \
               'PyFuncemeClimateTools/shp/brazil'.format(home)
    
    maps.readshapefile(shapedir, 'brazil', drawbounds=True, linewidth=.5, color='k')
			    
    bar = fig.colorbar(ff, spacing='uniform', ticks=levsn, extendfrac='auto', extend='both', drawedges=True)
    bar.set_ticklabels(levsn)
    
    # Save figure 
    fig_name = 'neb_pr_anom_{0}_{1}_{2}_eta40.png'.format(now.strftime('%Y%'
			    'm%d'), now.strftime('%Y%m%d'),
			     dates[-2].strftime('%Y%m%d'))
    fig_save_name = path_fig + fig_name
    plt.savefig(fig_save_name, bbox_inches='tight')
    plt.close
# 
