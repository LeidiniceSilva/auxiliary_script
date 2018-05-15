# -*- coding: utf-8 -*-

"""
This script creates images of accumulated (0-7 and 7-15 days) 
    precipitation simulated by GFS0p25 model.
"""


import matplotlib as mpl
# mpl.use('Agg')

from matplotlib.colors import BoundaryNorm
from mpl_toolkits.basemap import Basemap
from datetime import datetime, timedelta
from matplotlib import colors as c
from netCDF4 import Dataset

import matplotlib.pyplot as plt
import numpy as np
import argparse
import glob
import os

from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter, landscape
from reportlab.lib.units import inch
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
    
    aux_dates = map(str, np.arange(args.date, 7, dtype='datetime64[D]'))
    dates1 = [datetime.strptime(d, '%Y-%m-%d') for d in aux_dates]

    aux_dates = map(str, np.arange((now+timedelta(7)).strftime('%Y-%m-%d'), 7, dtype='datetime64[D]'))
    dates2 = [datetime.strptime(d, '%Y-%m-%d') for d in aux_dates]
    
    aux_dates = map(str, np.arange(args.date, 15, dtype='datetime64[D]'))
    dates3 = [datetime.strptime(d, '%Y-%m-%d') for d in aux_dates] 
    
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
    
    # Open gfs025 file 
    path_data = '{0}/io/gfs025/data/{1}/{2}00p25/'.format(home, now.year, now.strftime('%Y%m%d'))

    ii = 0
    pr = []
    nc = []
    
    for i in range(1,93):
        
    	if i <= 80:
            ii = i*3
        else:
            ii += 12
	    
        name = 'gfs.sam.0p25.f{0:03d}.nc'.format(ii)
        nc = path_data + name
        aux = Dataset(nc)	
        pr.append(aux.variables['apcpsfc'][0][:])
        aux.close()
    
    first_fcst = np.array(pr[:56])
    second_fcst = np.array(pr[56:89])
    total_fcst = np.array(pr[:91])
    del pr

    aux = Dataset(nc)
    lat = aux.variables['latitude'][:]
    lon = aux.variables['longitude'][:] - 360
    aux.close()
    
    first_fcst = first_fcst.sum(axis=0)   
    second_fcst = second_fcst.sum(axis=0)
    total_fcst = total_fcst.sum(axis=0)

    # First week 
    z = [0, 0]
    aux_month = dates1[0].month

    for d in dates1:
        if d.month == aux_month:
            z[0] += 1
        else:
            z[1] += 1

    cmap_clim = clim[dates1[0].month-1,:,:] * z[0]
    
    if z[1] != 0:
        aux = [cmap_clim, clim[dates1[-1].month - 1, :, :] * z[1]]
        cmap_clim = np.array(aux).sum(axis=0)
    else:
        cmap_clim = np.array(cmap_clim)

    x, y = np.meshgrid(lon, lat)
    
    cmap_clim1 = np.array(interp(cmap_clim, clim_lon, clim_lat, x, y, order=1))
    del cmap_clim, aux_month, z

    # Second week 
    z = [0, 0]
    aux_month = dates2[0].month

    for d in dates2:
        if d.month == aux_month:
            z[0] += 1
        else:
            z[1] += 1

    cmap_clim = clim[dates2[0].month - 1, :, :] * z[0]

    if z[1] != 0:
        aux = [cmap_clim, clim[dates2[-1].month - 1, :, :] * z[1]]
        cmap_clim = np.array(aux).sum(axis=0)
    else:
        cmap_clim = np.array(cmap_clim)

    cmap_clim2 = np.array(interp(cmap_clim, clim_lon, clim_lat, x, y, order=1))
    del cmap_clim, aux_month, z
    
    # Total week 
    z = [0, 0]
    aux_month = dates3[0].month

    for d in dates3:
        if d.month == aux_month:
            z[0] += 1
        else:
            z[1] += 1

    cmap_clim = clim[dates3[0].month - 1, :, :] * z[0]

    if z[1] != 0:
        aux = [cmap_clim, clim[dates3[-1].month - 1, :, :] * z[1]]
        cmap_clim = np.array(aux).sum(axis=0)
    else:
        cmap_clim = np.array(cmap_clim)

    cmap_clim3 = np.array(interp(cmap_clim, clim_lon, clim_lat, x, y, order=1))
    del cmap_clim, aux_month, z
    
    # Calculate anomaly 
    first_anom_fcst = np.array(first_fcst - cmap_clim1)
    second_anom_fcst = np.array(second_fcst - cmap_clim2)
    total_anom_fcst = np.array(total_fcst - cmap_clim3)
            
    
    
    # Anomaly levs
    levsn  = (-170, -140, -110, -80, -50, -20, 20, 50, 80, 110, 140, 170)

    my_colorsn = ('#7F0000', '#EA2A2A', '#FF4500', '#FF8000', '#FCD17D', '#FFFF00', '#FFFFFF',
		  '#CCFFFF', '#6EFFFF', '#00FFFF', '#1199FF', '#2A2AEA', '#00007F')		  

    # Accumuled levs
    levs  = (10, 20, 30, 60, 90, 120, 150, 180, 210)

    my_colors = ('#FFFFFF', '#4734F7', '#1E69A1', '#009B4A', '#58BD2E', '#AFDF10', '#FFF700',
		 '#FF8200', '#FF0300', '#9A0000')		  
		  
		  
		  
    print "--------------- Calculate first week pr anomaly from CE --------------"
    
    title = u'Anomalia de Precipitação Acumulada (mm) \n{0} 00Z -- {1} 00Z'.format(now.strftime('%d/%m/%Y'),
                                    (now + timedelta(7)).strftime('%d/%m/%Y'))      
    
    s1 = u'Modelo: GFS0p25'
    s2 = u'Anomalia de precipitação ' \
         u'derivada da climatologia mensal do CMAP (1981-2010)'
    s3 = u'Inicialização da previsão: {0} 00Z'.format(now.strftime('%d/%m/%Y'))

    # Plot figures 
    fig = plt.figure(figsize=(10,10))    
    plt.title(title, size=14)
    
    plt.text(-43., -9.6, s1, fontsize=9)
    plt.text(-43., -9.8, s2, fontsize=7)
    plt.text(-43., -10., s3, fontsize=7)
    
    maps = Basemap(projection='cyl', llcrnrlat=-9., urcrnrlat=-1., llcrnrlon=-43., urcrnrlon=-36.)    

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

    ff1 = plt.contourf(x, y, first_anom_fcst, cmap=my_cmap, norm=norml, levels=levsn, extend='both')

    shapedir = '{0}/apps/PyFuncemeClimateTools-0.1.12/' \
               'PyFuncemeClimateTools/shp/brazil'.format(home)
    
    maps.readshapefile(shapedir, 'brazil', drawbounds=True, linewidth=.5, color='k')
			    
    bar = fig.colorbar(ff1, spacing='uniform', ticks=levsn, extendfrac='auto', extend='both', drawedges=True)
    bar.set_ticklabels(levsn)
    
    # Save figure 
    fig_name = 'ce_pr_anom_{0}_{0}_{1}_gfs025.png'.format(now.strftime('%Y%m%d'),
							 (now + timedelta(6)).strftime('%Y%m%d'))
			    
    path_fig = './figs/{0}/{1}/'.format(now.year, now.strftime('%Y%m%d'))
    fig_save_name = path_fig + fig_name

    if not os.path.exists(path_fig):
        os.makedirs(path_fig)

    plt.savefig(fig_save_name, bbox_inches='tight')
    plt.close
    
    
    print "-------------- Calculate second week pr anomaly from CE --------------"
    
    title =  u'Anomalia de Precipitação Acumulada (mm) \n{0} 00Z -- {1} 00Z'.format((now + timedelta(7)).strftime('%d/%m/%Y'),
										   (now + timedelta(14)).strftime('%d/%m/%Y'))     
    
    s1 = u'Modelo: GFS0p25'
    s2 = u'Anomalia de precipitação ' \
         u'derivada da climatologia mensal do CMAP (1981-2010)'
    s3 = u'Inicialização da previsão: {0} 00Z'.format(now.strftime('%d/%m/%Y'))

    fig = plt.figure(figsize=(10,10))    
    plt.title(title, size=14)
    
    plt.text(-43., -9.6, s1, fontsize=9)
    plt.text(-43., -9.8, s2, fontsize=7)
    plt.text(-43., -10., s3, fontsize=7)    
    
    maps = Basemap(projection='cyl', llcrnrlat=-9., urcrnrlat=-1., llcrnrlon=-43., urcrnrlon=-36.)    

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
        		        
    ff2 = plt.contourf(x, y, second_anom_fcst, cmap=my_cmap, norm=norml, levels=levsn, extend='both',)

    shapedir = '{0}/apps/PyFuncemeClimateTools-0.1.12/' \
               'PyFuncemeClimateTools/shp/brazil'.format(home)
    
    maps.readshapefile(shapedir, 'brazil', drawbounds=True, linewidth=.5, color='k')
			    
    bar = fig.colorbar(ff2, spacing='uniform', ticks=levsn, extendfrac='auto', extend='both', drawedges=True)
    bar.set_ticklabels(levsn)
    
    # Save figure 
    fig_name = 'ce_pr_anom_{0}_{1}_{2}_gfs025.png'.format(now.strftime('%Y%'
			    'm%d'), (now + timedelta(7)).strftime('%Y%m%d'),
			    (now + timedelta(13)).strftime('%Y%m%d'))		    
    
    fig_save_name = path_fig + fig_name
    plt.savefig(fig_save_name, bbox_inches='tight')
    plt.close
    
    
    print "-------------- Calculate total week pr anomaly from CE --------------"
    
    title =  u'Anomalia de Precipitação Acumulada (mm) \n{0} 00Z -- {1} 00Z'.format(now.strftime('%d/%m/%Y'),
                                                                                   (now + timedelta(15)).strftime('%d/%m/%Y'))    
    s1 = u'Modelo: GFS0p25'
    s2 = u'Anomalia de precipitação ' \
         u'derivada da climatologia mensal do CMAP (1981-2010)'
    s3 = u'Inicialização da previsão: {0} 00Z'.format(now.strftime('%d/%m/%Y'))

    fig = plt.figure(figsize=(10,10))    
    plt.title(title, size=14)
    
    plt.text(-43., -9.6, s1, fontsize=9)
    plt.text(-43., -9.8, s2, fontsize=7)
    plt.text(-43., -10., s3, fontsize=7)    
    
    maps = Basemap(projection='cyl', llcrnrlat=-9., urcrnrlat=-1., llcrnrlon=-43., urcrnrlon=-36.)    

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
        		        
    ff3 = plt.contourf(x, y, total_anom_fcst, cmap=my_cmap, norm=norml, levels=levsn, extend='both',)

    shapedir = '{0}/apps/PyFuncemeClimateTools-0.1.12/' \
               'PyFuncemeClimateTools/shp/brazil'.format(home)
    
    maps.readshapefile(shapedir, 'brazil', drawbounds=True, linewidth=.5, color='k')
			    
    bar = fig.colorbar(ff3, spacing='uniform', ticks=levsn, extendfrac='auto', extend='both', drawedges=True)
    bar.set_ticklabels(levsn)
    
    # Save figure 
    fig_name = 'ce_pr_anom_{0}_{1}_{2}_gfs025.png'.format(now.strftime('%Y%m%d'), (now.strftime('%Y%m%d')), 
							 (now + timedelta(14)).strftime('%Y%m%d'))
    fig_save_name = path_fig + fig_name
    plt.savefig(fig_save_name, bbox_inches='tight')
    plt.close
    
    
    print "--------------- Calculate accum pr first week from CE --------------"
    
    title = u'Precipitação Acumulada (mm) \n{0} 00Z -- {1} 00Z'.format(now.strftime('%d/%m/%Y'),
                                                                      (now + timedelta(7)).strftime('%d/%m/%Y'))      
    
    s1 = u'Modelo: GFS0p25'
    s2 = u'Inicialização da previsão: {0} 00Z'.format(now.strftime('%d/%m/%Y'))

    # Plot figures 
    fig = plt.figure(figsize=(10,10))    
    plt.title(title, size=14)
    
    plt.text(-43., -9.6, s1, fontsize=9)
    plt.text(-43., -9.8, s2, fontsize=7)
    
    maps = Basemap(projection='cyl', llcrnrlat=-9., urcrnrlat=-1., llcrnrlon=-43., urcrnrlon=-36.)    

    maps.drawmeridians(np.arange(maps.lonmin,maps.lonmax+3,2),labels=[0,0,0,1], linewidth=0.0)  
    maps.drawparallels(np.arange(maps.latmin,maps.latmax+3,2),labels=[1,0,0,0], linewidth=0.0)
    
    x, y = maps(lon, lat)
    cpalunder = my_colors[0]
    cpalover = my_colors[-1]
    barcolor = my_colors[1:-1]
    my_cmap = c.ListedColormap(barcolor)
    my_cmap.set_under(cpalunder)
    my_cmap.set_over(cpalover)
    norml = BoundaryNorm(levs, ncolors=my_cmap.N, clip=True)

    ff4 = plt.contourf(x, y, first_fcst, cmap=my_cmap, norm=norml, levels=levs, extend='both')

    shapedir = '{0}/apps/PyFuncemeClimateTools-0.1.12/' \
               'PyFuncemeClimateTools/shp/brazil'.format(home)
    
    maps.readshapefile(shapedir, 'brazil', drawbounds=True, linewidth=.5, color='k')
			    
    bar = fig.colorbar(ff4, spacing='uniform', ticks=levs, extendfrac='auto', extend='both', drawedges=True)
    bar.set_ticklabels(levs)
    
    # Save figure 
    fig_name = 'ce_pr_accum_{0}_{0}_{1}_gfs025.png'.format(now.strftime('%Y%m%d'),
							  (now + timedelta(6)).strftime('%Y%m%d'))
			    
    fig_save_name = path_fig + fig_name
    plt.savefig(fig_save_name, bbox_inches='tight')
    plt.close
    
    
    print "-------------- Calculate accum pr second week from CE --------------"
    
    title =  u'Precipitação Acumulada (mm) \n{0} 00Z -- {1} 00Z'.format((now + timedelta(7)).strftime('%d/%m/%Y'),
								       (now + timedelta(14)).strftime('%d/%m/%Y'))     
    
    s1 = u'Modelo: GFS0p25'
    s2 = u'Inicialização da previsão: {0} 00Z'.format(now.strftime('%d/%m/%Y'))

    fig = plt.figure(figsize=(10,10))    
    plt.title(title, size=14)
    
    plt.text(-43., -9.6, s1, fontsize=9)
    plt.text(-43., -9.8, s2, fontsize=7)
    
    maps = Basemap(projection='cyl', llcrnrlat=-9., urcrnrlat=-1., llcrnrlon=-43., urcrnrlon=-36.)    

    maps.drawmeridians(np.arange(maps.lonmin,maps.lonmax+3,2),labels=[0,0,0,1], linewidth=0.0)  
    maps.drawparallels(np.arange(maps.latmin,maps.latmax+3,2),labels=[1,0,0,0], linewidth=0.0)
    
    x, y = maps(lon, lat)
    cpalunder = my_colors[0]
    cpalover = my_colors[-1]
    barcolor = my_colors[1:-1]
    my_cmap = c.ListedColormap(barcolor)
    my_cmap.set_under(cpalunder)
    my_cmap.set_over(cpalover)
    norml = BoundaryNorm(levs, ncolors=my_cmap.N, clip=True)
    
    ff5 = plt.contourf(x, y, second_fcst, cmap=my_cmap, norm=norml, levels=levs, extend='both',)

    shapedir = '{0}/apps/PyFuncemeClimateTools-0.1.12/' \
               'PyFuncemeClimateTools/shp/brazil'.format(home)
    
    maps.readshapefile(shapedir, 'brazil', drawbounds=True, linewidth=.5, color='k')
			    
    bar = fig.colorbar(ff5, spacing='uniform', ticks=levs, extendfrac='auto', extend='both', drawedges=True)
    bar.set_ticklabels(levs)
    
    # Save figure 
    fig_name = 'ce_pr_accum_{0}_{1}_{2}_gfs025.png'.format(now.strftime('%Y%'
			    'm%d'), (now + timedelta(7)).strftime('%Y%m%d'),
			    (now + timedelta(13)).strftime('%Y%m%d'))		    
    fig_save_name = path_fig + fig_name
    plt.savefig(fig_save_name, bbox_inches='tight')
    plt.close
    
    
    print "-------------- Calculate accum pr total week from CE --------------"
    
    title =  u'Precipitação Acumulada (mm) \n{0} 00Z -- {1} 00Z'.format(now.strftime('%d/%m/%Y'),
                                                                                    (now + timedelta(15)).strftime('%d/%m/%Y'))    
    
    s1 = u'Modelo: GFS0p25'
    s2 = u'Inicialização da previsão: {0} 00Z'.format(now.strftime('%d/%m/%Y'))

    fig = plt.figure(figsize=(10,10))    
    plt.title(title, size=14)
    
    plt.text(-43., -9.6, s1, fontsize=9)
    plt.text(-43., -9.8, s2, fontsize=7)
    
    maps = Basemap(projection='cyl', llcrnrlat=-9., urcrnrlat=-1., llcrnrlon=-43., urcrnrlon=-36.)    

    maps.drawmeridians(np.arange(maps.lonmin,maps.lonmax+3,2),labels=[0,0,0,1], linewidth=0.0)  
    maps.drawparallels(np.arange(maps.latmin,maps.latmax+3,2),labels=[1,0,0,0], linewidth=0.0)
    
    x, y = maps(lon, lat)
    cpalunder = my_colors[0]
    cpalover = my_colors[-1]
    barcolor = my_colors[1:-1]
    my_cmap = c.ListedColormap(barcolor)
    my_cmap.set_under(cpalunder)
    my_cmap.set_over(cpalover)
    norml = BoundaryNorm(levs, ncolors=my_cmap.N, clip=True)
        		        
    ff6 = plt.contourf(x, y, total_fcst, cmap=my_cmap, norm=norml,levels=levs, extend='both',)

    shapedir = '{0}/apps/PyFuncemeClimateTools-0.1.12/' \
               'PyFuncemeClimateTools/shp/brazil'.format(home)
    
    maps.readshapefile(shapedir, 'brazil', drawbounds=True, linewidth=.5, color='k')
			    
    bar = fig.colorbar(ff6, spacing='uniform', ticks=levs, extendfrac='auto', extend='both', drawedges=True)
    bar.set_ticklabels(levs)
    
    # Save figure 
    fig_name = 'ce_pr_accum_{0}_{1}_{2}_gfs025.png'.format(now.strftime('%Y%m%d'), (now.strftime('%Y%m%d')),
							  (now + timedelta(14)).strftime('%Y%m%d'))
    fig_save_name = path_fig + fig_name
    plt.savefig(fig_save_name, bbox_inches='tight')
    plt.close
    

    






