# -*- coding: utf-8 -*-

import os
from PIL import Image
from netCDF4 import Dataset
from numpy import *
import numpy as np
import math
import matplotlib as mpl
mpl.use('Agg')
from spi import spi
from PyFuncemeClimateTools import CreateNetCDF as cn
from dateutil.relativedelta import relativedelta
from PyFuncemeClimateTools import PlotMaps as pm
from datetime import datetime
from scipy.stats import gamma, norm
import numpy.ma as ma
from concatmes import concatmes

# Script para calcular o spi do periodo selecionado

datatual = "201608"

# Número de meses que será calculado o spi. (ex: nmes = 3, será calculado o spi
# [1, 3, 4, 6, 12, 18, 24] para o mês atual e os dois meses anteriores)
mesatual = int(datatual.split('-', 1)[1])

# Periodo climatologico
data_inic = "200301"
data_final = "201512"

ispi = [1, 3, 4, 6, 12, 18, 24]

stspi = ('SPI_1', 'SPI_3', 'SPI_4', 'SPI_6', 'SPI_12', 'SPI_18', 'SPI_24')

d1 = datetime.strptime(data_inic, "%Y%m")
d2 = datetime.strptime(data_final, "%Y%m")
d3 = datetime.strptime(datatual, "%Y%m")

# logo da funceme
dlogo = "{0}/figs/logo".format(os.getcwd())

# Paleta do spi
levs = (-2., -1.6, -1.3, -0.8, -0.5, 0.5, 0.8, 1.3, 1.6, 2.)
my_colors = ('#750000', '#ff0000', '#ff8000', '#fcd17d', '#ffff00',
             '#ffffff', '#00ffff', '#7dd1fa', '#0080ff', '#0000ff', '#000075')

ano1 = d1.year
ano2 = d2.year
nanos = d3.year - d1.year + 1
nmes = 12*(d3.year - d1.year) + (d3.month - d1.month) + 1

[pcp, lat, lon] = concatmes(data_inic, datatual)

[nt, ny, nx] = pcp.shape

mes = []
anos = []
idx = []
anoclim = []
mesclim = []
pmon = pcp

# selecionar periodo climatológico
cc = 0
for ii in range(pcp.shape[0]):
    jj = mod(ii, 12)+1
    pmon[ii, :, :] = pcp[ii, :, :]
    mes.append(jj)
    ano = fix(ano1 + (ii)/12)
    anos.append(int(ano))

for index, s in enumerate(anos):
    if s > ano1 - 1  and s < ano2 + 2:
        idx.append(index)
        anoclim.append(s)
phind = pmon[idx, :, :]

for index, s in enumerate(mes):
    if index >= idx[0] and index <= idx[-1]:
        mesclim.append(s)

parm = np.empty((12, 2, ny, nx))*np.nan

# Cálcudo dos parametros
for ig in range(len(ispi)):

    # Salvar as figuras
    directory = "{0}/figs/spi/spi{1}".format(os.getcwd(), ispi[ig])

    print "Calculando SPI", ispi[ig], "..."

    if not os.path.exists(directory):
        os.makedirs(directory)

    spi_matrix = np.empty((nanos*12, ny, nx))*np.nan

    for i in range(nx):
        for j in range(ny):

            ned = (ano2 - ano1 + 1)*12 + ispi[ig] - 1

            squeeze_ts = np.array(phind[0:ned, j, i])
            ts = squeeze_ts

            if np.isnan(ts).any() == True:
                ts = np.empty(0)
                parm[:, :, j, i] = NaN
            else:

                [var, phat] = spi(ts, ispi[ig], 12)
                parm[:, :, j, i] = phat

    parm[np.where(np.isnan(parm))] = -999

    parm = ma.masked_where(parm == -999, parm)

    # condicoes do número do spi para inicio do somatorio
    beg = ispi[ig] - 1

    # Cálculo do spi
    for it in range(nmes - 1, nmes):
        print it
        for i in range(nx):
            for j in range(ny):
                squeeze_ts = np.ma.sum(pmon[it-beg:it+1, j, i], 0)
                ts = squeeze_ts
                ss = np.mod(it, 12)
                parm1 = parm[ss, 0, j, i]
                parm2 = parm[ss, 1, j, i]
                if ts == 0.:
                    ts = ts + 0.1

                y = gamma.cdf(ts, parm1, scale=parm2)
                spi_matrix[it, j, i] = norm.ppf(y)

    spire = spi_matrix.reshape(-1, 12, ny, nx)

    spire[np.where(np.isnan(spire))] = -999

    spire = ma.masked_where(spire == -999, spire)

    sspi = stspi[ig]

    meses = ['Jan', 'Fev', 'Mar', 'Abr', 'Mai', 'Jun', 'Jul', 'Ago',
             'Set', 'Out', 'Nov', 'Dez']

    yr = (nmes-1)/12
    for year in range(yr, yr+1):

        for month in range(mesatual-1, mesatual):

            i = ispi[ig]
            ano = year + ano1

            # Gerar nomes das figuras
            dt = "{0}-{1:02d}".format(d1.year + year, month+1)
            d = datetime.strptime(dt, "%Y-%m")
            dspi = d - relativedelta(months=+i-1)

            if i != 1:
                strn = u'{0}/{1} - {2}/{3}'.format(meses[dspi.month-1], dspi.year,
                                                   meses[d.month-1], d.year)
            else:
                strn = u'{0}/{1}'.format(meses[dspi.month-1], dspi.year)

            figtitle = u'MERGE Monitoramento {0}\nÍndice de Precipitação ' \
                       u'Padronizado (SPI-{1})'.format(strn, i)

            figname = u'{0}/spi{1}-monit-merge-{2:02d}-{3}-{4}'\
                      .format(directory, i, month + 1, meses[month].lower(), ano)

            tifname = u'{0}/SPI_{1}{2:02d}_m{3}-MERGE.tif'\
                      .format(directory, ano, month + 1, i)

            print sspi
            print figtitle
            print figname
            print tifname

            pm.plotmap(spire[year, month, :, :], lat, lon, latsouthpoint=-20,
                       latnorthpoint=0, lonwestpoint=-50, loneastpoint=-34,
                       maptype='fill', resol='h', fig_name=figname, barloc='right',
                       barcolor=my_colors, barlevs=levs, fig_title=figtitle,
                       barinf='both', ocean_mask=1,
                       parallels=np.arange(-90, 90, 2),
                       meridians=np.arange(-180, 180, 2))

            # background = Image.open(figname + '.png')
            # lg = "{0}/FUNCEME_LOGO.png".format(dlogo)
            # foreground = Image.open(lg)
            # foreground = foreground.resize((90, 70), Image.ANTIALIAS)
            # background.paste(foreground, (334, 440), foreground)
            # background.save(figname + '.png', optimize=True, quality=95)

            # gerar nc
            fstr = "{0}.nc".format(figname)

            cn.create_netcdf(spire[year, month, :, :], lat, lon, fileout=fstr,
                             varname='spi', varunits='mm', ntime=1,
                             fillvalue=-999, varlongname='SPI')

            # Gerar Tif
            cmd = "gdal_translate {0}  {1}".format(fstr, tifname)
            os.system(cmd)
