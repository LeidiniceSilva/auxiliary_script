# -*- coding: utf-8 -*-

__author__      = "Alexandre Xavier"
__email__       = "alexandre.candido.xavier.ufes@gmail.com"
__date__        = "03/20/2020"
__description__ = "This script plot climatology graphics and maps from OBS basedata"

import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy as cart
from cartopy.feature import NaturalEarthFeature, BORDERS, COASTLINE

# pegando variavel a ser trabalhada
ncfile = 'prec_daily_UT_Brazil_v2.2_19800101_19891231.nc'
prec = xr.open_mfdataset(ncfile, combine='by_coords')['prec']

# criando mascara para o continente e mar
mask_ocean = 2 * np.ones(prec.shape[1:]) * np.isnan(prec.isel(time=0))
mask_land = 1 * np.ones(prec.shape[1:]) * ~np.isnan(prec.isel(time=0))
mask_array = mask_ocean + mask_land

# definindo intervalo da seria historica para os calculos e
# reamostrando para a media anual e mensal (sao 30 anos)
date_start, date_end = '1980-01-01', '1989-12-31'
prec_normal_anual = prec.loc[dict(time=slice(date_start, date_end))].resample(time='Y').sum('time').mean('time')
prec_normal_mensal = prec.loc[dict(time=slice(date_start, date_end))].resample(time='M').sum('time').\
                                groupby('time.month').mean(dim='time')

# incorporando mascara
prec_normal_anual.coords['mask'] = (('latitude', 'longitude'), mask_array)

# definindo limites estaduais
states = NaturalEarthFeature(category='cultural', scale='50m',
                             facecolor='none',
                             name='admin_1_states_provinces_shp')

# plotando
fig = plt.figure(figsize=(7, 6))

proj = cart.crs.PlateCarree()

ax1 = fig.add_subplot(1, 1, 1, projection=proj, aspect='auto')

ax1.outline_patch.set_visible(False)

# mapa do Brasil normal anual
im = prec_normal_anual.where(prec_normal_anual.mask == 1).plot.imshow(x='longitude', y='latitude',
                               vmin=200, vmax=3400, add_colorbar=False,
                               ax=ax1, transform=proj, cmap=plt.cm.Spectral)

cbar_ax = fig.add_axes([0.2, .05, .6, .02])

cbar_ax.set_title('Precipitação (mm)')

fig.colorbar(im, cax=cbar_ax, orientation="horizontal", extend='both')

ax1.add_feature(states, edgecolor='gray', facecolor='none')

ax1.add_feature(BORDERS)

ax1.add_feature(COASTLINE)

ax1.outline_patch.set_visible(False)

ax1.set_extent([-75, -33, -34, 6])

ax1.set_position([.2, .2, .6, .6])

# plotando as normais mensais
# tamanho do grafico das normais mensais
length, width = .25, .15

# posicao para santos
lat, lon = -23.9, -46.8

# posicao absoluta do grafixo
pos_x, pos_y = .7, .15

ax2 = fig.add_subplot(1,1,1)

ax2.set_position([pos_x, pos_y, length, width])

prec_normal_mensal.sel(latitude=lat, longitude=lon, method='nearest').plot(ax=ax2)

ax2.tick_params(axis='both', which='minor', labelsize=6)

ax2.set_xticks(np.arange(1, 13, 1).tolist())

ax2.set_title('')

# anotacao
transform = cart.crs.PlateCarree()._as_mpl_transform(ax1)

ax1.annotate('Santos', xy=(lon, lat), xycoords=transform, ha='right', va='top')

# Tem que encontrar manualmente a posicao do grafico em lat o lon, neste caso lon=-40, lat=-28
ax1.annotate('', xy=(lon, lat), xytext=(-40, -28),
            arrowprops=dict(facecolor='gray',
                            arrowstyle="simple",
                            connectionstyle="arc3,rad=-0.2",
                            alpha=0.5),
            xycoords=transform,
            ha='right', va='top')

# Path out to save figure
name_out = 'test.png'

plt.savefig(name_out, dpi=200, bbox_inches='tight')

plt.show()

plt.close()
