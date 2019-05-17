# -*- coding: utf-8 -*-

__author__ = "Leidinice Silva"
__copyright__ = "Copyright 2016, Funceme Hydropy Project"
__credits__ = ["Francisco Vasconcelos Junior", "Marcelo Rodrigues", "Enzo Pinheiro"]
__license__ = "GPL"
__version__ = "1.0.1"
__email__ = "leidinice.silvae@funceme.br"
__date__ = 07/25/2016


# Description 
# Importar arquivos .nc para .tif

# Import Datas
import os

nc_name = 'esi_12WK_{0:4d}_SA.nc'
tif_name = "esi_12WK_2016183_SA.tif"

cmd = "gdal_translate {0} {1}".format(nc_name, tif_name)
os.system(cmd)
