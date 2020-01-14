# -*- coding: utf-8 -*-
"""
by: Emily Sturdivant, emilysturdivant@gmail.com
requires: python 3.6+

OVERVIEW: Try spatial operations to extract backscatter values at plot polygons
"""
# Import packages
import pandas as pd
pd.__version__
import numpy as np
import matplotlib.pyplot as plt
import os
import unicodedata
import json
from platform import python_version
import geopandas as gpd
import rasterio
from rasterio.mask import mask
from shapely.geometry import mapping

print(python_version())

#%%

#%% Initialize
# Set working directory
home = r'/Users/emilysturdivant/GitHub/biomass-espanola'

# Filenames
cwd_fname = r'/Users/emilysturdivant/Documents/CIGA/data/N19W073_18_sl_HV_F02DAR.tif'
shp_fname = os.path.join(home, 'data', 'Augier.shp')

#%%
shapefile = gpd.read_file(shp_fname, crs ='EPSG:32618')
# import fiona; help(fiona.open)

# extract the geometry in GeoJSON format
geoms = shapefile.geometry.values # list of shapely geometries
geometry = geoms[0] # shapely geometry
# transform to GeJSON format
geoms = [mapping(geoms[0])]
geoms

# extract the raster values values within the polygon
# help(rasterio.open)
with rasterio.open(cwd_fname, crs='EPSG:4326') as src:
     out_image, out_transform = mask(src, geoms, crop=True)
