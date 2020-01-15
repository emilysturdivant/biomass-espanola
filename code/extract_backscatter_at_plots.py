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
home = r'/home/esturdivant/code/biomass-espanola'

#%% Filenames

#%% Import values from output of QGIS zonal stats
plots_zstats = pd.read_csv(os.path.join(home, 'plots_zstats_07gamma0_qgis.csv'))




plots_zstats.plot.scatter(x='2007_count', y='2007_mean')






















#%% Try extraction in python...
# Filenames
cwd_fname = r'/home/esturdivant/Documents/biota_out/Gamma0_2010_haiti_HVdb/Gamma0_2010_N19W073.tif'
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
with rasterio.open(cwd_fname, crs='EPSG:4326') as src:
    geom = rasterio.warp.transform_geom(dataset.crs, 'EPSG:4326', geom, precision=6)
    out_image, out_transform = mask(src, geoms, crop=True)
