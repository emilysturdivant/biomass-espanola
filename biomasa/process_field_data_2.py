# -*- coding: utf-8 -*-
"""
by: Emily Sturdivant, emilysturdivant@gmail.com
requires: python 3.6+

Pre-process and QC field inventory data.

OUTPUT: write CSVs of field data (haiti_biomass_v2_stacked.csv and haiti_biomass_v2_mplots.csv)

Python environment kernel: using py3_geo on Ubuntu and Python 3 on Mac
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
import re
import geopandas as gpd
import sys
try:
    proj_dir = os.path.dirname(os.path.realpath(__file__))
except:
    proj_dir = os.path.dirname(os.path.realpath('process_field_data_2.py'))
sys.path.append(proj_dir) # Add the script location to the system path just to make sure this works.
from biomasa.core import *

print(python_version())

#%%
# Set working directory
home = proj_dir #r'/Users/emilysturdivant/GitHub/biomass-espanola'
# home = r'/home/esturdivant/code/biomass-espanola' # work desktop

#%% Filenames
json_fname = os.path.join(home, 'data', 'standardize_creole.json')
lookup_fname = os.path.join(home, 'data', 'exploded_specieslookup.csv')
creole_wds_fname = os.path.join(home, 'data', 'creole_wooddensity_GWDBYavg_allSDs.csv')
field_data_fname = os.path.join(home, 'data', 'haiti_data_filled.csv')
plots_overview_fname = os.path.join(home, 'data', 'haiti_plots_meta.csv')
out_data = os.path.join(home, 'data', 'haiti_data_wds2.csv')
shp_allplots = os.path.join(home, 'data', 'AllPlots.shp')

#%% Load datasets created in other scripts
# Load creole-species lookup table
lookup_df = pd.read_csv(lookup_fname)
# Load alt_to_name dict - dict to fix typos in common names
alt_to_name = json.load(open(json_fname))
# Aggregated wood densities - creole to wood density lookup
creole_wds = pd.read_csv(creole_wds_fname)
# Field data
df_filled = pd.read_csv(field_data_fname)

#%% get plot centroids
# read file
plots_gdf.columns
plots_gdf = gpd.read_file(shp_allplots)[['plot_no', 'biomass_tf', 'geometry']] # WGS84 UTM 18N
# get area in hectares
plots_gdf['area_ha'] = plots_gdf.area*0.0001
# get lat/lon coordinates in separate columns
plots_gdf = plots_gdf.to_crs('EPSG:4326')
plots_gdf['lon'] = plots_gdf.centroid.x
plots_gdf['lat'] = plots_gdf.centroid.y

# join to mplots
mplots = pd.read_csv(plots_overview_fname)
mplots = mplots.join(plots_gdf.set_index('plot_no'), on='plot_no')\
            .drop(columns=['geometry'])
mplots.to_csv(os.path.join(home, 'data', 'mplots_geoms.csv'), index=False)

#%% Add plot data to field dataset
df_filled = df_filled.join(mplots.set_index('plot_no')[['lon', 'lat']], on='plot_no')

#%% Add wds to data
wds = creole_wds[['all_names', 'mean_gn', 'sd_gn']].rename(columns={'all_names':'creole', 'mean_gn':'meanWD', 'sd_gn':'sdWD'}).set_index('creole')
mstems = df_filled.join(wds, on='sp_creole')

# Fill Null values in meanWD with plot means
plot_wds = mstems.groupby('plot_no')['meanWD'].transform('median')
mstems['meanWD'] = mstems['meanWD'].fillna(plot_wds)
# Are there still NaN wood densities?
mstems['meanWD'].isna().sum()

# Export
mstems.to_csv(out_data, index=False)

#%% Look at specific sites for QC

mstems[mstems['plot_no']==16, :]
deluger2 = mstems.set_index('plot_no').loc[16, ['sp_creole', 'dbh_cm', 'ht_m', 'meanWD', 'sdWD']]
deluger2.describe()
deluger2[deluger2['dbh_cm']>200]
