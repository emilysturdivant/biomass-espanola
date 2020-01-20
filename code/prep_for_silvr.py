# -*- coding: utf-8 -*-
"""
by: Emily Sturdivant, emilysturdivant@gmail.com
requires: python 3.6+

OVERVIEW: Getting started with processing Haiti field data
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

print(python_version())

#%%
# Set working directory
home = r'/Users/emilysturdivant/GitHub/biomass-espanola'
home = r'/home/esturdivant/code/biomass-espanola' # work desktop

|#%% Load pre-processed field data
df = pd.read_csv(os.path.join(home, 'data', 'haiti_biomass_v2_stacked.csv'))
df.head()

#%% Create creole to species_name lookup and species_name to bwayo_wd lookup
field_species = pd.read_csv(os.path.join(home, 'data', 'master_lookup.csv'))

# Simplistic: take the first genus matching the creole name
lookup_genus = field_species[['genus', 'creole']].groupby('creole').first()
lookup_genus

'''
TODO: investigate sensitivity to different creole to species matches
'''

#%% Join genus to field data DF
df = df.join(lookup_genus, on='sp_creole', how='left')
df.loc[60:70, :]
df.to_csv(os.path.join(home, 'mstems_genus_rough.csv'), index=False)
