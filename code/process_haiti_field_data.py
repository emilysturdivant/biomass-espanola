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
import geopandas as gpd

print(python_version())

#%% Functions
'''
FUNCTIONS
'''
def strip_accents(text):
    """
    Strip accents from input String.

    :param text: The input string.
    :type text: String.

    :returns: The processed String.
    :rtype: String.
    from: https://stackoverflow.com/a/31607735
    """
    try:
        text = unicode(text, 'utf-8')
    except (TypeError, NameError): # unicode is a default on python 3
        pass
    text = unicodedata.normalize('NFD', text)
    text = text.encode('ascii', 'ignore')
    text = text.decode("utf-8")
    return str(text)

def get_plot_data(data_fname, col_ints, verbose=True):
    # Get data for plot indicated by col_ints
    df = pd.read_excel(data_fname, 'Plots', skiprows=[0,1,2], header=0,
        usecols=col_ints, names=['sp_creole','dbh_cm','ht_m','ba_m2'], #dtype={'ht_m':'float'},
        skip_blank_lines = True,
        converters={'sp_creole':lambda x : strip_accents(x.strip().lower()), 'ht_m': lambda x: str(x).strip("''")})
    df = df.dropna(axis=0, how='all').astype({'ht_m':'float'})
    # Create padding row in plots without trees
    if not len(df):
        df = df.reindex(labels=[1], fill_value=np.nan)
    # Get plot number, shapefile name, and area for given plot
    plot = pd.read_excel(data_fname, 'Plots', header=0, usecols=col_ints, nrows=1)
    plot_no=plot.columns[0].split('#')[1]
    if verbose:
        print(f'Plot number: {plot_no} | {plot.iloc[0,0]} | Area: {plot.iloc[0,2]} ha')

    df = df.assign(plot_no=plot_no,
                plot_shp=plot.iloc[0,0],
                plot_area=plot.iloc[0,2])
    return(df)

#%% Initialize
'''
Initialize
'''
# Set working directory
home = r'/Users/emilysturdivant/GitHub/biomass-espanola' # laptop
# home = r'/home/esturdivant/code/biomass-espanola' # work desktop

# Filenames
data_fname = os.path.join(home, 'data', 'haiti_biomass_v2.xlsx')
json_fname = os.path.join(home, 'standardize_creole.json')
lookup_fname = os.path.join(home, 'lookup_famAvgWD_byCreole.xlsx')
zstats_csv = os.path.join(home, 'plots_zstats_07gamma0_qgis.csv')

# Load creole name replacement dictionary
alt_to_name = json.load(open(json_fname))

#%% Import values from output of QGIS zonal stats
plots_gam0 = pd.read_csv(os.path.join(home, zstats_csv))
plots_gam0.loc[:, ['Name']]
plots_gam0.Name



#%% Load field data and convert to table of all plots stacked
'''
Convert field data to table of all plots stacked
'''
# Iteratively load each plot in turn and concatenate
col_init = [0,1,2,3]
df = get_plot_data(data_fname, col_init)
for i in range(1,36):
    col_ints = np.add(col_init, 4*i)
    df2 = get_plot_data(data_fname, col_ints, verbose=False)
    df = pd.concat([df, df2], ignore_index=True)

# Standardize creole names
df.loc[:, 'sp_creole'] = df.sp_creole.replace(alt_to_name)
df.columns

#%% Perform QC on data

#%% Check DBH measurements - look for outliers
df.dbh_cm
fig, ax = plt.subplots(figsize=(16,8))
ax.scatter(df['dbh_cm'], df['plot_no'])
ax.set_xlabel('DBH')
ax.set_ylabel('Plot Number')
plt.show()
df[df['dbh_cm'] < 1]

df.boxplot(column='dbh_cm', by='plot_no', figsize=(16,8))


fig, ax = plt.subplots(figsize=(16,8))
ax.scatter(df['dbh_cm'], df['ht_m'])
ax.set_xlabel('DBH')
ax.set_ylabel('Height')
plt.show()


fig, ax = plt.subplots(figsize=(16,8))
ax.scatter(df['plot_no'], df['ht_m'])
ax.set_xlabel('Plot')
ax.set_ylabel('Height')
plt.show()
df[df['ht_m'] > 60]

#%% Identify duplicate records
# Look at consecutive duplicated entries
cols = ["sp_creole", "dbh_cm", "plot_no"]
dups1 = (df[cols].shift() == df[cols]).all(axis=1).replace({False: np.nan})
de_dup = df[cols].loc[dups1.fillna(dups1.shift(-1)).fillna(False)]
de_dup

# All duplicates
df[df.duplicated(subset=['sp_creole', 'dbh_cm', 'ht_m', 'plot_no'], keep=False)]
# Duplicates within a plot (plot 9)
df_plt9 = df[df['plot_no'] == '9']
df_plt9[df_plt9.duplicated(subset=['sp_creole', 'dbh_cm', 'ht_m'], keep=False)]
df_plt9_mango = df_plt9[df_plt9['sp_creole'] == 'mango']
df_plt9_mango[df_plt9_mango.duplicated(subset=['dbh_cm', 'ht_m'], keep=False)]

#%% Look at DBH and height ranges by species
df.boxplot(column='dbh_cm', by='sp_creole', figsize=(20,8), rot=90)
df.boxplot(column='ht_m', by='sp_creole', figsize=(20,8), rot=90)

# Convert unknown species to something standardized


# Write to CSV
df.to_csv(os.path.join(home, 'haiti_biomass_v2_stacked.csv'), index=False)

#%% Basal Area: the area of land that is occupied by the cross-section of a tree. (silvR)
# BA = π(DBH/2)^2

#%% Stem volume: an estimate of the ammount of space that the tree bole occupies. (silvR)
# Calculate with equation from unknown source


#%%
'''

'''
# Read in lookup table (relates common creole name to wood density)
dens_lookup = pd.read_excel(lookup_fname)
dfdens = df.join(dens_lookup.set_index('creole'), on='sp_creole', how='left')
dfdens.dtypes

dfdens

#%% Perform allometric calculation for plot

# Use Model 4 from Chave et al. (2014) when we have height values
# select rows with height values
trees_wHt = dfdens.loc[~dfdens['ht_m'].isna(), :]
agb_wHt = lambda x: 0.0673 * (x['gwd_density'] * x['dbh_cm']**2 * x['ht_m'])**0.976
# trees_wHt.loc[:, ['biomass']] = trees_wHt.apply(agb_wHt, axis=1)
biomass_wHt = trees_wHt.assign(agb=agb_wHt)
biomass_wHt


trees_noHt = dfdens[dfdens['ht_m'].isna()]
agb_noHt =
biomass_noHt = trees_noHt.assign(biomass=agb_noHt)
biomass_noHt












#%% Get value for E (combination of TS, CWD, and PS) for Hispaniola
# Downloaded CWD from http://chave.ups-tlse.fr/pantropical_allometry.htm
# Downloaded TS and PS from http://worldclim.org/version2
# Working through this https://www.datacamp.com/community/tutorials/geospatial-data-python

# Load rasters
cwd_fname = os.path.join(home, 'data', 'CWD.tif')
ts_fname = os.path.join(home, 'data', 'wc2.0_bio_2.5m_04.tif')
ps_fname = os.path.join(home, 'data', 'wc2.0_bio_2.5m_15.tif')
cwd_world = cwd_fname
ts_world = ts_fname
ps_world = ps_fname

# Isolate to hispaniola
# Load and merge haiti and DR shapefiles

# filenames
haiti_fname = os.path.join(home, 'data', 'HTI_adm0.shp')
dr_fname = os.path.join(home, 'data', 'DOM_adm0.shp')

# import fiona
#
# # Load and check data
# haiti = fiona.open(haiti_fname)
# print(haiti.schema)
# print(haiti.next()) # (GeoJSON format) {'geometry': {'type': 'LineString', 'coordinates': [(0.0, 0.0), (25.0, 10.0), (50.0, 50.0)]}, 'type': 'Feature', 'id': '0', 'properties': OrderedDict([(u'FID', 0.0)])}
# dr = fiona.open(dr_fname)
# print(dr.schema) # {'geometry': 'LineString', 'properties': OrderedDict([(u'FID', 'float:11')])}
# #first feature of the shapefile
# print(dr.next())



import fiona

haiti = gpd.read_file(haiti_fname)
haiti.head()
haiti.columns
haiti.plot()
dr = gpd.read_file(dr_fname)
print(dr.head())
dr.describe

# Merge
hispaniola = gpd.sjoin(haiti, dr, how='inner', op='intersects')
hispaniola
hispaniola.plot()

hispaniola.dissolve(by='EU')


# Mask
cwd_hisp =
ts_hisp =
ps_hisp =

# Calculate E
E = 1.e-3 * (0.178*ts_hisp - 0.938*cwd_hisp - 6.61*ps_hisp)


agb_noHt = exp(-1.803 - 0.976*E + 0.976*ln(wood_dens) + 2.673 ln(dbh) - 0.0299 * ln(dbh)**2)
