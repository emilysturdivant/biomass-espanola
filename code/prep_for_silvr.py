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

#%% Functions
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

def split_species_binomial(df, binomial_fld='by_binomial'):
    def get_second_word(x):
        ser = []
        for lst in x:
            if len(lst) > 1:
                ser += [lst[1].strip()]
            else:
                ser += ['']
        return(ser)
    binomials_split_ser = df[binomial_fld].str.split(' ')
    df = df.assign(
        genus=[i[0].capitalize() for i in binomials_split_ser],
        species=[' '.join(i[1:]).strip() for i in binomials_split_ser],
        species_abbr=get_second_word(binomials_split_ser)
        )
    return(df)

def get_plot_data(data_fname, col_ints, verbose=True):
    # Get data for plot indicated by col_ints
    df = pd.read_excel(data_fname, 'Plots', skiprows=[0,1,2], header=0,
        usecols=col_ints, names=['sp_creole','dbh_cm','ht_m','ba_m2'], #dtype={'ht_m':'float'},
        skip_blank_lines = True, na_values='0',
        converters={'sp_creole':lambda x : strip_accents(x.strip().lower()), 'ht_m': lambda x: str(x).strip("''")})
    df = df.dropna(axis=0, how='all').astype({'ht_m':'float'})
    # Create padding row in plots without trees
    if not len(df):
        df = df.reindex(labels=[1], fill_value=np.nan)
    # Get plot number, shapefile name, and area for given plot
    plot = pd.read_excel(data_fname, 'Plots', header=0, usecols=col_ints, nrows=1)
    df['plot_no'] = plot.columns[0].split('#')[1]
    if verbose:
        print(f'Plot number: {df['plot_no']} | {plot.iloc[0,0]} | Area: {plot.iloc[0,2]} ha')
    df = df.assign(plot_shp=plot.iloc[0,0], plot_area=plot.iloc[0,2])
    return(df)

#%%
# Set working directory
home = r'/Users/emilysturdivant/GitHub/biomass-espanola'
home = r'/home/esturdivant/code/biomass-espanola' # work desktop

#%% Work with field data - Look at species in all plots
data_fname = os.path.join(home, 'data', 'haiti_biomass_v2.xlsx')
gwd_fname = os.path.join(home, 'data', 'GlobalWoodDensityDatabase.xlsx')
by_fname = os.path.join(home, 'data', 'bwayo_species.xlsx')
json_fname = os.path.join(home, 'standardize_creole.json')
zstats_csv = os.path.join(home, 'plots_zstats_07gamma0_qgis.csv')

#%% Load species table digitized from Bwa Yo and split binomial into genus and species
by_df = pd.read_excel(by_fname, header=0, usecols=[0,1,2,3],
    names=['by_binomial', 'creole', 'BY_spec_grav', 'family'],
    converters={'by_binomial':lambda x : x.lower(),
        'family':lambda x : x.split(' (')[0].lower().capitalize()})
by_df.replace({'Capparaceae': 'Brassicaceae', 'Sterculiaceae': 'Malvaceae'}, inplace=True)

# Split species binomial into genus, species, and species_abbr
by_df = split_species_binomial(by_df, binomial_fld='by_binomial')

#%% Create supplemental wood density from all Bwa Yo values
bwayo_wd = by_df.loc[~by_df['BY_spec_grav'].isna(),
    ['genus', 'species', 'species_abbr', 'BY_spec_grav', 'family']]
# Replace spp. with NaN
bwayo_wd = bwayo_wd.assign(species_abbr=bwayo_wd['species_abbr'].replace('spp.', np.nan))
# convert WD range to mean
bwayo_wd = bwayo_wd.assign(wd_avg=[np.mean([float(i) for i in
    re.findall(r"[0-9.]+", str(s))]) for s in bwayo_wd['BY_spec_grav']])

# Export to CSV
bwayo_wd_for_export = bwayo_wd[['genus', 'species', 'wd_avg']].rename(columns={'wd_avg':'wd'})
bwayo_wd_for_export.to_csv(os.path.join(home, 'bwayo_densities.csv'), index=False)

# Join wood density averages to the Bwa Yo DF.
by_df = by_df.join(bwayo_wd['wd_avg'])

#%% Explode DF by the creole names column. For every row with multiple creole names, duplicate species row.
creole_df = (by_df
            .assign(creole=by_df.creole.str.split(','))
            .explode('creole')
            .reset_index(drop=True)
        )
creole_df['creole'] = creole_df.creole.str.strip()
by_creole_names = creole_df.creole.unique().tolist()

#%% Extract species in field data from BY df
# Load creole name replacement dictionary
alt_to_name = json.load(open(json_fname))

# Create series of all species columns (labeled 'sp')
spec_df = pd.read_excel(data_fname, 'Plots', skiprows=[0,1,2],
                        usecols=lambda x : x.startswith('sp'))
spec_ser = spec_df.stack().apply(lambda x : strip_accents(x).strip().lower()).replace(alt_to_name)

# Extract species in field data from BY df
field_species_uniq = pd.Series(spec_ser.unique())
field_species = by_df.loc[by_df['creole'].isin(field_species_uniq)].reset_index(drop=True)

#%% Create creole to species_name lookup and species_name to bwayo_wd lookup
lookup_genus = field_species[['genus', 'creole']].groupby('creole').first()
lookup_genus
|#%% Load field data and convert to table of all plots stacked
# Iteratively load each plot in turn and concatenate
col_init = [0,1,2,3]
df = get_plot_data(data_fname, col_init, verbose=False)
for i in range(1,36):
    col_ints = np.add(col_init, 4*i)
    df2 = get_plot_data(data_fname, col_ints, verbose=False)
    df = pd.concat([df, df2], ignore_index=True)

# Standardize creole names
df.loc[:, 'sp_creole'] = df.sp_creole.replace(alt_to_name)
df.loc[60:70, :]

# Replace np.nan with 0 in dbh_cm column
df.replace({np.nan:0})

#%% Make and export plots DF
mplots = df[['plot_no', 'plot_shp', 'plot_area']].groupby('plot_no').first()
mplots.to_csv(os.path.join(home, 'data', 'haiti_biomass_v2_mplots.csv'), index=True)
mplots

#%% Join genus to field data DF
df = df.join(lookup_genus, on='sp_creole', how='left')
df.loc[60:70, :]
df.to_csv(os.path.join(home, 'mstems_genus_rough.csv'), index=False)
