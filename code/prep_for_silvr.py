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
df.loc[60:70, :]

#%% Create creole to species_name lookup and species_name to bwayo_wd lookup
field_species = pd.read_csv(os.path.join(home, 'data', 'master_lookup.csv'))

# GroupBy creole name
grouped = field_species.groupby('creole')#, as_index=False)

# Simplistic: take the first genus matching the creole name
grouped.first().head()
grouped.first()[['genus']].describe()
grouped.first()[['family', 'genus', 'species']].describe()
# Attempts at modes...
# First, get mode. But if there is no mode, then a list is returned.
lookup_modes = grouped.agg(lambda x: x.mode())
lookup_modes.loc['bwa savann', :]
# If there is no top genus and there is only one family, set genus to NaN to use the family.
for row in lookup_modes.iterrows():
    if not isinstance(row[1].genus, str) and isinstance(row[1].family, str):
        print(f'{row[0]}: {row[1].genus} --> {row[1].family}')
        lookup_modes.loc[row[0], 'genus'] = np.nan
lookup_modes.loc['bwa dom', :]

# There are 7 creole names for which there is no top genus.
for row in lookup_modes.iterrows():
    if not isinstance(row[1].genus, str):
        print(f'{row[0]}: {row[1].genus}')
multi_geni = [row[0] for row in lookup_modes.iterrows() if not isinstance(row[1].genus, str)]
lookup_modes.loc[multi_geni, :]
lookup_modes.loc['bwa blan', :]
field_species[field_species['creole'] == 'bwa blan']

np.median(row[1].wd_avg)


grouped['genus'].agg(mode=lambda x: x.mode()).head()
grouped.get_group('bwa savann')

grouped.genus.agg(mode=lambda x: x.mode())
grouped.family.agg(mode=lambda x: x.mode())
grouped.by_binomial.agg(mode=lambda x: x.mode()).head()
grouped.agg(
    top_genus=('genus', lambda x: x.mode()),
    top_binomial=('by_binomial', lambda x: x.mode()),
    top_species=('species', lambda x: x.mode())
)


# Below here is not yet working:
grouped.agg(mode=lambda x: x.mode()).head()
grouped[['genus', 'creole']].agg(mode=lambda x: x.mode())

animals.groupby("kind").agg(
min_height=pd.NamedAgg(column='height', aggfunc='min'),
max_height=pd.NamedAgg(column='height', aggfunc='max'),
average_weight=pd.NamedAgg(column='weight', aggfunc=np.mean),
)


'''
TODO: investigate sensitivity to different creole to species matches
'''

#%% Join genus to field data DF
df = df.join(lookup_genus, on='sp_creole', how='left')
df.loc[60:70, :]
df.to_csv(os.path.join(home, 'mstems_FamGenSpec_lookupfirst.csv'), index=False)




#%% Parse Global Wood Density
gwd_fname = os.path.join(home, 'data', 'GlobalWoodDensityDatabase.xlsx')
gwd_df = pd.read_excel(gwd_fname, sheet_name='Data', header=0,
    names=['gwd_num', 'family', 'binomial', 'wd',
        'region', 'gwd_ref_no'],
    index_col='gwd_num',
    converters={'binomial':lambda x : x.lower()})

# Extract rows from GWD with species present in field data
field_binoms = field_species['by_binomial']
gwd_lookup = gwd_df.loc[gwd_df['binomial'].isin(field_binoms)]

# Print QC info
print(f'Unique species in field data: {len(field_binoms)}')
gwd_unique_matches = gwd_lookup['gwd_binomial'].unique()
print(f'Matching species in GWD: {len(gwd_unique_matches)}')
unlisted = field_binoms[~field_binoms.isin(gwd_df['gwd_binomial'])]
print(f'Field species missing from GWD: {len(unlisted)}')


#%% Get mean Bwa Yo wood density at the three levels.
creole_wds = field_species[['creole', 'family', 'genus', 'by_binomial', 'wd_avg']]
BY_fam_mean = creole_wds.groupby('family')['wd_avg'].mean().rename('by_fam_mean')
creole_wds = creole_wds.join(BY_fam_mean, on='family')
BY_gen_mean = creole_wds.groupby('genus')['wd_avg'].mean().rename('by_gen_mean')
creole_wds = creole_wds.join(BY_gen_mean, on='genus')
BY_spec_mean = creole_wds.groupby('by_binomial')['wd_avg'].mean().rename('by_spec_mean')
creole_wds = creole_wds.join(BY_spec_mean, on='by_binomial')
creole_wds = creole_wds.groupby('creole').mean()
creole_wds
any(creole_wds['by_spec_mean'].isna())
any(creole_wds['by_gen_mean'].isna())
any(creole_wds['by_fam_mean'].isna())

creole_wds.plot.scatter('by_spec_mean', 'by_fam_mean')
creole_wds.plot.scatter('by_spec_mean', 'by_gen_mean')
ax1 = creole_wds.hist('by_spec_mean', alpha=0.4)
creole_wds.hist('by_gen_mean', ax=ax1, alpha=0.4)
creole_wds.hist('by_fam_mean', ax=ax1, alpha=0.4)
plt.show()
#%% Import wood densities from R
creole_wds2 = field_species[['creole', 'family', 'genus', 'species']]

# Get unique values and join to creole in field_species
wds = pd.read_csv(os.path.join(home, 'getWoodDensity_SpecGenFam_BY.csv'))
wds = wds.groupby(['family', 'genus', 'species'])['meanWD'].first().rename('gwdby_spgnfm_mean')
creole_wds2 = creole_wds2.join(wds, on=['family', 'genus', 'species'])
# Get unique values and join to creole in field_species
wds = pd.read_csv(os.path.join(home, 'getWoodDensity_GenFam_BY.csv'))
wds = wds.groupby(['family', 'genus'])['meanWD'].first().rename('gwdby_gnfm_mean')
creole_wds2 = creole_wds2.join(wds, on=['family', 'genus'])
# Get unique values and join to creole in field_species
wds = pd.read_csv(os.path.join(home, 'getWoodDensity_Gen_BY.csv'))
wds = wds.groupby(['family', 'genus'])['meanWD'].first().rename('gwdby_gn_mean')
creole_wds2 = creole_wds2.join(wds, on=['family', 'genus'])
# Get mean for each creole name
creole_wds2 = creole_wds2.groupby('creole').mean()
creole_wds = creole_wds.join(creole_wds2, on='creole')

#
creole_wds.plot.scatter('by_spec_mean', 'gwdby_gnfm_mean')
creole_wds.plot.scatter('gwdby_spgnfm_mean', 'gwdby_gnfm_mean')
creole_wds.plot.scatter('gwdby_spgnfm_mean', 'gwdby_gn_mean')

ax1 = mean_creole_wds.hist('wd_avg', alpha=0.5)
mean_creole_wds.hist('meanWD', ax=ax1, alpha=0.5)
plt.show()


#%% Replicate getWoodDensity from BIOMASS - using World as region and including the Bwa Yo wood densities
# def getWoodDensity(gwd_df, genus, species, stand=None, family=None, xtra_wd_data=None, verbose=True):
xtra_wd_data = field_species
xtra_wd_data = xtra_wd_data[['creole', 'by_binomial', 'genus', 'species_abbr', 'family', 'wd_avg']].rename(columns={'wd_avg':'wd', 'by_binomial': 'binomial', 'species_abbr':'species'})
# load GWD - we're not going to subset by region because there's not enough Central America Tropics

# merge Bwa Yo WD with GWD, joining by family, genus, and species, (outer join?)
# Drop rows without wood density values.
xtra_wd_data.dropna(subset=['wd'], inplace=True)
wd_data = pd.concat([gwd_df, xtra_wd_data], sort=True)

# Match creole to GWD values by binomial
xtra_wd_data[['binomial', 'creole']]
gwd_df[gwd_df['binomial']=='prosopis juliflora']
gwd_df2 = gwd_df.join(xtra_wd_data[['binomial', 'creole']].set_index('binomial'), on='binomial')
gwd_df2[gwd_df2['binomial']=='simarouba spp.']

# Load input data (mstems) - get genus, species, family, stand - only genus is required
# get a list of taxa, get all the unique combinations of family, genus, species
# coalesce replaces NaN values in one column with values in other column

# Select only the taxa that match taxa in field data.
# Check for matches between family and WD data.

# Compute mean and standard deviation at subsequent levels...
# Species level
meanSP = wd_data.groupby(['family', 'genus', 'species'])['wd'].agg(
    meanWDsp='mean',
    nIndsp='count',
    sdWDsp=lambda x: x.std() if x.count() > 10 else np.nan)
# Add values [meanWD, nInd, sdWD] to inputData
# Genus level
meanGN = wd_data.groupby(['family', 'genus'])['wd'].agg(
    meanWDgn='mean',
    nIndgn='count',
    sdWDgn=lambda x: x.std() if x.count() > 10 else np.nan)
# Add values to inputData where species left NaNs. - But all species have meanWD values...
# Family level
meanFM = wd_data.groupby(['family'])['wd'].agg(
    meanWDfm='mean',
    nIndfm='count',
    sdWDfm=lambda x: x.std() if x.count() > 10 else np.nan)
# Add values to inputData where species and genus left NaNs.
# Plot level
wd_data.join()
fam_means = field_species.join(fam_mean_density, on='family', how='left', lsuffix='_gwd')

wdmeans = meanSP.rename(columns={'meanWDsp': 'meanWD', 'nIndsp': 'nInd', 'sdWDsp': 'sdWD'})
wdmeans = wdmeans.reset_index()
wdmeans[~wdmeans['meanWD'].isna()]
wdmeans['meanWD'].isna().count()
wdmeans.iloc[40:45, :]
wdmeans.meanWD.fillna(meanGN.meanWDgn)


#%% Create lookup with family averaged values
fam_mean_density = gwd_df.groupby('gwd_family')['gwd_density'].mean()
fam_means = field_species.join(fam_mean_density, on='family', how='left', lsuffix='_gwd')
# Write to lookup_famAvgWD_allRegions.xlsx: mean WD for each family (all regions)
# fieldspec_byfam.to_excel(os.path.join(home, 'lookup_famAvgWD_allRegions.xlsx'), index=False)

# For each creole name, get mean WD for all matching species
creole_dens = fam_means.groupby('creole')['gwd_density'].mean()
# Write lookup_famAvgWD_byCreole.xlsx: mean family WD (all regions) of all families matching the creole name (usually only one match)
creole_dens.to_excel(os.path.join(home, 'lookup_famAvgWD_byCreole.xlsx'), index=True)
