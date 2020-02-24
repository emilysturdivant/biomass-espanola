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
# home = r'/home/esturdivant/code/biomass-espanola' # work desktop

#%% Filenames
out_wd_fname = os.path.join(home, 'data', 'bwayo_densities_2.csv')
plots_overview_fname = os.path.join(home, 'data', 'haiti_plots_meta.csv')
by_table_fname = os.path.join(home, 'data', 'exploded_specieslookup.csv')
out_filled_data_fname = os.path.join(home, 'data', 'haiti_data_filled.csv')
master_lookup = os.path.join(home, 'data', 'master_lookup_2.csv')

#%% Parse Global Wood Density and do cross-check
gwd_fname = os.path.join(home, 'data', 'GlobalWoodDensityDatabase.xlsx')
gwd_df = pd.read_excel(gwd_fname, sheet_name='Data', header=0,
    names=['gwd_num', 'family', 'binomial', 'wd',
        'region', 'gwd_ref_no'],
    index_col='gwd_num',
    converters={'binomial':lambda x : x.lower()})
gwd_df = split_species_binomial(gwd_df, binomial_fld='binomial') # in process_field_data_2

# Extract rows from GWD with species present in field data
lookup_field_species = pd.read_csv(master_lookup)
field_binoms = lookup_field_species[binomial_fld]
gwd_lookup = gwd_df.loc[gwd_df['binomial'].isin(field_binoms)]

# Print QC info
field_binoms = lookup_field_species[binomial_fld]
print(f'Unique species in field data: {len(field_binoms)}')
gwd_lookup = gwd_df.loc[gwd_df['binomial'].isin(field_binoms)]
gwd_unique_matches = gwd_lookup['binomial'].unique()
print(f'Matching species in GWD: {len(gwd_unique_matches)}')
unlisted = field_binoms[~field_binoms.isin(gwd_df['binomial'])]
print(f'Field species missing from GWD: {len(unlisted)}')

field_genus = pd.Series(lookup_field_species['genus'].unique())
print(f'Unique genus in field data: {len(field_genus)}')
gwd_lookup_genus = gwd_df.loc[gwd_df['genus'].isin(field_genus)]
gwd_unique_matches = gwd_lookup_genus['genus'].unique()
print(f'Matching genus in GWD: {len(gwd_unique_matches)}')
unlisted = field_genus[~field_genus.isin(gwd_df['genus'])]
print(f'Field genus missing from GWD: {len(unlisted)}')
unlisted

#%% Get mean Bwa Yo wood density at the three levels.
field_species = pd.read_csv(by_table_fname)
binomial_fld = 'species_binomial'
comm_name_fld = 'all_names'
by_wds = field_species[[comm_name_fld, 'family', 'genus', binomial_fld, 'wd_avg']]

# Species means
wd_agg = by_wds.groupby(binomial_fld)['wd_avg'].agg(mean='mean', sd='std')
# wd_agg = creole_wds.groupby(binomial_fld)['wd_avg'].mean().rename('by_spec_mean')
by_wds = by_wds.join(wd_agg, on=binomial_fld, rsuffix='_BYsp')
# Genus means
wd_agg = by_wds.groupby('genus')['wd_avg'].agg(mean='mean', sd='std', med_BYgn='median')
# wd_agg = creole_wds.groupby('genus')['wd_avg'].mean().rename('_by_gen')
by_wds = by_wds.join(wd_agg, on='genus', rsuffix='_BYgn')
# Family means
wd_agg = by_wds.groupby('family')['wd_avg'].agg(mean='mean', sd='std')
# wd_agg = creole_wds.groupby('family')['wd_avg'].mean().rename('by_fam_mean')
by_wds = by_wds.join(wd_agg, on='family', rsuffix='_BYfm')
creole_wds = by_wds.groupby(comm_name_fld).mean()

creole_wds = creole_wds.drop(['wd_avg'], axis=1).rename(columns={'mean':'mean_BYsp', 'sd':'sd_BYsp'})
creole_wds.loc[['pwa valye'],:]
creole_wds.describe()

#%% Import wood densities from R
# Initialize with my lookup table.
wds_fromR = field_species[[comm_name_fld, 'family', 'genus', 'species']]

# GWD+BY species, genus, and family
wds = pd.read_csv(os.path.join(home, 'getWoodDensity_SpecGenFam_BY.csv'))
wds = wds.groupby(['family', 'genus', 'species'])['meanWD', 'sdWD'].first().rename(columns={'meanWD':'mean_GWDBYspgnfm', 'sdWD': 'sd_GWDBYspgnfm'})
wds_fromR = wds_fromR.join(wds, on=['family', 'genus', 'species'])
# GWD+BY Genus and Family
wds = pd.read_csv(os.path.join(home, 'getWoodDensity_GenFam_BY.csv'))
wds = wds.groupby(['family', 'genus'])['meanWD', 'sdWD'].first().rename(columns={'meanWD':'mean_GWDBYgnfm', 'sdWD': 'sd_GWDBYgnfm'})
wds_fromR = wds_fromR.join(wds, on=['family', 'genus'])
# GWD+BY Genus
wds = pd.read_csv(os.path.join(home, 'getWoodDensity_Gen_BY.csv'))
wds = wds.groupby(['family', 'genus'])['meanWD', 'sdWD'].first().rename(columns={'meanWD':'mean_GWDBYgn', 'sdWD': 'sd_GWDBYgn'})
wds_fromR = wds_fromR.join(wds, on=['family', 'genus'])
# GWD species, genus, and family
wds = pd.read_csv(os.path.join(home, 'getWoodDensity_SpecGenFam.csv'))
wds = wds.groupby(['family', 'genus', 'species'])['meanWD', 'sdWD'].first().rename(columns={'meanWD':'mean_GWDspgnfm', 'sdWD': 'sd_GWDspgnfm'})
wds_fromR = wds_fromR.join(wds, on=['family', 'genus', 'species'])
# GWD genus, and family
wds = pd.read_csv(os.path.join(home, 'getWoodDensity_GenFam.csv'))
wds = wds.groupby(['family', 'genus'])['meanWD', 'sdWD'].first().rename(columns={'meanWD':'mean_GWDgnfm', 'sdWD': 'sd_GWDgnfm'})
wds_fromR = wds_fromR.join(wds, on=['family', 'genus'])
# GWD+BY Genus and Family with CentralAmericaTrop
wds = pd.read_csv(os.path.join(home, 'getWoodDensity_GenFam_BY_CAT.csv'))
wds = wds.groupby(['family', 'genus'])['meanWD', 'sdWD'].first().rename(columns={'meanWD':'mean_GWDBYgnfm_CAT', 'sdWD': 'sd_GWDBYgnfm_CAT'})
wds_fromR = wds_fromR.join(wds, on=['family', 'genus'])
# Get mean for each creole name
creole_wds = creole_wds.join(wds_fromR.groupby(comm_name_fld).mean(), on=comm_name_fld)
creole_wds.to_csv(os.path.join(home, 'data', 'options_for_creole_wooddensity.csv'))

#%% Create my own GWD-BY combo by

#%% Look at relationships between columns
creole_wds_means = creole_wds.filter(like='mean')

# Correlations
creole_wds_means.corr()

# Scatter matrix
fig, ax = plt.subplots(figsize=(12,12))
pd.plotting.scatter_matrix(creole_wds_means, alpha=1, ax=ax)
fig.savefig(os.path.join(home, 'qc_plots', 'wooddensities_scatter_matrix.png'))

# Correlations visualized
fig = plt.figure(figsize=(9, 9))
plt.matshow(creole_wds_means.corr('pearson'), fignum=fig.number)
plt.xticks(range(creole_wds_means.shape[1]), creole_wds_means.columns, fontsize=14, rotation=45)
plt.yticks(range(creole_wds_means.shape[1]), creole_wds_means.columns, fontsize=14)
cb = plt.colorbar()
cb.ax.tick_params(labelsize=14);
fig.savefig(os.path.join(home, 'qc_plots', 'wooddensities_corrmatrix.png'))

# Correlations visualized
fig = plt.figure(figsize=(9, 9))
plt.matshow(creole_wds_means.corr('spearman'), fignum=fig.number)
plt.xticks(range(creole_wds_means.shape[1]), creole_wds_means.columns, fontsize=14, rotation=45)
plt.yticks(range(creole_wds_means.shape[1]), creole_wds_means.columns, fontsize=14)
cb = plt.colorbar()
cb.ax.tick_params(labelsize=14);
fig.savefig(os.path.join(home, 'qc_plots', 'wooddensities_corrmatrix_spearman.png'))

#%% QC
# Print QC info
creole_wds_means = creole_wds.filter(like='mean')

# How many NaNs in each column?
creole_wds.isna().sum()
len(creole_wds)

# Look at means (and numbers of NaNs)
nanct = creole_wds_means.isna().sum()
nanct.name = 'NaNs'
desc = creole_wds_means.describe(percentiles=[]).append(nanct).transpose()
desc

# Look at standard deviations
creole_wds_sds = creole_wds.filter(like='sd')
nanct = creole_wds_sds.isna().sum()
nanct.name = 'nan_count'
desc = creole_wds_sds.describe(percentiles=[]).append(nanct).transpose()
desc
desc4fig = desc.drop(['count', 'nan_count'], axis=1)

fig = plt.figure(figsize=(5, 4))
plt.matshow(desc4fig, fignum=fig.number)
plt.xticks(range(desc4fig.shape[1]), desc4fig.columns, fontsize=14)
plt.yticks(range(len(desc4fig)), desc4fig.index, fontsize=14)
cb = plt.colorbar()
cb.ax.tick_params(labelsize=14);
fig.savefig(os.path.join(home, 'qc_plots', 'wooddensities_sd_summary.png'))

creole_wds[creole_wds['mean_BYsp'].isna()].index
creole_wds[creole_wds['mean_BYgn'].isna()].index
creole_wds[creole_wds['mean_GWDBYspgnfm'].isna()].index
creole_wds[creole_wds['mean_GWDBYgnfm'].isna()].index
creole_wds[creole_wds['mean_GWDBYgn'].isna()].index
creole_wds[creole_wds['mean_GWDspgnfm'].isna()].index

creole_wds_means[creole_wds_means['mean_GWDBYgn'].isna()].index
# Why do damari and dalmari have Null values when the binomial has a BY wood density?
'dalmari' in creole_wds_means[creole_wds_means['mean_BYgn'].isna()].index
'abbe marron' in creole_wds_means[creole_wds_means['mean_BYgn'].isna()].index
# Maybe because of what I export for use in R
# I import bwayo_densities_2.csv as the supplemental wood density data in getWoodDensity.
# Maybe it doesn't attempt to aggregate by genus for the supplemental table. 



#%% Make wood density array (join one to many) of same form as data df for use in computeAGB()
# Load pre-processed field data
df = pd.read_csv(out_filled_data_fname)
# Join
creole_wds_means = creole_wds.filter(like='mean')
df = df.join(creole_wds_means, on='sp_creole')
# Export
df.to_csv(os.path.join(home, 'data', 'mstems_with_wooddensities.csv'), index=False)

# QC
df.columns
# Isolate to just species and wood density values
df2 = df.dropna(subset=['sp_creole'], how='all')\
    .drop(['dbh_cm', 'ht_m', 'plot_fname', 'date', 'plot_id', 'plot_loc', 'data_scribe', 'plot_no', 'plot_shp', 'plot_area'], axis=1)\
    .drop_duplicates()
len(df2)

df2[df2['mean_BYsp'].isna()]
df2[df2['mean_BYsp'].isna()]['sp_creole']
df2[df2['mean_BYgn'].isna()]['sp_creole']
df2[df2['mean_BYfm'].isna()]['sp_creole']
df2[df2['mean_GWDBYspgnfm'].isna()]['sp_creole']
df2[df2['mean_GWDBYgnfm'].isna()]['sp_creole']
df2[df2['mean_GWDgnfm'].isna()]['sp_creole']
df2[df2['mean_GWDBYgn'].isna()]['sp_creole']



#%% computeAGB() from BIOMASS https://github.com/AMAP-dev/BIOMASS/blob/master/R/computeAGB.R
D = diam
genus
species
family = None
plot = None
height
latitude
longitude

# arguments
D = df['dbh_cm']
WD = df['mean_GWDBYgnfm']
H = None
coord = (19, -72)
Dlim = None

# Parameters verification

# Compute E
E =

# Modified Eq. 7 from Chave et al. 2014 Global Change Biology
AGB = exp(-2.023977 - 0.89563505 * E + 0.92023559 * log(WD) + 2.79495823 * log(D) - 0.04606298 * (log(D)^2)) / 1000







#%% Parse Global Wood Density - copied to above and possibly modified
gwd_fname = os.path.join(home, 'data', 'GlobalWoodDensityDatabase.xlsx')
gwd_df = pd.read_excel(gwd_fname, sheet_name='Data', header=0,
    names=['gwd_num', 'family', 'binomial', 'wd',
        'region', 'gwd_ref_no'],
    index_col='gwd_num',
    converters={'binomial':lambda x : x.lower()})
gwd_df = split_species_binomial(gwd_df, binomial_fld='binomial') # in process_field_data_2

# Extract rows from GWD with species present in field data
lookup_field_species = pd.read_csv(master_lookup)
field_binoms = lookup_field_species[binomial_fld]
gwd_lookup = gwd_df.loc[gwd_df['binomial'].isin(field_binoms)]

#%% Replicate getWoodDensity from BIOMASS - using World as region and including the Bwa Yo wood densities
# def getWoodDensity(gwd_df, genus, species, stand=None, family=None, xtra_wd_data=None, verbose=True):
xtra_wd_data = field_species
xtra_wd_data = xtra_wd_data[['creole', binomial_fld, 'genus', 'species_abbr', 'family', 'wd_avg']].rename(columns={'wd_avg':'wd', binomial_fld: 'binomial', 'species_abbr':'species'})
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