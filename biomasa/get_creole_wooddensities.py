# -*- coding: utf-8 -*-
"""
by: Emily Sturdivant, emilysturdivant@gmail.com
requires: python 3.6+

OVERVIEW: Calculate wood densities for creole names based on Global Wood Density Database (GWD) and Bwa Yo (BY). Includes multiple methods for aggregating the wood densities.

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
home = r'/Users/emilysturdivant/GitHub/biomass-espanola'
# home = r'/home/esturdivant/code/biomass-espanola' # work desktop

#%% Filenames
by_wd_fname = os.path.join(home, 'data', 'bwayo_densities_wFam.csv')
plots_overview_fname = os.path.join(home, 'data', 'haiti_plots_meta.csv')
by_table_fname = os.path.join(home, 'data', 'exploded_specieslookup.csv')
out_filled_data_fname = os.path.join(home, 'data', 'haiti_data_filled.csv')
master_lookup = os.path.join(home, 'data', 'master_lookup_2.csv')
out_creole_wds_GWDBYavg = os.path.join(home, 'data', 'creole_wooddensity_GWDBYavg_allSDs.csv')

#%% Import
field_species = pd.read_csv(by_table_fname).rename(columns={'species_binomial': 'binomial'})
comm_name_fld = 'all_names'
creole_lookup = field_species[[comm_name_fld, 'binomial', 'genus', 'family']]

#%% Parse Global Wood Density and do cross-check
gwd_fname = os.path.join(home, 'data', 'GlobalWoodDensityDatabase.xlsx')
gwd_df = pd.read_excel(gwd_fname, sheet_name='Data', header=0,
    names=['gwd_num', 'family', 'binomial', 'wd',
        'region', 'gwd_ref_no'],
    index_col='gwd_num',
    converters={'binomial':lambda x : x.lower()})
gwd_df = split_species_binomial(gwd_df, binomial_fld='binomial') # in process_field_data_2
gwd_df['wd'].describe()

# Extract rows from GWD with species present in field data
lookup_field_species = pd.read_csv(master_lookup).rename(columns={'species_binomial': 'binomial'})
binom_fld = 'binomial'
field_binoms = lookup_field_species[binom_fld]
gwd_lookup = gwd_df.loc[gwd_df['binomial'].isin(field_binoms)]

#%% Print QC info for GWD
field_binoms = lookup_field_species[binom_fld]
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

#%% Append BwaYo to GWD and get mean WDs
# Get BwaYo densities
by_df = pd.read_csv(by_wd_fname).rename(columns={'species_binomial': 'binomial'})

# Concatenate and calculate means
gwdby_df = pd.concat([gwd_df, by_df], sort=False)\
            .drop(['region', 'gwd_ref_no', 'species_extd'], axis=1)

# Creole wood densities from BY
by_wds = get_mean_WDs_2(df=by_df, lookup_df=creole_lookup, binom_fld='binomial', wd_fld='wd', comm_name_fld='all_names')
# Creole wood densities from GWD
gwd_wds = get_mean_WDs_2(df=gwd_df, lookup_df=creole_lookup, binom_fld='binomial', wd_fld='wd', comm_name_fld='all_names')
# Combined
gwdby_wds = get_mean_WDs_2(df=gwdby_df, lookup_df=creole_lookup, binom_fld='binomial', wd_fld='wd', comm_name_fld='all_names')
# Combined with BY weighted equal to all GWD
gwdby_wds2 = get_mean_wds_3(gwd_df, by_df, creole_lookup, binom_fld, comm_name_fld)

# QC
# Look at means (and numbers of NaNs)
nanct = by_wds.isna().sum()
nanct.name = 'NaNs'
desc = by_wds.describe(percentiles=[]).append(nanct).transpose()
desc

# Look at means (and numbers of NaNs)
nanct = gwd_wds.isna().sum()
nanct.name = 'NaNs'
desc = gwd_wds.describe(percentiles=[]).append(nanct).transpose()
desc

# Look at means (and numbers of NaNs)
nanct = gwdby_wds.isna().sum()
nanct.name = 'NaNs'
desc = gwdby_wds.describe(percentiles=[]).append(nanct).transpose()
desc

# Look at means (and numbers of NaNs)
nanct = gwdby_wds2.isna().sum()
nanct.name = 'NaNs'
desc = gwdby_wds2.describe(percentiles=[]).append(nanct).transpose()
desc

# Export
gwdby_wds2.to_csv(out_creole_wds_GWDBYavg)
out_creole_wds_GWDBYavg_gn = os.path.join(home, 'data', 'creole_wooddensity_GWDBYavg_gn.csv')
gwdby_wds2['mean_gn'].to_csv(out_creole_wds_GWDBYavg_gn)

#%% Troubleshoot get_mean_wds_3 to better estimate standard deviations
def get_mean_wds_3(df1, df2, lookup_df, binom_fld='binomial', comm_name_fld='all_names'):
    # Combine WD statistics from GWD and BY to species, genus, and family level
    wd_sp = agg_wd_stats_2dfs(df1, df2, group_fld = binom_fld, suffix = '_sp')
    wd_gn = agg_wd_stats_2dfs(df1, df2, group_fld = 'genus', suffix = '_gn')
    wd_fm = agg_wd_stats_2dfs(df1, df2, group_fld = 'family', suffix = '_fm')
    # Aggregate by creole
    wd_sp_cr = lookup_df.join(wd_sp, on=binom_fld).groupby(comm_name_fld).median()
    wd_gn_cr = lookup_df.join(wd_gn, on='genus').groupby(comm_name_fld).median()
    wd_fm_cr = lookup_df.join(wd_fm, on='family').groupby(comm_name_fld).median()
    # Join
    creole_wds = wd_sp_cr.join(wd_gn_cr).join(wd_fm_cr)
    creole_wds = fillNAs_highertaxonlevel(creole_wds)
    print(f'NaNs in WDs aggregated by creole and joined: \n{creole_wds.isna().sum()}')
    return(creole_wds)

df1 = gwd_df
df2 = by_df
lookup_df = creole_lookup
# Combine WD statistics from GWD and BY to species, genus, and family level
wd_sp = agg_wd_stats_2dfs(df1, df2, group_fld = binom_fld, suffix = '_sp')
wd_gn = agg_wd_stats_2dfs(df1, df2, group_fld = 'genus', suffix = '_gn')
wd_fm = agg_wd_stats_2dfs(df1, df2, group_fld = 'family', suffix = '_fm')
# Aggregate by creole
wd_sp_cr = lookup_df.join(wd_sp, on=binom_fld).groupby(comm_name_fld).median()
wd_gn_cr = lookup_df.join(wd_gn, on='genus').groupby(comm_name_fld).median()
wd_fm_cr = lookup_df.join(wd_fm, on='family').groupby(comm_name_fld).median()
# Join
creole_wds = wd_sp_cr.join(wd_gn_cr).join(wd_fm_cr)
creole_wds = fillNAs_highertaxonlevel(creole_wds)

len(wd_gn)
wd_gn.isna().sum()
len(wd_gn_cr)
wd_gn_cr.isna().sum()
creole_wds.isna().sum()

# when I change sd_10 to std, I get many more standard deviation measurements
def agg_wd_stats_2dfs(df1, df2, group_fld = 'family', agg_fld = 'wd', suffix = '_fm'):
    wd_agg1 = df1.groupby(group_fld)[agg_fld].agg(mean='mean', sd='std', med='median')
    wd_agg2 = df2.groupby(group_fld)[agg_fld].agg(mean='mean', sd='std', med='median')
    # Aggregate stats for the two groups
    groups = pd.concat([wd_agg1, wd_agg2], sort=False).groupby(group_fld)
    # perform separate aggregation on each column
    means = groups['mean'].mean()
    sds = groups['sd'].agg(sd=sd_pooled)
    meds = groups['med'].median()
    # Join the three genus-level statistics
    wd_agg = sds.join(meds).join(means)
    # Rename
    wd_agg.rename(columns={'mean':'mean'+suffix, 'sd':'sd'+suffix, 'med':'med'+suffix}, inplace=True)
    return(wd_agg)

group_fld = 'genus'
wd_agg1 = df1.groupby(group_fld)['wd'].agg(mean='mean', sd=sd_10, med='median')
wd_agg2 = df2.groupby(group_fld)['wd'].agg(mean='mean', sd='std', med='median')
len(wd_agg1)
wd_agg1.isna().sum()
len(wd_agg2)
wd_agg2.isna().sum()

#%% Compare/QC
gwdby_wds2.sort_index(axis=1).loc[['pwa valye'], :]
gwdby_wds.sort_index(axis=1).loc[['pwa valye'],:]

diff = gwdby_wds2['mean_gn'] - gwdby_wds['mean_gn']
diff.describe()
diff[diff > 0.2]
by_df
gwd_wds.loc[['bwa mit'], :].sort_index(axis=1)
by_wds.loc[['bwa mit'], :].sort_index(axis=1)
gwdby_wds2.loc[['bwa mit'], :].sort_index(axis=1)
gwdby_wds.loc[['bwa mit'], :].sort_index(axis=1)
field_species[field_species[comm_name_fld]=='bwa mit'].genus
gwd_df[gwd_df.genus == 'Eugenia'].wd

#%% Import wood densities from R - these tend to use dataset means more often than makes sense
# Initialize with my lookup table.
wds_fromR = field_species[[comm_name_fld, 'family', 'genus', 'species']]

# GWD+BY species, genus, and family
wds = pd.read_csv(os.path.join(home, 'getWoodDensity_GWDBYspgnfm.csv'))
wds = wds.groupby(['family', 'genus', 'species'])['meanWD', 'sdWD'].first().rename(columns={'meanWD':'mean_GWDBYspgnfm', 'sdWD': 'sd_GWDBYspgnfm'})
wds_fromR = wds_fromR.join(wds, on=['family', 'genus', 'species'])
# GWD+BY Genus and Family
wds = pd.read_csv(os.path.join(home, 'getWoodDensity_GWDBYgnfm.csv'))
wds = wds.groupby(['family', 'genus'])['meanWD', 'sdWD'].first().rename(columns={'meanWD':'mean_GWDBYgnfm', 'sdWD': 'sd_GWDBYgnfm'})
wds_fromR = wds_fromR.join(wds, on=['family', 'genus'])
# GWD+BY Genus
wds = pd.read_csv(os.path.join(home, 'getWoodDensity_GWDBYgn.csv'))
wds = wds.groupby(['family', 'genus'])['meanWD', 'sdWD'].first().rename(columns={'meanWD':'mean_GWDBYgn', 'sdWD': 'sd_GWDBYgn'})
wds_fromR = wds_fromR.join(wds, on=['family', 'genus'])
# GWD species, genus, and family
wds = pd.read_csv(os.path.join(home, 'getWoodDensity_GWDspgnfm.csv'))
wds = wds.groupby(['family', 'genus', 'species'])['meanWD', 'sdWD'].first().rename(columns={'meanWD':'mean_GWDspgnfm', 'sdWD': 'sd_GWDspgnfm'})
wds_fromR = wds_fromR.join(wds, on=['family', 'genus', 'species'])
# GWD genus, and family
wds = pd.read_csv(os.path.join(home, 'getWoodDensity_GWDgnfm.csv'))
wds = wds.groupby(['family', 'genus'])['meanWD', 'sdWD'].first().rename(columns={'meanWD':'mean_GWDgnfm', 'sdWD': 'sd_GWDgnfm'})
wds_fromR = wds_fromR.join(wds, on=['family', 'genus'])
# GWD+BY Genus and Family with CentralAmericaTrop
wds = pd.read_csv(os.path.join(home, 'getWoodDensity_GWDBYgnfm_CAT.csv'))
wds = wds.groupby(['family', 'genus'])['meanWD', 'sdWD'].first().rename(columns={'meanWD':'mean_GWDBYgnfm_CAT', 'sdWD': 'sd_GWDBYgnfm_CAT'})
wds_fromR = wds_fromR.join(wds, on=['family', 'genus'])
# Get mean for each creole name
creole_wds_fromR = wds_fromR.groupby(comm_name_fld).median()

#QC
wds_fromR[wds_fromR['all_names']=='bayawonn']
wds_fromR[wds_fromR['genus']=='Ochroma'].sd_GWDBYgnfm.unique()
creole_wds_fromR[creole_wds_fromR.index=='zakasya'].mean_GWDBYgnfm
creole_wds_fromR[creole_wds_fromR.index=='bwa madam'].sd_GWDBYgnfm
gwdby_wds[gwdby_wds.index=='bayawonn']
gwd_wds[gwd_wds.index=='bayawonn']
by_wds[by_wds.index=='bayawonn']
def sd_10(x): return(x.std() if x.count() > 10 else np.nan)
wd_agg = gwdby_df.groupby('genus')['wd'].agg(mean_gn='mean', sd_gn=sd_10, med_gn='median', count='count')
wd_agg = gwd_df.groupby('genus')['wd'].agg(mean_gn='mean', sd_gn=sd_10, med_gn='median', count='count')
wd_agg = by_df.groupby('genus')['wd'].agg(mean_gn='mean', sd_gn=sd_10, med_gn='median', count='count')
wd_agg[wd_agg.index=='Prosopis']

wds_fromR[wds_fromR['genus']=='Acacia'].mean_GWDgnfm.unique()
wd_agg[wd_agg.index=='Acacia']

#%% Compare/QC
creole_wds_BY.loc[['pwa valye'], :]
creole_wds_pyGWD.loc[['pwa valye'],:]
creole_wds_GWDBYavg.loc[['pwa valye'], :]
creole_wds_fromR.loc[['pwa valye'],['mean_GWDspgnfm', 'mean_GWDgnfm']]
creole_wds_fromR.columns

diff = creole_wds_fromR['mean_GWDBYgnfm'] - gwdby_wds['mean_gn']
diff.describe()
diff[diff < -0.16]
creole_wds_fromR.loc[['bois cotelette'], ['mean_GWDBYgn']]
gwdby_wds.loc[['bois cotelette'], :]
wds_fromR[wds_fromR[comm_name_fld]=='bois cotelette'].genus
wds_fromR[wds_fromR.genus == 'Drypetes'].mean_GWDBYgnfm
wd_agg[wd_agg.index == 'Drypetes']
gwdby_wds[gwdby_wds.index=='bois cotelette']
#%%
creole_wds_concat = creole_wds_BY.join(gwdby_wds, rsuffix='_GWDBY')\
                .join(gwdby_wds2, lsuffix='_BY', rsuffix='_GWDBYavg')\
                .join(creole_wds_fromR)


#%% Look at relationships between columns
creole_wds_means = creole_wds_concat.filter(like='mean')
creole_wds_gnfm_means = creole_wds_concat[['mean_gn_BY', 'mean_gn_GWDBY', 'mean_gn_GWDBYavg', 'mean_GWDBYgn']]

# Correlations
creole_wds_means.corr()

# Scatter matrix
fig, ax = plt.subplots(figsize=(12,12))
pd.plotting.scatter_matrix(creole_wds_means, alpha=1, ax=ax)
fig.savefig(os.path.join(home, 'qc_plots', 'wds_scattmatrix_all.png'))

# Correlations visualized
fig = plt.figure(figsize=(9, 9))
plt.matshow(creole_wds_means.corr('pearson'), fignum=fig.number)
plt.xticks(range(creole_wds_means.shape[1]), creole_wds_means.columns, fontsize=14, rotation=45)
plt.yticks(range(creole_wds_means.shape[1]), creole_wds_means.columns, fontsize=14)
cb = plt.colorbar()
cb.ax.tick_params(labelsize=14);
fig.savefig(os.path.join(home, 'qc_plots', 'wds_corrmatrix_all.png'))

# Scatter matrix
fig, ax = plt.subplots(figsize=(12,12))
pd.plotting.scatter_matrix(creole_wds_gnfm_means, alpha=1, ax=ax)
fig.savefig(os.path.join(home, 'qc_plots', 'wds_scattmatrix_all_gnfm.png'))

# Correlations visualized
fig = plt.figure(figsize=(9, 9))
plt.matshow(creole_wds_gnfm_means.corr('spearman'), fignum=fig.number)
plt.xticks(range(creole_wds_gnfm_means.shape[1]), creole_wds_gnfm_means.columns, fontsize=14, rotation=45)
plt.yticks(range(creole_wds_gnfm_means.shape[1]), creole_wds_gnfm_means.columns, fontsize=14)
cb = plt.colorbar()
cb.ax.tick_params(labelsize=14);
fig.savefig(os.path.join(home, 'qc_plots', 'wds_corrmatrix_all_gnfm_spearman.png'))

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
creole_wds_means = creole_wds_concat.filter(like='mean')

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

# Filter WD means to only species in field data
# Load pre-processed field data
df = pd.read_csv(out_filled_data_fname)
sp_creole_in_data = df.sp_creole.unique()
data_wd_means = creole_wds_means[creole_wds_means.index.isin(sp_creole_in_data)]
# Look at Nulls
data_wd_means[data_wd_means['mean_sp_BY'].isna()].index
data_wd_means[data_wd_means['mean_gn_BY'].isna()].index
data_wd_means[data_wd_means['mean_gn_GWDBY'].isna()].index
data_wd_means[data_wd_means['mean_gn_GWDBYavg'].isna()].index
data_wd_means[data_wd_means['mean_sp_GWDBYavg'].isna()].index
data_wd_means[data_wd_means['mean_gnfm'].isna()].index

# dalmari is an example where the binomial has a BY wood density, but the mean_GWDBYgn does not.
'dalmari' in data_wd_means[data_wd_means['mean_gn_BY'].isna()].index
# mean_BYgn does have a value for dalmari, as it should
'abbe marron' in data_wd_means[data_wd_means['mean_gn_BY'].isna()].index
# I import bwayo_densities_2.csv as the supplemental wood density data in getWoodDensity.
# Maybe getWoodDensity doesn't attempt to aggregate by genus for the supplemental table.


by_wds = pd.read_csv(out_wd_fname, index_col=['genus', 'species'])
binom = field_species.loc[field_species.all_names == 'dalmari', ['genus', 'species']]
by_wds.loc[binom.values.tolist()[0]]
by_wds.loc[binom.values.tolist()[0][0]]
# bwayo_densities_2 has a WD value for the dalmari binomial. That genus only appears with that binomial.

binom = field_species.loc[field_species.all_names == 'trompette', ['genus', 'species']]
by_wds.loc[binom.values.tolist()[0]]
by_wds.loc[binom.values.tolist()[0][0]]
# Trompette is the same situation.

binom = field_species.loc[field_species.all_names == 'bande', ['genus', 'species']]
by_wds.loc[binom.values.tolist()[0]]
by_wds.loc[binom.values.tolist()[0][0]]
# The bande genus (no species taxa specified) has no WD value

by_wds.loc['Garcinia']
field_species.loc[field_species.genus == 'Garcinia']

#%% Make wood density array (join one to many) of same form as data df for use in computeAGB()
# Load pre-processed field data
df = pd.read_csv(out_filled_data_fname)
sp_creole_in_data = df.sp_creole.unique()
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
field_binoms = lookup_field_species[binom_fld]
gwd_lookup = gwd_df.loc[gwd_df['binomial'].isin(field_binoms)]

#%% Replicate getWoodDensity from BIOMASS - using World as region and including the Bwa Yo wood densities
# def getWoodDensity(gwd_df, genus, species, stand=None, family=None, xtra_wd_data=None, verbose=True):
xtra_wd_data = field_species
xtra_wd_data = xtra_wd_data[['creole', binom_fld, 'genus', 'species_abbr', 'family', 'wd_avg']].rename(columns={'wd_avg':'wd', binomial_fld: 'binomial', 'species_abbr':'species'})
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
