# -*- coding: utf-8 -*-
"""
by: Emily Sturdivant, emilysturdivant@gmail.com
requires: python 3.6+

OVERVIEW: Getting started with processing Haiti field data

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

print(python_version())

#%% Function copied from process_field_data_2
def split_species_binomial(df, binomial_fld='binomial'):
    def get_second_word(x):
        ser = []
        for lst in x:
            if len(lst) > 1:
                ser += [lst[1].strip()]
            else:
                ser += ['']
        return(ser)
    binomials_split = df[binomial_fld].str.split(' ',expand=True)
    df = df.assign(
        genus=binomials_split[0].str.capitalize(),
        species_extd=binomials_split.iloc[:, 1:].apply(lambda x: ' '.join(x.dropna().astype(str).values), axis=1),
        species=binomials_split[1].str.lower()
        )
    # Remove 'spp.' from species column
    df = df.assign(species=df['species'].replace('spp.', np.nan))
    return(df)

def get_mean_WDs(df, binom_fld='binomial', wd_fld='wd'):
    def sd_10(x): return(x.std() if x.count() > 10 else np.nan)
    # Species means
    wd_agg = df.groupby(binom_fld)[wd_fld].agg(mean='mean', sd=sd_10)
    df = df.join(wd_agg, on=binom_fld, rsuffix='_sp')
    # Genus means
    wd_agg = df.groupby('genus')[wd_fld].agg(mean='mean', sd='std', med_gn='median')
    df = df.join(wd_agg, on='genus', rsuffix='_gn')
    # Family means
    wd_agg = df.groupby('family')[wd_fld].agg(mean='mean', sd=sd_10, med_fm='median')
    df = df.join(wd_agg, on='family', rsuffix='_fm')
    # Fill NAs in mean_BYgn with mean_BYfm
    df['mean_gnfm'] = df['mean_gn'].fillna(df['mean_fm'])
    df['med_gnfm'] = df['med_gn'].fillna(df['med_fm'])
    df['sd_gnfm'] = df['sd_gn'].fillna(df['sd_fm'])
    # Drop wd and rename mean and sd
    df = df.drop(wd_fld, axis=1).rename(columns={'mean':'mean_sp', 'sd':'sd_sp'})
    return(df)

def get_mean_WDs_2(df, lookup_df, binom_fld='binomial', wd_fld='wd', comm_name_fld='all_names'):
    def sd_10(x): return(x.std() if x.count() > 10 else np.nan)
    # Species means
    wd_agg = df.groupby(binom_fld)[wd_fld].agg(mean_sp='mean', sd_sp=sd_10, med_sp='median')
    wd_cr_sp = lookup_df.join(wd_agg, on=binom_fld).groupby(comm_name_fld).median()
    # print(f'NaNs in WDs aggregated by species: \n{wd_agg.isna().sum()}')
    # print(f'NaNs in WDs aggregated by species and creole: \n{wd_cr_sp.isna().sum()}')
    # Genus means
    wd_agg = df.groupby('genus')[wd_fld].agg(mean_gn='mean', sd_gn=sd_10, med_gn='median')
    wd_cr_gn = lookup_df.join(wd_agg, on='genus').groupby(comm_name_fld).median()
    # print(f'NaNs in WDs aggregated by genus: \n{wd_agg.isna().sum()}')
    # print(f'NaNs in WDs aggregated by genus and creole: \n{wd_cr_gn.isna().sum()}')
    # Family means
    wd_agg = df.groupby('family')[wd_fld].agg(mean_fm='mean', sd_fm=sd_10, med_fm='median')
    wd_cr_fm = lookup_df.join(wd_agg, on='family').groupby(comm_name_fld).median()
    # print(f'NaNs in WDs aggregated by family: \n{wd_agg.isna().sum()}')
    # print(f'NaNs in WDs aggregated by family and creole: \n{wd_cr_fm.isna().sum()}')
    # Join
    creole_wds = wd_cr_sp.join(wd_cr_gn).join(wd_cr_fm)
    # print(f'NaNs in WDs aggregated by creole and joined: \n{creole_wds.isna().sum()}')
    # Fill NAs in genus
    creole_wds['med_gn'] = creole_wds['med_gn'].fillna(creole_wds['med_fm'])
    na_idx = creole_wds[creole_wds['mean_gn'].isna()].index
    # Fill NAs in genus means
    mngn = creole_wds['mean_gn'].fillna(creole_wds['mean_fm'])
    creole_wds.loc[na_idx, 'mean_gn'] = mngn[na_idx]
    # Fill NAs in genus SDs. Only fill NAs in sd_10 where mean was filled.
    sdgn = creole_wds['sd_gn'].fillna(creole_wds['sd_fm'])
    creole_wds.loc[na_idx, 'sd_gn'] = sdgn[na_idx]

    # Fill NAs in species
    creole_wds['med_sp'] = creole_wds['med_sp'].fillna(creole_wds['med_gn'])
    na_idx = creole_wds[creole_wds['mean_sp'].isna()].index
    # Fill NAs in genus means
    mnsp = creole_wds['mean_sp'].fillna(creole_wds['mean_gn'])
    creole_wds.loc[na_idx, 'mean_sp'] = mnsp[na_idx]
    # Fill NAs in genus SDs. Only fill NAs in sd_10 where mean was filled.
    sdsp = creole_wds['sd_sp'].fillna(creole_wds['sd_gn'])
    creole_wds.loc[na_idx, 'sd_sp'] = sdsp[na_idx]
    # print(f'NaNs in WDs after filling: \n{creole_wds.isna().sum()}')
    # return
    return(creole_wds)

#%%
# Set working directory
home = r'/Users/emilysturdivant/GitHub/biomass-espanola'
home = r'/home/esturdivant/code/biomass-espanola' # work desktop

#%% Filenames
by_wd_fname = os.path.join(home, 'data', 'bwayo_densities_wFam.csv')
plots_overview_fname = os.path.join(home, 'data', 'haiti_plots_meta.csv')
by_table_fname = os.path.join(home, 'data', 'exploded_specieslookup.csv')
out_filled_data_fname = os.path.join(home, 'data', 'haiti_data_filled.csv')
master_lookup = os.path.join(home, 'data', 'master_lookup_2.csv')

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

# Print QC info
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

# Combined
gwdby_wds = get_mean_WDs_2(df=gwdby_df, lookup_df=creole_lookup, binom_fld='binomial', wd_fld='wd', comm_name_fld='all_names')

by_wds = get_mean_WDs_2(df=by_df, lookup_df=creole_lookup, binom_fld='binomial', wd_fld='wd', comm_name_fld='all_names')

gwd_wds = get_mean_WDs_2(df=gwd_df, lookup_df=creole_lookup, binom_fld='binomial', wd_fld='wd', comm_name_fld='all_names')

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

# gwdby_wds1 = get_mean_WDs(gwdby_df, 'binomial')
# # Get mean WDs for each creole name
# creole_wds_pyGWDBY = creole_lookup\
#                     .join(wds.set_index('binomial'), on='binomial')\
#                     .groupby(comm_name_fld).mean()
# # QC
# creole_wds_pyGWDBY.describe()

#%% Get mean Bwa Yo wood density at the three levels.
# by_df = field_species[[comm_name_fld, 'family', 'genus', 'binomial', 'wd_avg']]
by_wds = get_mean_WDs_2(by_df, binom_fld='binomial', wd_fld='wd')
# Get mean WDs for each creole name
creole_wds_BY = by_wds.groupby(comm_name_fld).mean()
# QC
creole_wds_BY.loc[['pwa valye'],:]
creole_wds_BY.filter(like='mean_gn').describe()
creole_wds_BY.filter(like='mean_gn').isna().sum()


by_wds2 = get_mean_WDs_2(df=by_df, lookup_df=creole_lookup, binom_fld='binomial', wd_fld='wd_avg', comm_name_fld=comm_name_fld)
by_wds2.loc[['pwa valye'],:]
by_wds2.filter(like='mean_gn').describe()
by_wds2.filter(like='mean_gn').isna().sum()


#%% Create function to get wood densities and join after each get
creole_lookup = creole_lookup.rename(columns={'species_binomial':'binomial'})
gwd_wds = get_mean_WDs_2(df=gwd_df, lookup_df=creole_lookup, binom_fld='binomial', wd_fld='wd', comm_name_fld='all_names')

# Get mean Bwa Yo wood density at the three levels.
gwd_wds1 = get_mean_WDs(gwd_df, binom_fld='binomial', wd_fld='wd')
# gwd_wds1.join(creole_lookup.set_index('binomial').drop(['family', 'genus'], axis=1), on='binomial')
creole_wds_GWD = creole_lookup.drop(['family', 'genus'], axis=1).join(gwd_wds1.set_index('binomial'), on='binomial').drop(['gwd_ref_no'], axis=1)
# Get mean WDs for each creole name
creole_wds_GWD = creole_wds_GWD.groupby(comm_name_fld).mean()

# Compare two versions
gwd_wds.filter(like='mean').join(creole_wds_GWD.filter(like='mean'), lsuffix='2', rsuffix='1')
# Look at means (and numbers of NaNs)
nanct = gwd_wds.isna().sum()
nanct.name = 'NaNs'
desc = gwd_wds.describe(percentiles=[]).append(nanct).transpose()
desc

# Look at means (and numbers of NaNs)
nanct = creole_wds_GWD.isna().sum()
nanct.name = 'NaNs'
desc = creole_wds_GWD.describe(percentiles=[]).append(nanct).transpose()
desc


#%% Append BwaYo to GWD and get mean WDs
# Get BwaYo densities
by_df = pd.read_csv(by_wd_fname).rename(columns={'species_binomial': 'binomial'})

# Concatenate and calculate means
gwdby_df = pd.concat([gwd_df, by_df], sort=False)\
            .drop(['region', 'gwd_ref_no', 'species_extd'], axis=1)
gwdby_wds = get_mean_WDs(gwdby_df, 'binomial')

# Get mean WDs for each creole name
creole_wds_pyGWDBY = creole_lookup\
                    .join(wds.set_index('binomial'), on='binomial')\
                    .groupby(comm_name_fld).mean()

# QC
creole_wds_pyGWDBY.describe()

#%% Just get genus-family means
by_wds2 = get_mean_WDs_2(gwd_df, lookup_df=creole_lookup, binom_fld='binomial', wd_fld='wd', comm_name_fld=comm_name_fld)

# GWD
df = gwd_df
wd_fld = 'wd'
def sd_10(x): return(x.std() if x.count() > 10 else np.nan)
# Get genus mean WDs for each creole name
wd_agg = df.groupby('genus')[wd_fld].agg(mean='mean', sd_10=sd_10, med='median')
creole_wds = creole_lookup.join(wd_agg, on='genus').groupby(comm_name_fld).median()
# Get family mean WDs for each creole name
wd_agg = df.groupby('family')[wd_fld].agg(mean='mean', sd_10=sd_10, med='median')
creole_wds_fm = creole_lookup.join(wd_agg, on='family').groupby(comm_name_fld).median()
# Fill NAs. Only fill NAs in sd_10 where mean was filled.
na_idx = creole_wds[creole_wds_gn['mean'].isna()].index
filled = creole_wds.fillna(creole_wds_fm)
creole_wds.loc[na_idx, 'mean'] = filled.loc[na_idx, 'mean']
creole_wds.loc[na_idx, 'sd_10'] = filled.loc[na_idx, 'sd_10']

creole_wds_gnfm2 = creole_wds


pd.concat([creole_wds_pyGWD['mean_gnfm'], creole_wds_gnfm2['mean']], axis=1)

# Look at means (and numbers of NaNs)
nanct = creole_gnfm_GWDBYavg1.isna().sum()
nanct.name = 'NaNs'
desc = creole_gnfm_GWDBYavg1.describe(percentiles=[]).append(nanct).transpose()
desc
#%% Take average of BY and GWD species means and then calculate creole means
gwd_gnfm = gwd_wds.set_index('binomial')['mean_gnfm'].to_frame(name = 'gwd')
by_gnfm = by_wds.set_index('binomial')['mean_gnfm'].to_frame(name='by')
gwdby_gnfm_avg = gwd_gnfm.join(by_gnfm).mean(axis=1).to_frame(name='gwdby_gnfm_avg')
creole_gnfm_GWDBYavg1 = creole_lookup\
                    .join(gwdby_gnfm_avg, on='binomial', how='outer')\
                    .groupby(comm_name_fld).mean()

# Look at means (and numbers of NaNs)
nanct = creole_gnfm_GWDBYavg1.isna().sum()
nanct.name = 'NaNs'
desc = creole_gnfm_GWDBYavg1.describe(percentiles=[]).append(nanct).transpose()
desc

df_concat = pd.concat([gwd_wds, by_wds], sort=False)
df_concat.head()
gwdby_avg = df_concat.groupby('binomial').mean()
gwdby_avg = df_concat.groupby('genus').mean()
# Get mean WDs for each creole name
creole_wds_GWDBYavg1 = field_species[[comm_name_fld, 'genus']]\
                    .join(gwdby_avg, on='genus', how='outer')\
                    .groupby(comm_name_fld).mean()

creole_wds_GWDBYavg1
# Look at means (and numbers of NaNs)
nanct = creole_wds_GWDBYavg1.isna().sum()
nanct.name = 'NaNs'
desc = creole_wds_GWDBYavg1.describe(percentiles=[]).append(nanct).transpose()
desc

#%% Take average of BY and GWD creole means

creole_wds_pyGWD.head()
creole_wds_BY.head()
creole_wds_GWDBYavg = pd.concat([creole_wds_pyGWD, creole_wds_BY])\
                        .groupby(df_concat.index).mean()

# Look at means (and numbers of NaNs)
nanct = creole_wds_GWDBYavg.isna().sum()
nanct.name = 'NaNs'
desc = creole_wds_GWDBYavg.describe(percentiles=[]).append(nanct).transpose()
desc
nans_idx = creole_wds_GWDBYavg['mean_gnfm'].isna()
nans_idx.sum()
creole_wds_pyGWD[nans_idx].isna().sum()
creole_wds_BY[nans_idx]
creole_wds_pyGWDBY[nans_idx].isna().sum()
creole_wds_GWDBYavg[creole_wds_GWDBYavg['mean_gnfm'].isna()]


#%% Import wood densities from R
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
creole_wds_fromR = wds_fromR.groupby(comm_name_fld).mean()

#QC
wds_fromR[wds_fromR['genus']=='Ochroma'].mean_GWDgnfm

#%% Compare/QC
creole_wds_BY.loc[['pwa valye'], :]
creole_wds_pyGWD.loc[['pwa valye'],:]
creole_wds_pyGWDBY.loc[['pwa valye'],:]
creole_wds_GWDBYavg.loc[['pwa valye'], :]
creole_wds_fromR.loc[['pwa valye'],['mean_GWDspgnfm', 'mean_GWDgnfm']]
creole_wds_fromR.columns

diff = creole_wds_fromR['mean_GWDgnfm'] - creole_wds_pyGWD['mean_gnfm']
diff.describe()
diff[diff > 0.4]
creole_wds_fromR.loc[['bois madame'], ['mean_GWDgnfm']]
creole_wds_pyGWD.loc[['bois madame'], :]
by_wds[by_wds[comm_name_fld]=='bois madame'].binomial
gwd_wds[gwd_wds.genus == 'Ochroma']

#%%
creole_wds_concat = creole_wds_BY.join(creole_wds_pyGWDBY, rsuffix='_GWDBY')\
                .join(creole_wds_GWDBYavg, lsuffix='_BY', rsuffix='_GWDBYavg')\
                .join(creole_wds_fromR)


#%% Look at relationships between columns
creole_wds_means = creole_wds_concat.filter(like='mean')
creole_wds_gnfm_means = creole_wds_concat[['mean_gnfm_BY', 'mean_gnfm_GWDBY', 'mean_gnfm_GWDBYavg', 'mean_GWDBYgnfm']]

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
data_wd_means = creole_wds_means[creole_wds_means.index.isin(sp_creole_in_data)]
# Look at Nulls
data_wd_means[data_wd_means['mean_BYsp'].isna()].index
data_wd_means[data_wd_means['mean_BYgn'].isna()].index
data_wd_means[data_wd_means['mean_GWDBYspgnfm'].isna()].index
data_wd_means[data_wd_means['mean_GWDBYgnfm'].isna()].index
data_wd_means[data_wd_means['mean_GWDBYgn'].isna()].index
data_wd_means[data_wd_means['mean_GWDspgnfm'].isna()].index
data_wd_means[data_wd_means['mean_GWDBYgn'].isna()].index
# dalmari is an example where the binomial has a BY wood density, but the mean_GWDBYgn does not.
'dalmari' in data_wd_means[data_wd_means['mean_BYgn'].isna()].index
# mean_BYgn does have a value for dalmari, as it should
'abbe marron' in data_wd_means[data_wd_means['mean_BYgn'].isna()].index
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
