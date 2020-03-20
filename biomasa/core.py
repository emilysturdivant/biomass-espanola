# -*- coding: utf-8 -*-
"""
by: Emily Sturdivant, emilysturdivant@gmail.com
requires: python 3.6+

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


#%% General use functions
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

def split_species_binomial(df, binomial_fld='species_binomial'):
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

#%% Functions to get field data from individual spreadsheets
def get_field_data_formulario(data_fname):
    df = pd.read_excel(data_fname, 'Hoja1', skiprows=[0,1,2,3,4,5,6,7,8], header=0,
        usecols=[1,2,3], names=['sp_creole','dbh_cm','ht_m'],
        skip_blank_lines = True, na_values='0',
        converters={'sp_creole':lambda x : strip_accents(x.strip().lower()), 'ht_m': lambda x: str(x).strip("''")})
    # Drop rows that have all Null values - usually parsed because of some merged cells in the excel sheet.
    df = df.dropna(subset=['sp_creole', 'dbh_cm', 'ht_m'], how='all')
    # get header information
    # Get plot number, shapefile name, and area for given plot
    head_df = pd.read_excel(data_fname, 'Hoja1', header=None, skiprows=[0],
            usecols=[0,1,2,3,4], nrows=3, skip_blank_lines = True)
    df = df.assign(
        plot_fname = os.path.basename(data_fname),
        date = head_df.loc[0,1],
        plot_id = head_df.loc[0,4],
        plot_loc = head_df.loc[1,2],
        data_scribe = head_df.loc[2,2]
        )
    return(df)

#%% Functions for processing field data
def load_supplemental_species_data(by_fname, binomial_fld='species_binomial',
    in_wd_fld='by_spec_grav', save_wooddensities='out_wd_fname.csv'):
    # Load species table as DF including genus, species_ext, and species means
    df = pd.read_excel(by_fname, header=0, usecols=[0,1,2,3,4,5],
        converters={binomial_fld:lambda x : x.lower(),
            'family':lambda x : x.split(' (')[0].lower().capitalize()})
    # Split species binomial into genus, species_ext, and species
    df = split_species_binomial(df, binomial_fld=binomial_fld)
    # Create supplemental wood density means. Convert WD range to mean
    wd_means = df.loc[~df[in_wd_fld].isna(), ['genus', 'species', in_wd_fld]]
    wd_means = wd_means.assign(wd_avg=[np.mean([float(i) for i in
        re.findall(r"[0-9.]+", str(s))]) for s in wd_means[in_wd_fld]])
    if save_wooddensities:
        # Export to CSV
        wd_means.drop(columns=in_wd_fld)\
            .rename(columns={'wd_avg':'wd'})\
            .to_csv(save_wooddensities, index=False)
    # Join wood density averages to the DF.
    df = df.join(wd_means['wd_avg'])
    return(df)

def explode_names_to_specieslookup(df, collist=['by_names', 'dot_names', 'name_guesses'],
    outnames_col='all_names', remove_patterns=[r'p.p.$', r'\?']):
    # Create dataframe that connects a local name to each species binomial that it might reference.
    # Concatenate name columns to all_names
    for col in collist:
        try:
            df[col] = df[col].str.split(',')
        except AttributeError:
            print(f'Column "{col}" may already be in list format.')
            pass
    df[outnames_col] = df[collist].stack().groupby(level=0).sum()
    # Explode
    creole_df = df.explode(outnames_col).reset_index(drop=True)
    # Tidy up - clean strings, drop columns, and drop duplicate rows
    creole_df[outnames_col] = creole_df[outnames_col].str.lower().str.strip()
    if outnames_col in collist:
        collist.remove(outnames_col)
    creole_df = creole_df.drop(collist, axis=1).drop_duplicates()
    # Remove optional characters
    for pat in remove_patterns:
        creole_df[outnames_col] = creole_df[outnames_col].replace(pat, '', regex=True)
    # Return
    return(creole_df)

#%% Functions for aggregating wood densities
def sd_10(x): return(x.std() if x.count() > 10 else np.nan)

def sd_pooled(x):
    if not any(x.isna()):
        return(np.sqrt(x.sum()/x.count()))
    else:
        return(np.nan)

def agg_wd_stats_2dfs(df1, df2, group_fld = 'family', agg_fld = 'wd', suffix = '_fm'):
    wd_agg1 = df1.groupby(group_fld)['wd'].agg(mean='mean', sd=sd_10, med='median')
    wd_agg2 = df2.groupby(group_fld)['wd'].agg(mean='mean', sd=sd_10, med='median')
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

def fillNAs_highertaxonlevel(creole_wds):
    # Genus filled from Family
    na_idx = creole_wds[creole_wds['mean_gn'].isna()].index
    # Median
    creole_wds['med_gn'] = creole_wds['med_gn'].fillna(creole_wds['med_fm'])
    # Mean
    creole_wds['mean_gn'] = creole_wds['mean_gn'].fillna(creole_wds['mean_fm'])
    # SDs. Only fill NAs in sd_10 where mean was filled.
    sdgn = creole_wds['sd_gn'].fillna(creole_wds['sd_fm'])
    creole_wds.loc[na_idx, 'sd_gn'] = sdgn[na_idx]

    # Species filled from genus
    na_idx = creole_wds[creole_wds['mean_sp'].isna()].index
    # Median
    creole_wds['med_sp'] = creole_wds['med_sp'].fillna(creole_wds['med_gn'])
    # Fill NAs in genus means
    creole_wds['mean_sp'] = creole_wds['mean_sp'].fillna(creole_wds['mean_gn'])
    # Fill NAs in genus SDs. Only fill NAs in sd_10 where mean was filled.
    sdsp = creole_wds['sd_sp'].fillna(creole_wds['sd_gn'])
    creole_wds.loc[na_idx, 'sd_sp'] = sdsp[na_idx]
    # Return
    return(creole_wds)

def get_mean_WDs(df, binom_fld='binomial', wd_fld='wd'):
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
    # Fill NAs
    creole_wds = fillNAs_highertaxonlevel(creole_wds)
    # print(f'NaNs in WDs after filling: \n{creole_wds.isna().sum()}')
    return(creole_wds)

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
