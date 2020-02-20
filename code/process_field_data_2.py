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

def split_species_binomial(df, binomial_fld='species_binomial'):
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
    # Remove 'spp.' from species_abbr column
    df = df.assign(species_abbr=df['species_abbr'].replace('spp.', np.nan))
    return(df)

def get_plot_data(data_fname, col_ints, verbose=True):
    # Get data for plot indicated by col_ints
    df = pd.read_excel(data_fname, 'Plots', skiprows=[0,1,2], header=0,
        usecols=col_ints, names=['sp_creole','dbh_cm','ht_m','ba_m2'], #dtype={'ht_m':'float'},
        skip_blank_lines = True, na_values='0',
        converters={'sp_creole':lambda x : strip_accents(x.strip().lower()), 'ht_m': lambda x: str(x).strip("''")})
    # Drop rows that have no data and convert height column to float.
    df = df.dropna(axis=0, how='all').astype({'ht_m':'float'})
    # Create padding row in plots without trees
    if not len(df):
        df = df.reindex(labels=[1], fill_value=np.nan)
    # Get plot number, shapefile name, and area for given plot
    plot = pd.read_excel(data_fname, 'Plots', header=0, usecols=col_ints, nrows=1)
    df['plot_no'] = plot.columns[0].split('#')[1]
    if verbose:
        print(f'Plot number: {df.plot_no} | {plot.iloc[0,0]} | Area: {plot.iloc[0,2]} ha')
    df = df.assign(plot_shp=plot.iloc[0,0], plot_area=plot.iloc[0,2])
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

#%%
# Set working directory
home = r'/Users/emilysturdivant/GitHub/biomass-espanola'
# home = r'/home/esturdivant/code/biomass-espanola' # work desktop

#%% Filenames
data_fname = os.path.join(home, 'data', 'haiti_biomass_v2.xlsx')
by_fname = os.path.join(home, 'data', 'bwayo_species_2.xlsx')
json_fname = os.path.join(home, 'standardize_creole.json')
out_wd_fname = os.path.join(home, 'data', 'bwayo_densities_2.csv')

#%% Load species table digitized from Bwa Yo and split binomial into genus and species
by_df = pd.read_excel(by_fname, header=0, usecols=[0,1,2,3,4,5],
    converters={'species_binomial':lambda x : x.lower(),
        'family':lambda x : x.split(' (')[0].lower().capitalize()})

# Split species binomial into genus, species, and species_abbr
by_df = split_species_binomial(by_df, binomial_fld='species_binomial')

# Create supplemental wood density means from all Bwa Yo values. Convert WD range to mean
bwayo_wd = by_df.loc[~by_df['by_spec_grav'].isna(),
    ['genus', 'species_abbr', 'by_spec_grav']]
bwayo_wd = bwayo_wd.assign(wd_avg=[np.mean([float(i) for i in
    re.findall(r"[0-9.]+", str(s))]) for s in bwayo_wd['by_spec_grav']])

# Export to CSV
bwayo_wd_for_export = bwayo_wd[['genus', 'species_abbr', 'wd_avg']].rename(columns={'species_abbr':'species', 'wd_avg':'wd'})
bwayo_wd_for_export.to_csv(out_wd_fname), index=False)

#%% Explode Bwa Yo DF by the creole names column. For every row with multiple creole names, duplicate species row.
# Join wood density averages to the Bwa Yo DF.
df = by_df.join(bwayo_wd['wd_avg'])
# Explode DF to name to species lookup.
creole_df = explode_names_to_specieslookup(df, collist=['by_names', 'dot_names', 'name_guesses'],
    outnames_col='all_names', remove_patterns=[r'p.p.$', r'\?'])
# Export to CSV
creole_df.to_csv(os.path.join(home, 'data', 'exploded_specieslookup.csv'), index=False)

#%% Extract species in field data from BY df - prep for lookup table
# Create series of all species columns (labeled 'sp')
spec_df = pd.read_excel(data_fname, 'Plots', skiprows=[0,1,2],
                        usecols=lambda x : x.startswith('sp'))
alt_to_name = json.load(open(json_fname)) # dict to standardize creole names
spec_ser = spec_df.stack().apply(lambda x : strip_accents(x).strip().lower()).replace(alt_to_name)

# Extract species in field data from exploded BY df
field_species_uniq = pd.Series(spec_ser.unique())
field_species = creole_df.loc[creole_df['all_names'].isin(field_species_uniq)].reset_index(drop=True)

# Export CSV
field_species.to_csv(os.path.join(home, 'data', 'master_lookup2.csv'), index=False)

#%% Convert unknown species to something standardized
unknowns = field_species_uniq[~field_species_uniq.isin(field_species['all_names'])]
ct = 0
for val in unknowns:
    ct += 1
    alt_to_name[val] = np.nan
    for key, value in alt_to_name.items():
         if val == value:
             alt_to_name[key] = np.nan
print(f'{ct} species in field data not identified.')

#%% Load field data and convert to table of all plots stacked
# Iteratively load each plot in turn and concatenate
col_init = [0,1,2,3]
df = get_plot_data(data_fname, col_init, verbose=False)
for i in range(1,36):
    col_ints = np.add(col_init, 4*i)
    df2 = get_plot_data(data_fname, col_ints, verbose=False)
    df = pd.concat([df, df2], ignore_index=True)

df
# Standardize creole names
df.loc[:, 'sp_creole'] = df.sp_creole.replace(alt_to_name)

# QC
df.loc[60:70, :]
df[df['dbh_cm'] == 0]
df[df['dbh_cm'].isna()]

# Replace np.nan with 0 in dbh_cm column
# df.replace({np.nan:0}, inplace=True)
df.to_csv(os.path.join(home, 'data', 'haiti_biomass_v2_stacked.csv'), index=False)
df = pd.read_csv(os.path.join(home, 'data', 'haiti_biomass_v2_stacked.csv'))

df.join(creole_df, )

#%% Make and export plots DF
mplots = df[['plot_no', 'plot_shp', 'plot_area']].groupby('plot_no').first()
mplots.to_csv(os.path.join(home, 'data', 'haiti_biomass_v2_mplots.csv'), index=True)
mplots


#%% Perform QC on data
#%% Identify duplicate records
# Look at consecutive duplicated entries
cols = ["sp_creole", "dbh_cm", 'ht_m', "plot_no"]
dups1 = (df[cols].shift() == df[cols]).all(axis=1).replace({False: np.nan})
de_dup = df[cols].loc[dups1.fillna(dups1.shift(-1)).fillna(False)]
de_dup
len(de_dup)

# All duplicates
df[df.duplicated(subset=['sp_creole', 'dbh_cm', 'ht_m', 'plot_no'], keep=False)]
# Duplicates within a plot (plot 9)
df_plt9 = df[df['plot_no'] == '9']
df_plt9[df_plt9.duplicated(subset=['sp_creole', 'dbh_cm', 'ht_m'], keep=False)]
df_plt9_mango = df_plt9[df_plt9['sp_creole'] == 'mango']
df_plt9_mango[df_plt9_mango.duplicated(subset=['dbh_cm', 'ht_m'], keep=False)]

#%% Check DBH measurements - look for outliers

fig, ax = plt.subplots(figsize=(16,8))
ax.scatter(df['dbh_cm'], df['plot_no'])
ax.set_xlabel('DBH')
ax.set_ylabel('Plot Number')
plt.show()
# plt.savefig('example.pdf')
#%%
plot = df[~df['dbh_cm'].isna()].boxplot(column='dbh_cm', figsize=(16, 2), whis=[1,99], notch=True, vert=False, showmeans=True, meanline=True)
fig = plot.get_figure()
fig.savefig(os.path.join(home, 'qc_plots', 'boxes_dbh_plotno.png'))

#%%
plot = df[~df['dbh_cm'].isna()].boxplot(column='dbh_cm', by='plot_no', figsize=(16,8), whis=[1,99], notch=True, showmeans=True)
fig = plot.get_figure()
fig.savefig(os.path.join(home, 'qc_plots', 'boxes_dbh_plotno.png'))

#%%
plot = df.plot.scatter(x='dbh_cm', y='ht_m', figsize=(16,8), alpha=0.05)
fig = plot.get_figure()
fig.tight_layout()
fig.savefig(os.path.join(home, 'qc_plots', 'scatter_dbh_height_suspicious_pattern.png'))

df['ht_m'].plot.hist(bins=40)
df[['dbh_cm', 'ht_m']].hist(bins=20)

#%%
plot = df[~df['ht_m'].isna()].boxplot(column='ht_m', by='plot_no', figsize=(16,8))
fig = plot.get_figure()
fig.savefig(os.path.join(home, 'qc_plots', 'boxes_height_plotno.png'))

df_lame = df[df.sp_creole=='figye']
df_lame[df_lame.ht_m > 18]


#%% Look at DBH and height distributions by species
df['sp_creole'] = df['sp_creole'].fillna('UNKNOWN')
plot = df.boxplot(column='dbh_cm', by='sp_creole', figsize=(8,16), vert=False)
fig = plot.get_figure()
fig.savefig(os.path.join(home, 'qc_plots', 'boxes_dbh_species.png'))

#%%
df_sub = df[df['sp_creole'].isin(['mango', 'latanye lame', 'figye', 'monben', 'kaliptis'])]
plot = df_sub.boxplot(column='dbh_cm', by='sp_creole', figsize=(8,4), vert=False)
fig = plot.get_figure()
fig.tight_layout()
fig.savefig(os.path.join(home, 'qc_plots', 'boxes_dbh_species_outliers.png'))

#%%
ax = df[df.sp_creole == 'mango'].plot.scatter(x='dbh_cm', y='ht_m', color='Blue', label='mango');
df[df.sp_creole == 'latanye lame'].plot.scatter(x='dbh_cm', y='ht_m', color='Green', label='latanye lame', ax=ax);
df[df.sp_creole == 'figye'].plot.scatter(x='dbh_cm', y='ht_m', color='Cyan', label='figye', ax=ax);
df[df.sp_creole == 'kampech'].plot.scatter(x='dbh_cm', y='ht_m', color='Yellow', label='kampech', ax=ax);
plot = df[df.sp_creole == 'kaliptis'].plot.scatter(x='dbh_cm', y='ht_m', color='Red', label='kaliptis', ax=ax);
fig = plot.get_figure()
fig.savefig(os.path.join(home, 'qc_plots', 'scatter_dbh_height_outliers.png'))

#%%
plot = df.boxplot(column='ht_m', by='sp_creole', figsize=(8,16), vert=False)
fig = plot.get_figure()
fig.savefig(os.path.join(home, 'qc_plots', 'boxes_ht_species.png'))
#%%

#%% Join genus to field data DF
# Simplistic: take the first genus matching the creole name
lookup_genus = field_species[['genus', 'family', 'species', 'creole']].groupby('creole').first()
dfjoined = df.join(lookup_genus, on='sp_creole', how='left')

plot = dfjoined.boxplot(column='dbh_cm', by='family', figsize=(8,16), vert=False)
fig = plot.get_figure()
fig.tight_layout()
fig.savefig(os.path.join(home, 'qc_plots', 'boxes_dbh_family.png'))
