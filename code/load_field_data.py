# -*- coding: utf-8 -*-
"""
by: Emily Sturdivant, emilysturdivant@gmail.com
requires: python 3.6+

Load field inventory data from individual plot tables

Python environment kernel: using 'py3_geo' on Ubuntu and 'Python 3' on Mac
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

#%%
# Set working directory
home = r'/Users/emilysturdivant/GitHub/biomass-espanola'
home = r'/home/esturdivant/code/biomass-espanola' # work desktop

#%% Filenames
data_folder = os.path.join(home, 'data', 'formularios')
json_fname = os.path.join(home, 'standardize_creole.json')
by_table_fname = os.path.join(home, 'data', 'exploded_specieslookup.csv')
plots_overview_fname = os.path.join(home, 'data', 'haiti_plots_meta.csv')
out_raw_data_fname = os.path.join(home, 'data', 'haiti_data_raw.csv')
out_cleaned_data_fname = os.path.join(home, 'data', 'haiti_data.csv')
out_filled_data_fname = os.path.join(home, 'data', 'haiti_data_filled.csv')
master_lookup = os.path.join(home, 'data', 'master_lookup_2.csv')

#%% Load datasets created in other scripts
# Load creole_df
creole_df = pd.read_csv(by_table_fname)
# Load alt_to_name dict - dict to fix typos in common names
alt_to_name = json.load(open(json_fname))

#%%
# Read data from all formularios into dataframe
df = None
for f in os.listdir(data_folder):
    data_fname = os.path.join(data_folder, f)
    if df is not None:
        df2 = get_field_data_formulario(data_fname)
        df = pd.concat([df, df2], ignore_index=True)
    else:
        df = get_field_data_formulario(data_fname)

# Join NoBiomass plots and plot areas
mplots = pd.read_csv(plots_overview_fname)
df = df.join(mplots.set_index('plot_id'), on='plot_id', how='outer')

# Save
df.to_csv(out_raw_data_fname, index=False)
df = pd.read_csv(out_raw_data_fname)

#%% Prep to convert unknown species names to NaN
# List accounted-for names
names1 = pd.Series(v for v in alt_to_name.values())
unknowns = names1[~names1.isin(creole_df['all_names'])].unique()
if len(unknowns) > 0:
    print(f'Standardized name(s) "{unknowns}" not identified. Double check entries in name_to_alts dictionary.')
namelist = list(alt_to_name.keys()) + list(creole_df['all_names'])

# Get series of all names entered in field data
spec_ser = pd.Series(df[['sp_creole']].groupby('sp_creole').first().index)

# Convert unknown species to NaN
unknowns = spec_ser[~spec_ser.isin(namelist)]
for val in unknowns:
    alt_to_name[val] = np.nan

#%% Perform replacement
df['sp_creole'] = df['sp_creole'].replace(alt_to_name)
field_sp_names = pd.Series(df[['sp_creole']].groupby('sp_creole').first().index)

#%% Save data
df.to_csv(out_cleaned_data_fname, index=False)
df = pd.read_csv(out_cleaned_data_fname)

# Quick QC
len(df)
df.sample(5)
len(df[df['dbh_cm'] == 0])
len(df[df['ht_m'] == 0])
len(df[df['dbh_cm'].isna()])
df[df['dbh_cm'].isna()]
nans = df[pd.concat([df['dbh_cm'].isna(), df['sp_creole'].isna(), df['ht_m'].isna()], axis=1).all(axis=1)]
if len(nans) > 8:
    print(f"WARNING: There are {len(nans)} rows with Null in all three sp, dbh, and ht, more than ther are NoBiomass plots.")

#%% Fill Null DBHs with plot median
# Check boxplot by plot with means to compare means and medians
plot = df[~df['dbh_cm'].isna()].boxplot(column='dbh_cm', by='plot_no', figsize=(16,8), whis=[1,99], notch=True, showmeans=True)
plot = df[~df['ht_m'].isna()].boxplot(column='ht_m', by='plot_no', figsize=(16,8), whis=[1,99], notch=True, showmeans=True)

df_filled = df.copy()
df_filled['dbh_cm'] = df.groupby('plot_id')['dbh_cm']\
    .transform(lambda x: x.fillna(x.median()))
df_filled[df['dbh_cm'].isna()].sample(5)

# Check boxplot by plot with means to compare means and medians
plot = df_filled[~df_filled['dbh_cm'].isna()]\
    .boxplot(column='dbh_cm', by='plot_no', figsize=(16,8),
        whis=[1,99], notch=True, showmeans=True)
sub_df = df[~df['dbh_cm'].isna()]
sub_df = sub_df[sub_df['sp_creole'].isin(['bois blanc', 'gommier', 'satanye'])]
sub_df = sub_df[df['plot_no']==33]
plot = sub_df.boxplot(column='dbh_cm', by='sp_creole', figsize=(16,8),
        whis=[1,99], notch=True, showmeans=True)

sub_df = df_filled[df_filled['plot_no']==33]
sub_df = sub_df[sub_df['sp_creole'].isin(['bois blanc', 'gommier', 'satanye'])]
plot = sub_df.boxplot(column='dbh_cm', by='sp_creole', figsize=(16,8),
        whis=[1,99], notch=True, showmeans=True)

#%% Save data
df_filled.to_csv(out_filled_data_fname, index=False)
df_filled = pd.read_csv(out_filled_data_fname)

#%% Extract species in field data from BY df - prep for lookup table

# Extract species in field data from exploded BY df
field_species_lookup = creole_df.loc[creole_df['all_names'].isin(field_sp_names)].reset_index(drop=True)

# Export CSV
field_species_lookup.to_csv(master_lookup, index=False)




#%% Perform QC on data
df = df_filled
#%% Identify duplicate records
# Look at consecutive duplicated entries
cols = ["sp_creole", "dbh_cm", 'ht_m', "plot_no"]
dups1 = (df[cols].shift() == df[cols]).all(axis=1).replace({False: np.nan})
de_dup = df[cols].loc[dups1.fillna(dups1.shift(-1)).fillna(False)]
de_dup
len(de_dup)

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

#%% Boxplot by plot with means
plot = df[~df['dbh_cm'].isna()].boxplot(column='dbh_cm', by='plot_no', figsize=(16,8), whis=[1,99], notch=True, showmeans=True)
fig = plot.get_figure()
fig.savefig(os.path.join(home, 'qc_plots', 'boxes_dbh_plotno.png'))

#%%
plot = df.plot.scatter(x='dbh_cm', y='ht_m', figsize=(16,8), alpha=0.05)
fig = plot.get_figure()
fig.tight_layout()
fig.savefig(os.path.join(home, 'qc_plots', 'scatter_dbh_height_suspicious_pattern.png'))

#%%
df['ht_m'].plot.hist(bins=40)
df[['dbh_cm', 'ht_m']].hist(bins=20)

#%%
plot = df[~df['ht_m'].isna()].boxplot(column='ht_m', by='plot_no', figsize=(16,8))
fig = plot.get_figure()
fig.savefig(os.path.join(home, 'qc_plots', 'boxes_height_plotno.png'))

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

#%% Scatterplot of DBH v. height for certain colors
ax = df.plot.scatter(x='dbh_cm', y='ht_m', figsize=(16,8), alpha=0.05);
df[df.sp_creole == 'mango'].plot.scatter(x='dbh_cm', y='ht_m', color='Blue', label='mango', ax=ax);
df[df.sp_creole == 'latanye'].plot.scatter(x='dbh_cm', y='ht_m', color='Maroon', label='latanye', ax=ax);
df[df.sp_creole == 'figye'].plot.scatter(x='dbh_cm', y='ht_m', color='Cyan', label='figye', ax=ax);
df[df.sp_creole == 'campeche'].plot.scatter(x='dbh_cm', y='ht_m', color='Yellow', label='campeche', ax=ax);
df[df.sp_creole == 'UNKNOWN'].plot.scatter(x='dbh_cm', y='ht_m', color='Green', label='UNKNOWN', ax=ax);
plot = df[df.sp_creole == 'kaliptis'].plot.scatter(x='dbh_cm', y='ht_m', color='Red', label='kaliptis', ax=ax);

fig = plot.get_figure()
fig.savefig(os.path.join(home, 'qc_plots', 'scatter_dbh_height_outliers.png'))

#%% Boxplots of height distribution by species name
plot = df.boxplot(column='ht_m', by='sp_creole', figsize=(8,16), vert=False)
fig = plot.get_figure()
fig.savefig(os.path.join(home, 'qc_plots', 'boxes_ht_species.png'))
#%%

#%% Join genus to field data DF
# Simplistic: take the first genus matching the creole name
lookup_genus = field_species_lookup[['genus', 'family', 'species', 'all_names']].groupby('all_names').first()
dfjoined = df.join(lookup_genus, on='sp_creole', how='left')

plot = dfjoined.boxplot(column='dbh_cm', by='family', figsize=(8,16), vert=False)
fig = plot.get_figure()
fig.tight_layout()
fig.savefig(os.path.join(home, 'qc_plots', 'boxes_dbh_family.png'))
