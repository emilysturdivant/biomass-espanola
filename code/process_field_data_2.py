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
    binomials_split = df[binomial_fld].str.split(' ',expand=True)
    df = df.assign(
        genus=binomials_split[0].str.capitalize(),
        species_extd=binomials_split.iloc[:, 1:].apply(lambda x: ' '.join(x.dropna().astype(str).values), axis=1),
        species=binomials_split[1].str.lower()
        )
    # Remove 'spp.' from species column
    df = df.assign(species=df['species'].replace('spp.', np.nan))
    return(df)

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

#%%
# Set working directory
home = r'/Users/emilysturdivant/GitHub/biomass-espanola'
# home = r'/home/esturdivant/code/biomass-espanola' # work desktop

#%% Filenames
by_fname = os.path.join(home, 'data', 'bwayo_species_2.xlsx')
json_fname = os.path.join(home, 'standardize_creole.json')
out_wd_fname = os.path.join(home, 'data', 'bwayo_densities_2.csv')
lookup_fname = os.path.join(home, 'data', 'exploded_specieslookup.csv')

#%% Create name-to-binomial-to-wd lookup table
# Load species table digitized from Bwa Yo and split binomial into genus and species
by_df = load_supplemental_species_data(by_fname, binomial_fld='species_binomial',
    in_wd_fld='by_spec_grav', save_wooddensities=out_wd_fname)
by_df.drop(columns=['by_spec_grav', 'by_names', 'dot_names', 'name_guesses', 'species_extd'])\
    .rename(columns={'wd_avg':'wd'})\
    .to_csv(os.path.join(home, 'data', 'bwayo_densities_wFam.csv'), index=False)

# Explode Bwa Yo DF by the creole names column. For every row with multiple creole names, duplicate species row.
creole_df = explode_names_to_specieslookup(by_df, collist=['by_names', 'dot_names', 'name_guesses'],
    outnames_col='all_names', remove_patterns=[r'p.p.$', r'\?'])
# Export to CSV
creole_df.to_csv(lookup_fname, index=False)
