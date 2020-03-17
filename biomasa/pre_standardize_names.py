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
import sys
try:
    proj_dir = os.path.dirname(os.path.realpath(__file__))
except:
    proj_dir = os.path.dirname(os.path.realpath('process_field_data_2.py'))
sys.path.append(proj_dir) # Add the script location to the system path just to make sure this works.
from biomasa.core import *

print(python_version())

#%% Set filenames
home = proj_dir #r'/Users/emilysturdivant/GitHub/biomass-espanola'
# home = r'/home/esturdivant/code/biomass-espanola' # work desktop
json_fname = os.path.join(home, 'data', 'standardize_creole.json')
by_fname = os.path.join(home, 'data', 'bwayo_species_2.xlsx')
lookup_fname = os.path.join(home, 'data', 'exploded_specieslookup.csv')

#%% Standardize species names
name_to_alts = {
    'bayawonn': ['buayawonn'],
    'bresillet': ['bouziyet', 'bwabusiet', 'brissiet', 'breziyet'],
    'bois blanc': ['biosblanc', 'boiblanc', 'boisblanc', 'bois blac', 'bos blanc', 'bios blanc'],
    'bois bleu': ['bwa bleu', 'bwableu', 'baubleu'],
    'bwa dom': ['buadom', 'boadon', 'bwa don'],
    "bois d'homme": ['boisdhomme'],
    'bwa doti': ['bwadolti'],
    'jambet': ['yambet'],
    'bwa kaka': ['bwacaca'],
    'kapab': ['kapalp'],
    'bois lait': ['boislet'],
    'bwa loray': ['lorai'],
    'bois major': ['mayor'],
    'bwa palmis': ['bwa palmia', 'bwapalmis'],
    'bwa petro': ['bopetro', 'boapreta', 'paupreta', 'parpreto'],
    'bwa pini':['bwapini', 'bwabwa pini', 'guapini', 'bwa bwa pini'],
    'poupe': ['pope'],
    'bwa santi': ['bwa senti'],
    'bwa savann': ['bwa saban'],
    'savann': ['saban', 'salbann'],
    'bois savane': ['bois savanne'],
    'bwa chenn': ['boachen'],
    'dalmari': ['delmari'],
    'delen': ['de lin', 'dele', 'delin'],
    'madlenn': ['madeleine'],
    'divi divi': ['divi'],
    'kaliptis': ['eucaliptos', 'eucalypto', 'eucalipto'],
    'figye': ['fieuier', 'ficus'],
    'flambwayan': ['framboyan'],
    'fwenn': ['fruen', 'fuen'],
    'frene': ['fren'],
    'gommier': ['gombier', 'gomier'],
    'gwayav': ['guayaba'],
    'gwenn': ['guen'],
    'kachiman': ['kashima', 'kashuma'],
    'cachiman': ['cachimem', 'cachemam', 'cachimam'],
    'calebasse': ['calabase', 'calbesse', 'calbasse'],
    'acajou': ['acayu', 'acayou'],
    'casse': ['cass'],
    'cassia': ['casseia'],
    'kaymit': ['cawimite', 'cayimit', 'kaymite', 'camite'],
    'kenep': ['queneb', 'kuinip', 'quinipi'],
    'kokoye': ['coco'],
    'koma':['bwa koma', 'lakoma'],
    'corossol': ['corosol', 'corosole', 'corolosol',  'collossol', 'corosore', 'corosorole'],
    'latanye': ['latanie', 'latani', 'latenier', 'letanier'],
    'lorie': ['lo rieue', 'loreiue'],
    'madam yas': ['madame jass'],
    'magerit': ['magritte'],
    'mango': ['mnago', 'mangos', 'mamgo'],
    'monben': ['mombe','momber','monbein','monbin', 'momben', 'monbenr'],
    'papay': ['papaya'],
    'pendoula': ['pandola'],
    'pistach': ['pistache'],
    'piyon': ['piyong'],
    'risin': ['arisin'],
    'satanye':['saiteyen', 'santayet', 'santeyet', 'santyet', 'satanyet', 'satayet'],
    'cedre': ['cede'],
    'sikren': ['sicre','sikre', 'sakrin', 'sacrin'],
    'citron': ['citroin'],
    'tamarenn': ['tamarindo', 'tamarin'],
    'tcha tcha': ['chacha', 'bwachacha'],
    'twompet': ['trompete'],
    'twa pawol': ['twa parol', 'twa palol'],
    'zaboka': ['aguacates', 'aguacate'],
    'acacia': ['acassia'],
    'zamann': ['zanmann'],
    'zoranj dous': ['naranja', 'naranaja'], # more general would be sitwon, zoranj
    'limon frans': ['limon'],
    np.nan: ['1', '2', '3', '4', 'aaa', 'bbb', 'crickcrack', 'kkk', 'kambala', 'tambala', 'wichinpit', 'z']
    }
# Create replacement dict to standardize inconsistencies
alt_to_name = dict((v,k) for k,vs in name_to_alts.items() for v in vs)
json.dump(alt_to_name, open(json_fname, 'w'))

# Save original as table
name_to_alts_df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in name_to_alts.items() ])).transpose()
name_to_alts_df.to_csv(os.path.join(home, 'data', 'standardize_creole_table.csv'))

#%% Create name-to-binomial-to-wd lookup table
# Load table digitized from Bwa Yo and split binomial to genus // species
by_df = load_supplemental_species_data(by_fname,
    binomial_fld='species_binomial',
    in_wd_fld='by_spec_grav',
    save_wooddensities=out_wd_fname)
by_df\
    .drop(columns=['by_spec_grav', 'by_names', 'dot_names', 'name_guesses', 'species_extd'])\
    .rename(columns={'wd_avg':'wd'})\
    .to_csv(os.path.join(home, 'data', 'bwayo_densities_wFam.csv'), index=False)

# Explode Bwa Yo DF by the creole names column.
# For every row with multiple creole names, duplicate species row.
lookup_df = explode_names_to_specieslookup(by_df,
    collist=['by_names', 'dot_names', 'name_guesses'],
    outnames_col='all_names', remove_patterns=[r'p.p.$', r'\?'])
# Export to CSV
lookup_df.to_csv(lookup_fname, index=False)


# #%% Extract values in field data from BY df
# field_species_uniq = pd.Series(spec_ser.unique())
# field_species = by_df.loc[by_df['creole'].isin(field_species_uniq)].reset_index(drop=True)
#
# #%% Create creole to species_name lookup and species_name to bwayo_wd lookup
# field_species.head(2)
#
# # What to do about creole names that match more than one genus?
# field_species[['genus', 'creole', 'wd_avg']].groupby('creole').agg(lambda x: x.unique())
#
# def most_common(List):
#     if not len(List):
#         return(np.nan)
#     elif len(List) == 1:
#         return(List[0])
#     else: # list is greater than 1
#         return(mode(List))
# from statistics import mode
# l = ['Mor','Ulm', 'Set', 'Ulm', 'Set', 'Set']
# max(set(l), key = l.count)
# mode(l)
# most_common(l)
#
# lookup_genus_by_creole = field_species[['genus', 'creole']].groupby('creole').first()
#
#
# # What to do about creole names that match more than one family?
# field_species[['family', 'creole', 'wd_avg']].groupby('creole').agg(lambda x: x.unique())
#
# # What to do about genus spp. entries? What does the splitGenusSpecies() R function do?
#
#
#
# #%% Write to excel
# # out_fname = os.path.join(home, 'wooddensity_lookup_20191218.xlsx')
# # field_species.to_excel(out_fname, index=False)
#
# '''---------------------------------------------------------------------------
# Create lookup table to find wood density for each creole name.
# Uses average family values from Global Wood Density database. (1/1/2020)
# ---------------------------------------------------------------------------'''
# #%% Parse Global Wood Density
# gwd_fname = os.path.join(home, 'data', 'GlobalWoodDensityDatabase.xlsx')
# gwd_df = pd.read_excel(gwd_fname, sheet_name='Data', header=0,
#     names=['gwd_num', 'gwd_family', 'gwd_binomial', 'gwd_density',
#         'gwd_region', 'gwd_ref_no'],
#     index_col='gwd_num',
#     converters={'gwd_binomial':lambda x : x.lower(),
#         'gwd_family':lambda x : x.lower()})
#
# # Extract rows from GWD with species present in field data
# field_binoms = field_species['by_binomial']
# gwd_lookup = gwd_df.loc[gwd_df['gwd_binomial'].isin(field_binoms)]
#
# # Print QC info
# print(f'Unique species in field data: {len(field_binoms)}')
# gwd_unique_matches = gwd_lookup['gwd_binomial'].unique()
# print(f'Matching species in GWD: {len(gwd_unique_matches)}')
# unlisted = field_binoms[~field_binoms.isin(gwd_df['gwd_binomial'])]
# print(f'Field species missing from GWD: {len(unlisted)}')
# unlisted_names = field_species_uniq[~field_species_uniq.isin(by_creole_names)]
# print(f'Unidentified species in field data: {len(unlisted_names)}')
#
# #%% Create lookup with family averaged values
# fam_mean_density = gwd_df.groupby('gwd_family')['gwd_density'].mean()
# fam_means = field_species.join(fam_mean_density, on='family', how='left', lsuffix='_gwd')
# # Write to lookup_famAvgWD_allRegions.xlsx: mean WD for each family (all regions)
# # fieldspec_byfam.to_excel(os.path.join(home, 'lookup_famAvgWD_allRegions.xlsx'), index=False)
#
# # For each creole name, get mean WD for all matching species
# creole_dens = fam_means.groupby('creole')['gwd_density'].mean()
# # Write lookup_famAvgWD_byCreole.xlsx: mean family WD (all regions) of all families matching the creole name (usually only one match)
# creole_dens.to_excel(os.path.join(home, 'lookup_famAvgWD_byCreole.xlsx'), index=True)
