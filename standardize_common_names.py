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
from platform import python_version

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

#%%
# Set working directory
home = r'/Users/emilysturdivant/GitHub/biomass-espanola'

#%% Work with field data - Look at species in all plots
data_fname = os.path.join(home, 'haiti_biomass_v2.xlsx')

# Standardize species names
name_to_alts = {
    'abey': ['abi'],
    'bayawonn': ['buayawonn', 'bayawonn'],
    'benzoliv': ['moringa', 'olivye'],
    'breziyet': ['bouziyet', 'bwabusiet', 'brissiet'],
    'bwa bande': ['bande'],
    'bwa blan': ['biosblanc', 'boiblanc', 'boisblanc', 'bois blac', 'bos blanc', 'bois blanc'],
    'bois bleu': ['bwa bleu', 'bwableu', 'bois bleu'],
    'bwa don': ['buadom', 'boisdhomme'],
    'bwa doti': ['bwadolti'],
    'bwa doule': ['noni'],
    'bwa jambet': ['jambet', 'mayi bouyi'],
    'bwa kaka': ['bwacaca', 'bois caca'],
    'bwa kapab': ['capable', 'kapab'],
    'bwa let': ['boislet', 'bois lait', 'letchie'],
    'bwa loray': ['bois loraille', 'bois loraille', 'lorai'],
    'bwa mabel': ['mabe'], # or mabi?
    'bwa majo': ['mayor', 'bois major'],
    'bwa nwa': ['bois noir', 'nwa', 'nua'],
    'bwa palmis': ['bwa palmia', 'bwapalmis'],
    'bwa pen': ['pin'],
    'bwa petro': ['bwapetro'], # unlocated
    'bwa pine':['bwapini', 'pini', 'bwabwa pini', 'bwa pini', 'bwa pine'],
    'bwa poupe': ['pope', 'poupe'],
    'bwa santi': ['bwa senti'],
    'bwa savann': ['bwa saban', 'bois savanne', 'saban', 'salbann', 'savann'],
    'chadek': ['chadeque'],
    'chenn': ['bwa chenn', 'chene', 'chenn'],
    'damari': ['dalmari', 'delmari'],
    'delen': ['de lin', 'dele', 'delin', 'madeleine', 'madlenn'],
    'divi divi': ['divi'],
    'kaliptis': ['eucaliptos', 'eucalypto', 'eucalyptus', 'eucalipto'],
    'figye': ['fieuier', 'ficus', 'figuier'],
    'flambwayan': ['framboyan'],
    'fwenn': ['fren', 'frene'],
    'gayak': ['gaiac'],
    'gommier': ['gombier', 'gomier'],
    'gomye': ['bursera'],
    'gwayav': ['guayaba'],
    'gwenn': ['guen'],
    'kachiman': ['cachimem', 'cachemam', 'kashima', 'cachimam', 'cachiman'],
    'kalbas': ['calabase', 'calbesse', 'calbasse', 'calebassier', 'calebasse'],
    'kampech': ['campeche'],
    'kajou': ['acayu', 'kajou', 'acajou', 'acayou'],
    'kas': ['cass', 'casse'],
    'kasya': ['cassia', 'casseia'],
    'kaymit': ['cawimite', 'cayimit', 'kaymit', 'kaymite', 'camite'],
    'kenep': ['quenepe', 'queneb'],
    'kokoye': ['coco'],
    'koma':['bwa koma', 'lakoma', 'akoma'],
    'kowosol': ['corosol', 'corosole', 'corolosol', 'corossol', 'collossol', 'corosore', 'corosorole', 'guanabana', 'corossolier'],
    'latanye lame': ['lame'],
    'latanye': ['latanye', 'latanie', 'latani', 'latenier', 'letanier', 'latanier'],
    'lisina': ['leucaena'],
    'lorie': ['lo rieue'],
    'madam yas': ['madame jass', 'madame yass'],
    'magerit': ['magritte'],
    'mango': ['mnago', 'mangos', 'mamgo'],
    'monben': ['mombe','momber','monbein','monbin', 'monben', 'momben', 'monbenr', 'mombin'],
    'nim': ['lila', 'neem/lila', 'neem'],
    'papay': ['papaya'],
    'pendoula': ['pandola'],
    'pistach': ['pistache'],
    'piyon': ['gliricidia', 'piyong'],
    'pwa valye': ['pois vallier'],
    'risin': ['arisin'],
    'satanye':['saiteyen', 'santayet', 'santeyet', 'santyet', 'satanyet', 'satayet'],
    'sed': ['cede'],
    'sikren': ['sicre','sikre', 'sakrin', 'sacrin'],
    'sitwon': ['citroin', 'limon'],
    'tamarenn': ['tamarindo', 'tamarenn', 'tamarin'],
    'tcha tcha': ['chacha', 'bwachacha'],
    'twompet': ['trompete', 'trompette'],
    'twa pawol': ['twa parol', 'twa palol','twa pabel','twa pable'],
    'zaboka': ['aguacates', 'aguacate', 'sambuca'],
    'zabriko': ['abricot'],
    'zakasya': ['acassia', 'acacia'],
    'zamann': ['amande', 'zanmann'],
    'zoranj dous': ['naranaja', 'naranja'], # more general would be sitwon, zoranj
    }

# name_to_alts = {}
# for k,v in alt_to_name.items():
#     if v in name_to_alts:
#         name_to_alts[v].append(k)
#     else:
#         name_to_alts[v] = [k]
# Create replacement dict for to standardize inconsistencies
alt_to_name = {}
for k, vals in name_to_alts.items():
    for v in vals:
        alt_to_name[v] = k

#%% Create series of all species columns (labeled 'sp')
spec_df = pd.read_excel(data_fname, 'Plots', skiprows=[0,1,2],
                        usecols=lambda x : x.startswith('sp'))
spec_ser = spec_df.stack().apply(lambda x : strip_accents(x).strip().lower()).replace(alt_to_name)

#%% List unique species with number of occurences and write to excel sheet.
spec_list_cnts = spec_ser.value_counts().rename_axis(['common_name'])
field_spec_list = pd.Series(spec_ser.unique())

# Look at different groupings
specs_mult = spec_list_cnts[spec_list_cnts > 2]
specs_mult.sort_index()
err_specs = spec_list_cnts[spec_list_cnts < 3]
err_specs.sort_index()

# Run replace again
# spec_ser.replace(alt_to_name, inplace=True)
# spec_list_cnts = spec_ser.value_counts().rename_axis(['common_name'])

#%% Write to excel
out_lookup = os.path.join(home, 'unique_species_standardized_20191218.xlsx')
spec_list_cnts.to_excel(out_lookup, 'round3', index=True)


#%% Look at species table digitized from Bwa Yo
by_fname = os.path.join(home, 'bwayo_species.xlsx')
by_df = pd.read_excel(by_fname, header=0, usecols=[0,1,2,3], names=['BY_binomial', 'creole', 'BY_spec_grav', 'family'], converters={'BY_binomial':lambda x : x.lower(), 'family':lambda x : x.split(' (')[0].lower()})

#%% Explode DF by the creole names column. Convert values to list and for every row with multiple names, create multiple species columns.
by_df = (
        by_df.assign(creole=by_df.creole.str.split(','))
                .explode('creole')
                .reset_index(drop=True))
by_df = by_df.assign(creole=by_df.creole.str.strip())
by_spec_list = by_df.creole.unique().tolist()

# Compare standardized common species names to species names in field data
unlisted_names = field_spec_list[~field_spec_list.isin(by_spec_list)]

#%% Extract values in field data from BY df, only need creole, BY_binomial, family, BY_spec_grav.
field_spec_list
fieldspec_lookup = by_df.loc[by_df['creole'].isin(field_spec_list)]

#%% Write to excel
out_fname = os.path.join(home, 'wooddensity_lookup_20191218.xlsx')
fieldspec_lookup.to_excel(out_fname, index=False)






#%% Load species table (Especies.xlsx by José Luis)
species_fname = os.path.join(home, 'Especies.xlsx')
spec_lookup = pd.read_excel(species_fname, 'Hoja1',
                            names=['common_name', 'common_alt', 'species', 'synonyms'])
common_syns = spec_lookup.loc[:,['common_name', 'common_alt']]
common_syns = common_syns.stack().apply(lambda x : x.strip().lower()).unstack()
comm_names = common_syns.stack().unique().tolist()


# Compare standardized common species names to species names in field data
unlisted_names = []
for sp in spec_list_cnts.index:
    if not sp in comm_names:
        unlisted_names += [sp]

unlisted_names
