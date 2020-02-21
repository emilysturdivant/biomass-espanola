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

def split_species_binomial(df, binomial_fld='by_binomial'):
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
    return(df)

#%%
# Set working directory
home = r'/Users/emilysturdivant/GitHub/biomass-espanola'
home = r'/home/esturdivant/code/biomass-espanola' # work desktop

#%% Work with field data - Look at species in all plots
data_fname = os.path.join(home, 'data', 'haiti_biomass_v2.xlsx')
gwd_fname = os.path.join(home, 'data', 'GlobalWoodDensityDatabase.xlsx')
by_fname = os.path.join(home, 'data', 'bwayo_species.xlsx')

# Standardize species names
# name_to_alts = {
#     'abbe': ['abi', 'abey'],
#     'bayawonn': ['buayawonn', 'bayawonn'],
#     'benzoliv': ['moringa', 'olivye'],
#     'breziyet': ['bouziyet', 'bwabusiet', 'brissiet'],
#     'bwa bande': ['bande'],
#     'bwa blan': ['biosblanc', 'boiblanc', 'boisblanc', 'bois blac', 'bos blanc', 'bois blanc'],
#     'bois bleu': ['bwa bleu', 'bwableu', 'bois bleu', 'baubleu'],
#     'bwa dom': ['buadom', 'boisdhomme', 'boadon', 'bwa don'],
#     'bwa doti': ['bwadolti'],
#     'bwa doule': ['noni'],
#     'bwa jambet': ['jambet', 'yambet'],
#     'bwa kaka': ['bwacaca', 'bois caca'],
#     'bwa kapab': ['capable', 'kapab', 'kapalp'],
#     'bwa let': ['boislet', 'bois lait'],
#     'bwa loray': ['bois loraille', 'bois loraille', 'lorai'],
#     'bwa mabel': ['mabe'], # or mabi? # or momben?
#     'bwa majo': ['mayor', 'bois major'],
#     'bwa nwa': ['bois noir', 'nwa', 'nua'],
#     'bwa palmis': ['bwa palmia', 'bwapalmis'],
#     'bwa pen': ['pin'],
#     'bwa pit': ['bopetro', 'boapreta', 'paupreta', 'parpreto', 'bwa petro']
#     'bwa pine':['bwapini', 'pini', 'bwabwa pini', 'bwa pini', 'bwa pine', 'guapini'],
#     'bwa poupe': ['pope', 'poupe'],
#     'bwa santi': ['bwa senti', 'bwa satsi'], # not sure...
#     'bwa savann': ['bwa saban', 'bois savanne', 'saban', 'salbann', 'savann'],
#     'chenn': ['bwa chenn', 'chene', 'boachen'],
#     'damari': ['dalmari', 'delmari'],
#     'delen': ['de lin', 'dele', 'delin', 'madeleine', 'madlenn'],
#     'divi divi': ['divi'],
#     'kaliptis': ['eucaliptos', 'eucalypto', 'eucalyptus', 'eucalipto'],
#     'figye': ['fieuier', 'ficus', 'figuier'],
#     'flambwayan': ['framboyan'],
#     'fwenn': ['fren', 'frene', 'fruen', 'fuen'],
#     'gayak': ['gaiac'],
#     'gommier': ['gombier', 'gomier'],
#     'gomye': ['bursera'],
#     'gwayav': ['guayaba'],
#     'gwenn': ['guen'],
#     'kachiman': ['cachimem', 'cachemam', 'kashima', 'cachimam', 'cachiman', 'kashuma'],
#     'kalbas': ['calabase', 'calbesse', 'calbasse', 'calebassier', 'calebasse'],
#     'kampech': ['campeche'],
#     'kajou': ['acayu', 'kajou', 'acajou', 'acayou'],
#     'kas': ['cass', 'casse'],
#     'kasya': ['cassia', 'casseia'],
#     'kaymit': ['cawimite', 'cayimit', 'kaymit', 'kaymite', 'camite'],
#     'kenep': ['quenepe', 'queneb', 'kuinip', 'quinipi'],
#     'kokoye': ['coco'],
#     'koma':['bwa koma', 'lakoma', 'akoma'],
#     'kowosol': ['corosol', 'corosole', 'corolosol', 'corossol', 'collossol', 'corosore', 'corosorole', 'guanabana', 'corossolier'],
#     'latanye lame': ['lame'],
#     'latanye': ['latanye', 'latanie', 'latani', 'latenier', 'letanier', 'latanier'],
#     'diversifolia': ['leucaena', 'diversifolia'],
#     'lorie': ['lo rieue', 'loreiue'],
#     'madam yas': ['madame jass', 'madame yass'],
#     'magerit': ['magritte'],
#     'mango': ['mnago', 'mangos', 'mamgo'],
#     'monben': ['mombe','momber','monbein','monbin', 'monben', 'momben', 'monbenr', 'mombin'],
#     'nim': ['lila', 'neem/lila', 'neem'],
#     'papay': ['papaya'],
#     'pendoula': ['pandola'],
#     'pistach': ['pistache'],
#     'piyon': ['gliricidia', 'piyong'],
#     'pwa valye': ['pois vallier'],
#     'risin': ['arisin'],
#     'satanye':['saiteyen', 'santayet', 'santeyet', 'santyet', 'satanyet', 'satayet'],
#     'sed': ['cede'],
#     'sikren': ['sicre','sikre', 'sakrin', 'sacrin'],
#     'sitwon': ['citroin', 'limon'],
#     'tamarenn': ['tamarindo', 'tamarenn', 'tamarin'],
#     'tcha tcha': ['chacha', 'bwachacha'],
#     'twompet': ['trompete', 'trompette'],
#     'twa pawol': ['twa parol', 'twa palol','twa pabel','twa pable'],
#     'zaboka': ['aguacates', 'aguacate', 'sambuca'],
#     'zabriko': ['abricot'],
#     'zakasya': ['acassia', 'acacia'],
#     'zamann': ['amande', 'zanmann'],
#     'zoranj dous': ['naranaja', 'naranja'], # more general would be sitwon, zoranj
#     }
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
json.dump(alt_to_name, open(os.path.join(home, 'standardize_creole.json'), 'w'))

# Save original as table
name_to_alts_df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in name_to_alts.items() ])).transpose()
name_to_alts_df.to_csv(os.path.join(home, 'data', 'standardize_creole_table.csv'))


#%% Extract values in field data from BY df
field_species_uniq = pd.Series(spec_ser.unique())
field_species = by_df.loc[by_df['creole'].isin(field_species_uniq)].reset_index(drop=True)

#%% Create creole to species_name lookup and species_name to bwayo_wd lookup
field_species.head(2)

# What to do about creole names that match more than one genus?
field_species[['genus', 'creole', 'wd_avg']].groupby('creole').agg(lambda x: x.unique())

def most_common(List):
    if not len(List):
        return(np.nan)
    elif len(List) == 1:
        return(List[0])
    else: # list is greater than 1
        return(mode(List))
from statistics import mode
l = ['Mor','Ulm', 'Set', 'Ulm', 'Set', 'Set']
max(set(l), key = l.count)
mode(l)
most_common(l)

lookup_genus_by_creole = field_species[['genus', 'creole']].groupby('creole').first()


# What to do about creole names that match more than one family?
field_species[['family', 'creole', 'wd_avg']].groupby('creole').agg(lambda x: x.unique())

# What to do about genus spp. entries? What does the splitGenusSpecies() R function do?



#%% Write to excel
# out_fname = os.path.join(home, 'wooddensity_lookup_20191218.xlsx')
# field_species.to_excel(out_fname, index=False)

'''---------------------------------------------------------------------------
Create lookup table to find wood density for each creole name.
Uses average family values from Global Wood Density database. (1/1/2020)
---------------------------------------------------------------------------'''
#%% Parse Global Wood Density
gwd_df = pd.read_excel(gwd_fname, sheet_name='Data', header=0,
    names=['gwd_num', 'gwd_family', 'gwd_binomial', 'gwd_density',
        'gwd_region', 'gwd_ref_no'],
    index_col='gwd_num',
    converters={'gwd_binomial':lambda x : x.lower(),
        'gwd_family':lambda x : x.lower()})

# Extract rows from GWD with species present in field data
field_binoms = field_species['by_binomial']
gwd_lookup = gwd_df.loc[gwd_df['gwd_binomial'].isin(field_binoms)]

# Print QC info
print(f'Unique species in field data: {len(field_binoms)}')
gwd_unique_matches = gwd_lookup['gwd_binomial'].unique()
print(f'Matching species in GWD: {len(gwd_unique_matches)}')
unlisted = field_binoms[~field_binoms.isin(gwd_df['gwd_binomial'])]
print(f'Field species missing from GWD: {len(unlisted)}')
unlisted_names = field_species_uniq[~field_species_uniq.isin(by_creole_names)]
print(f'Unidentified species in field data: {len(unlisted_names)}')

#%% Create lookup with family averaged values
fam_mean_density = gwd_df.groupby('gwd_family')['gwd_density'].mean()
fam_means = field_species.join(fam_mean_density, on='family', how='left', lsuffix='_gwd')
# Write to lookup_famAvgWD_allRegions.xlsx: mean WD for each family (all regions)
# fieldspec_byfam.to_excel(os.path.join(home, 'lookup_famAvgWD_allRegions.xlsx'), index=False)

# For each creole name, get mean WD for all matching species
creole_dens = fam_means.groupby('creole')['gwd_density'].mean()
# Write lookup_famAvgWD_byCreole.xlsx: mean family WD (all regions) of all families matching the creole name (usually only one match)
creole_dens.to_excel(os.path.join(home, 'lookup_famAvgWD_byCreole.xlsx'), index=True)
