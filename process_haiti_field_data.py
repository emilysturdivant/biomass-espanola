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

#%% Read in first 4 columns
col_ints = [0,1,2,3]
col_ints = np.array([32, 33, 34, 35]) # first plot with trees
col_ints = np.add(col_ints, 4)

# Get plot number, shapefile name, and area for given plot
plot = pd.read_excel(data_fname, 'Plots', header=0, usecols=col_ints, nrows=1)
plot_no = plot.columns[0].split('#')[1]
plot_shp = plot.iloc[0,0]
plot_area = plot.iloc[0,2]
print(f'Plot number: {plot_no} | {plot_shp} | Area: {plot_area} ha')

# Get plot data
df = pd.read_excel(data_fname, 'Plots', skiprows=[0,1,2], header=0, usecols=col_ints, names=['sp','dap_cm','ht_m','ba_m2'], converters={'sp':lambda x : strip_accents(x.strip().lower())})
df.head()








# Standardize species names
name_to_alts = {
    'abey': ['abi'],
    'acacia': ['acassia'],
    'acajou': ['acayu', 'kajou'],
    'amande': ['almond', 'amand', 'amod'],
    'aguacate': ['aguacates'],
    'bande': ['bendi', 'bombe'],
    'bayawonn': ['buayawonn'],
    'bois blanc': ['biosblanc', 'boiblanc', 'boisblanc', 'boisblanc2'],
    'bois bleu': ['bwa bleu', 'bwableu'],
    'bois caca': ['bwacaca'],
    'bois lait': ['boislet'],
    'bois loraille': [],
    'bois major': ['mayor'],
    'bwa don': ['buadom', 'boisdhomme'],
    'bwa doti': ['bwadolti'],
    'bwa nwa': ['nua', 'nwa'],
    'bwa palmis': ['bwa palmia', 'bwapalmis', 'palmis'],
    'bwa petro': ['bwapetro'],
    'bwa pini':['bwapini', 'pini'],
    'bwa poupe': ['pope', 'poupe'],
    'cachiman': ['cachimem', 'cachemam', 'kashima'],
    'casse': ['cass'],
    'cayoux': ['cayaoux'],
    'calebasse': ['calabase', 'calbesse', 'calbasse', 'calebassier'],
    'damari': ['dalmari', 'delmari'],
    'delin': ['de lin', 'dele', 'delen'],
    'eucalipto': ['eucaliptos', 'eucalypto', 'eucalyptus'],
    'figuier': ['fieuier'],
    'fwenn': ['fren'],
    'gommier': ['gombier', 'gomier'],
    'kapab': ['capable'],
    'kasya': ['cassia'],
    'kenep': ['quenepe', 'queneb'],
    'kowosol': ['corosol', 'corosole', 'corolosol', 'corossol', 'collossol', 'corosore', 'corosorole'],
    'lame': ['lamp'],
    'latanier': ['latanye', 'latanie', 'latani'],
    'lorie': ['lo rieue', 'lorai'],
    'madame yass': ['madame jass'],
    'madlenn': ['madeleine'],
    'mango': ['mnago', 'mango 26', 'mangos', 'mamgo'],
    'monben': ['mombe','momber','monbein','monbin'],
    'naranja': ['naranaja'],
    'neem': ['lila', 'neem/lila'],

    'satanye':['saiteyen', 'santayet', 'santeyet', 'santyet', 'satanyet', 'satayet'],
    'savann': ['bwa saban', 'bois savanne', 'saban'],
    'sed': ['cede'],
    'sikren': ['sicre','sikre', 'sakrin', 'sacrin'],
    'tamarindo': ['tamarenn', 'tamarin'],
    'tcha tcha': ['bois noir', 'nwa', 'nua', 'chacha', 'bwachacha'],
    'trompette': ['trompete'],
    'twa pawol': ['twa parol', 'twa palol','twa pabel','twa pable'],
    'zamann': ['zanmann'],

    'akoma':['bwa koma', 'lakoma'],
    'bande': ['bendi', 'bombe'],
    'camite': ['cawimite', 'cayimit', 'kaymit'],
    'chene': ['bwa chenn', 'chene', 'chenn'],
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
spec_ser

#%% List unique species with number of occurences and write to excel sheet.
spec_list_cnts = spec_ser.value_counts().rename_axis(['common_name'])
specs_mult = spec_list_cnts[spec_list_cnts > 2]
specs_mult.sort_index()
err_specs = spec_list_cnts[spec_list_cnts < 3]

# Run replace again
spec_ser.replace(alt_to_name, inplace=True)
spec_list_cnts = spec_ser.value_counts().rename_axis(['common_name'])

#%% Write to excel
out_lookup = os.path.join(home, 'unique_species_v2_standardized.xlsx')
spec_list_cnts.to_excel(out_lookup, 'round3', index=True)






#%% Load species table
species_fname = os.path.join(home, 'Especies.xlsx')
spec_lookup = pd.read_excel(species_fname, 'Hoja1',
                            names=['common_name', 'common_alt', 'blank', 'genus', 'synonyms'])
common_syns = spec_lookup.loc[:,['common_name', 'common_alt']]
common_syns = common_syns.stack().apply(lambda x : x.strip().lower()).unstack()
comm_names = common_syns.stack().unique().tolist()


# Compare standardized common species names to species names in field data
unlisted_names = []
for sp in spec_list_cnts.index:
    if not sp in comm_names:
        unlisted_names += [sp]

unlisted_names


#%% Read Global Wood Density table
gwd_fpath = r'/Users/emilysturdivant/Documents/Literature/GlobalWoodDensityDatabase.xls'

# Read data, skipping the first three rows with plot level data
gwd_df = pd.read_excel(gwd_fpath, 'Data', index_col=0)
gwd_df.columns
gwd_df.index
com_name = 'bayawonn'
binomial = 'Prosopis juliflora'
spec_gwdnum = 12857
gwd_speci = gwd_df[gwd_df['Binomial'] == binomial]
if len(gwd_speci) > 1:
    print(gwd_speci)


gwd_speci[gwd_speci['Region'] == 'South America (tropical)']
gwd_speci[gwd_speci['Reference Number'] == 110]
gwd_speci[gwd_speci['Reference Number'] == 111]



gwd_speci = gwd_df.loc[spec_gwdnum]
gwd_speci

# Create a table with species, binomial, region, reference number, and wood density
# May have multiple rows for one species (e.g. both central america and south america (tropical) are good options). In that case, we can use .groupby().
# If entries exist for the priority references, use them.


best_refs = [110, 112, 111, 113, 116, 141, 180, 186]
