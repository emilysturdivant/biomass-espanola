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

# Get data for plot indicated by col_ints
df = pd.read_excel(data_fname, 'Plots', skiprows=[0,1,2], header=0, usecols=col_ints, names=['sp','dap_cm','ht_m','ba_m2'], converters={'sp':lambda x : strip_accents(x.strip().lower())})
df.head()


#%% Look at species table digitized from Bwa Yo
by_fname = os.path.join(home, 'bwayo_species.xlsx')
by_df = pd.read_excel(by_fname, header=0, usecols=[0,1,2,3,5], names=['BY_binomial', 'creole', 'BY_spec_grav', 'family', 'synonyms'], converters={'BY_binomial':lambda x : x.lower(), 'family':lambda x : x.split(' (')[0].lower()})

#%% Explode DF by the creole names column. Convert values to list and for every row with multiple names, create multiple species columns.
by_df = by_df.assign(creole=by_df.creole.str.split(','))
by_df = by_df.explode('creole')







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
