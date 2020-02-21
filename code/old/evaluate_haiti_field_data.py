# -*- coding: utf-8 -*-
"""
by: Emily Sturdivant, emilysturdivant@gmail.com

OVERVIEW: Getting started with processing Haiti field data
"""
# Import packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Set working directory
home = r'/Users/emilysturdivant/GitHub/biomass-espanola'

#%% Work with field data - Look at species in all plots
data_fname = os.path.join(home, 'haiti_biomass_v2.xlsx')

#%% Create series of all species columns (labeled 'sp')
spec_df = pd.read_excel(data_fname, 'Plots', skiprows=[0,1,2],
                        usecols=lambda x : x.startswith('sp'))
spec_ser = spec_df.stack().apply(lambda x : x.strip().lower())

#%% List unique species with number of occurences and write to excel sheet.
spec_list_cnts = spec_ser.value_counts().rename_axis(['common_name'])
specs_mult = spec_list_cnts[spec_list_cnts > 2]
specs_mult.sort_index()
err_specs = spec_list_cnts[spec_list_cnts < 3]

#%% Write to excel
out_lookup = os.path.join(os.path.dirname(data_fname), 'species_inPlots_v2.xlsx')
spec_list_cnts.to_excel(out_lookup, 'v2', index=True)

# Standardize species names
# spec_ser.replace({'mnago':'mango', 'mango 26':'mango', 'mangos':'mango',
#                 'aguacates':'aguacate',
#                 'acassia':'acacia',
#                 'biosblanc':'boisblanc','boiblanc':'boisblanc','boisblanc2':'boisblanc',
#                 'buadom':'bwa don'
#                 'naranaja':'naranja',
#                 'eucaliptos':'eucalipto',
#                 'cayaoux':'cayoux',
#                 'calabase':'calebasse','calbesse':'calebasse', 'calbasse':'calebasse',
#                 'buayawonn':'bayawonn',
#                 'latanie':'latanye','latani':'latanye'
#                 'lamp':'lame',
#                 'corosol':'corossol','corosole':'corossol','corolosol':'corossol'
#                 'capable':'kapab',
#                 'fren':'fwenn',
#                 'tamarenn':'tamarindo'}, inplace=True)




#%% Compare v1 and v2
# File name
data_fname1 = os.path.join(home, 'haiti_biomass_v1.xlsx')
# Read plot data, skipping the first three rows with plot level data
plotsdf2 = pd.read_excel(data_fname, 'Plots', skiprows=[0,1,2], mangle_dupe_cols=True)
plotsdf1 = pd.read_excel(data_fname1, 'Plots', skiprows=[0,1,2], mangle_dupe_cols=True)
# Return differences
df = pd.concat([plotsdf1, plotsdf2])
df = df.reset_index(drop=True)
df_gpby = df.groupby(list(df.columns))
idx = [x[0] for x in df_gpby.groups.values() if len(x) == 1]
df.reindex(idx)

# Compare
spec_list_v2 = spec_list_cnts
spec_df_v1 = pd.read_excel(data_fname1, 'Plots', skiprows=[0,1,2],
                        usecols=lambda x : x.startswith('sp'))
spec_ser_v1 = spec_df_v1.stack().apply(lambda x : x.strip().lower())
spec_list_cnts_v1 = spec_ser_v1.value_counts().rename_axis(['common_name'])

#%% Write to excel
out_lookup = os.path.join(home, 'species_inPlots_v1.xlsx')
spec_list_cnts_v1.to_excel(out_lookup, 'v1', index=True)

spec_list_cnts_v1.sort_index()
spec_list_cnts_v1

spec_list_v2.sort_index()


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
