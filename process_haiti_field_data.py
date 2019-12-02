import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

#%% Work with field data - Look at species in all plots
xlsf = r'/Users/emilysturdivant/Downloads/haiti_biomass_v1.xlsx'

# Read plot data, skipping the first three rows with plot level data
plotsdf = pd.read_excel(xlsf, 'Plots', skiprows=[0,1,2], mangle_dupe_cols=True)
spcols = [c for c in plotsdf.columns if c.startswith('sp')]
spec_df = plotsdf.loc[:, spcols]

# List species columns ('sp')
spec_df = pd.read_excel(xlsf, 'Plots', skiprows=[0,1,2], usecols=lambda x : x.startswith('sp'))
spec_ser = spec_df.stack().apply(lambda x : x.strip().lower())

# Standardize species names
spec_ser.replace({'mnago':'mango', 'mango 26':'mango', 'mangos':'mango',
                'aguacates':'aguacate',
                'acassia':'acacia',
                'biosblanc':'boisblanc','boiblanc':'boisblanc','boisblanc2':'boisblanc',
                'buadom':'bwa don'
                'naranaja':'naranja',
                'eucaliptos':'eucalipto',
                'cayaoux':'cayoux',
                'calabase':'calebasse','calbesse':'calebasse', 'calbasse':'calebasse',
                'buayawonn':'bayawonn',
                'latanie':'latanye','latani':'latanye'
                'lamp':'lame',
                'corosol':'corossol','corosole':'corossol','corolosol':'corossol'
                'capable':'kapab',
                'fren':'fwenn',
                'tamarenn':'tamarindo'}, inplace=True)

#%% List unique species with number of occurences and write to excel sheet.
spec_list_cnts = spec_ser.value_counts()
specs_mult = spec_list_cnts[spec_list_cnts > 2]
specs_mult
specs_mult.sort_index()
err_specs = spec_list_cnts[spec_list_cnts < 3]
err_specs
type(spec_list_cnts)

#%% Write to excel
spec_list_cnts = spec_list_cnts.rename_axis(['common_name'])
out_lookup = os.path.join(os.path.dirname(xlsf), 'species_lookuptable.xlsx')
spec_list_cnts.to_excel(out_lookup, 'species_lookup', index=True)




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
