'''
Thoughts on how to perform lookup...
previously at end of standardize_common_names.py
'''
# Options for density lookup table:
# average by family
# manually-selected most relevant density for species (based on region)
# average by species
# Bwa Yo density, if available and if not, most reasonable density from GWD.

# If a binomial has multiple entries, look for


fieldspec_lookup.groupby('creole').agg(uniq)

'''
GLOBAL WOOD DENSITY TABLE
'''
#%% Read Global Wood Density table
gwd_speci[gwd_speci['Region'] == 'South America (tropical)']
gwd_speci[gwd_speci['Reference Number'] == 110]
gwd_speci[gwd_speci['Reference Number'] == 111]
gwd_speci = gwd_df.loc[spec_gwdnum]

# Create a table with species, binomial, region, reference number, and wood density
# May have multiple rows for one species (e.g. both central america and south america (tropical) are good options). In that case, we can use .groupby().
# If entries exist for the priority references, use them.

best_refs = [110, 112, 111, 113, 116, 141, 180, 186]

#%% Load species table (Especies.xlsx by José Luis)
species_fname = os.path.join(home, 'Especies.xlsx')
spec_lookup = pd.read_excel(species_fname, 'Hoja1',
                            names=['common_name', 'common_alt', 'species', 'synonyms'])
common_syns = spec_lookup.loc[:,['common_name', 'common_alt']]
common_syns = common_syns.stack().apply(lambda x : x.strip().lower()).unstack()
comm_names = common_syns.stack().unique().tolist()


# Extract rows from GWD with species present in field data
field_binoms = fieldspec_lookup['by_binomial']
gwd_lookup = gwd_df.loc[gwd_df['gwd_binomial'].isin(field_binoms)]

# %% Experiment with groupby on GWD
count = lambda x: x.count()
gwd_lookup.groupby('gwd_binomial')['gwd_density'].agg(
        mean=np.mean,
        count=count,
        single_entry=single)

# for each field species, get GWD values if there is only one matching GWD entry
single = lambda x: x if x.count() == 1 else np.nan
gwd_lookup.groupby('gwd_binomial').agg(single)
# for each field species, list all unique values in each column
uniq = lambda x: x.unique()
gwd_lookup.groupby('gwd_binomial').agg(uniq)

#%% Join the GWD values to the species lookup table.
gwd_uniq_bybinom = gwd_lookup.groupby('gwd_binomial').agg(uniq)
lookup = fieldspec_lookup.join(gwd_uniq_bybinom, on='by_binomial', how='outer')

#%% Write to excel
out_fname = os.path.join(home, 'lookup_combo_20191220.xlsx')
fieldspec_lookup.to_excel(out_fname, index=False)
'''
'''
#%%
'''
Load data for one plot. Now used as get_plot_data()
'''
# Read in first 4 columns
col_ints = [0,1,2,3]
# col_ints = np.add(col_ints, 4*8) # first plot with trees

# Get plot number, shapefile name, and area for given plot
plot = pd.read_excel(data_fname, 'Plots', header=0, usecols=col_ints, nrows=1)
plot_no = plot.columns[0].split('#')[1]
plot_shp = plot.iloc[0,0]
plot_area = plot.iloc[0,2]
print(f'Plot number: {plot_no} | {plot.iloc[0,0]} | Area: {plot.iloc[0,2]} ha')

# Get data for plot indicated by col_ints
df = pd.read_excel(data_fname, 'Plots', skiprows=[0,1,2], header=0,
    usecols=col_ints, names=['sp_creole','dbh_cm','ht_m','ba_m2'], #dtype={'ht_m':'float'},
    skip_blank_lines = True,
    converters={'sp_creole':lambda x : strip_accents(x.strip().lower()), 'ht_m': lambda x: str(x).strip("''")})
df = df.dropna(axis=0, how='all').astype({'ht_m':'float'})
df.loc[:, 'sp_creole'] = df.sp_creole.replace(alt_to_name)
