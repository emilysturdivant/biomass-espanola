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



#%% Basal Area: the area of land that is occupied by the cross-section of a tree. (silvR)
# BA = π(DBH/2)^2

#%% Stem volume: an estimate of the amount of space that the tree bole occupies. (silvR)
# Calculate with equation from unknown source


#%%
# Read in lookup table (relates common creole name to wood density)
dens_lookup = pd.read_excel(lookup_fname)
dfdens = df.join(dens_lookup.set_index('creole'), on='sp_creole', how='left')
dfdens.dtypes

dfdens

#%% Perform allometric calculation for plot

# Use Model 4 from Chave et al. (2014) when we have height values
# select rows with height values
trees_wHt = dfdens.loc[~dfdens['ht_m'].isna(), :]
agb_wHt = lambda x: 0.0673 * (x['gwd_density'] * x['dbh_cm']**2 * x['ht_m'])**0.976
# trees_wHt.loc[:, ['biomass']] = trees_wHt.apply(agb_wHt, axis=1)
biomass_wHt = trees_wHt.assign(agb=agb_wHt)
biomass_wHt


trees_noHt = dfdens[dfdens['ht_m'].isna()]
agb_noHt =
biomass_noHt = trees_noHt.assign(biomass=agb_noHt)
biomass_noHt












#%% Get value for E (combination of TS, CWD, and PS) for Hispaniola
# Downloaded CWD from http://chave.ups-tlse.fr/pantropical_allometry.htm
# Downloaded TS and PS from http://worldclim.org/version2
# Working through this https://www.datacamp.com/community/tutorials/geospatial-data-python

# Load rasters
cwd_fname = os.path.join(home, 'data', 'CWD.tif')
ts_fname = os.path.join(home, 'data', 'wc2.0_bio_2.5m_04.tif')
ps_fname = os.path.join(home, 'data', 'wc2.0_bio_2.5m_15.tif')
cwd_world = cwd_fname
ts_world = ts_fname
ps_world = ps_fname

# Isolate to hispaniola
# Load and merge haiti and DR shapefiles

# filenames
haiti_fname = os.path.join(home, 'data', 'HTI_adm0.shp')
dr_fname = os.path.join(home, 'data', 'DOM_adm0.shp')

# import fiona
#
# # Load and check data
# haiti = fiona.open(haiti_fname)
# print(haiti.schema)
# print(haiti.next()) # (GeoJSON format) {'geometry': {'type': 'LineString', 'coordinates': [(0.0, 0.0), (25.0, 10.0), (50.0, 50.0)]}, 'type': 'Feature', 'id': '0', 'properties': OrderedDict([(u'FID', 0.0)])}
# dr = fiona.open(dr_fname)
# print(dr.schema) # {'geometry': 'LineString', 'properties': OrderedDict([(u'FID', 'float:11')])}
# #first feature of the shapefile
# print(dr.next())



import fiona

haiti = gpd.read_file(haiti_fname)
haiti.head()
haiti.columns
haiti.plot()
dr = gpd.read_file(dr_fname)
print(dr.head())
dr.describe

# Merge
hispaniola = gpd.sjoin(haiti, dr, how='inner', op='intersects')
hispaniola
hispaniola.plot()

hispaniola.dissolve(by='EU')


# Mask
cwd_hisp =
ts_hisp =
ps_hisp =

# Calculate E
E = 1.e-3 * (0.178*ts_hisp - 0.938*cwd_hisp - 6.61*ps_hisp)


agb_noHt = exp(-1.803 - 0.976*E + 0.976*ln(wood_dens) + 2.673 ln(dbh) - 0.0299 * ln(dbh)**2)
