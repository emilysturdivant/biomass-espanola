# Test to download and process data for Haiti
# Import biota module
import biota
import os

'''
# Download using CLI:
Activate environment with biota installed ('biota')
# years 2007-2016
biota download -lon -73 -lat 20 -r -o /home/esturdivant/Documents/ALOS

# 2018 (ALOS-2)
for VAR in -68 -69 -70 -71 -72 -73 -74
do
    biota download -lon $VAR -lat 20 -y 2018 -r -o /home/esturdivant/Documents/ALOS
    biota download -lon $VAR -lat 19 -y 2018 -r -o /home/esturdivant/Documents/ALOS
done
biota download -lon -72 -lat 18 -y 2018 -r -o /home/esturdivant/Documents/ALOS
biota download -lon -75 -lat 19 -y 2018 -r -o /home/esturdivant/Documents/ALOS

# Mac:
python download.py -lon -68 -lat 20 -y 2017 -r -o /Users/emilysturdivant/Documents/CIGA/ALOS

'''
#%% functions
def getGamma0_HV_nofilt(data_dir, y1, output_dir=None, coord_list=None, filter=False, polarization='HV'):
    if not output_dir:
        if filter:
            output_dir = data_dir+f'_g0nu_{polarization}lee_{y1}'
        output_dir = data_dir+f'_g0nu_{polarization}_{y1}'
        if os.path.exists(output_dir):
            output_dir += '_2'
    os.makedirs(output_dir, exist_ok=True)
    if not coord_list:
        # Create list of tile coordinates
        coord_list = []
        for latitude in range(19, 21):
            for longitude in range(-74, -67):
                coord_list += [[latitude, longitude]]
        coord_list += [[18, -72], [19, -75]]
    for lat, lon in coord_list:
        # Print progress
        print('Doing latitude: {}, longitude: {}'.format(str(lat), str(lon)))
        # Load the ALOS tile with specified options
        try:
            tile = biota.LoadTile(data_dir, lat, lon, y1, lee_filter = filter, output_dir = output_dir)
        except:
            print('error')
            continue
        # Calculate gamma0 and output to GeoTiff
        gamma0 = tile.getGamma0(polarisation=polarization, output=True)
        # Calculate AGB using slope and intercept from R
        # agb = tile.getAGB(slope=slope, intercept=intercept, output = True)
    print('ALL DONE.')

def getAGB_forRegion(data_dir, y1, slope, intercept, output_dir=None, coord_list=None, filter=False):
    if not output_dir:
        if filter:
            output_dir = data_dir+f'_AGB_{y1}'
        output_dir = data_dir+f'_AGB_{y1}'
    ct = 1
    while os.path.exists(output_dir):
        ct += 1
        output_dir += f'_{ct}'
    os.makedirs(output_dir, exist_ok=False)
    if not coord_list:
        # Create list of tile coordinates
        coord_list = []
        for latitude in range(19, 21):
            for longitude in range(-74, -67):
                coord_list += [[latitude, longitude]]
        coord_list += [[18, -72], [19, -75]]
    for lat, lon in coord_list:
        # Print progress
        print('Doing latitude: {}, longitude: {}'.format(str(lat), str(lon)))
        # Load the ALOS tile with specified options
        try:
            tile = biota.LoadTile(data_dir, lat, lon, y1, lee_filter = filter, output_dir = output_dir)
        except:
            print('error')
            continue
        # Calculate AGB using slope and intercept from R
        agb = tile.getAGB(slope=slope, intercept=intercept, output = True)
    print('ALL DONE.')

#%% Initialize file paths
data_dir = r'/Users/emilysturdivant/Documents/CIGA/ALOS'

#%% Create list of tile coordinates
coord_list = []
for latitude in range(19, 21):
    for longitude in range(-74, -67):
        coord_list += [[latitude, longitude]]
coord_list += [[18, -72], [19, -75]]

#%% Run processing loop
y1 = 2018
output_dir = r'/Users/emilysturdivant/Documents/CIGA/biota_out/AGB_2018'
getGamma0_HV_nofilt(data_dir, y1, output_dir, coord_list)

#%% Set slope and intercept of AGB-backscatter regression
slope = 1032
intercept= 0.006948
y1 = 2018
output_dir = r'/Users/emilysturdivant/Documents/CIGA/biota_out/AGB_2018'
getAGB_forRegion(data_dir, y1, slope, intercept, output_dir, coord_list)

#%%
'''
# In QGIS:
gdal_merge.py -ot Float32 -of GTiff -o /home/esturdivant/Documents/biota_out/AGB_2018_v2/AGB_2018_v2.tif --optfile /tmp/processing_c054cfe79c5c4e7b946b1b6603d98e7b/4e8652a2142740759b28c1fcffea58da/mergeInputFiles.txt

Raster > Miscellaneous > Merge
Input layers: directory X

Display output...
Transparency: No Data value: 0
Symbology: Singlebnad pseudocolor; Cumulative count cut: 2-98%; Color ramp: Viridis. Classify. Apply.
Histogram analysis
'''

#%% Try to work with rasters
import numpy as np
from matplotlib import pyplot as plt
plt.gray()
# Load raster values into a numpy array.
fpath = r'/Users/emilysturdivant/PROJECTS/Haiti_biomass/biota_out/g0nu_2018_nofilt_HV.tif'
g0 = biota.IO.loadArray(fpath)
g0
plt.imshow(g0, interpolation='nearest')
plt.show()


mask_path = r'/Users/emilysturdivant/PROJECTS/Haiti_biomass/LULC/Haiti2017_water.tif'
mask = biota.IO.loadArray(mask_path)
plt.imshow(mask, interpolation='nearest')
plt.show()
mask = biota.IO.loadArray(mask_path) != 255
g0maskd = np.ma.array(g0, mask = mask)
g0filt = biota.filter.enhanced_lee_filter(g0)

#%% Experimental...
tile = biota.LoadTile(data_dir, latitude, longitude, y1, lee_filter = True, forest_threshold = 15., area_threshold = 1, output_dir = output_dir)

# Add river lines to the mask with a 250 m buffer
# tile_2007.updateMask('auxillary_data/TZA_water_lines_dcw.shp', buffer_size = 250)
# Calculate gamma0 and output to GeoTiff
gamma0 = tile.getGamma0(polarisation='HV', output=True)
print('Complete.')

year = 2018
shp = r'/home/esturdivant/code/biomass-espanola/data/AllPlots.shp'
dataloc =

data_dict = extractGamma0(dataloc, year, shp, plot_field='plot_no', agb_field='AGB', buffer_size = 0, verbose = True, units = 'natural')
    # Args:
    #     shp: A shapefile containing plot data
    #     year: Year to load
    #     plot_field: The shapefile field containing plot names
    #     agb_field: The shapefile field containing AGB estimates
    # Returns:
    #     A dictionary containing plot names, AGB and gamma0 values

(slp, intcpt) = fitLinearModel(data_dict)
'''
Fit a linear model to relate AGB to gamma0 backscatter.

Args:
    data_dict: Dictionary output from extractGamma0() function.
Returns:
    model slope, model intercept
'''

# Calculate AGB and output to GeoTiff
# gamma0_2007 = tile.getAGB(output = True)
#
# # Calculate Woody cover and output to GeoTiff
# gamma0_2007 = tile.getWoodyCover(output = True)
