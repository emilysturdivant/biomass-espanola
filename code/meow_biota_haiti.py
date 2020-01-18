# Test to download and process data for Haiti
# Import biota module
import biota


'''
# Download using CLI:
Activate environment with biota installed ('biota')
# years 2007-2016
biota download -lon -73 -lat 20 -r -o /home/esturdivant/Documents/ALOS

# years 2017
biota download -lon -73 -lat 20 -y 2018 -r -o /home/esturdivant/Documents/ALOS

# years 2018 and 2019 result in error
biota download -lon -73 -lat 20 -y 2018 -r -o /home/esturdivant/Documents/ALOS
biota download -lon -74 -lat 20 -y 2018 -r -o /home/esturdivant/Documents/ALOS
biota download -lon -72 -lat 20 -y 2018 -r -o /home/esturdivant/Documents/ALOS
biota download -lon -73 -lat 19 -y 2018 -r -o /home/esturdivant/Documents/ALOS

for VAR in -69 -73 -74 -75
do
    biota download -lon $VAR -lat 20 -r -o /home/esturdivant/Documents/ALOS
    biota download -lon $VAR -lat 20 -y 2017 -r -o /home/esturdivant/Documents/ALOS
done

for VAR in -69 -73 -74 -75
do
    biota download -lon $VAR -lat 20 -y 2018 -r -o /home/esturdivant/Documents/ALOS
done
biota download -lon -72 -lat 18 -y 2018 -r -o /Users/emilysturdivant/Documents/CIGA/ALOS
'''


# Process data
data_dir = r'/home/esturdivant/Documents/ALOS'
output_dir = r'/home/esturdivant/Documents/biota_out/g0nu_2018_HV_lee'

data_dir = r'/Users/emilysturdivant/Documents/CIGA/ALOS'
output_dir = r'/Users/emilysturdivant/Documents/CIGA/biota_out/AGB_2017_v1'

# Set slope and intercept of AGB-backscatter regression
slope = 2426.26
intercept= 10.21

# Create list of tile coordinates
coord_list = []
for latitude in range(19, 21):
    for longitude in range(-74, -67):
        coord_list += [[latitude, longitude]]
coord_list += [[18, -72], [19, -75]]

y1 = 2017
for lat, lon in coord_list:
    # Print progress
    print('Doing latitude: {}, longitude: {}'.format(str(lat), str(lon)))
    # Load the ALOS tile with specified options
    try:
        tile = biota.LoadTile(data_dir, lat, lon, y1, lee_filter = True, output_dir = output_dir)
    except:
        print('error')
        continue
    # Calculate gamma0 and output to GeoTiff
    # gamma0 = tile.getGamma0(polarisation='HV', output=True)
    agb = tile.getAGB(slope= 2426.26, intercept= 10.21, output = True)
    print('Complete.')



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
