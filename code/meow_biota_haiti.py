# Test to download and process data for Haiti
# Import biota module
import biota

# Download data using GUI
    # Input_folder = /home/esturdivant/code/meows/haiti/DATA
    # Output_folder = /home/esturdivant/code/meows/haiti/outputs
    # Latitude = -5 deg
    # Longitude = 35 deg
    # Tile_size = 1x1
    # Forest_properties = True
    # Forest_change = False
    # Year 1 = 2007
    # Year 2 = 2010
    # Gamma0 = False
    # AGB = True
    # Forest_cover = False
    # Biomass_change = False
    # Forest_cover_change = False
    # Deforestation_risk = False
    # Filter = True
    # Resample_factor = 1
    # Polarisation = HV
    # Area_threshold = 1
    # Biomass_threshold = 10
    # Area_change_threshold = 0
    # Biomass_change_treshold = 0
    # Percentage_change_threshold = 0
'''
# Download using CLI:
# years 2007-2016
biota download -lon -73 -lat 20 -r -o /home/esturdivant/Documents/ALOS

# years 2017
biota download -lon -73 -lat 20 -y 2017 -r -o /home/esturdivant/Documents/ALOS

# years 2018 and 2019 result in error
'''

# Process data
data_dir = r'/home/esturdivant/Documents/ALOS'
output_dir = r'/home/esturdivant/Documents/biota_out/N20W073_gamma0_leefilter'
latitude = 20 # maxlat
longitude = -73 # minlon
y1 = 2010
# y2 = 2018

# Load tile from downloaded data
tile_y1 = biota.LoadTile(data_dir, latitude, longitude, y1, lee_filter=True, output_dir=output_dir)

# Look at tile properties
tile_y1.year
tile_y1.lat
tile_y1.lon
tile_y1.directory
tile_y1.satellite
tile_y1.xSize, tile_y1.ySize
tile_y1.xRes, tile_y1.yRes
tile_y1.extent # minlon, minlat, maxlon, maxlat
tile_y1.geo_t # A geo_transform QMetaObject
tile_y1.proj # Projection wkt

# Extract backscatter information and save to GeoTiff
gamma0_y1 = tile_y1.getGamma0(polarisation='HV', units='decibels', output=True, show=True)


output_dir = r'/home/esturdivant/Documents/biota_out/gamma0_2015_HV_lee'
latrange = range(20,21)
lonrange = range(-74, -67)
y1 = 2015
for latitude in latrange:
    for longitude in lonrange:
        # Print progress
        print('Doing latitude: {}, longitude: {}'.format(str(latitude), str(longitude)))
        # Load the ALOS tile with specified options
        try:
            tile = biota.LoadTile(data_dir, latitude, longitude, y1, lee_filter = True, forest_threshold = 15., area_threshold = 1, output_dir = output_dir)
        except:
            print('error')
            continue
        # Add river lines to the mask with a 250 m buffer
        # tile_2007.updateMask('auxillary_data/TZA_water_lines_dcw.shp', buffer_size = 250)
        # Calculate gamma0 and output to GeoTiff
        gamma0 = tile.getGamma0(polarisation='HV', units='decibels', output=True)
        print('Complete.')

        # Calculate AGB and output to GeoTiff
        # gamma0_2007 = tile.getAGB(output = True)
        #
        # # Calculate Woody cover and output to GeoTiff
        # gamma0_2007 = tile.getWoodyCover(output = True)

tile = biota.LoadTile(data_dir, latitude, longitude, y1, lee_filter = True, forest_threshold = 15., area_threshold = 1, output_dir = output_dir)


# for VAR in -69 -73 -74 -75
# do
#     biota download -lon $VAR -lat 20 -r -o /home/esturdivant/Documents/ALOS
#     biota download -lon $VAR -lat 20 -y 2017 -r -o /home/esturdivant/Documents/ALOS
# done
