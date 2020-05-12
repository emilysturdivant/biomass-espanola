# Test to download and process data for Haiti
# Import biota module
import biota
import os

'''
# Download using CLI:
Activate environment with biota installed ('biota')
cd path/to/biota/cli

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

python download.py -lon -72 -lat 21 -y 2018 -r -o /Users/emilysturdivant/PROJECTS/Haiti_biomass/ALOS

python download.py -lon -72 -lat 21 -y 2019 -r -o /Users/emilysturdivant/PROJECTS/Haiti_biomass/ALOS

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
data_dir = r'/Users/emilysturdivant/PROJECTS/biomass-espanola/data/ALOS'
output_dir = r'/Users/emilysturdivant/PROJECTS/biomass-espanola/results/g0nu_HV'

#%% Create list of tile coordinates
coord_list = []
for latitude in range(19, 21):
    for longitude in range(-74, -67):
        coord_list += [[latitude, longitude]]
coord_list += [[18, -72], [19, -75], [21, -73]]

#%% Run processing loop
y1 = 2018
getGamma0_HV_nofilt(data_dir, y1, output_dir, coord_list)

#%% Set slope and intercept of AGB-backscatter regression
slope = 1032
intercept= 0.006948
y1 = 2018
output_dir = r'/Users/emilysturdivant/Documents/CIGA/biota_out/AGB_2018'
getAGB_forRegion(data_dir, y1, slope, intercept, output_dir, coord_list)

#%% Merge and Crop with rasterio. - TRY BELOW FIRST (not tested)
# Rasterio package isn't installing properly in 'biota' env. Installed in 'python3'
import fiona
import rasterio
import rasterio.mask
from rasterio.merge import merge
from rasterio.plot import show
import glob
import os
import numpy as np
from shapely.geometry import shape
from shapely.geometry import box
import geopandas as gpd

# Set inputs
y1= 2018
output_dir = r'/Users/emilysturdivant/PROJECTS/Haiti_biomass/biota_out'
shp_fp = r"/Users/emilysturdivant/PROJECTS/Haiti_biomass/contextual_data/HTI_adm/HTI_adm0.shp"
mosaic_fp = os.path.join(output_dir, f'g0nu_{y1}_HV_haiti.tif')
mosaiccrop_fp = os.path.join(output_dir, f'g0nu_{y1}_HV_haiti2.tif')

# MERGE tiles output by biota
# Open "Gamma0_[year]..." TIFF files in output dir
fps = glob.iglob(os.path.join(output_dir, f'Gamma0_{y1}*.tif'))
src_files_to_mosaic = []
for fp in fps:
    with rasterio.open(fp) as src:
        src_files_to_mosaic.append(src)
# Merge open files and set metadata
mosaic, out_trans = merge(src_files_to_mosaic)

# Set metadata
out_meta = src.meta.copy()
out_meta.update({"driver": "GTiff",
    "height": mosaic.shape[1],
    "width": mosaic.shape[2],
    "transform": out_trans})
# Write mosaic to file
with rasterio.open(mosaic_fp, "w", **out_meta) as dest:
    dest.write(mosaic)

#%% Crop to shapefile bbox
# Get polygon bounding box
with fiona.open(shp_fp, "r") as shapefile:
    shapes = [feature["geometry"] for feature in shapefile]
bb = shape(shapes[0]).bounds
geo = gpd.GeoDataFrame({'geometry': box(bb[0], bb[1], bb[2], bb[3])},
                        index=[0],
                        crs=fiona.crs.from_epsg(4326))
def getFeatures(gdf):
    # Parse features from GeoDataFrame for rasterio
    import json
    return [json.loads(gdf.to_json())['features'][0]['geometry']]
# Crop image to polygons
with rasterio.open(mosaic_fp) as src:
    coords = getFeatures(geo.to_crs(crs=src.crs.data))
    out_img, out_transform = rasterio.mask.mask(dataset=src, shapes=coords, crop=True)
    out_meta = src.meta.copy()
    out_meta.update({"driver": "GTiff",
                    "height": out_img.shape[1],
                    "width": out_img.shape[2],
                    "transform": out_transform}
                    )
    # Write cropped mosaic to file
    with rasterio.open(mosaiccrop_fp, "w", **out_meta) as dest:
        dest.write(out_img)
show(out_img, vmin=0.001, vmax=0.3)

#%% Merge and Crop with rasterio. - NOT TESTED
# Rasterio package isn't installing properly in 'biota' env. Installed in 'python3'
import fiona
import rasterio
import rasterio.mask
from rasterio.merge import merge
from rasterio.plot import show
import glob
import os
import numpy as np
from shapely.geometry import shape
from shapely.geometry import box
import geopandas as gpd
import json

# Set inputs
y1= 2018
output_dir = r'/Users/emilysturdivant/PROJECTS/Haiti_biomass/biota_out'
shp_fp = r"/Users/emilysturdivant/PROJECTS/Haiti_biomass/contextual_data/HTI_adm/HTI_adm0.shp"
mosaic_fp = os.path.join(output_dir, f'g0nu_{y1}_HV_haiti.tif')
mosaiccrop_fp = os.path.join(output_dir, f'g0nu_{y1}_HV_haiti3.tif')

# Get polygon bounding box
with fiona.open(shp_fp, "r") as shapefile:
    shapes = [feature["geometry"] for feature in shapefile]
bb = shape(shapes[0]).bounds
geo = gpd.GeoDataFrame({'geometry': box(bb[0], bb[1], bb[2], bb[3])},
                        index=[0],
                        crs=fiona.crs.from_epsg(4326))

# MERGE tiles output by biota
# Open "Gamma0_[year]..." TIFF files in output dir
fps = glob.iglob(os.path.join(output_dir, f'Gamma0_{y1}*.tif'))
src_files_to_mosaic = []
for fp in fps:
    with rasterio.open(fp) as src:
        src_files_to_mosaic.append(src)
# Merge open files
mosaic, out_trans = merge(src_files_to_mosaic)
# Crop mosaic
geo = geo.to_crs(crs=mosaic.crs.data)
coords = [json.loads(geo.to_json())['features'][0]['geometry']]
out_img, out_transform = rasterio.mask.mask(dataset=mosaic, shapes=coords, crop=True)
out_meta = mosaic.meta.copy()
out_meta.update({"driver": "GTiff",
                "height": out_img.shape[1],
                "width": out_img.shape[2],
                "transform": out_transform})
# Write cropped mosaic to file
with rasterio.open(mosaiccrop_fp, "w", **out_meta) as dest:
    dest.write(out_img)

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

#%% Apply the Radar Enhanced Lee Filter in biota to the whole mosaic.
import numpy as np
from matplotlib import pyplot as plt
import biota
from osgeo import gdal
import os

y1= 2018
output_dir = r'/Users/emilysturdivant/PROJECTS/Haiti_biomass/biota_out'

# Load raster values into a numpy array.
mosaiccrop_fp = os.path.join(output_dir, f'g0nu_{y1}_HV.tif')
g0 = biota.IO.loadArray(mosaiccrop_fp)
plt.imshow(g0, interpolation='nearest')
plt.show()

# Mask 0 values
g0maskd = np.ma.masked_where(g0 <= 0, g0, copy=False)
g0filt = biota.filter.enhanced_lee_filter(g0maskd)
plt.imshow(g0filt, interpolation='nearest')
plt.show()

# Write output
ds = gdal.Open(nofilt_fp, 0)
geo_t = ds.GetGeoTransform()
proj = ds.GetProjection()
filt_fn = f'{os.path.splitext(os.path.basename(nofilt_fp))[0]}_leeBiota2.tif'
biota.IO.outputGeoTiff(g0filt, filt_fn, geo_t, proj, output_dir=output_dir,
    nodata=g0maskd.fill_value)


#%% Convert polygon to raster mask, multiply by mosaic mask => valid area => subtract/inverse
import rasterio
import rasterio.plot
# from rasterio.enums import Resampling
import pyproj
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import shutil

# Load shapefile
with fiona.open(r"/Users/emilysturdivant/PROJECTS/Haiti_biomass/contextual_data/Hispaniola/Hisp_adm0.shp", "r") as shapefile:
    shapes = [feature["geometry"] for feature in shapefile]

# Build overviews
out_fn = f'/Users/emilysturdivant/PROJECTS/Haiti_biomass/biota_out/py_testing/g0nu_2018_HV_leeBiota.tif'
path = shutil.copy(filt_fn, out_fn)
factors = [2, 4, 8, 16]
with rasterio.open(path, 'r+') as dst:
    dst.build_overviews(factors, rasterio.enums.Resampling.average)
    dst.update_tags(ns='rio_overview', resampling='average')

# Try both downsampling and windowing
with rasterio.open(path) as src:
    profile = src.profile.copy()
    oviews = src.overviews(1)
    oview = oviews[-1]
    thumbnail = src.read(1, out_shape=(1, int(src.height // oview), int(src.width // oview)))
    subset = src.read(1, window=rasterio.windows.Window(9000, 12000, 280, 280))

# Plot thumbnail (downsampled) with colorbar and labels
plt.imshow(thumbnail, vmin=0.001, vmax=0.3)
plt.colorbar(shrink=0.5)
plt.title('Overview {}'.format(thumbnail.shape))
plt.xlabel('Column #')
plt.ylabel('Row #')
# Plot Subset (window)
plt.imshow(subset, vmin=0.001, vmax=0.3)
plt.colorbar()
plt.title('Subset {}'.format(thumbnail.shape))
plt.xlabel('Column #')
plt.ylabel('Row #')

# Try both downsampling and windowing
thumb_fn = f'/Users/emilysturdivant/PROJECTS/Haiti_biomass/biota_out/py_testing/g0nu_2018_HV_leeBiota_thumbnail.tif'
with rasterio.open(path) as src:
    profile = src.profile.copy()
    oview = 16
    thumbnail = src.read(1, out_shape=(1, int(src.height // oview), int(src.width // oview)))
    aff = src.transform
    newaff = rasterio.Affine(aff.a * oview, aff.b, aff.c,
                             aff.d, aff.e * oview, aff.f)
    profile.update({
            'dtype': 'float32',
            'height': thumbnail.shape[0],
            'width': thumbnail.shape[1],
            'transform': newaff})
    with rasterio.open(out_fn, 'w', **profile) as dst:
        dst.write_band(1, thumbnail)

out_fn = f'/Users/emilysturdivant/PROJECTS/Haiti_biomass/biota_out/py_testing/mask1.tif'
# Reclassify raster to masked Boolean
with rasterio.open(thumb_fn) as src:
    profile = src.profile.copy()
    # Read as numpy array
    data = src.read(1)
    # Reclassify
    dat2 = data.copy()
    dat2[np.where(dat2 >=0)] = 1
    mask_rst = np.ma.masked_where(data == src.nodata, dat2, copy=False)
    # Export
    with rasterio.open(out_fn, 'w', **profile) as dst:
        dst.write_band(1, thumbnail)
out_fn = f'/Users/emilysturdivant/PROJECTS/Haiti_biomass/biota_out/py_testing/g0nu_2018_HV_leeBiota_maskpoly.tif'
# Mask with inverse of polygon
with rasterio.open(thumb_fn) as src:
    mask_poly, out_transform = rasterio.mask.mask(src, shapes, invert=True)
    out_meta = src.meta
    with rasterio.open(out_fn, 'w', **profile) as dst:
        dst.write(mask_poly)
with rasterio.open(out_fn) as src:
    mask_poly = src.read(1)

# Make plot with colorbar and labels
plt.imshow(mask_rst, vmin=0.001, vmax=0.3)
plt.colorbar()
mask_rst.shape
plt.imshow(mask_poly)
plt.colorbar()
mask_poly.shape

# Testing
src = rasterio.open(filt_fn)
mask_rst = src.read(1, window=rasterio.windows.Window(9000, 12000, 280, 280))
mask_rst = np.ma.masked_where(mask_rst == src.nodata, mask_rst, copy=False)
# mask_rst[np.where(mask_rst <= 0.3)] = -99
# mask_rst[mask_rst > 0] = 1
# Make plot with colorbar and labels
plt.imshow(mask_rst, vmin=0.001, vmax=0.3)
plt.colorbar()

type(mask_rst)
type(np.array(mask_rst).astype(rasterio.float32))
mask_rst.max()
mask_rst.min()
src.nodata

mask_rst.astype(rasterio.float32)


inland_nodata = mask_rst + mask_poly
# 0 = inside poly, but no data in raster
inland_nodata.nodatavals
inland_nodata.nodata
len(inland_nodata.nodata)

# [1,0,1]+[0,0,1] = [1,0,2]
# 0 = inside poly, but no data in raster
# if poly==True and raster==NoData: +1
out_fn = f'/Users/emilysturdivant/PROJECTS/Haiti_biomass/biota_out/py_testing/mask1.tif'
with rasterio.open(out_fn, 'w', **profile) as dst:
    dst.write(array)





#%% Experimental...
latitude=20
longitude=-73
y1=2018
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
