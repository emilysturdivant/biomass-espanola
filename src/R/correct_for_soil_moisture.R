# ---------------------------------------------------------------------------------------------
# Script to:
#     * Correct AGB/g0 for soil moisture
# Proceeds:
#     * 
# Requires:
#     * ESA soil moisture product from http://www.esa-soilmoisture-cci.org, downloaded via FTP
# ---------------------------------------------------------------------------------------------

# Load libraries
library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(ggplot2) # package for plotting
library(tidyverse)
library(tools)
library(stars)
library(rhdf5)
library(reshape)
# library(tmap)
# tmap_mode('view')
library(rasterVis)
library(smapr)

# Process Surface/Rootzone Soil Moisture Analysis Update (SPL4SMAU) using smapr -----------------
# https://nsidc.org/data/SPL4SMAU/versions/4
# date accessed: 2020-06-19
get_smap_rasters <- function(dates, dir='data/SoilMoisture/SMAP', id='SPL4SMAU', 
                             version=4, name='/Analysis_Data/sm_surface_analysis'){
  # Download files
  files <- find_smap(id, dates = dates, version = version)
  local_files <- download_smap(files, overwrite = FALSE, verbose = T)
  sm_raster <- extract_smap(local_files, name)
  names(sm_raster) <- names(sm_raster) %>% str_extract("T[0-9]{4,}")
  
  # Crop to Hispaniola
  ext <- raster('results/g0nu_HV/g0nu_2018_HV.tif') %>% 
    projectExtent(crs(sm_raster[[1]])) %>% 
    extent()
  sm_raster <- crop(sm_raster, ext)
  
  # Save
  d <- dates %>% str_replace_all('-','')
  n <- name %>% basename()
  out_fp <- file.path(dir, id, str_glue('{n}_{d}.grd'))
  newr <- sm_raster %>% writeRaster(out_fp, 
                                    format='raster', overwrite=T)
  hdr(newr, format = "ENVI")
  return(newr)
}
get_stack_means <- function(in_fp){
  # Load stack
  smr <- stack(in_fp)
  # Calculate mean
  mean_sm <- calc(smr, fun=mean)
  d <- in_fp %>% str_extract('[0-9]{8,}')
  names(mean_sm) <- str_glue('D{d}')
  return(mean_sm)
}

dates <- c('2015-08-30', '2016-09-02', '2016-09-08', '2016-09-09', 
           '2016-09-13', '2016-09-18', '2016-09-30', '2017-06-14', 
           '2017-09-14', '2017-09-19', '2017-09-24', '2018-11-13')
dates %>% lapply(get_smap_rasters)
# Load
smr <- stack('data/SoilMoisture/SMAP/SPL4SMAU/sm_surface_analysis_20160908.grd')
levelplot(smr)
citation("smapr")

# Get daily means
fps <- list.files(path='data/SoilMoisture/SMAP/SPL4SMAU', pattern='grd$',
                  full.names=T, recursive=F, include.dirs=F)
sm_means <- fps %>% lapply(get_stack_means)
sm_means <- stack(sm_means)
# Save raster stack 
out_fp <- file.path('data/SoilMoisture/SMAP/SPL4SMAU', str_glue('daily_means.grd'))
newr <- sm_means %>% writeRaster(out_fp, format='raster', overwrite=T)
hdr(newr, format = "ENVI")

# Look
levelplot(newr)

sm_means <- stack(file.path('data/SoilMoisture/SMAP/SPL4SMAU', str_glue('daily_means.grd')))
levelplot(sm_means)











# ESA soil moisture data ------------------------------------------------------------
convert_sm_netcdf_to_gtif <- function(in_fp, extent, var='sm'){
  # Set output filename
  b <- in_fp %>% basename() %>% str_split_fixed('-', n=6)
  b <- b[[1,6]] %>% str_sub(1, 8)
  out_fp <- file.path(dirname(in_fp), str_glue('ESA_SM_COMBINED_{b}_{var}_hisp.tif'))
  # Open file
  nc_data <- ncdf4::nc_open(in_fp)
  # Read the array 
  sm.array <- ncdf4::ncvar_get(nc_data, var) # store the data in a 3-dimensional array
  # Read the lon, lat, time, and fill value
  lon <- ncdf4::ncvar_get(nc_data, "lon")
  lat <- ncdf4::ncvar_get(nc_data, "lat")
  t <- ncdf4::ncvar_get(nc_data, "time")
  # Replace fill with NA
  fillvalue <- ncdf4::ncatt_get(nc_data, var, "_FillValue")
  sm.array[sm.array == fillvalue$value] <- NA
  # Convert array to raster
  r <- raster(t(sm.array), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), 
              crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  # Crop to Hispaniola
  rc <- crop(r, extent)
  writeRaster(rc, out_fp, "GTiff", overwrite=TRUE)
}
convert_netcdf_to_geotiff <- function(in_fp, extent, var){
  array <- ncdf4::ncvar_get(nc_data, var) # store the data in a 3-dimensional array
  # Read the lon, lat, time, and fill value
  lon <- ncdf4::ncvar_get(nc_data, "lon")
  lat <- ncdf4::ncvar_get(nc_data, "lat")
  # Convert array to raster
  r <- raster(t(array), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), 
              crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  # Set output filename
  fn <- file_path_sans_ext(basename(in_fp))
  out_fp <- file.path(dirname(in_fp), str_glue('{fn}_{var}_hisp.tif'))
  # Crop to Hispaniola
  rc <- crop(r, extent)
  # Write raster
  writeRaster(rc, out_fp, "GTiff", overwrite=TRUE)
}

# Check out the data -------------------------------------------------------------------------
# Downloaded netCDF files of soil moisture and ancillary rasters by registering with the ESA.
# host: ftp.geo.tuwien.ac.at; username: esacci_sm_v047; password: zN1xr3apWhEE; port: 22
# Downloaded using FileZilla
nc_data <- ncdf4::nc_open('data/SoilMoisture/ESACCI-SOILMOISTURE-L3S-SSMV-COMBINED-20160902000000-fv04.7.nc')
# Save the print(nc) dump to a text file
{
  sink('data/SoilMoisture/ESACCI-SOILMOISTURE-L3S-SSMV-COMBINED-20160902000000-fv04.7.txt')
  print(nc_data)
  sink()
}
# Verify data, dimensions, attributes
lon <- ncdf4::ncvar_get(nc_data, "lon")
head(lon)
sm.array <- ncdf4::ncvar_get(nc_data, "sm") # store the data in a 3-dimensional array
dim(sm.array)  # lon 1440 x lat 720
ncdf4::nc_close(nc_data) 

# Convert NetCDFs to cropped GeoTiffs -------------------------------------------------------------------------------
ext <- raster('results/g0nu_HV/g0nu_2018_HV.tif') %>% extent()
fps <- list.files(path='data/SoilMoisture',
                  pattern='nc$',
                  full.names=T,
                  recursive=F,
                  include.dirs=F)
rs <- fps %>% lapply(convert_sm_netcdf_to_gtif, extent=ext,  var='sm')
rs <- stack(rs)

# Ancillary datasets
rs_daynight <- fps %>% lapply(convert_sm_netcdf_to_gtif, extent=ext,  var='dnflag')
rs_uncertainty <- fps %>% lapply(convert_sm_netcdf_to_gtif, extent=ext,  var='sm_uncertainty')
rs_flag <- fps %>% lapply(convert_sm_netcdf_to_gtif, extent=ext,  var='flag') # I don't see any pixels flagged as just dense vegetation. 

r1 <- rs[[3]] - rs[[2]]
plot(rs_flag[[7]]) 
plot(rs[[2]])

tmap_mode('view')
tm_shape(rs[[1]]) + tm_raster() +
  tm_shape(rs[[2]]) + tm_raster() +
  tm_shape(rs[[3]]) + tm_raster() +
  tm_shape(rs[[4]]) + tm_raster() +
  tm_shape(rs[[5]]) + tm_raster() +
  tm_shape(rs[[6]]) + tm_raster() +
  tm_shape(rs[[7]]) + tm_raster()


# Ancillary files ---------------------------------------------------------------------------------------
# Land and Rainforest mask
in_fp <- 'data/SoilMoisture/ancillary/ESA-CCI-SOILMOISTURE-LAND_AND_RAINFOREST_MASK-fv04.2.nc'
# Open file
nc_data <- ncdf4::nc_open(in_fp)
{
  sink(str_c(file_path_sans_ext(in_fp),'txt', sep='.'))
  print(nc_data)
  sink()
}
rc <- convert_netcdf_to_geotiff(in_fp, ext, var = 'land')
plot(rc)
rc <- convert_netcdf_to_geotiff(in_fp, ext, var = 'rainforest')
plot(rc)

# Topographic complexity
in_fp <- 'data/SoilMoisture/ancillary/ESACCI-SOILMOISTURE-TOPOGRAPHIC_COMPLEXITY_V01.1.nc'
# Open file
nc_data <- ncdf4::nc_open(in_fp)
{
  sink(str_c(file_path_sans_ext(in_fp),'txt', sep='.'))
  print(nc_data)
  sink()
}
rc <- convert_netcdf_to_geotiff(in_fp, ext, var = 'topographic_complexity')
plot(rc)

# Wetland Fraction
in_fp <- 'data/SoilMoisture/ancillary/ESACCI-SOILMOISTURE-WETLAND_FRACTION_V01.1.nc'
# Open file
nc_data <- ncdf4::nc_open(in_fp)
{
  sink(str_c(file_path_sans_ext(in_fp),'txt', sep='.'))
  print(nc_data)
  sink()
}
rc <- convert_netcdf_to_geotiff(in_fp, ext, var = 'wetland_fraction')
plot(rc)

# Porosity
in_fp <- 'data/SoilMoisture/ESA/ancillary/ESACCI-SOILMOISTURE-POROSITY_V01.1.nc'
# Open file
nc_data <- ncdf4::nc_open(in_fp)
{
  sink(str_c(file_path_sans_ext(in_fp),'txt', sep='.'))
  print(nc_data)
  sink()
}
rc <- convert_netcdf_to_geotiff(in_fp, ext, var = 'porosity')
plot(rc)

# Perform interpolation  ---------------------------------------------------------------------------------------
# If files were already created:
fps <- list.files(path='data/SoilMoisture/ESA',
                  pattern='sm_hisp.tif$',
                  full.names=T,
                  recursive=F,
                  include.dirs=F)
rs <- stack(fps)


rs_df <- as.data.frame(rs, xy = TRUE) %>%
  melt(id.vars = c('x','y'))

ggplot() +
  geom_raster(data = rs_df , aes(x = x, y = y, fill = value)) +
  facet_wrap(~ variable)


# SMAP data  ---------------------------------------------------------------------------------------
convert_sm_hdf5_to_gtif <- function(in_fp){
  # Set output filename
  b <- in_fp %>% basename() %>% str_split('_')
  d <- b[[1]][[6]] %>% str_sub(1,8)
  tile <- b[[1]][[8]]
  out_fp <- file.path(dirname(in_fp), str_glue('SMAP_sm1km_{d}_{tile}.tif'))
  # Get data extent
  met <- h5readAttributes(in_fp, "/Metadata/Extent")
  ext <- extent(met$polygonPosList[2], met$polygonPosList[4],
                met$polygonPosList[5], met$polygonPosList[1])
  # Get data, replace fills with NA, transpose and convert to raster
  sm <- h5read(in_fp,"/Soil_Moisture_Retrieval_Data_1km/soil_moisture_1km")
  sm[sm==-9999] <- NA
  smr <- raster(t(sm), crs=crs(paste0("+init=epsg:",4326)))
  extent(smr) <- ext
  # Save as GeoTIFF
  writeRaster(smr, out_fp, overwrite=T)
}

# SPL2SMAP_S: SMAP/Sentinel-1 L2 Radiometer/Radar 30-Second Scene 3 km EASE-Grid Soil Moisture, Version 2
fps <- list.files(path='data/SoilMoisture/SMAP/3km_v2',
                  pattern='h5$',
                  full.names=T,
                  recursive=F,
                  include.dirs=F)
rs <- fps %>% lapply(convert_sm_hdf5_to_gtif)
fps <- list.files(path='data/SoilMoisture/SMAP/3km_v2/2017',
                  pattern='h5$',
                  full.names=T,
                  recursive=F,
                  include.dirs=F)
rs <- fps %>% lapply(convert_sm_hdf5_to_gtif)

# Take a look
tm_shape(rs[[1]]) + tm_raster()+
  tm_shape(rs[[2]]) + tm_raster()+
  tm_shape(rs[[3]]) + tm_raster()+
  tm_shape(rs[[4]]) + tm_raster()+
  tm_shape(rs[[5]]) + tm_raster()+
  tm_shape(rs[[6]]) + tm_raster()+
  tm_shape(rs[[7]]) + tm_raster()+
  tm_shape(rs[[8]]) + tm_raster()+
  tm_shape(rs[[9]]) + tm_raster()+
  tm_shape(rs[[10]]) + tm_raster()



