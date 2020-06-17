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

# Function ------------------------------------------------------------
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

library(tmap)
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
in_fp <- 'data/SoilMoisture/ancillary/ESACCI-SOILMOISTURE-POROSITY_V01.1.nc'
# Open file
nc_data <- ncdf4::nc_open(in_fp)
{
  sink(str_c(file_path_sans_ext(in_fp),'txt', sep='.'))
  print(nc_data)
  sink()
}
rc <- convert_netcdf_to_geotiff(in_fp, ext, var = 'porosity')
plot(rc)
