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
# SPL4SMAU: Download, crop, and average Surface/Rootzone Soil Moisture Analysis Update -----------------
id <- 'SPL4SMAU'
dates <- c('2015-08-30', '2016-09-02', '2016-09-08', '2016-09-09', 
           '2016-09-13', '2016-09-18', '2016-09-30', '2017-06-14', 
           '2017-09-14', '2017-09-19', '2017-09-24', '2018-11-13')

# Download, crop, and save as multiband GRD rasters
dates %>% lapply(get_smap_rasters, id=id, name='/Analysis_Data/sm_surface_analysis')

# Load
date <- '20160902'
smr <- stack(file.path('data/SoilMoisture/SMAP', id, str_glue('sm_surface_analysis_{date}.grd')))
levelplot(smr, main=id)
citation("smapr")

# Get daily means
fps <- list.files(path='data/SoilMoisture/SMAP/SPL4SMAU', pattern='grd$',
                  full.names=T, recursive=F, include.dirs=F)
sm_means <- fps %>% lapply(get_stack_means)
sm_means <- stack(sm_means)
levelplot(sm_means, main=id)

# Save raster stack 
out_fp <- file.path('data/SoilMoisture/SMAP', id, str_glue('daily_means.grd'))
newr <- sm_means %>% writeRaster(out_fp, format='raster', overwrite=T)
hdr(newr, format = "ENVI")

# Load and look
sm_means <- stack(file.path('data/SoilMoisture/SMAP', id, str_glue('daily_means.grd')))
levelplot(sm_means, main=id)

# SPL4SMGP: Download, crop, and average Surface/Rootzone Soil Moisture Analysis Update -----------------
id <- 'SPL4SMGP'
dates <- c('2015-08-30', '2016-09-02', '2016-09-08', '2016-09-09', 
           '2016-09-13', '2016-09-18', '2016-09-30', '2017-06-14', 
           '2017-09-14', '2017-09-19', '2017-09-24', '2018-11-13')
dates %>% lapply(get_smap_rasters, id=id, name='/Geophysical_Data/sm_surface')

# Load
date <- '20160908'
smr <- stack(file.path('data/SoilMoisture/SMAP', id, str_glue('sm_surface_{date}.grd')))
levelplot(smr, main=id)
citation("smapr")

# Get daily means
fps <- list.files(path=str_glue('data/SoilMoisture/SMAP/{id}'), pattern='grd$',
                  full.names=T, recursive=F, include.dirs=F)
sm_means <- fps %>% lapply(get_stack_means)
sm_means <- stack(sm_means)
levelplot(sm_means, main=id)

# Save raster stack 
out_fp <- file.path('data/SoilMoisture/SMAP', id, str_glue('{id}_daily_means.grd'))
sm_means <- sm_means %>% writeRaster(out_fp, format='raster', overwrite=T)
hdr(sm_means, format = "ENVI")

# Look
sm_means <- stack(file.path('data/SoilMoisture/SMAP', id, str_glue('{id}_daily_means.grd')))
levelplot(sm_means, main=id)










