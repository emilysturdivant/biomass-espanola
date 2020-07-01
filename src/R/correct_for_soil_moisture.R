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
library(tmap)
tmap_mode('view')
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
download_crop_average_smap_sm <- function(id, name, dates, folder='data/SoilMoisture/SMAP/'){
  # Download, crop, and save as multiband GRD rasters
  dates %>% lapply(get_smap_rasters, id=id, name=name)
  
  # Get daily means
  fps <- list.files(path=file.path(folder, id), pattern='grd$',
                    full.names=T, recursive=F, include.dirs=F)
  sm_means <- fps %>% lapply(get_stack_means)
  sm_means <- stack(sm_means)
  
  # Save raster stack 
  out_fp <- file.path(folder, id, str_glue('{id}_daily_means.grd'))
  sm_means <- sm_means %>% writeRaster(out_fp, format='raster', overwrite=T)
  hdr(sm_means, format = "ENVI")
  
  # Look
  sm_means <- stack(file.path(folder, id, str_glue('{id}_daily_means.grd')))
  levelplot(sm_means, main=id)
  return(sm_means)
}

# Download, crop, and average SMAP soil moisture products for every date -----------------
dates <- c('2015-08-30', '2016-09-02', '2016-09-08', '2016-09-09', 
           '2016-09-13', '2016-09-18', '2016-09-30', '2017-06-14', 
           '2017-09-14', '2017-09-19', '2017-09-24', '2018-11-13')

# SPL4SMGP: Surface/Rootzone Soil Moisture Analysis Update
id <- 'SPL4SMGP'
name <- '/Geophysical_Data/sm_surface'
geophys_means <- download_crop_average_smap_sm(id, name, dates)
levelplot(geophys_means, main=id)

smr <- stack(file.path(folder, id, str_glue('sm_surface_20160908.grd')))
levelplot(smr, main=id)

# SPL4SMAU: Surface/Rootzone Soil Moisture Analysis Update
id <- 'SPL4SMAU'
name <- '/Analysis_Data/sm_surface_analysis'
analysis_means <- download_crop_average_smap_sm(id, name, dates)
levelplot(analysis_means, main=id)

smr <- stack(file.path(folder, id, str_glue('sm_surface_20160908.grd')))
levelplot(smr, main=id)

# Compare SMAP images to ALOS mosaic -----------------
# Load AGB 
agb_fp <- 'results/tifs_by_R/agb18_v1_l1_mask_Ap3WUw25.tif'
sm_fp <- file.path('data/SoilMoisture/SMAP/', 'SPL4SMAU', 'SPL4SMAU_daily_means.grd')
date_fp <- 'results/tifs_by_R/hisp18_date.tif'
date_fp2 <- 'results/tifs_by_R/hisp18_date_r-smap.tif'

a_date <- raster(date_fp)
agb <- raster(agb_fp)


library(gdalUtils)
# Resample dates mosaic to SMAP
gdalwarp(srcfile=date_fp, dstfile=date_fp2, s_srs='EPSG::4326', t_srs=crs(sm),
         te=extent(sm), tr=c(xres(sm), yres(sm)), tap=T, overwrite=T)
library(rgdal)
OGRSpatialReference.SetFromUserInput('EPSG:4326')

# Create mosaic of SMAP by mapping to dates



# - Map date codes to dates
d <- c(463, 832, 836, 837, 841, 846, 860, # Unique date codes
       1117, 1209, 1214, 1219, 1634)
as.Date(d, '2014-05-24')
as.numeric(as.Date('20140524', "%Y%m%d")) - as.numeric(as.Date('2014-05-24'))





sm_means <- stack(sm_fp)
sm <- sm_means$D20160902



crs(sm)
crs(agb)
sm_proj <- projectRaster(sm, crs=crs(agb))
crs(sm_proj)
extent(sm_proj)

test_bb <- extent(c(xmin=-72.44, xmax=-72.42, ymin=18.4, ymax=18.42))
agbc <- crop(agb, test_bb)
levelplot(agbc)
adc <- crop(a_date, test_bb)
levelplot(adc)
smc <- crop(sm_proj, test_bb)
levelplot(smc)

tm_shape(agbc) + 
  tm_raster(breaks=seq(0, 300, 10), palette=palette(hcl.colors(8, "viridis"))) +
  tm_shape(smc) + 
  tm_raster()








