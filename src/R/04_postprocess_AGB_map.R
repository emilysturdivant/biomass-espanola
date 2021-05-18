# *************************************************************************************************
# Script to:
#     * Perform post-process on AGB map 
#       - count pixels meeting various criteria
# Proceeds:
#     * regression_AGB-g0.R - creates AGB map
#     * calculate_AGB.R - calculates AGB by plot from field data
#     * process_ALOS_tiles.R
# Requires:
#     * AGB map (agbr)
#     * pre-created masks
#
# *************************************************************************************************

# Load libraries ----
# library(geobgu) 
# library(broom)
# library(gdalUtils)
# library(tmap)
# require(graphics)
# library(rasterVis)
library("terra")
library("tidyverse")
# library(clipr)

# Set variables ----
year <- '2019'
code <- 'sl_HV'
suffix <- g0_variant <- 'lee11s10'
p3_agb_val = 310

# raw_dir <- file.path('data/raw/ALOS', year)
modeling_dir <- 'data/modeling'

tidy_dir <- 'data/tidy'
palsar_dir <- file.path(tidy_dir, str_c('palsar_', year))
masks_dir <- file.path(palsar_dir, 'masks')
landmask_fp <- file.path(masks_dir, 'hti_land_palsar.tif')

# Functions ------------------------------------------------------------------------------------------
get_counts <- function(valuesraster, msk_zone, msk_0, threshold=132){
  # For given area of interest (msk_zone!=NA), 
  # count non-masked values (msk_0!=NA) and values>threshold after masking 
  #  - msk_zone - defines zone of interest ("total")
  #  - msk_0 - which pixels to count inside msk_zone ("masked")
  #  - threshold - count values > threshold after masking
  # 1. Convert rasters to DF
  df <- 
    tibble(zone=as.vector(msk_zone[[1]]),
           agb=as.vector(valuesraster[[1]]) %>% round(4), 
           mask=as.vector(msk_0[[1]])) %>% 
    filter(!is.na(zone)) 
  total = length(df[[1]])
  # Filter by mask
  df1 <- df %>% filter(!is.na(mask))
  total1 = length(df1[[1]])
  # Filter by (upper) threshold value
  df2 <- df1 %>% filter(!agb > threshold)
  total2 = length(df2[[1]])
  # Tally it up
  cts_df <- tibble(
    masked = total-total1,
    over_thresh = total1-total2,
    valid=total2
  ) %>% 
    pivot_longer(everything()) %>% 
    mutate(pct = value/sum(value))
}

filter_by_thresh_and_tally <- function(df1, total, threshold){
  total1 = nrow(df1)
  
  # Filter by (upper) threshold value
  df2 <- df1 %>% filter(!agb > threshold)
  total2 = nrow(df2)
  
  # Tally it up
  cts_df <- tibble(
    masked = total-total1,
    over_thresh = total1-total2,
    valid = total2
  ) %>% 
    pivot_longer(everything()) %>% 
    mutate(pct = value/sum(value))
}

get_counts_raster <- function(valuesraster, msk_zone, msk_0, threshold=132){
  # For given area of interest (msk_zone!=NA), 
  # count non-masked values (msk_0!=NA) and values>threshold after masking 
  #  - msk_zone - defines zone of interest ("total")
  #  - msk_0 - which pixels to count inside msk_zone ("masked")
  #  - threshold - count values > threshold after masking
  # 1. Convert rasters to DF
  df <- 
    tibble(zone = as.vector(msk_zone),
           agb = as.vector(valuesraster) %>% round(4), 
           mask = as.vector(msk_0)) %>% 
    filter(!is.na(zone)) 
  total = nrow(df)
  
  # Filter by mask
  df1 <- df %>% filter(!is.na(mask))
  total1 = nrow(df1)
  
  # Filter by (upper) threshold value
  cts_df <- filter_by_thresh_and_tally(df1, total, threshold)
}

get_mask_counts_terra <- function(fn_suff, p3_agb_val, agb_fp, dir=file.path(results_dir, 'tifs_by_R')){
  
  # Load AGB
  agb.ras <- terra::rast(agb_fp)
  
  # Load default mask (ALOS, water, and urban)
  msk_land <- # 1== ALOS 2018 land (i.e. Normal, Layover, and Shadowing)- HAITI
    terra::rast(file.path(results_dir, "masks/hti18_maskLand.tif") )
  msk_AWUwb <- # NA==ALOS mask, LC17 Water and Urban, and OSM water with 25 m buffer
    terra::rast(file.path(results_dir, "masks/mask_ALOS_WaterUrban_water25.tif"))
  mskinv_T <-  # 1== LC17 tree cover
    terra::rast(file.path(results_dir, "masks/mask_inv_TreeCover.tif"))
  
  # set extent and crop
  msk_land <- msk_land %>% terra::crop(agb.ras)
  msk_AWUwb <- msk_AWUwb %>% terra::crop(agb.ras)
  mskinv_T <- mskinv_T %>% terra::crop(agb.ras)
  
  # Get masks
  # AGB <20 mask
  msk_u20_fp <- file.path(results_dir, str_c("R_out/mask_AGB_under20",fn_suff,"_hti.tif"))
  if(file.exists(msk_u20_fp)){
    msk_u20 <- terra::rast(msk_u20_fp)
  } else {
    msk_u20 <- agb.ras %>% 
      terra::classify(rbind(c(-9999,20,NA), c(20,9999,1)),
                      filename=msk_u20_fp, overwrite=T, 
                      wopt=list(datatype='INT1U', gdal='COMPRESS=LZW', NAflag=0))
  }
  
  msk_AWUwb_u20_fp <- file.path(results_dir, str_c("R_out/mask_AWUwb_u20",fn_suff,"_hti.tif"))
  if(file.exists(msk_AWUwb_u20_fp)){
    msk_AWUwb_u20 <- terra::rast(msk_AWUwb_u20_fp)
  } else {
    msk_AWUwb_u20 <- msk_AWUwb %>% 
      terra::mask(msk_u20, filename=msk_AWUwb_u20_fp, overwrite=T, 
                  wopt=list(datatype='INT1U', gdal='COMPRESS=LZW', NAflag=0))
  }
  
  # g0 >0.3 mask (AGB>316.76)
  msk_p3_fp <- file.path(results_dir, str_c("R_out/mask_ALOSoverpt3",fn_suff,"_hti.tif"))
  if(file.exists(msk_p3_fp)) {
    msk_p3 <- terra::rast(msk_p3_fp)
  } else {
    msk_p3 <- agb.ras %>% 
      terra::classify(rbind(c(p3_agb_val, 9999, NA), 
                            c(-9999, p3_agb_val, 1)),
                      filename=msk_p3_fp, overwrite=T, 
                      wopt=list(datatype='INT1U', gdal='COMPRESS=LZW', NAflag=0))
  }
  
  msk_Ap3WUwb_u20_fp <- file.path(results_dir, str_c("R_out/mask_msk_Ap3WUwb_u20",fn_suff,"_hti.tif"))
  if(file.exists(msk_Ap3WUwb_u20_fp)){
    msk_Ap3WUwb_u20 <- terra::rast(msk_Ap3WUwb_u20_fp) 
  } else {
    msk_Ap3WUwb_u20 <- msk_AWUwb_u20 %>% 
      terra::mask(msk_p3, filename=msk_Ap3WUwb_u20_fp, overwrite=T, 
                  wopt=list(datatype='INT1U', gdal='COMPRESS=LZW', NAflag=0))
  }
  
  # Apply non-negotiable mask (msk_AWUwb_u20) to AGB from SAGA filter
  if(file.exists(agb_out_fp)){
    # read
    agb.r <- terra::rast(agb_out_fp)
  } else {
    agb.r <- agb.ras %>% 
      terra::mask(msk_Ap3WUwb_u20, filename=agb_out_fp, overwrite=T, 
                  wopt=list(datatype='FLT4S', gdal='COMPRESS=LZW'))
  }

  # Get mean and count
  mn <- global(agb.r, 'mean', na.rm=T)
  ct_valid_pixels <- global(!is.na(agb.r), 'sum')
  est_area_ha <- (25*25*0.0001) * ct_valid_pixels
  est_total_AGB <- mn * est_area_ha
  
  # Count NAs as percent of all land
  Ap3WUwb_u20 <- get_counts_raster(agb.ras, msk_land, msk_Ap3WUwb_u20, 310) %>% 
    filter(name=='masked') %>% select('pct')
  
  # 1. Convert rasters to DF
  df <- tibble(zone = as.vector(msk_land),
               agb = as.vector(agb.ras) %>% round(4), 
               mask = as.vector(msk_AWUwb_u20)) %>% 
    filter(!is.na(zone)) 
  total <- nrow(df)
  
  # Filter by mask
  df1 <- df %>% filter(!is.na(mask))
  
  df <- filter_by_thresh_and_tally(df1, total, 310)
  AWUwb_u20 <- df %>% filter(name=='masked') %>% select('pct')
  over310 <- df %>% filter(name=='over_thresh') %>% select('pct')
  over300 <- filter_by_thresh_and_tally(df1, total, 300) %>% filter(name=='over_thresh') %>% select('pct')
  over250 <- filter_by_thresh_and_tally(df1, total, 250) %>% filter(name=='over_thresh') %>% select('pct')
  over200 <- filter_by_thresh_and_tally(df1, total, 200) %>% filter(name=='over_thresh') %>% select('pct')
  over132 <- filter_by_thresh_and_tally(df1, total, 132) %>% filter(name=='over_thresh') %>% select('pct')
  
  # get Tree Cover
  df <- get_counts_raster(agb.ras, mskinv_T, msk_u20, 132)
  (TCu20 <- df %>% filter(name=='masked') %>% select('pct'))
  (TCu20o132 <- df %>% filter(name=='over_thresh') %>% select('pct'))
  df <- get_counts_raster(agb.ras, mskinv_T, msk_AWUwb_u20, 132)
  (TC_AWUwb_u20 <- df %>% filter(name=='masked') %>% select('pct'))
  (TCover132 <- df %>% filter(name=='over_thresh') %>% select('pct'))
  df <- get_counts_raster(agb.ras, mskinv_T, mskinv_T, 132)
  (TCo132 <- df %>% filter(name=='over_thresh') %>% select('pct'))
  df <- get_counts_raster(agb.ras, mskinv_T, msk_AWUwb, 132)
  (TC_AWUwb <- df %>% filter(name=='masked') %>% select('pct'))
  
  # Get pct of TC class that intersects water buffer
  wb_fp <- file.path(results_dir, 'masks/vector/osm_water_buff25m.shp')
  wb <- vect(wb_fp)
  ct_TC_pixels <- global(!is.na(mskinv_T), 'sum')
  mskinv_T_wb <- mskinv_T %>% mask(wb)
  ct_TC_pixels_wb <- global(!is.na(mskinv_T_wb), 'sum')
  TC_wb <- ct_TC_pixels_wb / ct_TC_pixels
  
  # Assemble report 
  df_out <- 
    tibble(
      Ap3WUwb_u20 = Ap3WUwb_u20 %>% deframe, 
      AWUwb_u20 = AWUwb_u20 %>% deframe, 
      over310 = over310 %>% deframe, 
      over300 = over300 %>% deframe, 
      over250 = over250 %>% deframe, 
      over200 = over200 %>% deframe, 
      over132 = over132 %>% deframe, 
      TC_AWUwb = TC_AWUwb %>% deframe, 
      TC_wb = TC_wb %>% deframe,
      TCu20 = TCu20 %>% deframe, 
      TC_AWUwb_u20 = TC_AWUwb_u20 %>% deframe, 
      TCover132 = TCover132 %>% deframe, 
      TCu20o132 = TCu20o132 %>% deframe, 
      TCo132 = TCo132 %>% deframe, 
      avg_AGB = mn %>% deframe, 
      est_total_AGB = est_total_AGB %>% deframe
      ) %>% 
    t %>% 
    as_tibble(rownames = 'name')
}

get_mask_counts <- function(fn_suff, p3_agb_val, agb_in_fn){
  
  # Filenames
  agb_in_fp <- # This file is only Haiti
    file.path(results_dir, str_c("tifs_by_R/", agb_in_fn) )
  agb_hti_fp <- file.path(results_dir, str_c("tifs_by_R/agb18_v3_l0_hti",fn_suff,".tif"))
  agb_out_fp <- file.path(results_dir, str_c("tifs_by_R/agb18_v3_l1_mask_Ap3WUw25_u20_hti",fn_suff,".tif"))
  
  # Load default mask (ALOS, water, and urban)
  msk_land <- # 1== ALOS 2018 land (i.e. Normal, Layover, and Shadowing)- HAITI
    raster(file.path(results_dir, "masks/hti18_maskLand_clip2border.tif") )
  msk_AWUwb <- # NA==ALOS mask, LC17 Water and Urban, and OSM water with 25 m buffer
    raster(file.path(results_dir, "masks/mask_ALOS_WaterUrban_water25.tif"))
  mskinv_T <-  # 1== LC17 tree cover
    raster(file.path(results_dir, "masks/mask_inv_TreeCover.tif"))
  
  # Get cropped AGB and extent for cropping
  if(file.exists(agb_hti_fp)){
    # Load cropped AGB
    agb.ras <- raster(agb_hti_fp)
    
    # Get extent
    bb <- agb.ras %>% st_bbox() %>% as.vector()
    bbex <- extent(bb[c(1,3,2,4)])
  } else {
    # Load input AGB
    agb.ras <- raster(agb_in_fp)

    # Reproject if the resolution is different
    if(res(agb.ras) != res(msk_land)){
      t.ras <- raster(file.path(results_dir, str_c("tifs_by_R/agb18_v1_l0_hti.tif")))
      agb.ras <- agb.ras %>% resample(t.ras, method="ngb")
    }

    # set extent and crop
    bb <- st_intersection(st_bbox(agb.ras) %>% st_as_sfc(), 
                          st_bbox(msk_land) %>% st_as_sfc()) %>% 
      st_bbox() %>% as.vector()
    bbex <- extent(bb[c(1,3,2,4)])
    agb.ras <- agb.ras %>% crop(bbex)
    # save
    agb.ras %>% writeRaster(agb_hti_fp)
  }
  
  # set extent and crop
  msk_land <- msk_land %>% crop(bbex)
  msk_AWUwb <- msk_AWUwb %>% crop(bbex)
  mskinv_T <- mskinv_T %>% crop(bbex)
  
  # Get masks
  # AGB <20 mask
  # msk_u20_fp <- file.path(results_dir, str_c("R_out/mask_AGB_under20_raster",fn_suff,"_hti.rds")
  msk_u20_fp <- file.path(results_dir, str_c("R_out/mask_AGB_under20",fn_suff,"_hti.tif"))
  if(file.exists(msk_u20_fp)){
    # msk_u20 <- readRDS(msk_u20_fp)
    msk_u20 <- raster(msk_u20_fp)
  } else {
    msk_u20 <- agb.ras
    msk_u20[msk_u20 < 20] <- NA
    msk_u20[msk_u20 >= 20] <- 1
    msk_u20 %>% writeRaster(msk_u20_fp)
    # msk_u20 %>% saveRDS(file.path(results_dir, str_c("R_out/mask_AGB_under20_raster",fn_suff,"_hti.rds"))
  }
  
  # msk_AWUwb_u20_fp <- file.path(results_dir, str_c("R_out/mask_AWUwb_u20_raster",fn_suff,"_hti.rds")
  msk_AWUwb_u20_fp <- file.path(results_dir, str_c("R_out/mask_AWUwb_u20",fn_suff,"_hti.tif"))
  if(file.exists(msk_AWUwb_u20_fp)){
    # msk_AWUwb_u20 <- readRDS(msk_AWUwb_u20_fp)
    msk_AWUwb_u20 <- raster(msk_AWUwb_u20_fp)
  } else {
    msk_AWUwb_u20 <- msk_AWUwb * msk_u20
    # msk_AWUwb_u20 %>% saveRDS(msk_AWUwb_u20_fp)
    msk_AWUwb_u20 %>% writeRaster(msk_AWUwb_u20_fp)
  }
  
  # g0 >0.3 mask (AGB>316.76)
  # msk_p3_fp <- file.path(results_dir, str_c("R_out/mask_ALOSoverpt3_raster",fn_suff,"_hti.rds")
  msk_p3_fp <- file.path(results_dir, str_c("R_out/mask_ALOSoverpt3",fn_suff,"_hti.tif"))
  if(file.exists(msk_p3_fp)) {
    # msk_p3 <- readRDS(msk_p3_fp)
    msk_p3 <- raster(msk_p3_fp)
  } else {
    msk_p3 <- agb.ras
    msk_p3[msk_p3 > p3_agb_val] <- NA
    msk_p3[msk_p3 <= p3_agb_val] <- 1
    # msk_p3 %>% saveRDS(msk_p3_fp)
    msk_p3 %>% writeRaster(msk_p3_fp)
  }
  
  # msk_Ap3WUwb_u20_fp <- file.path(results_dir, str_c("R_out/mask_msk_Ap3WUwb_u20_raster",fn_suff,"_hti.rds")
  msk_Ap3WUwb_u20_fp <- file.path(results_dir, str_c("R_out/mask_msk_Ap3WUwb_u20",fn_suff,"_hti.tif"))
  if(file.exists(msk_Ap3WUwb_u20_fp)){
    # msk_Ap3WUwb_u20 <- readRDS(msk_Ap3WUwb_u20_fp)
    msk_Ap3WUwb_u20 <- raster(msk_Ap3WUwb_u20_fp)
  } else {
    msk_Ap3WUwb_u20 <-  msk_AWUwb_u20 * msk_p3
    # msk_Ap3WUwb_u20 %>% saveRDS(msk_Ap3WUwb_u20_fp)
    msk_Ap3WUwb_u20 %>% writeRaster(msk_Ap3WUwb_u20_fp)
  }
  
  # Apply non-negotiable mask (msk_AWUwb_u20) to AGB from SAGA filter
  if(file.exists(agb_out_fp)){
    # read
    agb.r <- raster(agb_out_fp)
  } else {
    agb.r <- agb.ras * msk_Ap3WUwb_u20
    agb.r %>% 
      writeRaster(agb_out_fp, overwrite=T)
  }
  
  # Get mean and count
  mn <- cellStats(agb.r, 'mean')
  ct_valid_pixels <- cellStats(!is.na(agb.r), 'sum')
  est_area_ha <- (25*25*0.0001) * ct_valid_pixels
  est_total_AGB <- mn * est_area_ha
  
  # Count NAs as percent of all land
  AWUwb <- get_counts_raster(agb.ras, msk_land, msk_Ap3WUwb_u20, 310) %>% filter(name=='masked') %>% select('pct')
  
  # 1. Convert rasters to DF
  df <- tibble(zone = values(msk_land),
               agb = values(agb.ras) %>% round(4), 
               mask = values(msk_AWUwb_u20)) %>% 
    filter(!is.na(zone)) 
  total = nrow(df)
  
  # Filter by mask
  df1 <- df %>% filter(!is.na(mask))
  
  df <- filter_by_thresh_and_tally(df1, total, 310)
  AWUwb_u20 <- df %>% filter(name=='masked') %>% select('pct')
  over310 <- df %>% filter(name=='over_thresh') %>% select('pct')
  over300 <- filter_by_thresh_and_tally(df1, total, 300) %>% filter(name=='over_thresh') %>% select('pct')
  over250 <- filter_by_thresh_and_tally(df1, total, 250) %>% filter(name=='over_thresh') %>% select('pct')
  over200 <- filter_by_thresh_and_tally(df1, total, 200) %>% filter(name=='over_thresh') %>% select('pct')
  over132 <- filter_by_thresh_and_tally(df1, total, 132) %>% filter(name=='over_thresh') %>% select('pct')
  
  # get Tree Cover
  df <- get_counts_raster(agb.ras, mskinv_T, msk_u20, 132)
  (TCu20 <- df %>% filter(name=='masked') %>% select('pct'))
  (TCo132 <- df %>% filter(name=='over_thresh') %>% select('pct'))
  df <- get_counts_raster(agb.ras, mskinv_T, msk_AWUwb_u20, 132)
  (TC_AWUwb_u20 <- df %>% filter(name=='masked') %>% select('pct'))
  (TCover132 <- df %>% filter(name=='over_thresh') %>% select('pct'))
  df <- get_counts_raster(agb.ras, mskinv_T, mskinv_T, 132)
  (TCo132 <- df %>% filter(name=='over_thresh') %>% select('pct'))
  df <- get_counts_raster(agb.ras, mskinv_T, msk_AWUwb, 132)
  (TC_AWUwb <- df %>% filter(name=='masked') %>% select('pct'))
  
  tibble(AWUwb=AWUwb, AWUwb_u20=AWUwb_u20, over310=over310, over300=over300, 
         over250=over250, over200=over200, over132=over132, 
         TCu20=TCu20, TCu20o132=TCo132, TC_AWUwb_u20=TC_AWUwb_u20, 
         TCover132=TCover132, TCo132=TCo132, TC_AWUwb=TC_AWUwb, 
         avg_AGB=mn, est_total_AGB=est_total_AGB) %>% 
    pivot_longer(cols=everything())
}

crop_to_intersecting_extents <- function(r1, r2, return_r1=T, return_r2=T, return_bb=F) {
  # Get intersection of the bounding boxes of the two rasters
  bb <- st_intersection(terra::ext(r1) %>% as.vector %>% st_bbox %>% st_as_sfc, 
                        terra::ext(r2) %>% as.vector %>% st_bbox %>% st_as_sfc) %>% 
    st_bbox %>% as.vector 
  
  # Convert to terra extent object
  bbex <- terra::ext(bb[c(1, 3, 2, 4)])
  
  if(return_bb){
    if(!return_r1 & !return_r2){
      return(bbex)
    } else if(return_r1 & !return_r2) {
      # Crop each raster
      r1 <- r1 %>% terra::crop(bbex)
      return(list(r1=r1, bb=bbex))
    } else if(!return_r1 & return_r2) {
      # Crop each raster
      r2 <- r2 %>% terra::crop(bbex)
      return(list(r2=r2, bb=bbex))
    } else if(return_r1 & return_r2) {
      # Crop each raster
      r1 <- r1 %>% terra::crop(bbex)
      r2 <- r2 %>% terra::crop(bbex)
      return(list(r1=r1, r2=r2, bb=bbex))
    }
  } else {
    if(!return_r1 & !return_r2){
      return()
    } else if(return_r1 & !return_r2) {
      # Crop each raster
      r1 <- r1 %>% terra::crop(bbex)
      return(r1)
    } else if(!return_r1 & return_r2) {
      # Crop each raster
      r2 <- r2 %>% terra::crop(bbex)
      return(r2)
    } else if(return_r1 & return_r2) {
      # Crop each raster
      r1 <- r1 %>% terra::crop(bbex)
      r2 <- r2 %>% terra::crop(bbex)
      return(list(r1=r1, r2=r2))
    }
  }
}


# Get mask counts --------------------------------------------------------------
suffix <- if(suffix == 'simple') '' else str_c('_', suffix)

mod_dir <- file.path(modeling_dir, g0_variant)
agb_fp <- file.path(mod_dir, str_c("agb_l0", suffix, ".tif"))

msk_WUwb_fp <- file.path(masks_dir, "mask_WaterUrban_water25.tif")

# ~ function to mask AGB ----
mask_agb <- function(agb_fp, p3_agb_val, msk_WUwb_fp, 
                     masks = c('alos_normal', 'WU', 'wb', 'p3', 'u20'))
  
# Load AGB
agb_ras <- terra::rast(agb_fp)

# AGB <20 mask
msk_u20 <- agb_ras %>% terra::classify(rbind(c(-9999, 20, NA),
                                             c(20, 9999, 1)))

# g0 >0.3 mask (AGB>316.76)
msk_p3 <- agb_ras %>% terra::classify(rbind(c(p3_agb_val, 9999, NA), 
                                            c(-9999, p3_agb_val, 1)))

# msk_out_fp <- file.path(dirname(agb_fp), str_c('mask_u20to', p3_agb_val, '.tif'))

# LC17 Water and Urban, and OSM water with 25 m buffer
msk_WUwb <- terra::rast(msk_WUwb_fp) %>% terra::crop(agb_ras)

masks2 <- c(masks, msk_WUwb)
prod_msks <- app(masks2, prod, 
                 filename = agb_out_fp, overwrite = T, 
                 wopt = list(datatype = agb_dtype, gdal = 'COMPRESS=LZW'))

# Check extents
if(ext(agb_ras) != ext(msk_WUwb)){
  print("Extents don't match.")
}

# Combine masks
msk_WUwb_u20 <- msk_WUwb %>% terra::mask(msk_u20)
msk_Ap3WUwb_u20 <- msk_WUwb_u20 %>% terra::mask(msk_p3)

# Apply non-negotiable mask (msk_AWUwb_u20) to AGB from SAGA filter
agb_out_fp <- str_c(tools::file_path_sans_ext(agb_fp), 
                    '_maskWUwb20to', p3_agb_val, '.tif')
agb_dtype <- raster::dataType(raster::raster(agb_fp))
agb_masked <- agb_ras %>% 
  terra::mask(msk_Ap3WUwb_u20, filename = agb_out_fp, overwrite = T, 
              wopt = list(datatype = agb_dtype, gdal = 'COMPRESS=LZW'))

# (cts_df_v3l3 <- get_mask_counts_terra(suffix, p3_agb_val, agb_fp))
# ~ Reworking function get_mask_counts_terra ----

# Load AGB
agb_ras <- terra::rast(agb_fp)

# Load masks and crop
# 1 == Intersection of PALSAR 2019 land and Haiti administrative boundary
msk_land <- terra::rast(landmask_fp) %>% 
  terra::crop(agb_ras)
# LC17 Water and Urban, and OSM water with 25 m buffer
msk_WUwb <- terra::rast(file.path(masks_dir, "mask_WaterUrban_water25.tif")) %>% 
  terra::crop(agb_ras)
# 1== LC17 tree cover
mskinv_T <- terra::rast(file.path(masks_dir, "mask_inv_TreeCover.tif")) %>% 
  terra::crop(agb_ras)

# Check extents
if(ext(agb_ras) != ext(msk_WUwb)){
  print("Extents don't match.")
}

# Get masks
# AGB <20 mask
msk_u20 <- agb_ras %>% terra::classify(rbind(c(-9999,20,NA), c(20,9999,1)))
msk_WUwb_u20 <- msk_WUwb %>% terra::mask(msk_u20)

# g0 >0.3 mask (AGB>316.76)
msk_p3 <- agb_ras %>% terra::classify(rbind(c(p3_agb_val, 9999, NA), c(-9999, p3_agb_val, 1)))
msk_Ap3WUwb_u20 <- msk_WUwb_u20 %>% terra::mask(msk_p3)

# Apply non-negotiable mask (msk_AWUwb_u20) to AGB from SAGA filter
agb_out_fp <- str_c(tools::file_path_sans_ext(agb_fp), 
                    '_maskWUwb20to', p3_agb_val, '.tif')
agb_dtype <- raster::dataType(raster::raster(agb_fp))
agb_masked <- agb_ras %>% 
  terra::mask(msk_Ap3WUwb_u20, filename = agb_out_fp, overwrite = T, 
              wopt = list(datatype = agb_dtype, gdal = 'COMPRESS=LZW'))

# if(file.exists(agb_out_fp)){
#   # read
#   agb_masked <- terra::rast(agb_out_fp)
# } else {
#   agb_masked <- agb_ras %>% 
#     terra::mask(msk_Ap3WUwb_u20, filename=agb_out_fp, overwrite=T, 
#                 wopt=list(datatype='FLT4S', gdal='COMPRESS=LZW'))
# }

# Get stats ----
# Get mean and count
mn <- global(agb_masked, 'mean', na.rm=T)
ct_valid_pixels <- global(!is.na(agb_masked), 'sum')
est_area_ha <- (25*25*0.0001) * ct_valid_pixels
est_total_AGB <- mn * est_area_ha

# Count NAs as percent of all land
Ap3WUwb_u20 <- get_counts_raster(agb_ras, msk_land, msk_Ap3WUwb_u20, 310) %>% 
  filter(name == 'masked') %>% select('pct')

# 1. Convert rasters to DF
df <- tibble(zone = as.vector(msk_land),
             agb = as.vector(agb_ras) %>% round(4), 
             mask = as.vector(msk_WUwb_u20)) %>% 
  filter(!is.na(zone)) 

# Count all pixels in land mask
total <- nrow(df)
ct_land_pix <- global(msk_land, 'sum', na.rm = TRUE)

# Filter by mask
df1 <- df %>% filter(!is.na(mask))

df <- filter_by_thresh_and_tally(df1, total, 310)
AWUwb_u20 <- df %>% filter(name=='masked') %>% select('pct')
over310 <- df %>% filter(name=='over_thresh') %>% select('pct')
over300 <- filter_by_thresh_and_tally(df1, total, 300) %>% filter(name=='over_thresh') %>% select('pct')
over250 <- filter_by_thresh_and_tally(df1, total, 250) %>% filter(name=='over_thresh') %>% select('pct')
over200 <- filter_by_thresh_and_tally(df1, total, 200) %>% filter(name=='over_thresh') %>% select('pct')
over132 <- filter_by_thresh_and_tally(df1, total, 132) %>% filter(name=='over_thresh') %>% select('pct')

# get Tree Cover
df <- get_counts_raster(agb_ras, mskinv_T, msk_u20, 132)
(TCu20 <- df %>% filter(name=='masked') %>% select('pct'))
(TCu20o132 <- df %>% filter(name=='over_thresh') %>% select('pct'))
df <- get_counts_raster(agb_ras, mskinv_T, msk_AWUwb_u20, 132)
(TC_AWUwb_u20 <- df %>% filter(name=='masked') %>% select('pct'))
(TCover132 <- df %>% filter(name=='over_thresh') %>% select('pct'))
df <- get_counts_raster(agb_ras, mskinv_T, mskinv_T, 132)
(TCo132 <- df %>% filter(name=='over_thresh') %>% select('pct'))
df <- get_counts_raster(agb_ras, mskinv_T, msk_WUwb, 132)
(TC_AWUwb <- df %>% filter(name=='masked') %>% select('pct'))

# Get pct of TC class that intersects water buffer
wb_fp <- file.path(results_dir, 'masks/vector/osm_water_buff25m.shp')
wb <- vect(wb_fp)
ct_TC_pixels <- global(!is.na(mskinv_T), 'sum')
mskinv_T_wb <- mskinv_T %>% mask(wb)
ct_TC_pixels_wb <- global(!is.na(mskinv_T_wb), 'sum')
TC_wb <- ct_TC_pixels_wb / ct_TC_pixels

# Assemble report 
df_out <- 
  tibble(
    Ap3WUwb_u20 = Ap3WUwb_u20 %>% deframe, 
    AWUwb_u20 = AWUwb_u20 %>% deframe, 
    over310 = over310 %>% deframe, 
    over300 = over300 %>% deframe, 
    over250 = over250 %>% deframe, 
    over200 = over200 %>% deframe, 
    over132 = over132 %>% deframe, 
    TC_AWUwb = TC_AWUwb %>% deframe, 
    TC_wb = TC_wb %>% deframe,
    TCu20 = TCu20 %>% deframe, 
    TC_AWUwb_u20 = TC_AWUwb_u20 %>% deframe, 
    TCover132 = TCover132 %>% deframe, 
    TCu20o132 = TCu20o132 %>% deframe, 
    TCo132 = TCo132 %>% deframe, 
    avg_AGB = mn %>% deframe, 
    est_total_AGB = est_total_AGB %>% deframe
  ) %>% 
  t %>% 
  as_tibble(rownames = 'name')
}





































# Get counts for v3l3 outside of function ----
# Filenames
dir <- file.path(results_dir, 'tifs_by_R')
agb_ras_fp <- file.path(dir, agb_fp) 
msk_land_fp <- file.path(results_dir, "masks/hti18_maskLand.tif")
wb_fp <- file.path(results_dir, 'masks/vector/osm_water_buff25m.shp')
mskinv_T_fp <- file.path(results_dir, "masks/mask_inv_TreeCover.tif")

# Crop to intersecting extents of AGB and land 
out <- crop_to_intersecting_extents(
  terra::rast(agb_ras_fp), 
  terra::rast(msk_land_fp), 
  return_r1=T, return_r2=T, return_bb=T)

# Get rasters
bbex <- out$bb
msk_land <- out$r2
agb.ras <- out$r1 %>% 
  terra::mask(msk_land)
mskinv_T <- terra::rast(mskinv_T_fp) %>% # 1== LC17 tree cover
  terra::crop(bbex)
mskinv_T_wb <- mskinv_T %>% 
  terra::mask(terra::vect(wb_fp))

plot(agb.ras)
plot(msk_land)

# Get mean and count
mn <- global(agb.ras, 'mean', na.rm=T)
sum_agb <- global(agb.ras, 'sum', na.rm=T)
ct_valid_pixels <- global(!is.na(agb.ras), 'sum')
est_area_ha <- (25*25*0.0001) * ct_valid_pixels
est_total_AGB <- mn * est_area_ha
ct_land_pixels <- global(!is.na(msk_land), 'sum')
est_area_ha_land <- (25*25*0.0001) * ct_land_pixels

# Convert to DF
df <- tibble(land = as.vector(msk_land), 
             agb = as.vector(agb.ras) %>% round(2), 
             zone = as.vector(mskinv_T),
             water_in_tc = as.vector(mskinv_T_wb)) %>% 
  filter(!is.na(land))



total <- nrow(df)

tc_df <- df %>% filter(zone == 1)
total_TC <- tc_df %>% nrow

# Tree cover class
(ct_TC_u20 <- tc_df %>% filter(agb < 20) %>% nrow)
ct_TC_o132 <- tc_df %>% filter(agb > 132) %>% nrow
ct_TC_o132 / total_TC
ct_TC_o300 <- tc_df %>% filter(agb > 300) %>% nrow
ct_TC_o300 / total_TC
ct_TC_o200 <- tc_df %>% filter(agb > 200) %>% nrow
ct_TC_o200 / total_TC
ct_TC_na <- tc_df %>% filter(is.na(agb)) %>% nrow
ct_TC_na / total_TC
ct_TC_wb <- tc_df %>% filter(water_in_tc == 1) %>% nrow
ct_TC_wb / total_TC


ct_na <- df %>% filter(is.na(agb)) %>% nrow
ct_na / total
ct_o132 <- df %>% filter(agb > 132) %>% nrow
ct_o132 / total
ct_o200 <- df %>% filter(agb > 200) %>% nrow
ct_o200 / total
ct_o250 <- df %>% filter(agb > 250) %>% nrow
ct_o250 / total
ct_o300 <- df %>% filter(agb > 300) %>% nrow
ct_o300 / total
ct_o310 <- df %>% filter(agb > 310) %>% nrow
ct_o310 / total
ct_u20 <- df %>% filter(agb < 20) %>% nrow
ct_u20 / total

ct_o20 <- df %>% filter(!is.na(agb), agb > 20) %>% nrow
ct_o20 / total

# Get pct of TC class that intersects water buffer

ct_TC_pixels_wb <- global(!is.na(mskinv_T_wb), 'sum')
TC_wb <- ct_TC_pixels_wb / ct_TC_pixels


# Visualize
ggplot(df, aes(x=agb)) +
  geom_histogram()

ggplot(df, aes(x=zone, y=agb)) +
  geom_bar(stat="identity")

# Get count for ALOS mask ----
msk_A <- terra::rast(file.path(results_dir, "masks/mask_ALOS18.tif")) %>% 
  terra::crop(bbex) %>% 
  terra::mask(msk_land)

# Convert to DF
df <- tibble(land = as.vector(msk_land), 
             zone = as.vector(msk_A)) %>% 
  filter(!is.na(land))

total <- nrow(df)
total_TC <- df %>% filter(zone == 1) %>% nrow







# Older code -------------------------------------------------------------------
# Apply masks to AGB -----------------------------------------------------------
agb.ras <- terra::rast(file.path(results_dir, "tifs_by_R/agb18_v1_l0.tif"))
msk_Ap3_WU_wb_u20 <- 
  readRDS(file.path(results_dir, 'R_out/mask_ALOS_pt3_WaterUrban_water25_AGBu20_stars.rds'))

agb_mskd <- agb.ras %>% terra::mask(msk_Ap3_WU_wb)
agb_mskd %>% as("Raster") %>% 
  writeRaster(file.path(results_dir, 'tifs_by_R/agb18_v1_l1_mask_Ap3WUw25.tif'))

agb_mskd <- agb.ras %>% terra::mask(msk_Ap3_WU_wb_u20)
agb_mskd %>% as("Raster") %>% 
  writeRaster(file.path(results_dir, 'tifs_by_R/agb18_v1_l1_mask_Ap3WUw25_u20.tif'))

agb_mskd2 <- agb_mskd %>% terra::mask(msk_o310)
agb_mskd2 %>% as("Raster") %>% 
  writeRaster(file.path(results_dir, 'tifs_by_R/agb18_v1_l1_mask_Ap3WUw25_u20o310.tif'))
agb_mskd2 <- read_stars(file.path(results_dir, 'tifs_by_R/agb18_v1_l1_mask_Ap3WUw25_u20o310.tif'))




# ~ change up masking -----
# SAGA filter (v3) 
fn_suff <- '_qLee' # filename suffix (from regression_AGB-g0.R)
p3_agb_val = 316.76
agb_in_fn <- 'agb18_v3_l0_qLee.tif'
dir=file.path(results_dir, 'tifs_by_R')

# Filenames
agb_in_fp <- # This file is only Haiti
  file.path(dir, agb_in_fn) 
agb_hti_fp <- file.path(dir, str_c("agb18_v3_l0_hti",fn_suff,".tif"))
agb_out_fp <- file.path(dir, str_c("agb18_v3_l1_mask_AWUw25_u20_hti",fn_suff,".tif"))

# Load default mask (ALOS, water, and urban)
msk_land <- # 1== ALOS 2018 land (i.e. Normal, Layover, and Shadowing)- HAITI
  terra::rast(file.path(results_dir, "masks/hti18_maskLand.tif") )
msk_AWUwb <- # NA==ALOS mask, LC17 Water and Urban, and OSM water with 25 m buffer
  terra::rast(file.path(results_dir, "masks/mask_ALOS_WaterUrban_water25.tif"))
mskinv_T <-  # 1== LC17 tree cover
  terra::rast(file.path(results_dir, "masks/mask_inv_TreeCover.tif"))

# Get cropped AGB and extent for cropping
agb.ras <- terra::rast(agb_hti_fp)

# set extent and crop
msk_land <- msk_land %>% terra::crop(agb.ras)
msk_AWUwb <- msk_AWUwb %>% terra::crop(agb.ras)
mskinv_T <- mskinv_T %>% terra::crop(agb.ras)

# Get masks
# AGB <20 mask
msk_u20_fp <- file.path(results_dir, str_c("R_out/mask_AGB_under20",fn_suff,"_hti.tif"))
if(file.exists(msk_u20_fp)){
  msk_u20 <- terra::rast(msk_u20_fp)
} else {
  msk_u20 <- agb.ras %>% 
    terra::classify(rbind(c(-9999,20,NA), c(20,9999,1)),
                    filename=msk_u20_fp, overwrite=T, 
                    wopt=list(datatype='INT1U', gdal='COMPRESS=LZW', NAflag=0))
}

msk_AWUwb_u20_fp <- file.path(results_dir, str_c("R_out/mask_AWUwb_u20",fn_suff,"_hti.tif"))
if(file.exists(msk_AWUwb_u20_fp)){
  msk_AWUwb_u20 <- terra::rast(msk_AWUwb_u20_fp)
} else {
  msk_AWUwb_u20 <- msk_AWUwb %>% 
    terra::mask(msk_u20, filename=msk_AWUwb_u20_fp, overwrite=T, 
                wopt=list(datatype='INT1U', gdal='COMPRESS=LZW', NAflag=0))
}

# Apply non-negotiable mask (msk_AWUwb_u20) to AGB from SAGA filter
if(file.exists(agb_out_fp)){
  # read
  agb.r <- terra::rast(agb_out_fp)
} else {
  agb.r <- agb.ras %>% 
    terra::mask(msk_AWUwb_u20, filename=agb_out_fp, overwrite=T, 
                wopt=list(datatype='FLT4S', gdal='COMPRESS=LZW'))
}


# Apply value cap to AGB ----#######################################################
saturation_pt <- 310
dir_out <- file.path(results_dir, 'tifs_by_R')
agb_fp <- file.path(dir_out, 'agb18_v3_l1_mask_AWUw25_u20_hti_qLee.tif')
# agb_fp <- file.path(dir_out, 'agb18_v3_l0_hti_qLee.tif')
agb_cap_out <- file.path(dir_out, str_c('agb18_v3_l2_mAWUw25u20_cap', 
                        saturation_pt, '.tif'))

# Load
agb_sat <- terra::rast(agb_fp)

# Apply saturation point to AGB and save
agb_sat[agb_sat > saturation_pt] <- saturation_pt
agb_sat %>% 
  terra::writeRaster(agb_cap_out, overwrite=T, 
                     wopt=list(datatype='FLT4U', gdal='COMPRESS=LZW'))

# Load
agb_sat <- raster(agb_cap_out)
ct_valid_agb = sum(!is.na(agb_sat[1]))


# PLOT using tmap (interactive) ----############################################
# Load contextual data for mapping
parks_hti <- # OSM water with 25 m buffer
  st_read('data/contextual_data/OSM_free/hti_nature_reserves_osm.shp')

# Make crop bounding box/extent for testing
# extent for raster
test_ext <- extent(-72.7, -72.5, 18.2, 18.35)
# bounding box for stars
test_bb <- st_bbox(c(xmin=-72.34, xmax=-72.17, ymin=18.3, ymax=18.38), crs=4326) %>% 
  st_as_sfc()

# Set tmap to interactive viewing (vs. "plot")
tmap_mode("view")

# View some masks
names(mskinv_o275) <- 'Mask'
tm_shape(agb_mskd2[test_bb]) + tm_raster() +
  tm_shape(mskinv_o275[test_bb]) + tm_raster(col="Mask", palette="red") +
  tm_shape(parks_hti) + tm_borders()

# Layers
lyr_g0 <- tm_shape(crop(g0, test_ext)) + 
  tm_raster(style="order",
            # breaks=seq(0, 0.1, 0.02), 
            palette=palette(hcl.colors(8, "viridis")))
lyr_water <- 
  tm_shape(water_polys) + 
  tm_borders() + 
  tm_fill(col="cyan") 
lyr_waterB <- tm_shape(water_polysb) + tm_borders() 

# Map
(map_g0w <- lyr_g0 + lyr_water + lyr_waterB)




