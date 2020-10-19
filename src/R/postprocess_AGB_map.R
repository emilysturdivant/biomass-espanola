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

# Load libraries
# library(stars)
library(geobgu) # not available for R 4.0
library(broom)
library(gdalUtils)
library(tmap)
require(graphics)
library(rasterVis)
# library(raster)
library(terra)
library(tidyverse)
library(clipr)

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


get_mask_counts_terra <- function(fn_suff, p3_agb_val, agb_in_fn, dir='results/tifs_by_R'){
  
  # Filenames
  agb_in_fp <- # This file is only Haiti
    file.path(dir, agb_in_fn) 
  agb_hti_fp <- file.path(dir, str_c("agb18_v3_l0_hti",fn_suff,".tif"))
  agb_out_fp <- file.path(dir, str_c("agb18_v3_l1_mask_Ap3WUw25_u20_hti",fn_suff,".tif"))
  
  # Load default mask (ALOS, water, and urban)
  msk_land <- # 1== ALOS 2018 land (i.e. Normal, Layover, and Shadowing)- HAITI
    terra::rast("results/masks/hti18_maskLand.tif") 
  msk_AWUwb <- # NA==ALOS mask, LC17 Water and Urban, and OSM water with 25 m buffer
    terra::rast("results/masks/mask_ALOS_WaterUrban_water25.tif")
  mskinv_T <-  # 1== LC17 tree cover
    terra::rast("results/masks/mask_inv_TreeCover.tif")
  
  # Get cropped AGB and extent for cropping
  if(file.exists(agb_hti_fp)){
    # Load cropped AGB
    agb.ras <- terra::rast(agb_hti_fp)
    
  } else {
    # Load input AGB
    agb.ras <- terra::rast(agb_in_fp)
    
    # Reproject if the resolution is different
    # if(any(res(agb.ras) == res(msk_land))){
    #   t.ras <- terra::rast(file.path(dir, str_c("agb18_v1_l0_hti.tif")))
    #   agb.ras <- agb.ras %>% terra::resample(t.ras, method="ngb")
    # }
    
    # set extent and crop
    bb <- st_intersection(ext(agb.ras) %>% as.vector %>% st_bbox %>% st_as_sfc, 
                          ext(msk_land) %>% as.vector %>% st_bbox %>% st_as_sfc) %>% 
      st_bbox %>% as.vector 
    bbex <- ext(bb[c(1, 3, 2, 4)])
    
    agb.ras <- agb.ras %>% terra::crop(bbex, filename=agb_hti_fp)
  }
  
  # set extent and crop
  msk_land <- msk_land %>% terra::crop(agb.ras)
  msk_AWUwb <- msk_AWUwb %>% terra::crop(agb.ras)
  mskinv_T <- mskinv_T %>% terra::crop(agb.ras)
  
  # Get masks
  # AGB <20 mask
  msk_u20_fp <- str_c("results/R_out/mask_AGB_under20",fn_suff,"_hti.tif")
  if(file.exists(msk_u20_fp)){
    msk_u20 <- terra::rast(msk_u20_fp)
  } else {
    msk_u20 <- agb.ras %>% 
      terra::classify(rbind(c(-9999,20,NA), c(20,9999,1)),
                      filename=msk_u20_fp, overwrite=T, 
                      wopt=list(datatype='INT1U', gdal='NBITS=1', NAflag=0))
  }
  
  msk_AWUwb_u20_fp <- str_c("results/R_out/mask_AWUwb_u20",fn_suff,"_hti.tif")
  if(file.exists(msk_AWUwb_u20_fp)){
    msk_AWUwb_u20 <- terra::rast(msk_AWUwb_u20_fp)
  } else {
    msk_AWUwb_u20 <- msk_AWUwb %>% 
      terra::mask(msk_u20, filename=msk_AWUwb_u20_fp, overwrite=T, 
                  wopt=list(datatype='INT1U', gdal='NBITS=1', NAflag=0))
  }
  
  # g0 >0.3 mask (AGB>316.76)
  msk_p3_fp <- str_c("results/R_out/mask_ALOSoverpt3",fn_suff,"_hti.tif")
  if(file.exists(msk_p3_fp)) {
    msk_p3 <- terra::rast(msk_p3_fp)
  } else {
    msk_p3 <- agb.ras %>% 
      terra::classify(rbind(c(p3_agb_val, 9999, NA), 
                            c(-9999, p3_agb_val, 1)),
                      filename=msk_p3_fp, overwrite=T, 
                      wopt=list(datatype='INT1U', gdal='NBITS=1', NAflag=0))
  }
  
  msk_Ap3WUwb_u20_fp <- str_c("results/R_out/mask_msk_Ap3WUwb_u20",fn_suff,"_hti.tif")
  if(file.exists(msk_Ap3WUwb_u20_fp)){
    msk_Ap3WUwb_u20 <- terra::rast(msk_Ap3WUwb_u20_fp) 
  } else {
    msk_Ap3WUwb_u20 <- msk_AWUwb_u20 %>% 
      terra::mask(msk_p3, filename=msk_Ap3WUwb_u20_fp, overwrite=T, 
                  wopt=list(datatype='INT1U', gdal='NBITS=1', NAflag=0))
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
  wb_fp <- 'results/masks/vector/osm_water_buff25m.shp'
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
    str_c("results/tifs_by_R/", agb_in_fn) 
  agb_hti_fp <- str_c("results/tifs_by_R/agb18_v3_l0_hti",fn_suff,".tif")
  agb_out_fp <- str_c("results/tifs_by_R/agb18_v3_l1_mask_Ap3WUw25_u20_hti",fn_suff,".tif")
  
  # Load default mask (ALOS, water, and urban)
  msk_land <- # 1== ALOS 2018 land (i.e. Normal, Layover, and Shadowing)- HAITI
    raster("results/masks/hti18_maskLand_clip2border.tif") 
  msk_AWUwb <- # NA==ALOS mask, LC17 Water and Urban, and OSM water with 25 m buffer
    raster("results/masks/mask_ALOS_WaterUrban_water25.tif")
  mskinv_T <-  # 1== LC17 tree cover
    raster("results/masks/mask_inv_TreeCover.tif")
  
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
      t.ras <- raster(str_c("results/tifs_by_R/agb18_v1_l0_hti.tif"))
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
  # msk_u20_fp <- str_c("results/R_out/mask_AGB_under20_raster",fn_suff,"_hti.rds")
  msk_u20_fp <- str_c("results/R_out/mask_AGB_under20",fn_suff,"_hti.tif")
  if(file.exists(msk_u20_fp)){
    # msk_u20 <- readRDS(msk_u20_fp)
    msk_u20 <- raster(msk_u20_fp)
  } else {
    msk_u20 <- agb.ras
    msk_u20[msk_u20 < 20] <- NA
    msk_u20[msk_u20 >= 20] <- 1
    msk_u20 %>% writeRaster(msk_u20_fp)
    # msk_u20 %>% saveRDS(str_c("results/R_out/mask_AGB_under20_raster",fn_suff,"_hti.rds"))
  }
  
  # msk_AWUwb_u20_fp <- str_c("results/R_out/mask_AWUwb_u20_raster",fn_suff,"_hti.rds")
  msk_AWUwb_u20_fp <- str_c("results/R_out/mask_AWUwb_u20",fn_suff,"_hti.tif")
  if(file.exists(msk_AWUwb_u20_fp)){
    # msk_AWUwb_u20 <- readRDS(msk_AWUwb_u20_fp)
    msk_AWUwb_u20 <- raster(msk_AWUwb_u20_fp)
  } else {
    msk_AWUwb_u20 <- msk_AWUwb * msk_u20
    # msk_AWUwb_u20 %>% saveRDS(msk_AWUwb_u20_fp)
    msk_AWUwb_u20 %>% writeRaster(msk_AWUwb_u20_fp)
  }
  
  # g0 >0.3 mask (AGB>316.76)
  # msk_p3_fp <- str_c("results/R_out/mask_ALOSoverpt3_raster",fn_suff,"_hti.rds")
  msk_p3_fp <- str_c("results/R_out/mask_ALOSoverpt3",fn_suff,"_hti.tif")
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
  
  # msk_Ap3WUwb_u20_fp <- str_c("results/R_out/mask_msk_Ap3WUwb_u20_raster",fn_suff,"_hti.rds")
  msk_Ap3WUwb_u20_fp <- str_c("results/R_out/mask_msk_Ap3WUwb_u20",fn_suff,"_hti.tif")
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


# Get mask counts --------------------------------------------------------------------
# Unfiltered AGB (v1) 
fn_suff <- '' # filename suffix (from regression_AGB-g0.R)
p3_agb_val = 310.02
agb_in_fn <- 'agb18_v1_l0_htiR.tif'

cts_df_v1 <- get_mask_counts(fn_suff, p3_agb_val, agb_in_fn)

# Aggregate 50 m (v2)
fn_suff <- '_agg50m' # filename suffix (from regression_AGB-g0.R)
p3_agb_val = 304.71
agb_in_fn <- 'agb18_v2_l0_hti_agg50m.tif'

cts_df_v2 <- get_mask_counts(fn_suff, p3_agb_val, agb_in_fn)

# SAGA filter (v3) 
fn_suff <- '_qLee' # filename suffix (from regression_AGB-g0.R)
p3_agb_val = 316.76
agb_in_fn <- 'agb18_v3_l0_qLee.tif'

cts_df_v3 <- get_mask_counts(fn_suff, p3_agb_val, agb_in_fn)



# SAGA filter (v3) after filling gaps by LC
fn_suff <- '_qLee_filledLCpatches' # filename suffix (from regression_AGB-g0.R)
p3_agb_val = 316.76
agb_in_fn <- 'agb18_v3_l0_qLee.tif'
agb_in_fn <- 'agb18_v3_l1_Ap3WUw25u20_hti_filled_LCpatches.tif'

(cts_df_v3l3 <- get_mask_counts_terra(fn_suff, p3_agb_val, agb_in_fn))
cts_df_v3l3$V1 %>% write_clip

# Get missing pixel counts for AGB with filled gaps
msk_land <- # 1== ALOS 2018 land (i.e. Normal, Layover, and Shadowing)- HAITI
  terra::rast("results/masks/hti18_maskLand.tif") 
mskinv_T <-  # 1== LC17 tree cover
  terra::rast("results/masks/mask_inv_TreeCover.tif")
wb_fp <- 'results/masks/vector/osm_water_buff25m.shp'
mskinv_T_wb <- mskinv_T %>% mask(vect(wb_fp))

dir <- 'results/tifs_by_R'
agb_in_fn <- 'agb18_v3_l1_Ap3WUw25u20_hti_filled_LCpatches.tif'
agb.ras <- file.path(dir, agb_in_fn) %>% rast

bb <- st_intersection(ext(agb.ras) %>% as.vector %>% st_bbox %>% st_as_sfc, 
                      ext(msk_land) %>% as.vector %>% st_bbox %>% st_as_sfc) %>% 
  st_bbox %>% as.vector 
bbex <- ext(bb[c(1, 3, 2, 4)])

agb.ras <- agb.ras %>% terra::crop(bbex) %>% mask(msk_land)
msk_land <- msk_land %>% terra::crop(bbex)
mskinv_T <- mskinv_T %>% terra::crop(bbex)
plot(agb.ras)
plot(msk_land)

# Convert to DF
df <- tibble(land = as.vector(msk_land), 
             agb = as.vector(agb.ras) %>% round(4), 
             zone = as.vector(mskinv_T),
             water_in_tc = as.vector(mskinv_T_wb)) %>% 
  filter(!is.na(land))
total <- nrow(df)

tc_df <- df %>% filter(zone == 1)
total_TC <- tc_df %>% nrow

# Tree cover class
(ct_TC_u20 <- tc_df %>% filter(agb < 20) %>% nrow)
ct_TC_u132 <- tc_df %>% filter(agb > 132) %>% nrow
ct_TC_u132 / total_TC
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

# Get pct of TC class that intersects water buffer

ct_TC_pixels_wb <- global(!is.na(mskinv_T_wb), 'sum')
TC_wb <- ct_TC_pixels_wb / ct_TC_pixels


# Visualize
ggplot(df, aes(x=agb)) +
  geom_histogram()

ggplot(df, aes(x=zone, y=agb)) +
  geom_bar(stat="identity")








# Older code -------------------------------------------------------------------
# Barplot
(bp <- ggplot(df, aes(x="", y=value, fill=group))+
    geom_bar(width = 1, stat = "identity")+ 
    scale_fill_manual(values=c("#E69F00", "#999999", "#999000")) +
    theme_minimal())
(pie <- bp + coord_polar("y", start=0))


# Apply masks to AGB --------------------------------------------------------------------------------------
agb.ras <- read_stars("results/tifs_by_R/agb18_v1_l0.tif")
msk_Ap3_WU_wb_u20 <- 
  readRDS('results/R_out/mask_ALOS_pt3_WaterUrban_water25_AGBu20_stars.rds')

agb_mskd <- agb.ras * msk_Ap3_WU_wb
agb_mskd %>% as("Raster") %>% 
  writeRaster('results/tifs_by_R/agb18_v1_l1_mask_Ap3WUw25.tif')

agb_mskd <- agb.ras * msk_Ap3_WU_wb_u20
agb_mskd %>% as("Raster") %>% 
  writeRaster('results/tifs_by_R/agb18_v1_l1_mask_Ap3WUw25_u20.tif')

agb_mskd2 <- agb_mskd * msk_o310
agb_mskd2 %>% as("Raster") %>% 
  writeRaster('results/tifs_by_R/agb18_v1_l1_mask_Ap3WUw25_u20o310.tif')
agb_mskd2 <- read_stars('results/tifs_by_R/agb18_v1_l1_mask_Ap3WUw25_u20o310.tif')

# Apply value cap to AGB ----#######################################################
agb_sat <- agb_mskd
agb_sat <- read_stars('results/tifs_by_R/agb18_v1_l1_mask_Ap3WUw25_u20.tif')
saturation_pt <- 132

# Apply saturation point to AGB and save
agb_sat[agb_sat > saturation_pt] <- saturation_pt
agb_sat %>% saveRDS('results/R_out/agb18_v1_l2_mask_Ap3WUw25_u20_cap132.rds')
agb_sat %>% as("Raster") %>% 
  writeRaster('results/tifs_by_R/agb18_v1_l2_mask_Ap3WUw25_u20_cap132.tif')

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

# Look at AGB distributions ---- ###################################################
agb.br <- brick(raster("results/agb/agb_2018_v6_mask2share.tif"),
                raster("results/agb/agb_2018_v6CI_2share.tif"))
names(agb.br) <- c('AGB', 'CI')
saveRDS(agb.br, file = "results/R_out/AGB_95ci.rds")

get_brick_stats <- function(lc.br){
  # Get selection of percentiles for each LC
  agb.qs <- data.frame(row.names=c('0%', '1%','2%','10%', '25%', '50%', '75%', '90%','98%', '99%','100%'))
  for (i in seq(1,nlayers(lc.br))){
    lc <- names(lc.br[[i]])
    agb.qs[[lc]] <- quantile(lc.br[[i]], probs=c(0, 0.01, 0.02, 0.1, 0.25, 0.5, 0.75, 0.9,0.98, 0.99, 1))
  }
  # Get means, SDs, and Skews
  stats <- list(mean=cellStats(lc.br, stat='mean', na.rm=TRUE),
                sd=cellStats(lc.br, stat='sd', na.rm=TRUE),
                skew=cellStats(lc.br, stat='skew', na.rm=TRUE)) %>%
    bind_cols()%>%
    t()
  colnames(stats) <- names(lc.br)
  agb.stats <- rbind(agb.qs, stats)
}
agb.stats <- get_brick_stats(agb.br)

# Graphing ---- ###################################################################
# Sample the distribution of values in the raster
agb.samp <- agb %>% 
  sampleRandom(100000, na.rm=TRUE) %>% 
  as.data.frame() %>% 
  rename(AGB='.')

# Plot density
p <- ggplot(agb.samp, aes(x=AGB)) + 
  geom_histogram(aes(y=..density..), fill="#69b3a2", color="#e9ecef", alpha=0.7, binwidth = 2.5)+
  geom_density(alpha=.2) + 
  scale_x_continuous(name = expression(paste("Aboveground biomass (Mg ha"^"-1", ")")),
                     breaks = seq(0, 150, 25),
                     limits=c(-5, 150)) +
  geom_vline(aes(xintercept=mean(AGB)), color="black", linetype="dashed", size=.5)+
  geom_vline(aes(xintercept=median(AGB)),
             color="black", linetype="dashed", size=.5)+
  theme_minimal()
p

mean(agb.samp$AGB)
median(agb.samp$AGB)
