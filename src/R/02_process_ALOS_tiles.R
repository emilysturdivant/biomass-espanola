# *************************************************************************************************
# Script to:
#     * Process ALOS mosaic tiles 
#       - load mosaic supplemental rasters, merge tiles, 
#       - make masks
# Preceeds:
#     * calculate_AGB.R - calculates AGB by plot from field data
# Requires:
#     * downloaded ALOS mosaic tiles (currently performed with biota in command line)
#     * study area polygons
#
# *************************************************************************************************
# Load libraries ----
library('gdalUtils')
library('terra')
library('sf')
library('whitebox')
library('tidyverse')

# Initialize ----
source('src/R/initialize_vars.R')

# polarization <- 'HV'
# units <- 'nu'
# code <- str_c(polarization, '_', units)
# year <- '2019'
# 
# suffix <- ''
# tidy_dir <- 'data/tidy'
# raw_dir <- file.path('data/raw/ALOS', year)
# g0_dir <- file.path(tidy_dir, str_c('palsar_', year))
# masks_dir <- file.path(g0_dir, 'masks')
# landmask_fp <- file.path(masks_dir, 'hti_land_palsar.tif')
# hti_poly_fp <- "data/tidy/contextual_data/HTI_adm/HTI_adm0_fix.shp"

# Output filepaths

# Bounding boxes
hisp_bb <- sf::st_bbox(c(xmin = -74.48133, ymax = 20.09044, 
                     xmax = -68.32267, ymin = 17.47022))
hti_bb <- sf::st_bbox(c(xmin = -74.48133, ymax = 20.09044, 
                    xmax = -71.61815, ymin = 18.02180))

# Functions ----

merge_tiles <- function(code, suffix, g0_dir, hti_bb, raw_dir) {
  # Set datatype
  dtype <- if(code == 'mask' | code == 'linci') 'Byte' else 'UInt16'
  nbits <- if(code == 'mask' | code == 'linci') '8' else '16'
  
  # List files to merge
  fps <- list.files(path = raw_dir, pattern = str_glue('{code}.*\\.tif'), full.names = TRUE)
  
  # Create mosaic
  gdalUtils::mosaic_rasters(fps, out_fp, 
                            projwin = c(hti_bb$xmin, hti_bb$ymax, 
                                        hti_bb$xmax, hti_bb$ymin),
                            ot = dtype,
                            co = list('COMPRESS=LZW', str_c('NBITS=', nbits)))
  
  return(out_fp)
}

mask_to_poly <- function(code, suffix, g0_dir, landmask_fp) { 
  # Set output filename
  out_fp <- file.path(g0_dir, str_glue('{code}{suffix}.tif'))
  
  # Mask
  dtype <- ifelse((code == 'mask' | code == 'linci'), 'INT1U', 'INT2U')
  lnd_msk <- terra::rast(landmask_fp)
  r <- terra::rast(out_fp)
  r %>% terra::mask(lnd_msk, 
                    filename = out_fp,
                    overwrite = TRUE,
                    wopt = list(datatype = dtype, gdal = 'COMPRESS=LZW'))
}

mask_with_options <- function(in_fp, masks_dir, code = 'HV_nu', 
                    masks = c('L', 'WU', 'wb'), overwrite = TRUE,
                    out_mask_dir) {
  
  # Get output filenames
  masks_str <- masks %>% str_c(collapse = '')
  masked_fp <- file.path(dirname(in_fp), 
                            str_c(code, '_mask', masks_str, '.tif'))
  msk_all_fp <- file.path(out_mask_dir, str_c('mask', masks_str, '.tif'))
  
  # Stop if file already exists
  if(file.exists(masked_fp) & !overwrite) {
    return(terra::rast(masked_fp))
  }
  
  # Use mask file if it already exists
  if(file.exists(msk_all_fp) & !overwrite) {
    ras <- terra::rast(in_fp)
    msk <- terra::rast(msk_all_fp)
    
    # Apply non-negotiable mask (msk_AWUwb_u20) to AGB from SAGA filter
    r_dtype <- raster::dataType(raster::raster(in_fp))
    r_masked <- ras %>% 
      terra::mask(msk, filename = masked_fp, overwrite = T, 
                  wopt = list(datatype = r_dtype, gdal = 'COMPRESS=LZW'))
    
    # Return
    return(r_masked)
  }
  
  # Load AGB
  ras <- terra::rast(in_fp)
  
  # Initialize mask raster (all 1's)
  msk <- ras %>% terra::classify(rbind(c(-Inf, Inf, 1)))
  
  # AGB <20 mask
  if('u20' %in% masks){
    msk_u20 <- ras %>% terra::classify(rbind(c(-Inf, 20, NA),
                                                 c(20, Inf, 1)))
    msk <- msk %>% terra::mask(msk_u20)
  }
  
  # LC17 Water and Urban, and OSM water with 25 m buffer 
  if('WU' %in% masks & 'wb' %in% masks){
    msk_WUwb_fp <- file.path(masks_dir, "mask_WaterUrban_water25.tif")
    msk_WUwb <- terra::rast(msk_WUwb_fp) %>% terra::crop(ras)
    
    # Check extents
    if(ext(ras) != ext(msk_WUwb)){
      print("Extents don't match.")
    }
    
    msk <- msk %>% terra::mask(msk_WUwb)
    
  } else if('WU' %in% masks) {
    msk_fp <- file.path(masks_dir, "mask_WaterUrban.tif")
    msk_WU <- terra::rast(msk_fp) %>% terra::crop(ras)
    
    # Check extents
    if(ext(ras) != ext(msk_WU)){
      print("Extents don't match.")
    }
    
    msk <- msk %>% terra::mask(msk_WU)
    
  } else if('U' %in% masks) {
    msk_fp <- file.path(masks_dir, "mask_Urban.tif")
    msk_U <- terra::rast(msk_fp) %>% terra::crop(ras)
    
    # Check extents
    if(ext(ras) != ext(msk_U)){
      print("Extents don't match.")
    }
    
    msk <- msk %>% terra::mask(msk_U)
    
  }
  
  # ALOS mask 
  if('A' %in% masks){
    msk_fp <- file.path(masks_dir, "mask_palsar_normal_2019.tif")
    msk_tmp <- terra::rast(msk_fp) %>% terra::crop(ras)
    
    msk <- msk %>% terra::mask(msk_tmp)
    
  } else if('L' %in% masks) {
    msk_fp <- file.path(masks_dir, "mask_palsar_layover_2019.tif")
    msk_tmp <- terra::rast(msk_fp) %>% terra::crop(ras)
    
    msk <- msk %>% terra::mask(msk_tmp)
  }
  
  # Save mask raster
  msk %>% terra::writeRaster(filename = msk_all_fp, 
                             overwrite = TRUE, 
                             datatype = 'INT1S', 
                             gdal = 'COMPRESS=LZW')
  
  # Apply non-negotiable mask (msk_AWUwb_u20) to AGB from SAGA filter
  r_dtype <- raster::dataType(raster::raster(in_fp))
  r_masked <- ras %>% 
    terra::mask(msk, 
                filename = masked_fp, 
                overwrite = T, 
                datatype = r_dtype, 
                gdal = 'COMPRESS=LZW')
  
  # Return
  return(r_masked)
}

filter_g0 <- function(in_fp, filter = 'median', params=list(filtsize=5), cap=NA) {
  
  r_dtype <- raster::dataType(raster::raster(in_fp))
  
  # cap <- 0.3
  if(!is.na(cap)) {
    out_fp <- str_c(tools::file_path_sans_ext(in_fp), 
                    str_glue('_cap{cap}.tif'))
    g0 <- terra::rast(in_fp) %>% 
      terra::classify(rbind(c(cap, Inf, cap)), 
                      filename = out_fp, 
                      overwrite = TRUE,
                      datatype = r_dtype, 
                      gdal = 'COMPRESS = DEFLATE')
    
    in_fp <- out_fp
  }
  
  # params <- list(filtsize = 11, sigma = 10)
  if(filter == 'lee') {
    # Lee sigma filter
    out_fp <- str_c(tools::file_path_sans_ext(in_fp), 
                    str_glue('_lee{params$filtsize}s{params$sigma}.tif'))
    wbt_lee_sigma_filter(in_fp, 
                         out_fp, 
                         filterx = params$filtsize, 
                         filtery = params$filtsize,
                         sigma = params$sigma,
                         verbose_mode = FALSE, 
                         compress_rasters = TRUE)
    
    in_fp <- out_fp
  } 
  
  # params <- list(filtsize = 5)
  if(filter == 'median') {
    # Median filter - good smoothing, preserves edges
    out_fp <- str_c(tools::file_path_sans_ext(in_fp), 
                    str_glue('_med{params$filtsize}.tif'))
    g0_filt <- terra::rast(in_fp) %>% 
      terra::focal(w = params$filtsize, 
                   fun = 'median', 
                   na.rm = TRUE, 
                   filename = out_fp, 
                   overwrite = TRUE,
                   datatype = r_dtype, 
                   gdal = 'COMPRESS = DEFLATE')
    
    in_fp <- out_fp
  }
  
  # params <- list(filtsize = 13)
  if(filter == 'conservative') {
    # Conservative smoothing - minimal changes, but reduces local maximum
    out_fp <- str_c(tools::file_path_sans_ext(in_fp), 
                    str_glue('_conserv{params$filtsize}.tif'))
    wbt_conservative_smoothing_filter(in_fp, 
                                      out_fp, 
                                      filterx = params$filtsize, 
                                      filtery = params$filtsize,
                                      verbose_mode = FALSE, 
                                      compress_rasters = TRUE)
    
    in_fp <- out_fp
  }
  
  if(filter == 'aggregate') {
    # recommended by Saatchi 2015 and performed by Michelakis et al. 2015
    out_fp <- str_c(tools::file_path_sans_ext(in_fp), '_agg50.tif')
    if(!file.exists(out_fp)) {
      g0 <- terra::rast(in_fp) %>% 
        terra::aggregate(fact = 2, 
                         fun = mean, 
                         na.rm = TRUE,
                         filename = out_fp,
                         overwrite = TRUE,
                         datatype = r_dtype, 
                         gdal = 'COMPRESS = DEFLATE')
    }
  }
  
  return(out_fp)
}

fill_gaps_from_polygons <- function(filt_ras, masked_fp, summary_pols_fp, out_fp, 
                                    lc_stat = 'median') {
  
  # Stop if file already exists
  if(file.exists(out_fp)) return(terra::rast(out_fp))
  
  # SF to SpatVector
  lc_all_vect <- terra::vect(summary_pols_fp)
  
  # Rasterize medians
  lcpatch_ras <- lc_all_vect %>% terra::rasterize(filt_ras, field = lc_stat)
  
  # Fill gaps and save
  ras_dtype <- raster::dataType(raster::raster(masked_fp))
  masked_ras <- terra::rast(masked_fp)
  agb_filled <- terra::cover(masked_ras, 
                             lcpatch_ras, 
                             filename = out_fp, 
                             overwrite = TRUE, 
                             datatype = ras_dtype)
  
  return(agb_filled)
}

summarize_raster_by_polygons <- function(lc_pols_fp, agb_ras, lc_pols_agb_fp) {
  
  if(file.exists(lc_pols_agb_fp)) return(lc_pols_agb_fp)
  
  # Load polygons
  lc_vect <- terra::vect(lc_pols_fp)
  
  # Don't use Water and Urban class
  lc_vect2 <- terra::subset(lc_vect, lc_vect$LC > 2) 
  
  # Extract AGB to polygons
  agb_ex <- terra::extract(agb_ras, lc_vect2) 
  names(agb_ex) <- c('ID', 'agb')
  
  # Convert matrix to tibble and get mean and count for each polygon
  agb <- agb_ex %>% 
    as_tibble() %>% 
    dplyr::group_by(ID) %>% 
    dplyr::summarise(median = median(agb, na.rm=T),
                     mean = mean(agb, na.rm=T),
                     ct = sum(!is.na(agb)), 
                     sd = sd(agb, na.rm = T))
  
  # Append columns to SF polys
  x <- terra::as.data.frame(lc_vect2, geom = 'WKT')
  lc_sf <- sf::st_as_sf(x, 
                        wkt="geometry", 
                        crs=crs(lc_vect2))
  
  # Append columns to SF polys
  lc_sf_all <- agb %>% 
    dplyr::select(-ID) %>% 
    cbind(lc_sf, .) %>% 
    dplyr::filter(ct > 1)
  
  # Save
  dir.create(dirname(lc_pols_agb_fp), showWarnings = FALSE, recursive = TRUE)
  lc_sf_all %>% st_write(lc_pols_agb_fp, append = FALSE)
  
  # Return
  return(lc_pols_agb_fp)
}

# Merge and mask PALSAR tiles -------------------------------------------------
# Downloaded from https://www.eorc.jaxa.jp/ALOS/en/palsar_fnf/data/2019/html/Grid24_fnf.htm
perform_merge = FALSE

if(!file.exists(hti_poly_fp)){
  hti_poly <- st_read("data/tidy/contextual_data/HTI_adm/HTI_adm0.shp") %>% 
    st_make_valid(hti_poly)
  hti_poly %>% st_write(hti_poly_fp, append=FALSE)
}

if(perform_merge) {
  
  # Create mosaic of mask file
  code <- 'mask'
  out_fp <- file.path(g0_dir, str_glue('{code}.tif'))
  msk_fp <- merge_tiles(code, out_fp, hti_bb, raw_dir)
  
  # Get Haiti land mask, combine PALSAR mosaic and Haiti boundary
  if(!file.exists(landmask_fp)){
    dn_mask <- terra::rast(msk_fp)
    dn_mask[dn_mask<75] <- NA
    dn_mask[!is.na(dn_mask)] <- 1
    
    # Crop
    msk_poly <- terra::vect(hti_poly_fp)
    dn_mask %>% 
      terra::crop(msk_poly) %>% 
      terra::mask(msk_poly, 
                  inverse = FALSE, 
                  filename = landmask_fp,
                  overwrite = TRUE,
                  wopt = list(datatype='INT1U', gdal='COMPRESS=LZW'))
  }
  
  # Merge and mask
  code_list <- c('sl_HV', 'sl_HH', 'linci', 'date')
  
  # Mosaic tiles with merge_tiles
  code_list %>% purrr::walk(merge_tiles, suffix, g0_dir, hti_bb, raw_dir)
  
  # Mask mosaics with the land mask created above
  code_list %>% purrr::walk(mask_to_poly, suffix, g0_dir, landmask_fp)
}

# Convert to g0 dB -------------------------------------------------------------
# The DN values can be converted to gamma-0 values in decibel unit (dB) using the following equation:
#   𝛾0 = 10log10〈𝐷𝑁2〉+ 𝐶𝐹
# where, CF is a calibration factor, and <> is the ensemble averaging. 
# The CF value is “-83.0 dB”for the PALSAR-2/PALSAR mosaic

if(!file.exists(g0_fp)) {

  # Load raw DN  
  sl_fp <- file.path(g0_dir, str_glue("sl_{polarization}.tif"))
  sl <- terra::rast(sl_fp)
  
  # Convert to decibels
  g0 = 10 * log10( sl^2 ) + -83 # agrees with biota getGamma0
  
  # Convert to natural units (from biota getGamma0)
  if(units == 'nu') {
    g0 <- 10 ^ (g0 / 10)
  }
  
  # Save
  g0 %>% terra::writeRaster(g0_fp)
  
}

# Mask ----
masks_str <- g0_mask %>% str_c(collapse = '')
g0_masked_fp <- file.path(dirname(g0_fp), str_c(code, '_mask', masks_str, '.tif'))
g0_masked <- mask_with_options(g0_fp, masks_dir, code = code, 
                     masks = g0_mask, overwrite = FALSE)

# De-noise filter -----
filt_fp <- filter_g0(g0_masked_fp, 'lee', params = list(filtsize = 11, sigma = 10))

if(go_variant == 'med5') {
  filt_fp <- filter_g0(g0_masked_fp, 'median', params = list(filtsize = 5))
}

# Interpolate by landcover ----
lc_stat <-  'median'
fn_prefix <- tools::file_path_sans_ext(filt_fp) %>% basename()
# g0_filled_fp <- file.path(mod_dir, str_c(fn_prefix, '_LCinterp.tif'))
g0_filled_fp <- str_c(tools::file_path_sans_ext(filt_fp), '_LCinterp.tif')
g0_variant <- g0_filled_fp %>% str_extract(str_glue('(?<={code}_).*(?=\\.tif)'))
mod_dir <- file.path('data', 'modeling', code, g0_variant)

summary_pols_fp <- file.path(mod_dir, 'by_landcover', str_glue('{fn_prefix}_{lc_stat}.gpkg'))
dir.create(dirname(summary_pols_fp), recursive = TRUE)

# Load g0 raster (filtered)
filt_ras <- terra::rast(filt_fp)

# Get zonal values - write as polygons with median g0 value
summarize_raster_by_polygons(lc_pols_fp, filt_ras, summary_pols_fp) 

# Use zonal value to fill gaps in g0 raster
fill_gaps_from_polygons(filt_ras, filt_fp, summary_pols_fp, g0_filled_fp)


# ~ Create masks (which are used by mask_with_options()) ----
# Masks: ALOS ----
# Replicate ALOS normal 
msk_norm_fp <- file.path(masks_dir, str_glue("mask_palsar_normal_{year}.tif"))
if(!file.exists(msk_norm_fp)) {
  msk_fp <- file.path(g0_dir, str_glue('mask{suffix}.tif'))
  msk <- terra::rast(msk_fp)
  msk_A <- msk %>% 
    # terra::crop()
    terra::classify(cbind(255,1), 
                    othersNA = TRUE,
                    filename = msk_norm_fp, 
                    overwrite = TRUE, 
                    wopt = list(datatype='INT1U', gdal='COMPRESS=LZW'))
}

# Layover mask
msk_layover_fp <- file.path(masks_dir, str_glue("mask_palsar_layover_{year}.tif"))
if(!file.exists(msk_layover_fp)) {
  msk_fp <- file.path(g0_dir, str_glue('mask{suffix}.tif'))
  msk <- terra::rast(msk_fp)
  msk_L <- msk %>% 
    terra::classify(rbind(cbind(-Inf, 99, 1), 
                    cbind(99, 101, NA), 
                    cbind(101, Inf, 1)), 
                    othersNA = TRUE,
                    filename = msk_layover_fp, 
                    overwrite = TRUE, 
                    wopt = list(datatype='INT1U', gdal='COMPRESS=LZW'))
}

# Masks: Land cover -------------------------------------------------------------------
raw_dir <- file.path('data/raw')
lc_fp <- file.path(raw_dir, 'landcover', 'Lemoiner', 'Haiti2017_Clip.tif')
lc_res_fp <- file.path(tidy_dir, 'landcover', 'Haiti2017_agbres.tif')
masks_dir <- file.path(g0_dir, 'masks')

# Resample landcover to PALSAR resolution
if(!file.exists(lc_res_fp)) {
  lc <- terra::rast(lc_fp)
  g0 <- terra::rast(g0_fp)
  lc <- terra::resample(lc, g0, method='near',
                        filename = lc_res_fp, 
                        overwrite = TRUE, 
                        wopt = list(datatype='INT1U', gdal='COMPRESS=LZW'))
}

lc <- terra::rast(lc_res_fp)

# Urban mask
msk_U_fp <- file.path(masks_dir, "mask_Urban.tif")
msk_U <- lc %>% 
  terra::classify(rbind(cbind(0,1.5,1), 
                        cbind(3,7,1)), 
                  include.lowest = TRUE,
                  othersNA = TRUE,
                  filename = msk_U_fp, 
                  overwrite = TRUE, 
                  wopt = list(datatype='INT1U', gdal='COMPRESS=LZW'))

# WaterUrban mask
msk_WU_fp <- file.path(masks_dir, "mask_WaterUrban.tif")
msk_WU <- lc %>% 
  terra::classify(cbind(3,7,1), 
                  include.lowest = TRUE,
                  othersNA = TRUE,
                  filename = msk_WU_fp, 
                  overwrite = TRUE, 
                  wopt = list(datatype='INT1U', gdal='COMPRESS=LZW'))

# Mask out all but forest
mskinv_T_fp <- file.path(masks_dir, "mask_inv_TreeCover.tif")
mskinv_T <- lc %>%
  terra::classify(cbind(4,1), 
                  include.lowest = TRUE,
                  othersNA = TRUE,
                  filename = mskinv_T_fp, 
                  overwrite = TRUE, 
                  wopt = list(datatype='INT1U', gdal='COMPRESS=LZW'))

# Combine ALOS and WaterUrban masks ----
# NA== ALOS mask: non-valid ALOS pixels; 1==Normal ALOS land pixels
msk_norm_fp <- file.path(masks_dir, str_glue("mask_palsar_normal_{year}.tif"))
msk_A <- terra::rast(msk_norm_fp)
msk_WU <- terra::rast(file.path(masks_dir, "mask_WaterUrban.tif")) # NA== WaterUrban and ALOS ocean; 1==all other land

msk_AWU <- msk_WU * msk_A
msk_AWU %>% 
  terra::writeRaster(file.path(masks_dir, "mask_AWU.tif"),
                     overwrite = TRUE, 
                     wopt = list(datatype='INT1U', gdal='COMPRESS=LZW'))

# Create water mask from OSM polygons with 25 m buffer
# OSM water with 25 m buffer
water_buff_fp <- file.path(masks_dir, 'vector', 'osm_water_buff25m.shp')
if(!file.exists(water_buff_fp)){
  st_read(file.path(tidy_dir, 'contextual_data', 'OSM_free', 
                    'gis_osm_water_a_free_1.shp')) %>%
    st_transform(32618) %>%
    st_buffer(dist = 25) %>%
    summarize() %>%
    st_transform(4326) %>%
    st_write(water_buff_fp, append = FALSE)
}

msk_waterbuff_fp <- file.path(masks_dir, "mask_water_buff25m.tif")
water_polysb <- terra::vect(water_buff_fp)

# Initialize raster
msk_wb = msk_AWU
values(msk_wb) <- 1
msk_wb <- msk_wb %>% 
  terra::mask(water_polysb, inverse=TRUE,
              filename = msk_waterbuff_fp, 
              overwrite = TRUE, 
              wopt = list(datatype='INT1U', gdal='COMPRESS=LZW'))
# names(msk_wb) <- 'Mask'
msk_wb %>% writeRaster(msk_waterbuff_fp,
                       overwrite = TRUE,
                       wopt = list(datatype='INT1U', gdal='COMPRESS=LZW'))

# Mask all water
msk_Ww_fp <- file.path(masks_dir, "mask_allwater.tif")
msk_waterbuff_fp <- file.path(masks_dir, "mask_water_buff25m.tif")

# Make Water mask from landcover
lc <- terra::rast(lc_res_fp)
msk_W_fp <- file.path(masks_dir, "mask_Water.tif")
msk_W <- lc %>% 
  terra::classify(cbind(2,7,1), 
                  include.lowest = TRUE,
                  othersNA = TRUE,
                  filename = msk_W_fp, 
                  overwrite = TRUE, 
                  wopt = list(datatype='INT1U', gdal='COMPRESS=LZW'))

# Combine Water and water buffer masks
msk_Ww <- msk_W * msk_wb
msk_Ww %>% writeRaster(msk_Ww_fp,
                       overwrite = TRUE,
                       wopt = list(datatype='INT1U', gdal='COMPRESS=LZW'))

# Make default mask - water and urban
msk_WUw_fp <- file.path(masks_dir, "mask_WaterUrban_water25.tif")
msk_WUw <- msk_WU * msk_wb
msk_WUw %>% writeRaster(msk_WUw_fp,
                       overwrite = TRUE,
                       wopt = list(datatype='INT1U', gdal='COMPRESS=LZW'))
