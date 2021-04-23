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

# Load libraries
#library(silvr)
# library(readr)
library(tidyverse)
# library(ggridges)
# library(rgdal)
library(tmap)
# library(here)
library(gdalUtils)
library(terra)
library(sf)
# library(stars)


results_dir <- 'data/results'
# masks_dir <- 'data/tidied/alos18_masks'
suffix <- ''
year <- '2019'
raw_dir <- file.path('data/raw/ALOS', year)
tidy_dir <- file.path('data/tidied', str_c('palsar', substr(year, 3, 4)))
masks_dir <- file.path(tidy_dir, 'masks')
landmask_fp <- file.path(masks_dir, 'hti_land_palsar.tif')
hti_poly_fp <- "data/contextual_data/HTI_adm/HTI_adm0_fix.shp"

hisp_bb <- st_bbox(c(xmin = -74.48133, ymax = 20.09044, 
                     xmax = -68.32267, ymin = 17.47022))
hti_bb <- st_bbox(c(xmin = -74.48133, ymax = 20.09044, 
                    xmax = -71.61815, ymin = 18.02180))

# Merge and mask PALSAR tiles -------------------------------------------------
# Downloaded from https://www.eorc.jaxa.jp/ALOS/en/palsar_fnf/data/2019/html/Grid24_fnf.htm

perform_merge = FALSE

merge_tiles <- function(code, suffix, tidy_dir, hti_bb, raw_dir) {
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

mask_raster <- function(code, suffix, tidy_dir, landmask_fp) { 
  # Set output filename
  out_fp <- file.path(tidy_dir, str_glue('{code}{suffix}.tif'))
  
  # Mask
  dtype <- if(code == 'mask' | code == 'linci') 'INT1U' else 'INT2U'
  lnd_msk <- terra::rast(landmask_fp)
  r <- terra::rast(out_fp)
  r %>% terra::mask(lnd_msk, 
                    filename = out_fp,
                    overwrite = TRUE,
                    wopt = list(datatype = dtype, gdal = 'COMPRESS=LZW'))
}

if(!file.exists(hti_poly_fp)){
  hti_poly <- st_read("data/contextual_data/HTI_adm/HTI_adm0.shp") %>% 
    st_make_valid(hti_poly)
  hti_poly %>% st_write(hti_poly_fp, append=FALSE)
}

if(perform_merge) {
  
  # Create mosaic of mask file
  code <- 'mask'
  out_fp <- file.path(tidy_dir, str_glue('{code}.tif'))
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
      terra::mask(msk_poly, inverse = FALSE, 
                  filename = landmask_fp,
                  overwrite = TRUE,
                  wopt = list(datatype='INT1U', gdal='COMPRESS=LZW'))
  }
  
  # Merge and mask
  code_list <- c('sl_HV', 'sl_HH', 'linci', 'date')
  
  # Mosaic tiles with merge_tiles
  code_list %>% purrr::walk(merge_tiles, suffix, tidy_dir, hti_bb, raw_dir)
  
  # Mask mosaics with the land mask created above
  code_list %>% purrr::walk(mask_raster, suffix, tidy_dir, landmask_fp)
}

# ~Look~ at proportion of values in each mask category----
report_fp <- file.path('data', 'reports', str_glue('palsar{year}_mask_pcts.csv'))
if(!file.exists(report_fp)) {
  
  # Load PALSAR mask
  msk_fp <- file.path(tidy_dir, str_glue('mask{suffix}.tif'))
  dn_mask <- terra::rast(msk_fp)
  
  # Get values
  rvals <- values(dn_mask)
  
  # Count pixels in each category and convert to percentage
  df <- tibble(group = c("Normal", "Layover", "Shadowing"), 
               count = c(sum(rvals == 255, na.rm=TRUE), 
                         sum(rvals == 100, na.rm=TRUE),
                         sum(rvals == 150, na.rm=TRUE)))
  
  # Convert counts to percentage of all land pixels
  landct <- sum(df$count)
  df <- df %>%
    add_row(group = 'Land', count = landct) %>% 
    mutate(pct = count / landct)
  
  # Barplot
  (bp <- ggplot(filter(df, group %in% c('Normal', 'Shadowing', 'Layover')), 
                aes(x="", y=pct, fill=group))+
      geom_bar(width = 1, stat = "identity")+ 
      theme_minimal() + coord_polar("y", start=0))
  
  # Save as CSV
  df %>% write_csv(report_fp)
}

# Replicate ALOS Normal mask ----
msk_norm_fp <- file.path(masks_dir, "mask_ALOS_normal.rds")
if(!file.exists(msk_norm_fp)) {
  msk_fp <- file.path(tidy_dir, str_glue('mask{suffix}.tif'))
  msk_A <- terra::rast(msk_fp)
  msk_A[msk_A < 254] <- NA
  msk_A[!is.na(msk_A)] <- 1
  msk_A %>% saveRDS(msk_norm_fp)
}

# Create inverse mask for > 7000 ----
g0_fp <- file.path(tidy_dir, "sl_HV.tif")

msk_norm_fp <- file.path(masks_dir, "mask_ALOS_gt7000.rds")
if(!file.exists(msk_norm_fp)) {
  g0 <- terra::rast(g0_fp)
  g0[g0 < 7000] <- NA
  g0[!is.na(g0)] <- 1
  g0 %>% saveRDS(msk_norm_fp)
}

# Aggregate to 50m ----
# recommended by Saatchi 2015 and performed by Michelakis et al. 2015
out_fp <- str_c(tools::file_path_sans_ext(g0_fp), '_agg50.tif')
if(!file.exists(out_fp)) {
  g0 <- rast(g0_fp)
  g0.agg <- aggregate(g0, fact = 2, fun = mean, na.rm = TRUE,
                      filename = out_fp,
                      overwrite = TRUE,
                      wopt = list(datatype='INT2U', gdal='COMPRESS=LZW'))
}

# De-speckle with Lee filter ----------------------------------------------------
# whitebox
# if (!require(devtools)) install.packages('devtools')
# devtools::install_github("giswqs/whiteboxR")
library('whitebox')

# Lee sigma filter
filtsize <- 11
sigma <- 10
out_fp <- str_c(tools::file_path_sans_ext(g0_fp), 
                str_glue('_smooth_f{filtsize}sig{sigma}.tif'))
wbt_lee_sigma_filter(g0_fp, 
                     out_fp, 
                     filterx = filtsize, filtery = filtsize,
                     sigma = sigma,
                     verbose_mode = FALSE, 
                     compress_rasters = TRUE)
g0_filt <- terra::rast(out_fp)
max(g0_filt)

# Median filter - good smoothing, preserves edges
filtsize <- 5
out_fp <- str_c(tools::file_path_sans_ext(g0_fp), 
                str_glue('_smooth_med{filtsize}.tif'))
wbt_median_filter(g0_fp, 
                  out_fp, 
                  filterx = filtsize, filtery = filtsize,
                  verbose_mode = FALSE, 
                  compress_rasters = TRUE)
g0_filt <- terra::rast(out_fp)
max(g0_filt)

# Saturate values at 7000 ----
capval <- 7000
cap_fp <- str_c(tools::file_path_sans_ext(g0_fp), 
                str_glue('_cap{capval}.tif'))
g0 <- terra::rast(g0_fp)
g0[g0 > capval] <- capval + 1
terra::writeRaster(g0,
                   cap_fp, 
                   overwrite = TRUE,
                   wopt = list(datatype='INT2U', gdal='COMPRESS=LZW'))

# Conservative smoothing - minimal changes, but reduces local maximum
filtsize <- 13
out_fp2 <- str_c(tools::file_path_sans_ext(cap_fp), 
                str_glue('_conserv{filtsize}.tif'))
wbt_conservative_smoothing_filter(cap_fp, 
                                  out_fp2, 
                                  filterx = filtsize, filtery = filtsize,
                                  verbose_mode = FALSE, 
                                  compress_rasters = TRUE)
g0_filt <- terra::rast(out_fp2)
max(g0_filt)

# Median filter
filtsize <- 3
out_fp3 <- str_c(tools::file_path_sans_ext(cap_fp), 
                 str_glue('_med{filtsize}.tif'))
wbt_median_filter(cap_fp, 
                  out_fp3, 
                  filterx = filtsize, filtery = filtsize,
                  verbose_mode = FALSE, 
                  compress_rasters = TRUE)
g0_filt <- terra::rast(out_fp3)
max(g0_filt)

# Lee filter
filtsize <- 11
sigma <- 10
out_fp3 <- str_c(tools::file_path_sans_ext(cap_fp), 
                 str_glue('_lee{filtsize}sig{sigma}.tif'))
wbt_lee_sigma_filter(cap_fp, 
                     out_fp3, 
                     filterx = filtsize, filtery = filtsize,
                     sigma = sigma,
                     verbose_mode = FALSE, 
                     compress_rasters = TRUE)
g0_filt <- terra::rast(out_fp3)
max(g0_filt)




# 
# 
# 
# 
# 
# # Make masks (Hispaniola extent) ------------------------------------------------------------
# # Load LC17 masked to ALOS land 
# lc <- readRDS(file.path(results_dir, "R_out/LC17_masked_to_ALOS_land_stars.rds"))
# lc <- raster("data/LULC/Hisp_2017_resALOS_mskLand.tif")
# 
# # WaterUrban mask
# msk_WU <- lc
# msk_WU[msk_WU<3] <- NA
# msk_WU[!is.na(msk_WU)] <- 1
# msk_WU %>% writeRaster(file.path(masks_dir, "mask_WaterUrban_raster.tif"))
# # msk_WU %>% saveRDS(file.path(results_dir, "R_out/mask_WaterUrban_raster.rds"))
# 
# # Water mask
# msk_W <- lc
# msk_W[msk_W==1] <- NA
# msk_W[!is.na(msk_W)] <- 1
# msk_W %>% writeRaster(file.path(masks_dir, "mask_Water_raster.tif"))
# msk_W %>% saveRDS(file.path(results_dir, "R_out/mask_WaterLC17_raster.rds"))
# 
# # Mask out all but forest
# mskinv_T <- lc
# mskinv_T[mskinv_T!=4] <- NA
# mskinv_T[mskinv_T==4] <- 1
# mskinv_T %>% writeRaster(file.path(masks_dir, "mask_inv_TreeCover.tif"))
# mskinv_T %>% saveRDS(file.path(results_dir, "R_out/mask_inv_TreeCover_raster.rds"))
# 
# # Mask out Bareland
# msk_B <- lc
# msk_B[msk_B!=3] <- 1
# msk_B[msk_B==3] <- NA
# msk_B %>% saveRDS(file.path(results_dir, "R_out/mask_Bareland_raster.rds"))
# 
# # Mask out all but grassland and shrubs
# mskinv_GS <- lc
# mskinv_GS[mskinv_GS<5] <- NA
# mskinv_GS[mskinv_GS>4] <- 1
# mskinv_GS %>% saveRDS(file.path(results_dir, "R_out/mask_inv_GrasslandShrubs_stars.rds"))
# 
# # Mask out all but tree cover, grassland and shrubs
# mskinv_GS <- lc
# mskinv_GS[mskinv_GS<4] <- NA
# mskinv_GS[mskinv_GS>3] <- 1
# mskinv_GS %>% saveRDS(file.path(results_dir, "R_out/mask_inv_LC17_vegTreeCGrasslandShrubs_stars.rds"))
# 
# # Presence (Inverse mask) of LC17 Water
# mskinv_W <- lc
# mskinv_W[mskinv_W!=1] <- NA
# mskinv_W %>% saveRDS(file.path(results_dir, "R_out/mask_inv_ALOS_Water_raster.rds"))
# 
# # Combine ALOS and WaterUrban masks
# msk_A <- readRDS(file.path(results_dir, "R_out/mask_ALOS_stars.rds"))
# msk_A %>% as("Raster") %>% writeRaster(file.path(masks_dir, "mask_ALOS18.tif"))
# 
# rm(list=ls()) 
# 
# msk_A <-  # NA== ALOS mask: non-valid ALOS pixels; 1==Normal ALOS land pixels # HISPANIOLA
#   raster(file.path(masks_dir, "mask_ALOS18.tif"))
# msk_WU <- # NA== WaterUrban and ALOS ocean; 1==all other land
#   raster(file.path(masks_dir, "mask_WaterUrban_raster.tif"))
# 
# msk_AWU <- msk_WU * msk_A
# msk_AWU %>% writeRaster(file.path(masks_dir, "mask_ALOS_WaterUrban.tif"))
# msk_AWU %>% saveRDS(file.path(results_dir, "R_out/mask_ALOS_WaterUrban_raster.rds"))
# 
# # Presence (Inverse mask) of LC17 Water/Urban with ALOS mask applied
# msk_WU <- lc
# msk_WU[msk_WU==2] <- 1 # Water is already 1 and now urban is as well
# msk_WU[msk_WU!=1] <- NA
# mskinv_WU <- msk_WU*msk_A
# mskinv_WU %>% # 1==where WaterUrban overlap valid ALOS values
#   saveRDS(file.path(results_dir, "R_out/mask_inv_ALOS_WU_raster.rds") )
# 
# # Create water mask from OSM polygons with 25 m buffer
# # Create OSM water with 25 m buffer
# st_read('data/contextual_data/OSM_free/gis_osm_water_a_free_1.shp') %>% 
#   st_transform(32618) %>% 
#   st_buffer(dist = 25) %>% 
#   summarize() %>% 
#   st_transform(4326) %>% 
#   st_write(file.path(masks_dir, 'vector/osm_water_buff25m.shp', append=FALSE))
# water_polysb <- st_read(file.path(masks_dir, 'vector/osm_water_buff25m.shp'))
# 
# water_polysb <- st_read(file.path(masks_dir, 'vector/osm_water_buff25m.shp'))
# msk_wb = msk_AWU
# values(msk_wb) <- 1
# msk_wb <- msk_wb %>% 
#   mask(water_polysb, inverse=TRUE)
# names(msk_wb) <- 'Mask'
# msk_wb %>% writeRaster(file.path(masks_dir, "mask_osm_water_buff25m.tif"))
# 
# # Mask OSM water 25 (add to ALOS mask)
# msk_A <-  # NA== ALOS mask: non-valid ALOS pixels; 1==Normal ALOS land pixels # HISPANIOLA
#   raster(file.path(masks_dir, "mask_ALOS18.tif"))
# msk_Aw <- msk_A %>% 
#   mask(water_polysb, inverse=TRUE)
# msk_Aw %>% saveRDS(file.path(results_dir, "R_out/mask_ALOS_OSMwater25_raster.rds"))
# 
# # Water from OSM polygons
# water_polys <- st_read('data/contextual_data/OSM_free/gis_osm_water_a_free_1.shp')
# mskinv_WP <- msk_A %>% 
#   mask(water_polys, inverse=FALSE)
# mskinv_WP %>% saveRDS(file.path(results_dir, "R_out/mask_inv_ALOS_OSMwater_raster.rds"))
# 
# # Water from OSM polygons with 25 m buffer
# mskinv_WPb <- msk_A %>% 
#   mask(water_polysb, inverse=FALSE)
# mskinv_WPb %>% saveRDS(file.path(results_dir, "R_out/mask_inv_OSMwater25mbuffer_raster.rds"))
# 
# # Create inverse mask of WU and OSM water 25 m
# mskinv_WU <- mskinv_WU %>% as("Raster")
# mskinv_WPb <- mskinv_WPb %>% st_as_stars()
# mskinv_WUWPb <- mskinv_WU + mskinv_WPb
# mskinv_WUWPb %>% saveRDS(file.path(results_dir, "R_out/mask_inv_WUWPb_raster.rds"))
# 
# # Mask all water
# msk_W <- readRDS(file.path(results_dir, "R_out/mask_WaterLC17_raster.rds"))
# msk_wb <- raster(file.path(masks_dir, "mask_osm_water_buff25m.tif"))
# msk_Ww <- msk_W * msk_wb 
# msk_Ww %>% saveRDS(file.path(results_dir, "R_out/mask_allwater_raster.rds"))
# 
# # Create mask of backscatter >0.3
# msk_p3 <- # Initialize mask
#   read_stars(file.path(results_dir, "g0nu_HV/g0nu_2018_HV.tif"))
# msk_p3[msk_p3>0.3] <- NA
# msk_p3[msk_p3<=0.3] <- 1
# msk_p3 %>% saveRDS(file.path(results_dir, "R_out/mask_ALOSoverpt3_stars.rds"))
# 
