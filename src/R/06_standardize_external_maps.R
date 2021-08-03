# ******************************************************************************
# Script to:
#     * Standardize external (AGB) maps
# Proceeds:
#     * 4. postprocess_AGB_map.R - performs post-processing such as masking
#     * 3. regression_AGB-g0.R - creates AGB map
#     * 2. calculate_AGB.R - calculates AGB by plot from field data
#     * 1. process_ALOS_tiles.R
# Preceeds:
#     * post_external_AGB_comparisons.R -  Compare output AGB map to other biomass estimates
# Requires:
#     * AGB map (agbr)
#     * external AGB maps
# ******************************************************************************
# Load libraries ----
library('gdalUtils')
library('terra')
library('tidyverse')

# Set variables ----
source('src/R/initialize_vars.R')

# g0_variant <- 'med5'
# year <- '2019'
# input_level <- 'l3'
input_level <- agb_input_level

# Load AGB map ----
agb_ras <- terra::rast(agb_fp)
e <- ext(agb_ras)
bb <- c(e$xmin, e$ymin, e$xmax, e$ymax)

land_poly <- terra::vect(hti_poly_fp)

# FUNCTIONS -----------------------------------------------------------------------------------
crop_and_set_NA <- function(in_fp, crop_fp, out_fp = NA, bb){
  # Crop to given extent and convert 0s to NA
  r <- terra::rast(in_fp)
  gdalUtils::gdalwarp(srcfile = in_fp, 
                      dstfile = crop_fp, 
                      te = bb, 
                      tr = c(xres(r), yres(r)), 
                      tap = TRUE, 
                      overwrite = TRUE)
  
  if(!is.na(out_fp)) {
    # Set 0 to NA
    in_dtype <- raster::dataType(raster::raster(crop_fp))
    r <- terra::rast(crop_fp)
    terra::NAflag(r) <- 0
    r %>% terra::writeRaster(out_fp, overwrite=T, datatype = in_dtype)
  }
}

crop_and_mask_to_polygon <- function(in_fp, msk_poly, out_fp, na = NA, dtype = 'INT2S'){
  # Function to crop and mask raster to polygon
  
  # Filepath for temp file
  tmp_fp <- tempfile(pattern = sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(in_fp)), 
                     tmpdir = tempdir(), 
                     fileext = "tif")
  
  # Crop (save to temp)
  r <- terra::rast(in_fp)
  e <- ext(msk_poly)
  bb <- c(e$xmin, e$ymin, e$xmax, e$ymax)
  gdalUtils::gdalwarp(srcfile=in_fp, 
                      dstfile = tmp_fp,
                      te = bb, 
                      tr = c(xres(r), yres(r)), 
                      tap = TRUE, 
                      overwrite = TRUE)
  
  # Get cropped file
  r <- terra::rast(tmp_fp)
  
  # Set 0 to NA
  if(!is.na(na)) {
    terra::NAflag(r) <- na
  }
  
  # Mask
  if (is.na(dtype)) dtype <- raster::dataType(raster::raster(tmp_fp))
  r %>% 
    terra::mask(msk_poly, 
                inverse=FALSE,
                filename = out_fp, 
                overwrite = TRUE, 
                datatype = dtype) 
}

# ~ Landcover: ================================================================
# CCI Landcover ----
# Download CCI LC 2015 v.2.0.7 using FileZilla 
# Crop to our Haiti extent and convert 0s to NA
in_fp <- file.path(raw_lc_dir, "ESA_CCI/ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7.tif")
lc_fp <- file.path(tidy_lc_dir, "ESA_CCI/CCI_LC_L4_2015_v207_hti.tif")
crop_and_mask_to_polygon(in_fp, land_poly, lc_fp, na = 0)

# Download CCI LC v2.1.1 using CDS API in Python
library('reticulate')
py_run_file('src/python/download_esa_landcover.py')

# Land cover v2.1.1 from netCDF file 
in_fp <- file.path(raw_lc_dir, "ESA_CCI/C3S-LC-L4-LCCS-Map-300m-P1Y-2019-v2.1.1.nc")
crop_fp <- file.path(tidy_lc_dir, "ESA_CCI/C3S-LC-L4-LCCS-Map-300m-P1Y-2019-v2.1.1_crop.tif")
lc_fp <- file.path(tidy_lc_dir, "ESA_CCI/CCI_LC_L4_2019_v211_hti.tif")
r <- terra::rast(str_glue('NETCDF:"{in_fp}":lccs_class'))
r %>% 
  terra::crop(agb_ras) %>% 
  as("Raster") %>%
  writeRaster(crop_fp, options=c("dstnodata=-99999"), overwrite=T)
crop_and_mask_to_polygon(crop_fp, land_poly, lc_fp, na = 0)

# Land cover v2.1.1 from netCDF file 
in_fp <- file.path(raw_lc_dir, "ESA_CCI/C3S-LC-L4-LCCS-Map-300m-P1Y-2017-v2.1.1.nc")
crop_fp <- file.path(tidy_lc_dir, "ESA_CCI/C3S-LC-L4-LCCS-Map-300m-P1Y-2017-v2.1.1_crop.tif")
lc_fp <- file.path(tidy_lc_dir, "ESA_CCI/CCI_LC_L4_2017_v211_hti.tif")
r <- terra::rast(str_glue('NETCDF:"{in_fp}":lccs_class'))
r %>% 
  terra::crop(agb_ras) %>% 
  as("Raster") %>%
  writeRaster(crop_fp, options=c("dstnodata=-99999"), overwrite=T)
crop_and_mask_to_polygon(crop_fp, land_poly, lc_fp, na = 0)

# Resample LC to ESA grid ----
in_fp <- lc_fps$haiti
lc_res_fp <- file.path(tidy_lc_dir, "Lemoiner", "Lemoiner_lc17_hti_resCCI.tif")
resample_to_raster(in_fp, agb_fps$esa$fp, lc_res_fp, method = 'near')

# Pct. pix in each LC class ----
lc <- terra::rast(lc_res_fp)
rvals <- values(lc)
df <-
  tibble(group = c("Water", "Urban", "Bareland", "Tree cover", "Grassland", "Shrubs"),
         count = c(sum(rvals==1, na.rm=TRUE),
                   sum(rvals==2, na.rm=TRUE),
                   sum(rvals==3, na.rm=TRUE),
                   sum(rvals==4, na.rm=TRUE),
                   sum(rvals==5, na.rm=TRUE),
                   sum(rvals==6, na.rm=TRUE))) %>%
  mutate(pct = count / sum(count))

df %>% write_csv(file.path('data/reports', '06_lc_pixel_counts_ALOSres.csv'))

# ~ Biomass: ====================================================================
# Resample our AGB to ESA grid for comparison ----
agb_res_fp <- str_c(tools::file_path_sans_ext(agb_fps$internal$fp),
                    '_resCCI.tif')
resample_to_raster(agb_fps$internal$fp, agb_fps$esa$fp, agb_res_fp)

# GlobBiomass ------
# Crop GlobBiomass to our Hispaniola extent and convert 0s to NA
in_fp <- file.path(raw_ext_maps_dir, "GlobBiomass/N40W100_agb.tif")
glob_fp <- file.path(tidy_maps_dir, "GlobBiomass/Glob_agb10_hti.tif")
crop_and_mask_to_polygon(in_fp, land_poly, glob_fp, na = 0)

# Do the same for error (per-pixel uncertainty expressed as standard error in m3/ha)
in_fp <- file.path(raw_ext_maps_dir, "GlobBiomass/N40W100_agb_err.tif")
glob_fp <- file.path(tidy_maps_dir, "GlobBiomass/Glob_agb10_err_hti.tif")
crop_and_mask_to_polygon(in_fp, land_poly, glob_fp, na = 0)

# ESA ------
# http://data.ceda.ac.uk/neodc/esacci/biomass/data/agb/maps/2017/v1.0/geotiff
# ftp://anon-ftp.ceda.ac.uk/neodc/esacci/biomass/data/agb/maps/2017/v1.0/geotiff/
# Crop ESA to our Hispaniola extent and convert 0s to NA
in_fp <- file.path(raw_ext_maps_dir, "ESA_CCI/N40W100_ESACCI-BIOMASS-L4-AGB-MERGED-100m-2017-fv1.0.tif")
esa_fp <- file.path(tidy_maps_dir, "ESA_CCI/ESA_agb17_hti.tif")
crop_and_mask_to_polygon(in_fp, land_poly, esa_fp, na = 0)

# Do the same for SD
in_fp <- file.path(raw_ext_maps_dir, "ESA_CCI/N40W100_ESACCI-BIOMASS-L4-AGB_SD-MERGED-100m-2017-fv1.0.tif")
esa_fp <- file.path(tidy_maps_dir, "ESA_CCI/ESA_agb17_SD_hti.tif")
crop_and_mask_to_polygon(in_fp, land_poly, esa_fp, na = 0)

# Avitabile  ------
in_fp <- file.path(raw_ext_maps_dir, "Avitabile/Avitabile_AGB_Map.tif")
avit_fp <- file.path(tidy_maps_dir, "Avitabile/Avitabile_AGB_hti.tif")
crop_and_mask_to_polygon(in_fp, land_poly, avit_fp, na = NA)

# Baccini -----
# https://data.globalforestwatch.org/datasets/8f93a6f94a414f9588ce4657a39c59ff_1
# http://gfw-data.s3.amazonaws.com/climate/WHRC_biomass/WHRC_V4/Processed/20N_080W_t_aboveground_biomass_ha_2000.tif
in_fp <- file.path(raw_ext_maps_dir, "Baccini/20N_080W_t_aboveground_biomass_ha_2000.tif")
bacc_fp <- file.path(tidy_maps_dir, "Baccini/Baccini_agb00_hti.tif")
crop_and_mask_to_polygon(in_fp, land_poly, bacc_fp, na = NA)

bacc_res_fp <- file.path(tidy_maps_dir, "Baccini/Baccini_agb00_hti_resCCI.tif")
resample_to_raster(bacc_fp, agb_fps$esa$fp, bacc_res_fp)

# Biomass loss from Baccini ------
in_fp <- file.path(raw_ext_maps_dir, "Baccini/Forest_Biomass_Loss/AGLB_Deforested_Tropical_America_2000_crop.tif")
bmloss_fp <- file.path(tidy_maps_dir, "Baccini/AGLB_Deforested_Tropical_America_2000_crop.tif")
crop_and_mask_to_polygon(in_fp, land_poly, bmloss_fp, na = NA)
