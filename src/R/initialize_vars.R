
# Load libraries ----
# library('BIOMASS')
library("patchwork")
library("gridExtra")
library("sf")
library("tidyverse")

# Codes
code <- 'HV_nu'
# suffix <- 'cappt2_conserv13'
year <- '2019'
agb_input_level <- 'l2'
agb_code <- 'l2_maskWUwb'
g0_mask <- c('L', 'U', 'W', 'wb')
saturation_pt <- 300
g0_variant <- 'med5'
g0_variant <- 'maskLUWwb_med5_LCinterp'
if(is.na(g0_variant) | g0_variant == '') g0_variant <- 'simple'

# Filepaths
results_dir <- 'data/results'
tidy_dir <- 'data/tidy'
raw_ext_maps_dir <- 'data/raw/biomass_maps'
raw_lc_dir <- "data/raw/landcover"
tidy_maps_dir <- 'data/tidy/biomass_maps'
tidy_lc_dir <- 'data/tidy/landcover'

# 01 ----
stems_fp_in <- "data/species_and_wds/haiti_data_wds2.csv"
plots_fp_in <- "data/species_and_wds/mplots_geoms.csv"
plots_shp <- file.path(tidy_dir, 'survey_plots', 'all_plots.shp')
# mstems_fp <- file.path(tidy_dir, 'survey_plots/mstems_agb_alldbh.rds')
# mstems_fp <- file.path(tidy_dir, 'survey_plots/mstems_agb_noXtrms.rds')
mstems_fp <- file.path(tidy_dir, 'survey_plots/mstems_agb.rds')
# field_agb_fp <- file.path(tidy_dir, 'survey_plots', 'plots_agb_noXtrms.rds')
field_agb_fp <- file.path(tidy_dir, 'survey_plots', 'plots_agb.rds')

# 02 ----
hti_poly_fp <- file.path(tidy_dir, "contextual_data/HTI_adm/HTI_adm0_fix.shp")
raw_dir <- file.path('data/raw/ALOS', year)
g0_dir <- file.path(tidy_dir, str_c('palsar_', year))
masks_dir <- file.path(g0_dir, 'masks')
landmask_fp <- file.path(masks_dir, 'hti_land_palsar.tif')
lc_pols_fp <- file.path(tidy_dir, "landcover/Lemoiner/Haiti2017_Clip_polys.gpkg")
lc_res_fp <- file.path(tidy_dir, 'landcover', 'Lemoiner', 'Haiti2017_agbres.tif')
g0_fp <- file.path(g0_dir, 'mosaic_variants', str_glue("{code}.tif"))

# 03
g0_dir <- file.path(tidy_dir, str_c('palsar_', year))
modeling_dir <- file.path('data/modeling', code)
# var_order <- c('simple', 'cappt2_conserv13', 'cappt2_conserv13_mean5', 'med5', 'lee11s10', 'maskLU_lee11s10_LCinterp')
agb_l0_fp <- file.path(modeling_dir, g0_variant, str_c("agb_l0_", g0_variant, ".tif"))

# 04 
masks_dir <- file.path('data', 'tidy', str_c('palsar_', year), 'masks')
landmask_fp <- file.path(masks_dir, 'hti_land_palsar.tif')

# 05 ----
# Set variables 
rmse_cv <- 21.2
lc_stat <-  'median'
input_level <- 'l2_mask'
mask_level <- 'LU'

# Input filepaths
agb_dir <- file.path('data', 'modeling', code, g0_variant)
# (agb_capped_fp <- list.files(agb_dir, str_c('agb_', input_level, '.*\\.tif'), full.names = TRUE))
(agb_masked_fp <- list.files(agb_dir, str_glue('agb_l2.*{mask_level}\\.tif'), 
                             full.names = TRUE))

lc_fps <- list(haiti = file.path(raw_lc_dir, "Lemoiner/Haiti2017_Clip.tif"), 
               dr = file.path(raw_lc_dir, "Lemoiner/DR_2017_clip.tif"))

# Output filepaths
agb_by_lc_prefix <- file.path(agb_dir, 'agb_by_landcover', 
                              str_glue('agb_{input_level}_{mask_level}_{lc_stat}_byLC'))
lc_pols_agb_fp <- str_c(agb_by_lc_prefix, '.gpkg')
agb_by_lc_fp <- str_c(agb_by_lc_prefix, '.tif') 
agb_by_lc_sd_fp <- str_c(agb_by_lc_prefix, '_sd.tif') 

agb_filled_fp <- file.path(agb_dir, str_glue('agb_l3_fillLC{lc_stat}_{input_level}_{mask_level}.tif'))
agb_filled_sd_fp <- file.path(agb_dir, str_glue('agb_l3_fillLC{lc_stat}_{input_level}_{mask_level}_sd.tif'))

# 06

# 07
agb_dir <- file.path(modeling_dir, g0_variant)

# 07 - External map filepaths ----
glob_fp <- file.path(tidy_maps_dir, "GlobBiomass/Glob_agb10_hti.tif")
esa_fp <- file.path(tidy_maps_dir, "ESA_CCI/ESA_agb17_hti.tif")
avit_fp <- file.path(tidy_maps_dir, "Avitabile/Avitabile_AGB_hti.tif")
bacc_fp <- file.path(tidy_maps_dir, "Baccini/Baccini_agb00_hti.tif")
bacc_res_fp <- file.path(tidy_maps_dir, "Baccini/Baccini_agb00_hti_resCCI.tif")

agb_fp <- file.path(agb_dir, str_glue('agb_{agb_code}.tif'))

# List of AGB maps
agb_fps <- list(internal = list(name = str_glue('This study ({agb_code})'),
                                fp = agb_fp),
                glob = list(name = 'GlobBiomass',
                            fp = glob_fp), 
                esa = list(name = 'CCI', 
                           fp = esa_fp), 
                avit = list(name = 'Avitabile',
                            fp = avit_fp), 
                bacc = list(name = 'Baccini',
                            fp = bacc_res_fp))

comparison_dir <- file.path(dirname(agb_fp), 'external_comparison', agb_code)
plot_ext_csv <- file.path(comparison_dir, str_c('field_plot_means_', agb_code, '.csv'))
ext_report_csv <- file.path('data/reports', str_glue('07_ext_comparison_metrics_{agb_code}.csv'))

# AGB palettes ----
agb1_palette <- c('#4c006f', '#8d4d00', '#f7e700', '#4ee43d', '#006016')
agb1b_palette <- c('#4c006f', '#9f4d28', '#f7e700', '#4ee43d', '#006016')
agb2_palette <- c('#4c006f', '#8d6639', '#f7e700', '#4ee43d', '#006016')
agb3_palette <- c('#8d4d00', '#f7e700', '#4ee43d', '#006016')
bouvet_palette <- c('#9f4d28', '#b67633', '#cc9e45', 
                    '#e7c754', '#feee5e', '#cbdc50', 
                    '#9ac545', '#65b438', '#33a029', '#34782d')
bouvet_palette <- c('#9f4d28', '#cc9e45', '#feee5e', '#65b438', '#34782d')

agb_pal <- list(colors = agb1_palette,
                min = 0, 
                max = 120)

# Functions ----
resample_to_raster <- function(agb_fp, ext_fp, agb_res_fp, method = 'bilinear') {
  
  in_dtype <- raster::dataType(raster::raster(agb_fp))
  
  # Crop to intersection of the two extents
  out <- crop_to_intersecting_extents(terra::rast(agb_fp), terra::rast(ext_fp))
  agb_ras <- out[[1]]
  agb_ext <- out[[2]]
  
  # Resample our AGB to external resolution
  agb_res <- agb_ras %>% terra::resample(agb_ext, method=method)
  
  # Crop internal again
  agb_res <- crop_to_intersecting_extents(agb_res, terra::rast(ext_fp), return_r2=F)
  
  # Save resampled AGB
  agb_res %>% terra::writeRaster(filename = agb_res_fp, 
                                 overwrite = TRUE, 
                                 datatype = in_dtype)
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