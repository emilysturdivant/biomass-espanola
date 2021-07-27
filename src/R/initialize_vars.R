
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
saturation_pt <- 300
g0_variant <- 'med5'
if(is.na(g0_variant) | g0_variant == '') g0_variant <- 'simple'

# Filepaths
results_dir <- 'data/results'
tidy_dir <- 'data/tidy'
raw_ext_maps_dir <- 'data/raw/biomass_maps'
raw_lc_dir <- "data/raw/landcover"
tidy_maps_dir <- 'data/tidy/biomass_maps'
tidy_lc_dir <- 'data/tidy/landcover'

# 01
stems_fp_in <- "data/species_and_wds/haiti_data_wds2.csv"
plots_fp_in <- "data/species_and_wds/mplots_geoms.csv"
plots_shp <- file.path(tidy_dir, 'survey_plots', 'all_plots.shp')
# mstems_fp <- file.path(tidy_dir, 'survey_plots/mstems_agb_alldbh.rds')
# mstems_fp <- file.path(tidy_dir, 'survey_plots/mstems_agb_noXtrms.rds')
mstems_fp <- file.path(tidy_dir, 'survey_plots/mstems_agb.rds')
# field_agb_fp <- file.path(tidy_dir, 'survey_plots', 'plots_agb_noXtrms.rds')
field_agb_fp <- file.path(tidy_dir, 'survey_plots', 'plots_agb.rds')

# 02
hti_poly_fp <- file.path(tidy_dir, "contextual_data/HTI_adm/HTI_adm0_fix.shp")
raw_dir <- file.path('data/raw/ALOS', year)
g0_dir <- file.path(tidy_dir, str_c('palsar_', year))
masks_dir <- file.path(g0_dir, 'masks')
landmask_fp <- file.path(masks_dir, 'hti_land_palsar.tif')
lc_pols_fp <- file.path(tidy_dir, "landcover/Lemoiner/Haiti2017_Clip_polys.gpkg")
lc_res_fp <- file.path(tidy_dir, 'landcover', 'Lemoiner', 'Haiti2017_agbres.tif')

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
