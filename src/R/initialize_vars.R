
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
agb_code <- 'l2_WU'
saturation_pt <- 300
g0_variant <- 'med5'
if(is.na(g0_variant) | g0_variant == '') g0_variant <- 'simple'

# Filepaths
results_dir <- 'data/results'
tidy_dir <- 'data/tidy'

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
lc_pols_fp <- file.path(tidy_dir, "landcover/Haiti2017_Clip_polys.gpkg")

# 03
g0_dir <- file.path(tidy_dir, str_c('palsar_', year))
modeling_dir <- file.path('data/modeling', code)
# var_order <- c('simple', 'cappt2_conserv13', 'cappt2_conserv13_mean5', 'med5', 'lee11s10', 'maskLU_lee11s10_LCinterp')
agb_l0_fp <- file.path(modeling_dir, g0_variant, str_c("agb_l0_", g0_variant, ".tif"))

# 04 
masks_dir <- file.path('data', 'tidy', str_c('palsar_', year), 'masks')
landmask_fp <- file.path(masks_dir, 'hti_land_palsar.tif')

# 05
# Set variables ----
rmse_cv <- 21.2
lc_stat <-  'median'
input_level <- 'l2_mask'
mask_level <- 'LU'

# Input filepaths
agb_dir <- file.path('data', 'modeling', code, g0_variant)
(agb_capped_fp <- list.files(agb_dir, str_c('agb_', input_level, '.*\\.tif'), full.names = TRUE))
(agb_masked_fp <- list.files(agb_dir, str_glue('agb_l2.*{mask_level}\\.tif'), 
                             full.names = TRUE))

lc_fps <- c("data/raw/landcover/Lemoiner/Haiti2017_Clip.tif", 
            "data/raw/landcover/Lemoiner/DR_2017_clip.tif")

# Output filepaths

agb_by_lc_prefix <- file.path(agb_dir, 'agb_by_landcover', 
                              str_glue('agb_{input_level}_{mask_level}_{lc_stat}_byLC'))
lc_pols_agb_fp <- str_c(agb_by_lc_prefix, '.gpkg')
agb_by_lc_fp <- str_c(agb_by_lc_prefix, '.tif') 
agb_by_lc_sd_fp <- str_c(agb_by_lc_prefix, '_sd.tif') 

agb_filled_fp <- file.path(agb_dir, str_glue('agb_l3_fillLC{lc_stat}_{input_level}_{mask_level}.tif'))
agb_filled_sd_fp <- file.path(agb_dir, str_glue('agb_l3_fillLC{lc_stat}_{input_level}_{mask_level}_sd.tif'))


# 07
agb_dir <- file.path(modeling_dir, g0_variant)

# 07 - External map filepaths
glob_fp <- file.path(tidy_dir, 'biomass_maps', "GlobBiomass/N40W100_agb_crop_hti.tif")
# esa_fp <- file.path(tidy_dir, 'biomass_maps', "ESA_CCI/ESA_agb17_crop_hti.tif")
# avit_fp <- file.path(tidy_dir, 'biomass_maps', "Avitabile/Avitabile_AGB_crop_hti.tif")
# bacc_fp <- file.path(tidy_dir, 'biomass_maps', "Baccini/20N_080W_t_aboveground_biomass_ha_2000_crop_hti.tif")

# (agb_fps <- list.files(agb_dir, str_c('agb_', agb_input_level, '.*[^(sd)]\\.tif'), 
#                        full.names = TRUE))
# agb_fp <- agb_fps[[2]]

