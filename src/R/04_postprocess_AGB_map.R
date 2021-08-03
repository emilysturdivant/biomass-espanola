# ******************************************************************************
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
# ******************************************************************************
# Load libraries ----
library("terra")
library("tidyverse")

# Set variables ----
source('src/R/initialize_vars.R')
# year <- '2019'
# code <- 'HV_nu'
# saturation_pt <- 300
# modeling_dir <- file.path('data/modeling', code)
# 
# masks_dir <- file.path('data', 'tidy', str_c('palsar_', year), 'masks')
# landmask_fp <- file.path(masks_dir, 'hti_land_palsar.tif')

# Get level 0 AGB
# (dirs <- list.dirs(modeling_dir, recursive = FALSE))
# (mod_dir <- dirs[[1]])
# (agb_fp <- list.files(mod_dir, 'agb_l0.*\\.tif$', full.names = TRUE))

# # Get variant code
# (items <- str_split(dirname(agb_l0_fp), '/') )
# g0_variant <- nth(items[[1]], -1)

# Apply value cap to AGB -------------------------------------------------------
# Apply saturation point to AGB and save
if(!is.na(saturation_pt)) {
  if(!file.exists(agb_cap_fp)){
    agb_dtype <- raster::dataType(raster::raster(agb_l0_fp))
    agb_sat <- terra::rast(agb_l0_fp) %>% 
      terra::classify(rbind(c(saturation_pt, Inf, saturation_pt),
                            c(-Inf, 0, 0)), 
                      filename = agb_cap_fp, 
                      overwrite=T, 
                      wopt = list(datatype = agb_dtype, gdal='COMPRESS=LZW'))
  }
}

# Mask agb --------------------------------------------------------------
# Apply masks to L2 capped AGB
if(!is.na(saturation_pt)) {
  agb_l2 <- mask_with_options(in_fp = agb_cap_fp, 
                              masks_dir, 
                              agb_fp,
                              masks = agb_mask_codes, 
                              overwrite = FALSE, 
                              agb_dir)
} else {
  # Apply masks to L0 capped AGB
  agb_l2 <- mask_with_options(in_fp = agb_l0_fp, 
                              masks_dir, 
                              agb_fp,
                              masks = agb_mask_codes, 
                              overwrite = FALSE, 
                              agb_dir)
}

# Histogram ----
p_hist <- plot_agb_hist_density(agb_fps$internal, agb_pal)
hist_fp <- file.path(agb_var_dir, str_c('agb_', agb_code, '_hist.png'))
dir.create(dirname(hist_fp), recursive = TRUE)
ggsave(hist_fp, p_hist, width = 5, height = 3, dpi = 120)

# Report (incomplete update): Get value frequencies within masks ---------------
get_frequencies <- function(agb_ras, msk_name){

  # Create reclass matrix
  rcl_mat <- rbind(c(-Inf, 20, -1),
                   c(20, 132, 0),
                   c(132, 200, 1), 
                   c(200, 250, 2), 
                   c(250, 300, 3), 
                   c(300, 310, 4), 
                   c(310, Inf, 5)) 
  colnames(rcl_mat) <- c('V1', 'V2', 'value')
  
  # Convert to tibble and create strings
  rcl_tbl <- rcl_mat %>% 
    as_tibble() %>% 
    mutate(name = str_c(V1, ' to ', V2))
  
  # Reclass masked raster
  agb_rcl <- agb_ras %>% 
    terra::classify(rcl_mat)
  
  # Get class frequencies
  cts <- agb_rcl %>% 
    freq(bylayer = FALSE) %>% 
    as_tibble() %>% 
    left_join(select(rcl_tbl, value, name)) %>% 
    mutate(mask = msk_name)
  
  # Return
  return(cts)
}

# Load
agb_ras <- terra::rast(agb_l0_fp)
msk_land <- terra::rast(landmask_fp)
agb_masked <- agb_l2

# Get name of mask
# g0_variant
mask_name <- agb_masked@ptr$filenames %>% basename() %>% 
  str_extract("(?<=mask).*(?=\\.tif)")
mask_name <- ifelse(mask_name == '', 'none', mask_name)

# Tree cover 
mskinv_T_fp <- file.path(masks_dir, "mask_inv_TreeCover.tif")
mskinv_T <- terra::rast(mskinv_T_fp) %>% 
  terra::crop(agb_ras) %>% 
  terra::mask(msk_land)

# Totals
total_land <- global(!is.na(msk_land), 'sum', na.rm = TRUE) %>% deframe()
total_tc <- global(!is.na(mskinv_T), 'sum', na.rm = TRUE) %>% deframe()
total_mskd <- global(!is.na(agb_masked), 'sum', na.rm = TRUE) %>% deframe()

# Default mask
cts <- get_frequencies(agb_masked, mask_name)
cts <- cts %>% 
  mutate(pct = count / total_mskd,
         pct_land = count / total_land)

# Tree cover, no mask
agb_tc <- agb_ras %>% terra::mask(mskinv_T)
tc_cts <- get_frequencies(agb_tc, 'tree cover unmasked') %>% 
  mutate(pct = count / total_tc,
         pct_land = count / total_land)

# Tree cover, with mask
agb_tc_mskd <- agb_masked %>% terra::mask(mskinv_T)
tc_mskd_cts <- get_frequencies(agb_tc_mskd, str_c('tree cover, ', mask_name)) %>% 
  mutate(pct = count / total_tc,
         pct_land = count / total_land)

# Get mean and count
(tc_mn <- global(agb_tc_mskd, 'mean', na.rm=T))
ct_tc_pixels <- global(!is.na(agb_tc_mskd), 'sum')
(est_tc_area_ha <- (25*25*0.0001) * ct_tc_pixels)
est_tc_area_ha <- 953800
(est_total_AGB <- tc_mn * est_tc_area_ha)

# Get pct of TC class inside water buffer
msk_wb_fp <- file.path(masks_dir, 'mask_water_buff25m.tif')
msk_wb <- terra::rast(msk_wb_fp)
mskinv_T_wb <- mskinv_T %>% terra::mask(msk_wb, inverse = TRUE)
total_tc_wb <- global(!is.na(mskinv_T_wb), 'sum', na.rm = TRUE)
tc_wb_cts <- total_tc_wb %>% 
  rename(count = sum) %>% 
  mutate(name = 'water', 
         mask = 'tree cover', 
         pct = count / total_tc,
         pct_land = count / total_land)

combo_df <- bind_rows(cts, tc_mskd_cts, tc_wb_cts) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 4)))

combo_df %>% write_csv(file.path(agb_dir, str_glue('04_pcts_masked_{agb_code}.csv')))
combo_df %>% write_csv(file.path('data/reports', str_glue('04_pcts_masked_{agb_code}.csv')))

combo <- combo_df %>% 
  mutate(
    name = str_replace_all(name, 
                           "20 to 132|132 to 200|200 to 250|250 to 300|300 to 310|310 to Inf", 
                           'over20')) %>% 
  dplyr::group_by(name, mask) %>% 
  dplyr::summarize(pct = sum(pct), 
                   pct_land = sum(pct_land), 
                   count = sum(count)) %>% 
  mutate(area_ha = round(count * (25*25*0.0001), 0),
         area_from_pct = round(pct_land * 2710675, 0))

combo %>% write_csv(file.path(agb_dir, str_glue('04_pcts_masked_agg_w_areas_{agb_code}.csv')))
combo %>% write_csv(file.path('data/reports', str_glue('04_pcts_masked_agg_w_areas_{agb_code}.csv')))



# Area of Haiti
hti_poly_fp <- "data/tidy/contextual_data/HTI_adm/HTI_adm0_fix.shp"
sf::st_area(st_read(hti_poly_fp))

# Get stats 
# Get mean and count
(mn <- global(agb_masked, 'mean', na.rm=T))
ct_valid_pixels <- global(!is.na(agb_masked), 'sum', na.rm=T)
(est_area_ha <- (25*25*0.0001) * ct_valid_pixels)
(est_total_AGB <- mn * est_area_ha)


