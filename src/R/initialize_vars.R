
# Load libraries ----
# library('BIOMASS')
library("patchwork")
library("gridExtra")
library("sf")
library('terra')
library("tidyverse")

# Codes
# Backscatter pre-processing
code <- 'HV_nu'
year <- '2019'
g0_filt <- 'med5'
g0_mask <- c('L', 'WU', 'wb') # c()
g0_interp <- 'LC' # ''
g0msk_code <- ifelse(length(g0_mask > 0), 
                     str_c('mask', str_c(g0_mask, collapse = ''), '_'), 
                     '')
g0interp_code <- ifelse(g0_interp == 'LC', '_LCinterp', '')
g0_variant <- str_c(g0msk_code, g0_filt, g0interp_code)
if(is.na(g0_variant) | g0_variant == '') g0_variant <- 'simple'

# AGB processing
saturation_pt <- '99.9 percentile'
# saturation_pt <- 300
agb_mask_codes <- c('WU', 'wb') # c('L', 'WU', 'U', 'wb', 'u20', 'o400')
sat_code <- ifelse(is.na(saturation_pt), '', 
                   ifelse(is.numeric(saturation_pt), str_c('cap', saturation_pt),
                          ifelse(saturation_pt == '99.9 percentile', 'cap999pctl', 
                                 'UNKNOWN')))
agbmsk_code <- ifelse(length(agb_mask_codes > 0), 
                     str_c('mask', str_c(agb_mask_codes, collapse = '')), 
                     '')
agb_code <- str_c(sat_code, agbmsk_code, sep = '_') %>% str_remove_all('^_')

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
lc_res_fp <- file.path(tidy_lc_dir, "Lemoiner", "Lemoiner_lc17_hti_resCCI.tif")
g0_fp <- file.path(g0_dir, 'mosaic_variants', str_glue("{code}.tif"))

# 03 ----
g0_dir <- file.path(tidy_dir, str_c('palsar_', year))
modeling_dir <- file.path('data/modeling', code)
agb_dir <- file.path(modeling_dir, g0_variant)
# var_order <- c('simple', 'cappt2_conserv13', 'cappt2_conserv13_mean5', 'med5', 'lee11s10', 'maskLU_lee11s10_LCinterp')
agb_l0_fp <- file.path(agb_dir, str_c("agb_l0_", g0_variant, ".tif"))

# 04 ----
masks_dir <- file.path('data', 'tidy', str_c('palsar_', year), 'masks')
landmask_fp <- file.path(masks_dir, 'hti_land_palsar.tif')
agb_cap_fp <- file.path(agb_dir, str_c('agb_', sat_code, '.tif'))

agb_fp <- file.path(agb_dir, str_glue('agb_{agb_code}.tif'))

# 05 ----
# Set variables 
# rmse_cv <- 21.2
# lc_stat <-  'median'
# input_level <- 'l2_mask'
# mask_level <- 'LU'

# Input filepaths
# (agb_capped_fp <- list.files(agb_dir, str_c('agb_', input_level, '.*\\.tif'), full.names = TRUE))
# (agb_masked_fp <- list.files(agb_dir, str_glue('agb_l2.*{mask_level}\\.tif'), 
#                              full.names = TRUE))

lc_fps <- list(haiti = file.path(raw_lc_dir, "Lemoiner/Haiti2017_Clip.tif"), 
               dr = file.path(raw_lc_dir, "Lemoiner/DR_2017_clip.tif"))
lc_fp <- lc_fps$haiti

# # Output filepaths
# agb_by_lc_prefix <- file.path(agb_dir, 'agb_by_landcover', 
#                               str_glue('agb_{input_level}_{mask_level}_{lc_stat}_byLC'))
# lc_pols_agb_fp <- str_c(agb_by_lc_prefix, '.gpkg')
# agb_by_lc_fp <- str_c(agb_by_lc_prefix, '.tif') 
# agb_by_lc_sd_fp <- str_c(agb_by_lc_prefix, '_sd.tif') 
# 
# agb_filled_fp <- file.path(agb_dir, str_glue('agb_l3_fillLC{lc_stat}_{input_level}_{mask_level}.tif'))
# agb_filled_sd_fp <- file.path(agb_dir, str_glue('agb_l3_fillLC{lc_stat}_{input_level}_{mask_level}_sd.tif'))

# 06

# 07 - External map filepaths ----
glob_fp <- file.path(tidy_maps_dir, "GlobBiomass/Glob_agb10_hti.tif")
esa_fp <- file.path(tidy_maps_dir, "ESA_CCI/ESA_agb17_hti.tif")
avit_fp <- file.path(tidy_maps_dir, "Avitabile/Avitabile_AGB_hti.tif")
bacc_fp <- file.path(tidy_maps_dir, "Baccini/Baccini_agb00_hti.tif")
bacc_res_fp <- file.path(tidy_maps_dir, "Baccini/Baccini_agb00_hti_resCCI.tif")

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

mskvals_list <- list(
  tc = list(values = 4, code = 'TC'),
  total = list(values = c(1,2,3,4,5,6), code = 'Total'),
  bgs = list(values = c(3,5,6), code = 'BGS'),
  wubgs = list(values = c(1,2,3,5,6), code = 'WUBGS')
)

agb_var_dir <- file.path(agb_dir, agb_code)
comparison_dir <- file.path(agb_var_dir, 'external_comparison')
plot_ext_csv <- file.path(comparison_dir, str_c('field_plot_means_', agb_code, '.csv'))
ground_compare_csv <- file.path(comparison_dir, str_glue('07_grounddata_comparison.csv'))
compare_by_lc_csv <- file.path(comparison_dir, str_c('07_comparison_by_LC.csv'))
sums_csv <- file.path(comparison_dir, 'by_LC', '07_sums_by_LC.csv')
ext_pcts_csv <- file.path(comparison_dir, 'by_LC', '07_pcts_by_lc.csv')
agg_tbl_csv <- file.path(comparison_dir, '07_aggregated_comparison.csv')

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


mask_with_options <- function(in_fp, masks_dir, masked_fp, 
                              masks = c('L', 'WU', 'wb'), overwrite = TRUE,
                              out_mask_dir) {
  
  # Get output filenames
  masks_str <- masks %>% str_c(collapse = '')
  # masked_fp <- file.path(dirname(in_fp), 
  #                        str_c(code, '_mask', masks_str, '.tif'))
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
    msk_u20 <- ras %>% terra::classify(rbind(c(-Inf, 20, NA), c(20, Inf, 1)))
    msk <- msk %>% terra::mask(msk_u20)
  }
  
  # AGB > X
  oX <- str_subset(masks, 'o\\d+')
  if(length(oX) > 0){

    if(str_detect(oX, '999pctl')){
      # if mask is on 99.9th percentile... 
      p999 <- global(ras, function(x) quantile(x, probs = 0.999, na.rm = TRUE))
      cap_val <- deframe(p999)
      
    } else {
      # if not, use the value specified
      cap_val <- as.numeric(str_extract(oX, '\\d+'))
      
    }
    
    # Apply mask
    msk_oX <- ras %>% terra::classify(rbind(c(-Inf, cap_val, 1), c(cap_val, Inf, NA)))
    msk <- msk %>% terra::mask(msk_oX)
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


crop_to_intersecting_extents <- function(r1, r2, return_r1=T, return_r2=T, return_bb=F) {
  # Get intersection of the bounding boxes of the two rasters
  r1_bbox <- terra::ext(r1) %>% as.vector() %>% sf::st_bbox() %>% sf::st_as_sfc()
  r2_bbox <- terra::ext(r2) %>% as.vector() %>% sf::st_bbox() %>% sf::st_as_sfc()
  bb <- sf::st_intersection(r1_bbox, r2_bbox)
  bb <- bb %>% sf::st_bbox()
  bb <- bb %>% as.vector()
  
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

plot_agb_hist_density <- function(x, agb_pal, bwidth=5, sample_size=100000, 
                                  xlim = c(min=0, max=300)) {
  
  # Get values
  lyr_name <- x$name
  ras <- terra::rast(x$fp)
  df <- terra::values(ras, dataframe = TRUE) %>% 
    rename(value = 1) %>% 
    drop_na()
  
  # Get statistics
  mu <- mean(df$value, na.rm=T)
  sd1 <- sd(df$value, na.rm=T)
  mn <- min(df$value, na.rm=T)
  mx <- max(df$value, na.rm=T)
  
  cuts <- quantile(df$value, probs = c(.1, .5, .9)) %>% as_tibble(rownames='ref')
  cuts_pts <- data.frame(x = c(mn, mx), 
                         y = c(0, 0))
  
  # Sample
  if (nrow(df) > sample_size) {
    print("Sampling...")
    df <- df %>% sample_n(sample_size)
  }
  
  # Plot
  p <- ggplot(df, aes(x=value)) + 
    geom_histogram(aes(y=..density.., fill = ..x..), binwidth = bwidth,
                   show.legend = FALSE)+
    scale_fill_gradientn(colors = agb_pal$colors, 
                         # breaks = c(-60, -30, 0, 30, 60),
                         # labels = c('-60 (internal < external)', '', '0', '', '60 (internal > external)'),
                         name = expression(paste("AGB (Mg ha"^"-1", ")")),
                         # na.value='lightgray', 
                         limits = c(agb_pal$min, agb_pal$max),
                         oob = scales::squish,
                         guide = guide_colorbar(barheight = 2)
    ) + 
    geom_density(alpha=.2,
                 show.legend = FALSE) + 
    scale_x_continuous(name = expression(paste("Aboveground biomass (Mg ha"^"-1", ")")),
                       # breaks = seq(min, max, 100),
                       # limits = c(min, max)
    ) +
    coord_cartesian(xlim = xlim) +
    geom_vline(mapping = aes(xintercept = value,
                             colour = ref),
               data = cuts,
               color="black", 
               # linetype="solid", 
               linetype="dashed", 
               size=.5,
               show.legend = FALSE) +
    geom_text(mapping = aes(x = value,
                            y = Inf,
                            label = ref,
                            hjust = -.1,
                            vjust = 1),
              data = cuts) +
    ggtitle(lyr_name) +
    theme_minimal() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(), 
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank())
}

# 07 functions ---

summarize_raster_differences <- function(ext_r, ext_name, int_r, out_fp = NULL){
  
  # Make DF
  r_stack <- rast(list(internal = int_r, external = ext_r))
  names(r_stack) <- c('internal', 'external')
  df <- terra::as.data.frame(r_stack, na.rm = TRUE)
  df <- df %>% mutate(value = internal - external)
  
  # Distribution of external values
  dist_ext <- summary(df$external) %>% 
    tibble(names=names(.), value=.) %>% 
    mutate(value = as.numeric(value))
  
  # Get metrics for differences
  diff_stats <- df %>% 
    # group_by(name) %>% 
    summarize(
      R2 = summary(lm(external ~ internal))$r.squared,
      adjR2 = summary(lm(external ~ internal))$adj.r.squared,
      N = sum(!is.na(value)), 
      RMSD = sqrt(mean(value^2, na.rm=T)),
      bias = abs(sum(value, na.rm=T) / sum(!is.na(value))),
      MBD = mean(value), 
      IQR = IQR(value), 
      SD = sd(value), 
      MAD = mean(abs(value)),
      cor.coef = cor.test(x=internal, y=external, method = 'pearson')$estimate,
      min_diff = min(value),
      max_diff = max(value),
      med_diff = median(value)
    ) %>% 
    mutate(across(where(is.numeric), ~ round(.x, 3)),
           name = ext_name)
  
  # Save
  if(is.character(out_fp)){
    diff_stats %>% write_csv(out_fp)
  }
  
  # Return
  return(list(diffs=df, summary=diff_stats))
}


scatter_differences <- function(df, ext_name = '', prefix, filename = NA,
                                sample_size = 5000000, xlim = c(0, 300)){
  
  # Sample
  if (length(df[[1]]) > sample_size) {
    print("Sampling...")
    df <- df %>% sample_n(sample_size)
  }
  
  # Calculate density surface outside of ggplot object
  f1 <- MASS::kde2d(df$internal, df$external, 
                    n = 300, lims = c(xlim, xlim))

  # Convert to rast and then to df
  r <- rast(f1$z) %>% flip(direction = 'vertical')
  ext(r) <- c(xmin = xlim[1], xmax= xlim[2],
              ymin = xlim[1], ymax= xlim[2])
  dfr <- as.data.frame(r, xy=TRUE)
  
  # Plot
  p <- ggplot(dfr) +
    geom_raster(aes(x=x, y=y, fill=lyr.1)) +
    scale_fill_gradientn(colors = c('white', '#F8735C', 
                                    '#772181', '#00000A'),
                         name = 'Density', 
                         limits = c(0, 0.0004),
                         expand = expansion()
                         ) +
    scale_y_continuous(name = str_c("External AGB estimate (t/ha)"),
                       breaks = seq(xlim[[1]], xlim[[2]], 100), 
                       limits = xlim,
                       expand = expansion()
    ) +
    scale_x_continuous(name = str_c("Internal AGB estimate (t/ha)"),
                       breaks = seq(xlim[[1]], xlim[[2]], 100), 
                       limits = xlim,
                       expand = expansion()
    ) +
    coord_fixed(ratio = 1, xlim = xlim) +
    theme_minimal() + 
    geom_abline(intercept = 0, slope = 1, col='black', size=.25) +
    # geom_abline(intercept = coef(ols)[1], slope = coef(ols)[2], 
    #             col='black', size=.25, linetype = 'dashed') +
    ggtitle(ext_name)
  
  # ols <- lm(external ~ internal, data=df)
  # 
  # # Table with slope and intercept
  # rep1 <- coef(ols) %>%
  #   as_tibble() %>% 
  #   mutate(value = round(value, 1), 
  #          name = c('intercept', 'slope')) %>% 
  #   column_to_rownames('name')
  # 
  # tab <- gridExtra::tableGrob(rep1, 
  #                             cols = NULL,
  #                             theme = gridExtra::ttheme_minimal(base_size = 9,
  #                                                               padding = unit(c(2,2), 'mm')))
  
  # # Overlay plot with regression table
  # p <- p + inset_element(tab, 
  #                        left = 0, bottom = 0.8,
  #                        right = 0.2, top = 1, 
  #                        on_top = TRUE) +
  #   theme(plot.background = NULL)
  
  # Save/Return
  if(!is.na(filename)) ggsave(filename, plot=p, width=6, height=6)
  
  # Return
  return(p)
  
}

compare_by_given_lc <- function(agb_x, agb_fp, lc_fp,
                                mskvals = list(values = 4, code = 'TC'), 
                                comparison_dir,
                                return_obj = 'summary', 
                                boundary_fp) {
  
  # External names
  ext_name <- agb_x$name
  ext_fp <- agb_x$fp
  ext_masked_fp <- str_c(tools::file_path_sans_ext(ext_fp), '_', mskvals$code, 'mask.tif')
  
  if (ext_name == 'Avitabile') {    
    agb_res_fp <- str_c(tools::file_path_sans_ext(agb_fp), '_res', ext_name, '.tif')
    lc_res_fp <- str_c(tools::file_path_sans_ext(lc_fp), '_res', ext_name, '.tif')
    
    # Resample LC and internal to external AGB res
    if(!file.exists(agb_res_fp)) resample_to_raster(agb_fp, ext_fp, agb_res_fp)
    if(!file.exists(lc_res_fp)) resample_to_raster(lc_fp, ext_fp, lc_res_fp, method = 'near')
    
  } else {
    agb_res_fp <- str_c(tools::file_path_sans_ext(agb_fp), '_resCCI.tif')
    lc_res_fp <- file.path(tidy_lc_dir, "Lemoiner", "Lemoiner_lc17_hti_resCCI.tif")
    
  }
  
  # Internal names
  agb_code <- str_remove(tools::file_path_sans_ext(agb_res_fp), '.*agb_')
  int_masked_fp <- str_c(tools::file_path_sans_ext(agb_res_fp), '_', mskvals$code, 'mask.tif')
  
  # Mask internal AGB to given LC class
  if (!file.exists(int_masked_fp)) {
    mask_to_classes(lc_res_fp, agb_res_fp, int_masked_fp, mskvals$values)
  } 
  
  # Mask external AGB to given LC class
  if (!file.exists(ext_masked_fp)) {
    mask_to_classes(lc_res_fp, ext_fp, ext_masked_fp, mskvals$values)
  }
  
  # Calculate metrics for landcover class
  out <- crop_to_intersecting_extents(terra::rast(int_masked_fp), 
                                      terra::rast(ext_masked_fp))
  agb_T <- out$r1
  ext_T <- out$r2
  
  if (return_obj == 'map') {
    # Subtract
    diff_ras <- agb_T - ext_T
    
    # Plot spatial map of differences 
    map_diff_fp <- file.path(comparison_dir, str_c('map_', agb_code, '_minus_', ext_name, '.png'))
    p_map <- plot_differences_map(diff_ras, ext_name, agb_code, boundary_fp, filename = map_diff_fp)
    
  }
  
  # Summary of differences ----
  smmry_diff_fp <- file.path(comparison_dir, 'by_LC', 'temp',
                             str_c(mskvals$code, '_summary_diff_', ext_name, '_v_', agb_code, '.csv'))
  dir.create(dirname(smmry_diff_fp), recursive = TRUE, showWarnings = FALSE) 
  df <- summarize_raster_differences(ext_T, ext_name, agb_T, out_fp=smmry_diff_fp)
  
  # Scatterplot ----
  if (return_obj == 'plot' | return_obj == 'all') {
    scatter_diff_fp <- file.path(comparison_dir, 'by_LC', 
                                 str_c(mskvals$code, '_scatt_', agb_code, '_v_', ext_name, '.png'))
    p_scatt <- scatter_differences(df$diffs, ext_name, agb_code, 
                                   filename = scatter_diff_fp,
                                   sample_size = 500000)
    
    # Overlay comparison report on plot
    metrics <- df$summary %>% 
      select(N, R2, RMSD, MBD) %>% 
      pivot_longer(everything())
    
    # Table with slope and intercept
    ols <- lm(external ~ internal, data=df$diffs)
    ols_met <- coef(ols) %>%
      as_tibble() %>% 
      mutate(name = c('intercept', 'slope'))
    
    # Create formatted strings
    rep1 <- bind_rows(metrics, ols_met) %>%
        mutate(
          valstr = signif(value, 3) %>% 
            format(big.mark = ',',
                   trim = TRUE,
                   digits = 2, 
                   scientific = FALSE,
                   drop0trailing = TRUE),
          string = str_c(name, ' = ', valstr))
    tab1 <- select(rep1, string) %>% deframe %>% str_c(collapse = '\n')
    grob <- grid::grobTree(grid::textGrob(tab1, x=0.05,  y=0.85, hjust=0),
                           gp=grid::gpar(col="black", fontsize=9))
    
    # Place on plot
    p_scatt <- p_scatt + 
      annotation_custom(grob)+
      geom_abline(intercept = coef(ols)[1], slope = coef(ols)[2], 
                  col='black', size=.25, linetype = 'dashed')
  }
  
  # Return ----
  out_list <- list()
  if ('plot' %in% return_obj) out_list <- c(out_list, list(plot = p_scatt))
  if ('diffs' %in% return_obj) out_list <- c(out_list, list(diffs = df$diffs))
  if ('summary' %in% return_obj) out_list <- c(out_list, list(summary = df$summary))
  if ('map' %in% return_obj) out_list <- c(out_list, list(map = p_map))
  return(out_list)
  
  # # Return ----
  # if (return_obj == 'plot') {
  #   return(p_scatt)
  #   
  # } else if (return_obj == 'diffs') {
  #   return(df$diffs)
  #   
  # } else if (return_obj == 'summary') {
  #   return(df$summary)
  #   
  # } else if (return_obj == 'map') {
  #   return(p_map)
  #   
  # } else if (return_obj == 'all') {
  #   return(list(plot = p_scatt,
  #               diffs = df$diffs,
  #               summary = df$summary,
  #               map = p_map))
  #   
  # } else {
  #   return()
  #   
  # }
}
