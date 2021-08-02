# ******************************************************************************
# Script to:
#     * Compare output AGB map to other biomass estimates
# Proceeds:
#     * post_fill_AGB_gaps.R - performs post-processing such as masking
#     * postprocess_AGB_map.R - performs post-processing such as masking
#     * regression_AGB-g0.R - creates AGB map
#     * calculate_AGB.R - calculates AGB by plot from field data
#     * process_ALOS_tiles.R
# Requires:
#     * AGB map (agbr)
#     * external AGB maps
# ******************************************************************************

# Load libraries ----
# library(raster)
# library(stars)
# library(tmap)
# library(units)
# library(tools)
# library(broom)
# library(geobgu)
# # library(gridExtra)
# library(patchwork)
# library(clipr)
# library('gdalUtils')
library('sf')
library('terra')
library('tidyverse')
library('patchwork')

source('src/R/initialize_vars.R')
lc_fp <- lc_fps$haiti

# Set variables ----------------------------------------------------------------
# Input filepaths
# (agb_fps <- list.files(agb_dir, str_c('agb_', agb_input_level, '.*[^(sd)]\\.tif'),
#                       full.names = TRUE))
# (agb_fp <- agb_fps[[2]])
# (agb_fp <- file.path(agb_dir, str_glue('agb_{agb_code}.tif')))

# FUNCTIONS -----------------------------------------------------------------------------------
crop_and_resample <- function(in_fp, out_fp, template, warp_method="near", na_value=NA, overwrite=F){
  
  fun <- function(in_fp, out_fp, template, na_value=NA){
    # Resample raster to template raster with crop from GDAL for speed. 
    
    if(class(template)=='character'){
      template <- read_stars(template)
    }
    
    # Crop using GDAL
    print('Cropping...')
    fp_crop <- tempfile(pattern = str_glue(file_path_sans_ext(basename(in_fp)),'_crop'), 
                        tmpdir = tempdir(), 
                        fileext = "tif")
    # fp_crop <- str_glue(file_path_sans_ext(basename(in_fp)),'_crop.',file_ext(in_fp))
    r <- raster(in_fp)
    gdalwarp(srcfile=in_fp, dstfile=fp_crop,
             te=st_bbox(template), tr=c(xres(r), yres(r)), tap=T, overwrite=T)
    
    # Load and resample
    print('Resampling...')
    r <- read_stars(fp_crop)
    
    if(!is.na(na_value)){
      r[r == na_value] <- NA
    }
    
    r <- r %>% st_warp(template, method=warp_method, use_gdal=TRUE)
    r %>% as("Raster") %>%
      writeRaster(out_fp, options=c("dstnodata=-99999"), overwrite=T)
    return(r) # output raster is stars format
  }
  
  # Only run if file doesn't already exist
  if(!file.exists(out_fp) | overwrite){
    
    r <- fun(in_fp, out_fp, template, na_value)
    
  } else {
    
    # If the file exists, load it
    print("File already exists. Loading.")
    r <- read_stars(out_fp)
    
  }
}

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

plot_hist_density <- function(df, min=-200, max=200, bwidth=50, sample_size=10000, 
                              ext_name, prefix, filename = NA){
  
  # Get statistics
  mu <- mean(df$value, na.rm=T)
  sd1 <- sd(df$value, na.rm=T)
  mn <- min(df$value, na.rm=T)
  mx <- max(df$value, na.rm=T)

  cuts <- quantile(df$value, probs = c(.1, .5, .9)) %>% as_tibble(rownames='ref')
  cuts_pts <- data.frame(x = c(mn, mx), 
                         y = c(0, 0))
  
  if (length(df[[1]]) > 1000000) {
    print("Sampling...")
    df <- df %>% sample_n(sample_size)
  }
  
  # Plot
  p <- ggplot(df, aes(x=value)) + 
    # geom_histogram(aes(y=..density..), fill="#69b3a2", color="#e9ecef", alpha=0.7, binwidth = bwidth)+
    geom_histogram(aes(y=..density.., fill = ..x..), binwidth = bwidth,
                              show.legend = FALSE)+
    scale_fill_gradientn(colors=c('#ca0020', '#f4a582', '#f7f7f7', '#92c5de', '#0571b0'), 
                         breaks = c(-60, -30, 0, 30, 60),
                         labels = c('-60 (internal < external)', '', '0', '', '60 (internal > external)'),
                         name = 'Difference in AGB (t/ha)',
                         # na.value='lightgray', 
                         limits = c(-60, 60), 
                         oob = scales::squish,
                         guide = guide_colorbar(barheight = 2)) + 
    geom_density(alpha=.2) + 
    scale_x_continuous(name = str_c("Difference (t/ha): our AGB (", prefix, ", Haiti) - ", ext_name),
                       breaks = seq(min, max, 100), 
                       limits = c(min, max)
    ) +
    coord_cartesian(xlim=c(min, max)) +
    geom_point(aes(x=mn, y=0))+
    geom_point(aes(x=mx, y=0))+
    geom_vline(mapping = aes(xintercept = value,
                             colour = ref),
               data = cuts,
               color="black", 
               # linetype="solid", 
               linetype="dashed", 
               size=.5,
               show.legend = FALSE) +
    # geom_vline(aes(xintercept=mu - sd1),
    #            color="black", linetype="dashed", size=.5)+
    # geom_vline(aes(xintercept=mu + sd1),
    #            color="black", linetype="dashed", size=.5)+
    geom_text(mapping = aes(x = value,
                            y = Inf,
                            label = ref,
                            hjust = -.1,
                            vjust = 1),
              data = cuts) +
    theme_minimal()
  
  # Save figure
  if(!is.na(filename)) {
    
    ggsave(filename, plot=p, width=8, height=4)
    
  }
  
  # Return   
  return(p)
}

scatter_differences <- function(df, ext_name, prefix, filename = NA,
                                sample_size = 5000000, xlim = c(0, 300)){
  
  # Sample
  if (length(df[[1]]) > sample_size) {
    print("Sampling...")
    df <- df %>% sample_n(sample_size)
  }
  
  ols <- lm(external ~ internal, data=df)

  # Table with slope and intercept
  rep1 <- coef(ols) %>%
    as_tibble() %>% 
    mutate(value = round(value, 2), 
           name = c('intercept', 'slope')) %>% 
    column_to_rownames('name')
  
  tab <- gridExtra::tableGrob(rep1, 
                              cols = NULL,
                              theme = gridExtra::ttheme_minimal(base_size = 9,
                                                                padding = unit(c(2,2), 'mm')))
  
  # Plot
  p <- ggplot(df, aes(x=internal, y=external)) + 
    # geom_point() +
    stat_density2d(geom="tile",
                   aes(fill = after_stat(density),
                       alpha = after_stat(density)^0.1),
                   contour_var = 'count',
                   contour=FALSE) +
    # geom_bin2d(bins = 75) +
    # geom_hex(bins = 150) +
    scale_fill_viridis_c(option = 'magma', direction = -1) +
    scale_alpha('alpha') +
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
    geom_abline(intercept = coef(ols)[1], slope = coef(ols)[2], 
                col='black', size=.25, linetype = 'dashed') +
    ggtitle(ext_name) +
    theme_minimal()

  # Overlay plot with regression table
  p <- p + inset_element(tab, 
                  left = 0, bottom = 0.8,
                  right = 0.2, top = 1, 
                  on_top = TRUE) +
    theme(plot.background = NULL)
  
  # Save/Return
  if(!is.na(filename)) ggsave(filename, plot=p, width=6, height=6)
 
  # Return
  return(p)

}

plot_differences_map <- function(diff_ras, ext_name, prefix, boundary_fp, filename = NA){
  
  # Convert raster to dataframe
  diff_df <- as.data.frame(diff_ras, xy=TRUE) %>% 
    drop_na() %>% 
    rename(diff = all_of(names(diff_ras)))
  
  # Load land polygon
  hti_poly <- st_read(boundary_fp)
  
  # Plot
  p <- ggplot(diff_df) +
    # Land poly
    geom_sf(data=hti_poly, lwd=0, fill='darkgray') +
    # Difference surface
    geom_tile(aes(x=x, y=y, fill=diff)) +
    scale_fill_gradientn(colors=c('#ca0020', '#f4a582', '#f7f7f7', '#92c5de', '#0571b0'), 
                         breaks = c(-60, -30, 0, 30, 60),
                         labels = c('-60 (internal < external)', '', '0', '', '60 (internal > external)'),
                         name = 'Difference in AGB (t/ha)',
                         # na.value='lightgray', 
                         limits = c(-60, 60), 
                         oob = scales::squish,
                         guide = guide_colorbar(barheight = 2)) + 
    coord_sf() +
    theme_minimal() + 
    theme(axis.title = element_blank(), 
          axis.text = element_text(color='gray'),
          legend.position = c(.16, .85))
  
  # Save figure
  if(!is.na(filename)) {
    
    ggsave(filename, plot=p, width=8, height=6)

  } 
  
  # Return
  return(p)
  
}

get_summary_wrapper <- function(agb_fp, ext_fp, ext_name, agb_code) {
  
  # Filepaths
  agb_res_fp <- str_c(tools::file_path_sans_ext(agb_fp), '_res', ext_name, '.tif')
  
  # Load and crop
  agb_res <- terra::rast(agb_res_fp)
  agb_ext <- crop_to_intersecting_extents(agb_res, terra::rast(ext_fp), return_r1=F)
  
  # Summarize differences
  out <- summarize_raster_differences(agb_ext, ext_name, agb_res)
  
}

# Extract AGB means at survey sites
extract_mean <- function(x, vec, touches = FALSE, method='simple'){
  
  # Load data in terra format
  lyr_name <- x$name
  ras <- terra::rast(x$fp)
  vec <- terra::vect(vec)
  
  # Convert to SpatVector and extract AGB values
  ex1 <- ras %>% terra::extract(vec, touches = touches, method = method) 
  ex2 <- ex1 %>% as_tibble()
  ex3 <- ex2 %>% rename(value = 2)
  
  # Convert matrix to tibble and get mean and count for each polygon
  ex4 <- ex3 %>% 
    group_by(ID) %>% 
    summarise(value = mean(value, na.rm=T) %>% round(2)) %>% 
    arrange(ID) %>% 
    select(value) %>% 
    deframe()
}

get_pcts <- function(x, hti_poly_fp) {
  
  # Load 
  if(class(x) == 'SpatRaster') {
    ras <- x
  } else if(is.list(x)) {
    ras <- terra::rast(x$fp) 
  } else if(is.character(x)) {
    ras <- terra::rast(x)
  } 

  msk_poly <- terra::vect(hti_poly_fp)
  
  # Reclassify
  ras_clas <- ras %>% 
    terra::classify(rbind(c(-Inf, Inf, 1),
                          c(NA, NA, 0)))
  
  # Crop
  ras_msk <- ras_clas %>% 
    terra::crop(msk_poly) %>% 
    terra::mask(msk_poly, 
                inverse = FALSE)
  
  # Get values
  rvals <- terra::values(ras_msk)
  
  # Count pixels in each category and convert to percentage
  df <- tibble(group = c("value", "NA"), 
               count = c(sum(rvals == 1, na.rm = TRUE), 
                         sum(rvals == 0, na.rm = TRUE)))
  
  # Convert counts to percentage of all land pixels
  landct <- sum(df$count)
  df <- df %>%
    dplyr::add_row(group = 'land', count = landct) %>% 
    dplyr::mutate(pct = count / landct)
  
  # Structure data
  df %>% 
    pivot_wider(names_from = group, values_from = count:pct) %>% 
    select(count_land, pct_value, pct_NA)
}

mask_to_classes <- function(mask_fp, ras_fp, out_fp, class_values) {
  # Load and pre-process layers
  out <- crop_to_intersecting_extents(terra::rast(mask_fp), terra::rast(ras_fp))
  lc <- out$r1
  agb_ext <- out$r2
  
  # Mask by each class in list and sum
  ext_T <- class_values %>% 
    purrr::map(~ mask(agb_ext, lc, inverse = TRUE, maskvalues = .x)) %>% 
    rast() %>% 
    app(sum, na.rm = TRUE, filename = out_fp, overwrite = TRUE)
}

get_pcts_masked <- function(agb_x, mskvals, lc_res_fp, hti_poly_fp) {
  # External names
  ext_name <- agb_x$name
  ext_fp <- agb_x$fp
  ext_masked_fp <- str_c(tools::file_path_sans_ext(ext_fp), '_', mskvals$code, 'mask.tif')
  
  if (!file.exists(ext_masked_fp)) {
    mask_to_classes(lc_res_fp, ext_fp, ext_masked_fp, mskvals$values)
  }
  
  df <- get_pcts(ext_masked_fp, hti_poly_fp)
  
  df %>% mutate(name = ext_name, 
                mask = mskvals$code)
  
  
  

  # Reclassify
  ras_clas <- ras %>% 
    terra::classify(rbind(c(-Inf, Inf, 1),
                          c(NA, NA, 0)))
  
  # Crop
  ras_msk <- ras_clas %>% 
    terra::crop(msk_poly) %>% 
    terra::mask(msk_poly, 
                inverse = FALSE)
  
  # Get values
  rvals <- terra::values(ras_msk)
  
  # Count pixels in each category and convert to percentage
  df <- tibble(group = c("value", "NA"), 
               count = c(sum(rvals == 1, na.rm = TRUE), 
                         sum(rvals == 0, na.rm = TRUE)))
  
  # Convert counts to percentage of all land pixels
  landct <- sum(df$count)
  df <- df %>%
    dplyr::add_row(group = 'land', count = landct) %>% 
    dplyr::mutate(pct = count / landct)
  
  # Structure data
  df %>% 
    pivot_wider(names_from = group, values_from = count:pct) %>% 
    select(count_land, pct_value, pct_NA)
  
}

get_pcts_by_lc <- function(agb_x, lc_res_fp) { 
  
  # External names
  ext_name <- agb_x$name
  ext_fp <- agb_x$fp
  
  if (ext_name == 'Avitabile') {    
    lc_res_fp <- str_c(tools::file_path_sans_ext(lc_fp), '_res', ext_name, '.tif')
    
    # Resample LC and internal to external AGB res
    if(!file.exists(lc_res_fp)) resample_to_raster(lc_fp, ext_fp, lc_res_fp, method = 'near')
    
  } 
  
  # Load 
  out <- crop_to_intersecting_extents(terra::rast(lc_res_fp), terra::rast(ext_fp))
  lc <- out$r1; ras <- out$r2
  
  # Reclassify masked raster to mask (1s and 0s)
  ras_clas <- ras %>% 
    terra::classify(rbind(c(-Inf, Inf, 1),
                          c(NA, NA, 0)))
  
  # Multiply
  ras_lc <- ras_clas * lc
  
  # Get count in each class
  lc_name = c("Null", "Water", "Urban", "Bareland", "Tree cover", "Grassland", "Shrubs")
  lc_names <- tibble(code = 0:6, lc_name = lc_name)
  
  # Count values present in each LC
  cnts_vals <- as.data.frame(ras_lc) %>% 
    rename(code = 1) %>% 
    count(code) %>% 
    left_join(lc_names) %>% 
    filter(code != 0) %>%
    mutate(pct_of_values = n / sum(n))
  
  # Count pixels in each LC
  cnts_lc <- as.data.frame(lc) %>% 
    rename(code = 1) %>% 
    count(code) %>% 
    left_join(lc_names)
  
  # Join to total LC counts and get percentages
  cnts_vals %>% 
    left_join(cnts_lc, by = c('code', 'lc_name'), suffix = c('_vals', '_lc')) %>% 
    mutate(pct_of_lc = n_vals / n_lc,
           name = ext_name)

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
    # diff_ras %>% terra::writeRaster(filename = diff_fp, overwrite = TRUE, datatype = in_dtype)
    
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
  }
  
  # Return ----
  if (return_obj == 'plot') {
    return(p_scatt)
    
  } else if (return_obj == 'diffs') {
    return(df$diffs)
    
  } else if (return_obj == 'summary') {
    return(df$summary)
    
  } else if (return_obj == 'map') {
    return(p_map)
    
  } else if (return_obj == 'all') {
    return(list(plot = p_scatt,
                diffs = df$diffs,
                summary = df$summary,
                map = p_map))
    
  } else {
    return()
    
  }
}

iterate_compare_by_lc <- function(mskvals, agb_fps, lc_fp, comparison_dir, 
                                  plot = TRUE,
                                  boundary_fp){
  
  agb_fp <- agb_fps$internal$fp
  agb_fps_sub <- agb_fps[c(2,3,4,5)]
  
  # Get summary table for all externals and row bind
  sum_tbl <- agb_fps_sub %>% 
    purrr::map_dfr(compare_by_given_lc, 
                   agb_fp, lc_fp,
                   mskvals = mskvals, 
                   comparison_dir,
                   return_obj = 'summary', 
                   boundary_fp = boundary_fp) %>% 
    mutate(mask = mskvals$code) %>% 
    mutate(map_code = names(agb_fps_sub))
  
  # Save table
  smmry_diff_fp <- file.path(comparison_dir, 'by_LC', 'temp',
                             str_c('summary_diffs_', mskvals$code, '.csv'))
  sum_tbl %>% write_csv(smmry_diff_fp)
  
  # Plots
  if (plot) {
    agb_fps_sub %>% 
      purrr::walk(compare_by_given_lc, agb_fp, lc_fp,
                  mskvals = mskvals, comparison_dir,
                  return_obj = 'plot', 
                  boundary_fp = hti_poly_fp)
  }
  
  return(sum_tbl)
}

# Compare to field data: r2, RMSE, bias for external maps -----------------------------
# Get plot backscatter mean to polygons
plots_agb <- readRDS(field_agb_fp) %>% arrange(plot_no)

# Extract plot means for all AGB maps
plots2 <- agb_fps %>% 
  purrr::map_dfc(extract_mean, plots_agb, TRUE) %>% 
  mutate(across(everything(), ~ replace_na(.x, 0)))

plot_means <- plots_agb %>% 
  st_drop_geometry() %>% 
  select(plot_no, AGB_ha) %>% 
  mutate(AGB_ha = round(AGB_ha, 2)) %>% 
  bind_cols(plots2) 

# Save
plot_ext_csv <- file.path(comparison_dir,
                     str_c('field_plot_means_', agb_code, '.csv'))
dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)
plot_means %>% write_csv(plot_ext_csv)

# Get differences
df <- plot_means %>% 
  pivot_longer(cols = internal:bacc, names_to = 'map_code', values_to = 'est') %>% 
  mutate(diff = est - AGB_ha) 

ground_metrics <- df %>% 
  group_by(map_code) %>% 
  summarize(
    R2 = summary(lm(est ~ AGB_ha))$r.squared,
    adjR2 = summary(lm(est ~ AGB_ha))$adj.r.squared,
    N = sum(!is.na(diff)), 
    RMSD = sqrt(mean(diff^2, na.rm=T)),
    bias = abs(sum(diff, na.rm=T) / sum(!is.na(diff))),
    MBD = mean(diff)
    ) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 3)))

# Save as CSV
ground_metrics %>% write_csv(ground_compare_csv)

# Pixel-to-pixel comparison ----
# testing ----
# mskvals <-  list(values = 4, code = 'TC')
# agb_x <- agb_fps$avit
# df <- compare_by_given_lc(agb_x, agb_fp, lc_fp = lc_fps$haiti,
#                     mskvals = mskvals, comparison_dir,
#                     return_obj = 'all', 
#                     boundary_fp = hti_poly_fp)
# 
# agb_x <- agb_fps$esa
# compare_by_given_lc(agb_x, agb_res_fp, lc_fp = lc_res_fp,
#                           mskvals = mskvals, comparison_dir,
#                           return_obj = 'plot', 
#                     boundary_fp = hti_poly_fp)
# 
# p_scatt <- scatter_differences(df, ext_name, agb_code, 
#                                sample_size = 500000)
# 
# # TC
# mskvals <- list(values = 4, code = 'TC')
# tc_sum_tbl <- agb_fps[c(2,3,4,5)] %>% 
#   purrr::map_dfr(compare_by_given_lc, agb_fp, lc_fp,
#                  mskvals = mskvals, comparison_dir,
#                  return_obj = 'summary', 
#                  boundary_fp = hti_poly_fp) %>% 
#   mutate(mask = mskvals$code) 
# 
# smmry_diff_fp <- file.path(comparison_dir, 'by_LC', 
#                            str_c('summary_diffs_', mskvals$code, '.csv'))
# tc_sum_tbl %>% write_csv(smmry_diff_fp)
# 
# # Plots
# agb_fps[c(2,3,4,5)] %>% 
#   purrr::walk(compare_by_given_lc, agb_fp, lc_fp,
#               mskvals = mskvals, comparison_dir,
#               return_obj = 'plot', 
#               boundary_fp = hti_poly_fp)
# 
# # Difference map
# mskvals <-  list(values = 4, code = 'TC')
# agb_x <- agb_fps$avit
# p_map <- compare_by_given_lc(agb_x, agb_fp, lc_fp = lc_fp,
#                           mskvals = mskvals, comparison_dir,
#                           return_obj = 'map', 
#                           boundary_fp = hti_poly_fp)

# ~ Summary metric tables and scatter density plots ----
mskvals_list <- list(
  list(values = 4, code = 'TC'),
  list(values = c(1,2,3,4,5,6), code = 'Total'),
  list(values = c(3,5,6), code = 'BGS'),
  list(values = c(1,2,3,5,6), code = 'WUBGS')
)

# Parallelize
# library('furrr')
# options(future.fork.enable = TRUE)
# plan(multicore, workers = 7)
# 
# # Fails with plot = TRUE
# sum_tbl <- mskvals_list %>%
#   furrr::future_map(iterate_compare_by_lc, agb_fps, lc_fp, comparison_dir, 
#                     plot = FALSE,
#                  boundary_fp = hti_poly_fp, 
#                  .options = furrr_options(seed = TRUE)) %>% 
#   bind_rows()

# Create plots 
sum_tbl <- mskvals_list %>%
  purrr::map_dfr(iterate_compare_by_lc, agb_fps, lc_fp, comparison_dir,
                 boundary_fp = hti_poly_fp)

# Combine
sum_tbl <- list.files(file.path(comparison_dir, 'by_LC', 'temp'),
           pattern = 'summary_diffs.*\\.csv', 
           full.names = TRUE) %>% 
  purrr::map_dfr(read_csv)

# Save
pix2pix_csv <- file.path(comparison_dir, str_c('07_pixel_comparison_by_LC.csv'))
sum_tbl %>% write_csv(pix2pix_csv)

# Get masked percentages by LC ----
res_fps <- agb_fps
res_fps$internal <- list(
  name = str_glue('Internal {agb_code}, 100m'),
  fp = str_c(tools::file_path_sans_ext(agb_fp), '_resCCI.tif')
)

map_names <- res_fps %>% map_chr('name') %>% as_tibble(rownames = 'map_code')
ext_pcts <- res_fps %>% 
  purrr::map_dfr(get_pcts_by_lc, lc_res_fp) %>% 
  left_join(map_names, by = c(name = 'value'))

# Save
ext_pcts_csv <- file.path(comparison_dir, str_c('07_pcts_by_lc_', agb_code, '.csv'))
ext_pcts %>% write_csv(ext_pcts_csv)

# Aggregate to mask categories
ext_pcts <- read_csv(ext_pcts_csv)
agg_fun <- function(mskvals, ext_pcts) {
  ext_pcts %>% 
    filter(code %in% mskvals$values) %>% 
    group_by(name, map_code) %>% 
    summarize(across(starts_with(c('n', 'pct')), sum)) %>% 
    mutate(pct_of_lc =  n_vals / n_lc, 
           mask = mskvals$code)
}

pcts_masked <- mskvals_list %>% purrr::map_dfr(agg_fun, ext_pcts)

# Combine pix-to-pix comparison ----
sum_tbl <- read_csv(pix2pix_csv)
p2p_tbl <- right_join(sum_tbl, pcts_masked, by = c('map_code', 'mask', 'name'))

# Save
p2p_tbl %>% write_csv(compare_by_lc_csv)

# Look
# p2p_tbl %>% 
#   filter(map_code == 'internal') %>% 
#   select(mask, pct_of_lc) 

# join pix2pix with field data comparison (for unmasked)
ground_metrics <- read_csv(ground_compare_csv) %>% mutate(mask = 'Total')
agg_tbl <- p2p_tbl %>%
  filter(mask == 'Total') %>% 
  right_join(ground_metrics, by = 'map_code', suffix = c('', '_ground'))

agg_tbl_csv <- file.path(comparison_dir, str_c('07_aggregated_comparison_', agb_code, '.csv'))
agg_tbl %>% write_csv(agg_tbl_csv)

# ~ Difference maps (only for all LCs) -----------------------------------------
mskvals <- list(values = c(1,2,3,4,5,6), code = 'Total')
agb_fps[c(2,3,4,5)] %>% 
  purrr::walk(compare_by_given_lc, agb_fp, lc_fp,
              mskvals = mskvals, comparison_dir,
              return_obj = 'map', 
              boundary_fp = hti_poly_fp)
