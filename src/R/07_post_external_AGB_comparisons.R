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
      template <- stars::read_stars(template)
    }
    
    # Crop using GDAL
    print('Cropping...')
    fp_crop <- tempfile(pattern = str_glue(tools::file_path_sans_ext(basename(in_fp)),'_crop'), 
                        tmpdir = tempdir(), 
                        fileext = "tif")
    # fp_crop <- str_glue(file_path_sans_ext(basename(in_fp)),'_crop.',file_ext(in_fp))
    r <- terra::rast(in_fp)
    gdalUtils::gdalwarp(srcfile=in_fp, dstfile=fp_crop,
             te=sf::st_bbox(template), tr=c(xres(r), yres(r)), tap=T, overwrite=T)
    
    # Load and resample
    print('Resampling...')
    r <- stars::read_stars(fp_crop)
    
    if(!is.na(na_value)){
      r[r == na_value] <- NA
    }
    
    r <- r %>% stars::st_warp(template, method=warp_method, use_gdal=TRUE)
    r %>% as("Raster") %>%
      raster::writeRaster(out_fp, options=c("dstnodata=-99999"), overwrite=T)
    return(r) # output raster is stars format
  }
  
  # Only run if file doesn't already exist
  if(!file.exists(out_fp) | overwrite){
    
    r <- fun(in_fp, out_fp, template, na_value)
    
  } else {
    
    # If the file exists, load it
    print("File already exists. Loading.")
    r <- stars::read_stars(out_fp)
    
  }
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
  lc_names <- tibble(LC = 0:6, lc_name = lc_name)
  
  # Count values present in each LC
  cnts_vals <- as.data.frame(ras_lc) %>% 
    rename(LC = 1) %>% 
    count(LC) %>% 
    left_join(lc_names) %>% 
    filter(LC != 0) %>%
    mutate(pct_of_values = n / sum(n))
  
  # Count pixels in each LC
  cnts_lc <- as.data.frame(lc) %>% 
    rename(LC = 1) %>% 
    count(LC) %>% 
    left_join(lc_names)
  
  # Join to total LC counts and get percentages
  cnts_vals %>% 
    left_join(cnts_lc, by = c('LC', 'lc_name'), suffix = c('_vals', '_lc')) %>% 
    mutate(pct_of_lc = n_vals / n_lc,
           name = ext_name)
  
}

get_masked_agb_totals <- function(agb_x, lc_fp, out_dir) {
  # External names
  ext_name <- agb_x$name
  ext_fp <- agb_x$fp

  # Conditionally resample LC
  if (str_detect(ext_name, '^This')) {
    pix_area_ha <- (25^2)/(100^2)
    lc_res_fp <- 'data/tidy/landcover/Lemoiner/Haiti2017_agbres.tif'
    
    if(!file.exists(lc_res_fp)) resample_to_raster(lc_fp, ext_fp, lc_res_fp, method = 'near')
    
  } else if (ext_name == 'Avitabile') {    
    pix_area_ha <- (1000^2)/(100^2)
    lc_res_fp <- str_c(tools::file_path_sans_ext(lc_fp), '_res', ext_name, '.tif')
    
    # Resample LC and internal to external AGB res
    if(!file.exists(lc_res_fp)) resample_to_raster(lc_fp, ext_fp, lc_res_fp, method = 'near')

  } else {
    pix_area_ha <- (100^2)/(100^2)
    lc_res_fp <- lc_fp
  }
  
  out <- crop_to_intersecting_extents(terra::rast(lc_res_fp), 
                                      terra::rast(ext_fp))

  # Convert LC and AGB to dataframe
  rstack <- c(out$r1, out$r2)
  r_df <- as.data.frame(rstack) %>% rename(LC = 1, agb = 2)
  
  # Convert from density to volume
  sum_tbl <- r_df %>% 
    group_by(LC) %>% 
    summarize(
      mean_agb_ha = mean(agb, na.rm = TRUE),
      sum_agb = sum(agb, na.rm = TRUE) * pix_area_ha,
      max_agb_ha = max(agb, na.rm = TRUE),
      name = ext_name)
  
  # Save
  out_fp <- file.path(out_dir, str_glue('agb_sum_by_LC_{ext_name}.csv'))
  sum_tbl %>% write_csv(out_fp)
  
  # Return
  return(sum_tbl)
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

# ~ Sums and percents by LC ----
# Sums by LC ----
map_names <- agb_fps %>% map_chr('name') %>% as_tibble(rownames = 'map_code')
out_dir <- file.path(comparison_dir, 'by_LC', 'temp')
dir.create(out_dir, recursive = TRUE)
sum_tbl <- agb_fps %>%
  purrr::map_dfr(get_masked_agb_totals, lc_res_fp, out_dir) %>% 
  left_join(map_names, by = c(name = 'value'))

# sum_tbl <- list.files(out_dir, 'agb_sum', full.names = TRUE) %>% 
#   purrr::map_dfr(read_csv) %>% 
#   left_join(map_names, by = c(name = 'value'))

# Save
sum_tbl %>% write_csv(sums_csv)

# Look at values for Tree Cover class
# sum_tbl <- read_csv(sums_csv)
# sum_tbl %>% filter(LC == 4)

# Pct masked by LC ----
# Resample our AGB to ESA grid for comparison
agb_res_fp <- str_c(tools::file_path_sans_ext(agb_fp), '_resCCI.tif')
if(!file.exists(agb_res_fp)) {
  resample_to_raster(agb_fps$internal$fp, agb_fps$esa$fp, agb_res_fp)
}

res_fps <- agb_fps
res_fps$internal <- list(
  name = str_glue('Internal {agb_code}, 100m'),
  fp = agb_res_fp
)

map_names <- res_fps %>% map_chr('name') %>% as_tibble(rownames = 'map_code')
ext_pcts <- res_fps %>% 
  purrr::map_dfr(get_pcts_by_lc, lc_res_fp) %>% 
  left_join(map_names, by = c(name = 'value'))

# Save
ext_pcts %>% write_csv(ext_pcts_csv)

# Combine sums and percents
sum_tbl <- read_csv(sums_csv)
ext_pcts <- read_csv(ext_pcts_csv) %>% 
  select(-name) %>%
  rename(
    n_vals_100m = 'n_vals',
         n_lc_100m = 'n_lc', 
         pct_of_vals_100m = 'pct_of_values', 
         pct_of_lc_100m = 'pct_of_lc')
out_by_lc <- full_join(sum_tbl, ext_pcts, 
                       by = c('LC', 'map_code'))

smmry_csv <- file.path(comparison_dir, '07_overview_by_LC.csv')
out_by_lc %>% write_csv(smmry_csv)

# Aggregate to mask categories ----
# Sums and percents
out_by_lc <- read_csv(smmry_csv)

agg_fun <- function(mskvals, x) {
  x %>% 
    filter(LC %in% mskvals$values) %>% 
    group_by(name, map_code) %>% 
    summarize(sum_agb = sum(sum_agb, na.rm = TRUE),
              mean_agb_ha = mean(mean_agb_ha, na.rm = TRUE),
              max_agb_ha = max(max_agb_ha, na.rm = TRUE),
              n_vals_100m = sum(n_vals_100m, na.rm = TRUE),
              n_lc_100m = sum(n_lc_100m, na.rm = TRUE),
              pct_of_vals_100m = sum(pct_of_vals_100m, na.rm = TRUE)) %>% 
    mutate(pct_of_lc_100m =  n_vals_100m / n_lc_100m, 
           mask = mskvals$code)
}

out_by_mask <- mskvals_list %>% purrr::map_dfr(agg_fun, out_by_lc)

# Save
smmry_masks_csv <- file.path(comparison_dir, '07_overview_by_mask.csv')
out_by_mask %>% write_csv(smmry_masks_csv)

out_by_mask %>% filter(mask == 'TC')

# ~ Pix-to-pix comparison ----

# Comparison metrics & scatter density plots ----
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

# Combine pix-to-pix comparison ----
sum_tbl <- read_csv(pix2pix_csv)
p2p_tbl <- right_join(sum_tbl, out_by_mask, by = c('map_code', 'mask', 'name'))

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

agg_tbl %>% write_csv(agg_tbl_csv)

# ~ Difference maps (only for all LCs) -----------------------------------------
mskvals <- list(values = c(1,2,3,4,5,6), code = 'Total')
agb_fps[c(2,3,4,5)] %>% 
  purrr::walk(compare_by_given_lc, agb_fp, lc_fp,
              mskvals = mskvals, comparison_dir,
              return_obj = 'map', 
              boundary_fp = hti_poly_fp)

r <- terra::rast(agb_fp)
min(r)
