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

# crop_and_mask_to_polygon <- function(in_fp, msk_poly, out_fp){
#   # Function to crop and mask raster to polygon
#   tmp_fp <- tempfile(pattern = sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(in_fp)), 
#                      tmpdir = tempdir(), 
#                      fileext = "tif")
#   # Crop
#   r <- raster(in_fp)
#   gdalwarp(srcfile=in_fp, dstfile=tmp_fp,
#            te=st_bbox(msk_poly), tr=c(xres(r), yres(r)), tap=T, overwrite=T)
#   
#   # Mask
#   msk_land <- raster(tmp_fp)
#   msk_land <- msk_land %>% 
#     mask(msk_poly, inverse=FALSE) 
#   msk_land %>% writeRaster(out_fp, overwrite=T)
# }
# 
# crop_and_mask_to_poly_terra <- function(in_fp, msk_poly, out_fp, out_dtype='FLT4S'){
#   # Function to crop and mask raster to polygon
#   msk_poly <- if(is.character(msk_poly)) terra::vect(msk_poly)
#   
#   # Crop
#   terra::rast(in_fp) %>% 
#     terra::crop(msk_poly) %>% 
#     terra::mask(msk_poly, inverse=FALSE, 
#                 filename = out_fp, overwrite=T, 
#                 wopt = list(datatype=out_dtype, gdal='COMPRESS=LZW'))
#   
#   return(out_fp)
# }

crop_to_intersecting_extents <- function(r1, r2, return_r1=T, return_r2=T, return_bb=F) {
  # Get intersection of the bounding boxes of the two rasters
  bb <- sf::st_intersection(terra::ext(r1) %>% as.vector() %>% sf::st_bbox() %>% sf::st_as_sfc(), 
                            terra::ext(r2) %>% as.vector() %>% sf::st_bbox() %>% sf::st_as_sfc()) %>% 
    sf::st_bbox() %>% as.vector()
  
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

summarize_raster_differences <- function(int_r, ext_r, out_fp = NULL){
  # Make DF
  df <- bind_cols(external=as.vector(ext_r[[1]]), 
                  internal=as.vector(int_r[[1]])) %>% 
    filter(!is.na(external), !is.na(internal)) %>% 
    mutate(value=internal-external)
  
  # Get summary stats
  pe <- cor.test(x=df$external, y=df$internal, method = 'pearson') 
  s <- summary(df$value) %>% 
    tibble(names=names(.), value=.) %>% 
    mutate(value = as.numeric(value))
  cs <- tibble(
    names=c('IQR', 'SD', 'MAD', 'cor', 'N'), 
    value=c(IQR(df$value), 
            sd(df$value), 
            mean(abs(df$value)),
            pe$estimate, 
            dim(df)[1]))
  s_tbl <- bind_rows(s, cs)
  
  if(is.character(out_fp)){
    s_tbl %>% write_csv(out_fp)
  }
  
  # Return
  return(list(diffs=df, summary=s_tbl))
}

plot_hist_density <- function(df, min=-200, max=200, bwidth=50, sample_size=10000, 
                              ext_name, prefix, filename = NA){
  
  # Get statistics
  mu <- mean(df$value, na.rm=T)
  sd1 <- sd(df$value, na.rm=T)
  mn <- min(df$value, na.rm=T)
  mx <- max(df$value, na.rm=T)
  # cuts <- data.frame(ref = c('mean', 'SD', 'SD'), 
  #                    value = c(mu, mu-sd1, mu+sd1),
  #                    stringsAsFactors = F)
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

scatter_AGB_differences_from_DF <- function(df, ext_name, prefix, filename = NA,
                                            sample_size = 5000000){
  
  # Sample
  if (length(df[[1]]) > 500000) {
    print("Sampling...")
    df <- df %>% sample_n(sample_size)
  }
  
  # Plot
  p <- ggplot(df, aes(x=external, y=internal)) + 
    # geom_point(alpha=0.1, size=0.1, fill="royalblue", color="royalblue") +
    # stat_density_2d(aes(fill = ..level..), geom = "polygon") +
    geom_bin2d(bins = 75) +
    # geom_hex(bins = 150) +
    scale_fill_viridis_c(option = 'magma', direction = -1) +
    # scale_fill_viridis_c(option = 'magma', direction = -1, trans = 'log') +
    labs(y = expression(paste("Internal AGB estimate")), 
         x = expression(paste("External AGB estimate"))) +
    coord_fixed(ratio = 1, xlim=c(0, 300), ylim=c(0, 300)) +
    theme_minimal() + 
    geom_smooth(method="lm", #formula=y~0+x, 
                se=TRUE, fullrange=TRUE, level=0.95, col='black', size=.25)
  
  # Save/Return
  if(!is.na(filename)) {
    
    ggsave(filename, plot=p, width=6, height=6)
    # Return
    return(p)
    
  } 

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

perform_comparison <- function(agb_fp, ext_fp, ext_name, boundary_fp, agb_code = 'l3'){
  # Compare our AGB to External AGB
  
  # Get filenames
  out_dir <- file.path(dirname(agb_fp), 'external_comparison', agb_code)
  tif_dir <- file.path(out_dir, 'diff_tifs')
  agb_res_fp <- file.path(tif_dir, str_c(agb_code, '_resamp_to', ext_name, '.tif'))
  diff_fp <- file.path(tif_dir, str_c('diff_', ext_name, '_v_', agb_code, '.tif'))
  dir.create(tif_dir, recursive = TRUE) 
  
  # Input datatype
  in_dtype <- raster::dataType(raster::raster(agb_fp))
  
  # Resample AGB to External
  if(!file.exists(agb_res_fp)){
    
    # Crop to intersection of the two extents
    out <- crop_to_intersecting_extents(terra::rast(agb_fp), terra::rast(ext_fp))
    agb_ras <- out[[1]]
    agb_ext <- out[[2]]
    
    # Resample our AGB to external resolution
    agb_res <- agb_ras %>% terra::resample(agb_ext, method="bilinear")
    
    # Crop internal again
    agb_res <- crop_to_intersecting_extents(agb_res, terra::rast(ext_fp), return_r2=F)
    
    # Save resampled AGB
    agb_res %>% terra::writeRaster(filename = agb_res_fp, 
                                   overwrite = TRUE, 
                                   datatype = in_dtype)
  } 
  
  # Difference map
  if(file.exists(diff_fp)){
    
    diff_ras <- terra::rast(diff_fp)
    
  } else {
    
    # Crop external to intersection of the two extents
    agb_res <- terra::rast(agb_res_fp)
    agb_ext <- crop_to_intersecting_extents(agb_res, terra::rast(ext_fp), return_r1=F)
    
    # Subtract
    diff_ras <- agb_res - agb_ext
    diff_ras %>% terra::writeRaster(filename = diff_fp, 
                                    overwrite = TRUE, 
                                    datatype = in_dtype)
  }
  
  # Get summary stats of differences and Pearson correlation
  # Crop external to intersection of the two extents
  agb_res <- terra::rast(agb_res_fp)
  agb_ext <- crop_to_intersecting_extents(agb_res, terra::rast(ext_fp), return_r1=F)
  
  # Summary of differences
  smmry_diff_fp <- file.path(out_dir, str_c('summary_diff_', ext_name, '_v_', agb_code, '.csv'))
  df <- summarize_raster_differences(agb_res, agb_ext, smmry_diff_fp)

  # Histogram
  hist_diff_fp <- file.path(out_dir, str_c('hist_diff_', ext_name, '_v_', agb_code, '.png'))
  p_hist <- plot_hist_density(df$diffs, -300, 300, bwidth=10, sample_size=500000,
                              ext_name, agb_code, filename = hist_diff_fp)
  
  # Scatterplot
  scatter_diff_fp <- file.path(out_dir, str_c('scatter_diff_', ext_name, '_v_', agb_code, '.png'))
  p_scatt <- scatter_AGB_differences_from_DF(df$diffs, ext_name, agb_code, 
                                             filename = scatter_diff_fp)
  
  # Spatial map of differences
  N <- df$summary %>% filter(names=='N') %>% dplyr::select(value) %>% deframe
  if (N < 3000000) {
    
    map_diff_fp <- file.path(out_dir, str_c('map_diff_', ext_name, '_v_', agb_code, '.png'))
    p_map <- plot_differences_map(diff_ras, ext_name, agb_code, boundary_fp, filename = map_diff_fp)
    
  } else {
    p_map <- NULL
  }
  
  # Return
  out <- list(p_hist = p_hist, 
       p_scatt = p_scatt, 
       p_map = p_map, 
       summary = df$summary)
  return(out)
}

get_summary_wrapper <- function(agb_fp, ext_fp, ext_name, agb_code) {
  
  # Filepaths
  out_dir <- file.path(dirname(agb_fp), 'external_comparison', agb_code)
  agb_res_fp <- file.path(out_dir, str_c(agb_code, '_resamp_to', ext_name, '.tif'))  
  
  # Load and crop
  agb_res <- terra::rast(agb_res_fp)
  agb_ext <- crop_to_intersecting_extents(agb_res, terra::rast(ext_fp), return_r1=F)
  
  # Summarize differences
  out <- summarize_raster_differences(agb_res, agb_ext)
  
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

get_pcts <- function(x) {
  
  lyr_name <- x$name
  
  # Load and reclassify
  ras_clas <- terra::rast(x$fp) %>% 
    terra::classify(rbind(c(-Inf, Inf, 1),
                          c(NA, NA, 0)))
  
  # Crop
  msk_poly <- terra::vect(hti_poly_fp)
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
  
  df %>% 
    pivot_wider(names_from = group, values_from = count:pct) %>% 
    select(count_land, pct_value, pct_NA)
}

# Run comparison and get graphics ----------------------------------------------
glob_out <- perform_comparison(agb_fp, glob_fp, 'GlobB', hti_poly_fp, agb_code)
bacc_out <- perform_comparison(agb_fp, bacc_fp, 'Baccini', hti_poly_fp, agb_code)
rm(bacc_out)
esa_out <- perform_comparison(agb_fp, esa_fp, 'ESA', hti_poly_fp, agb_code)
rm(esa_out)
avit_out <- perform_comparison(agb_fp, avit_fp, 'Avitabile', hti_poly_fp, agb_code)

# Summary table ----
glob_out <- get_summary_wrapper(agb_fp, glob_fp, 'GlobB', agb_code)
bacc_out <- get_summary_wrapper(agb_fp, bacc_fp, 'Baccini', agb_code)
esa_out <- get_summary_wrapper(agb_fp, esa_fp, 'ESA', agb_code)
avit_out <- get_summary_wrapper(agb_fp, avit_fp, 'Avitabile', agb_code)

# Histogram
avit_out <- get_summary_wrapper(agb_fp, avit_fp, 'Avitabile', agb_code)

ext_name <- 'Avitabile'
out_dir <- file.path(dirname(agb_fp), 'external_comparison', agb_code)
hist_diff_fp <- file.path(out_dir, str_c('hist_diff_', ext_name, '_v_', agb_code, '_2.png'))
p_hist <- plot_hist_density(avit_out$diffs, -300, 300, bwidth=5, sample_size=500000,
                            ext_name, agb_code, filename = hist_diff_fp)

# Join
summary_df <- plyr::join_all(
  list(setNames(glob_out$summary, c('name', 'GlobBiomass')),
       setNames(esa_out$summary, c('name', 'ESA')),
       setNames(avit_out$summary, c('name', 'Avitabile')),
       setNames(bacc_out$summary, c('name', 'Baccini'))),
  by = 'name', type = 'left')

# Save summary table
summary_csv <- file.path(dirname(agb_fp), 'external_comparison', 
                         str_c('comparison_summaries_', agb_code, '.csv'))
summary_df %>% mutate(across(where(is.numeric), ~ round(.x, 2))) %>% 
                        write_csv(summary_csv)

rm(bacc_out, esa_out, avit_out, glob_out)

# Get RMSE of external maps against our field data -----------------------------
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

df <- plot_means %>% 
  pivot_longer(cols = internal:bacc, values_to = 'est') %>% 
  mutate(diff = est - AGB_ha) 

diff_stats <- df %>% 
  group_by(name) %>% 
  summarize(
    r.squared = summary(lm(est ~ AGB_ha))$r.squared,
    adj.r.squared = summary(lm(est ~ AGB_ha))$adj.r.squared,
    N = sum(!is.na(diff)), 
    RMSD = sqrt(mean(diff^2, na.rm=T)),
    bias = abs(sum(diff, na.rm=T) / sum(!is.na(diff))),
    MBD = mean(diff)
    ) %>% 
  mutate(across(where(is.numeric), ~ signif(.x, 4)))

# Save
tbl_csv <- file.path(dirname(agb_fp), 'external_comparison', agb_code,
                     str_c('07_comparison_rmse_', agb_code, '.csv'))
diff_stats %>% write_csv(tbl_csv)

# Percent of land area included in map ----
ext_pcts <- agb_fps %>% purrr::map_dfr(get_pcts) %>% 
  mutate(map = names(agb_fps),
         map_name = map_chr(agb_fps, "name"))

# Save as CSV
ext_pcts_csv <- file.path('data/reports', str_glue('07_external_pcts_masked.csv'))
ext_pcts %>% write_csv(ext_pcts_csv)

# Combine ----
ext_metrics <- diff_stats %>% 
  left_join(ext_pcts, by = c(name = 'map')) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 3)))

# Save as CSV
ext_report_csv <- file.path('data/reports', str_glue('07_ext_comparison_metrics.csv'))
ext_metrics %>% write_csv(ext_report_csv)

# Pie chart ----
ext_pcts <- read_csv(ext_pcts_csv)
df <- ext_pcts %>% 
  pivot_longer(cols = starts_with('pct'), 
               names_prefix = 'pct_', 
               names_to = 'group', 
               values_to = 'pct') %>% 
  filter(map != 'internal')

(bp <- ggplot(df, 
              aes(x=pct, y=map_name, fill=group)) +
    geom_bar(stat = "identity", position = 'fill') +
    theme_minimal() + 
    labs(fill="",x=NULL,y=NULL,title="",caption=""))

# Table
(d2 <- df %>% 
    mutate(pct_coverage = str_c(round(pct*100, digits = 1), '%')) %>% 
    arrange(desc(pct)) %>% 
    filter(!group %in% c('NA')) %>% 
    select(map_name, pct_coverage))

tab <- gridExtra::tableGrob(d2, rows = NULL, cols = c('Map', 'Percent of land'))

# Display side-by-side
wrap_plots(tab, bp)

# Scatterplot - AGB against backscatter ----------------------------------------
plot_scatter_ext <- function(plot_means, y_var, name){
  
  ggplot(plot_means, aes(x=AGB_ha, y=.data[[y_var]])) + geom_point() +
    labs(x = expression(paste("Observed AGB (Mg ha"^"-1", ")")), 
         y = expression(paste("Estimated AGB (Mg ha"^"-1", ")")), 
         subtitle = str_c(name, ", mean difference: ", 
                          round(mean(df[[y_var]], na.rm=T), 1), " +/- ", 
                          round(sd(df[[y_var]], na.rm=T), 1), " t/ha")) +
    coord_cartesian(xlim=c(0,160), ylim=c(0,160)) +
    geom_abline(intercept = 0, slope = 1, linetype='dashed', alpha=.5) +
    geom_text(label=rownames(plot_means), hjust=0, vjust=0, size=3.5) + 
    theme_minimal()
  
}

p0 <- plot_scatter_ext(plot_means, 'internal', name='Our AGB')
p1 <- plot_scatter_ext(plot_means, 'GlobB', name='GlobBiomass')
p2 <- plot_scatter_ext(plot_means, 'ESA', name='ESA')
p3 <- plot_scatter_ext(plot_means, 'Avitabile', name='Avitabile')
p4 <- plot_scatter_ext(plot_means, 'Baccini', name='Baccini')

library('patchwork')
(p1 | p2 )/ (p3 | p4)
plot_fp <- file.path(dirname(agb_fp), 'external_comparison', agb_code,
                     str_c('comparison_scatter_', agb_code, '.png'))
ggsave(plot_fp, width=20, height=20, units='cm')

plot_fp <- file.path(dirname(agb_fp), 'external_comparison', agb_code,
                     str_c('comparison_scatter_', agb_code, '_ourAGB.png'))
ggsave(plot_fp, p0, width=10, height=10, units='cm')

# Histogram of external map ----
plot_agb_hist_density <- function(x, agb_pal, bwidth=5, sample_size=1000000, 
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

bacc_hist <- plot_agb_hist_density(agb_fps$bacc, agb_pal)
hist_fp <- file.path('data/reports/external_agb', str_glue('agb_hist_Baccini.png'))
dir.create(dirname(hist_fp))
ggsave(hist_fp, width = 5, height = 3, dpi = 120)

(esa_hist <- plot_agb_hist_density(agb_fps$esa, agb_pal))
hist_fp <- file.path('data/reports/external_agb', str_glue('agb_hist_ESA.png'))
ggsave(hist_fp, width = 5, height = 3, dpi = 120)

(glob_hist <- plot_agb_hist_density(agb_fps$glob, agb_pal, xlim=NULL))
hist_fp <- file.path('data/reports/external_agb', str_glue('agb_hist_GlobB.png'))
ggsave(hist_fp, width = 5, height = 3, dpi = 120)

(avit_hist <- plot_agb_hist_density(agb_fps$avit, agb_pal, xlim=NULL))
hist_fp <- file.path('data/reports/external_agb', str_glue('agb_hist_Avit.png'))
ggsave(hist_fp, width = 5, height = 3, dpi = 120)

