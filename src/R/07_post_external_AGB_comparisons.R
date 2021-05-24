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

# Set variables ----------------------------------------------------------------
g0_variant <- 'med5'
year <- '2019'
input_level <- 'l2'
agb_code <- 'l2_WU'

# Input filepaths
agb_dir <- file.path('data', 'modeling', g0_variant)
(agb_fp <- list.files(agb_dir, str_c('agb_', input_level, '.*[^(sd)]\\.tif'), 
                      full.names = TRUE)[[2]])
hti_poly_fp <- "data/tidy/contextual_data/HTI_adm/HTI_adm0_fix.shp"

# External map filepaths
tidy_dir <- 'data/tidy'
glob_fp <- file.path(tidy_dir, 'biomass_maps', "GlobBiomass/N40W100_agb_crop_hti.tif")
esa_fp <- file.path(tidy_dir, 'biomass_maps', "ESA_CCI/ESA_agb17_crop_hti.tif")
avit_fp <- file.path(tidy_dir, 'biomass_maps', "Avitabile/Avitabile_AGB_crop_hti.tif")
bacc_fp <- file.path(tidy_dir, 'biomass_maps', "Baccini/20N_080W_t_aboveground_biomass_ha_2000_crop_hti.tif")

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

summarize_raster_differences <- function(int_r, ext_r){
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

scatter_AGB_differences_from_DF <- function(df, ext_name, prefix, filename = NA){
  
  # Sample
  if (length(df[[1]]) > 100000) {
    print("Sampling...")
    df <- df %>% sample_n(100000)
  }
  
  # Plot
  p <- ggplot(df, aes(x=external, y=internal)) + 
    geom_point(alpha=0.1, size=0.1, fill="royalblue", color="royalblue") +
    labs(y = expression(paste("Internal AGB estimate")), 
         x = expression(paste("External AGB estimate"))) +
    coord_fixed(ratio = 1, xlim=c(0, 320), ylim=c(0, 320)) +
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

perform_comparison <- function(agb_fp, ext_fp, ext_name, boundary_fp, input_level = 'l3'){
  # Compare our AGB to External AGB
  
  # Get filenames
  out_dir <- file.path(dirname(agb_fp), 'external_comparison', input_level)
  agb_res_fp <- file.path(out_dir, str_c(input_level, '_resamp_to', ext_name, '.tif'))
  diff_fp <- file.path(out_dir, str_c('diff_', ext_name, '_v_', input_level, '.tif'))
  dir.create(out_dir, recursive = TRUE) 
  
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
  
  df <- summarize_raster_differences(agb_res, agb_ext)

  # Histogram
  hist_diff_fp <- file.path(out_dir, str_c('hist_diff_', ext_name, '_v_', input_level, '.png'))
  p_hist <- plot_hist_density(df$diffs, -300, 300, bwidth=10, sample_size=500000,
                              ext_name, input_level, filename = hist_diff_fp)
  
  # Scatterplot
  scatter_diff_fp <- file.path(out_dir, str_c('scatter_diff_', ext_name, '_v_', input_level, '.png'))
  p_scatt <- scatter_AGB_differences_from_DF(df$diffs, ext_name, input_level, 
                                             filename = scatter_diff_fp)
  
  # Spatial map of differences
  N <- df$summary %>% filter(names=='N') %>% dplyr::select(value) %>% deframe
  if (N < 3000000) {
    
    map_diff_fp <- file.path(out_dir, str_c('map_diff_', ext_name, '_v_', input_level, '.png'))
    p_map <- plot_differences_map(diff_ras, ext_name, input_level, boundary_fp, filename = map_diff_fp)
    
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

get_summary_wrapper <- function(agb_fp, ext_fp, ext_name, input_level) {
  
  out_dir <- file.path(dirname(agb_fp), 'external_comparison', input_level)
  agb_res_fp <- file.path(out_dir, str_c(input_level, '_resamp_to', ext_name, '.tif'))  
  agb_res <- terra::rast(agb_res_fp)
  agb_ext <- crop_to_intersecting_extents(agb_res, terra::rast(ext_fp), return_r1=F)
  out <- summarize_raster_differences(agb_res, agb_ext)
  
}

# Run comparison and get graphics ----------------------------------------------
bacc_out <- perform_comparison(agb_fp, bacc_fp, 'Baccini', hti_poly_fp, agb_code)
glob_out <- perform_comparison(agb_fp, glob_fp, 'GlobB', hti_poly_fp, agb_code)
esa_out <- perform_comparison(agb_fp, esa_fp, 'ESA', hti_poly_fp, agb_code)
avit_out <- perform_comparison(agb_fp, avit_fp, 'Avitabile', hti_poly_fp, agb_code)

# Summary table ----
bacc_out <- get_summary_wrapper(agb_fp, bacc_fp, 'Baccini', agb_code)
glob_out <- get_summary_wrapper(agb_fp, glob_fp, 'GlobB', agb_code)
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
# agb_fp <- file.path(results_dir, 'tifs_by_R', str_c('agb18_v3_l3_hti_qLee_masked_filledLCpatches.tif'))
# esa_fp <- "data/ext_AGB_maps/ESA_CCI_Biomass/ESA_agb17_cropNA.tif"
agb.ras <- terra::rast(agb_fp)
glob <- terra::rast(glob_fp)
esa <- terra::rast(esa_fp)
avit <- terra::rast(avit_fp)
bacc <- terra::rast(bacc_fp)

# Add plot backscatter mean to polygons
field_agb_fp <- file.path(tidy_dir, 'survey_plots', 'plots_agb.rds')
plots_agb <- readRDS(field_agb_fp)

# Extract AGB means at survey sites
extract_stats_terra <- function(x, y, lyr_name, touches = FALSE, method='simple'){
  
  # Convert to SpatVector and extract AGB values
  x %>% 
    terra::extract(terra::vect(y), touches = touches, method = method) %>% 
    as_tibble() %>% 
    rename(value = 2) %>% 
    
    # Convert matrix to tibble and get mean and count for each polygon
    group_by(ID) %>% 
    summarise({{lyr_name}} := mean(value, na.rm=T) %>% 
                round(2)#,
              # ct = sum(!is.na(value)), 
              # sd = sd(value, na.rm = T)
              )
}

(agb_plots <- extract_stats_terra(agb.ras, plots_agb, lyr_name='agb', TRUE))
(glob_plots <- extract_stats_terra(glob, plots_agb, 'glob', TRUE))
(esa_plots <- extract_stats_terra(esa, plots_agb, 'esa', TRUE))
(avit_plots <- extract_stats_terra(avit, plots_agb, 'avit', TRUE))
(bacc_plots <- extract_stats_terra(bacc, plots_agb, 'bacc', TRUE))

# Join outputs and replace NAs with 0s
plots2 <- list(plots_agb, agb_plots, glob_plots, esa_plots, avit_plots, bacc_plots) %>% 
  purrr::reduce(left_join, by=c('plot_no' = 'ID')) %>% 
  mutate(across(everything(), ~ replace_na(.x, 0)))

# Extract just AGB and backscatter as dataframe
g0.agb <- plots2[c('AGB_ha', 'agb', 'glob', 'esa', 'avit', 'bacc')] %>% 
  st_drop_geometry %>%
  as.data.frame() %>% 
  mutate(AGB = as.vector(AGB_ha), 
         agb = na_if(agb, 0))

# Look at values
(backscatter <- c(
  all_mean=mean(g0.agb$agb, na.rm=T),
  all_sd=sd(g0.agb$agb, na.rm=T),
  bio_mean=mean(g0.agb[9:36, ]$agb, na.rm=T),
  bio_sd=sd(g0.agb[9:36, ]$agb, na.rm=T),
  min=min(g0.agb$agb, na.rm=T),
  max=max(g0.agb$agb, na.rm=T)))

# Get differences
df <- g0.agb %>% 
  transmute(
    agb = agb - AGB, 
    glob = glob - AGB, 
    esa = esa - AGB, 
    avit = avit - AGB,
    bacc = bacc - AGB
    )

df %>% 
  summarize(across(agb:bacc, ~ sum(!is.na(.x))))

rmse <- df %>% 
  summarize(across(agb:bacc, ~ mean(.x^2, na.rm=T))) %>% 
  pivot_longer(cols=everything(), values_to='mse') %>% 
  mutate(rmse = sqrt(mse))

bias <- df %>% 
  summarize(across(agb:bacc, ~ abs(sum(.x, na.rm=T) / sum(!is.na(.x)))))  %>% 
  pivot_longer(cols=everything(), values_to='bias')

me <- df %>% 
  summarize(across(agb:bacc, ~ mean(.x, na.rm=T)))  %>% 
  pivot_longer(cols=everything(), values_to='mean_diff')

diff_stats <- list(rmse, bias, me) %>% 
  purrr::reduce(left_join) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 2)))

diff_stats %>% write_csv(file.path(dirname(agb_fp), 'external_comparison', 
                            str_c('comparison_rmse_', input_level, '.csv')))

# Scatterplot - AGB against backscatter ----------------------------------------
plot_scatter_ext <- function(g0.agb, y_var, name){
  ggplot(g0.agb, aes(x=AGB, y=.data[[y_var]])) + geom_point() +
    labs(x = expression(paste("Observed AGB (Mg ha"^"-1", ")")), 
         y = expression(paste("Estimated AGB (Mg ha"^"-1", ")")), 
         subtitle = str_c(name, ", mean difference: ", 
                          round(mean(df[[y_var]], na.rm=T), 1), " +/- ", 
                          round(sd(df[[y_var]], na.rm=T), 1), " t/ha")) +
    coord_cartesian(xlim=c(0,160), ylim=c(0,160)) +
    geom_abline(intercept = 0, slope = 1, linetype='dashed', alpha=.5) +
    geom_text(label=rownames(g0.agb), hjust=0, vjust=0, size=3.5) + 
    theme_minimal()
}

p0 <- plot_scatter_ext(g0.agb, 'agb', name='Our AGB')
p1 <- plot_scatter_ext(g0.agb, 'glob', name='GlobBiomass')
p2 <- plot_scatter_ext(g0.agb, 'esa', name='ESA')
p3 <- plot_scatter_ext(g0.agb, 'avit', name='Avitabile')
p4 <- plot_scatter_ext(g0.agb, 'bacc', name='Baccini')

library('patchwork')
(p1 | p2 )/ (p3 | p4)
plot_fp <- file.path(dirname(agb_fp), 'external_comparison', str_c('comparison_scatter_', input_level, '.png'))
ggsave(plot_fp, width=20, height=20, units='cm')

plot_fp <- file.path(dirname(agb_fp), 'external_comparison', str_c('comparison_scatter_', input_level, '_ourAGB.png'))
ggsave(plot_fp, p0, width=10, height=10, units='cm')


