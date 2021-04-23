# ---------------------------------------------------------------------------------------------
# Script to:
#     * Correct AGB/g0 for soil moisture
# Proceeds:
#     * 
# Requires:
#     * Folder for SMAP products
#     * JAXA dates raster for given backscatter mosaic
#     * ESA soil moisture product from http://www.esa-soilmoisture-cci.org, downloaded via FTP
# ---------------------------------------------------------------------------------------------

# Load libraries
# library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(ggplot2) # package for plotting
library(lubridate)
library(geobgu)
library(units)
library(tidyverse)
library(tools)
library(stars)
# library(rhdf5)
library(reshape)
library(tmap)
tmap_mode('view')
library(rasterVis)
# library(smapr)
library(patchwork)

results_dir <- 'data/results'

# Process Surface/Rootzone Soil Moisture Analysis Update (SPL4SMAU) using smapr -----------------
# https://nsidc.org/data/SPL4SMAU/versions/4
# date accessed: 2020-06-19
# Functions ----
get_smap_rasters <- function(date, dir='data/raw/soil_moisture/SMAP', id='SPL4SMAU', 
                             version=4, name='/Analysis_Data/sm_surface_analysis', 
                             return_raster = FALSE){
  # Output filename
  d <- date %>% str_replace_all('-','')
  n <- name %>% basename()
  out_fp <- file.path(dir, id, str_glue('{n}_{d}.grd'))
  
  if(!file.exists(out_fp)){
    
    # Download files
    files <- smapr::find_smap(id, dates = date, version = version)
    local_files <- smapr::download_smap(files, overwrite = FALSE, verbose = T)
    sm_raster <- smapr::extract_smap(local_files, name)
    names(sm_raster) <- names(sm_raster) %>% str_extract("T[0-9]{4,}")
    
    # Crop to Hispaniola
    ext <- raster(file.path(results_dir, 'g0nu_HV/g0nu_2018_HV.tif')) %>% 
      projectExtent(crs(sm_raster[[1]])) %>% 
      extent()
    sm_raster <- crop(sm_raster, ext)
    
    # Save
    newr <- sm_raster %>% writeRaster(out_fp, format='raster', 
                                      overwrite=T)
    hdr(newr, format = "ENVI")
    
    # If return_raster == TRUE, return new raster
    if(return_raster) return(newr)
    
  } else {
    
    # If file already exists and return_raster == TRUE, load file and return
    if(return_raster) return(raster(out_fp))
    
  }
  
  # Return filename if return_raster == FALSE
  return(out_fp)
  
}

get_stack_means <- function(in_fp){
  # Load stack
  smr <- stack(in_fp)
  
  # Calculate mean
  mean_sm <- calc(smr, fun=mean)
  
  # Set layer name
  d <- in_fp %>% str_extract('[0-9]{8,}')
  names(mean_sm) <- str_glue('D{d}')
  
  # Return
  return(mean_sm)
}

download_crop_average_smap_sm <- function(id, name, dates, folder='data/raw/soil_moisture/SMAP'){
  # Download, crop, and save as multiband GRD rasters
  dates %>% lapply(get_smap_rasters, id=id, name=name)
  
  # Get daily means from multi-band rasters
  fps <- list.files(path=file.path(folder, id), pattern='grd$',
                    full.names=T, recursive=F, include.dirs=F)
  sm_means <- fps %>% lapply(get_stack_means)
  sm_means <- stack(sm_means)
  
  # Subset to only the requested dates
  sm_lyrs <- dates %>% map_chr(
    function(x) str_c('D', str_replace_all(x, '-', ''))
  )
  lyr_length <- length(sm_lyrs)
  sm_means <- sm_means[[sm_lyrs]]
  
  # Save raster stack - daily means from all multiband rasters in the folder
  out_fp <- file.path(folder, id, str_glue('{id}_daily_means_{lyr_length}.grd'))
  sm_means <- sm_means %>% writeRaster(out_fp, format='raster', 
                                       overwrite=T)
  hdr(sm_means, format = "ENVI")
  
  # Look
  # sm_means <- stack(out_fp)
  # levelplot(sm_means, main=id)
  return(sm_means)
}

download_crop_average_smap_7daymean <- function(id, name, dates, folder='data/raw/soil_moisture/SMAP'){
  # For each date, get the 7-day average
  # Get date +/- 3 days
  # for(day in dates)
  get_weekly_avg <- function(day, id, name, folder){
    # Get list of 7 days
    d1 <- seq(ymd(day) - days(3), by='days', length.out = 7)
    ras_names <- d1 %>% map_chr(
      function(x) str_c('D', str_replace_all(x, '-', ''))
      )
    
    # Get stack of daily means and subset to the 7-day period
    sm_means <- download_crop_average_smap_sm(id, name, d1, folder)
    week_sm <- sm_means[[ras_names]]
    
    # Get 7-day mean
    mean_sm <- calc(week_sm, fun=mean)
    
    # Set name
    d <- str_replace_all(day, '-', '')
    names(mean_sm) <- str_glue('D{d}')
    
    return(mean_sm)
  }
  # Download, crop, and save as multiband GRD rasters
  sm_means <- dates %>% lapply(get_weekly_avg, id=id, name=name, folder=folder)
  sm_means <- stack(sm_means)
  
  # Save raster stack - daily means from all multiband rasters in the folder
  out_fp <- file.path(folder, id, str_glue('{id}_weekly_means.grd'))
  sm_means <- sm_means %>% writeRaster(out_fp, format='raster', overwrite=T)
  hdr(sm_means, format = "ENVI")
  
  return(sm_means)
}

resamp_to_template <- function(src_fp, template_fp, dst_fp, align_coords=F, 
                               resamp_method='average', overwrite=T){
  src <- raster(src_fp)
  tmplt <- raster(template_fp)
  bb <- extent(tmplt)
  gdalUtils::gdalwarp(srcfile = src_fp, 
                      dstfile = dst_fp, 
                      s_srs = crs(src), 
                      # target spatial reference, extent, and resolution from SM
                      t_srs = crs(tmplt),
                      te = c(bb[1], bb[3], bb[2], bb[4]), 
                      tr = c(xres(tmplt), yres(tmplt)), 
                      tap = align_coords, # align the coords of the output extent to the target resolution
                      r = resamp_method, # average resampling, weighted avg of non-NA pixels
                      overwrite = overwrite)
}

mosaic_sm_to_match_dates <- function(sms, date_codes_df, date_mosaic_resamp, sm_mosaic_fp){
  # Initialize output brick
  out <- brick(nrow = nrow(sms), ncol = ncol(sms), 
               nl = nrow(date_codes_df), crs = crs(sms))
  values(out) <- NA
  names(out) <- names(sms)
  extent(out) <- extent(sms)
  
  # Iterate through dates
  for (i in seq(nrow(date_codes_df))) {
    
    # Get date codes
    sm_l <- date_codes_df %>% slice(i:i) %>% select(sm_lyr) %>% deframe
    g0_d <- date_codes_df %>% slice(i:i) %>% select(g0_date) %>% deframe
    
    # Mask SM to date zone
    smr <- sms[[sm_l]]
    smr[date_mosaic_resamp != as.numeric(g0_d)] <- NA
    
    # Add masked layer to output brick
    out[[sm_l]] = smr
  }
  
  # Merge into one SM layer
  sm <- merge(out)
  
  # Save
  sms <- sm %>% writeRaster(sm_mosaic_fp, overwrite=T)
}

scatter_AGB_vs_SM <- function(df){
  
  # Sample if there are more than 100,000 rows
  if (nrow(df) > 100000) {
    print("Sampling...")
    df <- df %>% sample_n(100000)
  }
  
  # Scatterplot
  p <- ggplot(df, aes(x=SM, y=AGB)) + 
    geom_point(alpha=0.1, fill="royalblue", color="blue") +
    labs(y = expression(paste("AGB estimate")), 
         x = expression(paste("Soil moisture"))) +
    theme_minimal() + 
    # Regression line with 95% CI band
    geom_smooth(method="lm", se=TRUE, fullrange=TRUE, level=0.95, col='black')
  
  # Return plot
  return(p)
  
}
# Filenames ----
agb_fp <- file.path(results_dir, 'tifs_by_R/agb18_v1_l1_mask_Ap3WUw25.tif')

date_fp <- file.path(results_dir, 'tifs_by_R/hisp18_date.tif')
date_resamp_fp <- file.path(results_dir, 'tifs_by_R/hisp18_date_SMAPres.tif')
# sm_fp <- file.path('data/SoilMoisture/SMAP', 'SPL4SMAU', 'SPL4SMAU_weekly_means.grd')
# sm_mosaic_fp <- file.path('data/SoilMoisture/SMAP', 'SPL4SMAU', 'SPL4SMAU_weekly_mosaic_jaxa18_dates.tif')
sm_fp <- file.path('data/raw/soil_moisture/SMAP', 'SPL4SMAU', 'SPL4SMAU_daily_means_12.grd')
sm_mosaic_fp <- file.path('data/raw/soil_moisture/SMAP', 'SPL4SMAU', 'SPL4SMAU_daily_mosaic_jaxa18_dates.tif')

# Resample dates mosaic to SMAP ----
# a_date <- raster(date_fp)
# sm <- raster(sm_fp, band = 1)
# bb <- extent(sm)
# gdalUtils::gdalwarp(srcfile = date_fp, 
#                     dstfile = date_resamp_fp, 
#                     s_srs = crs(a_date), 
#                     t_srs = crs(sm),
#                     te = c(bb[1], bb[3], bb[2], bb[4]), 
#                     tr = c(xres(sm), yres(sm)), 
#                     tap = T, 
#                     overwrite = T)
resamp_to_template(date_fp, sm, date_resamp_fp, align_coords = T, resamp_method = 'near')

# Download, crop, average soil moisture products for every date -----------------
# Load dates mosaic resampled to SMAP and SMAP stack
date_mosaic <- stack(date_resamp_fp)

# - Map date codes to date as SM layer names
date_codes_df <- date_mosaic %>% 
  # Unique date codes
  as.data.frame() %>% distinct %>% drop_na %>% 
  # Convert to SM layer names
  transmute(g0_date = layer, 
            date_str = g0_date %>% as_date('2014-05-24') %>% as.character,
            sm_lyr = date_str %>% str_replace_all('-', '') %>% str_c('D', .))
date_codes_df

# SPL4SMGP: Surface/Rootzone Soil Moisture Geophysical Data ----
id <- 'SPL4SMGP'
name <- '/Geophysical_Data/sm_surface'
sm_means <- download_crop_average_smap_sm(id, name, date_codes_df$date_str)
levelplot(sm_means, main=id)

smr <- stack(file.path(folder, id, str_glue('sm_surface_20160908.grd')))
levelplot(smr, main=id)

# SPL4SMAU: Surface/Rootzone Soil Moisture Analysis Update ----
id <- 'SPL4SMAU'
name <- '/Analysis_Data/sm_surface_analysis'
folder <- 'data/SoilMoisture/SMAP'
sm_means <- download_crop_average_smap_sm(id, name, date_codes_df$date_str)
levelplot(sm_means, main=id)

# Look at daily means for given dates
date_ct <- nrow(date_codes_df)
daily_fp <- file.path(folder, id, str_c(id, '_daily_means_', date_ct, '.grd'))
daily <- stack(daily_fp)
daily_means <- daily[[date_codes_df$sm_lyr]]
levelplot(daily_means, main=id)

# Look at image for one date - throws error
smr <- stack(file.path(folder, id, str_glue('sm_surface_20160908.grd')))
levelplot(smr, main=id)

# Weekly averages ----
sm_means <- download_crop_average_smap_7daymean(id, name, date_codes_df$date_str)
levelplot(sm_means, main=str_c(id, ' weekly'))

fps <- list.files(file.path(folder, id), pattern='sm_surface_analysis_2015.*\\.grd')

weekly_fp <- file.path(folder, id, str_c(id, '_weekly_means.grd'))
sms <- stack(weekly_fp)
levelplot(sms, main=id)

# Compare SMAP images to ALOS mosaic -------------------------------------------

# Create mosaic of SMAP by mapping to dates ------------------------------------

# Load dates mosaic resampled to SMAP
date_mosaic_resamp <- stack(date_resamp_fp)

# Load SMAP rasters
sms <- stack(sm_fp)

# Create mosaic matching dates
sm_mosaic <- mosaic_sm_to_match_dates(sms, date_codes_df, date_mosaic_resamp, sm_mosaic_fp)

# Load
fps <- list.files(file.path('data/SoilMoisture/SMAP', 'SPL4SMAU'), 
                  'SPL4SMAU_.*_mosaic_jaxa18_dates.tif', full.names = T)
sm_merges <- stack(fps)

# Look
levelplot(sm_merges)

# Regress SM vs. AGB -----------------------------------------------------------

# Resample AGB to SMAP ----
agb_resamp_fp <- file.path(results_dir, "tifs_by_R/agb18_v1_l1_m1_SMAPres_avg.tif")
if(!file.exists(agb_resamp_fp)) {
  resamp_to_template(agb_fp, sm_mosaic_fp, agb_resamp_fp)
}

# Load new raster
agb_res <- raster(agb_resamp_fp)
sm <- raster(sm_mosaic_fp)

# Convert rasters to data.frame
df_smres <- bind_cols(AGB=as.vector(agb_res), 
                SM=as.vector(sm)) %>% 
  filter(!is.na(AGB), !is.na(SM))

# Get Pearson's correlation coefficient
(pe <- cor.test(x = df_smres$AGB, y = df_smres$SM, method = 'pearson') )

# Look
tm_shape(agb_res) + tm_raster() +
  tm_shape(sm) + tm_raster()

# Plot scatterplot
cor_est <- pe$estimate %>% round(2) %>% as.numeric
(p_smres <- scatter_AGB_vs_SM(df_smres) + 
    ggtitle(str_c('SMAP resolution (AGB averaged): p = ', cor_est)))

# Resample AGB to SMAP (testing other summary stats) ----
agb_resamp_fp <- file.path(results_dir, "tifs_by_R/agb18_v1_l1_m1_SMAPres_bilinear.tif")
if(!file.exists(agb_resamp_fp)) {
  resamp_to_template(agb_fp, sm_mosaic_fp, agb_resamp_fp, resamp_method = 'bilinear')
}

# Load new raster
agb_res <- raster(agb_resamp_fp)
sm <- raster(sm_mosaic_fp)

# Convert rasters to data.frame
df_smres <- bind_cols(AGB=as.vector(agb_res), 
                      SM=as.vector(sm)) %>% 
  filter(!is.na(AGB), !is.na(SM))

# Get Pearson's correlation coefficient
(pe <- cor.test(x = df_smres$AGB, y = df_smres$SM, method = 'pearson') )
lm(AGB ~ SM, df_smres)
# Plot scatterplot
cor_est <- pe$estimate %>% round(2) %>% as.numeric
(p_smres <- scatter_AGB_vs_SM(df_smres) + 
    ggtitle(str_c('SMAP resolution (AGB bilinear resample): p = ', cor_est)))

# Filter to 20 t/ha bins
int <- 0
df_sub <- df_smres %>% filter(AGB > int & AGB <= int + 20)
(pe <- cor.test(x = df_sub$AGB, y = df_sub$SM, method = 'pearson') )
cor_est <- pe$estimate %>% round(2) %>% as.numeric
# Scatterplot
(p <- ggplot(df_sub, aes(x=AGB, y=SM)) + 
  geom_point(alpha=0.1, fill="royalblue", color="blue") +
  theme_minimal() + 
  # Regression line with 95% CI band
  geom_smooth(method="lm", se=TRUE, fullrange=TRUE, level=0.95, col='black'))

# Bin AGB to 20 t/ha intervals and get SM means
agb_breaks <- seq(0, 340, 20)
agb_labels <- head(agb_breaks, -1) + 10
interval <- 0.05
sm_breaks <- seq(0, 0.4, interval)
sm_labels <- head(sm_breaks, -1) + interval/2
df_out <- df_smres %>% 
  mutate(AGB_binned = cut(AGB, breaks=agb_breaks, labels = agb_labels),
         SM_binned = cut(SM, breaks=sm_breaks, labels = sm_labels))

(p_box_agb <- ggplot(df_out, aes(x=SM, y=AGB_binned)) +
    geom_boxplot(outlier.shape=1) + 
    # coord_flip() +
    labs(x ='Soil moisture',
         y = 'AGB +/- 10 t/ha') +
    theme_minimal() +
    theme(
      legend.position = 'bottom'
    ) +
    # Add mean point
    stat_summary(fun=mean, geom="point", aes(group=AGB_binned), 
                 position=position_dodge(.9), size=4, shape=20) +
    geom_point(aes(shape = "media"), alpha = 0, show.legend = F))

# Bin SM to intervals 
(p_box_sm <- ggplot(df_out, aes(x=AGB, y=SM_binned)) +
    geom_boxplot(outlier.shape=1) + 
    coord_flip() +
    labs(x ='AGB',
         y = str_glue('Soil moisture +/- {interval}')) +
    theme_minimal() +
    theme(legend.position = 'bottom' ) +
    # Add mean point
    stat_summary(fun=mean, geom="point", aes(group=SM_binned), 
                 position=position_dodge(.9), size=4, shape=20) +
    geom_point(aes(shape = "media"), alpha = 0, show.legend = F) +
    guides(shape=guide_legend(title=NULL, override.aes = list(alpha = 1), 
                              order=2)))

p_box_agb | p_box_sm

# Resample soil moisture to AGB ------------------------------------------------
sm_mosaic_resamp_fp <- file.path('data/SoilMoisture/SMAP', 'SPL4SMAU', 
                                 'SPL4SMAU_daily_mosaic_jaxa18_resampALOS.tif')
if(!file.exists(sm_mosaic_resamp_fp)) {
  resamp_to_template(sm_mosaic_fp, agb_fp, sm_mosaic_resamp_fp)
}

df_agbagg_res_fp <- file.path(results_dir, 'R_out/df_smap_v_agb_agbresx4.rds')
if(!file.exists(df_agbagg_res_fp)){
  # Load new raster
  agb <- raster(agb_fp)
  sm_res <- raster(sm_mosaic_resamp_fp)
  
  agb_agg <- aggregate(agb, 4)
  sm_agg <- aggregate(sm_res, 4)
  
  # Convert rasters to data.frame
  df_agbagg_res <- bind_cols(AGB=as.vector(agb_agg), 
                             SM=as.vector(sm_agg)) %>% 
    filter(!is.na(AGB), !is.na(SM))
  
  # Save
  df_agbagg_res %>% saveRDS(file.path(results_dir, 'R_out/df_smap_v_agb_agbresx4.rds'))
  
} else {
  df_agbagg_res <- readRDS(file.path(results_dir, 'R_out/df_smap_v_agb_agbresx4.rds'))
}

# Get Pearson's correlation coefficient
(pe <- cor.test(x = df_agbagg_res$AGB, y = df_agbagg_res$SM, method = 'pearson'))

cor_est <- round(pe$estimate, 2)
(p_agbx4res <- scatter_AGB_vs_SM(df_agbagg_res) +
    ggtitle(str_c('4x AGB resolution: p = ', cor_est)))

p_smres / p_agbx4res

# Filter to 20 t/ha bins
int <- 0
df_sub <- df_agbagg_res %>% filter(AGB > int & AGB <= int + 20)
(pe <- cor.test(x = df_sub$AGB, y = df_sub$SM, method = 'pearson') )
cor_est <- pe$estimate %>% round(2) %>% as.numeric
# Scatterplot
(p <- ggplot(df_sub, aes(x=AGB, y=SM)) + 
    geom_point(alpha=0.1, fill="royalblue", color="blue") +
    theme_minimal() + 
    # Regression line with 95% CI band
    geom_smooth(method="lm", se=TRUE, fullrange=TRUE, level=0.95, col='black'))


# Bin AGB to 20 t/ha intervals and get SM means
agb_breaks <- seq(0, 340, 20)
agb_labels <- head(agb_breaks, -1) + 10
interval <- 0.05
sm_breaks <- seq(0, 0.4, interval)
sm_labels <- head(sm_breaks, -1) + interval/2
df_out <- df_agbagg_res %>% 
  mutate(AGB_binned = cut(AGB, breaks=agb_breaks, labels = agb_labels),
         SM_binned = cut(SM, breaks=sm_breaks, labels = sm_labels))

(p_box_agb <- ggplot(df_out, aes(x=SM, y=AGB_binned)) +
    geom_boxplot(outlier.shape=1) + 
    coord_flip() +
    labs(x ='Soil moisture',
         y = 'AGB +/- 10 t/ha') +
    theme_minimal() +
    theme(
      legend.position = 'bottom'
    ) +
    # Add mean point
    stat_summary(fun=mean, geom="point", aes(group=AGB_binned), 
                 position=position_dodge(.9), size=4, shape=20) +
    geom_point(aes(shape = "media"), alpha = 0, show.legend = F))

# Bin SM to intervals 
(p_box_sm <- ggplot(df_out, aes(x=AGB, y=SM_binned)) +
    geom_boxplot(outlier.shape=1) + 
    coord_flip() +
    labs(x ='AGB',
         y = str_glue('Soil moisture +/- {interval}')) +
    theme_minimal() +
    theme(legend.position = 'bottom' ) +
    # Add mean point
    stat_summary(fun=mean, geom="point", aes(group=SM_binned), 
                 position=position_dodge(.9), size=4, shape=20) +
    geom_point(aes(shape = "media"), alpha = 0, show.legend = F) +
    guides(shape=guide_legend(title=NULL, override.aes = list(alpha = 1), 
                              order=2)))

p_box_agb | p_box_sm


# Compare to backscatter and AGB at field plots ================================
fn_suff <- '_qLee'
g0_fp <- file.path(results_dir, "g0nu_HV/g0nu_2018_HV_haitiR.tif"))
sm_mosaic_fp <- file.path('data/SoilMoisture/SMAP', 'SPL4SMAU', 'SPL4SMAU_daily_mosaic_jaxa18_dates.tif')

g0 <- read_stars(g0_fp)
sm <- read_stars(sm_mosaic_fp)

prefix <- 'plots_SMg0agb'

# Add plot backscatter mean to polygons
plots_agb <- readRDS(file.path(results_dir, 'R_out/plots_agb.rds'))
plots_agb %>% mutate(
  g0_mean = geobgu::raster_extract(g0, plots_agb, fun = mean, na.rm = TRUE),
  sm_mean = geobgu::raster_extract(sm, plots_agb, fun = mean, na.rm = TRUE)
) %>% 
  saveRDS(file.path(results_dir, str_c('R_out/', prefix, fn_suff,'.rds')))
g0_AGB <- readRDS(file.path(results_dir, str_c('R_out/', prefix, fn_suff,'.rds')))
g0_AGB %>% 
  st_write(file.path(results_dir, str_c("plots_values/", prefix, fn_suff, ".shp")), append=FALSE)

# Scatterplot - AGB against backscatter ----------------------------------------
g0_AGB <- g0_AGB %>% 
  mutate(area = set_units(area, NULL), 
         AGB_ha = set_units(AGB_ha, NULL))
(p <- ggplot(g0_AGB, aes(y=g0_mean, x=AGB_ha, color=sm_mean)) + 
    geom_point() +
    labs(x = expression(paste("Aboveground biomass (Mg ha"^"-1", ")")), 
         y = expression(paste("Radar backscatter, ",sigma['HV']^0," (m"^2, "/m"^2, ")"))) + 
    geom_smooth(method='lm', se=TRUE, fullrange=TRUE, level=0.95, col='black', size=0.2) + 
    geom_rug() + 
    theme_minimal())
(p <- ggplot(g0_AGB, aes(color=g0_mean, x=AGB_ha, y=sm_mean)) + 
    geom_point() +
    labs(x = expression(paste("Aboveground biomass (Mg ha"^"-1", ")")), 
         y = expression(paste("Soil moisture"))) + 
    geom_smooth(method='lm', se=TRUE, fullrange=TRUE, level=0.95, col='black', size=0.2) + 
    geom_rug() + 
    theme_minimal())
(p <- ggplot(g0_AGB, aes(y=g0_mean, color=AGB_ha, x=sm_mean)) + 
    geom_point() +
    labs(x = expression(paste("Soil moisture")), 
         y = expression(paste("Radar backscatter, ",sigma['HV']^0," (m"^2, "/m"^2, ")"))) + 
    geom_smooth(method='lm', se=TRUE, fullrange=TRUE, level=0.95, col='black', size=0.2) + 
    geom_rug() + 
    theme_minimal())
ggsave(str_c("figures/qc_plots/scatter_g0_sm_agb",fn_suff,".png"), width=15, height=13, units='cm')

# Multiple regression ----
ols <- lm(AGB ~ g0_mean + sm_mean, data=g0_AGB)
ols %>% saveRDS(file.path(results_dir, str_c("R_out/ols_AGB_g0_sm",fn_suff,".rds")))
summary(ols)
coefficients(ols)
confint(ols, level=0.95)
layout(matrix(c(1,2,3,4),2,2)) # optional 4 graphs/page 
plot(ols)
ols1 <- lm(AGB ~ g0_mean, data=g0_AGB)

anova(ols1, ols)








