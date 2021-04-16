# ---------------------------------------------------------------------------------------------
# Script to:
#     * Standardize external (AGB) maps
# Proceeds:
#     * 4. postprocess_AGB_map.R - performs post-processing such as masking
#     * 3. regression_AGB-g0.R - creates AGB map
#     * 2. calculate_AGB.R - calculates AGB by plot from field data
#     * 1. process_ALOS_tiles.R
# Preceeds:
#     * post_external_AGB_comparisons.R -  Compare output AGB map to other biomass estimates
# Requires:
#     * AGB map (agbr)
#     * external AGB maps
# ---------------------------------------------------------------------------------------------

# Load libraries
library(raster)
library(tidyverse)
library(stars)
library(tmap)
library(gdalUtils)
library(units)
library(tools)
library(broom)
library(geobgu)
# library(gridExtra)
library(patchwork)
library(clipr)

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

crop_and_mask_to_polygon <- function(in_fp, msk_poly, out_fp){
  # Function to crop and mask raster to polygon
  tmp_fp <- tempfile(pattern = sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(in_fp)), 
                     tmpdir = tempdir(), 
                     fileext = "tif")
  # Crop
  r <- raster(in_fp)
  gdalwarp(srcfile=in_fp, dstfile=tmp_fp,
           te=st_bbox(msk_poly), tr=c(xres(r), yres(r)), tap=T, overwrite=T)
  
  # Mask
  msk_land <- raster(tmp_fp)
  msk_land <- msk_land %>% 
    mask(msk_poly, inverse=FALSE) 
  msk_land %>% writeRaster(out_fp, overwrite=T)
}

plot_hist_density <- function(df, min=-200, max=200, bwidth=50, sample_size=10000, 
                              ext_name, prefix, save_fig=TRUE){
  
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
    geom_histogram(aes(y=..density..), fill="#69b3a2", color="#e9ecef", alpha=0.7, binwidth = bwidth)+
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
               color="black", linetype="solid", size=.5,
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
  if(save_fig){
    
    hist_diff_fp <- file.path('figures/ext_comparisons', str_c('hist_diff_', ext_name, '_v_', prefix, '.png'))
    ggsave(hist_diff_fp, plot=p, width=8, height=4)
    return(p)
    
  } else {
    
    return(p)
    
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

scatter_AGB_differences_from_DF <- function(df, ext_name, prefix, save_fig=TRUE){
  
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
  if(save_fig){
    
    scatter_diff_fp <- file.path('figures/ext_comparisons', str_c('scatter_diff_', ext_name, '_v_', prefix, '.png'))
    ggsave(scatter_diff_fp, plot=p, width=6, height=6)
    # Return
    return(p)
    
  } else {
    
    return(p)
  }
  
}

# External map filepaths -------------------------------------------------------
glob_fp <- "data/ext_AGB_maps/GlobBiomass/N40W100_agb_cropNA.tif"
esa_fp <- "data/ext_AGB_maps/ESA_CCI_Biomass/ESA_agb17_cropNA.tif"
avit_fp <- "data/ext_AGB_maps/Avitabile_AGB_Map/Avitabile_AGB_crop.tif"
bacc_fp <- "data/ext_AGB_maps/Baccini/20N_080W_t_aboveground_biomass_ha_2000_crop.tif"
agb_fp <- 'results/tifs_by_R/agb18_v3_l1_mask_Ap3WUw25_u20_hti_qLee.tif'
agb_fp <- 'results/tifs_by_R/agb18_v3_l3_hti_qLee_masked_filledLCpatches.tif'

# Standardize external maps ----------------------------------------------------
# Load AGB map
agb.ras <- read_stars(agb_fp)

# Land cover from zipped file ------
# Crop GlobBiomass to our Hispaniola extent and convert 0s to NA
# **File too big to unzip and haven't been able to access without unzipping.**
in_fp <- "~/Downloads/dataset-satellite-land-cover-9bc57a16-b204-4a4a-a478-6c4a3d5510cc.zip/C3S-LC-L4-LCCS-Map-300m-P1Y-2018-v2.1.1.nc"
crop_fp <- "data/ext_AGB_maps/Landcover/C3S-LC-L4-LCCS-Map-300m-P1Y-2018-v2.1.1.tif"
r <- raster(in_fp)
gdalwarp(srcfile=in_fp, dstfile=crop_fp, te=st_bbox(agb.ras), tr=c(xres(r), yres(r)), tap=T, overwrite=T)
r <- read_stars(crop_fp)
r[r == 0] <- NA
r %>% as("Raster") %>%
  writeRaster("", 
              options=c("dstnodata=-99999"), overwrite=T)
lc_fp <- ""

# GlobBiomass ------
# Crop GlobBiomass to our Hispaniola extent and convert 0s to NA
in_fp <- "data/ext_AGB_maps/GlobBiomass/N40W100_agb.tif"
crop_fp <- "data/ext_AGB_maps/GlobBiomass/N40W100_agb_crop.tif"
r <- raster(in_fp)
gdalwarp(srcfile=in_fp, dstfile=crop_fp, te=st_bbox(agb.ras), tr=c(xres(r), yres(r)), tap=T, overwrite=T)
r <- read_stars(crop_fp)
r[r == 0] <- NA
r %>% as("Raster") %>%
  writeRaster(glob_fp, 
              options=c("dstnodata=-99999"), overwrite=T)

# Do the same for error (per-pixel uncertainty expressed as standard error in m3/ha)
in_fp <- "data/ext_AGB_maps/GlobBiomass/N40W100_agb_err.tif"
crop_fp <- "data/ext_AGB_maps/GlobBiomass/N40W100_agb_err_crop.tif"
r <- raster(in_fp)
gdalwarp(srcfile=in_fp, dstfile=crop_fp, te=st_bbox(agb.ras), tr=c(xres(r), yres(r)), tap=T, overwrite=T)
r <- read_stars(crop_fp)
r[r == 0] <- NA
r %>% as("Raster") %>%
  writeRaster("data/ext_AGB_maps/GlobBiomass/N40W100_agb_err_cropNA.tif", 
              options=c("dstnodata=-99999"), overwrite=T)

# ESA ------
# http://data.ceda.ac.uk/neodc/esacci/biomass/data/agb/maps/2017/v1.0/geotiff
# ftp://anon-ftp.ceda.ac.uk/neodc/esacci/biomass/data/agb/maps/2017/v1.0/geotiff/
# Crop ESA to our Hispaniola extent and convert 0s to NA
in_fp <- "data/ext_AGB_maps/ESA_CCI_Biomass/N40W100_ESACCI-BIOMASS-L4-AGB-MERGED-100m-2017-fv1.0.tif"
esa_fp <- "data/ext_AGB_maps/ESA_CCI_Biomass/ESA_agb17_crop.tif"
r <- raster(in_fp)
gdalwarp(srcfile=in_fp, dstfile=esa_fp, te=st_bbox(agb.ras), tr=c(xres(r), yres(r)), tap=T, overwrite=T)
# Convert 0s to NA and save in new file
agb_esa <- read_stars(esa_fp)
agb_esa[agb_esa == 0] <- NA
agb_esa %>% as("Raster") %>%
  writeRaster("data/ext_AGB_maps/ESA_CCI_Biomass/ESA_agb17_cropNA.tif", 
              options=c("dstnodata=-99999"), overwrite=T)
# Avitabile  ------
in_fp <- "data/ext_AGB_maps/Avitabile_AGB_Map/Avitabile_AGB_Map.tif"
avit_fp <- "data/ext_AGB_maps/Avitabile_AGB_Map/Avitabile_AGB_crop.tif"
r <- raster(in_fp)
gdalwarp(srcfile=in_fp, dstfile=avit_fp, te=st_bbox(agb.ras), tr=c(xres(r), yres(r)), tap=T, overwrite=T)

# Baccini -----
in_fp <- "data/ext_AGB_maps/Baccini/20N_080W_t_aboveground_biomass_ha_2000.tif"
bacc_fp <- "data/ext_AGB_maps/Baccini/20N_080W_t_aboveground_biomass_ha_2000_crop.tif"
r <- raster(in_fp)
gdalwarp(srcfile=in_fp, dstfile=bacc_fp, te=st_bbox(agb.ras), tr=c(xres(r), yres(r)), tap=T, overwrite=T)

# Biomass loss from Baccini ------
in_fp <- "data/ext_AGB_maps/Baccini/Forest_Biomass_Loss/AGLB_Deforested_Tropical_America_2000.tif"
bmloss_fp <- "data/ext_AGB_maps/Baccini/Forest_Biomass_Loss/AGLB_Deforested_Tropical_America_2000_crop.tif"
r <- raster(in_fp)
gdalwarp(srcfile=in_fp, dstfile=bmloss_fp, te=st_bbox(agb.ras), tr=c(xres(r), yres(r)), tap=T, overwrite=T)
