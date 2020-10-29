# ---------------------------------------------------------------------------------------------
# Script to:
#     * Compare output AGB map to other biomass estimates
# Proceeds:
#     * postprocess_AGB_map.R - performs post-processing such as masking
#     * regression_AGB-g0.R - creates AGB map
#     * calculate_AGB.R - calculates AGB by plot from field data
#     * process_ALOS_tiles.R
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



plot_differences_map <- function(diff_ras, ext_name, prefix, save_fig=TRUE){
  
  # Convert raster to dataframe
  diff_df <- as.data.frame(diff_ras, xy=TRUE) %>% 
    drop_na() %>% 
    rename(diff = all_of(names(diff_ras)))
  
  # Load land polygon
  hti_poly <- st_read("data/contextual_data/HTI_adm/HTI_adm0_fix.shp")
  
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
  if(save_fig){
    
    map_diff_fp <- file.path('figures/ext_comparisons', str_c('map_diff_', ext_name, '_v_', prefix, '.png'))
    ggsave(map_diff_fp, plot=p, width=8, height=6)
    return(p)
    
  } else {
    
    return(p)
    
  }
}

crop_to_intersecting_extents <- function(r1, r2, return_r1=T, return_r2=T, return_bb=F) {
  # Get intersection of the bounding boxes of the two rasters
  bb <- st_intersection(terra::ext(r1) %>% as.vector %>% st_bbox %>% st_as_sfc, 
                        terra::ext(r2) %>% as.vector %>% st_bbox %>% st_as_sfc) %>% 
    st_bbox %>% as.vector 
  
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

perform_comparison <- function(agb_fp, ext_fp, ext_name){
  # Compare our AGB to External AGB
  
  # Get filenames
  prefix <- agb_fp %>% basename %>% str_split_fixed('_', 4) %>% .[,1:3] %>% str_c(collapse = '')
  
  agb_res_fp <- file.path('results/tifs_by_R', str_c(prefix, '_resamp_to', ext_name, '_bl_hti.tif'))
  diff_fp <- file.path('results/ext_comparisons', str_c('diff_', ext_name, '_v_', prefix, '.tif'))
  
  # Resample AGB to External
  if(!file.exists(agb_res_fp)){
    
    # Crop to intersection of the two extents
    out <- crop_to_intersecting_extents(terra::rast(agb_fp), terra::rast(ext_fp))
    agb.ras <- out[[1]]
    agb_glob <- out[[2]]
    
    agb_res <- agb.ras %>% terra::resample(agb_glob, method="bilinear")
    
    # Crop again
    agb_res <- crop_to_intersecting_extents(agb_res, terra::rast(ext_fp), return_r2=F)
    
    agb_res %>% terra::writeRaster(filename=agb_res_fp, 
                                   overwrite=TRUE, 
                                   wopt=list(gdal='COMPRESS=LZW'))
    
  } 
  
  # Difference map
  if(file.exists(diff_fp)){
    
    diff_ras <- terra::rast(diff_fp)
    
  } else {
    
    # Crop to intersection of the two extents
    agb_res <- terra::rast(agb_res_fp)
    agb_glob <- crop_to_intersecting_extents(agb_res, terra::rast(ext_fp), return_r1=F)
    
    diff_ras <- agb_res - agb_glob
    diff_ras %>% terra::writeRaster(filename=diff_fp, 
                                     overwrite=TRUE, 
                                     wopt=list(gdal='COMPRESS=LZW'))
  }
  
  # Get summary stats of differences and Pearson correlation
  # Crop to intersection of the two extents
  agb_res <- terra::rast(agb_res_fp)
  agb_glob <- crop_to_intersecting_extents(agb_res, terra::rast(ext_fp), return_r1=F)
  
  df <- summarize_raster_differences(agb_res, agb_glob)
  diffs <- df$diffs
  
  # Plots
  p_hist <- plot_hist_density(diffs, -300, 300, bwidth=10, sample_size=500000,
                              ext_name, prefix, save_fig=T)
  p_scatt <- scatter_AGB_differences_from_DF(diffs, ext_name, prefix, save_fig=T)
  
  N <- df$summary %>% filter(names=='N') %>% select(value) %>% deframe
  if (N < 3000000) {
    p_map <- plot_differences_map(diff_ras, ext_name, prefix, save_fig=T)
  } else {
      p_map <- NULL
    }
  
  # Return
  return(list(p_hist=p_hist, p_scatt=p_scatt, p_map=p_map, summary=df$summary))
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

################################################################################
# Run comparison and get graphics ------------------------------------------------------------------
# GlobBiomass
glob_out <- perform_comparison(agb_fp, glob_fp, 'GlobB')
glob_out$summary %>% View()

glob_out$summary %>% write_clip()

# ESA 
esa_out <- perform_comparison(agb_fp, esa_fp, 'ESA')
esa_out$summary %>% write_clip()

# Avitabile 
avit_out <- perform_comparison(agb_fp, avit_fp, 'Avitabile')
avit_out$summary %>% write_clip()

# Baccini 
bacc_out <- perform_comparison(agb_fp, bacc_fp, 'Baccini')
bacc_out$summary %>% write_clip()

################################################################################
# ESA --------------------------------------------------------------------------
# Difference map: L1 mask
agb_res <- read_stars('results/tifs_by_R/resamp_toESA_agb18v3l1a_blrasterresample.tif')
agb_esa <- read_stars(esa_fp)
diff_esa <- agb_res - agb_esa
diff_esa %>% as("Raster") %>% 
  writeRaster('results/ext_comparisons/diff_ESA_v_agb18v3l1a.tif', overwrite=TRUE)

# View
tmap_mode("view")
# bounding box for stars
test_bb <- st_bbox(c(xmin=-72.34, xmax=-72.17, ymin=18.3, ymax=18.38), crs=4326) %>% 
  st_as_sfc()
test_bb <- extent(c(xmin=-72.34, xmax=-72.3, ymin=18.4, ymax=18.44)[c(1,3,2,4)])
tm_shape(crop(agb_res, test_bb)) + 
  tm_raster(breaks=seq(0, 300, 10), palette=palette(hcl.colors(8, "viridis"))) +
  tm_shape(crop(agb_esa, test_bb)) + 
  tm_raster(breaks=seq(0, 300, 10), palette=palette(hcl.colors(8, "viridis"))) +
  tm_shape(crop(diff_esa, test_bb)) + 
  tm_raster(breaks=seq(-200, 200, 10), palette=palette(hcl.colors(8, "viridis")))
# View
tmap_mode("view")
# bounding box for stars
test_bb <- st_bbox(c(xmin=-72.34, xmax=-72.3, ymin=18.4, ymax=18.4), crs=4326) %>% 
  st_as_sfc()
# test_bb <- extent(c(xmin=-72.34, xmax=-72.3, ymin=18.4, ymax=18.44)[c(1,3,2,4)])
tm_shape(agb_res[test_bb]) + 
  tm_raster(breaks=seq(0, 300, 10), palette=palette(hcl.colors(8, "viridis")))

+
  tm_shape(agb_esa[test_bb]) + 
  tm_raster(breaks=seq(0, 300, 10), palette=palette(hcl.colors(8, "viridis"))) +
  tm_shape(diff_esa[test_bb]) + 
  tm_raster(breaks=seq(-200, 200, 10), palette=palette(hcl.colors(8, "viridis")))

# Cap values
df_bl <- df_bl %>% 
  mutate(AGBcap132=ifelse(AGB > 132, 132, AGB)) %>% 
  mutate(value=AGBcap132-ESA)
cs <- c(df_bl$value %>% IQR(), df_bl$value %>% sd()) %>% as.vector()
pe <- cor.test(x=df_bl$ESA, y=df_bl$AGBcap132, method = 'pearson') 
s <- summary(df_bl$value) %>% as.vector()
(s_bl <- c(s, cs, dim(df_bl)[1], as.vector(pe$estimate)) %>% as_tibble())
hist(df_bl$value, xlim=c(-200, 200), breaks = 50)
# Plot histogram and density
(p <- plot_hist_density(df_bl, -300, 300, 10, sample_size=500000) + 
    scale_x_continuous(breaks = seq(-400, 400, 100)) + 
    xlim(c(-310, 310)) +
    xlab(expression(paste("Difference (t/ha): our AGB - ESA 2017"))))


# Resample ESA
agb_esa <- crop_and_resample(in_fp="data/ext_AGB_maps/ESA_CCI_Biomass/N40W100_ESACCI-BIOMASS-L4-AGB-MERGED-100m-2017-fv1.0.tif", 
                             out_fp='results/tifs_by_R/ext_ESAagb17_resampledBLna.tif', 
                             template=agb.ras, na_value=0, overwrite=F)

# Difference map: L1 mask
agb.ras <- read_stars('results/tifs_by_R/agb18_v1_l1_mask_Ap3WUw25.tif')
diff_esa <- agb.ras - agb_esa
diff_esa %>% saveRDS('results/R_out/diffmap_ESA_stars.rds')
diff_esa %>% as("Raster") %>% 
  writeRaster('results/ext_comparisons/diff_ESAagb17_mask_Ap3WUw25.tif', overwrite=TRUE)
diff_esa <- readRDS('results/R_out/diffmap_ESA_stars.rds')

# Difference map: L2 capped values
agb.ras <- read_stars('results/tifs_by_R/agb18_v1_l2_mask_Ap3WUw25_u20_cap132.tif')
diff_esa <- agb.ras - agb_esa
diff_esa %>% saveRDS('results/R_out/diffmap_ESA_u20cap132_stars.rds')
diff_esa %>% as("Raster") %>% 
  writeRaster('results/ext_comparisons/diff_ESAagb17_mask-Ap3WUw25u20-cap132.tif', overwrite=TRUE)
diff_esa <- readRDS('results/R_out/diffmap_ESA_u20cap132_stars.rds')

# Histogram of differences
hist(diff_esa[[1]])
hist(diff_esa[[1]], xlim=c(-200, 200), breaks = 50)
diff_df <- diff_esa[[1]] %>% 
  as.vector() %>% 
  as_tibble() %>% 
  filter(!is.na(value)) 
summary(diff_df)
diff_df$value %>% IQR()
diff_df$value %>% sd()
length(diff_df$value)

# Plot histogram and density
(p <- plot_hist_density(diff_df, -300, 300, 10, sample_size=500000) + 
  scale_x_continuous(breaks = seq(-400, 400, 100)) + 
  xlim(c(-310, 310)) +
  xlab(expression(paste("Difference (t/ha): our AGB - ESA 2017"))))

ggsave('figures/ext_comparisons/hist_diff_ESA.png', width=150, height=100, units='mm')

# ESA resampled, AGB > 20----------------------------------------------------------------------------------------
diff_esa <- readRDS('results/R_out/diffmap_ESA_stars.rds')
msk_u20 <- readRDS("results/R_out/mask_AGB_under20_stars.rds")
diff_esa_masked <- diff_esa * msk_u20
# # Extract values as DF and summarize
# meow_msk_df_esam <- msk_u20[[1]] %>% 
#   as.vector() %>% 
#   as_tibble() %>% 
#   filter(!is.na(value)) 
# summary(meow_msk_df_esam)

# Histogram of differences
hist(diff_esa_masked[[1]], xlim=c(-200, 200), breaks = 50)
mu <- mean(diff_esa_masked[[1]], na.rm=T)
sd1 <- sd(diff_esa_masked[[1]], na.rm=T)
abline(v=c(mu-sd1, mu, mu+sd1), col='black')

# Extract values as DF and summarize
diff_df_esam <- diff_esa_masked[[1]] %>% 
  as.vector() %>% 
  as_tibble() %>% 
  filter(!is.na(value)) 

# Look at summary values
summary(diff_df_esam)
diff_df_esam$value %>% IQR()
diff_df_esam$value %>% sd()
length(diff_df_esam$value)

# Plot histogram and density
(p <- plot_hist_density(diff_df_esam, -200, 200, 10, sample_size=500000) + 
    scale_x_continuous(breaks = seq(-400, 400, 100)) + 
    xlim(c(-310, 310)) +
    xlab(expression(paste("Difference (t/ha): our AGB (resampled) - Avitabile"))))
ggsave('figures/ext_comparisons/hist_diff_Avitabile_AGB_mask.png', width=150, height=100, units='mm')

# Avitabile  -----------------------------------------------------------------------------------------------
# Resample our AGB to Avitabile 
agb_avit <- raster(avit_fp)
agb.ras <- raster(agb_fp)
agb_res <- agb.ras %>% raster::resample(agb_avit, method="bilinear")
agb_res %>% as("Raster") %>%
  writeRaster('results/tifs_by_R/resamp_toAvit_agb18v3l1a_bl.tif', 
              options=c("dstnodata=-99999"), overwrite=T)

# View
tmap_mode("view")
# bounding box for stars
test_bb <- st_bbox(c(xmin=-72.34, xmax=-72.17, ymin=18.3, ymax=18.38), crs=4326) %>% 
  st_as_sfc()
test_bb <- extent(c(xmin=-72.34, xmax=-72.27, ymin=18.3, ymax=18.38)[c(1,3,2,4)])
tm_shape(crop(agb_res, test_bb)) + 
  tm_raster(breaks=seq(0, 300, 10), palette=palette(hcl.colors(8, "viridis"))) +
  tm_shape(crop(agb_avit, test_bb)) + 
  tm_raster(breaks=seq(0, 300, 10), palette=palette(hcl.colors(8, "viridis")))

# Get summary stats of differences and Pearson correlation
dfs <- summarize_raster_differences(agb_res, agb_avit)

# Plot scatterplot
diffs <- dfs$diffs
(p <- scatter_AGB_differences_from_DF(diffs))

# Difference map: L1 mask
agb_res <- read_stars('results/tifs_by_R/resamp_toAvit_agb18v3l1a_bl.tif')
agb_avit <- read_stars(avit_fp)
diff_avit <- agb_res - agb_avit
diff_avit %>% as("Raster") %>% 
  writeRaster('results/ext_comparisons/diff_Avitabile_v_agb18v3l1a.tif', overwrite=TRUE)


# resampled to ALOS ----------------------------------------------------------------------------------------
agb_avit <- crop_and_resample(in_fp="data/ext_AGB_maps/Avitabile_AGB_Map/Avitabile_AGB_Map.tif", 
                              out_fp="results/tifs_by_R/ext_Avitabile_AGB_warpBL.tif", 
                              template=agb.ras, na_value=NA, overwrite=F)

# Make difference map
diff_avit <- agb.ras - agb_avit
diff_avit %>% saveRDS('results/R_out/diffmap_Avitabile_AGBv1_mask_Ap3WUw25.rds')
diff_avit %>% as("Raster") %>% 
  writeRaster('results/ext_comparisons/diffmap_Avitabile_AGBv1_mask_Ap3WUw25.tif')
diff_avit <- readRDS('results/R_out/diffmap_Avitabile_AGBv1_mask_Ap3WUw25.rds')

# Histogram of differences
hist(diff_avit[[1]])
hist(diff_avit[[1]], xlim=c(-200, 200), breaks = 50)

# Extract values to DF and summarize
diff_df_avit <- diff_avit[[1]] %>% 
  as.vector() %>% 
  as_tibble() %>% 
  filter(!is.na(value)) 
summary(diff_df_avit)
diff_df_avit$value %>% IQR()
diff_df_avit$value %>% sd()
length(diff_df_avit$value)

# Plot histogram and density
(p <- plot_hist_density(diff_df, -200, 200, 10, sample_size=500000) + 
    scale_x_continuous(breaks = seq(-400, 400, 100)) + 
    xlim(c(-310, 310)) +
    xlab(expression(paste("Difference (t/ha): our AGB - Avitabile, resampled"))))
ggsave('figures/ext_comparisons/hist_diff_Avitabile_AGB.png', width=150, height=100, units='mm')

# Avitabile resampled, AGB > 20----------------------------------------------------------------------------------------
diff_avit <- readRDS('results/R_out/diffmap_Avitabile_AGBv1_mask_Ap3WUw25.rds')
msk_u20 <- readRDS("results/R_out/mask_AGB_under20_stars.rds")
diff_avit_masked <- diff_avit * msk_u20

# Histogram of differences
hist(diff_avit_masked[[1]], xlim=c(-200, 200), breaks = 50)
mu <- mean(diff_avit_masked[[1]], na.rm=T)
sd1 <- sd(diff_avit_masked[[1]], na.rm=T)
abline(v=c(mu-sd1, mu, mu+sd1), col='black')
diff_df_avitm <- diff_avit_masked[[1]] %>% 
  as.vector() %>% 
  as_tibble() %>% 
  filter(!is.na(value)) 

summary(diff_df_avitm)
diff_df_avitm$value %>% IQR()
diff_df_avitm$value %>% sd()
length(diff_df_avitm$value)

# Plot histogram and density
(p <- plot_hist_density(diff_df_avitm, -200, 200, 10, sample_size=500000) + 
    scale_x_continuous(breaks = seq(-400, 400, 100)) + 
    xlim(c(-310, 310)) +
    xlab(expression(paste("Difference (t/ha): our AGB (resampled) - Avitabile"))))
ggsave('figures/ext_comparisons/hist_diff_Avitabile_AGB_mask.png', width=150, height=100, units='mm')

# Avitabile, AGB resampled ------------------------------------------------------------------------
# Load and resample biomass map
agb_avit <- read_stars(fp_crop)
names(agb_avit) <- 'Avitabile 2000–2013'
agb_rasA <- agb.ras %>% st_warp(agb_avit, method="bilinear", use_gdal=TRUE)
agb_rasA %>% saveRDS("results/R_out/ext_AGBresamp_AvitabileBL_stars.rds")
agb_rasA %>% as("Raster") %>% 
  writeRaster("results/tifs_by_R/ext_AGBresamp_AvitabileBL_stars.tif")
agb_rasA <- read_stars("results/tifs_by_R/ext_AGBresamp_AvitabileBL_stars.tif")

# Make difference map
diff_avit2 <- agb_rasA - agb_avit
diff_avit2 %>% saveRDS('results/R_out/diffmap_Avitabile_resAGBv1_mask_Ap3WUw25.rds')
diff_avit2 %>% as("Raster") %>% 
  writeRaster('results/ext_comparisons/diffmap_Avitabile_resAGBv1_mask_Ap3WUw25.tif')
diff_avit2 <- readRDS('results/R_out/diffmap_Avitabile_resAGBv1_mask_Ap3WUw25.rds')

# Histogram of differences
diff_df2 <- diff_avit2[[1]] %>% 
  as.vector() %>% 
  as_tibble() %>% 
  filter(!is.na(value)) 
summary(diff_df2)
diff_df2$value %>% IQR()
diff_df2$value %>% sd()
length(diff_df2$value)
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
getmode(round(diff_df2$value, 4))
getmode(diff_df2$value)

# Plot density
mu <- mean(diff_df2$value)
(p <- ggplot(diff_df2, aes(x=value)) + 
    geom_histogram(aes(y=..density..), fill="#69b3a2", color="#e9ecef", alpha=0.7, binwidth = 10)+
    geom_density(alpha=.2) + 
    scale_x_continuous(name = expression(paste("Difference: AGB resampled - Avitabile (t/ha)")),
                       breaks = seq(-200, 300, 50),
                       limits=c(-200, 300)
                       ) +
    geom_vline(aes(xintercept=mu), color="black", linetype="solid", size=.5)+
    geom_vline(aes(xintercept=mu-sd(value)),
               color="black", linetype="dashed", size=.5)+
    geom_vline(aes(xintercept=mu+sd(value)),
               color="black", linetype="dashed", size=.5)+
    theme_minimal())
ggsave('figures/ext_comparisons/hist_diff_Avitabile_resAGB.png', width=150, height=100, units='mm')

# Baccini ----------------------------------------------------------------------------------------
in_fp <- "data/ext_AGB_maps/Baccini/20N_080W_t_aboveground_biomass_ha_2000.tif"
bacc_fp <- "data/ext_AGB_maps/Baccini/20N_080W_t_aboveground_biomass_ha_2000_crop.tif"
r <- raster(in_fp)
gdalwarp(srcfile=in_fp, dstfile=bacc_fp, te=st_bbox(agb.ras), tr=c(xres(r), yres(r)), tap=T, overwrite=T)

agb_bacc <- crop_and_resample(in_fp="data/ext_AGB_maps/20N_080W_t_aboveground_biomass_ha_2000.tif", 
                              out_fp="results/tifs_by_R/ext_Baccini2000_AGB_warpBL.tif", 
                              template=agb.ras, na_value=NA, overwrite=F)

# Make difference map
diff_bacc <- agb.ras - agb_bacc
diff_bacc %>% saveRDS('results/R_out/diffmap_Baccini2000_AGBv1mask1.rds')
diff_bacc %>% as("Raster") %>% 
  writeRaster('results/ext_comparisons/diffmap_Baccini2000_AGBv1mask1.tif')

# Look at maps ---------------------------------------------------------------------------------
p_esa <- plot(agb_esa)
library(pryr)
object_size(p_esa)
p_esa

# Make crop bounding box/extent for testing
# bounding box for stars
test_bb <- st_bbox(c(xmin=-72.34, xmax=-72.17, ymin=18.3, ymax=18.38), crs=4326) %>% 
  st_as_sfc()
parks_hti <- # Parks for context
  st_read('data/contextual_data/OSM_free/hti_nature_reserves_osm.shp')
lyr_context <- tm_shape(parks_hti) + tm_borders()

# Set tmap to interactive viewing (vs. "plot")
tmap_mode("view")
# tmaptools::palette_explorer()

# View some masks
tm_shape(agb.ras[test_bb]) + 
  tm_raster(breaks=seq(0, 300, 10), palette=palette(hcl.colors(8, "viridis"))) +
tm_shape(agb_esa[test_bb]) + 
  tm_raster(breaks=seq(0, 300, 10), palette=palette(hcl.colors(8, "viridis"))) +
tm_shape(agb_avit[test_bb]) + 
  tm_raster(breaks=seq(0, 300, 10), palette=palette(hcl.colors(8, "viridis"))) +
tm_shape(agb_bacc[test_bb]) + 
  tm_raster(breaks=seq(0, 300, 10), palette=palette(hcl.colors(8, "viridis"))) +
tm_shape(diff_esa[test_bb]) + 
  tm_raster(breaks=seq(-200, 200, 10), palette="RdBu", n=9) +
tm_shape(diff_avit[test_bb]) + 
  tm_raster(breaks=seq(-200, 200, 10), palette="RdBu", n=9) +
tm_shape(diff_bacc[test_bb]) + 
  tm_raster(breaks=seq(-200, 200, 10), palette="RdBu", n=9) 

# View some masks
lyr_agb <- agb.ras[test_bb] %>% 
  tm_shape() + 
  tm_raster(breaks=seq(0, 300, 10), palette=palette(hcl.colors(8, "viridis")))
lyr_esa <- agb_esa[test_bb] %>% 
  tm_shape() + 
  tm_raster(breaks=seq(0, 300, 10), palette=palette(hcl.colors(8, "viridis")))
lyr_avit <- agb_avit[test_bb] %>% 
  tm_shape() + 
  tm_raster(breaks=seq(0, 300, 10), palette=palette(hcl.colors(8, "viridis")))
lyr_esa_diff <- diff_esa[test_bb] %>% 
  tm_shape() + 
  tm_raster(breaks=seq(-200, 200, 10), palette="RdBu", n=9)
lyr_avit_diff <- diff_avit[test_bb] %>% 
  tm_shape() + 
  tm_raster(breaks=seq(-200, 200, 10), palette="RdBu", n=9)

(lyr_agb + lyr_esa + lyr_avit + lyr_esa_diff + lyr_avit_diff + lyr_context)

# Get forest area ---------------------------------------------------------------------------------


# AGB: Crop and mask to Haiti - # THIS DIDN'T MASK ********
hti_poly <- st_read("data/contextual_data/HTI_adm/HTI_adm0_fix.shp")
crop_and_mask_to_polygon('results/tifs_by_R/agb18_v1_l1_mask_Ap3WUw25_u20.tif', 
                         hti_poly, 
                         'results/tifs_by_R/agb18_hti_v1_l1_mask_Ap3WUw25_u20.tif')

# Land mask: Crop and mask to Haiti
hti_poly <- st_read("data/contextual_data/HTI_adm/HTI_adm0_fix.shp")
crop_and_mask_to_polygon('results/masks/hisp18_maskLand.tif', 
                         hti_poly, 
                         'results/masks/hti18_maskLand_clip2border.tif')

# AGB: Crop and mask to DR
dr_poly <- st_read("data/contextual_data/DOM_adm/DOM_adm0.shp")
crop_and_mask_to_polygon('results/tifs_by_R/agb18_v1_l1_mask_Ap3WUw25_u20.tif', 
                         dr_poly, 
                         'results/tifs_by_R/agb18_dr_v1_l1_mask_Ap3WUw25_u20.tif')

# Land mask: Crop and mask to DR
dr_poly <- st_read("data/contextual_data/DOM_adm/DOM_adm0.shp")
crop_and_mask_to_polygon('results/masks/hisp18_maskLand.tif', 
                         dr_poly, 
                         'results/masks/dr18_maskLand_clip2border.tif')

# Get area and number of pixels in Haiti
land_area <- sum(st_area(hti_poly))
units(land_area) <- with(ud_units, ha)

msk_land <- # 1== ALOS 2018 land (i.e. Normal, Layover, and Shadowing)
  read_stars("results/masks/hti18_maskLand_clip2border.tif") 
ct_land_hti = sum(!is.na(msk_land[[1]]))

# Get number of 'forest' pixels and calculate percent of land and approximate area
agb_sat <- read_stars('results/tifs_by_R/agb18_hti_v1_l2_mask_Ap3WUw25_u20_cap132.tif')
ct_valid_agb = sum(!is.na(agb_sat[[1]]))
# Calculations
pct_of_land <- ct_valid_agb / ct_land_hti
est_area <- ct_valid_agb / ct_land_hti * land_area

mean(agb_sat[[1]], na.rm=T)


# Get area and number of pixels in DR
land_area_dr <- sum(st_area(dr_poly))
units(land_area_dr) <- with(ud_units, ha)

msk_land_dr <- # 1== ALOS 2018 land (i.e. Normal, Layover, and Shadowing)
  read_stars("results/masks/dr18_maskLand_clip2border.tif") 
ct_land_dr = sum(!is.na(msk_land_dr[[1]]))

# Get number of 'forest' pixels and calculate percent of land and approximate area
agb_sat <- read_stars('results/tifs_by_R/agb18_dr_v1_l1_mask_Ap3WUw25_u20.tif')
ct_valid_agb = sum(!is.na(agb_sat[[1]]))
# Calculations
land_area_dr
pct_of_land <- ct_valid_agb / ct_land_dr
est_area <- ct_valid_agb / ct_land_dr * land_area_dr



# Apply saturation point to AGB and save
agb_sat[agb_sat > saturation_pt] <- saturation_pt
agb_sat %>% saveRDS('results/R_out/agb18_hti_v1_l2_mask_Ap3WUw25_u20_cap132.rds')
agb_sat %>% as("Raster") %>% 
  writeRaster('results/tifs_by_R/agb18_hti_v1_l2_mask_Ap3WUw25_u20_cap132.tif')

agb_sat <- raster('results/tifs_by_R/agb18_hti_v1_l2_mask_Ap3WUw25_u20_cap132.tif')
ct_valid_agb = sum(!is.na(agb_sat1[[1]]))

# Get RMSE of external maps against our field data --------------------------------------------------
agb.ras <- raster('results/tifs_by_R/agb18_v1_l1_mask_Ap3WUw25_u20.tif')
esa <- raster("data/ext_AGB_maps/ESA_CCI_Biomass/ESA_agb17_cropNA.tif")
avit <- raster(avit_fp)
bacc <- raster(bacc_fp)

# Add plot backscatter mean to polygons
plots_agb <- readRDS('results/R_out/plots_agb.rds')
g0_AGB <- plots_agb %>% mutate(
  agb_est = raster_extract(agb.ras, plots_agb, fun = mean, na.rm = TRUE),
  esa_est = raster_extract(esa, plots_agb, fun = mean, na.rm = TRUE),
  avit_est = raster_extract(avit, plots_agb, fun = mean, na.rm = TRUE),
  bacc_est = raster_extract(bacc, plots_agb, fun = mean, na.rm = TRUE), 
  AGB = as.vector(AGB_ha))
  # ,
  # agb_dev = agb_est-AGB_ha,
  # esa_dev = esa_est-AGB_ha,
  # avit_dev = avit_est-AGB_ha,
  # bacc_dev = bacc_est-AGB_ha
  # )

# Extract just AGB and backscatter as dataframe
g0.agb <- g0_AGB[c('AGB', 'agb_est', 'esa_est', 'avit_est', 'bacc_est')] %>% 
  st_set_geometry(NULL) %>%
  as.data.frame()

# Look at values
(backscatter <- c(
  all_mean=mean(g0.agb$agb_est, na.rm=T),
  all_sd=sd(g0.agb$agb_est, na.rm=T),
  bio_mean=mean(g0.agb[9:36, ]$agb_est, na.rm=T),
  bio_sd=sd(g0.agb[9:36, ]$agb_est, na.rm=T),
  min=min(g0.agb$agb_est, na.rm=T),
  max=max(g0.agb$agb_est, na.rm=T)))

# Get differences
df <- g0.agb %>% transmute(
  agb_est = agb_est - AGB, 
  esa_est = esa_est - AGB, 
  avit_est = avit_est - AGB,
  bacc_est = bacc_est - AGB)
mean(df$agb_est, na.rm=T)
sd(df$agb_est, na.rm=T)

# Scatterplot - AGB against backscatter ---- ##############################################
g0.agb$AGB <- as.vector(g0.agb$AGB)
p1 <- ggplot(g0.agb, aes(x=AGB, y=agb_est)) + geom_point() +
  labs(x = expression(paste("Observed AGB (Mg ha"^"-1", ")")), 
       y = expression(paste("Estimated AGB (Mg ha"^"-1", ")"))) +
  coord_cartesian(xlim=c(0,160), ylim=c(0,160)) +
  geom_abline(intercept = 0, slope = 1, linetype='dashed', alpha=.5) +
  geom_text(label=rownames(g0.agb), hjust=0, vjust=0, size=3.5) + 
  ggtitle(str_c("Our AGB | MD: ", 
                round(mean(df$agb_est, na.rm=T), 1), " +/- ", 
                round(sd(df$agb_est, na.rm=T), 1), " t/ha"))
p2 <- ggplot(g0.agb, aes(x=AGB, y=esa_est)) + geom_point() +
  labs(x = expression(paste("Observed AGB (Mg ha"^"-1", ")")), 
       y = expression(paste("Estimated AGB (Mg ha"^"-1", ")"))) +
  coord_cartesian(xlim=c(0,160), ylim=c(0,160)) +
  geom_abline(intercept = 0, slope = 1, linetype='dashed', alpha=.5) +
  geom_text(label=rownames(g0.agb), hjust=0, vjust=0, size=3.5) + 
  ggtitle(str_c("ESA | MD: ", 
          round(mean(df$esa_est, na.rm=T), 1), " +/- ", 
          round(sd(df$esa_est, na.rm=T), 1), " t/ha"))
p3 <- ggplot(g0.agb, aes(x=AGB, y=avit_est)) + geom_point() +
  labs(x = expression(paste("Observed AGB (Mg ha"^"-1", ")")), 
       y = expression(paste("Estimated AGB (Mg ha"^"-1", ")"))) +
  coord_cartesian(xlim=c(0,160), ylim=c(0,160)) +
  geom_abline(intercept = 0, slope = 1, linetype='dashed', alpha=.5) +
  geom_text(label=rownames(g0.agb), hjust=0, vjust=0, size=3.5) + 
  ggtitle(str_c("Avitabile | MD: ", 
                round(mean(df$avit_est, na.rm=T), 1), " +/- ", 
                round(sd(df$avit_est, na.rm=T), 1), " t/ha"))
p4 <- ggplot(g0.agb, aes(x=AGB, y=bacc_est)) + geom_point() +
  labs(x = expression(paste("Observed AGB (Mg ha"^"-1", ")")), 
       y = expression(paste("Estimated AGB (Mg ha"^"-1", ")"))) +
  coord_cartesian(xlim=c(0,160), ylim=c(0,160)) +
  geom_abline(intercept = 0, slope = 1, linetype='dashed', alpha=.5) +
  geom_text(label=rownames(g0.agb), hjust=0, vjust=0, size=3.5) + 
  ggtitle(str_c("Baccini | MD: ", 
                round(mean(df$bacc_est, na.rm=T), 1), " +/- ", 
                round(sd(df$bacc_est, na.rm=T), 1), " t/ha"))
grid.arrange(p1, p2, p3, p4, nrow = 2, ncol=2)

