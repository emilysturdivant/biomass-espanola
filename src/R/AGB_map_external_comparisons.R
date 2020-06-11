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
library(gridExtra)

# FUNCTIONS -----------------------------------------------------------------------------------
crop_and_resample <- function(in_fp, out_fp, template, warp_method="near", na_value=NA, overwrite=F){
  # Resample raster to template raster with crop from GDAL for speed. 
  fun <- function(in_fp, out_fp, template, na_value=NA){
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
  if(!file.exists(esa_res_fp) | overwrite){
    r <- fun(in_fp, out_fp, template, na_value)
  }else{
    print("File already exists. Loading.")
    r <- read_stars(esa_res_fp)
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
plot_hist_density <- function(df, min=-200, max=200, bwidth=50, sample_size=10000){
  mu <- mean(df$value, na.rm=T)
  sd1 <- sd(df$value, na.rm=T)
  mn <- min(df$value, na.rm=T)
  mx <- max(df$value, na.rm=T)
  if (length(df[[1]]) > 1000000) {
    print("Sampling...")
    df <- df %>% sample_n(sample_size)
  }
  p <- ggplot(df, aes(x=value)) + 
    geom_histogram(aes(y=..density..), fill="#69b3a2", color="#e9ecef", alpha=0.7, binwidth = bwidth)+
    geom_density(alpha=.2) + 
    scale_x_continuous(name = expression(paste("Difference (t/ha): our AGB - Pan-tropical dataset")),
                       breaks = seq(min, max, 50)
                       # ,
                       # limits=c(min, max)
    ) +
    coord_cartesian(xlim=c(min, max)) +
    geom_vline(aes(xintercept=mu), color="black", linetype="solid", size=.5) +
    geom_vline(aes(xintercept=mu - sd1),
               color="black", linetype="dashed", size=.5)+
    geom_vline(aes(xintercept=mu + sd1),
               color="black", linetype="dashed", size=.5)+
    geom_vline(aes(xintercept=mn),
               color="black", linetype="solid", size=.5)+
    geom_vline(aes(xintercept=mx),
               color="black", linetype="solid", size=.5)+
    theme_minimal()
  return(p)
}

# Load AGB map -------------------------------------------------------------------------------
# agb.ras <- read_stars("results/tifs_by_R/agb18_v1_l0.tif")agb18_v1_l1_mask_Ap3WUw25.tif
# agb.ras <- read_stars('results/tifs_by_R/agb18_v1_l1_mask_Ap3WUw25_u20o310.tif')
agb.ras <- read_stars('results/tifs_by_R/agb18_v1_l2_mask_Ap3WUw25_u20_cap132.tif')
names(agb.ras) <- 'Our AGB 2018'

# Biomass loss from Baccini ------------------------------------------------------------------
in_fp <- "data/ext_AGB_maps/Baccini/Forest_Biomass_Loss/AGLB_Deforested_Tropical_America_2000.tif"
bmloss_fp <- "data/ext_AGB_maps/Baccini/Forest_Biomass_Loss/AGLB_Deforested_Tropical_America_2000_crop.tif"
r <- raster(in_fp)
gdalwarp(srcfile=in_fp, dstfile=bmloss_fp, te=st_bbox(agb.ras), tr=c(xres(r), yres(r)), tap=T, overwrite=T)


# GlobBiomass ----------------------------------------------------------------------------------------
# Crop GlobBiomass to our Hispaniola extent and convert 0s to NA
in_fp <- "data/ext_AGB_maps/GlobBiomass/N40W100_agb.tif"
crop_fp <- "data/ext_AGB_maps/GlobBiomass/N40W100_agb_crop.tif"
r <- raster(in_fp)
gdalwarp(srcfile=in_fp, dstfile=crop_fp, te=st_bbox(agb.ras), tr=c(xres(r), yres(r)), tap=T, overwrite=T)
r <- read_stars(crop_fp)
r[r == 0] <- NA
r %>% as("Raster") %>%
  writeRaster("data/ext_AGB_maps/GlobBiomass/N40W100_agb_cropNA.tif", 
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

# ESA ----------------------------------------------------------------------------------------
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
# Resample our AGB to ESA 
agb.ras <- read_stars('results/tifs_by_R/agb18_v1_l1_mask_Ap3WUw25_u20.tif')
agb_res <- agb.ras %>% st_warp(agb_esa, method="bilinear", use_gdal=TRUE)
agb_res %>% as("Raster") %>%
  writeRaster('results/tifs_by_R/resESAbl_agb18_v1_l1_mask_Ap3WUw25_u20.tif', 
              options=c("dstnodata=-99999"), overwrite=T)

# Compare
df_bl <- bind_cols(ESA=as.vector(agb_esa[[1]]), 
                AGB=as.vector(agb_res[[1]])) %>% 
  filter(!is.na(ESA), !is.na(AGB)) %>% 
  mutate(value=AGB-ESA)
s <- summary(df_bl$value) %>% as.vector()
cs <- c(df_bl$value %>% IQR(), df_bl$value %>% sd()) %>% as.vector()
pe <- cor.test(x=df_bl$ESA, y=df_bl$AGB, method = 'pearson') 
s_bl <- c(s, cs, dim(df_bl)[1], as.vector(pe$estimate)) %>% as_tibble()
hist(df_bl$value, xlim=c(-200, 200), breaks = 50)
# Plot histogram and density
(p <- plot_hist_density(df_bl, -300, 300, 10, sample_size=500000) + 
    scale_x_continuous(breaks = seq(-400, 400, 100)) + 
    xlim(c(-310, 310)) +
    xlab(expression(paste("Difference (t/ha): our AGB - ESA 2017"))))

# Resample our AGB to ESA 
agb.ras <- read_stars('results/tifs_by_R/agb18_v1_l1_mask_Ap3WUw25_u20.tif')
agb_resNN <- agb.ras %>% st_warp(agb_esa, method="near", use_gdal=TRUE)
agb_resNN %>% as("Raster") %>%
  writeRaster('results/tifs_by_R/resESAnn_agb18_v1_l1_mask_Ap3WUw25_u20.tif', 
              options=c("dstnodata=-99999"), overwrite=T)
# Compare
df_nn <- bind_cols(ESA=as.vector(agb_esa[[1]]), 
                AGB=as.vector(agb_resNN[[1]])) %>% 
  filter(!is.na(ESA), !is.na(AGB)) %>% 
  mutate(value=AGB-ESA)
s <- summary(df_nn$value) %>% as.vector()
cs <- c(df_nn$value %>% IQR(), df_nn$value %>% sd()) %>% as.vector()
pe <- cor.test(x=df_nn$ESA, y=df_nn$AGB, method = 'pearson') 
s_nn <- c(s, cs, dim(df_nn)[1], as.vector(pe$estimate)) %>% as_tibble()
hist(df_nn$value, xlim=c(-200, 200), breaks = 50)
# Plot histogram and density
(p <- plot_hist_density(df_nn, -300, 300, 10, sample_size=500000) + 
    scale_x_continuous(breaks = seq(-400, 400, 100)) + 
    xlim(c(-310, 310)) +
    xlab(expression(paste("Difference (t/ha): our AGB - ESA 2017"))))

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
# Plot scatterplot
df <- df_bl
if (length(df[[1]]) > 100000) {
  print("Sampling...")
  df <- df %>% sample_n(100000)
}
(p <- ggplot(df, aes(x=ESA, y=AGB)) + 
  geom_point(alpha=0.05, fill="royalblue", color="blue") +
  labs(y = expression(paste("Our AGB")), 
       x = expression(paste("ESA 2017 AGB"))) +
  theme_minimal() + 
  geom_smooth(method="lm", se=TRUE, fullrange=TRUE, level=0.95, col='black'))



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
in_fp <- "data/ext_AGB_maps/Avitabile_AGB_Map/Avitabile_AGB_Map.tif"
avit_fp <- "data/ext_AGB_maps/Avitabile_AGB_Map/Avitabile_AGB_crop.tif"
r <- raster(in_fp)
gdalwarp(srcfile=in_fp, dstfile=avit_fp, te=st_bbox(agb.ras), tr=c(xres(r), yres(r)), tap=T, overwrite=T)

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

