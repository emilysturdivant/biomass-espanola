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


# FUNCTIONS -----------------------------------------------------------------------------------
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
                       breaks = seq(min, max, 50),
                       limits=c(min, max)
    ) +
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
agb.ras <- read_stars("results/tifs_by_R/agb18_v1_l0.tif")agb18_v1_l1_mask_Ap3WUw25.tif
agb.ras <- read_stars('results/tifs_by_R/agb18_v1_l1_mask_Ap3WUw25_u20o310.tif')
agb.ras <- read_stars('results/tifs_by_R/agb18_v1_l1_mask_Ap3WUw25.tif')
names(agb.ras) <- 'Our AGB 2018'

# Filtering IN PROGRESS -----------------------------------------------------------------------------------
# Load EBImage up
require(EBImage)
# Read in your image
im = readImage('/path/to/your/image')
# Apply a median filter with window size of 7
im.med = medianFilter(im, size=7)
# Display image
display(im.med)

# ESA ----------------------------------------------------------------------------------------
fp_crop <- "data/ext_AGB_maps/ESA_CCI_Biomass/ESA_AGB_100m_2017_crop.tif"
bb <- st_bbox(agb.ras) #te=c(xmin,ymin,xmax,ymax)
gdalwarp(srcfile="data/ext_AGB_maps/ESA_CCI_Biomass/N40W100_ESACCI-BIOMASS-L4-AGB-MERGED-100m-2017-fv1.0.tif", 
         dstfile=fp_crop, 
         te=bb, overwrite=T)

# Load and resample ESA biomass map
agb_esa <- read_stars(fp_crop)
agb_esa[agb_esa == 0] <- NA
names(agb_esa) <- 'ESA 2017'
agb_esa <- agb_esa %>% st_warp(agb.ras, method="bilinear", use_gdal=TRUE)
agb_esa %>% saveRDS('results/R_out/ext_ESA_AGB17_resampledbl_stars.rds')
agb_esa %>% as("Raster") %>% 
  writeRaster('results/tifs_by_R/ext_ESAagb17_resampledBLna.tif', options=c("dstnodata=-99999"))

# Make difference map
diff_esa <- agb.ras - agb_esa
diff_esa %>% saveRDS('results/R_out/diffmap_ESA_stars.rds')
diff_esa %>% as("Raster") %>% 
  writeRaster('results/ext_comparisons/diff_ESAagb17_mask_Ap3WUw25.tif', overwrite=TRUE)
diff_esa <- readRDS('results/R_out/diffmap_ESA_stars.rds')

# Clip (mask and crop) to Haiti
# hti_poly <- st_read("data/contextual_data/HTI_adm/HTI_adm0_fix.shp")
# bb <- st_bbox(hti_poly) #te=c(xmin,ymin,xmax,ymax)
# gdalwarp(srcfile='results/ext_comparisons/diff_ESAagb17_mask_Ap3WUw25.tif', 
#          dstfile='results/ext_comparisons/diff_ESAagb17_hti_mask_Ap3WUw25.tif', 
#          cutline="data/contextual_data/HTI_adm/HTI_adm0_fix.shp",
#          te=bb, overwrite=T)
# diff_esa <- read_stars('results/ext_comparisons/diff_ESAagb17_hti_mask_Ap3WUw25.tif')

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

# Avitabile, resampled to ALOS ----------------------------------------------------------------------------------------
fp_crop <- "data/ext_AGB_maps/Avitabile_AGB_Map/Avitabile_AGB_crop.tif"
bb <- st_bbox(agb.ras) #te=c(xmin,ymin,xmax,ymax)
gdalwarp(srcfile="data/ext_AGB_maps/Avitabile_AGB_Map/Avitabile_AGB_Map.tif", 
         dstfile=fp_crop, 
         te=bb, overwrite=T)

# Load and resample biomass map
agb_avit <- read_stars(fp_crop)
names(agb_avit) <- 'Avitabile 2000–2013'
agb_avit <- agb_avit %>% st_warp(agb.ras, method="bilinear", use_gdal=TRUE)
agb_avit %>% saveRDS("results/R_out/ext_Avitabile_AGB_warpBL_stars.rds")
agb_avit %>% as("Raster") %>% 
  writeRaster("results/tifs_by_R/ext_Avitabile_AGB_warpBL.tif")
agb_avit <- read_stars("results/tifs_by_R/ext_Avitabile_AGB_warpBL.tif")

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
fp_crop <- "data/ext_AGB_maps/Baccini/20N_080W_t_agb_ha_2000_crop.tif"
bb <- st_bbox(agb.ras) #te=c(xmin,ymin,xmax,ymax)
gdalwarp(srcfile="data/ext_AGB_maps/20N_080W_t_aboveground_biomass_ha_2000.tif", 
         dstfile=fp_crop, 
         te=bb, overwrite=T)

# Load and resample biomass map
agb_bacc <- read_stars(fp_crop)
agb_bacc <- agb_bacc %>% st_warp(agb.ras, method="bilinear", use_gdal=TRUE)
names(agb_bacc) <- 'Baccini 2000'
agb_bacc %>% saveRDS("results/R_out/ext_Baccini2000_AGB_warpBL_stars.rds")
agb_bacc %>% as("Raster") %>% 
  writeRaster("results/tifs_by_R/ext_Baccini2000_AGB_warpBL.tif")

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
# Function to crop and mask raster to polygon
crop_and_mask_to_polygon <- function(in_fp, msk_poly, out_fp){
  tmp_fp <- tempfile(pattern = sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(in_fp)), 
                     tmpdir = tempdir(), 
                     fileext = "tif")
  # Crop
  gdalwarp(srcfile=in_fp, 
           dstfile=tmp_fp, 
           te=st_bbox(msk_poly), overwrite=T)
  # Mask
  msk_land <- raster(tmp_fp)
  msk_land <- msk_land %>% 
    mask(msk_poly, inverse=FALSE) 
  msk_land %>% writeRaster(out_fp, overwrite=T)
}

# AGB: Crop and mask to Haiti - 
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
           
           
           