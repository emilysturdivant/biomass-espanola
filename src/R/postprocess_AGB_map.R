# *************************************************************************************************
# Script to:
#     * Perform post-process on AGB map 
# Proceeds:
#     * regression_AGB-g0.R - creates AGB map
#     * calculate_AGB.R - calculates AGB by plot from field data
#     * process_ALOS_tiles.R
# Requires:
#     * AGB map (agbr)
#     * pre-created masks
#
# *************************************************************************************************

# Load libraries
library(raster)
library(tidyverse)
library(stars)
library(geobgu)
library(broom)
library(gdalUtils)
library(tmap)
require(graphics)

# Functions ------------------------------------------------------------------------------------------
get_counts <- function(valuesraster, msk_zone, msk_0, threshold=132){
  # For given area of interest (msk_zone!=NA), 
  # count non-masked values (msk_0!=NA) and values>threshold after masking 
  #  - msk_zone - defines zone of interest ("total")
  #  - msk_0 - which pixels to count inside msk_zone ("masked")
  #  - threshold - count values > threshold after masking
  # 1. Convert rasters to DF
  df <- 
    tibble(zone=as.vector(msk_zone[[1]]),
           agb=as.vector(valuesraster[[1]]) %>% round(4), 
           mask=as.vector(msk_0[[1]])) %>% 
    filter(!is.na(zone)) 
  total = length(df[[1]])
  # Filter by mask
  df1 <- df %>% filter(!is.na(mask))
  total1 = length(df1[[1]])
  # Filter by (upper) threshold value
  df2 <- df1 %>% filter(!agb > threshold)
  total2 = length(df2[[1]])
  # Tally it up
  cts_df <- tibble(
    masked = total-total1,
    str_c('over_',threshold) = total1-total2,
    valid=total2
  ) %>% 
    pivot_longer(everything()) %>% 
    mutate(pct = value/sum(value))
}

# Work with different versions of AGB  --------------------------------------------------------------------
# Load AGB map
fn_suff <- '_qLee' # filename suffix (from regression_AGB-g0.R)
agb.ras <- read_stars(str_c("results/tifs_by_R/agb18_v3_l0",fn_suff,".tif"))

# Load default mask (ALOS, water, and urban)
msk_AWUwb <- # NA==ALOS mask, LC17 Water and Urban, and OSM water with 25 m buffer
  readRDS(str_c("results/R_out/mask_AGB_ALOS_WaterUrban_water25_stars",fn_suff,".rds"))
msk_land <- # 1== ALOS 2018 land (i.e. Normal, Layover, and Shadowing)
  read_stars("results/masks/hti18_maskLand_clip2border.tif") 

# set extent and crop
bb <- st_intersection(st_bbox(agb.ras) %>% st_as_sfc(), 
                      st_bbox(msk_land) %>% st_as_sfc()) %>% 
  st_bbox() %>% as.vector()
bbex <- extent(bb[c(1,3,2,4)])
agb.ras <- agb.ras %>% as("Raster") %>% crop(bbex)
agb.ras %>% writeRaster(str_c("results/tifs_by_R/agb18_v3_l0_hti",fn_suff,".tif"))
msk_land <- msk_land %>% as("Raster") %>% crop(bbex)
msk_AWUwb <- msk_AWUwb %>% as("Raster") %>% crop(bbex)

agb.ras <- agb.ras %>% st_as_stars()
msk_land <- msk_land %>% st_as_stars()
msk_AWUwb <- msk_AWUwb %>% st_as_stars()

# AGB <20 mask
msk_u20 <- agb.ras
msk_u20[msk_u20 < 20] <- NA
msk_u20[msk_u20 >= 20] <- 1
msk_u20 %>% saveRDS(str_c("results/R_out/mask_AGB_under20_stars",fn_suff,"_hti.rds"))

# Count NAs
df <- get_counts(agb.ras, msk_land, msk_AWUwb, 310)
msk_AWUwb_u20 <- msk_AWUwb * msk_u20
df <- get_counts(agb.ras, msk_land, msk_AWUwb_u20, 300)
df <- get_counts(agb.ras, msk_land, msk_u20, 300)
df <- get_counts(agb.ras, msk_land, msk_AWUwb_u20, 132)

# get Tree Cover
df <- get_counts(agb.ras, msk_land, msk_u20, 300)



# Original  ---------------------------------------------------------------------------------------

# Load AGB map
agb.ras <- read_stars(str_c("results/tifs_by_R/agb18_v3_l0",fn_suff,".tif"))

# Make masks based on AGB --------------------------------------------------------------------------
# AGB <20 mask
msk_u20 <- agb.ras
msk_u20[msk_u20 < 20] <- NA
msk_u20[msk_u20 >= 20] <- 1
msk_u20 %>% saveRDS(str_c("results/R_out/mask_AGB_under20_stars",fn_suff,".rds"))

# Inverse <20 mask
mskinv_u20 <- agb.ras
mskinv_u20[mskinv_u20 < 20] <- 1
mskinv_u20[mskinv_u20 >= 20] <- 0
mskinv_u20 %>% saveRDS(str_c("results/R_out/mask_inv_AGB_under20_stars",fn_suff,".rds"))

# AGB >310 mask
val <- 310
msk_o310 <- agb.ras
msk_o310[msk_o310 > val] <- NA
msk_o310[msk_o310 <= val] <- 1
msk_o310 %>% saveRDS(str_c("results/R_out/mask_AGB_over",val,"_stars",fn_suff,".rds"))

# Inverse >310 mask
mskinv_o310 <- agb.ras
mskinv_o310[msk_o310 == 1] <- 0
mskinv_o310[mskinv_o310 > val] <- 1
mskinv_o310 %>% saveRDS(str_c("results/R_out/mask_inv_AGB_over",val,"_stars",fn_suff,".rds"))

# AGB >300 mask
val <- 300
msk_o300 <- agb.ras
msk_o300[msk_o300 > val] <- NA
msk_o300[msk_o300 <= val] <- 1
msk_o300 %>% saveRDS(str_c("results/R_out/mask_AGB_over",val,"_stars",fn_suff,".rds"))

# Inverse >300 mask
mskinv_o300 <- agb.ras
mskinv_o300[msk_o300 == 1] <- NA
mskinv_o300[mskinv_o310 > val] <- 1
mskinv_o300 %>% saveRDS(str_c("results/R_out/mask_inv_AGB_over",val,"_stars",fn_suff,".rds"))

# Inverse >275 mask
val <- 275
mskinv_o275 <- agb.ras
mskinv_o275[mskinv_o275 <= val] <- NA
mskinv_o275[mskinv_o275 > val] <- 1
names(mskinv_o275) <- 'Mask'
mskinv_o275 %>% saveRDS(str_c("results/R_out/mask_inv_AGB_over",val,"_stars",fn_suff,".rds"))
mskinv_o275 <- readRDS(str_c("results/R_out/mask_inv_AGB_over",val,"_stars",fn_suff,".rds"))

# AGB >132 mask
val <- 132
msk_o132 <- agb.ras
msk_o132[msk_o132 > val] <- NA
msk_o132[msk_o132 <= val] <- 1
msk_o132 %>% saveRDS(str_c("results/R_out/mask_AGB_over",val,"_stars",fn_suff,".rds"))

# Load masks ---------------------------------------------------------------------------------
# Previously-created from backscatter, land cover, and water features
msk_land <- # 1== ALOS 2018 land (i.e. Normal, Layover, and Shadowing)
  readRDS("results/R_out/mask_landALOS_hisp18_stars.rds") 
msk_A <- # NA== ALOS mask: non-valid ALOS pixels; 1==Normal ALOS land pixels
  readRDS("results/R_out/mask_ALOS_stars.rds") 
msk_p3 <- # NA== g0>0.3
  readRDS("results/R_out/mask_ALOSoverpt3_stars.rds")
msk_WU <- # NA== WaterUrban and ALOS ocean; 1==all other land
  readRDS("results/R_out/mask_WaterUrban_stars.rds")
msk_Aw <- # NA== ALOS mask and OSM water with 25 m buffer
  readRDS("results/R_out/mask_ALOS_OSMwater25_stars.rds")
msk_B <- # NA== LC17 Bareland
  readRDS("results/R_out/mask_Bareland_stars.rds")
mskinv_T <- # 1== LC17 tree cover
  readRDS("results/R_out/mask_inv_TreeCover_stars.rds")

# AGB-based
msk_u20 <- # NA==AGB<20 and ALOS mask
  readRDS("results/R_out/mask_AGB_under20_stars.rds")
mskinv_u20 <- # 1== where AGB<20 in valid ALOS land; 0==AGB>=20
  readRDS("results/R_out/mask_inv_AGB_under20_stars.rds")
msk_o310 <- # NA==AGB>310 and ALOS mask
  readRDS("results/R_out/mask_AGB_over310_stars.rds")
mskinv_o310 <- # 1== where AGB>310 in valid ALOS land; 0==AGB<=310
  readRDS("results/R_out/mask_inv_AGB_over310_stars.rds")

# Combined masks
msk_AWUwb <- # NA==ALOS mask, LC17 Water and Urban, and OSM water with 25 m buffer
  readRDS(str_c("results/R_out/mask_AGB_ALOS_WaterUrban_water25_stars",fn_suff,".rds"))
msk_AWU <- # 1==Valid ALOS land without WaterUrban
  readRDS("results/R_out/mask_ALOS_WaterUrban_stars.rds")
msk_Ww <- # NA==ALOS mask, LC17 water, and OSM water with 25 m buffer
  readRDS("results/R_out/mask_allwater_stars.rds")
msk_default <- # NA==ALOS mask, LC17 water and urban, OSM water 25 m buffer, AGB<20
  readRDS("results/R_out/mask_ALOS_pt3_WaterUrban_water25_AGBu20_stars.rds")

water_polysb <- # OSM water with 25 m buffer
  st_read('results/masks/vector/osm_water_buff25m.shp')
sum(units::set_units(st_area(water_polysb), ha))

# Inverse masks
# In general, inverse masks are created by multipling inverse mask by msk_A
mskinv_p3 <- # 1== g0>0.3; 0==g0<=0.3; NA=ALOS mask
  readRDS("results/R_out/mask_inv_ALOSoverpt3_stars.rds")
mskinv_p2 <- # 1== g0>0.2; 0==g0<=0.2; NA=ALOS mask
  readRDS("results/R_out/mask_inv_ALOSoverpt2_stars.rds")
mskinv_WU <- # 1==where WaterUrban overlap valid ALOS values
  readRDS("results/R_out/mask_inv_ALOS_WU_raster.rds") 
mskinv_WP <- # 1==where OSM water (no buffer) overlap valid ALOS values
  readRDS("results/R_out/mask_inv_ALOS_OSMwater_raster.rds")
mskinv_WPb <- # 1==where OSM water w/ 25 m buffer overlap valid ALOS values
  readRDS("results/R_out/mask_inv_OSMwater25mbuffer_raster.rds")
# mskinv_WUWPb <- # >0 == where WaterUrban and OSM water 25m overlap valid ALOS land
#   readRDS("results/R_out/mask_inv_WUWPb_raster.rds")

# Report starting number of NAs
df_mskA <- readRDS('results/R_out/mask_pcts.rds')
df_LCcounts <- readRDS('results/R_out/lc_ALOS_pixel_counts_tbl.rds')

# Get counts from mostly inverse masks --------------------------------------------------------
# Combine some masks
msk_WUwb <- msk_WU * msk_Ww # Mask of all water and urban
msk_AWUwb <- msk_WUwb * msk_A
msk_AWUwb %>% saveRDS(str_c("results/R_out/mask_AGB_ALOS_WaterUrban_water25_stars",fn_suff,".rds"))

msk_default <- msk_Ww * msk_WU * msk_p3 * msk_u20

# Get counts and percents of masked and unmasked pixels
cell_ct <- length(msk_A[[1]]) 
cts_df <- tibble(
    total = cell_ct,
    alos = sum(is.na(msk_A[[1]])),
    land = sum(is.na(msk_land[[1]])),
    g0_lt_p3 = sum(is.na(msk_p3[[1]])),
    water_urban = sum(is.na(msk_WU[[1]])),
    all_water = sum(is.na(msk_Ww[[1]])),
    WUwb = sum(is.na(msk_WUwb[[1]])),
    alos_WUwb = sum(is.na(msk_AWUwb[[1]])),
    agb_lt_20 = sum(is.na(msk_u20[[1]])),
    agb_gt_310 = sum(is.na(msk_o310[[1]]))
  ) %>% 
  pivot_longer(everything()) %>% 
  rename(masked = value) %>% 
  mutate(unmasked = cell_ct - masked)
land_ct <- cts_df %>% filter(name=='land') %>% select(unmasked) %>% as.numeric()
ocean_ct <- cts_df %>% filter(name=='land') %>% select(masked) %>% as.numeric()
cts_df <- cts_df %>% 
  mutate(masked_in_land = masked - ocean_ct) %>% 
  mutate(pct_of_land=masked_in_land/land_ct)
cts_df %>% saveRDS('results/R_out/df_counts_masks.rds')
cts_df <- readRDS('results/R_out/df_counts_masks.rds')

# Get percents
sum(is.na(msk_default[[1]]))
(pct = (sum(is.na(msk_default[[1]])) - ocean_ct) / land_ct)
(pct1 = (sum(is.na(msk_all1[[1]])) - ocean_ct) / land_ct)

# Look at AGB > 132
(o132 <- sum(is.na(msk_o132[[1]])))
(pct = (o132 - ocean_ct) / land_ct)

msk_def_o132 <- msk_o132 * msk_default
(pct_o132 = (sum(is.na(msk_def_o132[[1]])) - sum(is.na(msk_default[[1]]))) / land_ct)

# Look at values in Tree Cover
mskinv_T <- # 1== LC17 tree cover
  readRDS("results/R_out/mask_inv_TreeCover_stars.rds")
msk_default <- # NA==ALOS mask, LC17 water and urban, OSM water 25 m buffer, AGB<20
  readRDS("results/R_out/mask_ALOS_pt3_WaterUrban_water25_AGBu20_stars.rds")

treecover_masked <- mskinv_T * msk_default
ct_treecover <- sum(mskinv_T[[1]]==1, na.rm=T)
nanct_treecover <- sum(is.na(mskinv_T[[1]]))
nanct2_treecover <- sum(is.na(treecover_masked[[1]]))
ct_treecover/land_ct
(pct_tc = (sum(is.na(treecover_masked[[1]])) - nanct_treecover) / ct_treecover)

tc_o132 <- treecover_masked * msk_o132
(pct_tc_o132 = (sum(is.na(tc_o132[[1]])) - nanct2_treecover) / ct_treecover)

# Look at values in Grassland, Shrubs
mskinv_GS <- # 1== LC17 grassland and shrubs
  readRDS("results/R_out/mask_inv_GrasslandShrubs_stars.rds")
msk_default <- # NA==ALOS mask, LC17 water and urban, OSM water 25 m buffer, AGB<20
  readRDS("results/R_out/mask_ALOS_pt3_WaterUrban_water25_AGBu20_stars.rds")
msk_o132 <- 
  readRDS("results/R_out/mask_AGB_over132_stars.rds")

GS_masked <- mskinv_GS * msk_default
ct_GS <- sum(mskinv_GS[[1]]==1, na.rm=T)
nanct_GS <- sum(is.na(mskinv_GS[[1]]))
nanct2_GS <- sum(is.na(GS_masked[[1]]))
ct_GS/land_ct
(pct_tc = (sum(is.na(GS_masked[[1]])) - nanct_GS) / ct_GS)

GS_o132 <- GS_masked * msk_o132
(pct_GS_o132 = (sum(is.na(GS_o132[[1]])) - nanct2_GS) / ct_GS)

# Look at values in Tree Cover, Grassland, Shrubs
mskinv_GS[is.na(mskinv_GS)] <- 0
mskinv_T[is.na(mskinv_T)] <- 0
mskinv_veg <- mskinv_T + mskinv_GS
mskinv_veg[mskinv_veg==0] <- NA
mskinv_veg %>% saveRDS("results/R_out/mask_inv_LC17_vegTGS_stars.rds")

ct1 <- sum(mskinv_veg[[1]]==1, na.rm=T)
nanct <- sum(is.na(mskinv_veg[[1]]))

veg_masked <- mskinv_veg * msk_default
nanct2 <- sum(is.na(veg_masked[[1]]))
ct1/land_ct
(pct_tc = (nanct2 - nanct) / ct1)

veg_o132 <- veg_masked * msk_o132
(pct_veg_o132 = (sum(is.na(veg_o132[[1]])) - nanct2) / ct1)

# Write get_counts: Count non-masked values (msk_0!=NA) and values>threshold within land category ----
get_counts <- function(valuesraster, msk_zone, msk_0, threshold=132){
  # For given area of interest (msk_zone!=NA), 
  # count non-masked values (msk_0!=NA) and values>threshold after masking 
  #  - msk_zone - defines zone of interest ("total")
  #  - msk_0 - which pixels to count inside msk_zone ("masked")
  #  - threshold - count values > threshold after masking
  # 1. Convert rasters to DF
  df <- 
    tibble(zone=as.vector(msk_zone[[1]]),
           agb=as.vector(valuesraster[[1]]) %>% round(4), 
           mask=as.vector(msk_0[[1]])) %>% 
    filter(!is.na(zone)) 
  total = length(df[[1]])
  # Filter by mask
  df1 <- df %>% filter(!is.na(mask))
  total1 = length(df1[[1]])
  # Filter by (upper) threshold value
  df2 <- df1 %>% filter(!agb > threshold)
  total2 = length(df2[[1]])
  # Tally it up
  cts_df <- tibble(
    masked = total-total1,
    over_thresh = total1-total2,
    valid=total2
    ) %>% 
    pivot_longer(everything()) %>% 
    mutate(pct = value/sum(value))
}
df <- get_counts(agb.ras, mskinv_GS, msk_default, 132)

# Barplot
(bp <- ggplot(df, aes(x="", y=value, fill=group))+
    geom_bar(width = 1, stat = "identity")+ 
    scale_fill_manual(values=c("#E69F00", "#999999", "#999000")) +
    theme_minimal())
(pie <- bp + coord_polar("y", start=0))


# Mask out WU
msk_u20_noWU <- msk_u20 * msk_WU
msk_u20_noWU %>% saveRDS("results/R_out/mask_inv_AGBunder20_WU_stars.rds")
mskinv_u20_noWU <- # 1== AGB<20 with ALOS and WaterUrban masks applied
  readRDS("results/R_out/mask_inv_AGBunder20_WU_stars.rds")

# Combine all masks
msk_Ap3_WU_wb <- msk_WU * msk_p3 * msk_Aw 
msk_Ap3_WU_wb %>% 
  saveRDS('results/R_out/mask_ALOS_pt3_WaterUrban_water25_stars.rds')
msk_Ap3_WU_wb <- 
  readRDS('results/R_out/mask_ALOS_pt3_WaterUrban_water25_stars.rds')

msk_Ap3_WU_wb_u20 <- msk_Ap3_WU_wb * msk_u20
msk_Ap3_WU_wb_u20 %>% 
  saveRDS('results/R_out/mask_ALOS_pt3_WaterUrban_water25_AGBu20_stars.rds')


# Apply masks to AGB --------------------------------------------------------------------------------------
agb.ras <- read_stars("results/tifs_by_R/agb18_v1_l0.tif")
msk_Ap3_WU_wb_u20 <- 
  readRDS('results/R_out/mask_ALOS_pt3_WaterUrban_water25_AGBu20_stars.rds')

agb_mskd <- agb.ras * msk_Ap3_WU_wb
agb_mskd %>% as("Raster") %>% 
  writeRaster('results/tifs_by_R/agb18_v1_l1_mask_Ap3WUw25.tif')

agb_mskd <- agb.ras * msk_Ap3_WU_wb_u20
agb_mskd %>% as("Raster") %>% 
  writeRaster('results/tifs_by_R/agb18_v1_l1_mask_Ap3WUw25_u20.tif')

agb_mskd2 <- agb_mskd * msk_o310
agb_mskd2 %>% as("Raster") %>% 
  writeRaster('results/tifs_by_R/agb18_v1_l1_mask_Ap3WUw25_u20o310.tif')
agb_mskd2 <- read_stars('results/tifs_by_R/agb18_v1_l1_mask_Ap3WUw25_u20o310.tif')

# Apply value cap to AGB ----#######################################################
agb_sat <- agb_mskd
agb_sat <- read_stars('results/tifs_by_R/agb18_v1_l1_mask_Ap3WUw25_u20.tif')
saturation_pt <- 132

# Apply saturation point to AGB and save
agb_sat[agb_sat > saturation_pt] <- saturation_pt
agb_sat %>% saveRDS('results/R_out/agb18_v1_l2_mask_Ap3WUw25_u20_cap132.rds')
agb_sat %>% as("Raster") %>% 
  writeRaster('results/tifs_by_R/agb18_v1_l2_mask_Ap3WUw25_u20_cap132.tif')

# PLOT using tmap (interactive) ----############################################
# Load contextual data for mapping
parks_hti <- # OSM water with 25 m buffer
  st_read('data/contextual_data/OSM_free/hti_nature_reserves_osm.shp')

# Make crop bounding box/extent for testing
# extent for raster
test_ext <- extent(-72.7, -72.5, 18.2, 18.35)
# bounding box for stars
test_bb <- st_bbox(c(xmin=-72.34, xmax=-72.17, ymin=18.3, ymax=18.38), crs=4326) %>% 
  st_as_sfc()

# Set tmap to interactive viewing (vs. "plot")
tmap_mode("view")

# View some masks
names(mskinv_o275) <- 'Mask'
tm_shape(agb_mskd2[test_bb]) + tm_raster() +
  tm_shape(mskinv_o275[test_bb]) + tm_raster(col="Mask", palette="red") +
  tm_shape(parks_hti) + tm_borders()

# Layers
lyr_g0 <- tm_shape(crop(g0, test_ext)) + 
  tm_raster(style="order",
            # breaks=seq(0, 0.1, 0.02), 
            palette=palette(hcl.colors(8, "viridis")))
lyr_water <- 
  tm_shape(water_polys) + 
  tm_borders() + 
  tm_fill(col="cyan") 
lyr_waterB <- tm_shape(water_polysb) + tm_borders() 

# Map
(map_g0w <- lyr_g0 + lyr_water + lyr_waterB)

# Look at AGB distributions ---- ###################################################
agb.br <- brick(raster("results/agb/agb_2018_v6_mask2share.tif"),
                raster("results/agb/agb_2018_v6CI_2share.tif"))
names(agb.br) <- c('AGB', 'CI')
saveRDS(agb.br, file = "results/R_out/AGB_95ci.rds")

get_brick_stats <- function(lc.br){
  # Get selection of percentiles for each LC
  agb.qs <- data.frame(row.names=c('0%', '1%','2%','10%', '25%', '50%', '75%', '90%','98%', '99%','100%'))
  for (i in seq(1,nlayers(lc.br))){
    lc <- names(lc.br[[i]])
    agb.qs[[lc]] <- quantile(lc.br[[i]], probs=c(0, 0.01, 0.02, 0.1, 0.25, 0.5, 0.75, 0.9,0.98, 0.99, 1))
  }
  # Get means, SDs, and Skews
  stats <- list(mean=cellStats(lc.br, stat='mean', na.rm=TRUE),
                sd=cellStats(lc.br, stat='sd', na.rm=TRUE),
                skew=cellStats(lc.br, stat='skew', na.rm=TRUE)) %>%
    bind_cols()%>%
    t()
  colnames(stats) <- names(lc.br)
  agb.stats <- rbind(agb.qs, stats)
}
agb.stats <- get_brick_stats(agb.br)

# Graphing ---- ###################################################################
# Sample the distribution of values in the raster
agb.samp <- agb %>% 
  sampleRandom(100000, na.rm=TRUE) %>% 
  as.data.frame() %>% 
  rename(AGB='.')

# Plot density
p <- ggplot(agb.samp, aes(x=AGB)) + 
  geom_histogram(aes(y=..density..), fill="#69b3a2", color="#e9ecef", alpha=0.7, binwidth = 2.5)+
  geom_density(alpha=.2) + 
  scale_x_continuous(name = expression(paste("Aboveground biomass (Mg ha"^"-1", ")")),
                     breaks = seq(0, 150, 25),
                     limits=c(-5, 150)) +
  geom_vline(aes(xintercept=mean(AGB)), color="black", linetype="dashed", size=.5)+
  geom_vline(aes(xintercept=median(AGB)),
             color="black", linetype="dashed", size=.5)+
  theme_minimal()
p

mean(agb.samp$AGB)
median(agb.samp$AGB)
