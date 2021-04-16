# *************************************************************************************************
# Script to:
#     * Process ALOS mosaic tiles 
#       - load mosaic supplemental rasters, merge tiles, 
#       - make masks
# Preceeds:
#     * calculate_AGB.R - calculates AGB by plot from field data
# Requires:
#     * downloaded ALOS mosaic tiles (currently performed with biota in command line)
#     * study area polygons
#
# *************************************************************************************************

# Load libraries
#library(silvr)
library(readr)
library(raster)
library(tidyverse)
library(ggridges)
library(rgdal)
library(tmap)
library(here)
library(gdalUtils)
library(stars)

# Merge ALOS backscatter tiles (Gamma0) produced by biota --------------------------------------------
year <- 2018
mosaic_rasters(list.files(path='results/g0nu_HV/Gamma0_lee', pattern='.tif', full.names = TRUE), 
               str_c('results/g0nu_HV/g0nu_', year, '_HV_lee.tif'), 
               projwin=c(-74.48133, 20.09044, -68.32267, 17.47022))

# Merge supplemental ALOS mosaic rasters ------------------------------------------------------------
merge.alos.tiles <- function(path, pattern, clip_poly, filename=FALSE){
  # List filenames of ALOS mask tiles. Tiles range from 18-21, 68-75. ENVI format.
  fps <- list.files(path=path,
                    pattern=pattern,
                    full.names=TRUE,
                    recursive=TRUE,
                    include.dirs=FALSE)
  fps <- fps[1:2]
  # Load and crop rasters and store in list
  load.crop.na <- function(fp, ext) try({
    r <- raster(fp); 
    r <- crop(r, ext)
  }, silent=TRUE)
  rs <- fps %>% lapply(load.crop.na, ext=clip_poly)
  # Remove any that threw error (don't overlap with the polygon)
  l <- lapply(function(r) if(class(r) == 'try-error') {l <- TRUE} else {l <- FALSE})
  # Merge the tiles
  dn <- do.call(merge, rs[!l])
  if(filename) writeRaster(dn, filename)
  return(dn)
}
# Load island polygons
isl_poly <- readOGR(dsn="data/contextual_data/Hispaniola", layer='Hisp_adm0')
isl_poly <- readOGR(dsn=here("data", "Hispaniola"), layer='Hisp_adm0')
isl_poly <- readOGR(dsn="results/masks/vector", layer='hisp18_maskLand_clean')
sum(area(isl_poly))/1000000 

# Load polygon
isl_poly <- sf::st_read("results/masks/vector/hisp18_maskLand_clean.shp")

# get polygon area
sum(units::set_units(st_area(isl_poly), ha))

hti_poly <- st_read("data/contextual_data/HTI_adm/HTI_adm0.shp") %>% 
  st_make_valid(hti_poly)
hti_poly %>% st_write("data/contextual_data/HTI_adm/HTI_adm0_fix.shp", append=F)





dn_date <- merge.alos.tiles('data/ALOS', '_18_date_F02DAR$', 
                            isl_poly, 'data/R_out/hisp18_date.tif')
# writeRaster(dn_date, 'data/R_out/hisp18_date.tif')
plot(dn_date)
dn_linci <- merge.alos.tiles('data/ALOS', '_18_linci_F02DAR$', 
                             isl_poly, 'data/R_out/hisp18_localincidence.tif')
# writeRaster(dn_linci, 'data/R_out/hisp18_localincidence.tif')
dn_mask <- merge.alos.tiles('data/ALOS', '_18_mask_F02DAR$', 
                            isl_poly, 'data/R_out/hisp18_mask.tif')
plot(dn_mask)
mask_land <- dn_mask
mask_land[mask_land<75] <- NA
mask_land[!is.na(mask_land)] <- 1
writeRaster(mask_land, 'results/masks/hisp18_maskLand.tif')
mask_land <- raster('results/masks/hisp18_maskLand.tif')

# Mask to Haiti extent
msk_land <- terra::rast('results/masks/hisp18_maskLand.tif')
hti_poly <- terra::vect("data/contextual_data/HTI_adm/HTI_adm0_fix.shp")
msk_land <- msk_land %>% 
  crop(hti_poly) %>% 
  mask(hti_poly, filename='results/masks/hti18_maskLand.tif', overwrite=T)

# Load previously created mosaics
dn_linci <- raster('results/tifs_by_R/hisp18_localincidence.tif')
dn_mask <- raster('results/tifs_by_R/hisp18_mask.tif')
dn_date <- raster('results/tifs_by_R/hisp18_date.tif')
plot(dn_linci)
plot(dn_mask)
plot(dn_date)
plot(isl_poly, add=TRUE)

# Trying to use tmap, but the files are too big.
lyr_linci <- tm_shape(dn_linci) + tm_raster()
lyr_mask <- tm_shape(dn_mask) + tm_raster()
lyr_date <- tm_shape(dn_date) + tm_raster()
lyr_isl <- tm_shape(isl_poly) + tm_borders() 
tmap_mode("plot")
(map_linci <- lyr_linci + lyr_isl)
map_all <- tmap_arrange(lyr_linci + lyr_isl, 
                        lyr_date + lyr_isl, 
                        lyr_date + lyr_isl)
map_all

# STARS - Look at proportion of values in each mask category 
dn_mask <- read_stars('results/tifs_by_R/hisp18_mask.tif')
rvals <- dn_mask[[1]]
df <- 
  tibble(group = c("Normal", "Layover", "Shadowing"), 
         count = c(sum(rvals==255, na.rm=TRUE), 
                   sum(rvals==100, na.rm=TRUE),
                   sum(rvals==150, na.rm=TRUE))) %>% 
  mutate(pct= count / sum(count)) %>%
  add_row(group='Land', count=sum(count))
saveRDS(df, 'results/R_out/mask_pcts.rds')
df <- readRDS('results/R_out/mask_pcts.rds')

# Replicate ALOS Normal mask
msk_A <- read_stars('results/tifs_by_R/hisp18_mask.tif')
msk_A[msk_A<254] <- NA
msk_A[!is.na(msk_A)] <- 1
msk_A %>% saveRDS("results/R_out/mask_ALOS_raster.rds")
msk_A %>% as("Raster") %>% writeRaster("results/masks/mask_ALOS18.tif")

# Barplot
(bp <- ggplot(df, aes(x="", y=value, fill=group))+
    geom_bar(width = 1, stat = "identity")+ 
    scale_fill_manual(values=c("#E69F00", "#999999")) +
    theme_minimal())
(pie <- bp + coord_polar("y", start=0))

# Crop backscatter to Hispaniola extent and replace 0 with NA ---- ############################
# Crop backscatter to the common Hispaniola extent
g0_fp <- "results/g0nu_HV/g0nu_2018_HV.tif"
gdal_translate("results/g0nu_HV/g0nu_2018_HV_uncropped.tif", 
               g0_fp, 
               projwin=st_bbox(msk_AWU)[c(1, 4, 3, 2)])
# Load backscatter and set 0 to NA
g0 <- read_stars("results/g0nu_HV/g0nu_2018_HV_crop.tif")
g0[g0==0] <- NA
g0 %>% write_stars(g0_fp)

# Aggregate to 50m, as recommended by Saatchi 2015 and performed by Michelakis et al. 2015
g0.nofilt <- raster(g0_fp)
sum(is.na(g0.nofilt))
g0.agg <- aggregate(g0.nofilt, fact=2, fun=mean, na.rm=TRUE,
                    filename="results/g0nu_HV/g0nu_2018_agg50m.tif",
                    overwrite=TRUE)

# # Perform SAGA Lee filter ------------------------------------------------------
# g0.nofilt <- terra::rast(g0_fp)
# 
# install.packages('RSAGA')
# library(RSAGA)
# 
# work_env <- rsaga.env(workspace = getwd())
# rsaga.get.modules('grid_filter')
# rsaga.get.usage('grid_filter', 3)
# 
# rsaga.import.gdal(g0_fp)

# Crop WorldClim rasters -------------------------------------------------------------------
bb <- st_bbox(g0) #te=c(xmin,ymin,xmax,ymax)
srclist <- list.files(path='data/WorldClim/monthly_precip_2010_2018', 
                      pattern='.tif', full.names = TRUE)
dstlist <- file.path('data', 'WorldClim', 'monthly_precip_hisp', basename(srclist))

for (i in seq_along(srclist)) {
  gdalwarp(srcfile=srclist[[i]], 
           dstfile=dstlist[[i]], 
           te=bb, 
           overwrite=T)
}

# Make masks (Hispaniola extent) ------------------------------------------------------------
# Load LC17 masked to ALOS land 
lc <- readRDS("results/R_out/LC17_masked_to_ALOS_land_stars.rds")
lc <- raster("data/LULC/Hisp_2017_resALOS_mskLand.tif")

# WaterUrban mask
msk_WU <- lc
msk_WU[msk_WU<3] <- NA
msk_WU[!is.na(msk_WU)] <- 1
msk_WU %>% writeRaster("results/masks/mask_WaterUrban_raster.tif")
msk_WU %>% saveRDS("results/R_out/mask_WaterUrban_raster.rds")

# Water mask
msk_W <- lc
msk_W[msk_W==1] <- NA
msk_W[!is.na(msk_W)] <- 1
msk_W %>% writeRaster("results/masks/mask_Water_raster.tif")
msk_W %>% saveRDS("results/R_out/mask_WaterLC17_raster.rds")

# Mask out all but forest
mskinv_T <- lc
mskinv_T[mskinv_T!=4] <- NA
mskinv_T[mskinv_T==4] <- 1
mskinv_T %>% writeRaster("results/masks/mask_inv_TreeCover.tif")
mskinv_T %>% saveRDS("results/R_out/mask_inv_TreeCover_raster.rds")

# Mask out Bareland
msk_B <- lc
msk_B[msk_B!=3] <- 1
msk_B[msk_B==3] <- NA
msk_B %>% saveRDS("results/R_out/mask_Bareland_raster.rds")

# Mask out all but grassland and shrubs
mskinv_GS <- lc
mskinv_GS[mskinv_GS<5] <- NA
mskinv_GS[mskinv_GS>4] <- 1
mskinv_GS %>% saveRDS("results/R_out/mask_inv_GrasslandShrubs_stars.rds")

# Mask out all but tree cover, grassland and shrubs
mskinv_GS <- lc
mskinv_GS[mskinv_GS<4] <- NA
mskinv_GS[mskinv_GS>3] <- 1
mskinv_GS %>% saveRDS("results/R_out/mask_inv_LC17_vegTreeCGrasslandShrubs_stars.rds")

# Presence (Inverse mask) of LC17 Water
mskinv_W <- lc
mskinv_W[mskinv_W!=1] <- NA
mskinv_W %>% saveRDS("results/R_out/mask_inv_ALOS_Water_raster.rds")

# Combine ALOS and WaterUrban masks
msk_A <- readRDS("results/R_out/mask_ALOS_stars.rds")
msk_A %>% as("Raster") %>% writeRaster("results/masks/mask_ALOS18.tif")

rm(list=ls()) 

msk_A <-  # NA== ALOS mask: non-valid ALOS pixels; 1==Normal ALOS land pixels # HISPANIOLA
  raster("results/masks/mask_ALOS18.tif")
msk_WU <- # NA== WaterUrban and ALOS ocean; 1==all other land
  raster("results/masks/mask_WaterUrban_raster.tif")

msk_AWU <- msk_WU * msk_A
msk_AWU %>% writeRaster("results/masks/mask_ALOS_WaterUrban.tif")
msk_AWU %>% saveRDS("results/R_out/mask_ALOS_WaterUrban_raster.rds")

# Presence (Inverse mask) of LC17 Water/Urban with ALOS mask applied
msk_WU <- lc
msk_WU[msk_WU==2] <- 1 # Water is already 1 and now urban is as well
msk_WU[msk_WU!=1] <- NA
mskinv_WU <- msk_WU*msk_A
mskinv_WU %>% # 1==where WaterUrban overlap valid ALOS values
  saveRDS("results/R_out/mask_inv_ALOS_WU_raster.rds") 

# Create water mask from OSM polygons with 25 m buffer
# Create OSM water with 25 m buffer
st_read('data/contextual_data/OSM_free/gis_osm_water_a_free_1.shp') %>% 
  st_transform(32618) %>% 
  st_buffer(dist = 25) %>% 
  summarize() %>% 
  st_transform(4326) %>% 
  st_write('results/masks/vector/osm_water_buff25m.shp', append=FALSE)
water_polysb <- st_read('results/masks/vector/osm_water_buff25m.shp')

water_polysb <- st_read('results/masks/vector/osm_water_buff25m.shp')
msk_wb = msk_AWU
values(msk_wb) <- 1
msk_wb <- msk_wb %>% 
  mask(water_polysb, inverse=TRUE)
names(msk_wb) <- 'Mask'
msk_wb %>% writeRaster("results/masks/mask_osm_water_buff25m.tif")

# Mask OSM water 25 (add to ALOS mask)
msk_A <-  # NA== ALOS mask: non-valid ALOS pixels; 1==Normal ALOS land pixels # HISPANIOLA
  raster("results/masks/mask_ALOS18.tif")
msk_Aw <- msk_A %>% 
  mask(water_polysb, inverse=TRUE)
msk_Aw %>% saveRDS("results/R_out/mask_ALOS_OSMwater25_raster.rds")

# Water from OSM polygons
water_polys <- st_read('data/contextual_data/OSM_free/gis_osm_water_a_free_1.shp')
mskinv_WP <- msk_A %>% 
  mask(water_polys, inverse=FALSE)
mskinv_WP %>% saveRDS("results/R_out/mask_inv_ALOS_OSMwater_raster.rds")

# Water from OSM polygons with 25 m buffer
mskinv_WPb <- msk_A %>% 
  mask(water_polysb, inverse=FALSE)
mskinv_WPb %>% saveRDS("results/R_out/mask_inv_OSMwater25mbuffer_raster.rds")


# View masks
tm_shape(msk_wb[test_bb]) + tm_raster(col='Mask', palette='red') +
  tm_shape(water_polysb) + tm_borders()


# Create inverse mask of WU and OSM water 25 m
mskinv_WU <- mskinv_WU %>% as("Raster")
mskinv_WPb <- mskinv_WPb %>% st_as_stars()
mskinv_WUWPb <- mskinv_WU + mskinv_WPb
mskinv_WUWPb %>% saveRDS("results/R_out/mask_inv_WUWPb_raster.rds")

# Mask all water
msk_W <- readRDS("results/R_out/mask_WaterLC17_raster.rds")
msk_wb <- raster("results/masks/mask_osm_water_buff25m.tif")
msk_Ww <- msk_W * msk_wb 
msk_Ww %>% saveRDS("results/R_out/mask_allwater_raster.rds")

# Crop for testing
bb <- st_bbox(c(xmin=-72.68, xmax=-72.51, ymin=18.22, ymax=18.32)) %>% 
  st_as_sfc()
st_crs(bb) <- 4326
msk_Ac <- msk_A[bb]

# Create mask of backscatter >0.3
msk_p3 <- # Initialize mask
  read_stars("results/g0nu_HV/g0nu_2018_HV.tif")
msk_p3[msk_p3>0.3] <- NA
msk_p3[msk_p3<=0.3] <- 1
msk_p3 %>% saveRDS("results/R_out/mask_ALOSoverpt3_stars.rds")

# Create inverse mask of backscatter >0.3
mskinv_p3 <- # Initialize mask
  read_stars("results/g0nu_HV/g0nu_2018_HV.tif")
mskinv_p3[mskinv_p3<=0.3] <- 0
mskinv_p3[mskinv_p3>0.3] <- 1
mskinv_p3 %>% saveRDS("results/R_out/mask_inv_ALOSoverpt3_stars.rds")
(p3 <- sum(mskinv_p3[[1]]==1, na.rm=TRUE))

# Create inverse mask of backscatter >0.29
mskinv_p29 <- # Initialize mask
  read_stars("results/g0nu_HV/g0nu_2018_HV.tif")
mskinv_p29[mskinv_p29<=0.29] <- NA
mskinv_p29[mskinv_p29>0.29] <- 1
mskinv_p29 %>% saveRDS("results/R_out/mask_inv_ALOSoverpt29_stars.rds")
(p29 <- sum(mskinv_p29[[1]]==1, na.rm=TRUE))

# Create inverse mask of >0.2
mskinv_p2 <- g0; rm(g0)
mskinv_p2[mskinv_p2<=0.2] <- 0
mskinv_p2[mskinv_p2>0.2] <- 1
mskinv_p2 %>% saveRDS("results/R_out/mask_inv_ALOSoverpt2_stars.rds")
(p2 <- sum(mskinv_p2[[1]]==1, na.rm=TRUE))

# PLOT using tmap (interactive) ----############################################
test_ext <- extent(-72.7, -72.5, 18.2, 18.35)
tmap_mode("view")

# View masks
tm_shape(crop(msk_W, test_ext)) + tm_raster() +
  tm_shape(crop(msk_Aw, test_ext)) + tm_raster() +
  tm_shape(water_polysb) + tm_borders()

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
lyr_waterlc <- tm_shape(msk_W) + tm_raster()

# Map
(map_g0w <- lyr_g0 + lyr_waterlc + lyr_water + lyr_waterB)

# PLOT with ggplot ----############################################
g0c %>% 
  rasterToPoints() %>% 
  data.frame() %>% 
  mutate(cuts = cut(backscatter, breaks=seq(0, 0.1, 0.02))) %>% 
  
  ggplot() + 
  geom_tile(aes(x=x,y=y,fill=cuts)) + 
  scale_fill_brewer(type = "seq", palette = "YlGn") +
  coord_equal() +
  theme_minimal() +
  theme(panel.grid.major = element_blank()) +
  xlab("Longitude") + ylab("Latitude")
# + geom_sf(data=water_polys, fill='blue') 

g0m %>% 
  # convert raster to point dataframe for ggplot
  rasterToPoints() %>% 
  data.frame() %>% 
  mutate(cuts = cut(backscatter, breaks=seq(0, 0.1, 0.02))) %>% 
  # use geom_tile to plot DF
  ggplot() + 
  geom_tile(aes(x=x, y=y, fill=cuts)) + 
  scale_fill_brewer(type = "seq", palette = "YlGn") +
  coord_equal() +
  # theme_bw() +
  theme(panel.grid.major = element_blank()) +
  xlab("Longitude") + ylab("Latitude")

# Look at biota Lee filter ----###########################################
g0 <- read_stars("results/g0nu_HV/g0nu_2018_HV.tif")
g0lee <- read_stars('results/g0nu_HV/Gamma0_lee/Gamma0_2018_N19W075.tif')
g0agg <- read_stars('results/g0nu_HV/g0nu_2018_haiti_agg50m.tif')

mskinv_p3 <- readRDS("results/R_out/mask_inv_ALOSoverpt3_stars.rds")
mskinv_p29 <- readRDS("results/R_out/mask_inv_ALOSoverpt29_stars.rds")
names(mskinv_p29) <- 'Mask'
names(mskinv_p3) <- 'Mask'

# bounding box to crop for testing
test_bb <- st_bbox(c(xmin=-74.05, xmax=-74, ymin=18.32, ymax=18.38)) %>% 
  st_as_sfc()
st_crs(test_bb) <- 4326
parks_hti <- # Parks for context
  st_read('data/contextual_data/OSM_free/hti_nature_reserves_osm.shp')

# Look at map (small area)
tmap_mode("view")
tm_shape(g0lee[test_bb]) + tm_raster(breaks=seq(0, 0.3, 0.01), palette=palette(hcl.colors(8, "viridis"))) +
  tm_shape(g0agg[test_bb]) + tm_raster(breaks=seq(0, 0.3, 0.01), palette=palette(hcl.colors(8, "viridis"))) +
  tm_shape(g0[test_bb]) + tm_raster(breaks=seq(0, 0.3, 0.01), palette=palette(hcl.colors(8, "viridis"))) +
  tm_shape(mskinv_p3[test_bb]) + tm_raster(col="Mask", palette=c("red", "white")) +
  tm_shape(mskinv_p29[test_bb]) + tm_raster(col="Mask", palette="red") +
  tm_shape(parks_hti) + tm_borders()
tm_shape(g0lee[test_bb]) + tm_raster(breaks=seq(0, 0.3, 0.01), palette=palette(hcl.colors(8, "viridis"))) +
  tm_shape(g0[test_bb]) + tm_raster(breaks=seq(0, 0.3, 0.01), palette=palette(hcl.colors(8, "viridis")))


# Crop backscatter to Haiti extent ----###########################################
hti_poly <- st_read("data/contextual_data/HTI_adm/HTI_adm0_fix.shp")
units::set_units(st_area(hti_poly), ha)
r <- raster("results/g0nu_HV/g0nu_2018_HV_haitiR.tif") %>% 
  mask(hti_poly, inverse=FALSE) 
r %>% writeRaster("results/g0nu_HV/g0nu_2018_HV_hticlip.tif", overwrite=T)
g0 <- raster("results/g0nu_HV/g0nu_2018_HV_hticlip.tif")

msk_land <- raster("results/masks/hti18_maskLand_clip2border.tif")
extent(msk_land)
extent(g0)

# View
test_ext <- extent(-71.8, -71.65, 18.22, 18.32)
tm_shape(crop(msk_land, test_ext)) + tm_raster() +
  tm_shape(crop(g0, test_ext)) + tm_raster() +
  tm_shape(hti_poly) + tm_borders()

# Load files
g0 <- read_stars("results/g0nu_HV/g0nu_2018_HV.tif")
g0 <- raster("results/g0nu_HV/g0nu_2018_HV.tif")

hti_poly <- readOGR(dsn="data/contextual_data/HTI_adm", layer='HTI_adm0')
isl_poly <- readOGR(dsn="data/contextual_data/Hispaniola", layer='Hisp_adm0')

# Crop backscatter and reclass 0s to NA
g0 <- crop(g0, hti_poly)
g0[g0==0] <- NA
writeRaster(g0, "results/g0nu_HV/g0nu_2018_HV_haitiR.tif", overwrite=TRUE)
g0 <- raster("results/g0nu_HV/g0nu_2018_HV_haitiR.tif")

# OLD: Count inland NA value in biota output ----
# Convert island to land mask
tmp_isl <- rasterize(isl_poly, g0, field=1) # land==1; sea==NA

# Multiply backscatter classes x polygon classes 
# reclass to (NA==NA; 1==land; 2==sea)
msk <- overlay(g0, tmp_isl, 
               fun=function(r1, r2){
                 r1[!is.na(r1)] <- 0 # values in land-->0
                 r1[is.na(r1)] <- 1 # NAs in land-->1
                 r2[is.na(r2)] <- 2 # sea-->2 (land==1)
                 r <- r1*r2 # create output (sea==2)
                 r[r==1] <- NA # NAs in land --> NA
                 r[r==0] <- 1 # valid land --> 1
                 return(r)}
)
plot(msk) 
# define as categorical variable
f <- as.factor(msk) # or ratify(msk) ?
x <- levels(f)[[1]]
x$code <- c("land", "sea")
levels(f) <- x
writeRaster(f, "results/masks/mask_NAinland.tif")

# Land mask
land <- msk
land[is.na(land)] <- 1
land[land==2] <- NA
writeRaster(land, "results/masks/mask_land18.tif")
land <- raster("results/masks/mask_land18.tif")
land_poly <- rasterToPolygons(land) # very slow, much faster with gdal_polygonize in QGIS 
land_poly <- readOGR(dsn="results/masks/vector", layer='mask_land18_hti')
land_area <- sum(area(land_poly)) / 1000000

# Get values and count NAs
vals <- getValues(msk)
df <- data.frame(
  group = c("Land", "Nulls"), 
  value = c(sum(vals==1, na.rm=TRUE), 
            sum(is.na(vals)))
)
na_pct <- df$value[2] / (df$value[1] + df$value[2])

# Convert to area
msk[is.na(msk)] <- 0
cell_size <- area(msk, na.rm=TRUE, weights=FALSE)
cell_size <- cell_size[!cell_size==2]
land_area<-length(cell_size)*median(cell_size)
cell_size <- cell_size[is.na(cell_size)]
NA_area<-length(cell_size)*median(cell_size)

# Get counts
vals <- getValues(mask_ras)
na_ct <- sum(is.na(vals))
nonan_ct <- sum(!is.na(vals))
sea_ct <- sum(vals==0, na.rm=TRUE)
valid_ct <- sum(vals==1, na.rm=TRUE)

sea_ct + valid_ct == nonan_ct
na_ct/valid_ct
