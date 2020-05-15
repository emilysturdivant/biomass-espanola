#library(silvr)
library(readr)
library(raster)
library(tidyverse)
library(ggridges)
library(rgdal)
library(tmap)
library(rgdal)
library(here)
library(gdalUtils)
library(here)

# Merge Gamma0 backscatter tiles produced by biota ----
year <- 2018
mosaic_rasters(list.files(path='results/g0nu_HV/Gamma0_lee', pattern='.tif', full.names = TRUE), 
               str_c('results/g0nu_HV/g0nu_', year, '_HV_lee.tif'), 
               projwin=c(-74.48133, 20.09044, -68.32267, 17.47022))

# Merge ALOS mosaic rasters ----
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
  # l <- vector(mode="logical", length(rs))
  # for(i in seq(1, length(rs))){
  #   if(class(rs[[i]])=='try-error') l[i] <- TRUE
  # }
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
isl_poly <- sf::st_read("results/masks/vector/hisp18_maskLand_clean.shp")
sum(area(isl_poly))/1000000



plot(dn_mask)
dn_date <- merge.alos.tiles('data/ALOS', '_18_date_F02DAR$', 
                            isl_poly, 'data/R_out/hisp18_date.tif')
# writeRaster(dn_date, 'data/R_out/hisp18_date.tif')
plot(dn_date)
dn_linci <- merge.alos.tiles('data/ALOS', '_18_linci_F02DAR$', 
                             isl_poly, 'data/R_out/hisp18_localincidence.tif')
# writeRaster(dn_linci, 'data/R_out/hisp18_localincidence.tif')
dn_mask <- merge.alos.tiles('data/ALOS', '_18_mask_F02DAR$', 
                            isl_poly, 'data/R_out/hisp18_mask.tif')
dn_mask[dn_mask<75] <- NA
dn_mask[!is.na(dn_mask)] <- 1
writeRaster(dn_mask, 'results/masks/hisp18_maskLand.tif')
mask_land <- raster('results/masks/hisp18_maskLand.tif')

# Load previously created mosaics
dn_linci <- raster('results/tifs_by_R/hisp18_localincidence.tif')
dn_mask <- raster('results/tifs_by_R/hisp18_mask.tif')
dn_date <- raster('results/tifs_by_R/hisp18_date.tif')
plot(dn_linci)
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

# Barplot
(bp <- ggplot(df, aes(x="", y=value, fill=group))+
    geom_bar(width = 1, stat = "identity")+ 
    scale_fill_manual(values=c("#E69F00", "#999999")) +
    theme_minimal())
(pie <- bp + coord_polar("y", start=0))

# Crop backscatter to Hispaniola extent and replace 0 with NA ---- ############################
# Crop backscatter to the common Hispaniola extent
gdal_translate("results/g0nu_HV/g0nu_2018_HV.tif", 
               "results/g0nu_HV/g0nu_2018_HV_crop.tif", 
               projwin=st_bbox(msk_AWU)[c(1, 4, 3, 2)])
# Load backscatter and set 0 to NA
g0 <- read_stars("results/g0nu_HV/g0nu_2018_HV_crop.tif")
g0[g0==0] <- NA
g0 %>% write_stars("results/g0nu_HV/g0nu_2018_HV.tif")


# Crop backscatter to Haiti extent ----
# Load files
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
