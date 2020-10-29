# ---------------------------------------------------------------------------------------------
# Script to:
#     * Fill gaps in AGB map by interpolating by land cover class 
# Proceeds:
#     * (possibly postprocess_AGB_map.R)
#     * regression_AGB-g0.R - creates AGB map
#     * calculate_AGB.R - calculates AGB by plot from field data
#     * process_ALOS_tiles.R
# Requires:
#     * AGB map
#     * LULC map
# ---------------------------------------------------------------------------------------------

library(sf)
library(terra)
library(tidyverse)
library(tmap)
tmap_mode('view')


lc_fps <- c("data/LULC/Haiti2017_Clip.tif", 
            "data/LULC/DR_2017_clip.tif")
lc_fp_out <- "data/LULC/Hisp_2017_resALOS_terra.tif"
lc_out <- "data/LULC/Hisp_2017_resALOS_mskLand.tif"
agb_fp <- 'results/tifs_by_R/agb18_v3_l1_mask_Ap3WUw25_u20_hti_qLee.tif'
msk_lnd_fp <- 'results/masks/hti18_maskLand.tif'
agb_from_lc_fp <- 'results/tifs_by_R/LCpatches_agb18_v3l1_Ap3WUw25u20_hti.tif'

# Load AGB ================================================================
agb_ras <- terra::rast(agb_fp)

# Merge and resample Haiti and DR land cover to AGB res ------------------------
rs <- sapply(lc_fps, 
             function(fp) {
               r <- terra::rast(fp)
               r <- terra::resample(r, agb_ras, method='near')
               return(r)}
)

# Merge, i.e. cover
lc2 <- terra::cover(rs$`data/LULC/Haiti2017_Clip.tif`, 
                    rs$`data/LULC/DR_2017_clip.tif`, 
                    filename=lc_fp_out, overwrite=T)

# Create AGB surface from mean AGB for each LC class (6 values) ================
# Load and crop LC
lc <- terra::rast(lc_fp_out) %>% 
  terra::crop(agb_ras) 

# Mask LC to land pixels
lc <- lc * terra::rast(msk_lnd_fp)

# Compute mean AGB for each LC (zonal statistics) 
zonal_stats <- terra::zonal(agb_ras, lc, 'mean', na.rm=T)

# Reclass LC to AGB values -> AGB by LC surface
agb_by_lc <- terra::classify(lc, zonal_stats)
plot(agb_by_lc)

# Fill gaps in AGB with AGB by LC ----------------------------------------------
agb_filled_fp <- 'results/tifs_by_R/agb18_v3_l1_Ap3WUw25u20_hti_filled_6zones.tif'
agb_filled <- terra::cover(agb_ras, agb_by_lc, filename=agb_filled_fp)

plot(agb_filled)
plot(agb_ras)

# Get mean AGB for each land cover patch (extract) =============================
# Polygonize LULC raster ----
# output filenames
lc_multipols_fp <- "data/LULC/Haiti2017_Clip_multipolys"
lc_pols_fp <- "data/LULC/Haiti2017_Clip_polys.geojson"

# Polygonize
lc_fps[[1]] %>% 
  terra::rast() %>% 
  terra::as.polygons() %>% 
  terra::writeVector(lc_multipols_fp)

# Un-dissolve: Convert multipolygons to polygons
lc_sf <- lc_multipols_fp %>% 
  sf::st_read() %>% 
  rename('LC' = 1) %>% 
  sf::st_cast('POLYGON')

# Save
lc_sf %>% st_write(lc_pols_fp, delete_dsn=T)
lc_sf <- st_read(lc_pols_fp)

# Extract AGB values ----
lc_sf <- lc_sf %>% 
  # Don't use Water and Urban class as they are already masked from AGB 
  filter(LC > 2)

agb_ex <- lc_sf %>% 
  # Convert to SpatVector and extract AGB values
  terra::vect() %>% 
  terra::extract(agb_ras, .) 

agb <- agb_ex %>% 
  # Convert matrix to tibble and get mean and count for each polygon
  as_tibble() %>% 
  group_by(ID) %>% 
  summarise(mean = mean(agb18_v3_l1_mask_Ap3WUw25_u20_hti_qLee, na.rm=T),
            ct = sum(!is.na(agb18_v3_l1_mask_Ap3WUw25_u20_hti_qLee)), 
            sd = sd(agb18_v3_l1_mask_Ap3WUw25_u20_hti_qLee, na.rm = T))

# Append columns to SF polys
lc_sf_all <- agb %>% 
  select(-ID) %>% 
  cbind(lc_sf, .) %>% 
  filter(ct > 1)

# Save
lc_pols_agb_fp <- "data/LULC/Haiti2017_Clip_polys_meanAGB.geojson"
lc_sf_all %>% st_write(lc_pols_agb_fp, delete_dsn = T)

# Look
agb_raster <- raster::raster(agb_fp)
tm_shape(lc_sf_all) + tm_fill(col='agb.mean') +
  tm_shape(agb_raster) + tm_raster()

# Fill missing AGB with patch means --------------------------------------------
agb_from_lc_sd_fp <- 'results/tifs_by_R/LCpatches_agb18_v3l1_Ap3WUw25u20_hti_sd.tif'
agb_filled_fp <- 'results/tifs_by_R/agb18_v3_l1_Ap3WUw25u20_hti_filled_LCpatches.tif'

# SF to SpatVector
lc_all_vect <- lc_sf_all %>% 
  terra::vect()

# Rasterize means
agb_by_lcpatch <- lc_all_vect %>% 
  terra::rasterize(agb_ras, field = 'mean', filename=agb_from_lc_fp, overwrite=T)

# Fill gaps and save
agb_filled <- terra::cover(agb_ras, agb_by_lcpatch, filename=agb_filled_fp, overwrite=T)

# Create uncertainty layer for agb_filled ----
# Rasterize SDs
agb_by_lcpatch_sd <- lc_all_vect %>% 
  terra::rasterize(agb_err, field = 'sd', filename=agb_from_lc_sd_fp, overwrite=T)

agb_filled_err <- terra::cover(agb_err, agb_by_lcpatch_sd, filename=agb_filled_fp, overwrite=T)


