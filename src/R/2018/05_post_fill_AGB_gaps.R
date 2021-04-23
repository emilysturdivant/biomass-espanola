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

results_dir <- 'data/results'

lc_fps <- c("data/raw/landcover/Lemoiner/Haiti2017_Clip.tif", 
            "data/raw/landcover/Lemoiner/DR_2017_clip.tif")
lc_fp_out <- "data/LULC/Hisp_2017_resALOS_terra.tif"
lc_out <- "data/LULC/Hisp_2017_resALOS_mskLand.tif"
agb_masked_fp <- file.path(results_dir, 'tifs_by_R', 
                           'agb18_v3_l1_mask_Ap3WUw25_u20_hti_qLee.tif')
msk_lnd_fp <- file.path(results_dir, 'masks', 
                        'hti18_maskLand.tif')
agb_from_lc_fp <- file.path(results_dir, 'tifs_by_R', 
                            'LCpatches_agb18_v3l1_Ap3WUw25u20_hti.tif')
agb_filled_fp <- file.path(results_dir, 'tifs_by_R', 
                           'agb18_v3_l3_Ap3WUw25u20_hti_filled_LCpatches.tif')

agb_capped_fp <- file.path(results_dir, 'tifs_by_R', 
                           'agb18_v3_l2_nomask_cap310.tif')

lc_pols_agb_fp <- str_glue("data/LULC/Haiti2017_Clip_polys_meanAGB.geojson")
lc_pols_agb_fp <- str_glue("data/LULC/Haiti2017_polys_AGBzonal.gpkg")

lc_stat <-  'median'
agb_filled_fp <- file.path(results_dir, 'tifs_by_R', str_glue('agb18_v3_l3_nomask_cap310_hti_LC{lc_stat}.tif'))
agb_from_lc_fp <- file.path(results_dir, 'tifs_by_R', str_glue('LCpatches_agb18_v3l2_nomask_cap310_{lc_stat}.tif'))
agb_from_lc_sd_fp <- file.path(results_dir, 'tifs_by_R', str_glue('LCpatches_agb18_v3l2_nomask_cap310_hti_{lc_stat}_sd.tif'))
agb_filled_sd_fp <- file.path(results_dir, 'tifs_by_R', str_glue('agb18_v3_l3_nomask_cap310_hti_LC{lc_stat}_sd.tif'))

mask_code <- 'mAWUw25u20'
cap_code <- 'cap310'
agb_masked_fp <- file.path(results_dir, 'tifs_by_R', 
                           str_glue('agb18_v3_l2_{mask_code}_{cap_code}.tif'))
# lc_pols_agb_fp <- file.path('data/LULC', 
#                             str_glue("Haiti2017_Clip_polys_{lc_stat}AGB_v3l2_{mask_code}_{cap_code}.gpkg"))
agb_filled_fp <- file.path(results_dir, 'tifs_by_R', 
                           str_glue('agb18_v3_l3_{mask_code}_{cap_code}_hti_LC{lc_stat}.tif'))
agb_from_lc_fp <- file.path(results_dir, 'tifs_by_R', 
                            str_glue('LCpatches_agb18_v3l2_{mask_code}_{cap_code}_{lc_stat}.tif'))
agb_from_lc_sd_fp <- file.path(results_dir, 'tifs_by_R', 
                               str_glue('LCpatches_agb18_v3l2_{mask_code}_{cap_code}_hti_{lc_stat}_sd.tif'))
agb_filled_sd_fp <- file.path(results_dir, 'tifs_by_R', 
                              str_glue('agb18_v3_l3_{mask_code}_{cap_code}_hti_LC{lc_stat}_sd.tif'))


# Load AGB ================================================================
agb_ras <- terra::rast(agb_capped_fp)

# ~ AGB mean by LC patch =========================================================
# Get mean AGB for each land cover patch (extract)
# Polygonize LULC raster ----
# output filenames
lc_multipols_fp <- "data/LULC/Haiti2017_Clip_multipolys"
lc_pols_fp <- "data/LULC/Haiti2017_Clip_polys.geojson"

if(!file.exists(lc_pols_fp)) {
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
}

# Extract AGB values ----
if(!file.exists(lc_pols_agb_fp)){
  # Load polygons
  lc_vect <- terra::vect(lc_pols_fp)
  
  # Don't use Water and Urban class
  lc_vect2 <- terra::subset(lc_vect, lc_vect$LC > 2) 
  
  agb_ex <- terra::extract(agb_ras, lc_vect2) 
  names(agb_ex) <- c('ID', 'agb')
  
  # Convert matrix to tibble and get mean and count for each polygon
  agb <- agb_ex %>% 
    as_tibble() %>% 
    group_by(ID) %>% 
    summarise(median = median(agb, na.rm=T),
              mean = mean(agb, na.rm=T),
              ct = sum(!is.na(agb)), 
              sd = sd(agb, na.rm = T))
  
  # Append columns to SF polys
  lc_sf <- sf::st_as_sf(as.data.frame(lc_vect2, geom=TRUE), 
                        wkt="geometry", 
                        crs=crs(lc_vect2))
  lc_sf_all <- agb %>% 
    select(-ID) %>% 
    cbind(lc_sf, .) %>% 
    filter(ct > 1)
  
  # Save
  lc_sf_all %>% st_write(lc_pols_agb_fp, delete_dsn = T)
  
} 

# Look
agb_raster <- raster::raster(agb_masked_fp)
lc_sf_all <- st_read(lc_pols_agb_fp)
tm_shape(lc_sf_all) + tm_fill(col='mean') +
  tm_shape(agb_raster) + tm_raster()

# Fill missing AGB with patch means --------------------------------------------
# agb_from_lc_sd_fp <- file.path(results_dir, 'tifs_by_R/LCpatches_agb18_v3l1_Ap3WUw25u20_hti_sd.tif')
# agb_filled_fp <- file.path(results_dir, 'tifs_by_R/agb18_v3_l3_Ap3WUw25u20_hti_filled_LCpatches.tif')
# agb_filled_sd_fp <- file.path(results_dir, 'tifs_by_R/agb18_v3_l3_Ap3WUw25u20_hti_filled_LCpatches_sd.tif')

# SF to SpatVector
lc_all_vect <- terra::vect(lc_pols_agb_fp)

if(!file.exists(agb_filled_fp)) {
  
  # Rasterize means
  agb_by_lcpatch <- lc_all_vect %>% 
    terra::rasterize(agb_ras, field = lc_stat, 
                     filename=agb_from_lc_fp, overwrite=T, 
                     wopt=list(datatype='FLT4S', gdal='COMPRESS=LZW'))
  
  # Fill gaps and save
  agb_ras <- terra::rast(agb_masked_fp)
  agb_filled <- terra::cover(agb_ras, agb_by_lcpatch, 
                             filename=agb_filled_fp, overwrite=T, 
                             wopt=list(datatype='FLT4S', gdal='COMPRESS=LZW'))
  
} else {
  
  agb_by_lcpatch <- terra::rast(agb_filled_fp)
  
}

# Create uncertainty layer for agb_filled ----
agb_err <- agb_ras <- terra::rast(agb_masked_fp)
agb_err[agb_err > 0] <- 23 # cross-validation RMSE

if(!file.exists(agb_filled_sd_fp)) {
  # Rasterize SDs
  agb_by_lcpatch_sd <- lc_all_vect %>% 
    terra::rasterize(agb_err, field = 'sd', filename=agb_from_lc_sd_fp, overwrite=T, 
                     wopt=list(datatype='FLT4S', gdal='COMPRESS=LZW'))
  
  agb_filled_err <- terra::cover(agb_err, agb_by_lcpatch_sd,
                                 filename=agb_filled_sd_fp, overwrite=T, 
                                 wopt=list(datatype='FLT4S', gdal='COMPRESS=LZW'))
  
}

agb_by_lcpatch_sd <- terra::rast(agb_from_lc_sd_fp)
plot(agb_by_lcpatch_sd)


# Pre-process land cover -------------------------------------------------------
if(!file.exists(lc_fp_out)){
  # Resample Haiti and DR land cover to AGB res 
  rs <- sapply(lc_fps, 
               function(fp) {
                 r <- terra::rast(fp)
                 r <- terra::resample(r, agb_ras, method='near')
                 return(r)}
  )
  
  # Merge, i.e. cover
  lc2 <- terra::cover(rs$`data/LULC/Haiti2017_Clip.tif`, 
                      rs$`data/LULC/DR_2017_clip.tif`, 
                      filename=lc_fp_out, overwrite=T, 
                      wopt=list(datatype='FLT4S', gdal='COMPRESS=LZW'))
}

# ~ AGB mean by LC (6 values) ----------------------------------------------------
agb_filled_fp <- file.path(results_dir, 'tifs_by_R/agb18_v3_l3_Ap3WUw25u20_hti_filled_6zones.tif')

# Create AGB surface from mean AGB for each LC class (6 values) 
# Load and crop LC
lc <- terra::rast(lc_fp_out) %>% 
  terra::crop(agb_ras) 

# Mask LC to land pixels
lc <- lc * terra::rast(msk_lnd_fp)

# Compute mean AGB for each LC (zonal statistics) 
zonal_stats <- terra::zonal(agb_ras, lc, 'mean', na.rm=T)

# Reclass LC to AGB values -> AGB by LC surface
agb_by_lc <- terra::classify(lc, zonal_stats)

# Fill gaps in AGB with AGB by LC ----------------------------------------------
agb_filled <- terra::cover(agb_ras, agb_by_lc, filename=agb_filled_fp, 
                           wopt=list(datatype='FLT4S', gdal='COMPRESS=LZW'))

if(file.exists(agb_filled_fp)) agb_filled <- terra::rast(agb_filled_fp)

plot(agb_filled)
plot(agb_ras)
