# ******************************************************************************
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
# ******************************************************************************

# Load libraries ----
library('sf')
library('terra')
library('tidyverse')
library('tmap')
tmap_mode('view')

# Set variables ----
g0_variant <- 'lee11s10'
rmse_cv <- 21.2
year <- '2019'
lc_stat <-  'median'
input_level <- 'l2_mask'
mask_level <- 'LU'
code <- 'HV_nu'

# Input filepaths
agb_dir <- file.path('data', 'modeling', code, g0_variant)
(agb_capped_fp <- list.files(agb_dir, str_c('agb_', input_level, '.*\\.tif'), full.names = TRUE))
(agb_masked_fp <- list.files(agb_dir, str_glue('agb_l2.*{mask_level}\\.tif'), 
                             full.names = TRUE))

lc_fps <- c("data/raw/landcover/Lemoiner/Haiti2017_Clip.tif", 
            "data/raw/landcover/Lemoiner/DR_2017_clip.tif")

# Output filepaths
lc_pols_fp <- "data/tidy/landcover/Haiti2017_Clip_polys.gpkg"

agb_by_lc_prefix <- file.path(agb_dir, 'agb_by_landcover', 
                              str_glue('agb_{input_level}_{mask_level}_{lc_stat}_byLC'))
lc_pols_agb_fp <- str_c(agb_by_lc_prefix, '.gpkg')
agb_by_lc_fp <- str_c(agb_by_lc_prefix, '.tif') 
agb_by_lc_sd_fp <- str_c(agb_by_lc_prefix, '_sd.tif') 

agb_filled_fp <- file.path(agb_dir, str_glue('agb_l3_fillLC{lc_stat}_{input_level}_{mask_level}.tif'))
agb_filled_sd_fp <- file.path(agb_dir, str_glue('agb_l3_fillLC{lc_stat}_{input_level}_{mask_level}_sd.tif'))

# Load AGB ================================================================
agb_ras <- terra::rast(agb_capped_fp)
agb_dtype <- raster::dataType(raster::raster(agb_capped_fp))

# ~ AGB mean by LC patch =========================================================
# Get mean AGB for each land cover patch (extract)

# Polygonize LULC raster ----
if(!file.exists(lc_pols_fp)) {
  # Temporary file
  lc_multipols_fp <- "data/tidy/landcover/Haiti2017_Clip_multipolys"
  
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
summarize_raster_by_polygons <- function(lc_pols_fp, agb_ras, lc_pols_agb_fp) {
  
  if(file.exists(lc_pols_agb_fp)) return(lc_pols_agb_fp)
  

  # Load polygons
  lc_vect <- terra::vect(lc_pols_fp)
  
  # Don't use Water and Urban class
  lc_vect2 <- terra::subset(lc_vect, lc_vect$LC > 2) 
  
  # Extract AGB to polygons
  agb_ex <- terra::extract(agb_ras, lc_vect2) 
  names(agb_ex) <- c('ID', 'agb')
  
  # Convert matrix to tibble and get mean and count for each polygon
  agb <- agb_ex %>% 
    as_tibble() %>% 
    dplyr::group_by(ID) %>% 
    dplyr::summarise(median = median(agb, na.rm=T),
              mean = mean(agb, na.rm=T),
              ct = sum(!is.na(agb)), 
              sd = sd(agb, na.rm = T))
  
  # Append columns to SF polys
  x <- terra::as.data.frame(lc_vect2, geom = 'WKT')
  lc_sf <- sf::st_as_sf(x, 
                        wkt="geometry", 
                        crs=crs(lc_vect2))
  
  # Append columns to SF polys
  lc_sf_all <- agb %>% 
    dplyr::select(-ID) %>% 
    cbind(lc_sf, .) %>% 
    dplyr::filter(ct > 1)
  
  # Save
  dir.create(dirname(lc_pols_agb_fp), showWarnings = FALSE, recursive = TRUE)
  lc_sf_all %>% st_write(lc_pols_agb_fp, append = FALSE)
  
  # Return
  return(lc_pols_agb_fp)
}

summarize_raster_by_polygons(lc_pols_fp, agb_ras, lc_pols_agb_fp) 

# Look
# agb_raster <- raster::raster(agb_masked_fp)
# lc_sf_all <- st_read(lc_pols_agb_fp)
# tm_shape(lc_sf_all) + tm_fill(col='mean') +
#   tm_shape(agb_raster) + tm_raster()

# Fill missing AGB with patch means --------------------------------------------
# agb_by_lc_sd_fp <- file.path(results_dir, 'tifs_by_R/LCpatches_agb18_v3l1_Ap3WUw25u20_hti_sd.tif')
# agb_filled_fp <- file.path(results_dir, 'tifs_by_R/agb18_v3_l3_Ap3WUw25u20_hti_filled_LCpatches.tif')
# agb_filled_sd_fp <- file.path(results_dir, 'tifs_by_R/agb18_v3_l3_Ap3WUw25u20_hti_filled_LCpatches_sd.tif')

# SF to SpatVector
lc_all_vect <- terra::vect(lc_pols_agb_fp)

if(!file.exists(agb_filled_fp)) {
  
  # Rasterize means
  agb_by_lcpatch <- lc_all_vect %>% 
    terra::rasterize(agb_ras, field = lc_stat)

  # Fill gaps and save
  agb_ras <- terra::rast(agb_masked_fp)
  agb_filled <- terra::cover(agb_ras, 
                             agb_by_lcpatch, 
                             filename = agb_filled_fp, 
                             overwrite = TRUE, 
                             wopt = list(datatype = agb_dtype, gdal='COMPRESS=LZW'))
  
}

# Create uncertainty layer for agb_filled ----
if(!file.exists(agb_filled_sd_fp)) {
  
  # RMSE from CV for all areas with final AGB
  agb_err <- terra::rast(agb_masked_fp) %>% 
    terra::classify(rbind(c(0, Inf, rmse_cv)))
  
  # Rasterize SDs
  agb_by_lcpatch_sd <- lc_all_vect %>% terra::rasterize(agb_err, field = 'sd')
  
  # Use SDs from landcover to fill gaps
  agb_filled_err <- terra::cover(agb_err, 
                                 agb_by_lcpatch_sd,
                                 filename = agb_filled_sd_fp, 
                                 overwrite = TRUE, 
                                 wopt = list(datatype= agb_dtype, gdal='COMPRESS=LZW'))
  
}

agb_by_lcpatch_sd <- terra::rast(agb_by_lc_sd_fp)
plot(agb_by_lcpatch_sd)


# Pre-process land cover -------------------------------------------------------
lc_fp_out <- "data/tidy/landcover/Hisp_2017_resALOS_terra.tif"

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
                      wopt=list(datatype='INT1U', gdal='COMPRESS=LZW'))
}
