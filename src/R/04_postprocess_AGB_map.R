# ******************************************************************************
# Script to:
#     * Perform post-process on AGB map 
#       - count pixels meeting various criteria
# Proceeds:
#     * regression_AGB-g0.R - creates AGB map
#     * calculate_AGB.R - calculates AGB by plot from field data
#     * process_ALOS_tiles.R
# Requires:
#     * AGB map (agbr)
#     * pre-created masks
#
# ******************************************************************************
# Load libraries ----
library("terra")
library("tidyverse")

# Set variables ----
suffix <- g0_variant <- 'lee11s10'
suffix <- g0_variant <- 'med5'
p3_agb_val = 310
year <- '2019'

masks_dir <- file.path('data', 'tidy', str_c('palsar_', year), 'masks')
landmask_fp <- file.path(masks_dir, 'hti_land_palsar.tif')

suffix <- if(suffix == 'simple') '' else str_c('_', suffix)
agb_fp <- file.path('data', 'modeling', g0_variant, str_c("agb_l0", suffix, ".tif"))

# Apply value cap to AGB -------------------------------------------------------
saturation_pt <- 310
agb_cap_fp <- file.path(dirname(agb_fp), 
                        str_c('agb_l1_cap', saturation_pt, '.tif'))

# Apply saturation point to AGB and save
if(!file.exists(agb_cap_fp)){
  agb_dtype <- raster::dataType(raster::raster(agb_fp))
  agb_sat <- terra::rast(agb_fp) %>% 
    terra::classify(rbind(c(saturation_pt, Inf, saturation_pt)), 
                    filename = agb_cap_fp, 
                    overwrite=T, 
                    wopt = list(datatype = agb_dtype, gdal='COMPRESS=LZW'))
}

# Mask agb --------------------------------------------------------------
mask_agb <- function(agb_fp, p3_agb_val, masks_dir, level_code = 'l1', 
                     masks = c('A', 'WU', 'wb', 'u20', 'p3'), overwrite = TRUE,
                     masked_agb_fp) {
  
  # Get output filename
  masks_str <- masks %>% 
    str_replace('p3', str_c('t', p3_agb_val)) %>% 
    str_c(collapse = '')
  masked_agb_fp <- file.path(dirname(agb_fp), 
                             str_c('agb_', level_code, '_mask', masks_str, '.tif'))
  msk_all_fp <- file.path(dirname(agb_fp), str_c('mask', masks_str, '.tif'))
  
  # Stop if file already exists
  if(file.exists(masked_agb_fp) & !overwrite) {
    return(terra::rast(masked_agb_fp))
  }
  
  # Use mask file if it already exists
  if(file.exists(msk_all_fp) & !overwrite) {
    agb_ras <- terra::rast(agb_fp)
    msk <- terra::rast(msk_all_fp)
    
    # Apply non-negotiable mask (msk_AWUwb_u20) to AGB from SAGA filter
    agb_dtype <- raster::dataType(raster::raster(agb_fp))
    agb_masked <- agb_ras %>% 
      terra::mask(msk, filename = masked_agb_fp, overwrite = T, 
                  wopt = list(datatype = agb_dtype, gdal = 'COMPRESS=LZW'))
    
    # Return
    return(agb_masked)
  }
  
  # Load AGB
  agb_ras <- terra::rast(agb_fp)
  
  # Initialize mask raster (all 1's)
  msk <- agb_ras %>% terra::classify(rbind(c(-Inf, Inf, 1)))
  
  # AGB <20 mask
  if('u20' %in% masks){
    msk_u20 <- agb_ras %>% terra::classify(rbind(c(-Inf, 20, NA),
                                                 c(20, Inf, 1)))
    msk <- msk %>% terra::mask(msk_u20)
  }
  
  # g0 >0.3 mask (AGB>316.76)
  if('p3' %in% masks){
    msk_p3 <- agb_ras %>% terra::classify(rbind(c(p3_agb_val, Inf, NA), 
                                                c(-Inf, p3_agb_val, 1)))
    msk <- msk %>% terra::mask(msk_p3)
  }
  
  # LC17 Water and Urban, and OSM water with 25 m buffer
  if('WU' %in% masks & 'wb' %in% masks){
    msk_WUwb_fp <- file.path(masks_dir, "mask_WaterUrban_water25.tif")
    msk_WUwb <- terra::rast(msk_WUwb_fp) %>% terra::crop(agb_ras)
    
    # Check extents
    if(ext(agb_ras) != ext(msk_WUwb)){
      print("Extents don't match.")
    }
    
    msk <- msk %>% terra::mask(msk_WUwb)
    
  } else if('WU' %in% masks) {
      msk_fp <- file.path(masks_dir, "mask_WaterUrban.tif")
      msk_WU <- terra::rast(msk_fp) %>% terra::crop(agb_ras)
      
      # Check extents
      if(ext(agb_ras) != ext(msk_WU)){
        print("Extents don't match.")
      }
      
      msk <- msk %>% terra::mask(msk_WU)
      
  } else if('U' %in% masks) {
    msk_fp <- file.path(masks_dir, "mask_Urban.tif")
    msk_U <- terra::rast(msk_fp) %>% terra::crop(agb_ras)
    
    # Check extents
    if(ext(agb_ras) != ext(msk_U)){
      print("Extents don't match.")
    }
    
    msk <- msk %>% terra::mask(msk_U)
    
  }
  
  # LC17 Water and Urban, and OSM water with 25 m buffer
  if('A' %in% masks){
    msk_A_fp <- file.path(masks_dir, "mask_palsar_normal_2019.tif")
    msk_A <- terra::rast(msk_A_fp) %>% terra::crop(agb_ras)
    
    msk <- msk %>% terra::mask(msk_A)
  }
  
  # Save mask raster
  msk %>% terra::writeRaster(filename = msk_all_fp, 
                             overwrite = T, 
                             datatype = 'INT1S', 
                             gdal = 'COMPRESS=LZW')
  
  # Apply non-negotiable mask (msk_AWUwb_u20) to AGB from SAGA filter
  agb_dtype <- raster::dataType(raster::raster(agb_fp))
  agb_masked <- agb_ras %>% 
    terra::mask(msk, 
                filename = masked_agb_fp, 
                overwrite = T, 
                datatype = agb_dtype, 
                gdal = 'COMPRESS=LZW')
  
  # Return
  return(agb_masked)
}

# Apply masks ----
# Apply WU, wb, u20, p3 mask to L0 agb
agb_masked <- mask_agb(agb_fp,
                       p3_agb_val, 
                       masks_dir, 
                       level_code = 'l1',
                       masks = c('WU', 'wb', 'u20', 'p3'), overwrite = FALSE)

# Apply WU, wb, u20 mask to L1 capped AGB
agb_l2 <- mask_agb(agb_fp = agb_cap_fp, 
                       p3_agb_val, 
                       masks_dir, 
                       level_code = 'l2',
                       masks = c('WU', 'wb', 'u20'), overwrite = FALSE)

# Apply WU, wb mask to L1 capped AGB
agb_l2 <- mask_agb(agb_fp = agb_cap_fp, 
                   p3_agb_val, 
                   masks_dir, 
                   level_code = 'l2',
                   masks = c('WU', 'wb'), overwrite = FALSE)

# Apply WU, wb mask to L1 capped AGB
agb_l2 <- mask_agb(agb_fp = agb_cap_fp, 
                   p3_agb_val, 
                   masks_dir, 
                   level_code = 'l2',
                   masks = c('WU', 'u20'), overwrite = FALSE)

# Apply WU, wb mask to L1 capped AGB
agb_l2 <- mask_agb(agb_fp = agb_cap_fp, 
                       p3_agb_val, 
                       masks_dir, 
                       level_code = 'l2',
                       masks = c('U'), overwrite = FALSE)



# Report (incomplete update): Get value frequencies within masks ---------------
get_frequencies <- function(agb_ras, msk_name = 'none'){
  # Reclass masked raster
  rcl_tbl <- rbind(c(-Inf, 20, -1),
                   c(20, 132, 0),
                   c(132, 200, 1), 
                   c(200, 250, 2), 
                   c(250, 300, 3), 
                   c(300, 310, 4), 
                   c(310, Inf, 5)) %>% 
    as_tibble() %>% 
    rename(value = V3) %>% 
    mutate(name = str_c(V1, ' to ', V2))
  
  # Convert to matrix
  rcl_mat <- rcl_tbl %>% select(-name) %>% as.matrix()
  
  # Reclass masked raster
  agb_rcl <- agb_ras %>% 
    terra::classify(rcl_mat)
  
  # Get class frequencies
  cts <- agb_rcl %>% 
    freq(bylayer = FALSE) %>% 
    as_tibble() %>% 
    left_join(select(rcl_tbl, value, name)) %>% 
    mutate(mask = msk_name)
  
  # Return
  return(cts)
}

# Load
agb_ras <- terra::rast(agb_fp)
msk_land <- terra::rast(landmask_fp)
agb_masked <- agb_l2
mask_name <- 'WUu20'

# Tree cover 
mskinv_T_fp <- file.path(masks_dir, "mask_inv_TreeCover.tif")
mskinv_T <- terra::rast(mskinv_T_fp) %>% 
  terra::crop(agb_ras) %>% 
  terra::mask(msk_land)

# Totals
total_land <- global(!is.na(msk_land), 'sum', na.rm = TRUE) %>% deframe()
total_tc <- global(!is.na(mskinv_T), 'sum', na.rm = TRUE) %>% deframe()
total_mskd <- global(!is.na(agb_masked), 'sum', na.rm = TRUE) %>% deframe()

# Default mask
cts <- get_frequencies(agb_masked, mask_name)
cts <- cts %>% 
  mutate(pct = count / total_mskd,
         pct_land = count / total_land)

# Tree cover, no mask
agb_tc <- agb_ras %>% terra::mask(mskinv_T)
tc_cts <- get_frequencies(agb_tc, 'tree cover unmasked') %>% 
  mutate(pct = count / total_tc,
         pct_land = count / total_land)

# Tree cover, with mask
agb_tc_mskd <- agb_masked %>% terra::mask(mskinv_T)
tc_mskd_cts <- get_frequencies(agb_tc_mskd, str_c('tree cover, ', mask_name)) %>% 
  mutate(pct = count / total_tc,
         pct_land = count / total_land)

# Get mean and count
(tc_mn <- global(agb_tc_mskd, 'mean', na.rm=T))
ct_tc_pixels <- global(!is.na(agb_tc_mskd), 'sum')
(est_tc_area_ha <- (25*25*0.0001) * ct_tc_pixels)
est_tc_area_ha <- 953800
(est_total_AGB <- tc_mn * est_tc_area_ha)

# Get pct of TC class inside water buffer
msk_wb_fp <- file.path(masks_dir, 'mask_water_buff25m.tif')
msk_wb <- terra::rast(msk_wb_fp)
mskinv_T_wb <- mskinv_T %>% terra::mask(msk_wb, inverse = TRUE)
total_tc_wb <- global(!is.na(mskinv_T_wb), 'sum', na.rm = TRUE)
tc_wb_cts <- total_tc_wb %>% 
  rename(count = sum) %>% 
  mutate(name = 'water', 
         mask = 'tree cover', 
         pct = count / total_tc,
         pct_land = count / total_land)

combo_df <- bind_rows(cts, tc_mskd_cts, tc_wb_cts)

combo <- combo_df %>% 
  mutate(
    name = str_replace_all(name, 
                           "20 to 132|132 to 200|200 to 250|250 to 300|300 to 310|310 to Inf", 
                           'over20'))

combo %>% 
  dplyr::group_by(name, mask) %>% 
  dplyr::summarize(pct = sum(pct), 
                   pct_land = sum(pct_land), 
                   count = sum(count)) %>% 
  mutate(area_ha = count * (25*25*0.0001),
         area_from_pct = pct_land * 2710675)

hti_poly_fp <- "data/tidy/contextual_data/HTI_adm/HTI_adm0_fix.shp"
sf::st_area(st_read(hti_poly_fp))

# Get stats 
# Get mean and count
(mn <- global(agb_masked, 'mean', na.rm=T))
ct_valid_pixels <- global(!is.na(agb_masked), 'sum', na.rm=T)
(est_area_ha <- (25*25*0.0001) * ct_valid_pixels)
(est_total_AGB <- mn * est_area_ha)






# Get counts for v3l3 outside of function ----
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

# Filenames
dir <- file.path(results_dir, 'tifs_by_R')
agb_ras_fp <- file.path(dir, agb_fp) 
msk_land_fp <- file.path(results_dir, "masks/hti18_maskLand.tif")
wb_fp <- file.path(results_dir, 'masks/vector/osm_water_buff25m.shp')
mskinv_T_fp <- file.path(results_dir, "masks/mask_inv_TreeCover.tif")

# Crop to intersecting extents of AGB and land 
out <- crop_to_intersecting_extents(
  terra::rast(agb_ras_fp), 
  terra::rast(msk_land_fp), 
  return_r1=T, return_r2=T, return_bb=T)

# Get rasters
bbex <- out$bb
msk_land <- out$r2
agb.ras <- out$r1 %>% 
  terra::mask(msk_land)
mskinv_T <- terra::rast(mskinv_T_fp) %>% # 1== LC17 tree cover
  terra::crop(bbex)
mskinv_T_wb <- mskinv_T %>% 
  terra::mask(terra::vect(wb_fp))

plot(agb.ras)
plot(msk_land)

# Get mean and count
mn <- global(agb.ras, 'mean', na.rm=T)
sum_agb <- global(agb.ras, 'sum', na.rm=T)
ct_valid_pixels <- global(!is.na(agb.ras), 'sum')
est_area_ha <- (25*25*0.0001) * ct_valid_pixels
est_total_AGB <- mn * est_area_ha
ct_land_pixels <- global(!is.na(msk_land), 'sum')
est_area_ha_land <- (25*25*0.0001) * ct_land_pixels

# Convert to DF
df <- tibble(land = as.vector(msk_land), 
             agb = as.vector(agb.ras) %>% round(2), 
             zone = as.vector(mskinv_T),
             water_in_tc = as.vector(mskinv_T_wb)) %>% 
  filter(!is.na(land))



total <- nrow(df)

tc_df <- df %>% filter(zone == 1)
total_TC <- tc_df %>% nrow

# Tree cover class
(ct_TC_u20 <- tc_df %>% filter(agb < 20) %>% nrow)
ct_TC_o132 <- tc_df %>% filter(agb > 132) %>% nrow
ct_TC_o132 / total_TC
ct_TC_o300 <- tc_df %>% filter(agb > 300) %>% nrow
ct_TC_o300 / total_TC
ct_TC_o200 <- tc_df %>% filter(agb > 200) %>% nrow
ct_TC_o200 / total_TC
ct_TC_na <- tc_df %>% filter(is.na(agb)) %>% nrow
ct_TC_na / total_TC
ct_TC_wb <- tc_df %>% filter(water_in_tc == 1) %>% nrow
ct_TC_wb / total_TC


ct_na <- df %>% filter(is.na(agb)) %>% nrow
ct_na / total
ct_o132 <- df %>% filter(agb > 132) %>% nrow
ct_o132 / total
ct_o200 <- df %>% filter(agb > 200) %>% nrow
ct_o200 / total
ct_o250 <- df %>% filter(agb > 250) %>% nrow
ct_o250 / total
ct_o300 <- df %>% filter(agb > 300) %>% nrow
ct_o300 / total
ct_o310 <- df %>% filter(agb > 310) %>% nrow
ct_o310 / total
ct_u20 <- df %>% filter(agb < 20) %>% nrow
ct_u20 / total

ct_o20 <- df %>% filter(!is.na(agb), agb > 20) %>% nrow
ct_o20 / total

# Get pct of TC class that intersects water buffer

ct_TC_pixels_wb <- global(!is.na(mskinv_T_wb), 'sum')
TC_wb <- ct_TC_pixels_wb / ct_TC_pixels


# Visualize
ggplot(df, aes(x=agb)) +
  geom_histogram()

ggplot(df, aes(x=zone, y=agb)) +
  geom_bar(stat="identity")



