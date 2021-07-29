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
source('src/R/initialize_vars.R')
# year <- '2019'
# code <- 'HV_nu'
# saturation_pt <- 300
# modeling_dir <- file.path('data/modeling', code)
# 
# masks_dir <- file.path('data', 'tidy', str_c('palsar_', year), 'masks')
# landmask_fp <- file.path(masks_dir, 'hti_land_palsar.tif')

# Get level 0 AGB
# (dirs <- list.dirs(modeling_dir, recursive = FALSE))
# (mod_dir <- dirs[[1]])
# (agb_fp <- list.files(mod_dir, 'agb_l0.*\\.tif$', full.names = TRUE))

# # Get variant code
# (items <- str_split(dirname(agb_l0_fp), '/') )
# g0_variant <- nth(items[[1]], -1)

# Apply value cap to AGB -------------------------------------------------------
agb_cap_fp <- file.path(dirname(agb_l0_fp), 
                        str_c('agb_l1_cap', saturation_pt, '.tif'))

# Apply saturation point to AGB and save
if(!file.exists(agb_cap_fp)){
  agb_dtype <- raster::dataType(raster::raster(agb_l0_fp))
  agb_sat <- terra::rast(agb_l0_fp) %>% 
    terra::classify(rbind(c(saturation_pt, Inf, saturation_pt),
                          c(-Inf, 0, 0)), 
                    filename = agb_cap_fp, 
                    overwrite=T, 
                    wopt = list(datatype = agb_dtype, gdal='COMPRESS=LZW'))
}

# Mask agb --------------------------------------------------------------
mask_agb <- function(agb_fp, masks_dir, level_code = 'l1', 
                     masks = c('L', 'WU', 'wb', 'u20'), overwrite = TRUE,
                     masked_agb_fp) {
  
  # Get output filename
  masks_str <- masks %>% 
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
  
  # ALOS mask
  if('A' %in% masks){
    msk_fp <- file.path(masks_dir, "mask_palsar_normal_2019.tif")
    msk_tmp <- terra::rast(msk_fp) %>% terra::crop(agb_ras)
    
    msk <- msk %>% terra::mask(msk_tmp)
    
  } else if('L' %in% masks) {
    msk_fp <- file.path(masks_dir, "mask_palsar_layover_2019.tif")
    msk_tmp <- terra::rast(msk_fp) %>% terra::crop(agb_ras)
    
    msk <- msk %>% terra::mask(msk_tmp)
  }
  
  # Save mask raster
  msk %>% terra::writeRaster(filename = msk_all_fp, 
                             overwrite = TRUE, 
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
# Apply layover and Urban mask to L1 capped AGB
agb_l2 <- mask_agb(agb_fp = agb_cap_fp, 
                   masks_dir, 
                   level_code = 'l2',
                   masks = c('L', 'U'), overwrite = FALSE)

# Apply WU, wb, u20, p3 mask to L0 agb
agb_masked <- mask_agb(agb_l0_fp,
                       masks_dir, 
                       level_code = 'l1',
                       masks = c('WU', 'wb', 'u20'), overwrite = FALSE)

# Apply WU, wb, u20 mask to L1 capped AGB
agb_l2 <- mask_agb(agb_fp = agb_cap_fp, 
                       masks_dir, 
                       level_code = 'l2',
                       masks = c('WU', 'wb', 'u20'), overwrite = FALSE)

# Apply WU, wb mask to L1 capped AGB
agb_l2 <- mask_agb(agb_fp = agb_cap_fp, 
                   masks_dir, 
                   level_code = 'l2',
                   masks = c('WU', 'wb'), overwrite = FALSE)

# Apply WU, u20 mask to L1 capped AGB
agb_l2 <- mask_agb(agb_fp = agb_cap_fp, 
                   masks_dir, 
                   level_code = 'l2',
                   masks = c('WU', 'u20'), overwrite = FALSE)

# Apply Urban mask to L1 capped AGB
agb_l2 <- mask_agb(agb_fp = agb_cap_fp, 
                       masks_dir, 
                       level_code = 'l2',
                       masks = c('U'), overwrite = FALSE)

# Apply L, WU, wb mask to L1 capped AGB
agb_l2 <- mask_agb(agb_fp = agb_cap_fp, 
                   masks_dir, 
                   level_code = 'l2',
                   masks = c('L', 'WU', 'wb'), overwrite = FALSE)

# Apply L, WU, wb mask to L1 capped AGB
agb_l2 <- mask_with_options(in_fp = agb_cap_fp, 
                            masks_dir, 
                            code = 'agb_l2',
                            masks = c('L', 'WU', 'wb'), overwrite = FALSE, 
                            dirname(agb_cap_fp))

# Apply WU, wb mask to L1 capped AGB
agb_l2 <- mask_with_options(in_fp = agb_cap_fp, 
                            masks_dir, 
                            code = 'agb_l2',
                            masks = c('WU', 'wb'), overwrite = FALSE, 
                            dirname(agb_cap_fp))

# Report (incomplete update): Get value frequencies within masks ---------------
get_frequencies <- function(agb_ras, msk_name){

  # Create reclass matrix
  rcl_mat <- rbind(c(-Inf, 20, -1),
                   c(20, 132, 0),
                   c(132, 200, 1), 
                   c(200, 250, 2), 
                   c(250, 300, 3), 
                   c(300, 310, 4), 
                   c(310, Inf, 5)) 
  colnames(rcl_mat) <- c('V1', 'V2', 'value')
  
  # Convert to tibble and create strings
  rcl_tbl <- rcl_mat %>% 
    as_tibble() %>% 
    mutate(name = str_c(V1, ' to ', V2))
  
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
agb_ras <- terra::rast(agb_l0_fp)
msk_land <- terra::rast(landmask_fp)
agb_masked <- agb_l2
# mask_name <- 'WUu20'

# Get name of mask
mask_name <- agb_ras@ptr$filenames %>% basename() %>% 
  str_extract("(?<=mask).*(?=\\.tif)")
mask_name <- agb_masked@ptr$filenames %>% basename() %>% 
  str_extract("(?<=mask).*(?=\\.tif)")
mask_name <- ifelse(mask_name == '', 'none', mask_name)

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

combo_df %>% write_csv(file.path(dirname(agb_cap_fp), '04_pcts_masked.csv'))
combo_df %>% write_csv(file.path('data/reports', '04_pcts_masked.csv'))

combo <- combo_df %>% 
  mutate(
    name = str_replace_all(name, 
                           "20 to 132|132 to 200|200 to 250|250 to 300|300 to 310|310 to Inf", 
                           'over20')) %>% 
  dplyr::group_by(name, mask) %>% 
  dplyr::summarize(pct = sum(pct), 
                   pct_land = sum(pct_land), 
                   count = sum(count)) %>% 
  mutate(area_ha = count * (25*25*0.0001),
         area_from_pct = pct_land * 2710675)

combo %>% write_csv(file.path(dirname(agb_cap_fp), '04_pcts_masked_agg_w_areas.csv'))
combo %>% write_csv(file.path('data/reports', '04_pcts_masked_agg_w_areas.csv'))




hti_poly_fp <- "data/tidy/contextual_data/HTI_adm/HTI_adm0_fix.shp"
sf::st_area(st_read(hti_poly_fp))

# Get stats 
# Get mean and count
(mn <- global(agb_masked, 'mean', na.rm=T))
ct_valid_pixels <- global(!is.na(agb_masked), 'sum', na.rm=T)
(est_area_ha <- (25*25*0.0001) * ct_valid_pixels)
(est_total_AGB <- mn * est_area_ha)



# Histogram ----
# Function from 07...Rmd
plot_agb_hist_density <- function(x, agb_pal, bwidth=5, sample_size=100000, 
                                  xlim = c(min=0, max=300)) {
  
  # Get values
  lyr_name <- x$name
  ras <- terra::rast(x$fp)
  df <- terra::values(ras, dataframe = TRUE) %>% 
    rename(value = 1) %>% 
    drop_na()
  
  # Get statistics
  mu <- mean(df$value, na.rm=T)
  sd1 <- sd(df$value, na.rm=T)
  mn <- min(df$value, na.rm=T)
  mx <- max(df$value, na.rm=T)
  
  cuts <- quantile(df$value, probs = c(.1, .5, .9)) %>% as_tibble(rownames='ref')
  cuts_pts <- data.frame(x = c(mn, mx), 
                         y = c(0, 0))
  
  # Sample
  if (nrow(df) > sample_size) {
    print("Sampling...")
    df <- df %>% sample_n(sample_size)
  }
  
  # Plot
  p <- ggplot(df, aes(x=value)) + 
    geom_histogram(aes(y=..density.., fill = ..x..), binwidth = bwidth,
                   show.legend = FALSE)+
    scale_fill_gradientn(colors = agb_pal$colors, 
                         # breaks = c(-60, -30, 0, 30, 60),
                         # labels = c('-60 (internal < external)', '', '0', '', '60 (internal > external)'),
                         name = expression(paste("AGB (Mg ha"^"-1", ")")),
                         # na.value='lightgray', 
                         limits = c(agb_pal$min, agb_pal$max),
                         oob = scales::squish,
                         guide = guide_colorbar(barheight = 2)
    ) + 
    geom_density(alpha=.2,
                 show.legend = FALSE) + 
    scale_x_continuous(name = expression(paste("Aboveground biomass (Mg ha"^"-1", ")")),
                       # breaks = seq(min, max, 100),
                       # limits = c(min, max)
    ) +
    coord_cartesian(xlim = xlim) +
    geom_vline(mapping = aes(xintercept = value,
                             colour = ref),
               data = cuts,
               color="black", 
               # linetype="solid", 
               linetype="dashed", 
               size=.5,
               show.legend = FALSE) +
    geom_text(mapping = aes(x = value,
                            y = Inf,
                            label = ref,
                            hjust = -.1,
                            vjust = 1),
              data = cuts) +
    ggtitle(lyr_name) +
    theme_minimal() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(), 
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank())
}

p_hist <- plot_agb_hist_density(agb_fps$internal, agb_pal)
hist_fp <- file.path(agb_dir, 'reports', str_c('agb_', agb_code, '_hist.png'))
dir.create(dirname(hist_fp))
ggsave(hist_fp, p_hist, width = 5, height = 3, dpi = 120)