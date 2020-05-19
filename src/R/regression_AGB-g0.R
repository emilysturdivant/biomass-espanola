# *************************************************************************************************
# Script to:
#     * Calibrate regression model between AGB and backscatter 
# Proceeds:
#     * calculate_AGB.R - calculates AGB by plot from field data
#     * process_ALOS_tiles.R
# Requires:
#     * backscatter image (g0)
#     * plot polygons with AGB by plot (plots_agb)
#
# *************************************************************************************************

# Load libraries
# library(gridExtra)
# library(rgdal)
# library(caret)
library(raster)
library(tidyverse)
library(stars)
library(geobgu)
library(broom)
library(gdalUtils)
library(tmap)
require(graphics)

# Get mean backscatter for each plot ---- #############################################
# Load raster and polygon data
# Raster was created by 
# 1) biota to download PALSAR mosaic tiles and convert HV backscatter to natural units
# 2) process_ALOS_tiles.R to merge the tiles 
g0 <- read_stars("results/g0nu_HV/g0nu_2018_HV_haitiR.tif")

# Aggregate to 50m, as recommended by Saatchi 2015 and performed by Michelakis et al. 2015
# g0.nofilt <- raster("results/g0nu_HV/g0nu_2018_nofilt_HV_haiti.tif")
# g0.nofilt[g0.nofilt == 0] <- NA
# g0.agg <- aggregate(g0.nofilt, fact=2, fun=mean, na.rm=TRUE, 
#                     filename="results/g0nu_HV/g0nu_2018_haiti_agg50m.tif", 
#                     overwrite=TRUE)

# Add plot backscatter mean to polygons
plots_agb <- readRDS('results/R_out/plots_agb.rds')
plots_agb %>% mutate(
    g0l_mean = raster_extract(g0, plots_agb, fun = mean, na.rm = TRUE)
  ) %>% 
  saveRDS('results/R_out/plots_g0agb.rds')
g0_AGB <- readRDS('results/R_out/plots_g0agb.rds')
g0_AGB %>% 
  st_write("results/plots_values/plots_g0agb.shp", append=FALSE)

# Extract just AGB and backscatter as dataframe
g0.agb <- g0_AGB[c('AGB_ha', 'g0l_mean')] %>% 
  st_set_geometry(NULL) %>%
  as.data.frame() %>% 
  rename(AGB = AGB_ha, backscatter = g0l_mean) 
g0.agb %>% 
  saveRDS("results/R_out/plots_g0agb_dfslim.rds")
g0.agb  <- readRDS("results/R_out/plots_g0agb_dfslim.rds")

# Look at values
(backscatter <- c(
  all_mean=mean(g0.agb$backscatter),
  all_sd=sd(g0.agb$backscatter),
  bio_mean=mean(g0.agb[9:36, ]$backscatter),
  bio_sd=sd(g0.agb[9:36, ]$backscatter),
  min=min(g0.agb$backscatter),
  max=max(g0.agb$backscatter)))

# Scatterplot - AGB against backscatter ---- ##############################################
g0.agb$AGB <- as.vector(g0.agb$AGB)
p <- ggplot(g0.agb, aes(x=backscatter, y=AGB)) + geom_point() +
  labs(y = expression(paste("Aboveground biomass (Mg ha"^"-1", ")")), 
       x = expression(paste("Radar backscatter, ",sigma['HV']^0," (m"^2, "/m"^2, ")")))
(p <- p + geom_smooth(method="lm", se=TRUE, fullrange=TRUE, level=0.95, col='black'))

# Histograms of plot AGB and backscatter ---- ############################################
# AGB
p1 <-ggplot(g0.agb, aes(x=AGB)) + 
  geom_histogram(bins=30) +
  labs(x = expression(paste("Aboveground biomass (MgC ha"^"-1", ")")), 
       y = "Number of plots (N = 36)")+ 
  ylim(0, 9) 
#+  geom_vline(aes(xintercept=mean(AGB)), color="black", linetype="dashed", size=.5)
# Backscatter
p2 <-ggplot(g0.agb, aes(x=backscatter)) + 
  geom_histogram() +
  labs(x = expression(
    paste("Radar backscatter, ", sigma['HV']^0, " (m"^2, "/m"^2, ")")), 
    y = "Number of plots (N = 36)")+ 
  ylim(0, 9)
#+ geom_vline(aes(xintercept=mean(backscatter)), color="black", linetype="dashed", size=.5)
grid.arrange(p1, p2, nrow = 2)

# Linear regression ---- ###########################################################
# See how correlation improves without plots 7 and 8, which are probably responding to soil moisture. 
# g0.agb <- g0.agb %>% filter(!(AGB==0 & backscatter>0.016))

# Basic OLS regression
# g0.agb$backscatter <- as.vector(g0.agb$backscatter)
ols <- lm(as.vector(AGB) ~ as.vector(backscatter), data=g0.agb, x=TRUE, y=TRUE)
ols %>% saveRDS("results/R_out/ols_AGBv1_g0v1.rds")

# Report results ---- ###########################################################
report_ols_results <- function(ols){
  # Regression Coefficients
  coefs <- tidy(ols) %>% cbind(confint(ols)) %>% 
    mutate(term = str_replace(term, '\\(Intercept\\)', 'int')) %>% 
    mutate(term = str_replace(term, 'backscatter', 'g0')) %>% 
    mutate(term = str_replace(term, 'as\\.vector\\(', '')) %>% 
    mutate(term = str_replace(term, '\\)', '')) %>% 
    pivot_wider(names_from=term, values_from=estimate:'97.5 %') %>% 
    cbind(glance(ols))
  # Accuracy metrics
  mse <- mean((residuals(ols))^2)
  rss <- sum(residuals(ols)^2)
  acc_metrics <- list(
    rmse = sqrt(mse), # RMSE
    mae = mean(abs(residuals(ols))),
    mse = mse,
    rss = rss,
    mss = sum(residuals(ols)^2)/ols$df.residual,
    rse = sqrt(rss / ols$df.residual),
    cov = cov2cor(vcov(ols))[1,2]
  ) %>% as.data.frame()
  # ANOVA
  ANOVA <- anova(ols) %>% tidy() %>% 
    mutate(term = str_replace(term, 'backscatter', 'g0anova')) %>% 
    mutate(term = str_replace(term, 'as\\.vector\\(', '')) %>% 
    mutate(term = str_replace(term, '\\)', '')) %>% 
    mutate(term = str_replace(term, 'Residuals', 'resid')) %>% 
    pivot_wider(names_from=term, values_from=df:p.value)
  # Correlation coefficients
  sp <- cor.test(x=g0.agb$backscatter, y=g0.agb$AGB, method = 'spearman') %>% 
    tidy() %>% 
    select(-c(alternative, method)) %>% 
    mutate(term='rhoSpear')
  pe <- cor.test(x=g0.agb$backscatter, y=g0.agb$AGB, method = 'pearson') %>% 
    tidy() %>% 
    select(-c(alternative, method)) %>% 
    mutate(term='Pearson')
  # Join correlation coefficients
  cors <- rbind(pivot_longer(sp, cols=estimate:p.value),
                pivot_longer(pe, cols=estimate:conf.high)
    ) %>%
    mutate(stat = str_c(name, term, sep="_")) %>% 
    select(-c(term, name)) %>% 
    column_to_rownames('stat') %>% 
    t() %>% as.data.frame()
  # Combine
  top <- cbind(
    coefs %>% select(starts_with('estimate'), adj.r.squared),
    ci_int=(confint(ols)[1,2]-confint(ols)[1,1])/2,
    ci_slope=(confint(ols)[2,2]-confint(ols)[2,1])/2,
    acc_metrics %>% as.data.frame() %>% select(rmse, mae, mse, rss, mss, rse),
    coefs %>% select(!starts_with('estimate'), -adj.r.squared),
    covariance=cov2cor(vcov(ols))[1,2],
    ANOVA,
    cors
  ) %>% t() %>% as.data.frame()
}
vals <- report_ols_results(ols)

# Plot error plots
opar <- par(mfrow = c(2,2), oma = c(0, 0, 1.1, 0))
plot(ols, las = 1)

# Repeated k-fold cross validation ---- ###########################################################
fxn.bias <- function(data, lev = NULL, model = NULL) {
  resids <- data$pred - data$obs
  rss <- sum(resids^2)
  n <- length(resids)
  df <- n-2
  mse <- rss / n
  c(RMSE=sqrt(mse),
    Rsquared=summary(lm(pred ~ obs, data))$r.squared,
    MAE=sum(abs(resids)) / n,
    MSE=mse,
    B=abs(sum(resids) / n),
    RSS=rss,
    MSS=rss/df,
    RSE=sqrt(rss / df))
}
set.seed(45)
model.10000x10 <- train(AGB ~ backscatter, data = g0.agb, method = "lm",
                        trControl = trainControl(method = "repeatedcv", 
                                                 number = 10, repeats = 10000,
                                                 summaryFunction = fxn.bias))
cv_results <- t(model.10000x10$results)
View(model.10000x10$results)
set.seed(45)
model.10000x5 <- train(AGB ~ backscatter, data = g0.agb, method = "lm",
                       trControl = trainControl(method = "repeatedcv", 
                                                number = 5, repeats = 10000,
                                                summaryFunction = fxn.bias))
model.10000x5 %>% saveRDS("results/R_out/CVmodel_g0nu_10000x5.rds")

cv_r <- model.10000x5$results
cv_r <- as_tibble(cbind(metric = names(cv_r), t(cv_r))) %>% 
  rename(value='1')
cv_r1 <- cv_r %>% 
  filter(!str_detect(metric, 'SD$'), !str_detect(metric, 'intercept'))
cv_r2 <- cv_r %>% 
  filter(str_detect(metric, 'SD$')) %>% 
  select(value) %>% 
  rename(SD=value)
cv_r <- bind_cols(cv_r1, cv_r2)
View(cv_r)
View(model.10000x5$results)
model.10000x5$results[-1]
mets <- cbind(vals[1], model.10000x5$results[-1]) %>% t()

# Create AGB raster ---- ###########################################################
# Load data 
g0 <- raster("results/g0nu_HV/g0nu_2018_HV.tif"); names(g0) <- 'backscatter'
ols <- readRDS("results/R_out/ols_AGBv1_g0v1.rds")

# Apply linear regression model to create AGB map
agb.ras <- raster::predict(g0, ols, na.rm=TRUE)
agb.ras %>% writeRaster("results/tifs_by_R/agb18_v1_l0.tif")

# Make masks based on AGB ----####################################################################
# AGB <20 mask
agb.ras <- read_stars("results/tifs_by_R/agb18_v1_l0.tif")
msk_u20 <- agb.ras
msk_u20[msk_u20 < 20] <- NA
msk_u20[msk_u20 >= 20] <- 1
msk_u20 %>% saveRDS("results/R_out/mask_AGB_under20_stars.rds")

# Inverse
mskinv_u20 <- agb.ras
mskinv_u20[mskinv_u20 < 20] <- 1
mskinv_u20[mskinv_u20 >= 20] <- 0
mskinv_u20 %>% saveRDS("results/R_out/mask_inv_AGB_under20_stars.rds")

# Load masks ----########################################################
# AGB-based
msk_u20 <- # NA==AGB<20 and ALOS mask
  readRDS("results/R_out/mask_AGB_under20_stars.rds")
mskinv_u20 <- # 1== where AGB<20 in valid ALOS land
  readRDS("results/R_out/mask_inv_AGB_under20_stars.rds")
mskinv_u20_noWU <- # 1== AGB<20 with ALOS and WaterUrban masks applied
  readRDS("results/R_out/mask_inv_AGBunder20_WU_stars.rds")
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

# Combined masks
msk_AWU <- # 1==Valid ALOS land without WaterUrban
  readRDS("results/R_out/mask_ALOS_WaterUrban_stars.rds")
msk_Ww <- # NA==ALOS mask, LC17 water, and OSM water with 25 m buffer
  readRDS("results/R_out/mask_ALOS_allwater_raster.rds")
msk_Ww %>% st_as_stars() %>% saveRDS("results/R_out/mask_ALOS_allwater_stars.rds")

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
mskinv_WUWPb <- # >0 == where WaterUrban and OSM water 25m overlap valid ALOS land
  readRDS("results/R_out/mask_inv_WUWPb_raster.rds")

# Report starting number of NAs
df_mskA <- readRDS('results/R_out/mask_pcts.rds')
df_LCcounts <- readRDS('results/R_out/lc_ALOS_pixel_counts_tbl.rds')

# Get counts from mostly inverse masks
(lnd <- sum(msk_land[[1]]==1, na.rm=TRUE))
(wu <- sum(mskinv_WU[[1]]==1, na.rm=TRUE))
(wp <- sum(getValues(mskinv_WP)==1, na.rm=TRUE))
(wpb <- sum(getValues(mskinv_WPb)==1, na.rm=TRUE)) 
(wuwpb <- sum(mskinv_WUWPb[[1]]>0, na.rm=TRUE))

(tr <- sum(mskinv_T[[1]]==1, na.rm=TRUE))



# Counts
(u20 <- sum(mskinv_u20[[1]]==1, na.rm=TRUE))
(u20_noWU <- sum(mskinv_u20_noWU[[1]]==1, na.rm=TRUE))
msk_u20_T <- msk_u20 * mskinv_T
(u20_T <- sum(msk_u20_T[[1]]==1, na.rm=TRUE))

# Mask out WU
msk_u20_noWU <- msk_u20*msk_WU
msk_u20_noWU %>% saveRDS("results/R_out/mask_inv_AGBunder20_WU_stars.rds")

# Combine all masks
msk_Ap3_WU_wb <- msk_WU * msk_p3 * msk_Aw 
msk_Ap3_WU_wb %>% saveRDS('results/R_out/mask_ALOS_pt3_WaterUrban_water25_stars.rds')
msk_Ap3_WU_wb_u20 <- msk_Ap3_WU_wb * msk_u20
msk_Ap3_WU_wb_u20 %>% saveRDS('results/R_out/mask_ALOS_pt3_WaterUrban_water25_AGBu20_stars.rds')



# PLOT using tmap (interactive) ----############################################
test_ext <- extent(-72.7, -72.5, 18.2, 18.35)
tmap_mode("view")

# Crop for testing
bb <- st_bbox(c(xmin=-72.68, xmax=-72.51, ymin=18.22, ymax=18.32)) %>% 
  st_as_sfc()
st_crs(bb) <- 4326

# View some masks
tm_shape(msk_Ap3_WU_wb[bb]) + tm_raster() +
  tm_shape(msk_Ap3_WU_wb_u20[bb]) + tm_raster() +
  tm_shape(water_polysb) + tm_borders()

# View some masks
tm_shape(crop(msk_Ap3_WU_wb, test_ext)) + tm_raster() +
  tm_shape(crop(msk_Ap3_WU_wb_u20, test_ext)) + tm_raster() +
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
