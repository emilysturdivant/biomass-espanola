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

# *** VARIABLE ***
g0_fname <- "results/g0nu_HV/g0nu_2018_HV_haitiR.tif"

# Get mean backscatter for each plot ---- #############################################
# Load raster and polygon data
# Raster was created by 
# 1) biota to download PALSAR mosaic tiles and convert HV backscatter to natural units
# 2) QGIS to merge the tiles 
# 3) biota to run Radar Enhanced Lee Filter 
# 3.b) Possibly mask out water and urban features and extreme values. Convert to No Data using Raster Calc (0/0)
# 4) Export Filtered Grid (No Data == -99999)
g0 <- read_stars(g0_fname)
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

# Masks ---- ####################################################################
# Report starting number of NAs
df_mskA <- readRDS('results/R_out/mask_pcts.rds')
df_LCcounts <- readRDS('results/R_out/lc_ALOS_pixel_counts_tbl.rds')

# ALOS Normal mask
msk_A <- readRDS("results/R_out/mask_ALOS_raster.rds")

# Load LC 
lc <- readRDS("results/R_out/LC17_masked_to_ALOS_land_raster.rds")

# Mask out WaterUrban from land
msk_WU <- lc; rm(lc)
msk_WU[msk_WU<3] <- NA
msk_WU[!is.na(msk_WU)] <- 1
msk_WU %>% saveRDS("results/R_out/mask_WaterUrban_raster.rds")
msk_WU <- readRDS("results/R_out/mask_WaterUrban_raster.rds")

# Combine ALOS and WaterUrban masks
msk_AWU <- msk_WU*msk_A
msk_AWU %>% saveRDS("results/R_out/mask_ALOS_WaterUrban_raster.rds")
msk_AWU <- readRDS("results/R_out/mask_ALOS_WaterUrban_raster.rds")

plot(msk_AWU[1, 11000:14000, 7000:9000])

# Create masks based on backscatter values
g0 <- read_stars("results/g0nu_HV/g0nu_2018_HV.tif")

# Create mask of >0.3
create_inverse_mask_at_threshold <- function(r, threshold, gt=TRUE, filename=NULL){
  if (gt) {
    r[r <= threshold] <- 0
    r[r > threshold] <- 1
  } else {
    r[r < threshold] <- 1
    r[r >= threshold] <- 0
  }
  if (filename) r %>% saveRDS(filename)
}
msk_p3 <- create_inverse_mask_at_threshold(
  g0, 0.3, "results/R_out/mask_ALOS_g0overp3inverse_raster.rds"
  )
msk_p3 <- g0
msk_p3[msk_p3<=0.3] <- 0
msk_p3[msk_p3>0.3] <- 1
msk_p3 %>% saveRDS("results/R_out/mask_ALOS_g0overp3inverse_raster.rds")
msk_p3 <- readRDS("results/R_out/mask_ALOS_g0overp3inverse_raster.rds")
(p3 <- sum(msk_p3[[1]]==1, na.rm=TRUE))

# Create mask of >0.3
msk_p2 <- g0; rm(g0)
msk_p2[msk_p2<=0.2] <- 0
msk_p2[msk_p2>0.2] <- 1
msk_p2 %>% saveRDS("results/R_out/mask_ALOS_g0overp2inverse_raster.rds")
msk_p2 <- readRDS("results/R_out/mask_ALOS_g0overp2inverse_raster.rds")
(p2 <- sum(msk_p2[[1]]==1, na.rm=TRUE))

# Create inverse mask for each new mask and then multiple by msk_A

# Get presence of water/urban and then mask with ALOS
msk_WU <- readRDS("results/R_out/LC17_masked_to_ALOS_land_raster.rds")
msk_WU[msk_WU==2] <- 1 # Water is already 1 and now urban is as well
msk_WU[msk_WU!=1] <- NA
(wu1 <- sum(msk_WU[[1]]==1, na.rm=TRUE))
mskinv_WU <- msk_WU*msk_A
(wu0 <- sum(mskinv_WU[[1]]==1, na.rm=TRUE))

plot(mskinv_WU[1, 11000:14000, 7000:9000])




# Crop g0 to match 
agb.ras <- read_stars("results/tifs_by_R/agb18_v1_l0.tif")
msk_u20 <- agb.ras
msk_u20[msk_u20 < 20] <- 1
msk_u20[msk_u20 >= 20] <- 0
msk_u20 %>% saveRDS("results/R_out/mask_AGB_under20_inverse_raster.rds")
(u20 <- sum(msk_u20[[1]]==1, na.rm=TRUE))

# Mask out WU
msk_WU <- readRDS("results/R_out/mask_WaterUrban_raster.rds")
msk_u20_noWU <- msk_u20*msk_WU
(u20_noWU <- sum(msk_u20_noWU[[1]]==1, na.rm=TRUE))

# Look at number of AGB<20 in the tree cover class
lc <- readRDS("results/R_out/LC17_masked_to_ALOS_land_raster.rds")
# Mask out all but forest
msk_Tinv <- lc; rm(lc)
msk_Tinv[msk_Tinv!=4] <- 0
msk_Tinv[msk_Tinv==4] <- 1

msk_u20_T <- msk_u20*msk_Tinv
plot(msk_Tinv[1, 11000:14000, 7000:9000])

(u20_T <- sum(msk_u20_T[[1]]==1, na.rm=TRUE))
(T <- sum(msk_Tinv[[1]]==1, na.rm=TRUE))




# agb.ras[agb.ras > 310] <- NA
# agb.ras[agb.ras < 20] <- NA
# agb.ras <- raster("results/tifs_by_R/agb18_haiti_v6_0to310.tif")
agb.20to310 <- agb.ras[agb.ras < 20] <- NA

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
