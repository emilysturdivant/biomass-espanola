# *************************************************************************************************
# Script to:
#     * Calibrate regression model between AGB and backscatter 
# Proceeds:
#     * calculate_AGB.R - calculates AGB by plot from field data
#     * process_ALOS_tiles.R
# Requires:
#     * backscatter image (g0)
#     * plot polygons with AGB by plot (plots_agb)
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
max(plots_agb$AGB_ha)
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
  # resids <- data$pred - data$obs
  resids <- residuals(lm(pred ~ obs, data))
  rss <- sum(resids^2)
  n <- length(resids)
  df <- n-2
  mse <- rss / n
  c(r.squared=summary(lm(pred ~ obs, data))$r.squared,
    adj.r.squared=summary(lm(pred ~ obs, data))$adj.r.squared,
    Bias=abs(sum(resids) / n),
    RMSE=sqrt(mse),
    MAE=sum(abs(resids)) / n,
    MSE=mse,
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
model.10000x5 <- readRDS("results/R_out/CVmodel_g0nu_10000x5.rds")

cv_r <- model.10000x5$results
cv_r <- as_tibble(cbind(metric = names(cv_r), t(cv_r))) %>% 
  rename(value='1')
cv_r1 <- cv_r %>% 
  filter(!str_detect(metric, 'SD$'), !str_detect(metric, 'intercept'))
cv_r2 <- cv_r %>% 
  filter(str_detect(metric, 'SD$')) %>% 
  select(value) %>% 
  rename(SD=value)
cv_r3 <- bind_cols(cv_r1, cv_r2)
View(cv_r3)
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
