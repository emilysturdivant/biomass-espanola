# Calibrate regression model between AGB and backscatter 
# Follows calculate_AGB.R, which calculates AGB by plot from field data
# Requires backscatter image, plot polygons, AGB by plot
library(readr)
library(gridExtra)
library(rgdal)
library(caret)
library(boot)
library(raster)
library(tidyverse)
library(sf)
library(stars)
library(geobgu)
library(broom)
# library(tmap)

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

# Scatterplot - AGB against backscatter ---- ##############################################
g0.agb$AGB <- as.vector(g0.agb$AGB)
p <- ggplot(g0.agb, aes(x=backscatter, y=AGB)) + geom_point() +
  labs(y = expression(paste("Aboveground biomass (Mg ha"^"-1", ")")), 
       x = expression(paste("Radar backscatter, ",sigma['HV']^0," (m"^2, "/m"^2, ")")))
(p <- p + geom_smooth(method="lm", se=TRUE, fullrange=TRUE, level=0.95, col='black'))

# Linear regression ---- ###########################################################
# See how correlation improves without plots 7 and 8, which are probably responding to soil moisture. 
# g0.agb <- g0.agb %>% filter(!(AGB==0 & backscatter>0.016))

# Basic OLS regression
ols <- lm(AGB ~ as.vector(backscatter), data=g0.agb, x=TRUE, y=TRUE)

# Report results ---- ###########################################################
report_ols_results <- function(ols){
  # Regression Coefficients
  coefs <- tidy(ols) %>% cbind(confint(ols)) %>% 
    mutate(term = str_replace(term, '\\(Intercept\\)', 'int')) %>% 
    mutate(term = str_replace(term, 'backscatter', 'g0')) %>% 
    pivot_wider(names_from=term, values_from=estimate:'97.5 %') %>% 
    cbind(glance(ols))
  coefs %>% t()
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
  )
  # ANOVA
  ANOVA <- anova(ols) %>% tidy() %>% 
    mutate(term = str_replace(term, 'backscatter', 'g0')) %>% 
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
    t()
  # Combine
  all_vals <- cbind(coefs, acc_metrics, ANOVA, cors) %>% t()
  select_vals <- cbind(
    coefs %>% select(starts_with('estimate')),
    ci_int=(confint(ols)[1,2]-confint(ols)[1,1])/2,
    ci_slope=(confint(ols)[2,2]-confint(ols)[2,1])/2,
    acc_metrics %>% as.data.frame() %>% select(rmse, mae, mse, rss, mss, rse)
  )
  return(list(select_vals, all_vals))
}
vals <- report_ols_results(ols)
vals[1]
vals[2]

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

# Pairs Bootstrap ---- ###########################################################
# OLS 
set.seed(45)
boot.ols.100k <- boot(g0.agb, function(data=g0.agb, index) {
  data <- data[index,] # we sample along rows of the data frame
  model.boot <- lm(AGB ~ backscatter, data=data)
  coef(model.boot)
}, R=100000)
# Results
boot.ols.100k
plot(boot.ols.100k, index=1)
cis <- list()
ci <- boot.ci(boot.ols.100k, conf=0.95, type=c("basic", "bca", "perc"), index=1)
cis[['b']] <- ci$bca[4:5]
ci <- boot.ci(boot.ols.100k, conf=0.95, type=c("basic", "bca", "perc"), index=2)
cis[['m']] <- ci$bca[4:5]
cis <- as.data.frame(cis, row.names = c('lwr', 'upr'))

# Save/Load outputs from before creating AGB raster
boot.ols.100k %>% saveRDS("results/R_out/boot_g0nu_100k.rds")

# Create AGB raster ---- ###########################################################
# Load backscatter using raster (vs. stars)
# g0 <- raster(g0_fname)
# names(g0) <- 'backscatter'
# agb.ras <- raster::predict(g0, ols, na.rm=TRUE)

# Apply linear regression model to create AGB map
names(g0) <- 'backscatter'
agb.ras <- raster::predict(g0, ols, na.rm=TRUE)
write_stars(agb.ras, "results/tifs_by_R/agb18_v1_l0.tif")

# Mask
# Report starting number of NAs
# Look at proportion of values in each mask category
dn_mask <- raster('results/tifs_by_R/hisp18_mask.tif')
dn_mask <- read_stars('results/tifs_by_R/hisp18_mask.tif')
vals <- getValues(dn_mask)
f <- as.factor(vals)
levels(f)
df <- data.frame(
  group = c("Normal", "Layover", "Shadowing"), 
  value = c(sum(vals==255, na.rm=TRUE), 
            sum(vals==100, na.rm=TRUE),
            sum(vals==150, na.rm=TRUE)))
df$value[2] / sum(df$value)
df$value[3] / sum(df$value)
saveRDS(df, 'results/R_out/mask_pcts.rds')
df <- readRDS('results/R_out/mask_pcts.rds')
msk_land_hisp <- read_stars('results/masks/hisp18_maskLand.tif')
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
