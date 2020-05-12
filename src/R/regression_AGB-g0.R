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
library(geobgu)
library(stars)
# library(tmap)

# g0_plots <- read_csv("results/plots_values/plots_g0nu_HV.csv")
g0_fname <- "results/g0nu_HV/g0nu_2018_HV_leeBiota.tif"
plots_shp <- "results/plots_values/all_plots.shp"

# Consolidate plot polygons ---- ######################################################
# Load shp of all Biomass plots
# plots_b <- st_read("data/plots_shps/BiomassPlots.shp")
# plots_b <- plots_b[c('Name','Area_ha')] %>% 
#   st_zm() %>% 
#   st_transform(4326)
mstems <- read_csv("data/species_and_wds/haiti_data_wds2.csv")
plots_no <- mstems %>% group_by(plot_no, plot_shp, plot_area) %>% 
  summarize()
plots_no$plot_shp <- plots_no$plot_shp %>% 
  str_replace('Pandiassou2', '2pandiassou') %>% 
  str_replace('Campeche', 'Camapeche')

# List shps
# fps <- list.files(path="data/plots_shps", pattern="^P.+\\.shp$", full.names=TRUE)
fps <- list.files(path="data/plots_shps", pattern="\\.shp$", full.names=TRUE)

# Load and standardize polygons
standardize.names <- function(fp) try({
  pol <- fp %>% 
    # Read shapefile
    sf::st_read(quiet=TRUE) %>% 
    # Remove Z dimension from those that have it
    st_zm() %>% 
    # Transform to lat long
    st_transform(4326) # %>% 
    # # Standardize names
    # plyr::rename(c('Shape_Area' = 'Area_ha', 
    #                'area_ha' = 'Area_ha'), warn_missing=FALSE)
  pol$area <- st_area(pol) %>% units::set_units(value = ha) 
  # nm <- fp %>% basename() %>% str_split('\\.', simplify=TRUE)
  # pol$Name <- nm[1]
  pol$Name <- fp %>% basename()
  pol <- pol[c('Name', 'area')]
}, silent=FALSE)
plots <- fps %>% lapply(standardize.names)

# Merge plots into one DF
plots <- do.call(rbind, plots)

# Join to plot_no from mstems
plot_polys <- plots %>% full_join(plots_no, by=c('Name' ='plot_shp')) %>% 
  filter(!is.na(plot_no)) 
plot_polys %>% 
  st_write(plots_shp)

# # Merge Biomass and NoBiomass polys
# plots <- do.call(rbind, plots_nb) %>% 
#   rbind(plots_b) %>% 
#   st_write("results/plots_values/")

# Get mean backscatter for each plot ---- #############################################
# Load raster and polygon data
# Raster was created by 
# 1) download PALSAR mosaic tiles and convert HV backscatter to natural units in biota
# 2) merge the tiles in QGIS, 
# 3) Run Radar Enhanced Lee Filter in python biota. 
# 3.b) Possibly mask out water and urban features and extreme values. Convert to No Data using Raster Calc (0/0)
# 4) Export Filtered Grid (No Data == -99999)
g0 <- read_stars(g0_fname)

# Extract backscatter values at plots using sf methods
plot_polys <- st_read(plots_shp)
plot_polys <-
  plot_polys %>% mutate(
    g0l_mean = raster_extract(g0, plot_polys, fun = mean, na.rm = TRUE)
  )
plot_polys %>%
  st_set_geometry(NULL) %>%
  knitr::kable()

# Aggregate to 50m, as recommended by Saatchi 2015 and performed by Michelakis et al. 2015
# g0.nofilt <- raster("results/g0nu_HV/g0nu_2018_nofilt_HV_haiti.tif")
# g0.nofilt[g0.nofilt == 0] <- NA
# g0.agg <- aggregate(g0.nofilt, fact=2, fun=mean, na.rm=TRUE, 
#                     filename="results/g0nu_HV/g0nu_2018_haiti_agg50m.tif", 
#                     overwrite=TRUE)

# Merge plot AGB and backscatter data
plots_agb <- readRDS('results/R_out/plots_agb.rds')
g0_AGB <- plot_polys[c('plot_no', 'g0l_mean')] %>% 
  merge(plots_agb, by='plot_no', all=TRUE) %>% 
  st_write("results/plots_values/plots_g0agb.shp")

# Extract just AGB and backscatter as dataframe
g0.agb <- g0_AGB[c('AGB_ha', 'g0l_mean')] %>% 
  as.data.frame() %>% 
  rename(AGB = AGB_ha, backscatter = g0l_mean) 
g0.agb %>% 
  saveRDS("results/R_out/plots_g0agb_dfslim.rds")

# Look at values
agb_summary <- c(all_mean=mean(g0.agb$AGB),
  all_sd=sd(g0.agb$AGB),
  bio_mean=mean(g0.agb[9:36, ]$AGB),
  bio_sd=sd(g0.agb[9:36, ]$AGB),
  min=min(g0.agb$AGB),
  max=max(g0.agb$AGB))
area_summary <- c(bio_mean=mean(g0_AGB[9:36, ]$area_ha),
  bio_sd=sd(g0_AGB[9:36, ]$area_ha),
  all_mean=mean(g0_AGB$area_ha),
  all_sd=sd(g0_AGB$area_ha),
  min=min(g0_AGB$area_ha),
  max=max(g0_AGB$area_ha))
summaries <- cbind(agb_summary, area_summary)

(g0.agb$AGB[g0.agb$AGB > 98])

# Histograms of plot AGB and backscatter ---- ############################################
# AGB
p1 <-ggplot(g0.agb, aes(x=AGB)) + 
  geom_histogram(bins=30) +
  labs(x = expression(paste("Aboveground biomass (MgC ha"^"-1", ")")), 
       y = "Number of plots (N = 36)")+ 
  ylim(0, 9) 
p1
#+  geom_vline(aes(xintercept=mean(AGB)), color="black", linetype="dashed", size=.5)
# Backscatter
p2 <-ggplot(g0.agb, aes(x=backscatter)) + 
  geom_histogram() +
  labs(x = expression(
    paste("Radar backscatter, ", sigma['HV']^0, " (m"^2, "/m"^2, ")")), 
    y = "Number of plots (N = 36)")+ 
  ylim(0, 9)
p2
#+ geom_vline(aes(xintercept=mean(backscatter)), color="black", linetype="dashed", size=.5)
grid.arrange(p1, p2, nrow = 2)

# Scatterplot - AGB against backscatter ---- ##############################################
p <- ggplot(g0.agb, aes(x=backscatter, y=AGB)) + geom_point() +
  labs(y = expression(paste("Aboveground biomass (Mg ha"^"-1", ")")), 
       x = expression(paste("Radar backscatter, ",sigma['HV']^0," (m"^2, "/m"^2, ")")))
(p <- p + geom_smooth(method="lm", se=TRUE, fullrange=TRUE, level=0.95, col='black'))

# Linear regression ---- ###########################################################
# See how correlation improves without plots 7 and 8, which are probably responding to soil moisture. 
g0.agb <- g0.agb %>% filter(!(AGB==0 & backscatter>0.016))
# Basic OLS regression
ols <- lm(AGB ~ backscatter, data=g0.agb, x=TRUE, y=TRUE)
summary(ols)
confint(ols)
mse <- mean((residuals(ols))^2)
rss <- sum(residuals(ols)^2)
acc_metrics <- list(
  rmse = sqrt(mse), # RMSE
  mae = mean(abs(residuals(ols))),
  mse = mse,
  rss = rss,
  mss = sum(residuals(ols)^2)/ols$df.residual,
  rse = sqrt(rss / ols$df.residual)
)
View(acc_metrics)
(acc_metrics %>% as_tibble() %>% t())

cov2cor(vcov(ols))
anova(ols)
coef(ols)
# Get Spearman's rank correlation coefficient
corr <- cor.test(x=g0.agb$backscatter, y=g0.agb$AGB, method = 'spearman')
corr$estimate
# Get Pearson's rank correlation coefficient
corr <- cor.test(x=g0.agb$backscatter, y=g0.agb$AGB, method = 'pearson')
corr$estimate

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
model.10000x5$finalModel
save(model.10000x10, file = "results/R_out/CVmodel_g0nuLee_10000x10.rds")

# Pairs Bootstrap ---- ###########################################################
set.seed(45)
# OLS 
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
save(mstems, mplots, polys, ex, g0.agb, model.10000x10, boot.ols.100k, file="calcAGB_outputs.RData")
load("calcAGB_outputs.RData")

# Create AGB raster ---- ###########################################################
# Load backscatter raster
g0 <- raster(g0_fname)
names(g0) <- 'backscatter'

# Apply linear regression model to create AGB map
agb.ras <- raster::predict(g0, ols, na.rm=TRUE)
writeRaster(agb.ras, "results/tifs_by_R/agb18_v8.tif")

# Mask further
# agb.ras[agb.ras > 310] <- NA
# agb.ras[agb.ras < 20] <- NA
agb.ras <- raster("results/tifs_by_R/agb18_haiti_v6_0to310.tif")
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
