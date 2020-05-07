# Script where most of the processing happens
#library(silvr)
library(readr)
library(BIOMASS)
library(gridExtra)
library(rgdal)
library(caret)
library(boot)
library(raster)
library(tidyverse)

# Load data - Desktop
mstems <- read_csv("~/code/biomass-espanola/data/mstems_with_wooddensities.csv", col_types = cols(plot_no = col_integer()))
mplots <- read_csv("~/code/biomass-espanola/data/haiti_plots_meta.csv", col_types = cols(plot_no = col_integer()))
bwayo_densities <- read_csv("~/code/biomass-espanola/data/bwayo_densities_2.csv", col_types = cols(wd = col_double()))
g0_plots <- read_csv("~/code/biomass-espanola/data/plots_g0nu2018_HV.csv")
creole_df <- read_csv("~/code/biomass-espanola/data/exploded_specieslookup.csv")

# Load data - Mac
mstems <- read_csv("~/GitHub/biomass-espanola/data/haiti_data_wds2.csv")
mplots <- read_csv("~/GitHub/biomass-espanola/data/mplots_geoms.csv", col_types = cols(plot_no = col_integer()))
g0_plots <- read_csv("~/GitHub/biomass-espanola/data/plots_g0nu_HV.csv")
creole_df <- read_csv("~/GitHub/biomass-espanola/data/exploded_specieslookup.csv")
g0_fname <- "~/PROJECTS/Haiti_biomass/biota_out/g0nu_2018_HV_biotaLee.tif"

# Look at data ----
summary(mstems$dbh_cm, na.rm=TRUE)
sd(mstems$dbh_cm, na.rm=TRUE)
# 
(p <-ggplot(mstems, aes(x=dbh_cm)) + 
  geom_histogram(binwidth=5) +
  labs(y = "")+
  scale_x_continuous(name = "Diameter at breast height (cm)") +
  ggtitle("Histogram of tree diameters (N = 6,256, bin width = 5 m)")+
  theme_minimal())
mstems %>% 
  mutate(dbh_new = ifelse(dbh_cm > 75, 75, dbh_cm)) %>% 
  ggplot(aes(dbh_new)) +
  geom_histogram(binwidth = 2) +
  labs(y = "")+
  scale_x_continuous(name = "Diameter at breast height (cm)",
                     breaks = seq(0, 75, 10),
                     limits=c(0, 75)) +
  ggtitle("Histogram of tree diameters (N = 6,256, bin width = 5 m, outliers grouped at 75 cm)")+
  theme_minimal()
(p <-ggplot(mstems, aes(x=ht_m)) + 
  geom_histogram(binwidth=1) +
  labs(x = expression(paste("Height (m)")), 
       y = "")+
  ggtitle("Histogram of tree heights (N = 2,843, bin width = 1 m)")+
  theme_minimal())

# Run Chave14 equation without computeAGB() ----
# HEIGHTS, input: mstems$dbh_cm, mstems$ht_m, output: mstems$H, mstems$Hrse
interp_heights <- function(.data){
  # Create H-D model and use to fill missing heights
  HDmodel <- modelHD(
    D = .data$dbh_cm,
    H = .data$ht_m,
    method="log1", # best model is log1 based on prior checking
    useWeight = TRUE,
    drawGraph = FALSE
  )
  # Retrieve height data
  .data$H <- .data$ht_m
  .data$Hrse <- 0.5 # recommended RSE for observed heights
  filt <- is.na(.data$H)
  .data$H[filt] <- retrieveH(D = .data$dbh_cm, model = HDmodel)$H[filt]
  .data$Hrse[filt] <- HDmodel$RSE
  return(.data)
}
mstems <- mstems %>% interp_heights()

# AGB, input: mstems[c(dbh_cm, meanWD, H, plot_no)]
mstems$agb <- computeAGB(
  D = mstems$dbh_cm,
  WD = mstems$meanWD,
  H = mstems$H
)
# compute AGB(Mg) per plot
AGBplot <- summaryByPlot(
  computeAGB(
    D = mstems$dbh_cm,
    WD = mstems$meanWD,
    H = mstems$H
  ), 
  mstems$plot_no
)

# Look at data ----
summary(mstems$meanWD)
sd(mstems$meanWD, na.rm=TRUE)
summary(mstems$agb)
sd(mstems$agb, na.rm=TRUE)

# Plot histograms and density plot
mstems.filt <- mstems %>% filter(agb < 10)
p <-ggplot(mstems.filt, aes(x=agb)) + 
  geom_histogram(binwidth=0.05) +
  labs(y = "")+
  scale_x_continuous(name = "AGB") +
  ggtitle("Histogram of AGB (N = 6,256, bin width = 0.05 m)")+
  theme_minimal()
p
limO <- 2
mstems.filt %>% 
  mutate(x_new = ifelse(agb > limO, limO, agb)) %>% 
  ggplot(aes(x_new)) +
  geom_histogram(binwidth = 0.002) +
  labs(y = "")+
  scale_x_continuous(name = "AGB",
                     breaks = seq(0, 5, 0.1),
                     limits=c(0, limO)) +
  ggtitle(str_c("Histogram of AGB (N = 6,256, bin width = 0.05 m, outliers grouped at ",limO,")"))+
  theme_minimal()
p.dens.agb <- ggplot(mstems.filt, aes(x=agb)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)+ 
  labs(x = expression(paste("Aboveground biomass (MgC)")), 
       y = "Density")
p.dens.agb
ggplot(mstems.filt, aes(x=plot_no, y=agb))+
  geom_boxplot()

# Look at per plot 
plot.means <- mstems %>% 
  group_by(plot_no) %>%
  summarise_at(vars(dbh_cm, H, meanWD, agb, plot_area), list(mean = mean))
summary(plot.means$agb_mean)
sd(plot.means$meanWD_mean, na.rm=TRUE)
plot.sums <- mstems %>% 
  group_by(plot_no) %>%
  summarise_at(vars(dbh_cm, H, meanWD, agb), list(sum = sum)) 
plot.sums <- merge(mplots, plot.sums, by='plot_no', all=TRUE)
plot.dens <- plot.sums %>% mutate_at(vars(dbh_cm_sum, H_sum, agb_sum), `/`, y = .$plot_area)
summary(plot.dens$meanWD_sum)
sd(plot.dens$dbh_cm_sum, na.rm=TRUE)

#----
# Convert AGB per plot to AGB per hectare
plots_agb <- merge(mplots, AGBplot, by.x='plot_no', by.y='plot', all=TRUE)
plots_agb$AGB_ha <- plots_agb$AGB / plots_agb$area_ha
plots_agb$AGB_ha[is.na(plots_agb$AGB_ha)] <- 0

# Look at data ----
summary(plots_agb$AGB_ha)
# Plots
ggplot(plots_agb, aes(x=plot_no, y=AGB_ha))+
  geom_boxplot()
p.dens.agb <- ggplot(plots_agb, aes(x=AGB_ha)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)+ 
  labs(x = expression(paste("Aboveground biomass (MgC)")), 
       y = "Density")
p.dens.agb
p <-ggplot(plots_agb, aes(x=AGB_ha)) + 
  geom_histogram(binwidth=5) +
  scale_y_continuous(name = "", breaks = seq(0, 8, 2)) +
  labs(y = "")+
  scale_x_continuous(name = "Aboveground Biomass (Mg/ha)") +
  ggtitle("Histogram of AGB (N = 36, bin width = 5)")+
  theme_minimal()
p
filt <- plots_agb %>% filter(AGB_ha > 0)
summary(filt$AGB_ha)

mean(mstems$dbh_cm, na.rm=TRUE)
sd(mstems$dbh_cm, na.rm=TRUE)
mean(mstems$H, na.rm=TRUE)
sd(mstems$H, na.rm=TRUE)
mean(mstems$meanWD, na.rm=TRUE)
sd(mstems$meanWD, na.rm=TRUE)
mean(mstems$sdWD, na.rm=TRUE)
mean(mstems$agb, na.rm=TRUE)
sd(mstems$agb, na.rm=TRUE)

# Save/load data ----
save.image(file = "~/GitHub/biomass-espanola/data/work_space.RData") 
load("~/GitHub/biomass-espanola/data/work_space.RData")

# Get mean backscatter for each plot ---- 
# Load raster and polygon data
# Raster was created by 
# 1) download PALSAR mosaic tiles and convert HV backscatter to natural units. 
# 2) merge the tiles in QGIS, 
# 3) Run SAGA::Multi Direction Lee Filter in QGIS with both est. noise==1 and method==1
# 3.b) Possibly mask out water and urban features and extreme values. Convert to No Data using Raster Calc (0/0)
# 4) Export Filtered Grid (No Data == -99999)
g0 <- raster("~/PROJECTS/Haiti_biomass/biota_out/g0nu_2018_haiti_qLee1.tif")
g0 <- raster("~/PROJECTS/Haiti_biomass/biota_out/g0nu_2018_nofilt_HV_haiti.tif")
g0 <- raster("~/PROJECTS/Haiti_biomass/biota_out/g0nu_2018_HV_biotaLee.tif")
polys <- readOGR(dsn="~/GitHub/biomass-espanola/data", layer='AllPlots')

# Aggregate to 50m, as recommended by Saatchi 2015 and performed by Michelakis et al. 2015
# g0.nofilt <- raster("~/PROJECTS/Haiti_biomass/biota_out/g0nu_2018_nofilt_HV_haiti.tif")
# g0.nofilt[g0.nofilt == 0] <- NA
# g0.agg <- aggregate(g0.nofilt, fact=2, fun=mean, na.rm=TRUE, 
#                     filename="~/PROJECTS/Haiti_biomass/biota_out/g0nu_2018_haiti_agg50m.tif", 
#                     overwrite=TRUE)

# Extract backscatter values at plots
ex <- raster::extract(g0, polys)
polys$g0l_mean <- ex %>% 
  lapply(function(x) if (!is.null(x)) mean(x, na.rm=TRUE) else NA ) %>% 
  unlist()
polys$g0l_count <- ex %>%
  lapply(function(x) length(x)) %>% 
  unlist()
g0_plots <- polys[c('plot_no', 'g0l_mean')]

# Merge plot AGB and backscatter data
g0_AGB <- g0_plots %>% 
  merge(plots_agb, by='plot_no', all=TRUE)
writeOGR(g0_AGB, dsn="~/GitHub/biomass-espanola/data", layer='plots_g0agb', driver="ESRI Shapefile")
g0.agb <- g0_AGB[c('AGB_ha', 'g0l_mean')] %>% 
  as.data.frame() %>% 
  rename(AGB = AGB_ha, backscatter = g0l_mean)

mean(g0.agb$AGB)
sd(g0.agb$AGB)
mean(g0.agb[9:36, ]$AGB)
sd(g0.agb[9:36, ]$AGB)
g0.agb$AGB[g0.agb$AGB > 98]

mean(g0_AGB[9:36, ]$area_ha)
sd(g0_AGB[9:36, ]$area_ha)
mean(g0_AGB$area_ha)
sd(g0_AGB$area_ha)
range(g0_AGB$area_ha)

# Save/load data ----
save.image(file = "~/GitHub/biomass-espanola/data/work_space.RData") 
load("~/GitHub/biomass-espanola/data/work_space.RData")

# Histograms of plot AGB and backscatter ---- 
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

#---- 
# Scatterplot - AGB against backscatter
p <- ggplot(g0.agb, aes(x=backscatter, y=AGB)) + geom_point() +
  labs(y = expression(paste("Aboveground biomass (Mg ha"^"-1", ")")), 
       x = expression(paste("Radar backscatter, ",sigma['HV']^0," (m"^2, "/m"^2, ")")))
(p <- p + geom_smooth(method="lm", se=TRUE, fullrange=TRUE, level=0.95, col='black'))

# Linear regression ----
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

# Repeated k-fold cross validation ----
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
save(model.10000x10, file = "~/PROJECTS/Haiti_biomass/R_out/CVmodel_g0nuLee_10000x10.rds")

# Save/load data ----
save.image(file = "~/GitHub/biomass-espanola/data/work_space.RData") 
load("~/GitHub/biomass-espanola/data/work_space.RData")

# Pairs Bootstrap ----
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

# Create AGB raster ----
# Load backscatter raster
g0 <- raster(g0_fname)
names(g0) <- 'backscatter'

# Apply linear regression model to create AGB map
agb.ras <- raster::predict(g0, ols, na.rm=TRUE)
writeRaster(agb.ras, "~/PROJECTS/Haiti_biomass/R_out/agb18_v8.tif")

# Mask further
# agb.ras[agb.ras > 310] <- NA
# agb.ras[agb.ras < 20] <- NA
agb.ras <- raster("~/PROJECTS/Haiti_biomass/R_out/agb18_haiti_v6_0to310.tif")
agb.20to310 <- agb.ras[agb.ras < 20] <- NA

# Look at AGB distributions ----
agb.br <- brick(raster("~/PROJECTS/Haiti_biomass/biota_out/agb_2018_v6_mask2share.tif"),
                 raster("~/PROJECTS/Haiti_biomass/biota_out/agb_2018_v6CI_2share.tif"))
names(agb.br) <- c('AGB', 'CI')
saveRDS(agb.br, file = "~/PROJECTS/Haiti_biomass/R_out/AGB_95ci.rds")

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

# Graphing ----
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

# Import ALOS mosaic rasters (ENVI files) ----
isl_poly <- readOGR(dsn="~/PROJECTS/Haiti_biomass/contextual_data/Hispaniola", layer='Hisp_adm0')
# List filenames of ALOS mask tiles. Tiles range from 18-21, 68-75
fps <- list.files(path='~/PROJECTS/Haiti_biomass/ALOS',
           pattern='_18_mask_F02DAR$',
           full.names=TRUE,
           recursive=TRUE,
           include.dirs=FALSE)

# Load and crop rasters and store in list
load.crop.na <- function(fp, ext) try({r <- raster(fp); r <- crop(r, ext)}, silent=TRUE)
rs <- lapply(fps, load.crop.na, ext=isl_poly)
# There are two that don't seem to overlap with the polygon so throw an error. Remove them.
l <- vector(mode="logical", length(rs))
for(i in seq(1, length(rs))){
  if(class(rs[[i]])=='try-error') l[i] <- TRUE
}
rs <- rs[!l]
rs

# Merge the tiles
dn <- do.call(merge, rs)
plot(dn)
vals <- getValues(dn)
df <- data.frame(
  group = c("Normal", "Layover", "Shadowing"), 
  value = c(sum(vals==255, na.rm=TRUE), 
            sum(vals==100, na.rm=TRUE),
            sum(vals==150, na.rm=TRUE)))
df$value[2] / sum(df$value)
df$value[3] / sum(df$value)
setwd("~/PROJECTS/Haiti_biomass/R_out")
save.image()
load(".RData")


# Count inland NA value in raster ----
# Load files
g0 <- raster("~/PROJECTS/Haiti_biomass/biota_out/g0nu_2018_HV.tif")
hti_poly <- readOGR(dsn="~/PROJECTS/Haiti_biomass/contextual_data/HTI_adm", layer='HTI_adm0')
isl_poly <- readOGR(dsn="~/PROJECTS/Haiti_biomass/contextual_data/Hispaniola", layer='Hisp_adm0')

# Crop backscatter and reclass 0s to NA
g0 <- crop(g0, hti_poly)
g0[g0==0] <- NA
writeRaster(g0, "~/PROJECTS/Haiti_biomass/biota_out/g0nu_2018_HV_haitiR.tif", overwrite=TRUE)
g0 <- raster("~/PROJECTS/Haiti_biomass/biota_out/g0nu_2018_HV_haitiR.tif")

# Convert island to land mask
tmp_isl <- rasterize(isl_poly, g0, field=1) # land==1; sea==NA

# Multiply backscatter classes x polygon classes 
# reclass to (NA==NA; 1==land; 2==sea)
msk <- overlay(g0, tmp_isl, 
               fun=function(r1, r2){
                 r1[!is.na(r1)] <- 0 # values in land-->0
                 r1[is.na(r1)] <- 1 # NAs in land-->1
                 r2[is.na(r2)] <- 2 # sea-->2 (land==1)
                 r <- r1*r2 # create output (sea==2)
                 r[r==1] <- NA # NAs in land --> NA
                 r[r==0] <- 1 # valid land --> 1
                 return(r)}
               )
plot(msk) 
# define as categorical variable
f <- as.factor(msk) # or ratify(msk) ?
x <- levels(f)[[1]]
x$code <- c("land", "sea")
levels(f) <- x
writeRaster(f, "~/PROJECTS/Haiti_biomass/biota_out/mask_NAinland.tif")

# Land mask
land <- msk
land[is.na(land)] <- 1
land[land==2] <- NA
writeRaster(land, "~/PROJECTS/Haiti_biomass/biota_out/mask_land18.tif")
land <- raster("~/PROJECTS/Haiti_biomass/biota_out/mask_land18.tif")
land_poly <- rasterToPolygons(land) # very slow, much faster with gdal_polygonize in QGIS 
land_poly <- readOGR(dsn="~/PROJECTS/Haiti_biomass/biota_out/vector", layer='mask_land18_hti')
land_area <- sum(area(land_poly)) / 1000000

# Get values and count NAs
vals <- getValues(msk)
df <- data.frame(
  group = c("Land", "Nulls"), 
  value = c(sum(vals==1, na.rm=TRUE), 
            sum(is.na(vals)))
)
na_pct <- df$value[2] / (df$value[1] + df$value[2])

# Convert to area
msk[is.na(msk)] <- 0
cell_size <- area(msk, na.rm=TRUE, weights=FALSE)
cell_size <- cell_size[!cell_size==2]
land_area<-length(cell_size)*median(cell_size)
cell_size <- cell_size[is.na(cell_size)]
NA_area<-length(cell_size)*median(cell_size)

# Barplot
(bp <- ggplot(df, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity")+ 
  scale_fill_manual(values=c("#E69F00", "#999999")) +
  theme_minimal())
(pie <- bp + coord_polar("y", start=0))


# Try ggplot plotting
test_spdf <- as(msk, "SpatialPixelsDataFrame")
test_df <- as.data.frame(test_spdf)
colnames(test_df) <- c("value", "x", "y")
ggplot() +  
  geom_tile(data=test_df, aes(x=x, y=y, fill=value), alpha=0.8) + 
  geom_polygon(data=OR, aes(x=long, y=lat, group=group), 
               fill=NA, color="grey50", size=0.25) +
  scale_fill_viridis() +
  coord_equal() +
  theme_map() +
  theme(legend.position="bottom") +
  theme(legend.key.width=unit(2, "cm"))



mask_ras <- raster("~/PROJECTS/Haiti_biomass/biota_out/mask_mosaic18_nosea.tif")
mask_ras
mask_vals <- getValues(mask_ras)

# Get counts
na_ct <- sum(is.na(mask_vals))
nonan_ct <- sum(!is.na(mask_vals))
sea_ct <- sum(mask_vals==0, na.rm=TRUE)
valid_ct <- sum(mask_vals==1, na.rm=TRUE)

sea_ct + valid_ct == nonan_ct
na_ct/valid_ct


