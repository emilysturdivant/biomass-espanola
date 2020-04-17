#library(silvr)
library(readr)
library(BIOMASS)
library(gridExtra)
library(rgdal)
library(caret)
library(boot)
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
g0_AGB <- read_csv("~/GitHub/biomass-espanola/data/plots_g0nu_AGB.csv")

#----
# Look at data
summary(mstems$dbh_cm, na.rm=TRUE)
sd(mstems$dbh_cm, na.rm=TRUE)
p <-ggplot(mstems, aes(x=dbh_cm)) + 
  geom_histogram(binwidth=5) +
  labs(y = "")+
  scale_x_continuous(name = "Diameter at breast height (cm)") +
  ggtitle("Histogram of tree diameters (N = 6,256, bin width = 5 m)")+
  theme_minimal()
p
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

p <-ggplot(mstems, aes(x=ht_m)) + 
  geom_histogram(binwidth=1) +
  labs(x = expression(paste("Height (m)")), 
       y = "")+
  ggtitle("Histogram of tree heights (N = 2,843, bin width = 1 m)")+
  theme_minimal()
p
# ------------------------
# Run Chave14 equation without computeAGB() 
# HEIGHTS, input: mstems$dbh_cm, mstems$ht_m, output: mstems$H, mstems$Hrse
interp_heights <- function(.data){
  # Create H-D model and use to fill missing heights
  HDmodel <- modelHD(
    D = .data$dbh_cm,
    H = .data$ht_m,
    method="log1", # best model is log1 based on prior checking
    useWeight = TRUE,
    drawGraph = TRUE
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

# Look at data
summary(mstems$meanWD)
sd(mstems$meanWD, na.rm=TRUE)
summary(mstems$agb)
sd(mstems$agb, na.rm=TRUE)
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

# ---- 
# Get mean backscatter for each plot
# Load raster and polygon data
# Raster was created by 
# 1) download PALSAR mosaic tiles and convert HV backscatter to natural units. 
# 2) merge the tiles in QGIS, 
# 3) Run SAGA::Multi Direction Lee Filter in QGIS with both est. noise==1 and method==1
# 4) Export Filtered Grid (No Data == -99999)
g0 <- raster("~/PROJECTS/Haiti_biomass/biota_out/g0nu_2018_haiti_qLee1.tif")
polys <- readOGR(dsn="~/GitHub/biomass-espanola/data", layer='AllPlots')

# Aggregate to 50m, as recommended by Saatchi 2015 and performed by Michelakis et al. 2015
# g0.nofilt <- raster("~/PROJECTS/Haiti_biomass/biota_out/g0nu_2018_nofilt_HV_haiti.tif")
# g0.nofilt[g0.nofilt == 0] <- NA
# g0.agg <- aggregate(g0.nofilt, fact=2, fun=mean, na.rm=TRUE, 
#                     filename="~/PROJECTS/Haiti_biomass/biota_out/g0nu_2018_haiti_agg50m.tif", 
#                     overwrite=TRUE)

# Extract backscatter values at plots
ex <- extract(g0, polys)
polys$g0l_mean <- unlist(lapply(ex, function(x) if (!is.null(x)) mean(x, na.rm=TRUE) else NA ))
polys$g0l_count <- unlist(lapply(ex, function(x) length(x)))
View(polys)
g0_plots <- polys[c('plot_no', 'g0l_mean')]

# Merge plot AGB and backscatter data
g0_AGB <- merge(plots_agb, g0_plots, by.x='plot_no', by.y='plot_no', all=TRUE)
write_csv(g0_AGB, "~/GitHub/biomass-espanola/data/plots_g0nuLee_AGB.csv")
g0.agb <- g0_AGB[c('AGB_ha', 'g0l_mean')] %>% 
  rename(AGB = AGB_ha, backscatter = g0l_mean)

save(mstems, mplots, polys, ex, g0.agb, model.10000x10, boot.ols.100k, file="calcAGB_outputs.RData")
load("calcAGB_outputs.RData")
writeOGR(polys, dsn="~/GitHub/biomass-espanola/data", layer='plots_agb', driver="ESRI Shapefile")
mplots2 <- cbind(mplots, g0.agb)
poly2 <- merge(polys, mplots2, by='plot_no')
names(poly2)
poly2 <- poly2[-seq(2,18)]
writeOGR(poly2, dsn="~/GitHub/biomass-espanola/data", layer='plots_agb', driver="ESRI Shapefile")


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

#---- 
# Scatterplot
# Plot AGB against backscatter
p <- ggplot(g0.agb, aes(x=backscatter, y=AGB)) + geom_point() +
  labs(y = expression(paste("Aboveground biomass (MgC ha"^"-1", ")")), 
       x = expression(paste("Radar backscatter, ",sigma['HV']^0," (m"^2, "/m"^2, ")")))
p <- p + geom_smooth(method="lm", se=TRUE, fullrange=TRUE, level=0.95, col='black')
p
p + labs(caption = 'OLS regression')

# Linear model
# Basic OLS regression
ols <- lm(AGB ~ backscatter, data=g0.agb, x=TRUE, y=TRUE)
summary(ols)
confint(ols)
mae <- mean(abs(residuals(ols)))
mse <- mean((residuals(ols))^2)
rmse <- sqrt(mse) # RMSE
rss <- sum(residuals(ols)^2)
mss <- sum(residuals(ols)^2)/ols$df.residual
rse <- sqrt(rss / ols$df.residual)
View(c(rmse, mae))
cov2cor(vcov(ols))
anova(ols)
coef(ols)

# Get Spearman's rank correlation coefficient
corr <- cor.test(x=g0.agb$backscatter, y=g0.agb$AGB, method = 'spearman')
corr$estimate
# Get Pearson's rank correlation coefficient
corr <- cor.test(x=g0.agb$backscatter, y=g0.agb$AGB, method = 'pearson')
corr$estimate

# 95% CI for data points
stats::predict(ols, newdata=(backscatter=0.5), interval="confidence", level = 0.95)
names(g0.lee) <- 'backscatter'
agb.95ci <- raster::predict(g0.lee, ols, 
                            fun='predict.se', 
                            #filename="~/PROJECTS/Haiti_biomass/R_out/agb_95ci.tif", 
                            na.rm=TRUE)
agb.95ci <- raster::predict(g0.lee, ols, 
                            fun=function(x){stats::predict(x, interval="confidence") if (!is.null(x)) x*m+b else NA}, 
                            #filename="~/PROJECTS/Haiti_biomass/R_out/agb_95ci.tif", 
                            na.rm=TRUE)

agb.95ci

# Plot error plots
opar <- par(mfrow = c(2,2), oma = c(0, 0, 1.1, 0))
plot(ols, las = 1)

# Run cross validation (10,000 x 10-fold)
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
    B=sum(resids) / n,
    RSS=rss,
    MSS=rss/df,
    RSE=sqrt(rss / df))
}
set.seed(45)
model.10000x10 <- train(AGB ~ backscatter, data = g0.agb, method = "lm",
                        trControl = trainControl(method = "repeatedcv", 
                                                 number = 10, repeats = 10000,
                                                 summaryFunction = fxn.bias))
View(t(model.10000x10$results))
save(model.10000x10, file = "~/PROJECTS/Haiti_biomass/R_out/CVmodel_g0nuLee_10000x10.rds")

# Pairs Bootstrap
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
boot.ci(boot.ols.100k, conf=0.95, type=c("basic", "bca", "perc"), index=1)
boot.ci(boot.ols.100k, conf=0.95, type=c("basic", "bca", "perc"), index=2)

#---- 
# Histogram
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
#+ geom_vline(aes(xintercept=mean(backscatter)), color="black", linetype="dashed", size=.5)
grid.arrange(p1, p2, nrow = 2)

#----
# Z-score histograms
Zscores_backscatter <- scale(g0_AGB$`2018mean`)
hist(Zscores_backscatter)
Zscores_AGB <- scale(g0_AGB$AGB_ha)
hist(Zscores_AGB)
# Interleaved histograms - not working
d1 <- cbind(Zscores_AGB, 'AGB')
d2 <- cbind(Zscores_backscatter, 'backscatter')
dat1 <- as.data.frame(rbind(d1, d2)) %>% rename(Zscore = V1, Variable =V2)
p<-ggplot(dat1, aes(x=Zscore, color=Variable)) +
  geom_histogram(fill="white", position="dodge", binwidth = 0.01)+
  theme(legend.position="top")
p


# Remove plot 16 and check correlation
g0_AGB2 <- g0_AGB[-c(16), ]
# Get Spearman's rank correlation coefficient
corr <- cor.test(x=g0_AGB2$`2018mean`, y=g0_AGB2$AGB_ha, method = 'spearman')
corr$estimate
# Plot AGB against backscatter
plot(g0_AGB2$`2018mean`, g0_AGB2$AGB_ha, xlab='2018 HV backscatter (plot median)', ylab='2019 AGB (tC/ha) from H-D model')
linreg <- lm(g0_AGB2$AGB ~ g0_AGB2$`2018mean`)
abline(lm(g0_AGB2$AGB ~ g0_AGB2$`2018mean`), col='blue')
linreg
abline(lm(g0_AGB$AGB ~ g0_AGB$`2018mean`), col='red')

#----
# Create AGB raster 
g0 <- raster("~/PROJECTS/Haiti_biomass/biota_out/g0nu_2018_haiti_qLee1.tif")
names(g0) <- 'backscatter'
agb.ras <- raster::predict(g0, ols, na.rm=TRUE)
agb.ras[agb.ras > 310] <- NA
agb.ras[agb.ras < 0] <- NA
writeRaster(agb.ras, "~/PROJECTS/Haiti_biomass/R_out/agb18_haiti_v6_0to310.tif")

# Run with cropped section
predict.ci.raster <- function(g0, model){
  names(g0) <- 'backscatter'
  df <- as.data.frame(g0)
  chunks <- split(df, (seq(nrow(df))-1) %/% 1000000) 
  agblist <- list()
  for (i in 1:length(chunks)){
    chunk <- chunks[[i]]
    agb1 <- stats::predict(model, newdata=chunk, 
                           interval="confidence", level = 0.95) %>%
      as.data.frame()
    agblist[[i]] <- agb1
  }
  agb <- bind_rows(agblist)
  agb.ci.v <- as.vector((agb$upr - agb$lwr)/2)
  ci.ras <- setValues(g0, values=agb.ci.v)
}

ras.agb.ci <- predict.ci.raster(g0, ols)
ras.agb.ci[agb.ras > 310] <- NA
ras.agb.ci[agb.ras < 0] <- NA
writeRaster(ras.agb.ci, "~/PROJECTS/Haiti_biomass/R_out/agb_CI_sub310.tif", overwrite=TRUE)


#----
# Look at AGB distributions
agb.ras <- raster("~/PROJECTS/Haiti_biomass/biota_out/agb_2018_v6_mask2share.tif")
ras.agb.ci <- raster("~/PROJECTS/Haiti_biomass/biota_out/agb_2018_v6CI_2share.tif")

#----
lc.br <- brick(raster("~/PROJECTS/Haiti_biomass/biota_out/agb_2018_v6_mask2share.tif"), 
               raster("~/PROJECTS/Haiti_biomass/biota_out/agb_2018_v6CI_2share.tif"),
               raster("~/PROJECTS/Haiti_biomass/LULC/Haiti2017_water_agb18v6.tif"), 
               raster("~/PROJECTS/Haiti_biomass/LULC/Haiti2017_urban_agb18v6.tif"), 
               raster("~/PROJECTS/Haiti_biomass/LULC/agb_lc3.tif"),
               raster("~/PROJECTS/Haiti_biomass/LULC/agb_lc4.tif"),
               raster("~/PROJECTS/Haiti_biomass/LULC/agb_lc5.tif"), 
               raster("~/PROJECTS/Haiti_biomass/LULC/agb_lc6.tif"),
               raster("~/PROJECTS/Haiti_biomass/LULC/agb_lc_over3.tif"))
names(lc.br) <- c('Haiti', 'CI', 'Water', 'Urban', 'Bareland', 'Tree cover', 'Grassland', 'Shrubs', 'Veg')
saveRDS(lc.br, file = "~/PROJECTS/Haiti_biomass/R_out/brick_AGBv6_withLC.rds")
lc.br <- readRDS(file = "~/PROJECTS/Haiti_biomass/R_out/rasterbrick_AGBv6byLC.rds")
lc.br

# Get selection of percentiles for each LC
agb.qs <- data.frame(row.names=c('0%', '1%','2%','10%', '25%', '50%', '75%', '90%','98%', '99%','100%'))
for (i in seq(1,nlayers(lc.br))){
  lc <- names(lc.br[[i]])
  agb.qs[[lc]] <- quantile(lc.br[[i]], probs=c(0, 0.01, 0.02, 0.1, 0.25, 0.5, 0.75, 0.9,0.98, 0.99, 1))
}
saveRDS(agb.qs, file = "~/PROJECTS/Haiti_biomass/R_out/agb_quantiles_v6.rds")
agb.qs <- readRDS(file = "~/PROJECTS/Haiti_biomass/R_out/agb_quantiles_v6.rds")
View(agb.qs)

# Get means, SDs, and Skews
stats <- list(mean=cellStats(lc.br, stat='mean', na.rm=TRUE),
     sd=cellStats(lc.br, stat='sd', na.rm=TRUE),
     skew=cellStats(lc.br, stat='skew', na.rm=TRUE)) %>%
  bind_cols()%>%
  t()
colnames(stats) <- names(lc.br)
agb.stats <- rbind(agb.qs, stats)

saveRDS(agb.stats, file = "~/PROJECTS/Haiti_biomass/R_out/agb_stats_v6.rds")
agb.stats <- readRDS(file = "~/PROJECTS/Haiti_biomass/R_out/agb_stats_v6.rds")



# Plot density
p <- ggplot(agb.samp, aes(x=AGB)) + 
  geom_histogram(aes(y=..density..), fill="#69b3a2", color="#e9ecef", alpha=0.7, binwidth = 2.5)+
  geom_density(alpha=.2) + 
  scale_x_continuous(name = expression(paste("Aboveground biomass (MgC ha"^"-1", ")")),
                     breaks = seq(0, 150, 25),
                     limits=c(-5, 150)) +
  geom_vline(aes(xintercept=mean(AGB)), color="black", linetype="dashed", size=.5)+
  geom_vline(aes(xintercept=median(AGB)),
             color="black", linetype="dashed", size=.5)+
  theme_minimal()
p


# Get mean and median
mu <- agb.samp.haiti %>% summarise(mean = mean(AGB), med=median(AGB))
dp <- ggplot(agb.samp.haiti, aes(x = AGB)) +
  geom_density(position="identity", alpha=0.2) +
  scale_x_continuous(name = "AGB",
                     breaks = seq(0, 150, 25),
                     limits=c(0, 150)) +
  scale_y_continuous(name = "Density", limits = c(0, 0.03)) +
  ggtitle("Density plot of AGB")+
  #geom_vline(data=mu, aes(xintercept=med, color=Category),
  #linetype="dashed")+
  geom_point(data=mu, aes(med, 0), size=2)+
  geom_point(data=mu, aes(mean, 0), size=2)
dp


# Histogram and density plots
p1 <- ggplot(agb.samp, aes(x=AGB)) + 
  geom_histogram(col='white', fill="gray", binwidth = 5)+
  labs(x = expression(paste("Aboveground biomass (MgC ha"^"-1", ")")), 
       y = "Pixel count (out of 100,000 pixel sample)")+
  geom_vline(aes(xintercept=mean(AGB)),
             color="black", linetype="dashed", size=.5)+
  geom_vline(aes(xintercept=median(AGB)),
             color="black", linetype="dashed", size=.5)
p1
p.dens.agb <- ggplot(agb.samp, aes(x=AGB)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.7)+ 
  labs(x = expression(paste("Aboveground biomass (MgC ha"^"-1", ")")), 
       y = "Density")+
  scale_x_continuous(name = expression(paste("Aboveground biomass (MgC ha"^"-1", ")")),
                     breaks = seq(0, 150, 25),
                     limits=c(-5, 150)) +
  geom_vline(aes(xintercept=mean(AGB)), color="black", linetype="dashed", size=.5)+
  geom_vline(aes(xintercept=median(AGB)),
             color="black", linetype="dashed", size=.5)+
  theme_minimal()
p.dens.agb

mean(agb.samp$AGB)
median(agb.samp$AGB)






# ------------------------
# Run Chave14 equation using computeAGB() in BIOMASS
# Get wood density by genus
species <- NULL
if (is.null(species)) species = rep('', length(mstems$genus))
family <- NULL
#region <- 'CentralAmericaTrop'
latitude <- 19
longitude <- -72

# Prepare coordinates, required without height
if (!is.null(latitude)) coord = data.frame(longitude = longitude, latitude = latitude) else coord = NULL
height=NULL

# Calculate AGB with Chave equation, return in Mg
mstems$AGB  <- computeAGB(mstems$dbh_cm, mstems$meanWD, H = height, 
                          coord = cbind(mstems$lon, mstems$lat), Dlim = 0)
summary(mstems$AGB)
mstems[which(mstems$AGB == max(mstems$AGB, na.rm=T)), ]
boxplot(AGB~plot_no, data=mstems)
boxplot(AGB~sp_creole, data=mstems)
plot(mstems$dbh_cm, mstems$AGB, xlab='diam', ylab='AGB')

# aggregate AGB by plot (MgC/ha)
agb_plot <- aggregate(AGB ~ plot_no, mstems, sum)
mplots <- merge(mplots, agb_plot, by='plot_no', all=TRUE)
mplots$AGB <- mplots$AGB / mplots$area_ha

# Add AGB to
g0_plots <- merge(mplots, g0_plots, by='plot_no', all=TRUE)
#write.csv(g0_plots, "~/GitHub/biomass-espanola/plots_g0nu2018_withAGB.csv")

# Plot AGB against backscatter
plot(g0_plots$`2018_mean`, g0_plots$AGB, xlab='2018 HV backscatter', ylab='2019 AGB (tC/ha)')
linreg <- lm(g0_plots$AGB ~ g0_plots$`2018_mean`)
abline(linreg)
linreg

# Get Spearman's rank correlation coefficient
corr <- cor.test(x=g0_plots$`2018_mean`, y=g0_plots$AGB, method = 'spearman')
corr$estimate

#-----------------------
library(silvr)
median(mplots$area_ha)

# aggregate Basal Area by plot (m^2/ha) - the area of land that is occupied by the cross-section of a tree.
# calculate basal area (m^2) from DBH (cm)
mstems$basal_area <- calculateBasalArea(mstems$dbh_cm)
ba_plot <- aggregate(basal_area ~ plot_no, mstems, sum)
mplots <- merge(mplots, ba_plot, by='plot_no', all=TRUE)
mplots$basal_area <- mplots$basal_area / mplots$area_ha

# calculate stem volume  (m^3) from DBH (cm) - an estimate of the ammount of space that the tree bole occupies.
mstems$volume <- calculateStemVolume(mstems$dbh_cm)
vol_plot <- aggregate(volume ~ plot_no, mstems, sum)
mplots <- merge(mplots, vol_plot, by='plot_no', all=TRUE)
mplots$volume <- mplots$volume / mplots$area_ha

# calculate stocking density () from DBH
mstems$stocking <- calculateStocking(mstems$dbh_cm)
mstems$stocking_10 <- calculateStocking(mstems$dbh_cm, min_diam = 5)
stocking_plot <- aggregate(stocking ~ plot_no, mstems, sum)
mplots <- merge(mplots, stocking_plot, by='plot_no', all=TRUE)
mplots$stocking <- mplots$stocking / mplots$area_ha

# Look at dominant species based on stocking densities. Add species_name (by_binomial)
dominant_species <- getDominantSpecies(mstems$sp_creole, mstems$plot_no)
dominant_species <- getDominantSpecies(mstems$sp_creole, mstems$plot_no, abundance = calculateBasalArea(mstems$dbh_cm))
dominant_species <- getDominantSpecies(mstems$sp_creole, mstems$plot_no, abundance = mstems$AGB_Chave14)

