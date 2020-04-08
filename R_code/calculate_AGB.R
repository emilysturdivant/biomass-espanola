#library(silvr)
library(readr)
library(BIOMASS)
library(gridExtra)
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
summary(mstems$agb)
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
plot.sums <- mstems %>% 
  group_by(plot_no) %>%
  summarise_at(vars(dbh_cm, H, meanWD, agb), list(sum = sum)) 
plot.sums <- merge(mplots, plot.sums, by='plot_no', all=TRUE)
plot.dens <- plot.sums %>% mutate_at(vars(dbh_cm_sum, H_sum, agb_sum), `/`, y = .$plot_area)
summary(plot.dens$meanWD_sum)

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

# Merge plot AGB and backscatter data
g0_AGB <- merge(plots_agb, g0_plots, by.x='plot_no', by.y='plot_no', all=TRUE)
g0_AGB$g0ha2018 <- g0_AGB$`2018sum` / g0_AGB$area_ha
write_csv(g0_AGB, "~/GitHub/biomass-espanola/data/plots_g0nu_AGB.csv")

# just get the two columns we care about
g0.agb <- as.data.frame(cbind(g0_AGB$AGB_ha, g0_AGB$'2018mean')) %>% 
  rename(AGB = V1, backscatter =V2)

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
ols <- lm(g0.agb$AGB ~ g0.agb$backscatter, x=TRUE, y=TRUE)
summary(ols)
cov2cor(vcov(ols))
anova(ols)
coef(ols)
mae <- mean(abs(residuals(ols)))
mse <- mean((residuals(ols))^2)
rmse <- sqrt(mse) # RMSE
rss <- sum(residuals(ols)^2)
mss <- sum(residuals(ols)^2)/ols$df.residual
rse <- sqrt(rss / ols$df.residual)

# Get Spearman's rank correlation coefficient
corr <- cor.test(x=g0.agb$backscatter, y=g0.agb$AGB, method = 'spearman')
corr$estimate
# Get Pearson's rank correlation coefficient
corr <- cor.test(x=g0.agb$backscatter, y=g0.agb$AGB, method = 'pearson')
corr$estimate

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
model.10000x10$results

#---- 
# Histogram
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

# ------------------------
# Load AGB statistics at plots from output AGB map
agb_plots <- read_csv("~/GitHub/biomass-espanola/data/qgis_out/plots_agb18v4.csv")[c('plot_no', 'agb18sum', 'agb18mean')]
agb_obs_est <- merge(plots_agb, agb_plots, by.x='plot_no', by.y='plot_no', all=TRUE)
agb_obs_est$diff <- agb_obs_est$AGB_ha - agb_obs_est$agb18mean
agb_obs_est$diff_sqr <- agb_obs_est$diff**2
sqrt(mean(agb_obs_est$diff_sqr))
mean(agb_obs_est$diff)
mean(abs(agb_obs_est$diff))

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

