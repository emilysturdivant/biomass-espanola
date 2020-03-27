#library(silvr)
library(readr)
library(BIOMASS)
library(ggplot2)

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

# ------------------------
# Run Chave14 equation without computeAGB() 
# best model is log1 so we use that to create the model
HDmodel <- modelHD(
  D = mstems$dbh_cm,
  H = mstems$ht_m,
  method="log1",
  useWeight = TRUE,
  drawGraph = TRUE
)
# Retrieve height data
mstems$Hmix <- mstems$ht_m
mstems$RSEmix <- 0.5
filt <- is.na(mstems$Hmix)
mstems$Hmix[filt] <- retrieveH(D = mstems$dbh_cm, model = HDmodel)$H[filt]
mstems$RSEmix[filt] <- HDmodel$RSE
# Apply AGBmonteCarlo --- not working
resultMC <- AGBmonteCarlo(
  D = mstems$dbh_cm, WD = mstems$meanWD, errWD = mstems$sdWD,
  H = mstems$Hmix, errH = mstems$RSEmix,
  Dpropag = "chave2004"
)
# Try performing the work of summaryByPlot() to troubleshoot/workaround the error.
library(data.table)
AGB_val <- resultMC$AGB_simu
is.list(AGB_val)
Plot <- data.table(plot = mstems$plot_no)
indice_tree <- Plot[is.na(plot), .I, by = plot]
# filter if there is there is NA in the AGB_val
filter <- rowSums(is.na(AGB_val)) > 0
# take the first tree in the database by plot
indice_first <- Plot[!is.na(plot) & !filter, .(indice_line = first(.I), plot = unique(plot)), by = plot]
nrow(indice_tree)
resAGB <- colSums(AGB_val[.I, ], na.rm = T)
AGB_summary <- (list(
  AGB = mean(resAGB),
  Cred_2.5 = quantile(resAGB, probs = 0.025),
  Cred_97.5 = quantile(resAGB, probs = 0.975)
))
Res <- as.data.frame(Plot[!is.na(plot), AGB_summary, by = plot])
View(AGB_summary)
View(resAGB)
AGB_val[.I, ]
colSums(AGB_val[.I, ], na.rm = F)
# drawPlot <- TRUE
# if (drawPlot) {
#   with(AGB[order(AGB)], {
#     plot(AGB,
#          pch = 20, xlab = "", ylab = "AGB (Mg/ha)", ylim = range(Cred_2.5, Cred_97.5),
#          las = 1, cex.lab = 1.3, xaxt = "n", main = "AGB by plot"
#     )
#     axis(1, at = seq(length(AGB)), labels = plot, las = 2)
#     segments(x0 = seq(length(AGB)), y0 = Cred_2.5, y1 = Cred_97.5, col = "red")
#   })
# }

# The rest of the section isn't working; contact Maxime (maxime.rejou@gmail.com)?
Res <- summaryByPlot(resultMC$AGB_simu, mstems$plot_no)
Res <- Res[order(Res$AGB), ]
plot(Res$AGB, pch = 20, xlab = "Plots", ylab = "AGB (Mg/ha)", ylim = c(0, max(Res$Cred_97.5)), las = 1, cex.lab = 1.3)
segments(1:nrow(Res), Res$Cred_2.5, 1:nrow(Res), Res$Cred_97.5, col = "red")

# compute AGB(Mg) per plot
AGBplot <- summaryByPlot(
  computeAGB(
    D = mstems$dbh_cm,
    WD = mstems$meanWD,
    H = mstems$Hmix
  ), 
  mstems$plot_no
)

# Convert AGB per plot to AGB per hectare
plots_agb <- merge(mplots, AGBplot, by.x='plot_no', by.y='plot', all=TRUE)
plots_agb$AGB_ha <- plots_agb$AGB / plots_agb$area_ha
plots_agb$AGB_ha[is.na(plots_agb$AGB_ha)] <- 0
summary(plots_agb$AGB_ha)

# Plot AGB against backscatter
g0_AGB <- merge(plots_agb, g0_plots, by.x='plot_no', by.y='plot_no', all=TRUE)
g0_AGB$g0ha2018 <- g0_AGB$`2018sum` / g0_AGB$area_ha
plot(g0_AGB$`2018mean`, g0_AGB$AGB_ha, xlab='2018 HV backscatter (g0 nu mean)', ylab='2019 AGB (tC/ha)')
linreg <- lm(g0_AGB$AGB_ha ~ g0_AGB$`2018mean`, x=TRUE, y=TRUE)
abline(linreg)
linreg$coefficients
linreg$residuals
summary(linreg)
anova(linreg)
confint(linreg, level=.95)
sqrt(var(linreg$residuals))
library(normwhn.test)
normality.test1(cbind(mstems$AGB_ha, mstems$`2018mean`))
# Use Model II regression
library(lmodel2)
lm2 <- lmodel2(AGB_ha ~ `2018mean`, g0_AGB, range.y='interval', range.x='relative')
summary(lm2)
lm2$regression.results
lm2$confidence.intervals
lm2$rsquare
lm2$H
lm2$r
lm2$P.param
lm2$eigenvalues
# 4 figures arranged in 2 rows and 2 columns
par(mfrow=c(2,2))
plot.lmodel2(lm2, 'OLS')
plot.lmodel2(lm2, 'MA')
plot.lmodel2(lm2, 'SMA')
plot.lmodel2(lm2, 'RMA')

# Get Spearman's rank correlation coefficient
corr <- cor.test(x=g0_AGB$`2018mean`, y=g0_AGB$AGB_ha, method = 'spearman')
corr$estimate

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
agb_obs_est$AGB_ha[is.na(agb_obs_est$AGB_ha)] <- 0
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
