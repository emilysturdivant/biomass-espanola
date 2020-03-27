#library(silvr)
library(readr)
library(BIOMASS)
library(ggplot2)
library(dplyr)
require(boot)

# Load data - Desktop
# Load data - Mac
g0_AGB <- read_csv("~/GitHub/biomass-espanola/data/plots_g0nu_HV.csv")

#----
# Scatterplot
plot(g0_AGB$`2018mean`, g0_AGB$AGB_ha, xlab='2018 HV backscatter (g0 nu mean)', ylab='2019 AGB (tC/ha)')

# Basic OLS regression
linreg <- lm(g0_AGB$AGB_ha ~ g0_AGB$`2018mean`, x=TRUE, y=TRUE)
abline(linreg)
summary(linreg)
confint(linreg)
cov2cor(vcov(linreg))
anova(linreg)
coef(linreg)

# just get the two columns we care about
dat1 <- as.data.frame(cbind(g0_AGB$AGB_ha, g0_AGB$'2018mean')) %>% 
  rename(AGB = V1, backscatter =V2)

#----
# Bootstrap options - manual
# from "Using the non-parametric bootstrap for regression models in R" by Ian Dworkin
# Non-parametric bootstrap: Pairs (Random x) approach
N = 10000 # Perform N bootstrap iterations
BootstrapRandomX <- function(dat=dat1, mod.formula=formula(AGB ~ backscatter)){
  dat.boot <- dat[sample(x = nrow(dat), size = nrow(dat), replace=T),] # samples along index
  boot.lm <- lm(mod.formula, data=dat.boot)
  coef(boot.lm)
}
vector.boot <- t(replicate(N, BootstrapRandomX()))
# standard error of the estimates via bootstrap 
# (standard deviations of those distributions)
apply(vector.boot, MARGIN = 2, sd)
# Percentile CIs (transposed to compare to simple confints)
t(apply(vector.boot, MARGIN = 2, quantile, probs=c(0.025, 0.975)))
# Histogram of distributions from Pairs bootstrap
par(mfrow=c(1,2))
MultipleHistograms <- function(X=vector.boot){
  for (i in 1:ncol(X)){
    hist(X[,i], freq=F,
         main=colnames(X)[i],
         xlab=colnames(X)[i])
  }
}
MultipleHistograms()
pairs(vector.boot)

# Use Bias-Corrected (BC) and accelerated (a) non-parametric bootstrap 
# confidence intervals (BCa) to adjust for biases in the Percentile Confidence 
# Intervals. We can calculate them with boot() in the boot library.

# Non-parametric bootstrap: Residual (Fixed effect / Experimental) approach
# Analogous analysis to Monte Carlo simulations to generate confidence intervals
# 1) fit model as normal and 2) get residuals, 
# 3) bootstrap the residuals from the model (r*)
# 4) add r* back onto fitted component of the model (i.e. b*x[i] + r*[i])
resid.model.1 <- resid(linreg)
par(mfrow=c(2,1))
plot(density(resid.model.1, bw=0.5))
plot(density(resid.model.1, bw=1))

par(mfrow=c(1,2))
plot(resid.model.1 ~ dat1$AGB)
plot(resid.model.1 ~ dat1$backscatter)

BootstrapFromResiduals <- function(mod.object = linreg, dat = dat1) {
  resids = mod.object$resid
  fittedValues = mod.object$fitted
  matr <- model.matrix(mod.object)
  # generating new values for each y[i], by adding bootstrapped resids to fitted values.
  Y <- fittedValues + sample(resids, length(resids), replace=T)
  # Using model.matrix for the predictors 
  model.boot <- lm(Y ~ 0 + matr, data=dat)
  coef(model.boot) # Extract coefficients
}
# Run and look at it.
residual.boot.N <- t(replicate(N, BootstrapFromResiduals()))
par(mfrow=c(1,2))
MultipleHistograms(X=residual.boot.N)
pairs(residual.boot.N)
apply(residual.boot.N, MARGIN = 2, sd)
t(apply(residual.boot.N, MARGIN=2, quantile, probs=c(0.025, 0.975)))

# Monte Carlo bootstrap (parametric) - simulating values in the response
SimulationUnderModel <- function(model = linreg) {
  # extract design matrix
  matr <- model.matrix(model)
  rse = summary(model)$sigma
  df = model$df
  # incorporate uncertainty in RSE
  rse.sim <- rse*sqrt(df/rchisq(1, df=df))
  # Simulate data (response) conditional on the simulated RSE.
  y.sim <- rnorm(n = nrow(matr), 
                 mean=matr%*%coef(model), sd=rse.sim)
  # 0 + design matrix (since the intercept is already in the design matrix)
  lm.sim <- lm(y.sim ~ 0 + matr) # fit model with simulated response
  coef(lm.sim)
}
# Run and look at it
sim.coef <- t(replicate(N, SimulationUnderModel()))
apply(sim.coef, MARGIN = 2, sd)
t(apply(sim.coef, MARGIN=2, quantile, probs=c(0.025, 0.975)))
par(mfrow=c(1,2))
MultipleHistograms(X=sim.coef)

# Compare
par(mfrow=c(2,1))
plot(density(residual.boot.N[,2], bw=10), 
     main="Comparing bootstrap methods for parameter uncertainty: backscatter", 
     lwd=2, lty=1)
lines(density(vector.boot[,2], bw=10), col='red', lwd=2, lty=1)
lines(density(sim.coef[,2], bw=10), col='purple', lwd=2, lty=1)
#legend('topright', legend=c("Residual Boot", "Pairs Boot", "Monte Carlo Normal"), 
#       col=c("black", "red", "purple"), lty=c(1,1,1), lwd=2, bg=NULL)
# Compare
#par(mfrow=c(1,1))
plot(density(residual.boot.N[,1], bw=0.5), 
     main="Comparing bootstrap methods for parameter uncertainty: AGB", 
     lwd=2, lty=1)
lines(density(vector.boot[,1], bw=0.5), col='red', lwd=2, lty=1)
lines(density(sim.coef[,1], bw=0.5), col='purple', lwd=2, lty=1)
#legend('topright', legend=c("Residual Boot", "Pairs Boot", "Monte Carlo Normal"), 
#       col=c("black", "red", "purple"), lty=c(1,1,1), lwd=2, bg=NULL)

# Use boot library
# I think this is set up to perform a Pairs Bootstrap
BootstrapFunctionRegression <- function(data=dat1, index) {
  data <- data[index,] # we sample along rows of the data frame
  model.boot <- lm(AGB ~ backscatter, data=data)
  coef(model.boot)
}
bootstrappedModel <- boot(dat1, BootstrapFunctionRegression, R=10000)
bootstrappedModel

plot(bootstrappedModel, index=1)
boot.ci(bootstrappedModel, conf=0.95, type=c("basic", "bca", "perc"), index=1)
boot.ci(bootstrappedModel, conf=0.95, type=c("basic", "bca", "perc"), index=2)

plot(bootstrappedModel$t[,1], bootstrappedModel$t[,2],
     xlab="t1", ylab="t2", pch=16)

# Change it to run a Residuals bootstrap
BootstrapFunctionRegression <- function(data=dat1, index) {
  mod.object <- lm(AGB ~ backscatter, data=data)
  resids = mod.object$resid
  fittedValues = mod.object$fitted
  matr <- model.matrix(mod.object)
  # generating new values for each y[i], by adding bootstrapped resids to fitted values.
  Y <- fittedValues + resids[index] # we sample along rows of the data frame
  # Using model.matrix for the predictors 
  model.boot <- lm(Y ~ 0 + matr, data=data)
  coef(model.boot)
}
bootstrappedModel <- boot(dat1, BootstrapFunctionRegression, R=10000)
bootstrappedModel

plot(bootstrappedModel, index=1)
boot.ci(bootstrappedModel, conf=0.95, type=c("basic", "bca", "perc"), index=1)
boot.ci(bootstrappedModel, conf=0.95, type=c("basic", "bca", "perc"), index=2)

# ----
sqrt(var(linreg$residuals))
library(normwhn.test)
DHtest <- normality.test1(cbind(g0_AGB$AGB_ha, g0_AGB$`2018mean`))
normality.test1(cbind(g0_AGB$AGB_ha))

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

Zscores_backscatter <- scale(g0_AGB$`2018mean`)
hist(Zscores_backscatter)
Zscores_AGB <- scale(g0_AGB$AGB_ha)
hist(Zscores_AGB)

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

#-----------------------
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
mstems$stocking_10 <- calculateStocking(mstems$dbh_cm, min_diam = 10)
stocking_plot <- aggregate(stocking ~ plot_no, mstems, sum)
mplots <- merge(mplots, stocking_plot, by='plot_no', all=TRUE)
mplots$stocking <- mplots$stocking / mplots$area_ha

# Look at dominant species based on stocking densities. Add species_name (by_binomial)
dominant_species <- getDominantSpecies(mstems$sp_creole, mstems$plot_no)
dominant_species <- getDominantSpecies(mstems$sp_creole, mstems$plot_no, abundance = calculateBasalArea(mstems$dbh_cm))
dominant_species <- getDominantSpecies(mstems$sp_creole, mstems$plot_no, abundance = mstems$AGB_Chave14)

