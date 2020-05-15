# Script where I experiment with cross-validation and bootstrapping to measure the regression fit
library(readr)
library(BIOMASS)
library(tidyverse)
require(boot)
require(MASS)

# Load data
g0.agb  <- readRDS("results/R_out/plots_g0agb_dfslim.rds")

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

# Experimentation, old ---- ###########################################################
# Load data - Desktop
# Load data - Mac
g0_AGB <- read_csv("~/GitHub/biomass-espanola/data/plots_g0nu_AGB.csv")
# just get the two columns we care about
g0.agb <- as.data.frame(cbind(g0_AGB$AGB_ha, g0_AGB$'2018mean')) %>% 
  rename(AGB = V1, backscatter =V2)

#----
# Scatterplot
p <- ggplot(g0.agb, aes(x=backscatter, y=AGB)) + geom_point() +
  labs(y = expression(paste("Aboveground biomass (MgC ha"^"-1", ")")), 
       x = expression(paste("Radar backscatter, ",sigma['HV']^0," (m"^2, "/m"^2, ")")))
p
p2 <- p + geom_smooth(method="lm", se=TRUE, fullrange=TRUE, level=0.95, col='black')
p2 + labs(caption = 'OLS regression')
p2
p1 + geom_smooth(method="rlm", col='red', se=TRUE, fullrange=TRUE, level=0.95) + 
  labs(caption = 'RLM regression')
p2 + geom_smooth(method="lm", col='red', fill='red', se=TRUE, fullrange=TRUE, level=0.95) + 
  geom_smooth(method="rlm", col='blue', fill='blue', se=TRUE, fullrange=TRUE, level=0.95)

# Manually construct confidence bands around OLS regression line
mm <- model.matrix(~ backscatter, data = g0.agb)
vars <- mm %*% vcov(ols) %*% t(mm)
sds <- sqrt(diag(vars))
t.val <- qt(1 - (1 - 0.95)/2, ols$df.residual)
t.val
g0.agb$LoCI.man <- ols$fitted.values - t.val * sds
g0.agb$HiCI.man <- ols$fitted.values + t.val * sds
p + geom_ribbon(aes(ymin=g0.agb$LoCI.man, ymax=g0.agb$HiCI.man), linetype=2, alpha=0.1)
ols.ci95 <- predict(ols, newdata = g0.agb, interval = 'confidence')
ols.pi95 <- predict(ols, newdata = g0.agb, interval = 'prediction')
p2 <- p + geom_ribbon(aes(ymin=ols.ci95[,2], ymax=ols.ci95[,3]), linetype=2, alpha=0.1)
p2 + geom_ribbon(aes(ymin=ols.pi95[,2], ymax=ols.pi95[,3]), linetype=2, alpha=0.1)



# plot BCa CI from bootstrapping - not sure if this is appropriate
g0.agb$loCI <- -6.96 + 799*g0.agb$backscatter
g0.agb$hiCI <- 9.63 + 1358*g0.agb$backscatter
p2 + geom_ribbon(aes(ymin=g0.agb$loCI, ymax=g0.agb$hiCI), linetype=2, alpha=0.1)

# plot parameter CI from OLS - not sure if this is appropriate
g0.agb$loCI <- -14.6454 + 713.18*g0.agb$backscatter
g0.agb$hiCI <- 14.64934 + 1349.83*g0.agb$backscatter
p2 + geom_ribbon(aes(ymin=g0.agb$loCI, ymax=g0.agb$hiCI), linetype=2, alpha=0.1)

# plot parameter CI from OLS - not sure if this is appropriate
g0.agb$loCI <- 14.64934 + 713.18*g0.agb$backscatter
g0.agb$hiCI <- -14.6454 + 1349.83*g0.agb$backscatter
p2 + geom_ribbon(aes(ymin=g0.agb$loCI, ymax=g0.agb$hiCI), linetype=2, alpha=0.1)

# 
boot.ols.100k$t0
confint(ols)[2,2] - ols$coefficients[2]
ols$model

# Bias
g0.agb$resids <- ols$residuals
p <- ggplot(g0.agb, aes(x=backscatter, y=resids)) + geom_point() +
  labs(y = expression(paste("Aboveground biomass (MgC ha"^"-1", ")")), 
       x = expression(paste("Radar backscatter, ",sigma['HV']^0," (m"^2, "/m"^2, ")")))
p
mean(abs(ols$residuals))
model.10000x5$results
model.10000x5


# Basic OLS regression
ols <- lm(g0.agb$AGB ~ g0.agb$backscatter, x=TRUE, y=TRUE)
summary(ols)
confint(ols)
cov2cor(vcov(ols))
anova(ols)
coef(ols)[2]*10000
coef(ols)[1]*100000000
opar <- par(mfrow = c(2,2), oma = c(0, 0, 1.1, 0))
plot(ols, las = 1)

ols$residuals
rmse <- sqrt(var(ols$residuals)) # RMSE
rmse <- sd(ols$residuals)
mse <- mean((residuals(ols))^2)
mse
rss <- sum(residuals(ols)^2)
rss
rse <- sqrt(rss / ols$df.residual)
rse
mean(abs(residuals(ols)))
summary(ols)$adj.r.squared
sigma(ols)

# OLS with intercept=0
ols.0 <- lm(g0.agb$AGB ~ 0+g0.agb$backscatter, x=TRUE, y=TRUE)
abline(ols.0, col='cyan')
summary(ols.0)
confint(ols.0)
cov2cor(vcov(ols))
anova(ols)
coef(ols)
opar <- par(mfrow = c(2,2), oma = c(0, 0, 1.1, 0))
plot(ols, las = 1)

#----
# Robust Linear Model regression, 
# with parameters set based on optimization performed by cross-validation below
rr <- rlm(AGB ~ backscatter, g0.agb, psi=psi.hampel)
rr$coefficients
abline(rr, col='red')
rr
anova(rr)
summary(rr)
opar <- par(mfrow = c(2,2), oma = c(0, 0, 1.1, 0))
plot(rr, las = 1)
View(rr$w)


rr.int0 <- rlm(AGB ~ 0 + backscatter, g0.agb, psi=psi.hampel)
rr.int0$coefficients
abline(rr, col='pink')
anova(rr.int0)
summary(rr.int0)
opar <- par(mfrow = c(2,2), oma = c(0, 0, 1.1, 0))
plot(rr.int0, las = 1)

View(rr.int0$w)

#----
# Cross-validation
set.seed(45)
# LOOCV
# Train the model and summarize results
model.loocv <- train(AGB ~ backscatter, data = g0.agb, method = "lm",
               trControl = trainControl(method = "LOOCV"))
print(model.loocv)
model.loocv$finalModel

# K-fold
# Train model and summarize results
model.5fold <- train(AGB ~ backscatter, data = g0.agb, method = "lm",
               trControl = trainControl(method = "cv", number = 5))
print(model.5fold)
model.5fold$finalModel

# Repeated K-fold
# Define training control, train model, and summarize results
# OLS w/ int=0, 10,000 x 10-fold
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
head(model.10000x10$resample)

set.seed(3)
model.2x5 <- train(AGB ~ backscatter, data = g0.agb, method = "lm",
                    trControl = trainControl(method = "repeatedcv", 
                                             number = 5, repeats = 2,
                                             summaryFunction = fxn.bias))
model.2x5$results
set.seed(3)
model.2x5.rmse <- train(AGB ~ backscatter, data = g0.agb, method = "lm",
                         trControl = trainControl(method = "repeatedcv", 
                                                  number = 5, repeats = 2))
model.2x5.rmse$results

model.10000x10 <- train(AGB ~ backscatter, data = g0.agb, method = "lm",
                       trControl = trainControl(method = "repeatedcv", 
                                                number = 10, repeats = 10000,
                                                returnResamp='all',
                                                trim=FALSE))
print(model.10000x10)
model.10000x10$results
model.10000x10$finalModel
head(model.10000x10$resample)
head(model.10000x10$metric)




# OLS w/ int=0, 10,000 x 10-fold
model.boot100k <- train(AGB ~ backscatter, data = g0.agb, method = "lm",
                             trControl = trainControl(method = "boot", 
                                                      number = 100000))
print(model.boot100k)
model.boot100k$results
model.boot100k$finalModel

# 5,000 x10-fold
model.5000x10.int0$results
model.5000x10.int0$finalModel

# Repeated K-fold with Robust Linear Regression
# RLM, 10,000 x 5-fold
model.5000x10.rlmI <- train(AGB ~ backscatter, data = g0.agb, method = "rlm", 
                          tuneGrid = expand.grid(intercept = TRUE, psi = 'psi.hampel'),
                          trControl = trainControl(method = "repeatedcv", 
                                                   number = 10, repeats = 5000))
print(model.5000x10.rlmI)
model.5000x10.rlmI$results
model.5000x10.rlmI$finalModel
abline(model.5000x10.rlmI$finalModel, col='red')


model.10000x10.rlmI
#----
# Pairs Bootstrap with the boot library
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

# OLS with int=0
boot.ols.int0.100k <- boot(g0.agb, function(data=g0.agb, index) {
    data <- data[index,] # we sample along rows of the data frame
    model.boot <- lm(AGB ~ 0 + backscatter, data=data)
    coef(model.boot)
  }, R=100000)
# Results
boot.ols.int0.100k
plot(boot.ols.int0.100k, index=1)
boot.ci(boot.ols.int0.100k, conf=0.95, type=c("basic", "bca", "perc"), index=1)
boot.ci(boot.ols.int0.100k, conf=0.95, type=c("basic", "bca", "perc"), index=2)

# Pairs bootstrap with RLM
boot.rlm.int0.100k <- boot(g0.agb, function(data=g0.agb, index) {
  data <- data[index,] # we sample along rows of the data frame
  model.boot <- rlm(AGB ~ 0 + backscatter, data=data, psi=psi.hampel)
  coef(model.boot)
}, R=100000)
# Results
boot.rlm.int0.100k
plot(boot.rlm.int0.100k, index=1)
boot.ci(boot.rlm.int0.100k, conf=0.95, type=c("basic", "bca", "perc"), index=1)
plot(boot.rlm.int0.100k$t[,1], boot.rlm.int0.100k$t[,2],
     xlab="t1", ylab="t2", pch=1)

# Pairs bootstrap with RLM w/ intercept
boot.rlm.100k <- boot(g0.agb, function(data=g0.agb, index) {
  data <- data[index,] # we sample along rows of the data frame
  model.boot <- rlm(AGB ~ backscatter, data=data, psi=psi.hampel)
  coef(model.boot)
}, R=100000)
# Results
boot.rlm.100k
plot(boot.rlm.100k, index=1)
boot.ci(boot.rlm.100k, conf=0.95, type=c("basic", "bca", "perc"), index=1)
boot.ci(boot.rlm.100k, conf=0.95, type=c("basic", "bca", "perc"), index=2)
plot(boot.rlm.100k$t[,1], boot.rlm.100k$t[,2],
     xlab="t1", ylab="t2", pch=1)

#----
# Save
save(boot.ols.100k, boot.ols.30k, boot.rlm.100k, boot.rlm.int0.100k, 
     model.10000x10, model.10000x10.rlm,  model.10000x10.rlmI, model.10000x10.int0, 
     model.10000x5, model.10000x5.boot, model.10000x5.rlm, model.10000x5.rlmI, 
     model.boot100k,
     file = "~/PROJECTS/Haiti_biomass/R_out/model_trains.RData")

#----
# Residuals bootstrap
BootstrapFunctionRegression <- function(data=g0.agb, index) {
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
bootstrappedModel <- boot(g0.agb, BootstrapFunctionRegression, R=10000)
bootstrappedModel

plot(bootstrappedModel, index=1)
boot.ci(bootstrappedModel, conf=0.95, type=c("basic", "bca", "perc"), index=1)
boot.ci(bootstrappedModel, conf=0.95, type=c("basic", "bca", "perc"), index=2)

# Bootstrap options - manual
# from "Using the non-parametric bootstrap for regression models in R" by Ian Dworkin
# Non-parametric bootstrap: Pairs (Random x) approach
N = 10000 # Perform N bootstrap iterations
BootstrapRandomX <- function(dat=g0.agb, mod.formula=formula(AGB ~ backscatter)){
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
plot(resid.model.1 ~ g0.agb$AGB)
plot(resid.model.1 ~ g0.agb$backscatter)

BootstrapFromResiduals <- function(mod.object = linreg, dat = g0.agb) {
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