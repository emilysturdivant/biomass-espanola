# *************************************************************************************************
# Script to:
#     * Calibrate regression model between AGB and backscatter; 
#     * Convert backscatter to AGB map
# Proceeds:
#     * calculate_AGB.R - calculates AGB by plot from field data
#     * process_ALOS_tiles.R
# Requires:
#     * backscatter image (g0)
#     * plot polygons with AGB by plot (plots_agb)
# *************************************************************************************************

# Load libraries
library(caret)
library(terra)
library(geobgu)
library(broom)
library(gdalUtils)
library(tmap)
require(graphics)
library(patchwork)
library(tidyverse)

year <- '2019'
g0_variant <- 'simple'
code <- 'sl_HV'

# raw_dir <- file.path('data/raw/ALOS', year)
tidy_dir <- 'data/tidy'
palsar_dir <- file.path(tidy_dir, str_c('palsar_', year))
masks_dir <- file.path(palsar_dir, 'masks')
landmask_fp <- file.path(masks_dir, 'hti_land_palsar.tif')
hti_poly_fp <- file.path(tidy_dir, "contextual_data/HTI_adm/HTI_adm0_fix.shp")
modeling_dir <- 'data/modeling'

hisp_bb <- st_bbox(c(xmin = -74.48133, ymax = 20.09044, 
                     xmax = -68.32267, ymin = 17.47022))
hti_bb <- st_bbox(c(xmin = -74.48133, ymax = 20.09044, 
                    xmax = -71.61815, ymin = 18.02180))

# results_dir <- 'data/results'
if(g0_variant == 'simple') {
  suffix <- ''
} else suffix <- g0_variant

list.files(palsar_dir, str_glue('{code}.*\\.tif$'))
g0_fp <- file.path(palsar_dir, str_glue("sl_HV{suffix}.tif"))

ex_shp <- file.path(modeling_dir, g0_variant, str_glue('plots_agb_g0{suffix}.gpkg'))
ex_csv <- file.path(modeling_dir, g0_variant, str_glue('plots_agb_g0{suffix}.csv'))

# Get mean backscatter for each plot --------------------------------------------------------------
# Load raster and polygon data
# Raster was pre-processed in 02_process_ALOS_tiles.R 
g0 <- rast(g0_fp)

# Add plot backscatter mean to polygons
field_agb_fp <- file.path(tidy_dir, 'survey_plots', 'plots_agb.rds')
plots_agb <- readRDS(field_agb_fp)

max(plots_agb$AGB_ha)

# Extract mean backscatter for each site polygon
ex <- terra::extract(g0, vect(plots_agb), mean, na.rm = TRUE)
names(ex) <- c('ID', 'g0_mean')

# Join mean to polygons
g0_agb <- bind_cols(plots_agb, ex)

# Save as polygons
g0_agb %>% st_write(ex_shp, delete_dsn = TRUE)

# Save as CSV
g0_agb %>% st_drop_geometry() %>% write_csv(ex_csv)

# Extract just AGB and backscatter as dataframe
# g0_agb <- st_read(ex_shp)
g0_agb <- st_drop_geometry(g0_agb) %>% 
  dplyr::select(AGB = AGB_ha, backscatter = g0_mean) %>% 
  mutate(AGB = as.numeric(AGB))

# ~Look~ at values -------------------------------------------------------------
# Scatterplot - AGB against backscatter ----
# (p <- ggplot(g0_agb, aes(x=backscatter, y=AGB)) + geom_point() +
#   labs(y = expression(paste("Aboveground biomass (Mg ha"^"-1", ")")), 
#        x = expression(paste("Radar backscatter, ",sigma['HV']^0," (m"^2, "/m"^2, ")"))) + 
#   geom_smooth(method='lm', se=TRUE, fullrange=TRUE, level=0.95, col='black', size=0.2) + 
#   # geom_rug() + 
#   theme_minimal())
# 
# ggsave(file.path('figures', year, 'modeling', str_glue('scatter_agb_g0{suffix}.png')),
#        width=14, height=12.5, units='cm')
# ggsave(file.path(modeling_dir, g0_variant, str_glue('scatter_agb_g0{suffix}.png')),
#        width=14, height=12.5, units='cm')

# Linear regression ---------------------------------------------------------------------------------
# See how correlation improves without plots 7 and 8, which are probably responding to soil moisture. 
# g0_agb <- g0_agb %>% filter(!(AGB==0 & backscatter>0.016))
g0_agb <- read_csv(ex_csv) %>% 
  dplyr::select(AGB = AGB_ha, backscatter = g0_mean) %>% 
  mutate(AGB = as.numeric(AGB))

# Basic OLS regression
ols_fp <- file.path(modeling_dir, g0_variant, str_c("ols", suffix, ".rds"))
ols <- lm(as.vector(AGB) ~ as.vector(backscatter), data=g0_agb, x=TRUE, y=TRUE)
ols %>% saveRDS(ols_fp)

# Report results -------------------------------------------------------------------------------------
#' Report OLS results
#' 
#' @param ols Object of class "lm"
#' @return Data frame with regression statistics for \code{ols}
#' @examples
#' vals <- report_ols_results(ols)
report_ols_results <- function(ols, df){
  
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
  sp <- cor.test(x=df$backscatter, y=df$AGB, method = 'spearman') %>% 
    tidy() %>% 
    select(-c(alternative, method)) %>% 
    mutate(term='rhoSpear')
  pe <- cor.test(x=df$backscatter, y=df$AGB, method = 'pearson') %>% 
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


ols_vals_fp <- file.path(modeling_dir, g0_variant, 
          str_c('ols_results', suffix, '.csv'))

# Convert training results to table
vals <- report_ols_results(ols, g0_agb)
vals %>% as_tibble(rownames = 'name') %>% 
  write_csv(ols_vals_fp)

# Plot error plots
png(file.path(modeling_dir, g0_variant, 
              str_c('ols_plots', suffix, '.png')))
opar <- par(mfrow = c(2,2), oma = c(0, 0, 1.1, 0))
plot(ols, las = 1)
dev.off()

# Repeated k-fold cross validation ---------------------------------------------
fxn.bias <- function(data, lev = NULL, model = NULL) {
  resids <- data$pred - data$obs
  # resids <- residuals(lm(pred ~ obs, data)) # produces very different results than the alternative above
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

repeats <- 1000
cv_num <- 5
cv_fp <- file.path(modeling_dir, g0_variant, 
                   str_c("crossval_", repeats, "x", cv_num, suffix, ".rds"))
if(!file.exists(cv_fp)) {
  set.seed(45)
  cv_mod <- train(AGB ~ backscatter, data = g0_agb, method = "lm",
                          trControl = trainControl(method = "repeatedcv", 
                                                   number = 10, repeats = 10000,
                                                   summaryFunction = fxn.bias))
  cv_mod %>% saveRDS(cv_fp)
}

# Load model
cv_mod <- readRDS(cv_fp)

# Format TEST values
cv_vals <- cv_mod$results[-1] %>% 
  tibble() %>% 
  pivot_longer(cols=everything()) %>% 
  mutate(type = 'test')

# Format TRAIN values
train_vals <- read_csv(ols_vals_fp) %>% 
  mutate(type = 'train')

train_vals <- vals[1] %>%
  as_tibble(rownames='name') %>%
  mutate(type = 'train')

# Join tables
mets <- rbind(train_vals, cv_vals)

# Save
full_results_csv <- file.path(modeling_dir, g0_variant, 
                              str_c('model_results', suffix, '.csv'))
mets %>% write_csv(full_results_csv)

# Create AGB raster --------------------------------------------------------------------------------
# Load data 
g0 <- rast(g0_fp); names(g0) <- 'backscatter'
ols <- readRDS(ols_fp)

# Apply linear regression model to create AGB map
agb_fp <- file.path(modeling_dir, g0_variant, str_c("agb_l0", suffix, ".tif"))
agb.ras <- terra::predict(g0, ols, na.rm=TRUE)
agb.ras %>% writeRaster(agb_fp)





# ~ OLD - Create 95% CI raster ---------------------------------------------------
# Function to create raster of confidence intervals
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
writeRaster(ras.agb.ci, file.path(results_dir, "agb/agb_CI_sub310.tif"), overwrite=TRUE)
