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
# Load libraries -----
library("caret")
library("terra")
library("sf")
library("broom")
library("patchwork")
library("tidyverse")

# Set variables ----
year <- '2019'
code <- 'HV_nu'
suffix <- g0_variant <- 'med5'

# raw_dir <- file.path('data/raw/ALOS', year)
tidy_dir <- 'data/tidy'
palsar_dir <- file.path(tidy_dir, str_c('palsar_', year))
modeling_dir <- file.path('data/modeling', code)
field_agb_fp <- file.path(tidy_dir, 'survey_plots', 'plots_agb_noXtrms.rds')
field_agb_fp <- file.path(tidy_dir, 'survey_plots', 'plots_agb.rds')

# Set output filepaths ----
# List palsar mosaic files
(fps <- list.files(file.path(palsar_dir, 'mosaic_variants'), 
                  str_glue('{code}.*\\.tif$'), 
                  full.names = TRUE))
(g0_fp <- fps %>% nth(16))

# (dirs <- list.dirs(modeling_dir, recursive = FALSE))
# mod_dir <- dirs[[4]]
# (g0_fp <- list.files(mod_dir, 'tif$', full.names = TRUE))

# Get variant code
(g0_variant <- suffix <- str_extract(g0_fp, str_glue("(?<={code}_).*(?=\\.tif)")))
if(is.na(suffix)) {
  suffix <- ''
  g0_variant <- 'simple'
} else {
  suffix <- str_c('_', suffix)
}

# code <- 'HV_nu_noXtrms'
# modeling_dir <- file.path('data/modeling', code)
# g0_fp <- file.path(palsar_dir, str_glue("{code}{suffix}.tif"))

# Get path for model outputs
mod_dir <- file.path(modeling_dir, g0_variant, 'calibration')
dir.create(mod_dir, recursive = TRUE)

# extracted backscatter values
ex_shp <- file.path(mod_dir, str_glue('plots_agb_g0{suffix}.gpkg'))
ex_csv <- file.path(mod_dir, str_glue('plots_agb_g0{suffix}.csv'))

# Linear regression
ols_fp <- file.path(mod_dir, str_c("ols", suffix, ".rds"))
ols_vals_fp <- file.path(mod_dir, str_c('ols_results', suffix, '.csv'))
full_results_csv <- file.path(mod_dir, str_c('model_results', suffix, '.csv'))

# AGB map
agb_fp <- file.path(modeling_dir, g0_variant, str_c("agb_l0", suffix, ".tif"))


# Functions ----
#' Report OLS results
#' 
#' @param ols Object of class "lm"
#' @return Data frame with regression statistics for \code{ols}
#' @examples
#' vals <- report_ols_results(ols)
report_ols_results <- function(ols, df){
  
  # Regression Coefficients
  coefs <- broom::tidy(ols) %>% 
    cbind(confint(ols)) %>% 
    mutate(term = str_replace(term, '\\(Intercept\\)', 'int')) %>% 
    mutate(term = str_replace(term, 'backscatter', 'g0')) %>% 
    mutate(term = str_replace(term, 'as\\.vector\\(', '')) %>% 
    mutate(term = str_replace(term, '\\)', '')) %>% 
    pivot_wider(names_from=term, values_from=estimate:'97.5 %') %>% 
    cbind(broom::glance(ols))
  
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
  ANOVA <- anova(ols) %>% 
    broom::tidy() %>% 
    mutate(term = str_replace(term, 'backscatter', 'g0anova')) %>% 
    mutate(term = str_replace(term, 'as\\.vector\\(', '')) %>% 
    mutate(term = str_replace(term, '\\)', '')) %>% 
    mutate(term = str_replace(term, 'Residuals', 'resid')) %>% 
    pivot_wider(names_from=term, values_from=df:p.value)
  
  # Correlation coefficients
  sp <- cor.test(x=df$backscatter, y=df$AGB, method = 'spearman') %>% 
    broom::tidy() %>% 
    select(-c(alternative, method)) %>% 
    mutate(term='rhoSpear')
  pe <- cor.test(x=df$backscatter, y=df$AGB, method = 'pearson') %>% 
    broom::tidy() %>% 
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
  
  top %>% as_tibble(rownames='name')
}

#' Values to summarize in cross-validation 
#' 
#' @param data Matrix of observed (obs) and predicted (pred) values
#' @return Named vector with accuracy metrics for residuals
#' @examples
#' # Run cross-validation
#' set.seed(45)
#' cv_mod <- caret::train(AGB ~ backscatter, data = g0_agb, method = "lm",
#'                        trControl = caret::trainControl(method = "repeatedcv", 
#'                                                        number = 10, repeats = 10000,
#'                                                        summaryFunction = fxn.bias))
fxn.bias <- function(data, lev = NULL, model = NULL) {
  resids <- data$pred - data$obs
  rss <- sum(resids^2)
  n <- length(resids)
  df <- n-2
  mse <- rss / n
  c(r.squared = summary(lm(pred ~ obs, data))$r.squared,
    adj.r.squared = summary(lm(pred ~ obs, data))$adj.r.squared,
    Bias = abs(sum(resids) / n),
    RMSE = sqrt(mse),
    MAE = sum(abs(resids)) / n,
    MSE = mse,
    RSS = rss,
    MSS = rss/df,
    RSE = sqrt(rss / df))
}

# Get mean backscatter for each plot --------------------------------------------------------------
# Load raster and polygon data
g0 <- terra::rast(g0_fp)             # pre-processed in 02_process_ALOS_tiles.R 

# Add plot backscatter mean to polygons
plots_agb <- readRDS(field_agb_fp)

# Extract mean backscatter for each site polygon
ex <- terra::extract(g0, vect(plots_agb), mean, na.rm = TRUE)
names(ex) <- c('ID', 'g0_mean')

# Join mean to polygons
g0_agb <- bind_cols(plots_agb, ex)

# Save as polygons and CSV
g0_agb %>% st_write(ex_shp, append = FALSE)
g0_agb %>% st_drop_geometry() %>% write_csv(ex_csv)

# Linear regression ------------------------------------------------------------
# Extract just AGB and backscatter as dataframe
g0_agb <- read_csv(ex_csv) %>% 
  dplyr::select(AGB = AGB_ha, backscatter = g0_mean) %>% 
  mutate(AGB = as.numeric(AGB))

# Does correlation improve without plots 7 and 8, (probably responding to soil moisture)? 
# g0_agb <- g0_agb %>% filter(!(AGB==0 & backscatter>0.016))

# Basic OLS regression
ols <- lm(as.vector(AGB) ~ as.vector(backscatter), data=g0_agb, x=TRUE, y=TRUE)
ols %>% saveRDS(ols_fp)

# OLS results ---------------------------------------------------------------
# Convert training results to table
train_vals <- report_ols_results(ols, g0_agb)
train_vals %>% write_csv(ols_vals_fp)

# Plot error plots
ols_error_png <- file.path(mod_dir, str_c('ols_plots', suffix, '.png'))
png(ols_error_png)
opar <- par(mfrow = c(2,2), oma = c(0, 0, 1.1, 0))
plot(ols, las = 1)
dev.off()

# Repeated k-fold cross validation ---------------------------------------------
repeats <- 1000 # number of complete sets of folds to compute
cv_num <- 5 # number of folds
cv_fp <- file.path(mod_dir, str_c("crossval_", repeats, "x", cv_num, suffix, ".rds"))
if(!file.exists(cv_fp)) {
  
  # Run cross-validation
  set.seed(45)
  cv_mod <- caret::train(AGB ~ backscatter, 
                         data = g0_agb, 
                         method = "lm",
                         trControl = caret::trainControl(method = "repeatedcv", 
                                                         number = cv_num, 
                                                         repeats = repeats,
                                                         summaryFunction = fxn.bias))
  
  # Save model
  cv_mod %>% saveRDS(cv_fp)
}

# Load model
cv_mod <- readRDS(cv_fp)

# Format TEST values
cv_vals <- cv_mod$results[-1] %>% 
  tibble() %>% 
  pivot_longer(cols=everything())

# Format TRAIN values
train_vals <- read_csv(ols_vals_fp)

# Join tables
full_results <- rbind(mutate(train_vals, type = 'train'), 
              mutate(cv_vals, type = 'test'))

# Save
full_results %>% write_csv(full_results_csv)

# ~Look~ at values ----
# Get AGB, backscatter, plot number dataframe
g0_agb <- read_csv(ex_csv) %>%
  dplyr::select(plot_no, AGB = AGB_ha, backscatter = g0_mean) %>% 
  mutate(AGB = as.numeric(AGB))

# Scatterplot - AGB against backscatter ----
(p <- ggplot(g0_agb, aes(x=backscatter, y=AGB)) + geom_point() +
  labs(y = expression(paste("Aboveground biomass (Mg ha"^"-1", ")")),
       x = expression(paste("Radar backscatter, ",sigma['HV']^0," (m"^2, "/m"^2, ")"))) +
  geom_smooth(method='lm', se=TRUE, fullrange=TRUE, level=0.95, col='black', size=0.2) +
  # annotate(geom = 'text', label = plot_no, x = -Inf, y = +Inf, hjust = 0, vjust = 1) +
  # geom_text(aes(label = plot_no), hjust = -.2, vjust = 1) +
  theme_minimal())

ggsave(file.path(mod_dir, str_glue('scatter_agb_g0{suffix}.png')),
       p, 
       width=14, height=12.5, units='cm')

# Inset table with regression summary ----
full_results <- read_csv(full_results_csv)

d1 <- full_results %>% 
  filter(type == 'train') %>% 
  pivot_wider() %>% 
  transmute(
    Intercept = signif(estimate_int, 3),
    Slope = signif(estimate_g0, 3), 
    P = {if(p.value < 0.001) '< 0.001' else if(p.value < 0.05) '< 0.05' else '>= 0.05'},
    Equation = str_glue("AGB = {Intercept} + {Slope}x")) %>% 
  dplyr::select(Equation, P)

d2 <- full_results %>% 
  filter(type == 'test') %>% 
  pivot_wider() %>% 
  transmute(
    Adj_R2 = str_c(format(adj.r.squared, digits = 2), " ± ", 
                           format(adj.r.squaredSD, digits = 2)),
    RMSE = str_c(format(RMSE, digits = 3), " ± ", format(RMSESD, digits = 2)),
    MAE = str_c(format(MAE, digits = 3), " ± ", format(MAESD, digits = 2)),
    Bias = str_c(format(Bias, digits = 3), " ± ", format(BiasSD, digits = 2))) 

d3 <- bind_cols(d1, d2) %>%
  pivot_longer(everything()) %>% 
  column_to_rownames('name')

tab <- gridExtra::tableGrob(d3, 
                            cols = NULL,
                            theme = gridExtra::ttheme_minimal(base_size = 9,
                                                              base_colour = 'violetred4',
                                                              padding = unit(c(2,2), 'mm')))

# Overlay plot with regression table
p +
  inset_element(tab, 
                  left = 0, bottom = 0.8,
                  right = 0.2, top = 1, 
                  on_top = TRUE) +
  theme(plot.background = NULL)

# Save
ggsave(file.path('figures', year, 'modeling', str_glue('scatter_agb_g0_{g0_variant}_table.png')),
       width=14, height=12.5, units='cm')

p + 
  geom_text(aes(label = plot_no), hjust = -.2, vjust = 1, 
            size = 3) +
  inset_element(tab, 
                left = 0, bottom = 0.8,
                right = 0.2, top = 1, 
                on_top = TRUE) +
  theme(plot.background = NULL)

ggsave(file.path(mod_dir, str_glue('scatter_agb_g0_{g0_variant}_table.png')),
       width=14, height=12.5, units='cm')

# Create AGB raster --------------------------------------------------------------------------------
# Load data 
g0 <- rast(g0_fp); names(g0) <- 'backscatter'
ols <- readRDS(ols_fp)

# Apply linear regression model to create AGB map
agb.ras <- terra::predict(g0, ols, na.rm=TRUE)
minmax(agb.ras)
agb.ras %>% writeRaster(agb_fp, 
                        overwrite = TRUE,
                        wopt = list(datatype='INT2S', gdal='COMPRESS=LZW'))

# Plot map ----
# Load as terra rast 
agb <- rast(agb_fp)

# Aggregate
agb_a <- terra::aggregate(agb, fact = 7)
agb_a_df <- as.data.frame(agb_a, xy = T) 

# Plot
(agb_map <- ggplot(agb_a_df, aes(x = x, y = y, fill = backscatter)) +
  # geom_raster() +
    geom_tile() +
  # geom_sf(data = mex, fill = "transparent", size = 0.2, color = "gray70") +
  colormap::scale_fill_colormap(expression(paste("Aboveground\nbiomass (Mg ha"^"-1", ")")), 
                                na.value = "transparent", 
                                colormap = colormap::colormaps$viridis,
                                limits = c(0, 200), 
                                oob = scales::squish
                                ) +
    coord_fixed() +
    # theme_minimal() +
    ggthemes::theme_map() +
  theme(legend.position=c(0.23, 0.9),
        legend.title.align=0,
        legend.justification = c(1,1),
        legend.box.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_blank()) 
  )

# Save
map_fp <- file.path('figures/agb_maps', code, g0_variant, str_c("agb_l0", suffix, ".png"))
dir.create(dirname(map_fp), recursive = TRUE)

ggsave(map_fp, agb_map, 
       width=7,
       height=5.6,
       dpi=120)
