# Check correlation with stem number
# Load libraries -----
library("sf")
library("patchwork")
library("tidyverse")

# Set variables ----
year <- '2019'
code <- 'HV_nu'
g0_variant <- 'med5'

modeling_dir <- file.path('data/modeling', code)
mod_dir <- file.path(modeling_dir, g0_variant, 'calibration')
ex_csv <- file.path(mod_dir, str_glue('plots_agb_g0_{g0_variant}.csv'))

tidy_dir <- 'data/tidy'
mstems_fp <- file.path(tidy_dir, 'survey_plots/mstems_agb.rds')

# Load data ----
mstems <- readRDS(mstems_fp) 
cts <- mstems %>% 
  count(plot_no, wt = !is.na(dbh_cm), name = 'n_stems')

g0_agb <- read_csv(ex_csv)
g0_agb_stems <- g0_agb %>% 
  left_join(cts) %>% 
  mutate(stem_dens = n_stems / area)

# Plot ----
(p <- ggplot(g0_agb_stems, aes(x=g0_mean, y=AGB_ha, color=stem_dens)) + 
    geom_point() +
   geom_text(aes(label = plot_no, hjust = -.1, vjust = 1)) + 
    labs(y = expression(paste("Aboveground biomass (Mg ha"^"-1", ")")), 
         x = expression(paste("Radar backscatter, ",sigma['HV']^0," (m"^2, "/m"^2, ")"))) + 
    geom_smooth(method='lm', se=TRUE, fullrange=TRUE, level=0.95, col='black', size=0.2) + 
    geom_rug() + 
    theme_minimal())
ggsave(str_c("figures/qc_plots/03c_scatter_agb",g0_variant,"_g0_stemdens.png"), width=15, height=13, units='cm')
# One of the outliers has a high stem density (>1250/ha) and one of the highest backscatters, but a mediocre AGB.

(p <- ggplot(g0_agb_stems, aes(color=g0_mean, y=AGB_ha, x=stem_dens)) + 
    geom_point() +
    geom_text(aes(label = plot_no, hjust = -.1, vjust = 1)) + 
    labs(y = expression(paste("Aboveground biomass (Mg ha"^"-1", ")")), 
         x = expression(paste("Stems per hectare"))) + 
    # geom_smooth(method='lm', se=TRUE, fullrange=TRUE, level=0.95, col='black', size=0.2) + 
    geom_rug() + 
    theme_minimal())
ggsave(str_c("figures/qc_plots/03c_scatter_stemdens_agb",g0_variant,"_g0.png"), width=15, height=13, units='cm')
# Stem density up to about 300 stems/ha tracks well with AGB. 
# There aren't plots with stem density b/w 300 and ~600. 
# More than 600 stems/ha has a negative pattern with AGB.

(p <- ggplot(g0_agb_stems, aes(y=g0_mean, color=AGB_ha, x=stem_dens)) + 
    geom_point() +
    geom_text(aes(label = plot_no, hjust = -.1, vjust = 1)) + 
    labs(x = expression(paste("Stems per hectare")), 
         y = expression(paste("Radar backscatter, ",sigma['HV']^0," (m"^2, "/m"^2, ")"))) + 
    geom_smooth(method='lm', se=TRUE, fullrange=TRUE, level=0.95, col='black', size=0.2) + 
    geom_rug() + 
    theme_minimal())
ggsave(str_c("figures/qc_plots/03c_scatter_stemdens_g0_agb",g0_variant,".png"), width=15, height=13, units='cm')


# AGB with stems > 10 cm ----
fig_dir <- 'figures/qc_plots/stems_gt10cm'

# Remove stems dbh < 10 cm from calculations
mstems <- readRDS(mstems_fp) %>% 
  filter(dbh_cm >= 10 | is.na(dbh_cm))

plot_agb <- mstems %>% 
  agb_by_plot(plot_polys) %>% 
  st_drop_geometry()

# Histogram
summary(plot_agb$AGB_ha)
ct <- plot_agb %>% filter(!is.na(AGB_ha)) %>% nrow
binwidth = 5
p <-ggplot(plot_agb, aes(x=AGB_ha)) + 
  geom_histogram(binwidth = binwidth) +
  scale_y_continuous(name = "", breaks = seq(0, 8, 2)) +
  labs(caption = str_glue("N = {ct}, bin width = {binwidth}")) +
  scale_x_continuous(name = "Aboveground biomass (Mg/ha)") +
  theme_minimal() +
  theme(axis.title.y = element_blank())
p

cts <- mstems %>% 
  count(plot_no, wt = !is.na(dbh_cm), name = 'n_stems')

g0_agb_stems <- read_csv(ex_csv) %>% 
  select(plot_no, area, g0_mean) %>% 
  left_join(cts) %>% 
  mutate(stem_dens = n_stems / area)

g0_agb_stems <- g0_agb_stems %>% 
  left_join(plot_agb, by = 'plot_no')

(p <- ggplot(g0_agb_stems, aes(x=g0_mean, y=AGB_ha, color=stem_dens)) + 
   geom_point() +
   geom_text(aes(label = plot_no, hjust = -.1, vjust = 1)) + 
   labs(y = expression(paste("Aboveground biomass (Mg ha"^"-1", ")")), 
        x = expression(paste("Radar backscatter, ",sigma['HV']^0," (m"^2, "/m"^2, ")"))) + 
   geom_smooth(method='lm', se=TRUE, fullrange=TRUE, level=0.95, col='black', size=0.2) + 
   geom_rug() + 
   theme_minimal())

ggsave(file.path(fig_dir, str_c("03c_scatter_agb",g0_variant,"_g0_stemdens.png")), width=15, height=13, units='cm')
# One of the outliers has a high stem density (>1250/ha) and one of the highest backscatters, but a mediocre AGB.

(p <- ggplot(g0_agb_stems, aes(color=g0_mean, y=AGB_ha, x=stem_dens)) + 
    geom_point() +
    geom_text(aes(label = plot_no, hjust = -.1, vjust = 1)) + 
    labs(y = expression(paste("Aboveground biomass (Mg ha"^"-1", ")")), 
         x = expression(paste("Stems per hectare"))) + 
    # geom_smooth(method='lm', se=TRUE, fullrange=TRUE, level=0.95, col='black', size=0.2) + 
    geom_rug() + 
    theme_minimal())
ggsave(file.path(fig_dir, str_c("03c_scatter_stemdens_agb",g0_variant,"_g0.png")), 
       width=15, height=13, units='cm')
# Stem density up to about 300 stems/ha tracks well with AGB. 
# There aren't plots with stem density b/w 300 and ~600. 
# More than 600 stems/ha has a negative pattern with AGB.

(p <- ggplot(g0_agb_stems, aes(y=g0_mean, color=AGB_ha, x=stem_dens)) + 
    geom_point() +
    geom_text(aes(label = plot_no, hjust = -.1, vjust = 1)) + 
    labs(x = expression(paste("Stems per hectare")), 
         y = expression(paste("Radar backscatter, ",sigma['HV']^0," (m"^2, "/m"^2, ")"))) + 
    geom_smooth(method='lm', se=TRUE, fullrange=TRUE, level=0.95, col='black', size=0.2) + 
    geom_rug() + 
    theme_minimal())
ggsave(file.path(fig_dir, str_c("03c_scatter_stemdens_g0_agb",g0_variant,".png")), 
       width=15, height=13, units='cm')

# Basic OLS regression
ols <- lm(as.vector(stem_dens) ~ as.vector(g0_mean), 
          data=g0_agb_stems, x=TRUE, y=TRUE)
broom::glance(ols)
ols <- lm(as.vector(AGB_ha) ~ as.vector(g0_mean), 
          data=g0_agb_stems, x=TRUE, y=TRUE)
broom::glance(ols)


# AGB with WD==1 ----
fig_dir <- 'figures/qc_plots/wd_pt7'

# Remove stems dbh < 10 cm from calculations
mstems <- readRDS(mstems_fp) %>% 
  filter(dbh_cm >= 5 | is.na(dbh_cm))

plot_agb <- mstems %>% 
  mutate(meanWD = 0.7) %>% 
  agb_by_plot(plot_polys) %>% 
  st_drop_geometry()

# Histogram
summary(plot_agb$AGB_ha)
ct <- plot_agb %>% filter(!is.na(AGB_ha)) %>% nrow
binwidth = 5
p <-ggplot(plot_agb, aes(x=AGB_ha)) + 
  geom_histogram(binwidth = binwidth) +
  scale_y_continuous(name = "", breaks = seq(0, 8, 2)) +
  labs(caption = str_glue("N = {ct}, bin width = {binwidth}")) +
  scale_x_continuous(name = "Aboveground biomass (Mg/ha)") +
  theme_minimal() +
  theme(axis.title.y = element_blank())
p

cts <- mstems %>% 
  count(plot_no, wt = !is.na(dbh_cm), name = 'n_stems')

g0_agb_stems <- read_csv(ex_csv) %>% 
  select(plot_no, area, g0_mean) %>% 
  left_join(cts) %>% 
  mutate(stem_dens = n_stems / area)

g0_agb_stems <- g0_agb_stems %>% 
  left_join(plot_agb, by = 'plot_no')

# Basic OLS regression
ols <- lm(as.vector(stem_dens) ~ as.vector(g0_mean), 
          data=g0_agb_stems, x=TRUE, y=TRUE)
broom::glance(ols)
ols <- lm(as.vector(AGB_ha) ~ as.vector(g0_mean), 
          data=g0_agb_stems, x=TRUE, y=TRUE)
broom::glance(ols)

# Plots
(p <- ggplot(g0_agb_stems, aes(x=g0_mean, y=AGB_ha, color=stem_dens)) + 
    geom_point() +
    geom_text(aes(label = plot_no, hjust = -.1, vjust = 1)) + 
    labs(y = expression(paste("Aboveground biomass (Mg ha"^"-1", ")")), 
         x = expression(paste("Radar backscatter, ",sigma['HV']^0," (m"^2, "/m"^2, ")"))) + 
    geom_smooth(method='lm', se=TRUE, fullrange=TRUE, level=0.95, col='black', size=0.2) + 
    geom_rug() + 
    theme_minimal())

ggsave(file.path(fig_dir, str_c("03c_scatter_agb",g0_variant,"_g0_stemdens.png")), width=15, height=13, units='cm')
# One of the outliers has a high stem density (>1250/ha) and one of the highest backscatters, but a mediocre AGB.

(p <- ggplot(g0_agb_stems, aes(color=g0_mean, y=AGB_ha, x=stem_dens)) + 
    geom_point() +
    geom_text(aes(label = plot_no, hjust = -.1, vjust = 1)) + 
    labs(y = expression(paste("Aboveground biomass (Mg ha"^"-1", ")")), 
         x = expression(paste("Stems per hectare"))) + 
    # geom_smooth(method='lm', se=TRUE, fullrange=TRUE, level=0.95, col='black', size=0.2) + 
    geom_rug() + 
    theme_minimal())
ggsave(file.path(fig_dir, str_c("03c_scatter_stemdens_agb",g0_variant,"_g0.png")), 
       width=15, height=13, units='cm')
# Stem density up to about 300 stems/ha tracks well with AGB. 
# There aren't plots with stem density b/w 300 and ~600. 
# More than 600 stems/ha has a negative pattern with AGB.

(p <- ggplot(g0_agb_stems, aes(y=g0_mean, color=AGB_ha, x=stem_dens)) + 
    geom_point() +
    geom_text(aes(label = plot_no, hjust = -.1, vjust = 1)) + 
    labs(x = expression(paste("Stems per hectare")), 
         y = expression(paste("Radar backscatter, ",sigma['HV']^0," (m"^2, "/m"^2, ")"))) + 
    geom_smooth(method='lm', se=TRUE, fullrange=TRUE, level=0.95, col='black', size=0.2) + 
    geom_rug() + 
    theme_minimal())
ggsave(file.path(fig_dir, str_c("03c_scatter_stemdens_g0_agb",g0_variant,".png")), 
       width=15, height=13, units='cm')

# AGB without giant Ficus ----
fig_dir <- 'figures/qc_plots/wo_xtrmFicus'

# Remove stems dbh < 10 cm from calculations
mstems <- readRDS(mstems_fp) %>% 
  filter((dbh_cm >= 5 & dbh_cm < 300) | is.na(dbh_cm))

plot_agb <- mstems %>% 
  agb_by_plot(plot_polys) %>% 
  st_drop_geometry()

# Histogram
summary(plot_agb$AGB_ha)
ct <- plot_agb %>% filter(!is.na(AGB_ha)) %>% nrow
binwidth = 5
p <-ggplot(plot_agb, aes(x=AGB_ha)) + 
  geom_histogram(binwidth = binwidth) +
  scale_y_continuous(name = "", breaks = seq(0, 8, 2)) +
  labs(caption = str_glue("N = {ct}, bin width = {binwidth}")) +
  scale_x_continuous(name = "Aboveground biomass (Mg/ha)") +
  theme_minimal() +
  theme(axis.title.y = element_blank())
p

cts <- mstems %>% 
  count(plot_no, wt = !is.na(dbh_cm), name = 'n_stems')

g0_agb_stems <- read_csv(ex_csv) %>% 
  select(plot_no, area, g0_mean) %>% 
  left_join(cts) %>% 
  mutate(stem_dens = n_stems / area)

g0_agb_stems <- g0_agb_stems %>% 
  left_join(plot_agb, by = 'plot_no')

# Basic OLS regression
ols <- lm(as.vector(AGB_ha) ~ as.vector(g0_mean), 
          data=g0_agb_stems, x=TRUE, y=TRUE)
broom::glance(ols)

# Plots
(p <- ggplot(g0_agb_stems, aes(x=g0_mean, y=AGB_ha, color=stem_dens)) + 
    geom_point() +
    geom_text(aes(label = plot_no, hjust = -.1, vjust = 1)) + 
    labs(y = expression(paste("Aboveground biomass (Mg ha"^"-1", ")")), 
         x = expression(paste("Radar backscatter, ",sigma['HV']^0," (m"^2, "/m"^2, ")"))) + 
    geom_smooth(method='lm', se=TRUE, fullrange=TRUE, level=0.95, col='black', size=0.2) + 
    geom_rug() + 
    theme_minimal())

ggsave(file.path(fig_dir, str_c("03c_scatter_agb",g0_variant,"_g0_stemdens.png")), width=15, height=13, units='cm')
# One of the outliers has a high stem density (>1250/ha) and one of the highest backscatters, but a mediocre AGB.

(p <- ggplot(g0_agb_stems, aes(color=g0_mean, y=AGB_ha, x=stem_dens)) + 
    geom_point() +
    geom_text(aes(label = plot_no, hjust = -.1, vjust = 1)) + 
    labs(y = expression(paste("Aboveground biomass (Mg ha"^"-1", ")")), 
         x = expression(paste("Stems per hectare"))) + 
    # geom_smooth(method='lm', se=TRUE, fullrange=TRUE, level=0.95, col='black', size=0.2) + 
    geom_rug() + 
    theme_minimal())
ggsave(file.path(fig_dir, str_c("03c_scatter_stemdens_agb",g0_variant,"_g0.png")), 
       width=15, height=13, units='cm')
# Stem density up to about 300 stems/ha tracks well with AGB. 
# There aren't plots with stem density b/w 300 and ~600. 
# More than 600 stems/ha has a negative pattern with AGB.

(p <- ggplot(g0_agb_stems, aes(y=g0_mean, color=AGB_ha, x=stem_dens)) + 
    geom_point() +
    geom_text(aes(label = plot_no, hjust = -.1, vjust = 1)) + 
    labs(x = expression(paste("Stems per hectare")), 
         y = expression(paste("Radar backscatter, ",sigma['HV']^0," (m"^2, "/m"^2, ")"))) + 
    geom_smooth(method='lm', se=TRUE, fullrange=TRUE, level=0.95, col='black', size=0.2) + 
    geom_rug() + 
    theme_minimal())
ggsave(file.path(fig_dir, str_c("03c_scatter_stemdens_g0_agb",g0_variant,".png")), 
       width=15, height=13, units='cm')


# g0 converted back to sigma0 ----
# Convert g0 back to sigma0 from linear power after getting arithmetic mean
fig_dir <- 'figures/qc_plots/03_sigma0'

# Get stems values
# Remove stems dbh < 10 cm from calculations
mstems <- readRDS(mstems_fp) %>% 
  filter((dbh_cm >= 5 & dbh_cm < 300) | is.na(dbh_cm))

plot_agb <- mstems %>% 
  agb_by_plot(plot_polys) %>% 
  st_drop_geometry()

cts <- mstems %>% 
  count(plot_no, wt = !is.na(dbh_cm), name = 'n_stems')

g0 <- 10 ^ (g0 / 10)
log10(g0_mean*10)

# Join to plot AGB data
g0_agb_stems <- read_csv(ex_csv) %>% 
  select(plot_no, area, g0_mean) %>% 
  mutate(g0_sig0a = log10(g0_mean*10),
         g0_sig0b = log10(g0_mean)*10) %>% 
  left_join(cts) %>% 
  mutate(stem_dens = n_stems / area)

g0_agb_stems <- g0_agb_stems %>% 
  left_join(plot_agb, by = 'plot_no')

# Basic OLS regression
ols <- lm(as.vector(AGB_ha) ~ as.vector(g0_sig0b), 
          data=g0_agb_stems, x=TRUE, y=TRUE)
broom::glance(ols)
ols <- lm(as.vector(AGB_ha) ~ as.vector(g0_sig0a), 
          data=g0_agb_stems, x=TRUE, y=TRUE)
broom::glance(ols)

# Plots
(p <- ggplot(g0_agb_stems, aes(x=g0_sig0b, y=AGB_ha, color=stem_dens)) + 
    geom_point() +
    geom_text(aes(label = plot_no, hjust = -.1, vjust = 1)) + 
    labs(y = expression(paste("Aboveground biomass (Mg ha"^"-1", ")")), 
         x = expression(paste("Radar backscatter, ",sigma['HV']^0," (m"^2, "/m"^2, ")"))) + 
    geom_smooth(method='lm', se=TRUE, fullrange=TRUE, level=0.95, col='black', size=0.2) + 
    geom_rug() + 
    theme_minimal())

ggsave(file.path(fig_dir, str_c("03c_scatter_agb",g0_variant,"_g0_stemdens.png")), width=15, height=13, units='cm')
# One of the outliers has a high stem density (>1250/ha) and one of the highest backscatters, but a mediocre AGB.

(p <- ggplot(g0_agb_stems, aes(color=g0_sig0b, y=AGB_ha, x=stem_dens)) + 
    geom_point() +
    geom_text(aes(label = plot_no, hjust = -.1, vjust = 1)) + 
    labs(y = expression(paste("Aboveground biomass (Mg ha"^"-1", ")")), 
         x = expression(paste("Stems per hectare"))) + 
    # geom_smooth(method='lm', se=TRUE, fullrange=TRUE, level=0.95, col='black', size=0.2) + 
    geom_rug() + 
    theme_minimal())
ggsave(file.path(fig_dir, str_c("03c_scatter_stemdens_agb",g0_variant,"_g0.png")), 
       width=15, height=13, units='cm')
# Stem density up to about 300 stems/ha tracks well with AGB. 
# There aren't plots with stem density b/w 300 and ~600. 
# More than 600 stems/ha has a negative pattern with AGB.

(p <- ggplot(g0_agb_stems, aes(y=g0_sig0b, color=AGB_ha, x=stem_dens)) + 
    geom_point() +
    geom_text(aes(label = plot_no, hjust = -.1, vjust = 1)) + 
    labs(x = expression(paste("Stems per hectare")), 
         y = expression(paste("Radar backscatter, ",sigma['HV']^0," (m"^2, "/m"^2, ")"))) + 
    geom_smooth(method='lm', se=TRUE, fullrange=TRUE, level=0.95, col='black', size=0.2) + 
    geom_rug() + 
    theme_minimal())
ggsave(file.path(fig_dir, str_c("03c_scatter_stemdens_g0_agb",g0_variant,".png")), 
       width=15, height=13, units='cm')

