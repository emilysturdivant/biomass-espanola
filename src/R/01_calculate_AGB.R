# ------------------------------------------------------------------------------
# Script to:
#     * Calculate plot AGB from field inventory data
# Proceeds:
#     * python scripts to preprocess field data
# Preceeds:
#     * process_ALOS_tiles.R - merges tiles and makes masks
#     * regression_AGB-g0.R - creates AGB map
# Requires:
#     * plot outline polygons 
#     * field data table with plot ID, tree name, DBH, H
#     * lookup table to match tree names with wood densities
# ------------------------------------------------------------------------------

# Load libraries ----
library(readr)
library(BIOMASS)
library(tidyverse)
library(sf)

# Filenames ----
results_dir <- 'data/results'
tidy_dir <- 'data/tidy'

stems_fp_in <- "data/species_and_wds/haiti_data_wds2.csv"
plots_fp_in <- "data/species_and_wds/mplots_geoms.csv"
plots_shp <- file.path(tidy_dir, 'survey_plots', 'all_plots.shp')

# Functions ----
# Load and standardize polygons
standardize.names <- function(fp) {
  try(
    {
      pol <- fp %>% 
        # Read shapefile
        sf::st_read(quiet=TRUE) %>% 
        # Remove Z dimension from those that have it
        st_zm() %>% 
        # Transform to lat long
        st_transform(4326)
      
      pol$area <- st_area(pol) %>% units::set_units(value = ha) 
      pol$Name <- basename(fp)
      
      pol$biomass_tf <- (if (str_detect(pol$Name, 'No_?biomass')) 0 else 1)
      
      pol <- pol[c('Name', 'area', 'biomass_tf')]
    }, 
    
    silent=FALSE
  )
}

# HEIGHTS, input: mstems$dbh_cm, mstems$ht_m, output: mstems$H, mstems$Hrse
interp_heights <- function(.data, method = 'log1'){
  # Create H-D model and use to fill missing heights
  
  # Create model
  HDmodel <- BIOMASS::modelHD(
    D = .data$dbh_cm,
    H = .data$ht_m,
    method = method, # best model is log1 based on prior checking
    useWeight = TRUE,
    drawGraph = FALSE
  )
  
  # Retrieve height data
  .data$H <- .data$ht_m
  .data$Hrse <- 0.5 # recommended RSE for observed heights
  filt <- is.na(.data$H)
  
  # Model missing heights 
  .data$H[filt] <- BIOMASS::retrieveH(D = .data$dbh_cm, model = HDmodel)$H[filt]
  .data$Hrse[filt] <- HDmodel$RSE
  
  # Return
  return(.data)
}

agb_by_plot <- function(mstems, plot_polys) {
  # compute AGB(Mg) per plot
  AGBplot <- BIOMASS::summaryByPlot(
    computeAGB(
      D = mstems$dbh_cm,
      WD = mstems$meanWD,
      H = mstems$H
    ), 
    mstems$plot_no
  )
  
  # Calculate AGB per hectare
  plots_agb <- #merge(plot_polys, AGBplot, by.x='plot_no', by.y='plot', all=TRUE) %>% 
    left_join(plot_polys, AGBplot, by = c(plot_no = 'plot')) %>% 
    mutate(AGB_ha = AGB / area,
           AGB_ha = replace_na(AGB_ha, 0)) 
  
  return(plots_agb)
}

# Load data ----
mstems <- read_csv(stems_fp_in)
mplots <- read_csv(plots_fp_in, col_types = cols(plot_no = col_integer()))

# Prep mstems ----
mstems <- mstems %>% 
  filter(is.na(dbh_cm) | (dbh_cm >= 5 & dbh_cm < 300))

# Consolidate plot polygons ---- ######################################################
# Extract plot values from stems
plots_no <- mstems %>% 
  group_by(plot_no, plot_shp) %>% 
  summarize()

# Revise shp names to match raw data
plots_no$plot_shp <- plots_no$plot_shp %>% 
  str_replace('Pandiassou2', '2pandiassou') %>% 
  str_replace('Campeche', 'Camapeche')

# Load and merge plot polygons into one DF
fps <- list.files(path="data/raw/survey_plot_shps", pattern="\\.shp$", full.names=TRUE)
plots <- fps %>% lapply(standardize.names)
plots <- do.call(rbind, plots)

# Join to plot_no from mstems
plot_polys <- plots %>% 
  full_join(plots_no, by=c('Name' ='plot_shp')) %>% 
  filter(!is.na(plot_no)) 

# Save
plot_polys %>% 
  st_write(plots_shp, append=FALSE)

# Load
plot_polys <- st_read(plots_shp)

# Compute AGB ---- 
mstems <- mstems %>% interp_heights('log1')

# AGB, input: mstems[c(dbh_cm, meanWD, H, plot_no)]
mstems$agb <- computeAGB(
  D = mstems$dbh_cm,
  WD = mstems$meanWD,
  H = mstems$H
)

# Save
mstems_fp <- file.path(tidy_dir, 'survey_plots/mstems_agb_alldbh.rds')
mstems_fp <- file.path(tidy_dir, 'survey_plots/mstems_agb.rds')
mstems_fp <- file.path(tidy_dir, 'survey_plots/mstems_agb_noXtrms.rds')
saveRDS(mstems, mstems_fp)
mstems <- readRDS(mstems_fp)

# Compute AGB per plot ---- 
mstems <- mstems %>% 
  filter(is.na(dbh_cm) | (dbh_cm >= 5 & dbh_cm < 300))

# compute AGB(Mg) per plot
plots_agb <- agb_by_plot(mstems, plot_polys)

plots_agb_fp <- file.path(tidy_dir, 'survey_plots', 'plots_agb.rds')
plots_agb_fp <- file.path(tidy_dir, 'survey_plots', 'plots_agb_noXtrms.rds')
saveRDS(plots_agb, plots_agb_fp)

plots_agb <- readRDS(plots_agb_fp)


# # Look at data 
# summary(mstems$dbh_cm, na.rm=TRUE)
# sd(mstems$dbh_cm, na.rm=TRUE)
# n_stems <- nrow(filter(mstems, !is.na(dbh_cm)))
# 
# # 
# (p <-ggplot(mstems, aes(x=dbh_cm)) + 
#   geom_histogram(binwidth=5) +
#   labs(y = "")+
#   scale_x_continuous(name = "Diameter at breast height (cm)") +
#     ggtitle(str_c("Histogram of tree diameters (N = ",
#                   format(n_stems, big.mark = ','), ", bin width = 5 m)")) +
#   theme_minimal())
# 
# 
# dbh_lim <- 140
# filter(mstems, dbh_cm > dbh_lim)
# mstems %>% 
#   mutate(dbh_new = ifelse(dbh_cm > dbh_lim, dbh_lim, dbh_cm)) %>% 
#   ggplot(aes(dbh_new)) +
#   geom_histogram(binwidth = 2) +
#   labs(y = "")+
#   scale_x_continuous(name = "Diameter at breast height (cm)",
#                      breaks = seq(0, dbh_lim, 10),
#                      limits=c(0, dbh_lim)) +
#   ggtitle(str_c("Histogram of tree diameters (N = ",
#                 format(n_stems, big.mark = ','), 
#                 ", bin width = 5 m, outliers grouped at ",
#                 dbh_lim, " cm)")
#           ) +
#   theme_minimal()
#   
# n_hts <- nrow(filter(mstems, !is.na(ht_m)))
# (p <-ggplot(mstems, aes(x=ht_m)) + 
#   geom_histogram(binwidth=1) +
#   labs(x = expression(paste("Height (m)")), 
#        y = "")+
#   ggtitle(str_c("Histogram of tree heights (N = ",
#                 format(n_hts, big.mark = ','), 
#                 ", bin width = 1 m)"))+
#   theme_minimal())
# 
# # Boxplots per plot 
# # DBH
# (p <-ggplot(mstems, aes(y=dbh_cm, x=plot_no, group=plot_no)) + 
#    geom_boxplot() +
#    labs(y = expression(paste("DBH (cm)")), 
#         x = "")+
#    ggtitle(str_c("Tree diameter (DBH; N = ",
#            format(n_stems, big.mark = ','), 
#            ")"))+
#    theme_minimal())
# ggsave('figures/qc_plots/boxes_byplot_dbh.png', width=6, height=2.5)
# 
# # Height
# (p <-ggplot(mstems, aes(y=ht_m, x=plot_no, group=plot_no)) + 
#     geom_boxplot() +
#     labs(y = expression(paste("Height (m)")), 
#          x = "")+
#     ggtitle(str_c("Tree heights (N = ",
#             format(n_hts, big.mark = ','), 
#             ")")) +
#     theme_minimal())
# ggsave('figures/qc_plots/boxes_byplot_ht.png', width=6, height=2.5)
# 
# # DBH
# (p <-ggplot(mstems, aes(x=dbh_cm, y=plot_no, group=plot_no)) + 
#     geom_boxplot() +
#     labs(x = expression(paste("DBH (cm)")), 
#          y = "")+
#     ggtitle(str_c("Tree diameter \n(DBH; N = ",
#             format(n_stems, big.mark = ','), 
#             ")")) +
#     theme_minimal()+
#     theme(axis.text.y=element_blank()))
# ggsave('figures/qc_plots/boxes_byplot_dbhX.png', width=2.5, height=6)
# 
# # Height
# (p <-ggplot(mstems, aes(x=ht_m, y=plot_no, group=plot_no)) + 
#     geom_boxplot() +
#     labs(x = expression(paste("Height (m)")), 
#          y = "")+
#     ggtitle(str_c("Tree heights (N = ",
#                   format(n_hts, big.mark = ','), 
#                   ")")) +
#     theme_minimal()+
#     theme(axis.text.y=element_blank()))
# ggsave('figures/qc_plots/boxes_byplot_htX.png', width=2.5, height=6)
# 
# mstems_temp <- mstems %>% group_by(plot_no) %>% 
#   summarize(mean_ht = median(ht_m, na.rm=T)) %>% 
#   arrange(desc(mean_ht)) %>% 
#   left_join(mstems) %>% 
#   mutate(plot_no = as.factor(plot_no))
# (p <-ggplot(mstems_temp, aes(y=ht_m, x=plot_no, group=plot_no)) + 
#     geom_boxplot() +
#     labs(y = expression(paste("Height (m)")), 
#          x = "")+
#     ggtitle(str_c("Tree heights (N = ",
#                   format(n_hts, big.mark = ','), 
#                   ")")) +
#     theme_minimal()+
#     theme(axis.text.x=element_blank()))

# # Look at data ---- 
# summary(mstems$meanWD)
# wd_sd=sd(mstems$meanWD, na.rm=TRUE)
# summary(mstems$agb)
# sd(mstems$agb, na.rm=TRUE)
# 
# # Plot histograms and density plot
# mstems.filt <- mstems %>% filter(agb < 10)
# (p <-ggplot(mstems.filt, aes(x=agb)) + 
#   geom_histogram(binwidth=0.05) +
#   labs(y = "")+
#   scale_x_continuous(name = "AGB") +
#   ggtitle("Histogram of AGB (bin width = 0.05 m)")+
#   theme_minimal())
# limO <- 2
# mstems.filt %>% 
#   mutate(x_new = ifelse(agb > limO, limO, agb)) %>% 
#   ggplot(aes(x_new)) +
#   geom_histogram(binwidth = 0.002) +
#   labs(y = "")+
#   scale_x_continuous(name = "AGB",
#                      breaks = seq(0, 5, 0.1),
#                      limits=c(0, limO)) +
#   ggtitle(str_c("Histogram of AGB (bin width = 0.05 m, outliers grouped at ",limO,")"))+
#   theme_minimal()
# (p.dens.agb <- ggplot(mstems.filt, aes(x=agb)) +
#   geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)+ 
#   labs(x = expression(paste("Aboveground biomass (MgC)")), 
#        y = "Density"))
# ggplot(mstems.filt, aes(x=plot_no, y=agb))+
#   geom_boxplot()
# 
# # Look at per plot 
# plot.means <- mstems %>% 
#   group_by(plot_no) %>%
#   summarise_at(vars(dbh_cm, H, meanWD, agb, plot_area), list(mean = mean))
# summary(plot.means$agb_mean)
# sd(plot.means$meanWD_mean, na.rm=TRUE)
# plot.sums <- mstems %>% 
#   group_by(plot_no) %>%
#   summarise_at(vars(dbh_cm, H, meanWD, agb), list(sum = sum)) 
# plot.sums <- merge(mplots, plot.sums, by='plot_no', all=TRUE)
# plot.dens <- plot.sums %>% mutate_at(vars(dbh_cm_sum, H_sum, agb_sum), `/`, y = .$plot_area)
# summary(plot.dens$meanWD_sum)
# sd(plot.dens$dbh_cm_sum, na.rm=TRUE)
# 
# 
# # Look at data ----
# summary(plots_agb$AGB_ha)
# # Plots
# ggplot(plots_agb, aes(x=plot_no, y=AGB_ha))+
#   geom_boxplot()
# p.dens.agb <- ggplot(plots_agb, aes(x=AGB_ha)) +
#   geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)+ 
#   labs(x = expression(paste("Aboveground biomass (MgC)")), 
#        y = "Density")
# p.dens.agb
# p <-ggplot(plots_agb, aes(x=AGB_ha)) + 
#   geom_histogram(binwidth=5) +
#   scale_y_continuous(name = "", breaks = seq(0, 8, 2)) +
#   labs(y = "")+
#   scale_x_continuous(name = "Aboveground Biomass (Mg/ha)") +
#   ggtitle("Histogram of AGB (N = 36, bin width = 5)")+
#   theme_minimal()
# p
# filt <- plots_agb %>% filter(AGB_ha > 0)
# summary(filt$AGB_ha)
# 
# 
# (summaries <- as_tibble(cbind(Mean = c(dbh=mean(mstems$dbh_cm, na.rm=TRUE),
#                               ht=mean(mstems$H, na.rm=TRUE),
#                               wd=mean(mstems$meanWD, na.rm=TRUE),
#                               sdWD=mean(mstems$sdWD, na.rm=TRUE),
#                               agb=mean(mstems$agb, na.rm=TRUE)), 
#                     SD = c(dbh=sd(mstems$dbh_cm, na.rm=TRUE), 
#                              ht=sd(mstems$H, na.rm=TRUE),
#                              wd=sd(mstems$meanWD, na.rm=TRUE),
#                              sdWD=NA, 
#                              agb=sd(mstems$agb, na.rm=TRUE))), 
#                     rownames = 'variable'))
# 
# summaries %>% write_csv('data/reports/01_field_plots.csv')
# 
# # Look at values
# plots_agb$AGB <- as.vector(plots_agb$AGB_ha)
# plots_agb$area <- as.vector(plots_agb$area)
# g0.agb <- plots_agb
# AGB <- c(
#   all_mean=mean(g0.agb$AGB),
#   all_sd=sd(g0.agb$AGB),
#   biomass_mean=mean(g0.agb[9:36, ]$AGB),
#   biomass_sd=sd(g0.agb[9:36, ]$AGB),
#   min=min(g0.agb$AGB),
#   max=max(g0.agb$AGB))
# g0_AGB <- plots_agb
# Area <- c(
#   all_mean=mean(g0_AGB$area),
#   all_sd=sd(g0_AGB$area),
#   biomass_mean=mean(g0_AGB[9:36, ]$area),
#   biomass_sd=sd(g0_AGB[9:36, ]$area),
#   min=min(g0_AGB$area),
#   max=max(g0_AGB$area))
# (summaries <- cbind(AGB, Area))
# (plots_agb %>% filter(AGB > 98))
