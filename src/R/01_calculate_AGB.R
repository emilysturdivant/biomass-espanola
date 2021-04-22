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

# Load libraries
library(readr)
library(BIOMASS)
library(tidyverse)
library(sf)

results_dir <- 'data/results'

# Filenames
plots_shp <- file.path(results_dir, "plots_values/all_plots.shp")

# Load CSVs
mstems <- read_csv("data/species_and_wds/haiti_data_wds2.csv")
mplots <- read_csv("data/species_and_wds/mplots_geoms.csv", col_types = cols(plot_no = col_integer()))
creole_df <- read_csv("data/species_and_wds/exploded_specieslookup.csv")

# Consolidate plot polygons ---- ######################################################
# Extract plot values from stems
plots_no <- mstems %>% group_by(plot_no, plot_shp) %>% 
  summarize()
# Revise shp names to match raw data
plots_no$plot_shp <- plots_no$plot_shp %>% 
  str_replace('Pandiassou2', '2pandiassou') %>% 
  str_replace('Campeche', 'Camapeche')

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

# Load and merge plot polygons into one DF
fps <- list.files(path="data/raw/survey_plot_shps", pattern="\\.shp$", full.names=TRUE)
plots <- fps %>% lapply(standardize.names)
plots <- do.call(rbind, plots)

# Join to plot_no from mstems
plot_polys <- plots %>% full_join(plots_no, by=c('Name' ='plot_shp')) %>% 
  filter(!is.na(plot_no)) 
plot_polys %>% 
  st_write(plots_shp, append=FALSE)

plot_polys <- st_read(plots_shp)

# Look at data ---- ####################################################################
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

# Boxplots per plot ----
(p <-ggplot(mstems, aes(y=dbh_cm, x=plot_no, group=plot_no)) + 
   geom_boxplot() +
   labs(y = expression(paste("DBH (cm)")), 
        x = "")+
   ggtitle("Tree diameter at breast height (DBH; N = 6,256)")+
   theme_minimal())
ggsave('figures/qc_plots/boxes_byplot_dbh.png', width=6, height=2.5)
(p <-ggplot(mstems, aes(y=ht_m, x=plot_no, group=plot_no)) + 
    geom_boxplot() +
    labs(y = expression(paste("Height (m)")), 
         x = "")+
    ggtitle("Tree heights (N = 2,843)")+
    theme_minimal())
ggsave('figures/qc_plots/boxes_byplot_ht.png', width=6, height=2.5)

(p <-ggplot(mstems, aes(x=dbh_cm, y=plot_no, group=plot_no)) + 
    geom_boxplot() +
    labs(x = expression(paste("DBH (cm)")), 
         y = "")+
    ggtitle("Tree diameter \n(DBH; N = 6,256)")+
    theme_minimal()+
    theme(axis.text.y=element_blank()))
ggsave('figures/qc_plots/boxes_byplot_dbhX.png', width=2.5, height=6)
(p <-ggplot(mstems, aes(x=ht_m, y=plot_no, group=plot_no)) + 
    geom_boxplot() +
    labs(x = expression(paste("Height (m)")), 
         y = "")+
    ggtitle("Tree heights (N = 2,843)")+
    theme_minimal()+
    theme(axis.text.y=element_blank()))
ggsave('figures/qc_plots/boxes_byplot_htX.png', width=2.5, height=6)

mstems_temp <- mstems %>% group_by(plot_no) %>% 
  summarize(mean_ht = median(ht_m, na.rm=T)) %>% 
  arrange(desc(mean_ht)) %>% 
  left_join(mstems) %>% 
  mutate(plot_no = as.factor(plot_no))
(p <-ggplot(mstems_temp, aes(y=ht_m, x=plot_no, group=plot_no)) + 
    geom_boxplot() +
    labs(y = expression(paste("Height (m)")), 
         x = "")+
    ggtitle("Tree heights (N = 2,843)")+
    theme_minimal()+
    theme(axis.text.x=element_blank()))

# Compute AGB ---- ##################################################################
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

# Save
saveRDS(mstems, file.path(results_dir, 'R_out/mstems_agb.rds'))
mstems <- readRDS(file.path(results_dir, 'R_out/mstems_agb.rds'))

# Compute AGB per plot and convert to AGB per hectare ---- #########################################
# compute AGB(Mg) per plot
AGBplot <- summaryByPlot(
  computeAGB(
    D = mstems$dbh_cm,
    WD = mstems$meanWD,
    H = mstems$H
  ), 
  mstems$plot_no
)

# Calculate AGB per hectare
plots_agb <- merge(plot_polys, AGBplot, by.x='plot_no', by.y='plot', all=TRUE)
plots_agb$AGB_ha <- plots_agb$AGB / plots_agb$area
plots_agb$AGB_ha[is.na(plots_agb$AGB_ha)] <- 0

saveRDS(plots_agb, file.path(results_dir, 'R_out/plots_agb.rds'))

# Look at data ---- ####################################################################
summary(mstems$meanWD)
wd_sd=sd(mstems$meanWD, na.rm=TRUE)
summary(mstems$agb)
sd(mstems$agb, na.rm=TRUE)

# Plot histograms and density plot
mstems.filt <- mstems %>% filter(agb < 10)
(p <-ggplot(mstems.filt, aes(x=agb)) + 
  geom_histogram(binwidth=0.05) +
  labs(y = "")+
  scale_x_continuous(name = "AGB") +
  ggtitle("Histogram of AGB (N = 6,256, bin width = 0.05 m)")+
  theme_minimal())
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
(p.dens.agb <- ggplot(mstems.filt, aes(x=agb)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)+ 
  labs(x = expression(paste("Aboveground biomass (MgC)")), 
       y = "Density"))
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


# Look at data ---- ####################################################################
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


(summaries <- cbind(Mean = c(dbh=mean(mstems$dbh_cm, na.rm=TRUE),
                              ht=mean(mstems$H, na.rm=TRUE),
                              wd=mean(mstems$meanWD, na.rm=TRUE),
                              sdWD=mean(mstems$sdWD, na.rm=TRUE),
                              agb=mean(mstems$agb, na.rm=TRUE)), 
                    SD = c(dbh=sd(mstems$dbh_cm, na.rm=TRUE), 
                             ht=sd(mstems$H, na.rm=TRUE),
                             wd=sd(mstems$meanWD, na.rm=TRUE),
                             sdWD=NA, 
                             agb=sd(mstems$agb, na.rm=TRUE))))

# Look at values
plots_agb$AGB <- as.vector(plots_agb$AGB_ha)
plots_agb$area <- as.vector(plots_agb$area)
AGB <- c(
  all_mean=mean(g0.agb$AGB),
  all_sd=sd(g0.agb$AGB),
  bio_mean=mean(g0.agb[9:36, ]$AGB),
  bio_sd=sd(g0.agb[9:36, ]$AGB),
  min=min(g0.agb$AGB),
  max=max(g0.agb$AGB))
Area <- c(
  all_mean=mean(g0_AGB$area),
  all_sd=sd(g0_AGB$area),
  bio_mean=mean(g0_AGB[9:36, ]$area),
  bio_sd=sd(g0_AGB[9:36, ]$area),
  min=min(g0_AGB$area),
  max=max(g0_AGB$area))
(summaries <- cbind(AGB, Area))
(plots_agb$AGB[plots_agb$AGB > 98])
