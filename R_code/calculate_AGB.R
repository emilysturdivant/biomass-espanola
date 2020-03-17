library(silvr)
library(readr)
library(BIOMASS)

# Load data - Desktop
mstems <- read_csv("~/code/biomass-espanola/data/mstems_with_wooddensities.csv", col_types = cols(plot_no = col_integer()))
mplots <- read_csv("~/code/biomass-espanola/data/haiti_plots_meta.csv", col_types = cols(plot_no = col_integer()))
bwayo_densities <- read_csv("~/code/biomass-espanola/data/bwayo_densities_2.csv", col_types = cols(wd = col_double()))
g0_plots <- read_csv("~/code/biomass-espanola/data/plots_g0nu2018_HV.csv")
creole_df <- read_csv("~/code/biomass-espanola/data/exploded_specieslookup.csv")

# Load data - Mac
mstems <- read_csv("~/GitHub/biomass-espanola/data/haiti_data_wds1.csv")
mplots <- read_csv("~/GitHub/biomass-espanola/data/mplots_geoms.csv", col_types = cols(plot_no = col_integer()))
g0_plots <- read_csv("~/GitHub/biomass-espanola/data/plots_g0nu2018_HV.csv")
creole_df <- read_csv("~/GitHub/biomass-espanola/data/exploded_specieslookup.csv")

# ------------------------
# Run Chave14 equation without computeAGB() 
# Compute height based on heigh-diameter model
result <- modelHD(
  D = mstems$dbh_cm,
  H = mstems$ht_m,
  useWeight = TRUE
)
result
# best model is log1 so we use that to create the model
HDmodel <- modelHD(
  D = mstems$dbh_cm,
  H = mstems$ht_m,
  method="log1",
  useWeight = TRUE
)
HDmodelPerPlot <- modelHD(
  D = mstems$dbh_cm,
  H = mstems$ht_m,
  method="log1",
  useWeight = TRUE,
  plot = mstems$plot_no
)
# too few height values in some plots
# Retrieve height data
dataHlocal <- retrieveH(
  D = mstems$dbh_cm,
  model = HDmodel
)
dataHchave <- retrieveH(
  D = mstems$dbh_cm,
  coord = mstems[, c("lon", "lat")]
)

# compute AGB(Mg) per tree
AGBtree <- computeAGB(
  D = mstems$dbh_cm,
  WD = mstems$meanWD,
  H = dataHlocal$H
)
AGBtreeChave <- computeAGB(
  D = mstems$dbh_cm, 
  WD = mstems$meanWD,
  coord = mstems[, c("lon", "lat")]
)
# compute AGB(Mg) per plot
AGBplot <- summaryByPlot(
  computeAGB(
    D = mstems$dbh_cm,
    WD = mstems$meanWD,
    H = dataHlocal$H
  ), 
  mstems$plot_no
)

AGBplotChave <- summaryByPlot(
  computeAGB(
    D = mstems$dbh_cm, 
    WD = mstems$meanWD,
    coord = mstems[, c("lon", "lat")]
  ),
  mstems$plot_no
)
plots_agb <- merge(mplots, AGBplot, by.x='plot_no', by.y='plot', all=TRUE)
plots_agb$AGB <- plots_agb$AGB / plots_agb$area_ha

# compare
summary(AGBtree - mstems$AGB)

# Plot AGB against backscatter
g0_AGB <- merge(plots_agb, g0_plots, by.x='plot_no', by.y='plot_no', all=TRUE)
plot(g0_AGB$`2018_mean`, g0_AGB$AGB, xlab='2018 HV backscatter', ylab='2019 AGB (tC/ha) from H-D model')
linreg <- lm(g0_AGB$AGB ~ g0_AGB$`2018_mean`)
abline(linreg)
linreg
# Get Spearman's rank correlation coefficient
corr <- cor.test(x=g0_AGB$`2018_mean`, y=g0_AGB$AGB, method = 'spearman')
corr$estimate

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

