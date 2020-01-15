library(silvr)
library(readr)
library(BIOMASS)

# Load data
mstems <- read_csv("~/code/biomass-espanola/mstems_genus_rough.csv", col_types = cols(plot_no = col_integer()))
mplots <- read_csv("~/code/biomass-espanola/data/haiti_biomass_v2_mplots.csv", col_types = cols(plot_no = col_integer()))
bwayo_densities <- read_csv("~/code/biomass-espanola/bwayo_densities.csv", col_types = cols(wd_avg = col_double()))
g0_plots <- read_csv("~/code/biomass-espanola/plots_zstats_07gamma0_qgis.csv")

# calculate basal area (m^2) from DBH (cm)
mstems$basal_area <- calculateBasalArea(mstems$dbh_cm)

# aggregate Basal Area by plot (m^2/ha) - the area of land that is occupied by the cross-section of a tree.
ba_plot <- aggregate(basal_area ~ plot_no, mstems, sum)
mplots <- merge(mplots, ba_plot, by='plot_no', all=TRUE)
mplots$basal_area <- mplots$basal_area / mplots$plot_area

# calculate stem volume  (m^3) from DBH (cm) - an estimate of the ammount of space that the tree bole occupies.
mstems$volume <- calculateStemVolume(mstems$dbh_cm)
vol_plot <- aggregate(volume ~ plot_no, mstems, sum)
mplots <- merge(mplots, vol_plot, by='plot_no', all=TRUE)
mplots$volume <- mplots$volume / mplots$plot_area

# calculate stocking density () from DBH
mstems$stocking <- calculateStocking(mstems$dbh_cm)
mstems$stocking_10 <- calculateStocking(mstems$dbh_cm, min_diam = 10)
stocking_plot <- aggregate(stocking ~ plot_no, mstems, sum)
mplots <- merge(mplots, stocking_plot, by='plot_no', all=TRUE)
mplots$stocking <- mplots$stocking / mplots$plot_area

# Run Chave14 equation using BIOMASS
# Get wood density by genus
genus <- mstems$genus
species <- NULL
plot <- mstems$plot_no
family <- NULL
region <- 'CentralAmericaTrop'
wd_data <- bwayo_densities
latitude <- 19
longitude <- -72
diam <- mstems$dbh_cm
height <- NULL # get this... 

# Load wood density
if (is.null(species)) species = rep('', length(genus))
wood_density_complete <- getWoodDensity(genus, species, stand = plot, family = family, region = region,
                               addWoodDensityData = wd_data, verbose = FALSE)
wood_density = wood_density_complete$meanWD

# Get height

# Prepare coordinates, if required
if (!is.null(latitude)) coord = data.frame(longitude = longitude, latitude = latitude) else coord = NULL

# Calculate AGB with Chave equation, return in kg
AGB <- computeAGB(diam, wood_density, H = height, coord = coord, Dlim = 0) * 1000
mstems$AGB_Chave14 <- AGB
summary(mstems$AGB_Chave14)

# aggregate AGB by plot (C tons/ha)
agb_plot <- aggregate(AGB_Chave14 ~ plot_no, mstems, sum)
mplots <- merge(mplots, agb_plot, by='plot_no', all=TRUE)
mplots$AGB_Chave14 <- mplots$AGB_Chave14 / mplots$plot_area

# Look at dominant species based on stocking densities. Add species_name (by_binomial)
dominant_species <- getDominantSpecies(mstems$sp_creole, mstems$plot_no)
dominant_species <- getDominantSpecies(mstems$sp_creole, mstems$plot_no, abundance = calculateBasalArea(mstems$dbh_cm))
dominant_species <- getDominantSpecies(mstems$sp_creole, mstems$plot_no, abundance = mstems$AGB_Chave14)

# Compare plot AGB to plot mean backscatter
g0_plots <- merge(mplots, g0_plots, by='plot_no', all=TRUE)
plot(g0_plots$`2007_mean`, g0_plots$AGB_Chave14)
