library(silvr)
library(readr)
library(BIOMASS)

# Load data
mstems <- read_csv("~/code/biomass-espanola/mstems_genus_rough.csv", col_types = cols(plot_no = col_integer()))
mplots <- read_csv("~/code/biomass-espanola/data/haiti_biomass_v2_mplots.csv", col_types = cols(plot_no = col_integer()))
bwayo_densities <- read_csv("~/code/biomass-espanola/bwayo_densities.csv", col_types = cols(wd_avg = col_double()))
g0_plots <- read_csv("~/code/biomass-espanola/data/plots_g0nu2018.csv")

# Load data
mstems <- read_csv("~/GitHub/biomass-espanola/mstems_genus_rough.csv", col_types = cols(plot_no = col_integer()))
mplots <- read_csv("~/GitHub/biomass-espanola/data/haiti_biomass_v2_mplots.csv", col_types = cols(plot_no = col_integer()))
bwayo_densities <- read_csv("~/GitHub/biomass-espanola/bwayo_densities.csv", col_types = cols(wd = col_double()))
g0_plots <- read_csv("~/GitHub/biomass-espanola/plots_zstats_07gamma0_qgis.csv")

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

# Prepare coordinates, required without height
if (!is.null(latitude)) coord = data.frame(longitude = longitude, latitude = latitude) else coord = NULL

# Calculate AGB with Chave equation, return in kg
AGB <- computeAGB(diam, wood_density, H = height, coord = coord, Dlim = 0) * 1000
mstems$AGB_BYwds <- AGB
summary(mstems$AGB_BYwds)

# Load wood density without Bwa Yo values
wood_density_complete <- getWoodDensity(genus, species, stand = plot, family = family, region = region, verbose = FALSE)
wood_density = wood_density_complete$meanWD

# Calculate AGB with Chave equation, return in kg
AGB <- computeAGB(diam, wood_density, H = height, coord = coord, Dlim = 0) * 1000
mstems$AGB_WDsFromCATrop <- AGB
summary(mstems$AGB_WDsFromCATrop)

# Load wood density without Bwa Yo values
wood_density_complete <- getWoodDensity(genus, species, stand = plot, family = family, region = 'World', verbose = FALSE)
wood_density = wood_density_complete$meanWD

# Calculate AGB with Chave equation, return in kg
AGB <- computeAGB(diam, wood_density, H = height, coord = coord, Dlim = 0) * 1000
mstems$AGB_WDsFromWorld <- AGB
summary(mstems$AGB_WDsFromWorld)

# aggregate AGB by plot (C kg/ha)
agb_plot <- aggregate(AGB_BYwds ~ plot_no, mstems, sum)
mplots <- merge(mplots, agb_plot, by='plot_no', all=TRUE)
mplots$AGB_BYwds <- mplots$AGB_BYwds / mplots$plot_area

# aggregate AGB by plot (C kg/ha)
agb_plot <- aggregate(AGB_WDsFromCATrop ~ plot_no, mstems, sum)
mplots <- merge(mplots, agb_plot, by='plot_no', all=TRUE)
mplots$AGB_WDsFromCATrop <- mplots$AGB_WDsFromCATrop / mplots$plot_area

# aggregate AGB by plot (C kg/ha)
agb_plot <- aggregate(AGB_WDsFromWorld ~ plot_no, mstems, sum)
mplots <- merge(mplots, agb_plot, by='plot_no', all=TRUE)
mplots$AGB_WDsFromWorld <- mplots$AGB_WDsFromWorld / mplots$plot_area

# Compare plot AGB to plot mean backscatter
g0_plots <- merge(mplots, g0_plots, by='plot_no', all=TRUE)
g0_plots$AGB_BYwds_tons <- g0_plots$AGB_BYwds/1000
plot(g0_plots$`2018_mean`, g0_plots$AGB_BYwds_tons, xlab='2018 HV backscatter', ylab='2019 AGB (tC/ha)')
scatter.smooth(x=g0_plots$`2018_mean`, y=g0_plots$AGB_BYwds_tons, main="Backscatter ~ Biomass", xlab='2018 HV backscatter', ylab='2019 AGB (tC/ha)') 
write.csv(g0_plots, "~/GitHub/biomass-espanola/plots_g0nu2018_withAGB.csv")

plot(mplots$plot_area, mplots$AGB_BYwds_tons, xlab='area', ylab='2019 AGB')


lm('2018_mean' ~ AGB_BYwds_tons, g0_plots)

# Get Spearman's rank correlation coefficient
corr <- cor.test(x=g0_plots$`2018_mean`, y=g0_plots$AGB_BYwds_tons, method = 'spearman')
corr$estimate


# Look at dominant species based on stocking densities. Add species_name (by_binomial)
dominant_species <- getDominantSpecies(mstems$sp_creole, mstems$plot_no)
dominant_species <- getDominantSpecies(mstems$sp_creole, mstems$plot_no, abundance = calculateBasalArea(mstems$dbh_cm))
dominant_species <- getDominantSpecies(mstems$sp_creole, mstems$plot_no, abundance = mstems$AGB_Chave14)


