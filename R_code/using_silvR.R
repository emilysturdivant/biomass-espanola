library(silvr)
library(readr)
library(BIOMASS)

# Load data - Desktop
mstems <- read_csv("~/code/biomass-espanola/mstems_genus_rough_nans.csv", col_types = cols(plot_no = col_integer()))
mplots <- read_csv("~/code/biomass-espanola/data/haiti_biomass_v2_mplots.csv", col_types = cols(plot_no = col_integer()))
bwayo_densities <- read_csv("~/code/biomass-espanola/bwayo_densities.csv", col_types = cols(wd = col_double()))
g0_plots <- read_csv("~/code/biomass-espanola/data/plots_g0nu2018.csv")

# Load data - Mac
mstems <- read_csv("~/GitHub/biomass-espanola/mstems_genus_rough_nans.csv", col_types = cols(plot_no = col_integer()))
mplots <- read_csv("~/GitHub/biomass-espanola/data/haiti_biomass_v2_mplots.csv", col_types = cols(plot_no = col_integer()))
bwayo_densities <- read_csv("~/GitHub/biomass-espanola/data/bwayo_densities.csv", col_types = cols(wd = col_double()))
g0_plots <- read_csv("~/GitHub/biomass-espanola/data/plots_g0nu2018.csv")

# aggregate Basal Area by plot (m^2/ha) - the area of land that is occupied by the cross-section of a tree.
# calculate basal area (m^2) from DBH (cm)
mstems$basal_area <- calculateBasalArea(mstems$dbh_cm)
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
species <- NULL
plot <- mstems$plot_no
family <- NULL
region <- 'CentralAmericaTrop'
wd_data <- bwayo_densities
latitude <- 19
longitude <- -72
height <- NULL # get this... 

mstems$dbh_cm <-

# Load wood density
if (is.null(species)) species = rep('', length(mstems$genus))
wood_densities <- getWoodDensity(mstems$genus, species, stand = plot, family = family, region = 'World',
                                  addWoodDensityData = wd_data, verbose = FALSE)
wood_density = wood_densities$meanWD

# Prepare coordinates, required without height
if (!is.null(latitude)) coord = data.frame(longitude = longitude, latitude = latitude) else coord = NULL

# Calculate AGB with Chave equation, return in Mg
AGB <- computeAGB(mstems$dbh_cm, wood_density, H = height, coord = coord, Dlim = 0)
mstems$AGB_BYwds <- AGB
summary(mstems$AGB_BYwds)
mstems[which(AGB == max(AGB, na.rm=T)), ]
boxplot(AGB_BYwds~plot_no, data=mstems)
boxplot(AGB_BYwds~plot_no, data=mstems)
plot(mstems$dbh_cm, mstems$AGB_BYwds, xlab='diam', ylab='AGB')

# aggregate AGB by plot (MgC/ha)
agb_plot <- aggregate(AGB_BYwds ~ plot_no, mstems, sum)

mplots <- merge(mplots, agb_plot, by='plot_no', all=TRUE)
mplots$AGB_BYwds <- mplots$AGB_BYwds / mplots$plot_area


# Add AGB to 
g0_plots <- merge(mplots, g0_plots, by='plot_no', all=TRUE)
#write.csv(g0_plots, "~/GitHub/biomass-espanola/plots_g0nu2018_withAGB.csv")

# Plot AGB against backscatter
plot(g0_plots$`2018_mean`, g0_plots$AGB_BYwds_tons, xlab='2018 HV backscatter', ylab='2019 AGB (tC/ha)')
linreg <- lm(g0_plots$AGB_BYwds ~ g0_plots$`2018_mean`)
abline(linreg)
linreg

# Get Spearman's rank correlation coefficient
corr <- cor.test(x=g0_plots$`2018_mean`, y=g0_plots$AGB_BYwds_tons, method = 'spearman')
corr$estimate


#-----------------------
# Look at dominant species based on stocking densities. Add species_name (by_binomial)
dominant_species <- getDominantSpecies(mstems$sp_creole, mstems$plot_no)
dominant_species <- getDominantSpecies(mstems$sp_creole, mstems$plot_no, abundance = calculateBasalArea(mstems$dbh_cm))
dominant_species <- getDominantSpecies(mstems$sp_creole, mstems$plot_no, abundance = mstems$AGB_Chave14)


