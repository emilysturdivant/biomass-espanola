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
mplots <- read_csv("~/GitHub/biomass-espanola/data/haiti_plots_meta.csv", col_types = cols(plot_no = col_integer()))
#bwayo_densities <- read_csv("~/GitHub/biomass-espanola/data/bwayo_densities_2.csv", col_types = cols(wd = col_double()))
g0_plots <- read_csv("~/GitHub/biomass-espanola/data/plots_g0nu2018_HV.csv")
creole_df <- read_csv("~/GitHub/biomass-espanola/data/exploded_specieslookup.csv")

# ------------------------
# Get mean wood densities for field data - aggregate creole means and plot means.
# Get wood density (parameters: GWD=World, Bwa Yo, genus, family)
species = rep('', length(creole_df$genus))
wood_densities <- getWoodDensity(creole_df$genus, species, family = creole_df$family, region = 'World',
                                 addWoodDensityData = bwayo_densities, verbose = FALSE)

# Replace dataset means with NA
wood_densities[which(wood_densities$levelWD == 'dataset'), c('meanWD', 'sdWD', 'levelWD', 'nInd')] <- NA

# Calculate creole means and join to mstems
creole_df$meanWD <- wood_densities$meanWD
creole_wds <- aggregate(meanWD ~ creole, creole_df, mean)
mstems <- merge(mstems, creole_wds, by.x='sp_creole', by.y='creole', all.x=TRUE, sort=FALSE)

# Calculate plot means and use them to replace NA values. 
plot_wds <- aggregate(meanWD ~ plot_no, mstems, mean)
mstems <- merge(mstems, plot_wds, by='plot_no', all.x=TRUE, sort=FALSE, suffixes=c("", "_plot"))
mstems$meanWD[is.na(mstems$meanWD)] <- mstems$meanWD_plot[is.na(mstems$meanWD)]
mstems <- subset(mstems, select = -meanWD_plot)


# ------------------------
# Run Chave14 equation 
# Prepare coordinates, required without height
long <- c(19.00)
lat <- c(-72.00)

# Compute E
computeE(cbind(long, lat))

# Calculate AGB with Chave equation, return in Mg
mstems$AGB_GWDBYgnfm <- computeAGB(mstems$dbh_cm, mstems$meanWD, H = NULL, coord = coord, Dlim = 0)

# aggregate AGB by plot (MgC/ha)
agb_plot <- aggregate(AGB_GWDBYgnfm ~ plot_no, mstems, sum)
mplots <- merge(mplots, agb_plot, by='plot_no', all=TRUE)
mplots$AGB_GWDBYgnfm <- mplots$AGB_GWDBYgnfm / mplots$plot_area

# ------------------------
# Run Chave14 equation using BIOMASS
# Get wood density by genus
species <- NULL
if (is.null(species)) species = rep('', length(mstems$genus))
family <- NULL
#region <- 'CentralAmericaTrop'
latitude <- 19
longitude <- -72

# convert Nulls in DBH to...?
mstems$dbh_cm <-

# Load wood density
wood_densities <- getWoodDensity(mstems$genus, species = mstems$species, stand = mstems$plot_no, family = mstems$family, region = 'World',
                                 addWoodDensityData = bwayo_densities, verbose = FALSE)

# Prepare coordinates, required without height
if (!is.null(latitude)) coord = data.frame(longitude = longitude, latitude = latitude) else coord = NULL

# Calculate AGB with Chave equation, return in Mg
AGB <- computeAGB(mstems$dbh_cm, wood_densities$meanWD, H = height, coord = coord, Dlim = 0)
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
plot(g0_plots$`2018_mean`, g0_plots$AGB_BYwds, xlab='2018 HV backscatter', ylab='2019 AGB (tC/ha)')
linreg <- lm(g0_plots$AGB_BYwds ~ g0_plots$`2018_mean`)
abline(linreg)
linreg

# Get Spearman's rank correlation coefficient
corr <- cor.test(x=g0_plots$`2018_mean`, y=g0_plots$AGB_BYwds, method = 'spearman')
corr$estimate


#-----------------------
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

# Look at dominant species based on stocking densities. Add species_name (by_binomial)
dominant_species <- getDominantSpecies(mstems$sp_creole, mstems$plot_no)
dominant_species <- getDominantSpecies(mstems$sp_creole, mstems$plot_no, abundance = calculateBasalArea(mstems$dbh_cm))
dominant_species <- getDominantSpecies(mstems$sp_creole, mstems$plot_no, abundance = mstems$AGB_Chave14)

