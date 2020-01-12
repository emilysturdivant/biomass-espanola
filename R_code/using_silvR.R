library(silvr)

# Load data
mstems
mplots

# calculate basal area (m^2) from DBH (cm)
mstems$basal_area <- calculateBasalArea(mstems$diam)

# aggregate Basal Area by plot (m^2/ha)
ba_plot <- aggregate(basal_area ~ plotcode, mstems, sum)
ba_plot$basal_area <- ba_plot$basal_area / mplots$plot_area # check this to make sure the plots match up.

# calculate stem volume  (m^3) from DBH (cm)
mstems$volume <- calculateStemVolume(mstems$diam)
vol_plot <- aggregate(volume ~ plotcode, mstems, sum)
vol_plot$volume <- vol_plot$volume / mplots$plot_area # check this to make sure the plots match up.

# calculate stocking density () from DBH
mstems$stocking <- calculateStocking(mstems$diam)
mstems$stocking_10 <- calculateStocking(mstems$diam, min_diam = 10)
stocking_plot <- aggregate(stocking ~ plotcode, mstems, sum)
stocking_plot$stocking <- stocking_plot$stocking / mplots$plot_area # check this to make sure the plots match up.

# 
mstems <- cbind(mstems, splitGenusSpecies(mstems$species_name))

# Run Chave14 equation using BIOMASS
# Get wood density by genus
genus <- mstems$genus
species <- NULL
plot <- NULL
family <- mstems$family
region <- 'CentralAmericaTrop'
wd_data <- 
if (is.null(species)) species = rep('', length(genus))
wood_density <- getWoodDensity(genus, species, stand = plot, family = family, region = region,
                               addWoodDensityData = wd_data, verbose = FALSE)

