library(silvr)
library(readr)
library(BIOMASS)

# Load data - Desktop
mstems <- read_csv("~/code/biomass-espanola/data/mstems_with_wooddensities.csv", col_types = cols(plot_no = col_integer()))
mplots <- read_csv("~/code/biomass-espanola/data/haiti_biomass_v2_mplots.csv", col_types = cols(plot_no = col_integer()))
bwayo_densities <- read_csv("~/code/biomass-espanola/data/bwayo_densities_2.csv", col_types = cols(wd = col_double()))
g0_plots <- read_csv("~/code/biomass-espanola/data/plots_g0nu2018_HV.csv")
creole_df <- read_csv("~/code/biomass-espanola/data/exploded_specieslookup.csv")

# Load data - Mac
mstems <- read_csv("~/GitHub/biomass-espanola/data/mstems_with_wooddensities.csv")
mplots <- read_csv("~/GitHub/biomass-espanola/data/haiti_biomass_v2_mplots.csv", col_types = cols(plot_no = col_integer()))
bwayo_densities <- read_csv("~/GitHub/biomass-espanola/data/bwayo_densities_2.csv", col_types = cols(wd = col_double()))
g0_plots <- read_csv("~/GitHub/biomass-espanola/data/plots_g0nu2018_HV.csv")
creole_df <- read_csv("~/GitHub/biomass-espanola/data/exploded_specieslookup.csv")

# ------------------------
# Calculate and export wood densities with different parameters for Haiti field species.
GWDBYspgnfm <- getWoodDensity(
                  genus = creole_df$genus, 
                  species = creole_df$species, 
                  family = creole_df$family, 
                  region = 'World',
                  addWoodDensityData = bwayo_densities, 
                  verbose = FALSE
                  )
write.csv(GWDBYspgnfm, "~/GitHub/biomass-espanola/getWoodDensity_GWDBYspgnfm.csv", row.names=FALSE)

# Without BwaYo wood densities
GWDspgnfm <- getWoodDensity(
                  genus=creole_df$genus, 
                  creole_df$species, 
                  family = creole_df$family, 
                  region = 'World', 
                  verbose = FALSE
                  )
write.csv(GWDspgnfm, "~/GitHub/biomass-espanola/getWoodDensity_GWDspgnfm.csv", row.names=FALSE)

# With BwaYo, without species
species = rep('', length(creole_df$genus))
GWDBYgnfm <- getWoodDensity(creole_df$genus, species, family = creole_df$family, region = 'World',
                                 addWoodDensityData = bwayo_densities, verbose = FALSE)
write.csv(GWDBYgnfm, "~/GitHub/biomass-espanola/getWoodDensity_GWDBYgnfm.csv", row.names=FALSE)

# Without BwaYo and without species
species = rep('', length(creole_df$genus))
GWDgnfm <- getWoodDensity(creole_df$genus, species, family = creole_df$family, region = 'World', verbose = FALSE)
write.csv(GWDgnfm, "~/GitHub/biomass-espanola/getWoodDensity_GWDgnfm.csv", row.names=FALSE)

# Null species and family, with BwaYo
species = rep('', length(creole_df$genus))
GWDBYgn <- getWoodDensity(creole_df$genus, species, family = NULL, region = 'World',
                                 addWoodDensityData = bwayo_densities, verbose = FALSE)
write.csv(GWDBYgn, "~/GitHub/biomass-espanola/getWoodDensity_GWDBYgn.csv", row.names=FALSE)

# Null species and only regional WDs CentralAmericaTrop, with BwaYo
species = rep('', length(creole_df$genus))
GWDBYgnfm_CAT <- getWoodDensity(creole_df$genus, species, family = creole_df$family, region = 'CentralAmericaTrop',
                                 addWoodDensityData = bwayo_densities, verbose = FALSE)
write.csv(GWDBYgnfm_CAT, "~/GitHub/biomass-espanola/getWoodDensity_GWDBYgnfm_CAT.csv", row.names=FALSE)

# Examine
sum(GWDBYspgnfm$levelWD == "family")
GWDBYspgnfm[which(GWDBYspgnfm$levelWD == 'family'), c('family', 'genus', 'species', 'meanWD', 'sdWD', 'levelWD', 'nInd')]

