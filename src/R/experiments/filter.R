# ---------------------------------------------------------------------------------------------
# Script to:
#     * Test use of RSAGA and RQGIS to then implement steps that are currently manual
# Proceeds:
#     * 
# Requires:
#     * 
# ---------------------------------------------------------------------------------------------

# Load libraries
library(RSAGA)
library(raster)
rsaga.env()

data(landslides)
l <- rsaga.get.libraries()
rsaga.get.modules(libs = "grid_filter")
rsaga.get.usage(lib = "grid_filter", module = "Multi Direction Lee Filter")
View(l)

g0 <- raster('results/g0nu_HV/g0nu_2018_HV_hticlip.tif')
g0 <- as.matrix(g0)
g0$header
write.sgrd(data = g0, file = file.path(tempdir(), "g0_hticlip"), write.header=F)

params = list(INPUT = file.path(tempdir(), "g0_hticlip.sgrd"),
              RESULT = file.path(tempdir(), "twi.sdat"),
              METHOD = 1)
rsaga.geoprocessor(lib = "grid_filter", module = "Multi Direction Lee Filter", 
                   param = params)


library(RQGIS)
install.packages('RQGIS3')
