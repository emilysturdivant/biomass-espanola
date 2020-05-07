library(raster)
library(tidyverse)

#----
slp.fname <- "~/PROJECTS/Haiti_biomass/SRTM/slope_1arc_utm18n.tif"
aspect.fname <- "~/PROJECTS/Haiti_biomass/SRTM/aspect_1arc_utm18n.tif"
agb.fname <- "~/PROJECTS/Haiti_biomass/R_out/agb18_v8.tif"

slp.resamp.out <- "~/PROJECTS/Haiti_biomass/SRTM/slope_1arc_resamptoAGB.tif"
aspect.resamp.out <- "~/PROJECTS/Haiti_biomass/SRTM/aspect_1arc_resamptoAGB.tif"

# Compare slope to AGB ----
# Compare slope/aspect to AGB to determine whether there is a correlation.
# Load rasters
slp <- raster(slp.fname)
agb.ras <- raster(agb.fname)
crs(slp)
crs(agb.ras)

# Resample slope raster to AGB resolution
slp.proj <- projectRaster(slp, agb.ras)
slp.res <- resample(slp.proj, agb.ras)
writeRaster(slp.res, slp.resamp.out)
slp.res <- raster(slp.resamp.out)

# Scatterplot
agb.ras[agb.ras > 300] <- NA
plot(slp.res, agb.ras, pch=5, col=rgb(red=0.2, green=0.2, blue=1.0, alpha=0.1))

# Compare aspect to AGB ----
# Load rasters
aspect <- raster(aspect.fname)
agb.ras <- raster(agb.fname)
crs(aspect)
crs(agb.ras)

# Resample aspect raster to AGB resolution
aspect.proj <- projectRaster(aspect, agb.ras)
aspect.res <- resample(aspect.proj, agb.ras)
writeRaster(aspect.res, aspect.resamp.out)

# Scatterplot
plot(aspect.res, agb.ras, pch=5, col=rgb(red=0.2, green=0.2, blue=1.0, alpha=0.1))

# ----
# Brick
br <- brick(agb.ras, slp.res)
names(br) <- c('AGB', 'Slope')
br
ggplot(br, aes(x=Slope, y=AGB)) +
  geom_point(alpha=0.2) 

cor(chunk3, agb.ras, method='spearman')

#---- 
# Mask out AGB in high slope areas and look at value distributions
agb.ras.m <- agb.ras
agb.ras.m[slp.res > 30] <- NA
br <- brick(agb.ras, agb.ras.m)
names(br) <- c('AGB', 'AGB slope sub30')
agb.ras.m[slp.res > 20] <- NA
st3 <- addLayer(br, agb.ras.m)
st3 <- br2
names(st3) <- c('AGB', 'AGB slope sub30', 'AGB slope sub20')

agb.stats <- get_brick_stats(st3)

# the distribution of values in the raster
agb.samp <- agb.ras.m %>% 
  sampleRandom(100000, na.rm=TRUE) %>% 
  as.data.frame() %>% 
  rename(AGB='.')

p.hist <- ggplot(agb.samp, aes(x=AGB)) + 
  geom_histogram(col='white', fill="gray", binwidth = 5)+
  labs(x = expression(paste("Aboveground biomass (Mg ha"^"-1", ")")), 
       y = "Pixel count (out of 100,000 pixel sample)")+
  geom_vline(aes(xintercept=mean(AGB)),
             color="black", linetype="dashed", size=.5)+
  geom_vline(aes(xintercept=median(AGB)),
             color="black", linetype="dashed", size=.5)
p.hist

# Ridge density with color-coded tail probability
library(ggridges)
# Sample
sample_rasters <- function(stack, sampSize=100){
  samps <- data.frame(AGB=numeric(0),Category=numeric(0))
  for (i in seq(1,nlayers(stack))){
    samp <- stack[[i]] %>% 
      sampleRandom(sampSize, na.rm=TRUE) %>%
      as.data.frame() %>%
      rename(AGB='.')
    # add condition to rerun if no samples returned
    samp$Category <- names(stack[[i]])
    samps <- rbind(samps, samp)
  }
  return(samps)
}
samp <- sample_rasters(st3, 1000000) %>% 
  na.omit()
View(head(samp))
samp$Category <- ordered(samp$Category,
                              levels = names(st3))
ggplot(samp, aes(x = AGB, y = Category, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", 
                      calc_ecdf = TRUE, 
                      scale=8) +
  scale_x_continuous(name = "Aboveground Biomass",
                     breaks = seq(0, 200, 25),
                     limits=c(-10, 175)) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1)+
  theme_ridges()+
  #xlab("Aboveground Biomass") +
  ylab("")+
  ggtitle("Density plots of AGB by land cover type")

dp <- ggplot(samp, aes(x = AGB, y = Category)) +
  geom_density(position="identity", alpha=0.2) +
  scale_x_continuous(name = "AGB",
                     breaks = seq(0, 150, 25),
                     limits=c(0, 150)) +
  scale_y_continuous(name = "Density", limits = c(0, 0.03)) +
  ggtitle("Density plot of AGB")+
  #geom_vline(data=mu, aes(xintercept=med, color=Category),
  #linetype="dashed")+
  geom_point(data=mu, aes(med, 0), size=2)
dp
