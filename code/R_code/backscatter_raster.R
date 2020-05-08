#library(silvr)
library(readr)
library(raster)
library(tidyverse)
library(ggridges)
library(rgdal)
library(tmap)
library(rgdal)
library(here)

# Merge ALOS mosaic rasters ----
merge.alos.tiles <- function(path, pattern, clip_poly, filename=FALSE){
  # List filenames of ALOS mask tiles. Tiles range from 18-21, 68-75. ENVI format.
  fps <- list.files(path=path,
                    pattern=pattern,
                    full.names=TRUE,
                    recursive=TRUE,
                    include.dirs=FALSE)
  fps <- fps[1:2]
  # Load and crop rasters and store in list
  load.crop.na <- function(fp, ext) try({
    r <- raster(fp); 
    r <- crop(r, ext)
  }, silent=TRUE)
  rs <- fps %>% lapply(load.crop.na, ext=clip_poly)
  # Remove any that threw error (don't overlap with the polygon)
  # l <- vector(mode="logical", length(rs))
  # for(i in seq(1, length(rs))){
  #   if(class(rs[[i]])=='try-error') l[i] <- TRUE
  # }
  l <- lapply(function(r) if(class(r) == 'try-error') {l <- TRUE} else {l <- FALSE})
  # Merge the tiles
  dn <- do.call(merge, rs[!l])
  if(filename) writeRaster(dn, filename)
  return(dn)
}
# Load island polygons
isl_poly <- readOGR(dsn="data/contextual_data/Hispaniola", layer='Hisp_adm0')
isl_poly <- readOGR(dsn=here("data", "Hispaniola"), layer='Hisp_adm0')
isl_poly <- readOGR(dsn="data/biota_out/vector", layer='hisp18_maskLand_clean')
isl_poly <- sf::st_read("data/biota_out/vector/hisp18_maskLand_clean.shp")
isl_poly <- sf::st_read(here("data", 'hisp18_maskLand_clean.shp'))
sum(area(isl_poly))/1000000



plot(dn_mask)
dn_date <- merge.alos.tiles('data/ALOS', '_18_date_F02DAR$', 
                            isl_poly, 'data/R_out/hisp18_date.tif')
# writeRaster(dn_date, 'data/R_out/hisp18_date.tif')
plot(dn_date)
dn_linci <- merge.alos.tiles('data/ALOS', '_18_linci_F02DAR$', 
                             isl_poly, 'data/R_out/hisp18_localincidence.tif')
# writeRaster(dn_linci, 'data/R_out/hisp18_localincidence.tif')
dn_mask <- merge.alos.tiles('data/ALOS', '_18_mask_F02DAR$', 
                            isl_poly, 'data/R_out/hisp18_mask.tif')
dn_mask[dn_mask<75] <- NA
dn_mask[!is.na(dn_mask)] <- 1
writeRaster(dn_mask, 'data/R_out/hisp18_maskLand.tif')
mask_land <- raster('data/R_out/hisp18_maskLand.tif')

# Load previously created mosaics
dn_linci <- raster('data/R_out/hisp18_localincidence.tif')
dn_mask <- raster('data/R_out/hisp18_mask.tif')
dn_date <- raster('data/R_out/hisp18_date.tif')
plot(dn_linci)
plot(isl_poly, add=TRUE)

# Trying to use tmap, but the files are too big.
lyr_linci <- tm_shape(dn_linci) + tm_raster()
lyr_mask <- tm_shape(dn_mask) + tm_raster()
lyr_date <- tm_shape(dn_date) + tm_raster()
lyr_isl <- tm_shape(isl_poly) + tm_borders() 
tmap_mode("plot")
(map_linci <- lyr_linci + lyr_isl)
map_all <- tmap_arrange(lyr_linci + lyr_isl, 
                        lyr_date + lyr_isl, 
                        lyr_date + lyr_isl)
map_all

# Look at proportion of values in each mask category
vals <- getValues(dn_mask)
f <- as.factor(vals)
levels(f)
df <- data.frame(
  group = c("Normal", "Layover", "Shadowing"), 
  value = c(sum(vals==255, na.rm=TRUE), 
            sum(vals==100, na.rm=TRUE),
            sum(vals==150, na.rm=TRUE)))
df$value[2] / sum(df$value)
df$value[3] / sum(df$value)
setwd("data/R_out")
save.image()
load(".RData")

# Barplot
(bp <- ggplot(df, aes(x="", y=value, fill=group))+
    geom_bar(width = 1, stat = "identity")+ 
    scale_fill_manual(values=c("#E69F00", "#999999")) +
    theme_minimal())
(pie <- bp + coord_polar("y", start=0))

# Count inland NA value in biota output ----
# Load files
g0 <- raster("data/biota_out/g0nu_2018_HV.tif")
hti_poly <- readOGR(dsn="data/contextual_data/HTI_adm", layer='HTI_adm0')
isl_poly <- readOGR(dsn="data/contextual_data/Hispaniola", layer='Hisp_adm0')

# Crop backscatter and reclass 0s to NA
g0 <- crop(g0, hti_poly)
g0[g0==0] <- NA
writeRaster(g0, "data/biota_out/g0nu_2018_HV_haitiR.tif", overwrite=TRUE)
g0 <- raster("data/biota_out/g0nu_2018_HV_haitiR.tif")

# Convert island to land mask
tmp_isl <- rasterize(isl_poly, g0, field=1) # land==1; sea==NA

# Multiply backscatter classes x polygon classes 
# reclass to (NA==NA; 1==land; 2==sea)
msk <- overlay(g0, tmp_isl, 
               fun=function(r1, r2){
                 r1[!is.na(r1)] <- 0 # values in land-->0
                 r1[is.na(r1)] <- 1 # NAs in land-->1
                 r2[is.na(r2)] <- 2 # sea-->2 (land==1)
                 r <- r1*r2 # create output (sea==2)
                 r[r==1] <- NA # NAs in land --> NA
                 r[r==0] <- 1 # valid land --> 1
                 return(r)}
)
plot(msk) 
# define as categorical variable
f <- as.factor(msk) # or ratify(msk) ?
x <- levels(f)[[1]]
x$code <- c("land", "sea")
levels(f) <- x
writeRaster(f, "data/biota_out/mask_NAinland.tif")

# Land mask
land <- msk
land[is.na(land)] <- 1
land[land==2] <- NA
writeRaster(land, "data/biota_out/mask_land18.tif")
land <- raster("data/biota_out/mask_land18.tif")
land_poly <- rasterToPolygons(land) # very slow, much faster with gdal_polygonize in QGIS 
land_poly <- readOGR(dsn="data/biota_out/vector", layer='mask_land18_hti')
land_area <- sum(area(land_poly)) / 1000000

# Get values and count NAs
vals <- getValues(msk)
df <- data.frame(
  group = c("Land", "Nulls"), 
  value = c(sum(vals==1, na.rm=TRUE), 
            sum(is.na(vals)))
)
na_pct <- df$value[2] / (df$value[1] + df$value[2])

# Convert to area
msk[is.na(msk)] <- 0
cell_size <- area(msk, na.rm=TRUE, weights=FALSE)
cell_size <- cell_size[!cell_size==2]
land_area<-length(cell_size)*median(cell_size)
cell_size <- cell_size[is.na(cell_size)]
NA_area<-length(cell_size)*median(cell_size)


mask_ras <- raster("data/biota_out/mask_mosaic18_nosea.tif")

# Get counts
vals <- getValues(mask_ras)
na_ct <- sum(is.na(vals))
nonan_ct <- sum(!is.na(vals))
sea_ct <- sum(vals==0, na.rm=TRUE)
valid_ct <- sum(vals==1, na.rm=TRUE)

sea_ct + valid_ct == nonan_ct
na_ct/valid_ct

# Merge and resample LULC for Hispaniola ----
fps <- c("data/LULC/Haiti2017_Clip.tif", 
        "data/LULC/DR_2017_clip.tif")
rs <- lapply(fps, 
             function(fp) {
               r <- raster(fp)
               r <- resample(r, dn_mask)
               system(paste0('gdalwarp {}'))
               return(r)}
             )
rs
# Merge the tiles
lc <- do.call(merge, rs)
writeRaster(lc, "data/LULC/Hisp_2017_resALOS.tif")
lc <- raster("data/LULC/Hisp_2017_resALOS.tif")

# Make water and urban mask
mask_land <- raster('data/R_out/hisp18_maskLand.tif') # 1 for land, NA for everything else
# mask_land <- crop(mask_land, extent(-73, -72, 19, 20))
# lc <- crop(lc, extent(-73, -72, 19, 20))
# Overlay
mask_land_wu <- overlay(mask_land, lc, 
               fun=function(r1, r2){
                 r2[r2<3] <- NA
                 r2[!is.na(r2)] <- 1
                 r <- r1*r2
                 return(r)}
)
# Check results - plot and see how NA count changed
plot(mask_land_wu)
newnas <- sum(is.na(getValues(mask_land_wu)), na.rm=TRUE) - 
  sum(is.na(getValues(mask_land)), na.rm=TRUE)
df <- data.frame(
  group = c("AllLand", "WaterUrban"), 
  value = c(sum(!is.na(getValues(mask_land)), na.rm=TRUE), 
            newnas)
  )
df$value[2] / df$value[1]


# Look at AGB map ----
# Load data - Desktop
# Load data - Mac
g0_AGB <- read_csv("~/GitHub/biomass-espanola/data/plots_g0nu_AGB.csv")
# just get the two columns we care about
g0.agb <- as.data.frame(cbind(g0_AGB$AGB_ha, g0_AGB$'2018mean')) %>% 
  rename(AGB = V1, backscatter =V2)
# Linear model
ols <- lm(g0.agb$AGB ~ g0.agb$backscatter, x=TRUE, y=TRUE)

# Extract backscatter values to polygons ----
g0 <- raster("data/biota_out/g0nu_2018_haiti_qLee1.tif")

polys <- readOGR(dsn="~/GitHub/biomass-espanola/data", layer='AllPlots')
ex <- extract(g0, polys)
polys$g0l_mean <- unlist(lapply(ex, function(x) if (!is.null(x)) mean(x, na.rm=TRUE) else NA ))
count <- unlist(lapply(ex, function(x) length(x)))
polys$g0l_count <- count
View(cbind(polys$X2018mean, polys$g0l_mean))
View(cbind(polys$X2018count, polys$g0l_count))

# Create 95% CI raster ----
# Function to create raster of confidence intervals
predict.ci.raster <- function(g0, model){
  names(g0) <- 'backscatter'
  df <- as.data.frame(g0)
  chunks <- split(df, (seq(nrow(df))-1) %/% 1000000) 
  agblist <- list()
  for (i in 1:length(chunks)){
    chunk <- chunks[[i]]
    agb1 <- stats::predict(model, newdata=chunk, 
                           interval="confidence", level = 0.95) %>%
      as.data.frame()
    agblist[[i]] <- agb1
  }
  agb <- bind_rows(agblist)
  agb.ci.v <- as.vector((agb$upr - agb$lwr)/2)
  ci.ras <- setValues(g0, values=agb.ci.v)
}
ras.agb.ci <- predict.ci.raster(g0, ols)
ras.agb.ci[agb.ras > 310] <- NA
ras.agb.ci[agb.ras < 0] <- NA
writeRaster(ras.agb.ci, "data/R_out/agb_CI_sub310.tif", overwrite=TRUE)

# Get AGB percentiles ----
agb <- raster("data/biota_out/AGB_2018_v5q_haiti.tif")
agb <- raster("data/biota_out/agb_2018_v6r.tif")

#agb[agb < =0] <- NA
qs.agb <- quantile(agb, probs=c(0, 0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999, 1))
View(qs.agb)
agb.q90s <- quantile(agb, probs=c(0.91, 0.92, 0.93, 0.94, 0.95))
View(agb.q90s)
# the distribution of values in the raster
agb.samp <- agb %>% 
  sampleRandom(100000, na.rm=TRUE) %>% 
  as.data.frame() %>% 
  rename(AGB='.')
agb.samp.haiti <- agb.samp

# Boxplot ----
agb.qs <- quantile(agb, probs=c(0, 0.02, 0.25,0.5,0.75, 0.98, 1))
agb.qs
iqr <- agb.qs[['75%']] - agb.qs[['25%']]
low.inner.fence <- agb.qs[['25%']] - iqr*1.5
low.outer.fence <- agb.qs[['25%']] - iqr*3
upper.inner.fence <- agb.qs[['75%']] + iqr*1.5
upper.outer.fence <- agb.qs[['75%']] + iqr*3
DF <- data.frame(x=c("AGB"), min=agb.qs[['2%']], 
                 low=agb.qs[['25%']], 
                 mid=agb.qs[['50%']], 
                 top=agb.qs[['75%']], 
                 max=agb.qs[['98%']])
ggplot(DF, aes(x=x, ymin = min, lower = low, middle = mid, upper = top, ymax = upper.inner.fence)) +
  geom_boxplot(stat = "identity")

# Histogram and density plots
p1 <- ggplot(agb.samp, aes(x=AGB)) + 
  geom_histogram(col='white', fill="gray", binwidth = 5)+
  labs(x = expression(paste("Aboveground biomass (MgC ha"^"-1", ")")), 
       y = "Pixel count (out of 100,000 pixel sample)")+
  geom_vline(aes(xintercept=mean(AGB)),
             color="black", linetype="dashed", size=.5)+
  geom_vline(aes(xintercept=median(AGB)),
             color="black", linetype="dashed", size=.5)
p1
p.dens.agb <- ggplot(agb.samp, aes(x=AGB)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)+ 
  labs(x = expression(paste("Aboveground biomass (MgC ha"^"-1", ")")), 
       y = "Density")+
  scale_x_continuous(breaks = seq(0, 300, by = 50))+
  geom_vline(aes(xintercept=mean(AGB)), color="black", linetype="dashed", size=.5)+
  geom_vline(aes(xintercept=median(AGB)),
             color="black", linetype="dashed", size=.5)
p.dens.agb
p <- ggplot(agb.samp, aes(x=AGB)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white", binwidth = 5)+
  geom_density(alpha=.2, fill="#FF6666") + 
  geom_vline(aes(xintercept=mean(AGB)),
              color="black", linetype="dashed", size=.5)
p
mean(agb.samp$AGB)
median(agb.samp$AGB)

#----
# Get zonal statistics in R. Produces different values than those from Q. 
lc <- raster("data/LULC/Haiti2017_Clip.tif")
# Resample LC raster to AGB resolution
lc.c <- resample(lc, agb.ras)
writeRaster(lc.c, filename="data/LULC/Haiti2017_Clip_AGBres.tif")
# Get means and SDs for LC class
lc.means <- zonal(agb.ras, lc.c)
lc.sds <- zonal(agb.ras, lc.c, fun='sd')
lc.zonal <- cbind(lc.means, lc.sds)
View(t(lc.zonal))

# Get AGB percentiles for each LULC class ----
get_brick_stats <- function(lc.br){
  # Get selection of percentiles for each LC
  agb.qs <- data.frame(row.names=c('0%', '1%','2%','10%', '25%', '50%', '75%', '90%','98%', '99%','100%'))
  for (i in seq(1,nlayers(lc.br))){
    lc <- names(lc.br[[i]])
    agb.qs[[lc]] <- quantile(lc.br[[i]], probs=c(0, 0.01, 0.02, 0.1, 0.25, 0.5, 0.75, 0.9,0.98, 0.99, 1))
  }
  # Get means, SDs, and Skews
  stats <- list(mean=cellStats(lc.br, stat='mean', na.rm=TRUE),
                sd=cellStats(lc.br, stat='sd', na.rm=TRUE),
                skew=cellStats(lc.br, stat='skew', na.rm=TRUE)) %>%
    bind_cols()%>%
    t()
  colnames(stats) <- names(lc.br)
  agb.stats <- rbind(agb.qs, stats)
}
lc.br <- brick(raster("data/biota_out/agb_2018_v6_mask2share.tif"), 
               raster("data/biota_out/agb_2018_v6CI_2share.tif"),
               raster("data/LULC/Haiti2017_water_agb18v6.tif"), 
               raster("data/LULC/Haiti2017_urban_agb18v6.tif"), 
               raster("data/LULC/agb_lc3.tif"),
               raster("data/LULC/agb_lc4.tif"),
               raster("data/LULC/agb_lc5.tif"), 
               raster("data/LULC/agb_lc6.tif"),
               raster("data/LULC/agb_lc_over3.tif"))
names(lc.br) <- c('Haiti', 'CI', 'Water', 'Urban', 'Bareland', 'Tree cover', 'Grassland', 'Shrubs', 'Veg')
saveRDS(lc.br, file = "data/R_out/brick_AGBv6_withLC.rds")
lc.br <- readRDS(file = "data/R_out/brick_AGBv6_withLC.rds")
lc.br

agb.stats <- get_brick_stats(lc.br)
# saveRDS(agb.qs, file = "data/R_out/agb_quantiles_v6.rds")
# agb.qs <- readRDS(file = "data/R_out/agb_quantiles_v6.rds")
# View(agb.qs)
saveRDS(agb.stats, file = "data/R_out/agb_stats_v6.rds")
agb.stats <- readRDS(file = "data/R_out/agb_stats_v6.rds")

# Sample AGB by LC rasters for plotting ----
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
lc.samp <- sample_rasters(lc.br, 1000000)
lc.sampNA <- lc.samp %>% na.omit()
lc.sampNA$Category <- ordered(lc.sampNA$Category,
                              levels = names(lc.br))
saveRDS(lc.sampNA, file = "data/R_out/lc_samp1mil_2share.rds")
lc.sampNA <- readRDS("data/R_out/lc_samp_1mil_v6.rds")

lc.samp.sub100 <- lc.sampNA[which(lc.sampNA$AGB <= 100),]
lc.samp.sub200 <- lc.sampNA[which(lc.sampNA$AGB <= 200),]

# ANOVA with multiple pair-wise comparison ----
sampLC <- lc.sampNA[-c('Haiti', 'CI', 'Veg', 'Water')]
sampLC <- lc.sampNA[lc.sampNA$Category != 'Haiti', , drop=FALSE]
sampLC <- sampLC[sampLC$Category != 'CI', , drop=FALSE]
sampLC <- sampLC[sampLC$Category != 'Veg', , drop=FALSE]
res.aov <- aov(AGB ~ Category, data = sampLC)
summary(res.aov)
tkhsd <- TukeyHSD(res.aov)
op <- par(mar = c(4,10,4,2) + 0.1)
plot(tkhsd, las=1)
pairwise.t.test(lc.sampNA$AGB, lc.sampNA$Category,
                p.adjust.method = "BH")
plot(res.aov, 1)
library(car)
leveneTest(AGB ~ Category, data = lc.sampNA)

lc.sampNA[283224,]
lc.sampNA[8052,]
lc.sampNA[9444,]

res.aov <- aov(AGB ~ Category, data = lc.sampNA)
summary(res.aov)
TukeyHSD(res.aov)
pairwise.t.test(lc.sampNA$AGB, lc.sampNA$Category,
                p.adjust.method = "BH")

library(ggpubr)
ggline(lc.sampNA, x = "Category", y = "AGB", 
       add = c("mean_se", "jitter"), 
       order = c("water", "urban", "bareland", "grassland", "shrubs", "tree.cover"),
       ylab = "AGB", xlab = "Land cover")

# Density plot
# Get group means and medians
mu <- lc.samp.sub200 %>%
  group_by(Category) %>%
  summarise(mean = mean(AGB), med=median(AGB))
# Density plot, overlaid on one axis and color-coded
custom.col <- c(water="#6699CC", urban="#888888", bareland="#661100", 
                grassland="#DDCC77", shrubs="#999933", tree.cover="#117733")
dp <- ggplot(lc.samp.sub200, aes(x = AGB, color = Category, fill = Category)) +
  geom_density(position="identity", alpha=0.2) +
  scale_x_continuous(name = "AGB",
                     breaks = seq(0, 150, 25),
                     limits=c(0, 150)) +
  scale_y_continuous(name = "Density", limits = c(0, 0.03)) +
  ggtitle("Density plot of AGB by land cover")+
  #geom_vline(data=mu, aes(xintercept=med, color=Category),
             #linetype="dashed")+
  geom_point(data=mu, aes(med, 0), size=2) +
  scale_color_manual(values=custom.col)
dp
# Ridge density plot with quantile lines
dp <- ggplot(lc.samp.sub200, aes(x = AGB, y = Category, fill = Category)) +
  #geom_density_ridges(scale = 8, alpha=0.7) +
  stat_density_ridges(quantile_lines = TRUE, scale = 8, alpha=0.7)+
  theme_ridges() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  ) +
  xlab("Aboveground Biomass") +
  ylab("")+
  ggtitle("Density of AGB by land cover type")
dp
# Ridge density with color-coded tail probability
ggplot(lc.sampNA, aes(x = AGB, y = Category, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
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

# Boxplot
# Get group means and medians
mu <- lc.sampNA %>%
  group_by(Category) %>%
  summarise(mean = mean(AGB), med=median(AGB))
custom.col <- c(water="#2980B9", urban="#99A3A4", bareland="#BA4A00", 
                grassland="#F8C471", shrubs="#45B39D", tree.cover="#1E8449")
custom.col <- c(water="#38A6A5", urban="#666666", bareland="#CC503E", 
                grassland="#EDAD08", shrubs="#73AF48", tree.cover="#0F8554")
custom.col <- c(water="#6699CC", urban="#888888", bareland="#661100", 
                grassland="#DDCC77", shrubs="#999933", tree.cover="#117733")
bp <- ggplot(lc.sampNA, aes(x=Category, y=AGB, fill=Category)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(alpha=0.9, notch=TRUE, outlier.alpha = 0.1) +
  stat_summary(fun=mean, geom="point", shape=5, size=4) +
  scale_fill_manual(values=custom.col) +
  scale_y_continuous(name = "Aboveground Biomass",
                     breaks = seq(0, 200, 25),
                     limits=c(-10, 200)) +
  coord_flip()+
  xlab("")+
  theme_minimal()+
  theme(legend.position="none")
bp


# Old LULC analysis----
lulc.samp2 <- lulc.samp %>% 
  pivot_longer(Water:Shrubs, 
               names_to='Land.cover', 
               values_to='AGB') %>% 
  drop_na()

# Boxplot
lulc.qs <- cbind(Water, Urban, Bareland, Treecover, Grassland, Shrubs) %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column('Land.cover')

lulc.qs <- lulc.qs %>% pivot_longer('0%':'100%', 'Percentile', 'Value')
lulc.qs$Percentile <- as.numeric(gsub("%", "", lulc.qs$Percentile))*0.01

ggplot(lulc.qs, aes(x=Percentile, y=Value, group=Land.cover)) +
  geom_line(aes(color=Land.cover))+
  geom_point(aes(color=Land.cover))

ggplot(lulc.qs, aes(x=Land.cover, 
                    ymin = '0%', 
                    lower = '25%', 
                    middle = '50%', 
                    upper = '75%', 
                    ymax = '100%')) + 
  geom_boxplot(stat = "identity")


# Backscatter values ----
g0 <- raster("data/biota_out/g0nu_2018_nofilt_HV.tif")
NAvalue(g0) <- 0
qs.g0 <- quantile(g0, probs=c(0, 0.1,0.5,0.9,0.99, 0.999, 1))
qs.g0[['100%']]
# There's a really high value
NAvalue(g0) <- qs.g0[['100%']]
qs.g0 <- quantile(g0, probs=c(0, 0.1,0.5,0.9,0.99, 0.999, 1))
qs.g0

# Get interquartile range etc.
g0.qs <- quantile(g0, probs=c(0, 0.02, 0.25,0.5,0.75, 0.98, 1))
g0.qs
iqr <- g0.qs[['75%']] - g0.qs[['25%']]
low.inner.fence <- g0.qs[['25%']] - iqr*1.5
low.outer.fence <- g0.qs[['25%']] - iqr*3
upper.inner.fence <- g0.qs[['75%']] + iqr*1.5
upper.outer.fence <- g0.qs[['75%']] + iqr*3
upper.outer.fence
g0.qs[['100%']]

# the distribution of values in the raster
g0.samp <- g0 %>% 
  sampleRandom(100000, na.rm=TRUE) %>% 
  as.data.frame(names=c('backscatter'))
g0.samp <- g0.samp2
p.dens.g0 <- ggplot(g0.samp, aes(x=backscatter)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)+ 
  labs(x = expression(paste("Radar backscatter, ", sigma['HV']^0, " (m"^2, "/m"^2, ")")), 
       y = "Density")+
  scale_x_continuous(breaks = seq(0, .3, by = .05))+
  xlim(0,0.3)+
  geom_vline(aes(xintercept=mean(backscatter)), color="black", linetype="dashed", size=.5)+
  geom_vline(aes(xintercept=median(backscatter)),
             color="black", linetype="dashed", size=.5)
p.dens.g0
