# ---------------------------------------------------------------------------------------------
# Script to:
#     * Evaluate influence of slope/aspect on AGB/backscatter
# Proceeds:
#     * regression_AGB-g0.R - creates AGB map
#     * calculate_AGB.R - calculates AGB by plot from field data
#     * process_ALOS_tiles.R
# Requires:
#     * AGB map (agbr)
#     * backscatter (g0)
#     * slope
#     * aspect
# ---------------------------------------------------------------------------------------------

# Load libraries
library(raster)
library(tidyverse)
library(stars)
library(tmap)

# Reproject and resample slope and aspect to match AGB ------------------------------------
agb.ras <- read_stars("results/tifs_by_R/agb18_v1_l0.tif")

# slope
slp1 <- read_stars("data/SRTM/slope_1arc_utm18n.tif")
slp <- slp1 %>% st_warp(agb.ras, method="bilinear", use_gdal=TRUE)
slp %>% saveRDS('results/R_out/slope_1arc_resamptoAGBbl.rds')

# aspect
asp1 <- read_stars("data/SRTM/aspect_1arc_utm18n.tif")
asp <- asp1 %>% st_warp(agb.ras, method="bilinear", use_gdal=TRUE)
asp %>% saveRDS('results/R_out/aspect_1arc_resamptoAGBbl.rds')

# Load rasters as stars
slp <- readRDS('results/R_out/slope_1arc_resamptoAGBbl.rds')
asp <- readRDS('results/R_out/aspect_1arc_resamptoAGBbl.rds')
g0 <- read_stars("results/g0nu_HV/g0nu_2018_HV.tif")

# Create DF of slope, aspect, and g0 -------------------------------------------------------
cell_numbers = g0
cell_numbers[[1]][] = 1:length(cell_numbers[[1]])
df <- 
  tibble(cell=as.vector(cell_numbers[[1]]),
         g0=as.vector(g0[[1]]) %>% round(4), 
         agb=as.vector(agb.ras[[1]]) %>% round(1),
         slope=as.vector(slp[[1]]) %>% round(1), 
         aspect=as.vector(asp[[1]]) %>% round(1)) %>% 
  filter(!is.na(g0)) 
df %>% saveRDS("results/R_out/df_srtm_g0agb.rds")
df  <- readRDS("results/R_out/df_srtm_g0agb.rds")

cor(df$agb, df$aspect, use="complete.obs")
cor(df$agb, df$slope)

# Compare slope to g0 ----------------------------------------------------------------------
df_o1 <- df %>% filter(g0>=1)
df_o1 %>% saveRDS("results/R_out/df_srtm_agb_g0_over1.rds")
df_u1 <- df %>% filter(g0<1)
df_u1 %>% saveRDS("results/R_out/df_srtm_agb_g0_under1.rds")
df_op4 <- df %>% filter(g0>0.4)
df_op4 %>% saveRDS("results/R_out/df_slope_v_g0_overpt4.rds")
df_op3 <- df %>% filter(g0>0.3)
df_op3 %>% saveRDS("results/R_out/df_slope_v_g0_overpt3.rds")
df_p3_1 <- df_u1 %>% filter(g0>0.3)
df_p3_1 %>% saveRDS("results/R_out/df_slope_v_g0_pt3to1.rds")
df_p2_3 <- df %>% filter(g0>0.2 & g0<=0.3)
df_p2_3 %>% saveRDS("results/R_out/df_slope_v_g0_pt2topt3.rds")
df_up2 <- df %>% filter(g0<=0.2)
df_up2 %>% saveRDS("results/R_out/df_slope_v_g0_underpt2.rds")
df_p1_2 <- df %>% filter(g0>0.1 & g0<=0.2)
df_p1_2 %>% saveRDS("results/R_out/df_slope_v_g0_pt1topt2.rds")
df_up1 <- df %>% filter(g0<=0.1)
df_up1 %>% saveRDS("results/R_out/df_slope_v_g0_underpt1.rds")

cor(df_o1$g0, df_o1$slope)
cor(df_op4$g0, df_op4$slope, use="complete.obs")
cor(df_p3$g0, df_p3$slope)
cor(df_p3_1$g0, df_p3_1$slope)
cor(df_p2_3$g0, df_p2_3$slope)
cor(df_up2$g0, df_up2$slope)
cor(df_p1_3$g0, df_p1_3$slope)
cor(df_up1$g0, df_up1$slope)

# Correlations between backscatter and aspect ----
cor(df_o1$g0, df_o1$aspect, use="complete.obs")
cor(df_op4$g0, df_op4$aspect, use="complete.obs")
cor(df_op3$g0, df_op3$aspect, use="complete.obs")
cor(df_p3_1$g0, df_p3_1$aspect, use="complete.obs")
cor(df_p2_3$g0, df_p2_3$aspect, use="complete.obs")
cor(df_up2$g0, df_up2$aspect, use="complete.obs")
cor(df_p1_2$g0, df_p1_2$aspect, use="complete.obs")
cor(df_up1$g0, df_up1$aspect, use="complete.obs")

# Look at correlations by slope ranges ----
df_sl_o30 <- df %>% filter(slope>=30)
df_sl_o30 %>% saveRDS("results/R_out/df_srtm_agbg0_sl_over30.rds")
cor(df_sl_o30$g0, df_sl_o30$slope)
length(df_sl_o30[[1]])
cor(df_sl_o30$g0, df_sl_o30$aspect, use="complete.obs")

df_sl_20_30 <- df %>% filter(slope<30 & slope>=20)
df_sl_20_30 %>% saveRDS("results/R_out/df_srtm_agbg0_sl_over30.rds")
cor(df_sl_20_30$g0, df_sl_20_30$slope)
length(df_sl_20_30[[1]])
cor(df_sl_20_30$g0, df_sl_20_30$aspect, use="complete.obs")

# Scatterplot
p_p3 <- ggplot(df_up1, aes(x=slope, y=g0)) + geom_point(alpha=0.05) +
  labs(y = expression(paste("Radar backscatter, ",sigma['HV']^0," (m"^2, "/m"^2, ")")), 
       x = expression(paste("Slope")))
(p_p3 <- p_p3 + geom_smooth(method="lm", se=TRUE, fullrange=TRUE, level=0.95, col='black'))
ggsave('figures/qc_plots/scatter_slope_g0_underpt1.png')

p_p2_3 <- ggplot(df_p2_3, aes(x=slope, y=g0)) + geom_point(alpha=0.05) +
  labs(y = expression(paste("Radar backscatter, ",sigma['HV']^0," (m"^2, "/m"^2, ")")), 
       x = expression(paste("Slope")))
(p_p2_3 <- p_p2_3 + geom_smooth(method="lm", se=TRUE, fullrange=TRUE, level=0.95, col='black'))

# Scatterplot
p <- ggplot(df, aes(x=slope, y=g0)) + geom_point(alpha=0.05) +
  labs(y = expression(paste("Radar backscatter, ",sigma['HV']^0," (m"^2, "/m"^2, ")")), 
     x = expression(paste("Slope")))
(p <- p + geom_smooth(method="lm", se=TRUE, fullrange=TRUE, level=0.95, col='black'))

# ----
# bounding box to crop for testing
test_bb <- st_bbox(c(xmin=-72.3, xmax=-72.23, ymin=18.32, ymax=18.38)) %>% 
  st_as_sfc()
st_crs(test_bb) <- 4326

# Look at map (small area)
tmap_mode("view")
tm_shape(g0t[test_bb]) + tm_raster() +
  tm_shape(slp[test_bb]) + tm_raster()

lm_opt3_slp <- lm(g0 ~ slope, data=df)
cov2cor(vcov(lm_opt3_slp))[1,2]

lm_opt3_slp <- lm(g0t[[1]] ~ slpt[[1]], na.exclude=TRUE)

# Use linear regression with step functions?




# Look at correlation between slope and backscatter > 0.3
agb.ras[agb.ras > 300] <- NA


# Compare slope/aspect to AGB to determine whether there is a correlation.
# Scatterplot

plot(slp, agb.ras, pch=5, col=rgb(red=0.2, green=0.2, blue=1.0, alpha=0.1))

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
br <- brick(agb.ras, slp)
names(br) <- c('AGB', 'Slope')
br
ggplot(br, aes(x=Slope, y=AGB)) +
  geom_point(alpha=0.2) 

cor(chunk3, agb.ras, method='spearman')

#---- 
# Mask out AGB in high slope areas and look at value distributions
agb.ras.m <- agb.ras
agb.ras.m[slp > 30] <- NA
br <- brick(agb.ras, agb.ras.m)
names(br) <- c('AGB', 'AGB slope sub30')
agb.ras.m[slp > 20] <- NA
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
