#library(silvr)
library(readr)
library(tidyverse)
library(raster)
library(ggridges)

# Load data - Desktop
# Load data - Mac
g0_AGB <- read_csv("~/GitHub/biomass-espanola/data/plots_g0nu_AGB.csv")
# just get the two columns we care about
g0.agb <- as.data.frame(cbind(g0_AGB$AGB_ha, g0_AGB$'2018mean')) %>% 
  rename(AGB = V1, backscatter =V2)
# Linear model
ols <- lm(g0.agb$AGB ~ g0.agb$backscatter, x=TRUE, y=TRUE)

#----
agb.hispaniola <- agb
agb.samp.hispaniola <- agb.samp

agb <- raster("~/PROJECTS/Haiti_biomass/biota_out/AGB_2018_v5q_haiti.tif")
agb
NAvalue(agb) <- 0
qs.agb <- quantile(agb, probs=c(0, 0.1,0.5,0.9,0.99, 0.999, 1))
qs.agb
agb.q90s <- quantile(agb, probs=c(0.91, 0.92, 0.93, 0.94, 0.95))
View(agb.q90s)
# the distribution of values in the raster
agb.samp <- agb %>% 
  sampleRandom(100000, na.rm=TRUE) %>% 
  as.data.frame()
agb.samp <- agb.samp %>% rename(AGB='.')
agb.samp.haiti <- agb.samp

# Boxplot
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
lc.br <- brick(raster("~/PROJECTS/Haiti_biomass/LULC/Haiti2017_water_agb18v5q.tif"), 
                    raster("~/PROJECTS/Haiti_biomass/LULC/Haiti2017_urban_agb18v5q.tif"), 
                    raster("~/PROJECTS/Haiti_biomass/LULC/Haiti2017_bareland_agb18v5q.tif"), 
                    raster("~/PROJECTS/Haiti_biomass/LULC/Haiti2017_treecover_agb18v5q.tif"), 
                    raster("~/PROJECTS/Haiti_biomass/LULC/Haiti2017_grassland_agb18v5q.tif"), 
                    raster("~/PROJECTS/Haiti_biomass/LULC/Haiti2017_shrubs_agb18v5q.tif"))
names(lc.br) <- c('water', 'urban', 'bareland', 'tree cover', 'grassland', 'shrubs')
NAvalue(lc.br) <- 0
saveRDS(lc.br, file = "~/PROJECTS/Haiti_biomass/R_out/rasterbrick_AGBbyLC.rds")
readRDS(file = "~/PROJECTS/Haiti_biomass/R_out/rasterbrick_AGBbyLC.rds")

# Get selection of percentiles for each LC
lulc.qs <- data.frame(row.names=c('0%', '1%','2%','10%', '25%', '50%', '75%', '90%','98%', '99%','100%'))
for (i in seq(1,nlayers(lc.br))){
  lc <- names(lc.br[[i]])
  lulc.qs[[lc]] <- quantile(lc.br[[i]], probs=c(0, 0.01, 0.02, 0.1, 0.25, 0.5, 0.75, 0.9,0.98, 0.99, 1))
}
saveRDS(lulc.qs, file = "~/PROJECTS/Haiti_biomass/R_out/lc_quantiles.rds")
means <- cellStats(lc.br, stat='mean', na.rm=TRUE)
sds <- cellStats(lc.br, stat='sd', na.rm=TRUE)
skews <- cellStats(lc.br, stat='skew', na.rm=TRUE)
lc.stats <- cbind(means, sds, skews)
saveRDS(lulc.qs, file = "~/PROJECTS/Haiti_biomass/R_out/lc_stats.rds")

# Sample AGB by LC rasters
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
lc.sampNA <- lc.samp %>% 
  mutate(AGB = na_if(AGB, 0)) %>% 
  na.omit()
saveRDS(lc.sampNA, file = "~/PROJECTS/Haiti_biomass/R_out/lc_samp_1mil.rds")
lc.sampNA <- readRDS("~/PROJECTS/Haiti_biomass/R_out/lc_samp_1mil.rds")

lc.sampNA$Category <- ordered(lc.sampNA$Category,
                         levels = c("water", "urban", "bareland", "grassland", "shrubs", "tree.cover"))
levels(lc.sampNA$Category)

# ANOVA with multiple pair-wise comparison
res.aov <- aov(AGB ~ Category, data = lc.sampNA)
summary(res.aov)
TukeyHSD(res.aov)
pairwise.t.test(lc.sampNA$AGB, lc.sampNA$Category,
                p.adjust.method = "BH")
plot(res.aov, 1)
library(car)
leveneTest(AGB ~ Category, data = lc.sampNA)

lc.sampNA[283224,]
lc.sampNA[8052,]
lc.sampNA[9444,]
lc.samp.sub100 <- lc.sampNA[which(lc.sampNA$AGB <= 100),]
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
mu <- lc.samp.sub100 %>%
  group_by(Category) %>%
  summarise(mean = mean(AGB), med=median(AGB))

dp <- ggplot(lc.samp.sub100, aes(x = AGB, color = Category, fill = Category)) +
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

dp <- ggplot(lc.samp.sub100, aes(x = AGB, y = Category, fill = Category)) +
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
  stat_summary(fun.y=mean, geom="point", shape=5, size=4) +
  scale_fill_manual(values=custom.col) +
  coord_flip()+
  ylab("Aboveground Biomass") +
  xlab("")+
  theme_minimal()+
  theme(legend.position="none")
bp


#----
# Old
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


# Backscatter values
g0 <- raster("~/PROJECTS/Haiti_biomass/biota_out/g0nu_2018_nofilt_HV.tif")
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
