# HEADER --------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: Plot descriptive statistics plots for cover and connectivity
#
# Script Description: Descriptive statistics plots for cover 
# and connectivity for each period, change between them, and a correlation
# plot for all variables

# LOAD LIBRARIES & INSTALL PACKAGES -----------------

# Change  library to C: (R: doesn't have enough space for packages):
.libPaths("R:/rsrch/cb751/lab/Charles/R/PackageLibrary")

# Load packages
library(tidyverse)
library(terra)
library(sf)
library(wesanderson)
library(corrplot)

# LOAD SPATIAL DATA -------------------------------------

# List all covariate SpatRasters
for (i in list.files("../Data/Spatial_data/DataForInlabru/spatRaster",
                     pattern =  "\\.tif$")) {
  
  assign(gsub(".tif", "", i),
         rast(paste0("../Data/Spatial_data/DataForInlabru/spatRaster/",
                     i))) 
}

# Boundaries
UK <- vect("../Data/Spatial_data/Boundaries_and_CRS/UK/GBR.shp")
Ireland <- vect("../Data/Spatial_data/Boundaries_and_CRS/Ireland/IRL.shp")
IsleOfMan <- vect("../Data/Spatial_data/Boundaries_and_CRS/IsleOfMan/IMN.shp")

# SET PARAMETERS ------------------------------------

# Set covariate objects
cov_R <- list(coverBF_scaled, coverCF_scaled, connW_scaled,
              coverBF_connW, coverCF_connW,
              GDD5_grp, WMIN_grp, tasCV_grp, soilM_grp, RAIN_grp)

# Set covariate object plot names
covNames <- c("Broadleaf Cover", "Coniferous Cover", "Connectivity",
              "Broadleaf-Connectivity interaction",
              "Coniferous-Connectivity interaction",
              "Growing degree days (GDD5)", "Minimum winter tempterature (WMIN)",
              "Temperature seasonality (tasCV)", "Soil moisture (soilM)",
              "Total annual rainfall (RAIN)")

# PROCESS DATA -------------------------------------

### CREATE STANDARDISED RASTERS FOR CORRELATION PLOT

# Separate spatRasts for each time period
cov_R_1990 <- lapply(cov_R, function(x) { x[[1]] }) %>% rast(.)
cov_R_2015 <- lapply(cov_R, function(x) { x[[2]] }) %>% rast(.)

# Rename 
names(cov_R_1990) <- names(cov_R_2015) <- covNames

### CREATE STANDARDISED RASTERS FOR COVER-CONNECTIVITY PLOTS

# Download BNG WKT string
download.file(url = "https://epsg.io/27700.wkt2?download=1",
              destfile = "../Data/Spatial_data/Boundaries_and_CRS/bng.prj")

# Read in BNG CRS
bng <- sf::st_crs("../Data/Spatial_data/Boundaries_and_CRS/bng.prj")$wkt

# Project back to bng in m (not km)
coverBF <- project(coverBF, bng)
coverCF <- project(coverCF, bng)
connW <- project(connW, bng)

# Standardise names
names(coverBF) <- names(coverCF) <- names(connW) <- c("year1990", "year2015")

# CORRELATION PLOT -------------------------------------------------

# 1990 correlation plot
cor1990 <- cov_R_1990 %>%
  as.data.frame %>%
  cor
png(filename = '../Writing/Plots/corr_all_1990.png',
    width = 30, height = 30, units = "cm", res = 300)
corrplot(cor1990, type = "upper", order = "original", 
         tl.col = "black", addCoef.col="black" ,tl.srt = 45)
dev.off()

# 2015 correlation plot
cor2015 <- cov_R_2015 %>%
  as.data.frame %>%
  cor
png(filename = '../Writing/Plots/corr_all_2015.png',
    width = 30, height = 30, units = "cm", res = 300)
corrplot(cor2015, type = "upper", order = "original" , 
         tl.col = "black",addCoef.col="black" ,tl.srt = 45)
dev.off()

# Combine into single plot
corBoth <- cowplot::ggdraw(clip = "on") +
  cowplot::draw_image('../Writing/Plots/corr_all_1990.png', -0.1, 0.4, 1.25, 0.625) +
  cowplot::draw_image('../Writing/Plots/corr_all_2015.png', -0.1, -0.1, 1.25, 0.625) +
  cowplot::draw_label("(a) 1990", 0.1, 0.95, size = 22) +
  cowplot::draw_label("(b) 2015", 0.1, 0.45, size = 22)

# Save
ggsave(filename = paste0("../Writing/Plots/", "allCorr.png"),
       corBoth,
       dpi = 600,
       units = "px", width = 4000, height = 8400)  

# Remove individual plots as not needed
unlink('../Writing/Plots/corr_all_1990.png')
unlink('../Writing/Plots/corr_all_2015.png') 

# 1KM STATIC PLOTS -----------------------------------------

# Calculate change from 1990 to 2015
coverBF[["change"]] <-  coverBF[["year2015"]] - coverBF[["year1990"]] 
coverCF[["change"]] <-  coverCF[["year2015"]] - coverCF[["year1990"]] 
connW[["change"]]  <-  connW[["year2015"]] - connW[["year1990"]] 

# Format to data frames and add type columns
coverBF_df <- as.data.frame(coverBF, xy = TRUE) %>% add_column(type = "Broadleaf")
coverCF_df <- as.data.frame(coverCF, xy = TRUE) %>% add_column(type = "Coniferous")
connW_df <- as.data.frame(connW, xy = TRUE) %>% add_column(type = "Connectivity")

# Join 1990 and 2015 LCM data frames together
plot_df <- rbind(coverBF_df, coverCF_df, connW_df)

# 10KM CHANGE PLOTS

# Aggregate to 10km for plot (mean)
coverBF_10k <- terra::aggregate(coverBF,
                                fact = 10,
                                sum, na.rm = TRUE) / 10^2
coverCF_10k <- terra::aggregate(coverCF,
                                fact = 10,
                                sum, na.rm = TRUE) / 10^2
connW_10k <- terra::aggregate(connW,
                              fact = 10,
                              sum, na.rm = TRUE) / 10^2

# Format to data frames and add type columns
coverBF_10k_df <- as.data.frame(coverBF_10k, xy = TRUE) %>% add_column(type = "Broadleaf")
coverCF_10k_df <- as.data.frame(coverCF_10k, xy = TRUE) %>% add_column(type = "Coniferous")
connW_10k_df <- as.data.frame(connW_10k, xy = TRUE) %>% add_column(type = "Connectivity")

# Join 1990 and 2015 LCM data frames together
plot_10k_df <- rbind(coverBF_10k_df, coverCF_10k_df, connW_10k_df)

# PLOT CORRELATION ------------------------------------

# Create 2015 data frame from plot data frame
corr_df <- pivot_wider( plot_df,
                         id_cols = c(x, y),
                         names_from = type,
                         values_from = c(year2015, change)) %>%
  # Create 2015 and change columns for total woodland
  mutate(year2015_Woodland = year2015_Broadleaf + year2015_Coniferous) %>%
  mutate(change_Woodland = change_Broadleaf + change_Coniferous)


### PLOT 2015 VALUES

# Create data frame
covConnCorr <-  pivot_longer(corr_df, # Gather into long format
                             cols = c(year2015_Woodland,
                                      year2015_Coniferous,
                                      year2015_Broadleaf)) %>%
  
  # Plot
  ggplot(aes(x = value ,
             y = year2015_Connectivity )) +
  
  # Add points
  geom_point(size = 0.5, shape = 16, alpha = 0.3) +
  
  # Facet wrap and strip labels
  facet_wrap( ~name ,
              labeller = as_labeller(
                c(year2015_Woodland = "All Woodland",
                  year2015_Coniferous = "Coniferous",
                  year2015_Broadleaf = "Broadleaf"))) +
  
  # Axis labels and theme adjustments
  labs(x = "Proportion woodland cover in 2015",
       y = "Connectivity (Amps)") + 
  theme_minimal() +
  theme(axis.title = element_text(size = 16), 
        axis.text = element_text(size = 14),  
        strip.text = element_text(size = 14)) 

# Save to .png file
ggsave(filename = paste0("../Writing/Plots/",
                         "Cover_connectivity_corellation.png"),
       covConnCorr,
       dpi = 600,
       units = "px", width = 8000, height = 3000)

### PLOT CHANGE VALUES

# Create data frame
covConnChangeCorr <- 
  pivot_longer(corr_df, # Gather into long format
               cols = c(change_Woodland,
                        change_Coniferous,
                        change_Broadleaf )) %>%
  # Plot
  ggplot(aes(x = value,
             y = change_Connectivity )) +
  
  # Add points
  geom_point(size = 0.5, shape = 16, alpha = 0.3) +
  
  # Facet wrap and strip labels
  facet_wrap( ~name ,
              labeller = as_labeller(
                c(change_Woodland = "All Woodland",
                  change_Coniferous = "Coniferous",
                  change_Broadleaf = "Broadleaf"))) +
  
  # Axis labels and theme
  labs(x = "Absolute proportion woodland cover change",
       y = "Change in connectivity (Amps)") + 
  theme_minimal() +
  theme(axis.title = element_text(size = 16), 
        axis.text = element_text(size = 14),  
        strip.text = element_text(size = 14)) 

# Save to .png file
ggsave(filename = paste0("../Writing/Plots/",
                         "Cover_connectivity_change_corellation.png"),
       covConnChangeCorr,
       dpi = 600,
       units = "px", width = 8000, height = 3000)

# PLOT COVER AND CONNECTIVITY ---------------------------------------

# Set shared theme parameters
sharedTheme <- theme(legend.justification = c(0.1, 0.1),
                     legend.position = c(0.1, 0.65),
                     plot.title = element_text(size = 16, hjust = 0.1, vjust = -5),
                     plot.margin = margin(-1, -1, -1, -1, "lines" ))

### CREATE STATIC PLOTS

# BROADLEAF COVER

# Plot - first set plot aesthetics 
staticBF_map <- ggplot(data = filter(plot_df, type == "Broadleaf"),
                       aes(x = x, y = y, 
                           colour = year2015 * 100,
                           fill = year2015 * 100)) +
  
  # Set equal coordinates
  coord_equal() +
  
  # Land cover raster
  geom_tile() +
  
  # Fill options
  scale_fill_gradient2( "Broadleaf cover",
                        midpoint = 75,
                        low = "grey90",
                        mid = "#238b45",
                        high = "#00441b",
                        aesthetics = c("fill", "colour"),
                        guide = guide_colourbar(direction = "vertical",
                                                ticks = TRUE,
                                                draw.ulim = FALSE,
                                                draw.llim = FALSE,
                                                title.position = "top",
                                                label.position = "left",
                                                label.theme = element_text(size = 8),
                                                title.theme = element_text(size = 10),
                                                barwidth = unit(1, "lines"),
                                                barheight = unit(5, "lines")),
                        labels = c("0%", "50%", "100%"),
                        breaks = c(0, 50, 100),
                        limits = c(0, 100) ) +
  
  # Country boundary polygons
  geom_sf(data = st_as_sf(UK), fill = NA, colour = "black", inherit.aes = FALSE) +
  geom_sf(data = st_as_sf(Ireland), fill = NA, colour = "black", inherit.aes = FALSE) +
  geom_sf(data = st_as_sf(IsleOfMan), fill = NA, colour = "black", inherit.aes = FALSE) +
  
  # Set theme parameters
  theme_void() +
  ggtitle("(a)") +
  sharedTheme

# CONIFEROUS COVER

# Plot - first set plot aesthetics 
staticCF_map <- ggplot(data = filter(plot_df, type == "Coniferous"),
                       aes(x = x, y = y,
                           colour = year2015 * 100,
                           fill = year2015 * 100)) +
  
  # Set equal coordinates
  coord_equal() +
  
  # Land cover raster
  geom_tile() +
  
  # Fill options
  scale_fill_gradient2( "Coniferous cover",
                        midpoint = 75,
                        low = "grey90",
                        mid = "#238b45",
                        high = "#00441b",
                        aesthetics = c("fill", "colour"),
                        guide = guide_colourbar(direction = "vertical",
                                                ticks = TRUE,
                                                draw.ulim = FALSE,
                                                draw.llim = FALSE,
                                                title.position = "top",
                                                label.position = "left",
                                                label.theme = element_text(size = 8),
                                                title.theme = element_text(size = 10),
                                                barwidth = unit(1, "lines"),
                                                barheight = unit(5, "lines")),
                        labels = c("0%", "50%", "100%"),
                        breaks = c(0, 50, 100),
                        limits = c(0, 100) ) +
  
  # Country boundary polygons
  geom_sf(data = st_as_sf(UK), fill = NA, colour = "black", inherit.aes = FALSE) +
  geom_sf(data = st_as_sf(Ireland), fill = NA, colour = "black", inherit.aes = FALSE) +
  geom_sf(data = st_as_sf(IsleOfMan), fill = NA, colour = "black", inherit.aes = FALSE) +
  
  # Set theme parameters
  theme_void() +
  ggtitle("(c)") +
  sharedTheme

### CONNECTIVITY

# Plot - first set plot aesthetics 
staticConn_map <- ggplot(data = filter(plot_df, type == "Connectivity"),
                       aes(x = x, y = y, 
                           colour = year2015,
                           fill = year2015)) +
  
  # Set equal coordinates
  coord_equal() +
  
  # Land cover raster
  geom_tile() +
  
  # Fill options
  scale_fill_gradient2( "Connectivity (Amps)",
                        midpoint = 150,
                        low = "grey90",
                        mid = "#238b45",
                        high = "#00441b",
                        aesthetics = c("fill", "colour"),
                        guide = guide_colourbar(direction = "vertical",
                                                ticks = TRUE,
                                                draw.ulim = FALSE,
                                                draw.llim = FALSE,
                                                title.position = "top",
                                                label.position = "left",
                                                label.theme = element_text(size = 8),
                                                title.theme = element_text(size = 10),
                                                barwidth = unit(1, "lines"),
                                                barheight = unit(5, "lines"))) +
  
  # Country boundary polygons
  geom_sf(data = st_as_sf(UK), fill = NA, colour = "black", inherit.aes = FALSE) +
  geom_sf(data = st_as_sf(Ireland), fill = NA, colour = "black", inherit.aes = FALSE) +
  geom_sf(data = st_as_sf(IsleOfMan), fill = NA, colour = "black", inherit.aes = FALSE) +
  
  # Set theme parameters
  theme_void() +
  ggtitle("(e)") +
  sharedTheme

### CREATE CHANGE PLOTS (BROADLEAF COVER, CONIFEROUS COVER, CONNECTIVITY)

# BROADLEAF COVER

# Plot - first set plot aesthetics 
changeBF_map <- ggplot(data = filter(plot_10k_df, type == "Broadleaf"),
                      aes(x = x, y = y,
                          colour = change,
                          fill = change)) +
  
  # Set equal coordinates
  coord_equal() +
  
  # Land cover raster
  geom_tile() +

  # Fill options
  scale_fill_gradient2( "Absolute proportion\nbroadleaf cover change",
                        midpoint = 0,
                        low = wes_palette("Zissou1")[5],
                        mid = "white",
                        high = wes_palette("Zissou1")[1],
                        aesthetics = c("fill", "colour"),
                        guide = guide_colourbar(direction = "vertical",
                                                ticks = TRUE,
                                                draw.ulim = FALSE,
                                                draw.llim = FALSE,
                                                title.position = "top",
                                                label.position = "left",
                                                label.theme = element_text(size = 8),
                                                title.theme = element_text(size = 10),
                                                barwidth = unit(1, "lines"),
                                                barheight = unit(5, "lines"))) +
  
  # Country boundary polygons
  geom_sf(data = sf::st_as_sf(UK), fill = NA, colour = "black", inherit.aes = FALSE) +
  geom_sf(data = sf::st_as_sf(Ireland), fill = NA, colour = "black", inherit.aes = FALSE) +
  geom_sf(data = sf::st_as_sf(IsleOfMan), fill = NA, colour = "black", inherit.aes = FALSE) +
  
  # Set theme parameters
  theme_void() +
  ggtitle("(b)") +
  sharedTheme


# CONIFEROUS COVER

# Plot - first set plot aesthetics 
changeCF_map <- ggplot(data = filter(plot_10k_df, type == "Coniferous"),
                      aes(x = x, y = y,
                          colour = change,
                          fill = change)) +
  
  # Set equal coordinates
  coord_equal() +
  
  # Land cover raster
  geom_tile() +
  
  # Fill options
  scale_fill_gradient2( "Absolute proportion\n coniferous cover change",
                        midpoint = 0,
                        low = wes_palette("Zissou1")[5],
                        mid = "white",
                        high = wes_palette("Zissou1")[1],
                        aesthetics = c("fill", "colour"),
                        guide = guide_colourbar(direction = "vertical",
                                                ticks = TRUE,
                                                draw.ulim = FALSE,
                                                draw.llim = FALSE,
                                                title.position = "top",
                                                label.position = "left",
                                                label.theme = element_text(size = 8),
                                                title.theme = element_text(size = 10),
                                                barwidth = unit(1, "lines"),
                                                barheight = unit(5, "lines"))) +
  
  # Country boundary polygons
  geom_sf(data = sf::st_as_sf(UK), fill = NA, colour = "black", inherit.aes = FALSE) +
  geom_sf(data = sf::st_as_sf(Ireland), fill = NA, colour = "black", inherit.aes = FALSE) +
  geom_sf(data = sf::st_as_sf(IsleOfMan), fill = NA, colour = "black", inherit.aes = FALSE) +
  
  # Set theme parameters
  theme_void() +
  ggtitle("(d)") +
  sharedTheme

# CONNECTIVITY

# Plot - first set plot aesthetics 
changeConn_map <- ggplot(data = filter(plot_10k_df, type == "Connectivity"),
                        aes(x = x, y = y,
                            colour = change,
                            fill = change)) +
  
  # Set equal coordinates
  coord_equal() +
  
  # Land cover raster
  geom_tile() +
  
  # Fill options
  scale_fill_gradient2( "Woodland connectivity\nchange (amps)",
                        midpoint = 0,
                        low = wes_palette("Zissou1")[5],
                        mid = "white",
                        high = wes_palette("Zissou1")[1],
                        aesthetics = c("fill", "colour"),
                        guide = guide_colourbar(direction = "vertical",
                                                ticks = TRUE,
                                                draw.ulim = FALSE,
                                                draw.llim = FALSE,
                                                title.position = "top",
                                                label.position = "left",
                                                label.theme = element_text(size = 8),
                                                title.theme = element_text(size = 10),
                                                barwidth = unit(1, "lines"),
                                                barheight = unit(5, "lines"))) +
  
  # Country boundary polygons
  geom_sf(data = sf::st_as_sf(UK), fill = NA, colour = "black", inherit.aes = FALSE) +
  geom_sf(data = sf::st_as_sf(Ireland), fill = NA, colour = "black", inherit.aes = FALSE) +
  geom_sf(data = sf::st_as_sf(IsleOfMan), fill = NA, colour = "black", inherit.aes = FALSE) +
  
  # Set theme parameters
  theme_void() +
  ggtitle("(f)") +
  sharedTheme

### COLLATE PLOTS AND SAVE

# Group into single grob
allChange <- gridExtra::arrangeGrob( staticBF_map, changeBF_map,
                                     staticCF_map, changeCF_map,
                                     staticConn_map, changeConn_map,
                                     nrow = 3, ncol = 2)

# Save to .png file
ggsave(filename = paste0("../Writing/Plots/",
                         "Cover_connectivity_change.png"),
       allChange,
       dpi = 600,
       units = "px", width = 4000, height = 8800) 
