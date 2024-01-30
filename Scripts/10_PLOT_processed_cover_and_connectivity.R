# HEADER --------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: Plot processed woodland cover and connectivity
#
# Script Description: Descriptive statistics and plots for cover 
# and connecitivity change summary

# LOAD LIBRARIES & INSTALL PACKAGES -----------------

# Change  library to C: (R: doesn't have enough space for packages):
.libPaths("R:/rsrch/cb751/lab/Charles/R/PackageLibrary")

# Load packages
library(tidyverse)
library(terra)

# LOAD SPATIAL DATA -------------------------------------

# SpatRasters
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

# CREATE STANDARDISED RASTERS --------------------------

# Download BNG WKT string
download.file(url = "https://epsg.io/27700.wkt2?download=1",
              destfile = "../Data/Spatial_data/Boundaries_and_CRS/bng.prj")

# Read in BNG CRS
bng <- sf::st_crs("../Data/Spatial_data/Boundaries_and_CRS/bng.prj")$wkt

# Convert SGDFs used in model to SpatRasts, and project back to bng in m (not km)
coverBF <- project(coverBF, bng)
coverCF <- project(coverCF, bng)
connW <- project(connW, bng)

# Standardise names
names(coverBF) <- names(coverCF) <- names(connW) <- c("year1990", "year2015")

# PLOT CORELLATION ------------------------------------

### PROCESS SPATRASTERS

# Create vectors of all 1x1km cells (remove NA) for cover and connectivity
coverBF_vector <- values(coverBF, na.rm = TRUE ) %>%
  as.vector
connW_vector <- values(connW, na.rm = TRUE ) %>%
  as.vector
# All woodland or broadleaf woodland
connW_df <- connW %>%
  as.data.frame(.) %>%
  gather("Year", "Connectivity")
coverW_df <- (coverBF + coverCF) %>%
  as.data.frame(.) %>%
  gather("Year", "Cover")

# Plot
ggplot(data.frame(coverBF_vector, connW_vector),
       aes(x = coverBF_vector ,
           y = connW_vector)) +
  geom_point(color = "blue", size = 1) +
  labs(title = "Scatter Plot Example",
       x = "X-axis",
       y = "Y-axis")

# Save to .png file
ggsave(filename = paste0("../Writing/Plots/",
                         "Cover_connectivity_corellation.png"),
       allChange,
       dpi = 600,
       units = "px", width = 6000, height = 3600)



# PLOT CHANGE ---------------------------------------

### PROCESS SPATRASTERS

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

# Calculate change from 1990 to 2015
coverBF_10k[["change"]] <-  coverBF_10k[["year2015"]] - coverBF_10k[["year1990"]] 
coverCF_10k[["change"]] <-  coverCF_10k[["year2015"]] - coverCF_10k[["year1990"]] 
connW_10k[["change"]]  <-  connW_10k[["year2015"]] - connW_10k[["year1990"]] 

# CREATE INDIVIDUAL PLOTS

# Shared theme
changeTheme <- theme(
  plot.background = element_rect(fill = "white",
                                 colour = "white"),
  legend.position = c(0.5, -0.04),
  strip.text.x = element_text(size = 10, face = "bold"),
  strip.text.y = element_text(size = 10, face = "bold"),
  panel.spacing.x = unit(-3, "lines"),
  panel.spacing.y = unit(-3, "lines"),
  plot.margin = margin(0, -0.5, 2, -0.5, "lines" ))

### Broadleaf

# Plot - first set plot aesthetics 
BFchangeMap <- ggplot(data = as.data.frame(coverBF_10k, xy = TRUE),
                    aes(x = x, y = y,
                        colour = change ,
                        fill = change)) +
  
  # Set equal coordinates
  coord_equal() +
  
  # Land cover raster
  geom_tile() +
  
  # Country boundary polygons
  geom_sf(data = sf::st_as_sf(UK), fill = NA, colour = "black", inherit.aes = FALSE) +
  geom_sf(data = sf::st_as_sf(Ireland), fill = NA, colour = "black", inherit.aes = FALSE) +
  geom_sf(data = sf::st_as_sf(IsleOfMan), fill = NA, colour = "black", inherit.aes = FALSE) +
  
  # facet_wrap("type",
  #            labeller = labeller(type = type.labs)) +
  scale_colour_gradient2(midpoint = 0,
                         low = "#a50026",
                         mid = "#ffffbf",
                         high = "#313695",
                         guide = NULL) +
  scale_fill_gradient2( "Absolute proportion broadleaf\ncover change 1990-2015     ",
                        midpoint = 0,
                        low = "#a50026",
                        mid = "#ffffbf",
                        high = "#313695",
                        guide = guide_colourbar(direction = "horizontal",
                                                ticks = TRUE,
                                                draw.ulim = FALSE,
                                                draw.llim = FALSE,
                                                title.position = "top",
                                                label.position = "bottom",
                                                label.theme = element_text(size = 10),
                                                title.theme = element_text(size = 10),
                                                barwidth = unit(8, "lines"),
                                                barheight = unit(0.5, "lines")) ) +
  
  # Set theme parameters
  theme_void() +
  changeTheme

### Coniferous

# Plot - first set plot aesthetics 
CFchangeMap <- ggplot(data = as.data.frame(coverCF_10k, xy = TRUE),
                      aes(x = x, y = y,
                          colour = change ,
                          fill = change)) +
  
  # Set equal coordinates
  coord_equal() +
  
  # Land cover raster
  geom_tile() +
  
  # Country boundary polygons
  geom_sf(data = sf::st_as_sf(UK), fill = NA, colour = "black", inherit.aes = FALSE) +
  geom_sf(data = sf::st_as_sf(Ireland), fill = NA, colour = "black", inherit.aes = FALSE) +
  geom_sf(data = sf::st_as_sf(IsleOfMan), fill = NA, colour = "black", inherit.aes = FALSE) +
  
  # facet_wrap("type",
  #            labeller = labeller(type = type.labs)) +
  scale_colour_gradient2(midpoint = 0,
                         low = "#a50026",
                         mid = "#ffffbf",
                         high = "#313695",
                         guide = NULL) +
  scale_fill_gradient2( "Absolute proportion coniferous\ncover change 1990-2015     ",
                        midpoint = 0,
                        low = "#a50026",
                        mid = "#ffffbf",
                        high = "#313695",
                        guide = guide_colourbar(direction = "horizontal",
                                                ticks = TRUE,
                                                draw.ulim = FALSE,
                                                draw.llim = FALSE,
                                                title.position = "top",
                                                label.position = "bottom",
                                                label.theme = element_text(size = 10),
                                                title.theme = element_text(size = 10),
                                                barwidth = unit(8, "lines"),
                                                barheight = unit(0.5, "lines")) ) +
  
  # Set theme parameters
  theme_void() +
  changeTheme

### Connectivity

# Plot - first set plot aesthetics 
connChangeMap <- ggplot(data = as.data.frame(connW_10k, xy = TRUE),
                      aes(x = x, y = y,
                          colour = change ,
                          fill = change)) +
  
  # Set equal coordinates
  coord_equal() +
  
  # Land cover raster
  geom_tile() +
  
  # Country boundary polygons
  geom_sf(data = sf::st_as_sf(UK), fill = NA, colour = "black", inherit.aes = FALSE) +
  geom_sf(data = sf::st_as_sf(Ireland), fill = NA, colour = "black", inherit.aes = FALSE) +
  geom_sf(data = sf::st_as_sf(IsleOfMan), fill = NA, colour = "black", inherit.aes = FALSE) +
  
  # facet_wrap("type",
  #            labeller = labeller(type = type.labs)) +
  scale_colour_gradient2(midpoint = 0,
                         low = "#a50026",
                         mid = "#ffffbf",
                         high = "#313695",
                         guide = NULL) +
  scale_fill_gradient2( "Woodland connectivity change\n1990-2015 (amps) ",
                        midpoint = 0,
                        low = "#a50026",
                        mid = "#ffffbf",
                        high = "#313695",
                        guide = guide_colourbar(direction = "horizontal",
                                                ticks = TRUE,
                                                draw.ulim = FALSE,
                                                draw.llim = FALSE,
                                                title.position = "top",
                                                label.position = "bottom",
                                                label.theme = element_text(size = 10),
                                                title.theme = element_text(size = 10),
                                                barwidth = unit(8, "lines"),
                                                barheight = unit(0.5, "lines")) ) +
  
  # Set theme parameters
  theme_void() +
  changeTheme

### COLLATE PLOTS AND SAVE

# Group into single grob
allChange <- gridExtra::arrangeGrob( BFchangeMap, CFchangeMap, connChangeMap,
                       nrow = 1, ncol = 3)

# Save to .png file
ggsave(filename = paste0("../Writing/Plots/",
                         "Cover_connectivity_change.png"),
       allChange,
       dpi = 600,
       units = "px", width = 6000, height = 3600)
