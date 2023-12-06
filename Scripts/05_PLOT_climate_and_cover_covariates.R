# HEADER -------------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: Plot unprocessed climate and land cover covariates
#
# Script Description: First create national boundary files to be 
# used for plots. Then source and neatly plot LCM woodland cover,
# as well as the five climate covariates for the UK.
# Can be used in SI if needed.

# LOAD LIBRARIES & INSTALL PACKAGES ----------------------

# Change  library to C: (R: doesn't have enough space for packages):
.libPaths("R:/rsrch/cb751/lab/Charles/R/PackageLibrary")

# Load packages
library(terra)
library(sf)
library(tidyverse)
library(ggplot2)

# DATA FILES -----------------------------------------

### BNG CRS

# Download BNG WKT string
download.file(url = "https://epsg.io/27700.wkt2?download=1",
              destfile = "../Data/Spatial_data/Boundaries_and_CRS/bng.prj")

# Read in BNG CRS
bng <- st_crs("../Data/Spatial_data/Boundaries_and_CRS/bng.prj")$wkt

### Land cover map files

# Land Cover Map of woodland in 1990
LCM1990file <- "../Data/Spatial_data/LCM/1990/UK1990all_1km.tif"

# Land Cover Map of woodland in 1990
LCM2015file <- "../Data/Spatial_data/LCM/2015/UK2015all_1km.tif"

### Climate covariate folders
GDD5dir <- "../Data/Spatial_data/Annual_climate_covar/GDD5"
WMINdir <- "../Data/Spatial_data/Annual_climate_covar/wmin"
soilMdir <- "../Data/Spatial_data/Annual_climate_covar/soilM"
tasCVdir <- "../Data/Spatial_data/Annual_climate_covar/tasCV"
RAINdir <- "../Data/Spatial_data/Annual_climate_covar/RAIN"

### CODE --------------------------------------------

# NATIONAL BOUNDARIES -------------------------------
# (UK, Ireland, and Isle of Man)

# Create data frame with name of country, and GADM code
boundaries <- data.frame(country  = c( "UK",  "Ireland", "IsleOfMan"),
                         GADMcode = c( "GBR", "IRL",     "IMN"))

# Loop through each country
for (i in 1:NROW(boundaries)) {
  
  # Assign year and resolution for row i
  country <- boundaries[i,"country"]
  GADMcode <- boundaries[i,"GADMcode"]
  
  # Download country boundary
  boundary <- geodata::gadm(country = GADMcode, level = 0,
                            path = paste0("../Data/Spatial_data/Boundaries_and_CRS/",
                                          country ))
  
  # Reproject country to BNG
  boundary <- project(boundary, bng)
  
  # Assign boundary to corresponding country
  assign(country, boundary)
  
  # If file not already saved, then write to file
  if(!file.exists(paste0("../Data/Spatial_data/Boundaries_and_CRS/",
                         country, "/", GADMcode, ".shp" ))) {
  
  # Save country as vector
  writeVector(boundary, 
              filename = paste0("../Data/Spatial_data/Boundaries_and_CRS/",
                                country ),
              layer = GADMcode,
              filetype = "ESRI Shapefile")
  }
}

# PLOT WOODLAND COVER -----------------------------------

# Read in 1990 and 2015 LCM
LCM1990r <- rast(LCM1990file)
LCM2015r <- rast(LCM2015file)

# Format to data frames and add year columns
LCM1990df <- as.data.frame(LCM1990r, xy = TRUE) %>%
  add_column("year" = "1990") 
LCM2015df <- as.data.frame(LCM2015r, xy = TRUE) %>%
  add_column("year" = "2015")

# Join 1990 and 2015 LCM data frames together
LCMdf <- rbind(LCM1990df, LCM2015df) %>% 
  gather(., type, cover, BF:W) %>%
  na.omit(.)

# Convert 'type' and 'year' columns to factors
LCMdf$year <- as.factor(LCMdf$year)
LCMdf$type <- as.factor(LCMdf$type)

# Specify out labels
type.labs <- c("Broadleaf", "Coniferous", "All woodland")
names(type.labs) <- c("BF", "CF", "W")

# Plot - first set plot aesthetics 
LCMmap <- ggplot(data = LCMdf, aes(x = x, y = y, 
                                   colour = cover * 100,
                                   fill = cover * 100)) +
  
  # Set equal coordinates
  coord_equal() +
  
  # Country boundary polygons
  geom_sf(data = st_as_sf(UK), fill = "grey90", colour = NA, inherit.aes = FALSE) +
  geom_sf(data = st_as_sf(Ireland), fill = "grey90", colour = NA, inherit.aes = FALSE) +
  geom_sf(data = st_as_sf(IsleOfMan), fill = "grey90", colour = NA, inherit.aes = FALSE) +
  
  # Land cover raster
  geom_tile() +
  facet_grid(vars(year),
             vars(type),
             switch = "y",
             labeller = labeller(type = type.labs)) +
  scale_colour_gradient2(midpoint = 75,
                         low = "grey80",
                         mid = "#238b45",
                         high = "#00441b",
                         guide = NULL) +
  scale_fill_gradient2( "Proportion cover  ",
                        midpoint = 75,
                        low = "grey80",
                        mid = "#238b45",
                        high = "#00441b",
                        guide = guide_colourbar(direction = "horizontal",
                                                ticks = TRUE,
                                                draw.ulim = FALSE,
                                                draw.llim = FALSE,
                                                title.position = "left",
                                                label.position = "bottom",
                                                label.theme = element_text(size = 12),
                                                title.theme = element_text(size = 12),
                                                barwidth = unit(8, "lines"),
                                                barheight = unit(0.5, "lines")),
                        labels = c("0%", "50%", "100%"),
                        breaks = c(0, 50, 100),
                        limits = c(0, 100) ) +
  
  # Set theme parameters
  theme_void() +
  theme(plot.background = element_rect(fill = "white",
                                       colour = "white"),
        legend.position = c(0.5, -0.01),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        panel.spacing.x = unit(-1, "lines"),
        panel.spacing.y = unit(-2, "lines"),
        plot.margin = margin(2, 2, 2, 2, "lines"))
  
# Save to .png file
ggsave(filename = paste0("../Writing/Plots/", "Woodland_LCM_cover.png"),
       LCMmap,
       dpi = 600,
       units = "px", width = 6000, height = 5000)
  

# PLOT WOODLAND COVER CHANGE --------------------------------

# Find change between 1990 and 2015
LCMchangeR <- LCM2015r - LCM1990r

# Aggregate to 10km for plot
LCMchangeR <- terra::aggregate(LCMchangeR,
                               fact = 5,
                               sum, na.rm = TRUE) / 5^2

# Format to data frame
LCMchangeDF <- as.data.frame(LCMchangeR, xy = TRUE) 

# Switch to long format to plot, an remove NA values
LCMchangeDF <- gather(LCMchangeDF, type, change, BF:W) %>%
  na.omit(.)

# Convert 'type' columns to factor
LCMchangeDF$type <- as.factor(LCMchangeDF$type)

# Specify out labels
type.labs <- c("Broadleaf", "Coniferous", "All woodland")
names(type.labs) <- c("BF", "CF", "W")

# Plot - first set plot aesthetics 
LCMchangeMap <- ggplot(data = LCMchangeDF, aes(x = x, y = y,
                                               colour = change * 100,
                                               fill = change * 100)) +
  
  # Set equal coordinates
  coord_equal() +
  
  # Country boundary polygons
  geom_sf(data = st_as_sf(UK), fill = "grey90", colour = NA, inherit.aes = FALSE) +
  geom_sf(data = st_as_sf(Ireland), fill = "grey90", colour = NA, inherit.aes = FALSE) +
  geom_sf(data = st_as_sf(IsleOfMan), fill = "grey90", colour = NA, inherit.aes = FALSE) +
  
  # Land cover raster
  geom_tile() +
  facet_wrap("type",
             labeller = labeller(type = type.labs)) +
  scale_colour_gradient2(midpoint = 0,
                         low = "#d73027",
                         mid = "#fee090",
                         high = "#4575b4",
                         guide = NULL) +
  scale_fill_gradient2( "Percentage cell change 1990-2015  ",
                        midpoint = 0,
                        low = "#d73027",
                        mid = "#fee090",
                        high = "#4575b4",
                        guide = guide_colourbar(direction = "horizontal",
                                                ticks = TRUE,
                                                draw.ulim = FALSE,
                                                draw.llim = FALSE,
                                                title.position = "left",
                                                label.position = "bottom",
                                                label.theme = element_text(size = 12),
                                                title.theme = element_text(size = 12),
                                                barwidth = unit(8, "lines"),
                                                barheight = unit(0.5, "lines")) ) +
  
  # Set theme parameters
  theme_void() +
  theme(plot.background = element_rect(fill = "white",
                                       colour = "white"),
        legend.position = c(0.5, -0.01),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14, face = "bold"),
        panel.spacing.x = unit(-1, "lines"),
        panel.spacing.y = unit(-2, "lines"),
        plot.margin = margin(2, 2, 2, 2, "lines"))

# Save to .png file
ggsave(filename = paste0("../Writing/Plots/", "Woodland_LCM_change.png"),
       LCMchangeMap,
       dpi = 600,
       units = "px", width = 6000, height = 4000)

# PLOT CLIMATE COVARIATES ---------------------------------

### Read in climate covariate data

GDD5df <- list.files(GDD5dir, # List .tif files in climate directory
                    full.names = TRUE, 
                    pattern = "\\.tif$") %>%
  rast(.) %>% # create spatRast of each .tif
  mean(.) %>% # Find the mean of every annual spatRast
  as.data.frame(., xy = TRUE) # Convert to spatial data frame

WMINdf <- list.files(WMINdir, # List .tif files in climate directory
                     full.names = TRUE, 
                     pattern = "\\.tif$") %>%
  rast(.) %>% # create spatRast of each .tif
  mean(.) %>% # Find the mean of every annual spatRast
  as.data.frame(., xy = TRUE) # Convert to spatial data frame

tasCVdf <- list.files(tasCVdir, # List .tif files in climate directory
                     full.names = TRUE, 
                     pattern = "\\.tif$") %>%
  rast(.) %>% # create spatRast of each .tif
  mean(.) %>% # Find the mean of every annual spatRast
  as.data.frame(., xy = TRUE) # Convert to spatial data frame

soilMdf <- list.files(soilMdir, # List .tif files in climate directory
                      full.names = TRUE, 
                      pattern = "\\.tif$") %>%
  rast(.) %>% # create spatRast of each .tif
  mean(.) %>% # Find the mean of every annual spatRast
  as.data.frame(., xy = TRUE) # Convert to spatial data frame

RAINdf <- list.files(RAINdir, # List .tif files in climate directory
                      full.names = TRUE, 
                      pattern = "\\.tif$") %>%
  rast(.) %>% # create spatRast of each .tif
  mean(.) %>% # Find the mean of every annual spatRast
  as.data.frame(., xy = TRUE) # Convert to spatial data frame

### Plot climate covariates seperately [bit verbose, potentially more succinct code possible]

# Plot GDD5
GDD5map <- ggplot(data = GDD5df, aes(x = x, y = y, 
                                     colour = mean, fill = mean)) +
  
  # Set equal coordinates
  coord_equal() +
  
  # Country boundary polygons
  geom_sf(data = st_as_sf(UK), fill = "grey90", colour = NA, inherit.aes = FALSE) +
  geom_sf(data = st_as_sf(Ireland), fill = "grey90", colour = NA, inherit.aes = FALSE) +
  geom_sf(data = st_as_sf(IsleOfMan), fill = "grey90", colour = NA, inherit.aes = FALSE) +
  
  # Land cover raster
  geom_tile() +
  scale_colour_gradient2(midpoint = 1500,
                         low = "#313695",
                         mid = "#fee090",
                         high = "#a50026",
                         guide = NULL) +
  scale_fill_gradient2( "",
                        midpoint = 1500,
                        low = "#313695",
                        mid = "#fee090",
                        high = "#a50026",
                        guide = guide_colourbar(ticks = TRUE,
                                                draw.ulim = FALSE,
                                                draw.llim = FALSE,
                                                title.position = "top",
                                                label.position = "right",
                                                label.theme = element_text(size = 12),
                                                title.theme = element_text(size = 12),
                                                barwidth = unit(2, "lines"),
                                                barheight = unit(8, "lines"))) +
  
  # Add title
  annotate(geom="text", 
           x=80000, y=1100000, 
           label = "Growing degree days\nabove 5\u00B0C\n[1980-2021 mean]", 
           size = 6) +
  
  # Set theme parameters
  theme_void() +
  theme(plot.background = element_rect(fill = "white",
                                       colour = "white"),
        legend.position = c(0.9,0.6),
        plot.margin = margin(-4,0,-4,0, "lines"))

# Plot WMIN
WMINmap <- ggplot(data = WMINdf, aes(x = x, y = y,
                                     colour = mean, fill = mean)) +
  
  # Set equal coordinates
  coord_equal() +
  
  # Country boundary polygons
  geom_sf(data = st_as_sf(UK), fill = "grey90", colour = NA, inherit.aes = FALSE) +
  geom_sf(data = st_as_sf(Ireland), fill = "grey90", colour = NA, inherit.aes = FALSE) +
  geom_sf(data = st_as_sf(IsleOfMan), fill = "grey90", colour = NA, inherit.aes = FALSE) +
  
  # Land cover raster
  geom_tile() +
  scale_colour_gradient2(midpoint = 0,
                         low = "#313695",
                         mid = "#fee090",
                         high = "#a50026",
                         guide = NULL) +
  scale_fill_gradient2( "",
                        midpoint = 0,
                        low = "#313695",
                        mid = "#fee090",
                        high = "#a50026",
                        guide = guide_colourbar(ticks = TRUE,
                                                draw.ulim = FALSE,
                                                draw.llim = FALSE,
                                                title.position = "top",
                                                label.position = "right",
                                                label.theme = element_text(size = 12),
                                                title.theme = element_text(size = 12),
                                                barwidth = unit(2, "lines"),
                                                barheight = unit(8, "lines"))) +
  
  # Add title
  annotate(geom="text", 
           x=80000, y=1100000, 
           label = "Mean minimum temperature\nof coldest month (\u00B0C)\n[1980-2021 mean]", 
           size = 6) +
  
  # Set theme parameters
  theme_void() +
  theme(plot.background = element_rect(fill = "white",
                                       colour = "white"),
        legend.position = c(0.9,0.6),
        plot.margin = margin(-4,0,-4,0, "lines"))

# Plot tasCV
tasCVmap <- ggplot(data = tasCVdf, aes(x = x, y = y, 
                                       colour = mean, fill = mean)) +
  
  # Set equal coordinates
  coord_equal() +
  
  # Country boundary polygons
  geom_sf(data = st_as_sf(UK), fill = "grey90", colour = NA, inherit.aes = FALSE) +
  geom_sf(data = st_as_sf(Ireland), fill = "grey90", colour = NA, inherit.aes = FALSE) +
  geom_sf(data = st_as_sf(IsleOfMan), fill = "grey90", colour = NA, inherit.aes = FALSE) +
  
  # Land cover raster
  geom_tile() +
  scale_colour_gradient2(midpoint = 0.016,
                         low = "#313695",
                         mid = "#fee090",
                         high = "#a50026",
                         guide = NULL) +
  scale_fill_gradient2( "",
                        midpoint = 0.016,
                        low = "#313695",
                        mid = "#fee090",
                        high = "#a50026",
                        guide = guide_colourbar(ticks = TRUE,
                                                draw.ulim = FALSE,
                                                draw.llim = FALSE,
                                                title.position = "top",
                                                label.position = "right",
                                                label.theme = element_text(size = 12),
                                                title.theme = element_text(size = 12),
                                                barwidth = unit(2, "lines"),
                                                barheight = unit(8, "lines"))) +
  
  # Add title
  annotate(geom="text", 
           x=80000, y=1100000, 
           label = "Temperature seasonality\n[1980-2021 mean]", 
           size = 6) +

  # Set theme parameters
  theme_void() +
  theme(plot.background = element_rect(fill = "white",
                                       colour = "white"),
        legend.position = c(0.9,0.6),
        plot.margin = margin(-4,0,-4,0, "lines"))

# Plot soil moisture
soilMmap <- ggplot(data = soilMdf[soilMdf$mean != 1,], aes(x = x, y = y,
                                                           colour = mean, fill = mean)) +
  
  # Set equal coordinates
  coord_equal() +
  
  # Country boundary polygons
  geom_sf(data = st_as_sf(UK), fill = "grey90", colour = NA, inherit.aes = FALSE) +
  geom_sf(data = st_as_sf(Ireland), fill = "grey90", colour = NA, inherit.aes = FALSE) +
  geom_sf(data = st_as_sf(IsleOfMan), fill = "grey90", colour = NA, inherit.aes = FALSE) +
  
  # Land cover raster excluding '1' values (lakes)
  geom_tile() +
  scale_colour_gradient2(midpoint = 0.2,
                         low = "#a50026",
                         mid = "#fee090",
                         high = "#313695",
                         guide = NULL) +
  scale_fill_gradient2(midpoint = 0.2,
                       low = "#a50026",
                       mid = "#fee090",
                       high = "#313695",
                       guide = guide_colourbar(ticks = TRUE,
                                                draw.ulim = FALSE,
                                                draw.llim = FALSE,
                                                title.position = "top",
                                                label.position = "right",
                                                label.theme = element_text(size = 12),
                                                title.theme = element_text(size = 12),
                                                barwidth = unit(2, "lines"),
                                                barheight = unit(8, "lines")),
                        limits = c(0.05, 0.55)) +

  # Add raster '1' values (lakes)
  geom_tile(data = soilMdf[soilMdf$mean == 1,], aes(x = x, y = y),
            colour = "#313695", fill = "#313695") +
  
  # Add title
  annotate(geom="text", 
           x=80000, y=1100000, 
           label="Soil Moisture\n[1981-2010 mean]",
           size = 6) +
  
  # Set theme parameters
  theme_void() +
  theme(plot.background = element_rect(fill = "white",
                                       colour = "white"),
        legend.position = c(0.9,0.6),
        plot.margin = margin(-4,0,-4,0, "lines"))

# Plot total precipitation
RAINmap <- ggplot(data = RAINdf, aes(x = x, y = y,
                                     colour = mean, fill = mean)) +
  
  # Set equal coordinates
  coord_equal() +
  
  # Country boundary polygons
  geom_sf(data = st_as_sf(UK), fill = "grey90", colour = NA, inherit.aes = FALSE) +
  geom_sf(data = st_as_sf(Ireland), fill = "grey90", colour = NA, inherit.aes = FALSE) +
  geom_sf(data = st_as_sf(IsleOfMan), fill = "grey90", colour = NA, inherit.aes = FALSE) +
  
  # Land cover raster
  geom_tile() +
  scale_colour_gradient2( "",
                        midpoint = 1500,
                        low = "#a50026",
                        mid = "#fee090",
                        high = "#313695",
                        guide = NULL) +
  scale_fill_gradient2( "",
                        midpoint = 1500,
                        low = "#a50026",
                        mid = "#fee090",
                        high = "#313695",
                        guide = guide_colourbar(ticks = TRUE,
                                                draw.ulim = FALSE,
                                                draw.llim = FALSE,
                                                title.position = "top",
                                                label.position = "right",
                                                label.theme = element_text(size = 12),
                                                title.theme = element_text(size = 12),
                                                barwidth = unit(2, "lines"),
                                                barheight = unit(8, "lines"))) +
  
  # Add title
  annotate(geom="text", 
           x=80000, y=1100000, 
           label="Total annual\nprecipitation (mm)\n[1980-2021 mean]",
           size = 6) +
  
  # Set theme parameters
  theme_void() +
  theme(plot.background = element_rect(fill = "white",
                                       colour = "white"),
        legend.position = c(0.9,0.6),
        plot.margin = margin(-4,0,-4,0, "lines"))

# Aggregate all climate covariate maps together
allClimateMap <- gridExtra::arrangeGrob(GDD5map, WMINmap, tasCVmap, soilMmap, RAINmap,
                                        nrow = 2, ncol = 3)

# Save to .png file
ggsave(filename = paste0("../Writing/Plots/", "Climate_covariate_means.png"),
       allClimateMap,
       dpi = 600,
       units = "px", width = 8000, height = 7000)
