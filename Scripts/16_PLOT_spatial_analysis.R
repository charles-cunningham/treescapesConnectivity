# HEADER --------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: SDM spatial analysis
#
# Script Description: Analysis of SDM inlabru model outputs to investigate 
# potential spatial priorities for targeting improving connectivity

# LOAD LIBRARIES & INSTALL PACKAGES -----------------

# Change  library to C: (R: doesn't have enough space for packages):
.libPaths("R:/rsrch/cb751/lab/Charles/R/PackageLibrary")

# Load packages
library(terra)
library(sf)
library(tidyverse)
library(RColorBrewer)

# SET PARAMETERS ------------------------------------

# DATA FILES ------------------------------------------

# Load SDM fixed effect summaries
load(file = "../Data/Species_data/SDM_fixed_effect_summaries.RData")

### Download BNG WKT string
download.file(url = "https://epsg.io/27700.wkt2?download=1",
              destfile = "../Data/Spatial_data/Boundaries_and_CRS/bng.prj")
bng <- sf::st_crs("../Data/Spatial_data/Boundaries_and_CRS/bng.prj")$wkt

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

# Load previous plots
tilePlot <- readRDS("../Data/Species_Data/tilePlotBF.RDS")

# DESCRIPTIVE POOLED STATS --------------------------------

# Group by taxa group and effect category
# (focus on just broadleaf species for spatial analysis)
group_df <- meta_df %>%
  group_by( broadleafAssociation, connectivitySig ) 

# Get basic numbers
summarise(group_df, length(species))

# GROUPED COVER AND CONNECTIVITY ANALYSIS -----------------------------------

# # Create raster for sum of every pooled group
# 
# # Loop through every group in group_df
# groupSumsR <- lapply(1:NROW(group_keys(group_df)) , function(i) {
#   
#   # Split group into separate df, get df i...
#   iGroup <- group_split(group_df)[[i]] %>%
#     dplyr::select(species, taxa) #... and select taxa and species
#   
#   # Create list of file names of occupancy plots
#   # for every row in iGroup
#   iGroupFiles <- apply(iGroup, 1, function(j) {
#     paste0("../Data/Species_data/SDMs/",
#            j["taxa"],
#            "/",
#            j["species"],
#            "/medianPred.tif")
#   })
#   
#   # Create spatRast (select second layer of each prediction [current])
#   iGroupRast <- lapply(iGroupFiles, 
#                        function(x) { rast(x, lyrs = 2) }) %>%
#     rast
#   
#   # Sum
#   iGroupSum <- sum(iGroupRast)
#   
#   # Convert projection back from km to m
#   iGroupSum <- project(iGroupSum, bng)
#   
#   # Change spatRast layer name to group key
#   names(iGroupSum) <- group_keys(group_df)[i, ] %>%
#     as.character  %>%
#     paste(., collapse = " ")
#   
#   return(iGroupSum)
#   
#   # Join all group sum spatRasts into single spatRast  
# }) %>% rast(.) 
# 
# # Calculate mean
# groupMeansR <- groupSumsR /
#   sapply(group_split(group_df), NROW)
# 
# # Save
# writeRaster(groupSumsR,
#             "../Data/Species_data/Summed_occupancy_grouped_by_effects.tif",
#             overwrite = TRUE)
# writeRaster(groupMeansR,
#             "../Data/Species_data/Mean_occupancy_grouped_by_effects.tif",
#             overwrite = TRUE)

# PLOT RICHNESS MAPS ----------------------------------------------------------

# Load spatRasters
groupSumsR <- rast("../Data/Species_data/Summed_occupancy_grouped_by_effects.tif")
groupMeansR <- rast("../Data/Species_data/Mean_occupancy_grouped_by_effects.tif")

# Create plot directory
dir.create("../Writing/Plots/GroupedOccPlots")

# Loop through each cover/connectivity group
for(i in 1:nlyr(groupMeansR)) {
  
  # Create data frame from grouped mean spatRast i
  iGroupMeans_df <- groupMeansR[[i]] %>%
    as.data.frame(., xy = TRUE)
  
  # Create ggplot
  occPlot <- ggplot(data = iGroupMeans_df,
                    aes(x = x, y = y,
                        colour = iGroupMeans_df[, 3],
                        fill = iGroupMeans_df[, 3])) +
    
    # Set equal coordinates
    coord_equal() +
    
    # Country boundary polygons
    geom_sf(data = st_as_sf(UK), fill = "grey90", colour = NA, inherit.aes = FALSE) +
    geom_sf(data = st_as_sf(Ireland), fill = "grey90", colour = NA, inherit.aes = FALSE) +
    geom_sf(data = st_as_sf(IsleOfMan), fill = "grey90", colour = NA, inherit.aes = FALSE) +
    
    # Raster
    geom_tile() +
    scale_colour_gradient(low = "#e5f5e0",
                          high = "#31a354",
                          guide = NULL) +
    scale_fill_gradient(
      "Mean relative\noccupancy probability",
      low = "#e5f5e0",
      high = "#31a354",
      #limits = c(0,0.5),
      guide = guide_colourbar(
        ticks = TRUE,
        draw.ulim = FALSE,
        draw.llim = FALSE,
        title.position = "top",
        label.position = "right",
        label.theme = element_text(size = 20),
        title.theme = element_text(size = 20, vjust = 0.5),
        barwidth = unit(2, "lines"),
        barheight = unit(8, "lines") )) +

    # Add title
    annotate(
      geom = "text",
      x = 20000,
      y = 100000,
      label = paste(NROW(group_split(group_df)[[i]]),
                    "species"),
      size = 12) +
    
    # Set theme parameters
    theme_void() +
    theme(plot.background = element_rect(fill = "white",
                                         colour = "white"),
          legend.position = c(0.3, 0.8),
          plot.margin = margin(-3, 0, -3, 0, "lines"))
  
  ggsave(filename = paste0("../Writing/Plots/GroupedOccPlots/",
                           "Occupancy_broadleaf",
                           strsplit(names(iGroupMeans_df[3]), " ")[[1]][1],
                           "_connectivity",
                           strsplit(names(iGroupMeans_df[3]), " ")[[1]][2],
                           ".png"),
         occPlot,
         dpi = 600,
         units = "px",
         width = 6000,
         height = 8000)
}

# PRIORITY MAP ---------------------------------------------

### PROCESS SPATIAL DATA

# Process BF cover
coverBF <- project (coverBF[["BF_2015"]],
                    gsub( "units=km", "units=m",
                          st_crs(groupSumsR)$proj4string ))
crs(coverBF) <- crs(groupMeansR)

# Process connectivity
connW <- project (connW[["conn_2015"]],
                  gsub( "units=km", "units=m",
                        st_crs(groupSumsR)$proj4string ))
crs(connW) <- crs(groupMeansR)

# What are priority landscapes for improving connectivity? 
# Based on previous plot from script 13, 
# they are landscapes with approximately <20% cover and
# <75 amps connectivity?
priorityData <- c(coverBF, connW) %>%
  as.data.frame(., xy = TRUE) %>%
  mutate(priority = case_when(BF_2015 < 0.2 &
                                conn_2015 < 75 ~ "TRUE",
                              TRUE ~ "FALSE"))

# Assign data (broadleaf species richness)
occ_df <- names(groupSumsR) %>% 
  grep("Y ", .) %>%
  groupSumsR[[.]] %>%
  sum %>%
  as.data.frame(., xy = TRUE)

### CREATE DATA FRAME
  
# Create bivariate data frame
bivariate_df <- full_join(occ_df,
                          priorityData[,c("x" ,"y","priority")]) %>%
  na.omit

# Change names (no spaces, intuitive)
names(bivariate_df) <- c("x", "y", "Occurrence", "Priority")
  
# Create 3 quantile buckets for occupancy
quantilesOcc <- bivariate_df %>%
  pull(Occurrence) %>%
  quantile(probs = 0:3/3, na.rm = TRUE)

# Cut into groups
bivariate_df <- bivariate_df %>%
  mutate(Occ_quantiles = cut(Occurrence, # based on quantiles
                             breaks = quantilesOcc,
                             include.lowest = TRUE) %>% 
           as.numeric) %>%
  mutate(plotValues = if_else(Priority == FALSE, 0, Occ_quantiles))

### PLOT

# Create richness map
richnessMap <- ggplot(data = occ_df) +
  
  # Add raster data
  geom_tile(aes( x = x, y = y,
                 fill = sum, colour = sum)) +
  scale_fill_distiller(aesthetics = c("fill", "colour"),
                    name = "Summed broadleaf-associated\nspecies occurrence probability",
                    palette = "OrRd",
                    direction = 1) +
  coord_fixed() +
  
  # Add country boundaries
  geom_sf(data = sf::st_as_sf(UK),
          fill = "NA",
          colour = "black",
          inherit.aes = FALSE) +
  geom_sf(data = sf::st_as_sf(Ireland),
          fill = "NA",
          colour = "black",
          inherit.aes = FALSE) +
  geom_sf(data = sf::st_as_sf(IsleOfMan),
          fill = "NA",
          colour = "black",
          inherit.aes = FALSE) +
  
  theme_void()+
  theme(legend.position = c(0.3,0.9),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))

# Create priority map
priorityMap <- ggplot(data = bivariate_df) +
  
  # Add raster data
  geom_tile(aes( x = x, y = y,
    fill = as.factor(plotValues),
    colour = as.factor(plotValues))) +
  coord_fixed() +
  
  # Set colours
  scale_fill_manual(aesthetics = c("fill", "colour"),
                    "",
                    breaks = c(3,2,1,0),
                    values = c("#31a354", "#addd8e", "#ffffcc","black"),
                    labels = c(
                      "Connectivity opportunity with high occurrence -\nprioritise resilience",
                      "Connectivity opportunity -\nbalance colonisation and resilience",
                      "Connectivity opportunity with low occurrence - \nprioritise colonisation",
                      "Outside connectivity opportunity space\n(high existing cover/connectivity)")) +
  
  # Add dotted lines around priority areas
  guides(color = guide_legend(override.aes = list(linetype = c(2, 2, 2, 0),
                                                  colour = "black",
                                                  linewidth = 0.5))) +
  
  # Add country boundaries
  geom_sf(data = sf::st_as_sf(UK),
          fill = "NA",
          colour = "black",
          inherit.aes = FALSE) +
  geom_sf(data = sf::st_as_sf(Ireland),
          fill = "NA",
          colour = "black",
          inherit.aes = FALSE) +
  geom_sf(data = sf::st_as_sf(IsleOfMan),
          fill = "NA",
          colour = "black",
          inherit.aes = FALSE) +
  
  theme_void()+
  theme(legend.position = c(0.3,0.9),
        legend.text = element_text(size = 16))

# Adjust tile plot for multi-pane plot
tilePlot <- tilePlot +
  theme(text = element_text(size = 20))

# Add everything together!
bivarPlot <- cowplot::ggdraw(clip = "on") +
  cowplot::draw_plot(tilePlot, 0, 0.55, 1, 0.45) +
  cowplot::draw_plot(richnessMap, 0, 0, 0.5, 0.55, vjust = 0.05) +
  cowplot::draw_plot(priorityMap, 0.5, 0, 0.5, 0.55, vjust = 0.05) +
  cowplot::draw_label("(a)", 0.015, 0.995, size = 22) +
  cowplot::draw_label("(b)", 0.015, 0.525, size = 22) +
  cowplot::draw_label("(c)", 0.515, 0.525, size = 22) +
  theme(plot.background = element_rect( fill = "white", colour = "white"))
  
# Save
ggsave(filename = paste0("../Writing/Plots/", "PrioritiesPlot.png"),
       bivarPlot,
       dpi = 600,
       units = "px", width = 8000, height = 11000)

