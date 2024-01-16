# HEADER --------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: SDM spatial analysis
#
# Script Description: Spatial analysis of SDM inlabru model outputs

# LOAD LIBRARIES & INSTALL PACKAGES -----------------

# Change  library to C: (R: doesn't have enough space for packages):
.libPaths("R:/rsrch/cb751/lab/Charles/R/PackageLibrary")

# Load packages
library(terra)
library(tidyverse)
library(tidybayes)
library(ggridges)

# SET PARAMETERS ------------------------------------

# DATA FILES ------------------------------------------

# Load SDM fixed effect summaries
load(file = "../Data/Species_data/SDM_fixed_effect_summaries.RData")

# DESCRIPTIVE POOLED STATS --------------------------------

# Group by taxa group and effect category
group_df <- meta_df %>%
  group_by( woodSp, connectivity_sig ) 

# Get basic numbers
summarise(group_df, length(species))

# GROUPED COVER AND CONNECTIVITY ANALYSIS -----------------------------------

# Create raster for sum of every pooled group

# Loop through every group in group_df
groupSumsR <- lapply(1:NROW(group_keys(group_df)) , function(i) {
  
  # Split group into separate df, get df i...
  iGroup <- group_split(group_df)[[i]] %>%
    dplyr::select(species, taxa) #... and select taxa and species
  
  # Create list of file names of occupancy plots (Script 11) 
  # for every row in iGroup
  iGroupFiles <- apply(iGroup, 1, function(j) {
    paste0("../Data/Species_data/SDMs/",
           j["taxa"],
           "/",
           j["species"],
           "/occPlot.tif")
  })
  
  # Create spatRast
  iGroupRast <- rast(iGroupFiles)
  
  # Sum
  iGroupSum <- sum(iGroupRast)
  
  # Change spatRast layer name to group key
  names(iGroupSum) <- group_keys(group_df)[i, ] %>%
    as.character  %>%
    paste(., collapse = " ")
  
  return(iGroupSum)
  
  # Join all group sum spatRasts into single spatRast  
}) %>% rast(.) 

# Calculate mean
groupMeansR <- groupSumsR /
  sapply(group_split(group_df), NROW)

# Save
writeRaster(groupSumsR,
            paste0("../Data/Species_data/Summed_occupancy_grouped_by_effects.tif"),
            overwrite = TRUE)
writeRaster(groupMeansR,
            paste0("../Data/Species_data/Mean_occupancy_grouped_by_effects.tif"),
            overwrite = TRUE)

### Plot

# Create plot directory
dir.create("../Writing/Plots/GroupedOccPlots")

# Loop through each cover/connectivity group
for(i in 1:nlyr(groupMeansR)) {
  
  # Create data frame from grouped mean spatRast i
  iGroupMeans_df <- groupMeansR[[i]] %>%
    as.data.frame(., xy = TRUE)
  
  # GGplot
  occPlot <- ggplot(data = iGroupMeans_df,
                    aes(x = x, y = y,
                        colour = iGroupMeans_df[, 3],
                        fill = iGroupMeans_df[, 3])) +
    
    # Set equal coordinates
    coord_equal() +
    
    # Country boundary polygons
    # geom_sf(data = st_as_sf(UK), fill = "grey90", colour = NA, inherit.aes = FALSE) +
    # geom_sf(data = st_as_sf(Ireland), fill = "grey90", colour = NA, inherit.aes = FALSE) +
    # geom_sf(data = st_as_sf(IsleOfMan), fill = "grey90", colour = NA, inherit.aes = FALSE) +
    
    # Rraster
    geom_tile() +
    scale_colour_gradient(low = "#e5f5e0",
                          high = "#31a354",
                          # limits = c(0,0.5),
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
        label.theme = element_text(size = 12),
        title.theme = element_text(size = 12),
        barwidth = unit(2, "lines"),
        barheight = unit(8, "lines") )) +
    
    # Add title
    annotate(
      geom = "text",
      x = 490000,
      y = 900000,
      label = names(iGroupMeans_df[3]),
      size = 6) +
    
    # Add title
    annotate(
      geom = "text",
      x = 90000,
      y = 300000,
      label = paste(NROW(group_split(group_df)[[i]]),
                    "species"),
      size = 6) +
    
    # Set theme parameters
    theme_void() +
    theme(plot.background = element_rect(fill = "white",
                                         colour = "white"),
          legend.position = c(0.8, 0.7),
          plot.margin = margin(-3, 0, -3, 0, "lines"))
  
  ggsave(filename = paste0("../Writing/Plots/GroupedOccPlots/",
                           "Occupancy_",
                           names(iGroupMeans_df[3]), ".png"),
         occPlot,
         dpi = 600,
         units = "px",
         width = 6000,
         height = 8000)
  
}

# CREATION PRIORITY MAPS ----------------------------- IN PROGRESS

groupSumsR <- rast("../Data/Species_data/Summed_occupancy_grouped_by_effects.tif")
groupMeansR <- rast("../Data/Species_data/Mean_occupancy_grouped_by_effects.tif")

test <- c(groupMeansR[[7]], groupMeansR[[8]], groupMeansR[[9]]) %>%
  mean %>%
  as.data.frame(., xy = TRUE)

# Occupancy only
ggplot(data = test) +
  geom_tile(aes(x=x,y=y, fill = mean, colour = mean)) +
  scale_colour_gradient(low = "#e5f5e0",
                        high = "#31a354",
                        # limits = c(0,0.5),
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
      label.theme = element_text(size = 12),
      title.theme = element_text(size = 12),
      barwidth = unit(2, "lines"),
      barheight = unit(8, "lines") )) +
  theme_void() + 
  coord_fixed() 


# Process BF cover
woodCover <- rast("../Data/Spatial_data/DataForInlabru/spatRaster/coverBF_scaled.tif")[["BF_2015"]]
woodCover <- project (woodCover, gsub( "units=km", "units=m",
                                       sf::st_crs(groupSumsR)$proj4string ))
crs(woodCover) <- crs(groupMeansR)

# Process connectivity
woodConn <- rast("../Data/Spatial_data/DataForInlabru/spatRaster/connW.tif")[["conn_2015"]]
woodConn <- project (woodConn, gsub( "units=km", "units=m",
                                     sf::st_crs(groupSumsR)$proj4string ))
crs(woodConn) <- crs(groupMeansR)

# Create bivariate data frame
bivariate_df <- c(woodConn, groupMeansR[["Prefer woodland Pos"]]) %>%
  as.data.frame(., xy = TRUE)

# Change names (no spaces, intuitive)
names(bivariate_df) <- c("x", "y", "Conn", "MeanOcc")

# Create 5 quantile buckets for occupancy
quantilesOcc <- bivariate_df %>%
  pull(MeanOcc) %>%
  quantile(probs = seq(0, 1, length.out = 5), na.rm = TRUE)
quantilesConn <- bivariate_df %>%
  pull(Conn) %>%
  quantile(probs = seq(0, 1, length.out = 5), na.rm = TRUE)

# Create 5x5 colour scale
bivariate_color_scale <- tibble(
  "1-1" = "#d3d3d3", # low x, low y
  "2-1" = "#b6cdcd",
  "3-1" = "#97c5c5",
  "4-1" = "#75bebe",
  "5-1" = "#52b6b6", # high x, low y
  "1-2" = "#cab6c5",
  "2-2" = "#aeb0bf",
  "3-2" = "#91aab9",
  "4-2" = "#70a4b2",
  "5-2" = "#4e9daa",
  "1-3" = "#c098b9",
  "2-3" = "#a593b3",
  "3-3" = "#898ead",
  "4-3" = "#6b89a6",
  "5-3" = "#4a839f",
  "1-4" = "#b77aab",
  "2-4" = "#9e76a6",
  "3-4" = "#8372a0",
  "4-4" = "#666e9a",
  "5-4" = "#476993",
  "1-5" = "#ad5b9c", # low x, high y
  "2-5" = "#955898",
  "3-5" = "#7c5592",
  "4-5" = "#60528d",
  "5-5" = "#434e87" # high x, high y
)%>%
  gather("group", "fill")

# Separate the groups
bivariate_color_scale <- bivariate_color_scale %>%
  separate(group, into = c("MeanOcc", "Conn"),
           sep = "-",
           remove = FALSE) %>%
  mutate(MeanOcc = as.integer(MeanOcc),
         Cover = as.integer(Conn))

# Cut into groups defined above and join fill
bivariate_df <- mutate(bivariate_df,
                       Occ_quantiles = cut(MeanOcc,
                                           breaks = quantilesOcc, # based on quantiles
                                           include.lowest = TRUE),
                       Conn_quantiles = cut(Conn,
                                            breaks = 5, # equal interval cuts 
                                            include.lowest = TRUE,
                                            na.rm = TRUE),
                       # by pasting the factors together as numbers we match the groups defined
                       # in the tibble bivariate_color_scale
                       group = paste0(as.numeric(Conn_quantiles), "-",
                                      as.numeric(Occ_quantiles))
) %>%
  # we now join the actual hex values per "group"
  # so each municipality knows its hex value based on the his gini and avg
  # income value
  left_join(bivariate_color_scale, by = "group")

### Create plot

# Create map
map <- ggplot(data = bivariate_df) +
  ggtitle(paste0("Modelled relative occurence proability of woodland species\n",
                 "with positive association with connectivity against woodland connectivity")) +
  geom_tile(aes(x=x,y=y, fill = fill, colour = fill),
            show.legend = FALSE) +
  scale_fill_identity() +
  scale_colour_identity() +
  theme_void() + 
  coord_fixed() 
geom_sf(data = st_as_sf(smoothUK), fill = "NA", colour = "black", inherit.aes = FALSE)

# Create legend
legend <- ggplot() +
  geom_tile(
    data = bivariate_color_scale,
    mapping = aes(
      x = MeanOcc,
      y = Conn,
      fill = fill)
  ) +
  scale_fill_identity() +
  labs(x = "Woodland connectivity ⟶️",
       y = "Mean Occurence ⟶️") +
  theme_void() +
  # make font small enough
  theme(
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10, angle = 90) ) +
  # quadratic tiles
  coord_fixed()

# Plot map and legend
bivarPlot <- cowplot::ggdraw() +
  cowplot::draw_plot(map, 0, 0, 1, 1) +
  cowplot::draw_plot(legend, 0.65, 0.55, 0.25, 0.25)

ggsave(filename = paste0("../Writing/Plots/", "Quick temp plots/", "conn_bivar.png"),
       bivarPlot,
       dpi = 600,
       units = "px", width = 6000, height = 5000)