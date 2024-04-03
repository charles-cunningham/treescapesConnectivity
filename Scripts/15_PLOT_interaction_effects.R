# HEADER --------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: Plot meta analysis
#
# Script Description: Plot figures showing results of brms meta-analysis

# LOAD LIBRARIES & INSTALL PACKAGES -----------------

# Change  library to C: (R: doesn't have enough space for packages):
.libPaths("R:/rsrch/cb751/lab/Charles/R/PackageLibrary")

# Load packages
library(brms)
library(terra)
library(tidyverse)
library(tidybayes)
library(wesanderson)

### RELOAD OBJECTS ---------------------------------------------------

# List of brms objects
brmsList <- c("connAll_brms", "connBF_brms", "connCF_brms", "connOpen_brms",
              "coverBF_brms", "coverCF_brms",
              "intBF_brms", "intCF_brms")

# Load brms objects
load(file = "../Data/Species_data/brms_objects.RData")

# Load SDM fixed effect summaries
load(file = "../Data/Species_data/SDM_fixed_effect_summaries.RData")

# CREATE INTERACTION PREDICTION DATAFRAME ------------------------------

# # Load scaling parameters
# load("../Data/Spatial_data/DataForInlabru/scalingParams.RData")
# 
# # Set number of prediction intervals
# nSteps <- 100
# 
# ### Create unscaled prediction grid for broadleaf and coniferous woodland
# 
# # Load cover/connectivity spatRasters
# for (i in list.files("../Data/Spatial_data/DataForInlabru/spatRaster",
#                      pattern =  "\\.tif$")) {
# 
#   assign(gsub(".tif", "", i),
#          rast(paste0("../Data/Spatial_data/DataForInlabru/spatRaster/",
#                      i)))
# }
# 
# # Create dataframes from spatRasters
# coverBF_df <- coverBF %>%
#   as.data.frame(.) %>%
#   gather("Year", "Cover")
# coverCF_df <- coverCF %>%
#   as.data.frame(.) %>%
#   gather("Year", "Cover")
# connW_df <- connW %>%
#   as.data.frame(.) %>%
#   gather("Year", "Connectivity")
# 
# # Calculate max connectivity value
# maxConnectivity <- max(connW_df$Connectivity)
# 
# # Create unscaled data frame of cover and connectivity values 
# # to predict over with 100 prediction steps over entire range
# # (will be finding values within intervals so use nSteps +1)
# BF_pred_df <-  expand.grid("BF_pred" = seq(0,
#                                            1,
#                                            length.out = nSteps + 1),
#                            "conn_pred" = seq(0,
#                                              maxConnectivity ,
#                                              length.out = nSteps + 1))
# CF_pred_df <-  expand.grid("CF_pred" = seq(0,
#                                            1,
#                                            length.out = nSteps + 1),
#                            "conn_pred" = seq(0,
#                                              maxConnectivity ,
#                                              length.out = nSteps + 1))
# 
# ### Remove prediction grid cells with no actual data points, we want predictions to
# ### be possible, i.e. high connectivity with very low cover is impossible
# 
# # Create data frame of observed cover and connectivity data points
# BF_data_pts <- data.frame("Cover" = coverBF_df$Cover,
#                           "Connectivity" = connW_df$Connectivity)
# CF_data_pts <- data.frame("Cover" = coverCF_df$Cover,
#                           "Connectivity" = connW_df$Connectivity)
# 
# # Convert cover data points to prediction interval bin number...
# BFcoverBins <- findInterval(BF_data_pts[, "Cover"],
#                           unique(BF_pred_df[, "BF_pred" ]),
#                           all.inside = TRUE) %>%
#   # ... and use to create associated interval mid-point
#   unique(BF_pred_df[, "BF_pred"] + 1/(nSteps*2))[.] 
# # Convert cover data points to prediction interval bin number...
# CFcoverBins <- findInterval(CF_data_pts[, "Cover"],
#                           unique(CF_pred_df[, "CF_pred" ]),
#                           all.inside = TRUE) %>%
#   # ...and use to create associated interval mid-point
#   unique(CF_pred_df[, "CF_pred"] + 1/(nSteps*2))[.]
# 
# # Convert cover data points to prediction interval bin number...
# connBins <- findInterval(BF_data_pts[,"Connectivity"],
#                           unique(BF_pred_df[, "conn_pred"]),
#                          all.inside = TRUE) %>%
#   # ...and use to create associated interval mid-point
#   unique(BF_pred_df[, "conn_pred"] + maxConnectivity/(nSteps*2))[.] # Can use BF as same number of rows
# 
# # Create a data frame of cross-occurrence of cover and connectivity bins
# BFpredBins <- table("Cover" = BFcoverBins,
#                     "Connectivity" = connBins) %>%
#   as.data.frame %>%
#   mutate_all(function(x) {as.numeric(as.character(x))} )
# CFpredBins <- table("Cover" = CFcoverBins,
#                     "Connectivity" = connBins) %>%
#   as.data.frame %>%
#   mutate_all(function(x) {as.numeric(as.character(x))} )
# 
# # Remove bins with no presences
# BFpredBins <- BFpredBins[BFpredBins$Freq > 0,]
# CFpredBins <- CFpredBins[CFpredBins$Freq > 0,]
# 
# ### Scale covariates
# 
# # Scale the prediction steps for broadleaf and coniferous woodland separately, and connectivity
# BFpredBins$coverScaled <- ( BFpredBins$Cover -
#                                scalingParams[scalingParams$variable == "coverBF", "variableMean"] ) /
#   scalingParams[scalingParams$variable == "coverBF", "variableSD"]
# CFpredBins$coverScaled <- ( CFpredBins$Cover -
#                               scalingParams[scalingParams$variable == "coverCF", "variableMean"] ) /
#   scalingParams[scalingParams$variable == "coverCF", "variableSD"]
# 
# BFpredBins$connScaled <-
#   (BFpredBins$Connectivity - scalingParams[scalingParams$variable == "connW", "variableMean"]) /
#   scalingParams[scalingParams$variable == "connW", "variableSD"]
# CFpredBins$connScaled <-
#   (CFpredBins$Connectivity - scalingParams[scalingParams$variable == "connW", "variableMean"]) /
#   scalingParams[scalingParams$variable == "connW", "variableSD"]
# 
# # Calculate scaled interaction terms for prediction
# BFpredBins$coverConnInt <- BFpredBins$coverScaled * BFpredBins$connScaled
# CFpredBins$coverConnInt <- CFpredBins$coverScaled * CFpredBins$connScaled
# 
# ### Calculate simplified prediction from coverScaled, connectivity, and their interaction
# 
# # Create empty comlumns to be populated
# BFpredBins$mean <- BFpredBins$q0.05 <- BFpredBins$q0.95 <-
#   CFpredBins$mean <- CFpredBins$q0.05 <- CFpredBins$q0.95 <-NA
# 
# # For each prediction combination...
# for (i in 1:NROW(BFpredBins)) {
# 
#   # Calculate interaction effect draws from combining random draw effect posteriors
#   # and scaled predication values
#   predDraws <- BFpredBins$coverScaled[i] * spread_draws(coverBF_brms, b_Intercept)$b_Intercept +
#     BFpredBins$connScaled[i] * spread_draws(connBF_brms, b_Intercept)$b_Intercept +
#     BFpredBins$coverConnInt[i] * spread_draws(intBF_brms, b_Intercept)$b_Intercept
# 
#   # Take mean and quantiles
#   BFpredBins$mean[i] <- mean( predDraws )
#   BFpredBins$q0.05[i] <- quantile( predDraws, 0.05 )
#   BFpredBins$q0.95[i] <- quantile( predDraws, 0.95 )
# }
# for (i in 1:NROW(CFpredBins)) {
# 
#   # Calculate interaction effect draws from combining random draw effect posteriors
#   # and scaled predication values
#   predDraws <- CFpredBins$coverScaled[i] * spread_draws(coverCF_brms, b_Intercept)$b_Intercept +
#     CFpredBins$connScaled[i] * spread_draws(connCF_brms, b_Intercept)$b_Intercept +
#     CFpredBins$coverConnInt[i] * spread_draws(intCF_brms, b_Intercept)$b_Intercept
# 
#   # Take mean and quantiles
#   CFpredBins$mean[i] <- mean( predDraws )
#   CFpredBins$q0.05[i] <- quantile( predDraws, 0.05 )
#   CFpredBins$q0.95[i] <- quantile( predDraws, 0.95 )
# }
# 
# # Save interaction prediction object
# save(list = c("BFpredBins", "CFpredBins"),
#      file = "../Data/Species_data/intPred.RData")

# LINE PLOT -----------------------------------------

# Load interaction prediction object
load(file = "../Data/Species_data/intPred.RData")

# Define the ggplot object for broadleaf and coniferous
for (i in c("BF", "CF")) {
  
  linePlot <- 
    get(paste0(i, "predBins")) %>%
    filter(Cover == 0.105  |
             Cover == 0.305 |
             Cover == 0.505) %>% 
    
    ggplot(aes(x = Connectivity,
               y = 1 - exp(-exp(mean)),
               group = as.factor(Cover),
               colour = as.factor(Cover))) +
    
    # Set plot title and axis labels
    xlab("Connectivity") +
    ylab("Relative Occurrence Probability") +
    
    # Define colour scale and legend label
    scale_colour_manual(aesthetics = c("colour", "fill"),
                        name = "",
                        values = wes_palette("Zissou1")[c(1,3,5)] %>% 
                          rev,
                        labels = c("Low cover (10%)",
                                   "Moderate cover (30%)",
                                   "High cover (50%)")) +
    
    
    # Set geom_line and geom_ribbon
    geom_line(linewidth = 1) +
    geom_ribbon(inherit.aes = TRUE,
                aes(ymin = 1 - exp(-exp(q0.05)),
                    ymax = 1 - exp(-exp(q0.95)),
                    fill = as.factor(Cover)),
                alpha = 0.1,
                colour = NA) + 
    
    
    # Set theme for a minimalistic appearance
    theme_classic(base_size = 18) +
    theme(legend.position = c(0.5, 0.14),
          axis.text.y = element_text(size = 16),
          axis.text.x = element_text(size = 16),
          legend.text = element_text(size = 18),
          legend.background = element_blank(),
          legend.title = element_text(size = 18))
  
  # Save
  ggsave(filename = paste0("../Writing/Plots/",
                           "linePlot", i, ".png"),
         linePlot,
         dpi = 600,
         units = "px",
         width = 5000,
         height = 5000)
}

# TILE PLOTS ------------------------------------------

# CONNECTIVITY-OCCUPANCY LINE PLOTS

# Load interaction prediction object
load(file = "../Data/Species_data/intPred.RData")

# Create vector of cover values from 5% to 95%
coverVals <- seq(0.055, 0.955, 0.1)

# Loop through coverVals
for (i in c("BF", "CF")) {
  for (j in coverVals) {
    
    # Create plot
    linePlot <- get(paste0(i, "predBins")) %>%
      filter(Cover == as.character(j)) %>%
      ggplot(aes(x = Connectivity,
                 y = 1 - exp(-exp(mean)))) +
      geom_ribbon(inherit.aes = TRUE,
                  aes(ymin = 1 - exp(-exp(q0.05)),
                      ymax = 1 - exp(-exp(q0.95))),
                  fill = "grey70") +
      xlab(NULL) +
      ylab(NULL) +
      theme_classic() +
      theme(axis.text = element_blank()) +
      ylim(0,1) +
      xlim(0,265)
    
    # Assign plot out of loop
    assign(paste0(i, "linePlot_", j), linePlot,
           envir = globalenv())
  }
}

# TILE PLOT WITHOUT TRANSPARENCY
# N.B. Following warning okay, connected to empty values removed before:
# Warning messages:
#   1: Raster pixels are placed at uneven horizontal intervals and will be shifted
# ℹ Consider using `geom_tile()` instead. 

# Define the ggplot object for broadleaf and coniferous
for (i in c("BF", "CF")) {
  
  # Define the ggplot object
  tilePlot <-
    # Get last two characters from i list
    substr(i, nchar(i) - 1, nchar(i)) %>%
    # Convert to name of connected object
    paste0(., "predBins") %>%
    # Get object
    get %>%
    
    # Start ggplot object
    ggplot(aes(x = Cover,
               y = Connectivity,
               colour = 1 - exp(-exp(mean)))) +
    
    # Set tile aesthetics
    geom_raster(aes(fill = 1 - exp(-exp(mean))), alpha = 0.8) +
    geom_contour(aes(z = 1 - exp(-exp(mean))),
                 colour = "black",
                 binwidth = 0.02, alpha = 0.3) +
    coord_cartesian(expand = FALSE) +
    
    # Define colour scale and legend label
    scale_fill_gradientn(aesthetics = c("fill", "colour"),
                         name = "Occurrence\nProbability",
                         colours = wes_palette("Zissou1",
                                               100,
                                               type = "continuous") %>%
                           rev,
                         guide = guide_colorbar(barwidth = 10)) +
    
    # Set plot title and axes labels
    ylab("Woodland connectivity (amps)") +
    
    # Set theme for a minimalistic appearance
    theme_classic(base_size = 16) +
    theme(legend.position = "bottom",
          plot.margin = margin(5, 20, 5, 5))
  
  #Add x-axis label
    if (i == "BF") {
    tilePlot <- tilePlot +
      xlab("Proportion broadleaf woodland cover")
  } else if (i == "CF") {
    tilePlot <- tilePlot +
      xlab("Proportion coniferous woodland cover")
  }

  # Add connectivity priority box and annotation for broadleaf
  if (i == "BF") {
    tilePlotPlusAnnotation <- tilePlot +
      # Add main connectivity priority box
      geom_rect(aes(xmin = 0,
                    xmax = 0.2,
                    ymin = 0,
                    ymax = 75),
                colour = "black",
                fill = NA, 
                linetype = 2,
                linewidth = 0.5) +
      # Add connectivity box key
      geom_rect(aes(xmin = 0.76,
                    xmax = 0.80,
                    ymin = 190,
                    ymax = 200),
                colour = "black",
                fill = NA,
                linetype = 2,
                linewidth = 0.5) +
      # Add key text
      annotate("text",
               x = 0.82, y = 195, size = 7,
               label = "Connectivity\nopportunity space",
               hjust = 0) +
      # Add arrow
      annotate("segment", x = 0, y = 0, xend = 0.32, yend = 120,
               arrow = arrow(length = unit(0.05, "npc"))) +
      # Add arrow text
      annotate("text",
               x = 0.26, y = 128, size = 6,
               label = "Increasing occurrence of invertebrates\nassociated with broadleaf woodland",
               angle = 42) +
      # Adjust legend position
      theme(legend.position = c(0.84,0.91)) +
      guides(fill = guide_colourbar(title.position="right",
                                    label.position = "left",
                                    title.hjust = 0))

    # Save object
    saveRDS(tilePlotPlusAnnotation,
            file = paste0("../Data/Species_Data/",
                          "tilePlot", i, ".RDS"))
  }
  
  # Save
  ggsave(filename = paste0("../Writing/Plots/",
                           "tilePlot", i, ".png"),
         tilePlot,
         dpi = 600,
         units = "px",
         width = 6000,
         height = 5000)
}

# TILE PLOT WITH TRANSPARENCY
# N.B. Following warning okay, connected to empty values removed before:
# Warning messages:
#   1: Raster pixels are placed at uneven horizontal intervals and will be shifted
# ℹ Consider using `geom_tile()` instead. 

# Define the ggplot object for broadleaf and coniferous
# (and for supplemental with additional info or not)
for (i in c("BF", "CF")) {
  
  # Define the ggplot object
  tilePlot <-
    # Get last two characters from i list
    substr(i, nchar(i) - 1, nchar(i)) %>%
    # Convert to name of connected object
    paste0(., "predBins") %>%
    # Get object
    get %>%
    
    # Start ggplot object
    ggplot(aes(x = Cover,
               y = Connectivity,
               z = 1 - exp(-exp(mean)))) +
    
    # Set tile aesthetics
    geom_raster(aes(fill = 1 - exp(-exp(mean)),
                  alpha = log(Freq))) +
    geom_contour(colour = "black", binwidth = 0.02, alpha = 0.3) +
    coord_cartesian(expand = FALSE) +
    
    # Define colour scale and legend label
    scale_fill_gradientn(aesthetics = c("fill", "colour"),
                         name = "Occurrence\nProbability",
                         colours = wes_palette("Zissou1",
                                               100,
                                               type = "continuous") %>%
                           rev,
                         guide = guide_colorbar(barwidth = 10)) +
    scale_alpha_continuous(name = " Log frequency\n of landscape") +

    # Set y-axis label
    ylab("Woodland connectivity (amps)") +
    
    # Set theme for a minimalistic appearance
    theme_classic(base_size = 16) +
    theme(legend.position = "bottom",
          plot.margin = margin(5, 20, 5, 5))
  
  #Add x-axis label
  if (i == "BF") {
    tilePlot <- tilePlot +
      xlab("Proportion broadleaf woodland cover")
  } else if (i == "CF") {
    tilePlot <- tilePlot +
      xlab("Proportion coniferous woodland cover")
  }
  
  # Save basic plots
  ggsave(filename = paste0("../Writing/Plots/",
                           "tilePlotAlpha", i, ".png"),
         tilePlot,
         dpi = 600,
         units = "px",
         width = 6000,
         height = 5000)
  
  ### Add small multiples for SI
  tilePlotSI <- tilePlot +
    
    # Raise y axis limit
    ylim(0,300) +
    
    # Add inset line plots
    annotation_custom(ggplotGrob(get(paste0(i, "linePlot_0.055"))),
                      xmin = 0,
                      xmax = 0.1,
                      ymin = 250,
                      ymax = 300) +
    annotation_custom(ggplotGrob(get(paste0(i, "linePlot_0.155"))),
                      xmin = 0.1,
                      xmax = 0.2,
                      ymin = 250,
                      ymax = 300) +
    annotation_custom(ggplotGrob(get(paste0(i, "linePlot_0.255"))),
                      xmin = 0.2,
                      xmax = 0.3,
                      ymin = 250,
                      ymax = 300) +
    annotation_custom(ggplotGrob(get(paste0(i, "linePlot_0.355"))),
                      xmin = 0.3,
                      xmax = 0.4,
                      ymin = 250,
                      ymax = 300) +
    annotation_custom(ggplotGrob(get(paste0(i, "linePlot_0.455"))),
                      xmin = 0.4,
                      xmax = 0.5,
                      ymin = 250,
                      ymax = 300) +
    annotation_custom(ggplotGrob(get(paste0(i, "linePlot_0.555"))),
                      xmin = 0.5,
                      xmax = 0.6,
                      ymin = 250,
                      ymax = 300) +
    annotation_custom(ggplotGrob(get(paste0(i, "linePlot_0.655"))),
                      xmin = 0.6,
                      xmax = 0.7,
                      ymin = 250,
                      ymax = 300) +
    annotation_custom(ggplotGrob(get(paste0(i, "linePlot_0.755"))),
                      xmin = 0.7,
                      xmax = 0.8,
                      ymin = 250,
                      ymax = 300) +
    annotation_custom(ggplotGrob(get(paste0(i, "linePlot_0.855"))),
                      xmin = 0.8,
                      xmax = 0.9,
                      ymin = 250,
                      ymax = 300) +
    annotation_custom( ggplotGrob(get(paste0(i, "linePlot_0.955"))),
                       xmin = 0.9,
                       xmax = 1,
                       ymin = 250,
                       ymax = 300) +
      # Add axes labels
      annotate("text",
               x = 0.5, y = 245, size = 4,
               label = "Woodland connectivity") +
      annotate("text",
               x = -0.03, y = 275, size = 4,
               label = "Occurrence\nProbability",
               angle = 90) +
      # Add white space around data to allong aex text to fit
      coord_cartesian(expand = TRUE)

  # Save complex plots
  ggsave(filename = paste0("../Writing/Plots/",
                           "tilePlotAlpha", i, "_SI.png"),
         tilePlotSI,
         dpi = 600,
         units = "px",
         width = 6000,
         height = 5000)
}
