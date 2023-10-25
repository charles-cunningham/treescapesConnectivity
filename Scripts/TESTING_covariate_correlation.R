# HEADER --------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: Check covariate correlation
#
# Script Description: 

# LOAD LIBRARIES & INSTALL PACKAGES -----------------

# Change  library to C: (R: doesn't have enough space for packages):
.libPaths("R:/rsrch/cb751/lab/Charles/R/PackageLibrary")

# Load packages
library(terra)
library(tidyverse)
library(corrplot)

# DATA FILES ------------------------------------------

load(file = "../Data/Spatial_data/DataForInlabru.RData")

# SET PARAMETERS ------------------------------------

# List all covariate spatial grid data frames
cov_SGDF <- c(coverBF_SGDF, coverCF_SGDF, connW_SGDF,
              coverBF_connW_SGDF, coverCF_connW_SGDF,
              GDD5_SGDF_grp, WMIN_SGDF_grp, tasCV_SGDF_grp, soilM_SGDF_grp, RAIN_SGDF_grp)

covNames <- c("Broadleaf Cover", "Coniferous Cover", "Connectivity",
               "Broadleaf-Connectivity interaction",
               "Coniferous-Connectivity interaction",
               "GDD5", "WMIN", "tasCV", "soilM", "RAIN")

# PROCESS DATA --------------------------------------

cov_R_1990 <- lapply(cov_SGDF, function(x) { rast(x[1]) }) %>% rast(.)
cov_R_2015 <- lapply(cov_SGDF, function(x) { rast(x[2]) }) %>% rast(.)

names(cov_R_1990) <- names(cov_R_2015) <- covNames

# CHECK CORRELATION ---------------------------------

full_df <- as.data.frame(cov_R_2015) %>% drop_na(.)
res <- cor(full_df)
#png(filename = 'figures/corr_plot.png', width = 45, height = 30, units = "cm", res = 300)
corrplot(res, type = "upper", order = "hclust",
         method = "color", addCoef.col="black", number.cex=0.75,
         tl.col = "black", tl.srt = 45, tl.cex = 0.8)
dev.off()
