
library(terra)

FragStats <- read.csv(paste0("G:/Shared drives/NERC Connected Treescapes project drive/",
                             "WP4 biodiversity, ecosystem function and nature recovery/Data/Land_cover/",
                             "Fragstats/LCM2020/Landscapes_fragstats.csv"))

FragStats <- FragStats %>%
  st_as_sf(., coords = c("X","Y"),
           na.fail = FALSE, # NOT SURE WHY THIS MAKES IT WORK, NEED TO CHECK
           crs = bng)

BLcover <- rast(paste0("G:/Shared drives/NERC Connected Treescapes project drive/",
                       "WP4 biodiversity, ecosystem function and nature recovery/Data/Land_cover/",
                       "Fragstats/LCM2020/Raster/Broadleaved_percentCover.tif"),
                crs = bng)

plot(FragStats)
