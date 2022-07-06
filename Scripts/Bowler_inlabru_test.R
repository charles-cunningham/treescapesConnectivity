#code partly from PhD student of Colin: 
#https://github.com/PhilipMostert/TZ_inla_spatial_temporal/blob/main/Scripts/TZ_sp_distribution_workflow_inlabru.R

library(INLA)
library(data.table)
library(tidyverse)
library(lubridate)
library(raster)
library(rgdal)
library(rlang)
library(inlabru)
library(BRCmap)
library(sf)
library(sp)
library(ggpolypath)
library(mapview)
library(leafsync)
library(RColorBrewer)

### get species data into visits #######

# sort species data
load("W:/PYWELL_SHARED/Pywell Projects/BRC/Charlie/1.c. New Model Rerun/1. Data/Cleaned Datasets/Dragonflies_170310_Cleaned_Data.rdata")

# choose taxa
myspecies <- "calopteryx_virgo"

# put into visits
taxa_data$visit <- paste(taxa_data$TO_GRIDREF, taxa_data$TO_STARTDATE, sep="_")

visitData <- taxa_data %>%
  add_column(Presence = 1) %>%
  pivot_wider(names_from = CONCEPT,
              values_from = Presence) %>%
  janitor::clean_names(case = "all_caps") %>%
  rename(Presence = toupper(myspecies)) %>%
  dplyr::select(TO_GRIDREF,TO_STARTDATE,
                YEAR, VISIT, Presence)


#change non-detections to zeros
visitData$Presence[is.na(visitData$Presence)] <- 0
table(visitData$Presence)

### add coordinates ########

allGrids <- sort(unique(visitData$TO_GRIDREF))
coords <- OSgrid2GB_EN(gridref = allGrids) %>%
  as_tibble() %>%
  mutate(TO_GRIDREF = allGrids)

#add on to taxa data
visitData <- visitData %>%
  left_join(.,coords,by = "TO_GRIDREF") 

#quick plot
ggplot() +
  geom_point(data = subset(visitData, Presence==0), aes(x=EASTING, y=NORTHING),
             color = "white") +
  geom_point(data = subset(visitData, Presence==1), aes(x=EASTING, y=NORTHING),
             color = "black")

### spatial subset ##############

# get map of wales
Wales <- UK_countries %>%
  st_as_sf() %>%
  filter(COUNTRY=="Wales") %>%
  st_transform(.,crs=27700)
plot(Wales)

# make data spatial
visitData_spatial <- visitData %>%
  st_as_sf(.,coords=c("EASTING","NORTHING"), crs=27700)

# filter non-spatial data
visitData <- visitData %>%
  filter(st_intersects(visitData_spatial, Wales, sparse = FALSE)[,1]) %>%
  filter(YEAR > 1999 & YEAR <2016)

# plot again
ggplot() +
  geom_point(data = subset(visitData, Presence==0), aes(x=EASTING, y=NORTHING),
             color = "white") +
  geom_point(data = subset(visitData, Presence==1), aes(x=EASTING, y=NORTHING),
             color = "black")

#make sp for INLA
visitData_sp <- visitData
coordinates(visitData_sp) <- c("EASTING", "NORTHING")
proj4string(visitData_sp) <- CRS("+init=epsg:27700")
Wales_sp <- as(Wales, 'Spatial')

### inla mesh ##########

#set up mesh
#meshbuilder() 

#rules of thumb
max.edge = diff(range(visitData$NORTHING))/8
bound.outer = diff(range(visitData$NORTHING))/15

mesh <- inla.mesh.2d(boundary = Wales_sp,
                     loc = visitData_sp,
                     max.edge = c(1,4)*max.edge/2,
                     cutoff = max.edge/5,
                     offset = c(max.edge, bound.outer),
                     crs = CRS("+init=epsg:27700"))

ggplot() + 
  gg(mesh)

ggplot() + 
  gg(mesh) + 
  gg(Wales_sp, col = "black") +
  gg(visitData_sp, col = "blue") + 
  coord_fixed() 

### fit spatial model ####

mySpace <- inla.spde2.pcmatern(
  mesh,
  prior.range = c(1, 0.01),   
  prior.sigma = c(4, 0.01))

comp_inlabru <- Presence ~ spatial(main = coordinates, model = mySpace) + Intercept(1)

inlabru1 <- bru(comp_inlabru, data = visitData_sp, 
                family = "binomial", Ntrials = 1)

### predict model ########################

preds = predict(inlabru1, pixels(mesh, mask = Wales_sp), 
                formula = ~ plogis(Intercept + spatial))

ggplot() + 
  gg(preds) +
  gg(Wales_sp) + 
  coord_fixed()+
  scale_fill_viridis_c(option="magma",direction=-1)

# interactive maps
target <- preds
mean.fn <- colorRampPalette(brewer.pal(9, 'YlOrRd'))
map.mean <- mapview(target, zcol = c("mean"), legend = TRUE, na.color=NA, 
                    col.regions = mean.fn(500), map.types = c("OpenStreetMap"))

### fit space-time model ##################

# create indices

#simply year to nearest decade
#floor_decade    = function(value){ return(value - value %% 10) }
#visitData_sp$iYear <- as.factor(floor_decade((visitData_sp$YEAR)))

visitData_sp$iYear <- as.factor(visitData$YEAR)
groups <- visitData_sp$iYear
visitData_sp$ind <- group_indices(visitData_sp@data, iYear)
nyear = length(unique(visitData_sp$iYear))

#define spatial pattern
mySpace <- inla.spde2.pcmatern(
  mesh,
  prior.range = c(1, 0.01),   
  prior.sigma = c(4, 0.01))

#fit = lgcp(cmp, 
#           df, 
#           samplers = sampleMTB_spatial,
#           domain=list(coordinates = mesh, iYear = 1:nyear))


comp_inlabru <- Presence ~ spatial(main = coordinates, 
                                   group = ind,
                                   ngroup = nyear,
                                   model = mySpace,
                                   control.group = list(model = "ar1"))+ 
  Intercept(1)

inlabru1 <- bru(comp_inlabru, data = visitData_sp, 
                family = "binomial", Ntrials = 1)

#predict#
ppxl <- pixels(mesh, mask = Wales_sp)
ppxl_all <- cprod(ppxl, data.frame(iYear = 1:nyear))

preds <- predict(inlabru1, ppxl_all, ~ data.frame(iYear = iYear, 
                                                  lambda = plogis(Intercept + spatial)))

ggplot() + 
  gg(preds, aes(fill=mean)) +
  gg(Wales_sp) + 
  coord_fixed()+
  scale_fill_viridis_c(option="magma",direction=-1) +
  facet_wrap(~iYear)

### include effort covariates ############

### include land cover ####################

#get land cover data

#urban cover
mtbsWGS$UrbanCover <- environData$urban[match(mtbsWGS$Value,environData$Value)]

#make as a spatial pixels data frame
landPixels <- SpatialPixelsDataFrame(points = coordinates(mtbsWGS),
                                     data = data.frame(urbancover = log10(mtbsWGS$UrbanCover+1)),
                                     proj4string = crs(mtbsWGS))
plot(landPixels)#WGS projection

#transform to utm projection
library(raster)
landRaster <- raster(landPixels)
landRaster <- projectRaster(landRaster,crs=(proj4string(mtbs)),method="bilinear")
landPixels <- as(landRaster, "SpatialPixelsDataFrame")
plot(landPixels)

#function to extract covariate info at random points
f.urban <- function(x, y) {
  # turn coordinates into SpatialPoints object:
  # with the appropriate coordinate reference system (CRS)
  spp <- SpatialPoints(data.frame(x = x, y = y), proj4string = fm_sp_get_crs(landPixels))
  proj4string(spp) <- fm_sp_get_crs(landPixels)
  # Extract elevation values at spp coords, from our elev SpatialGridDataFrame
  v <- over(spp, landPixels)
  if (any(is.na(v$urbancover))) {
    v$urbancover <- inlabru:::bru_fill_missing(urbancover, spp, v$urbancover)
  }
  return(v$elevation)
}

#might be needed:
#cmp <- coordinates ~ landPixels(f.urban(x, y), model = "linear") +
#                        mySmooth(coordinates, model = matern) + Intercept(1)

#or simpler option:
cmp <- coordinates ~ urbancover(landPixels, model = "linear") +
  mySmooth(coordinates, model = matern) + Intercept

fit = lgcp(cmp, 
           df, 
           samplers = sampleMTB_spatial,
           domain=list(coordinates = mesh))

fit$summary.fixed

lambda <- predict(fit, pixels(mesh, mask = Germany), ~ exp(mySmooth + urbancover + Intercept))
lambda <- predict(fit, pixels(mesh, mask = Germany), ~ exp(mySmooth + Intercept))
lambda <- predict(fit, pixels(mesh, mask = Germany), ~ exp(urbancover + Intercept))

ggplot() + 
  gg(lambda) + 
  gg(Germany) + 
  coord_fixed()+
  scale_fill_viridis_c(option="magma",direction=-1)

### spatio temporal model #################

### end ####################################