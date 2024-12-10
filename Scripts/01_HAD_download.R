# HEADER ---------------------------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: HAD data download (current)
#
# Script Description: Download files quickly from HAD database 
# through webscraping approach for years 1979 to present
#
# N.B. Original data downloaded for modelling was HaDUK-Grid v1.1 
# (https://data.ceda.ac.uk/badc/ukmo-hadobs/data/insitu/MOHC/HadOBS/HadUK-Grid/v1.1.0.0)
# but have updated link to most recent version for script
# Also, may timeout so have to run several times

### LOAD LIBRARIES & INSTALL PACKAGES ---------------------------------

# Change  library to R: (not enough space on C:)
.libPaths("R:/rsrch/cb751/lab/Charles/R/PackageLibrary")

# Load packages
library(tidyverse)
library(rvest)

### SET URLs, FILE LOCATIONS, PARAMETERS -----------------------------

# Login page
loginPageCEDA <- "https://auth.ceda.ac.uk/account/signin"

# Dataset link 
dataCEDA <- "https://data.ceda.ac.uk/badc/ukmo-hadobs/data/insitu/MOHC/HadOBS/HadUK-Grid/v1.2.0.ceda"

# Username and password (ENTER HERE)
username <- "*****"
password <- "*****"
  
# Set years to ignore
yearCutoff <- 1978 # >1980 but need 1979 for winter cold [MTCO] calculation

### OPEN SESSION AND LOGIN  ----------------------------------------

# Open session using login page
pgsession <- session(loginPageCEDA)

# Extract login form for session (in this case the submit is the 1st form)
pgform <- html_form(pgsession)[[1]]

# Fill in login form using username/password
filled_form <- html_form_set(pgform, username = username, password = password)

# Submit form
session_submit(pgsession, filled_form)

### TOTAL PRECIPITATION AMOUNT --------------------------
# N.B. HadUK data in mm, but CHESS-SCAPE data in kg/m^2/s,
# to convert multiply kg/m^2/s by 86400

# NAVIGATE TO DOWNLOAD PAGE & EXTRACT DOWNLOAD LINKS

# CEH download page
linkCEDA <- paste0(dataCEDA, "/1km/rainfall/mon/latest")

# Navigate to URL
page <- session_jump_to(pgsession, linkCEDA)

# Extract all links
links <- html_attr(html_nodes(page, "a"), "href")

# Subset to only download links 
links <- links[grep("download",links)]

# Extract year of data
linkYear <- str_extract(links, "(?<=rainfall_hadukgrid_uk_1km_mon_).{4}") %>% 
  as.integer()

# Subset to years above cutoff year 
links <- links[linkYear > yearCutoff]

# DOWNLOAD FILES

# Set folder in which to save files
destFolder <- "../Data/Spatial_data/HAD/rainfall"

# Create destination folder if it doesn't exist
if(!file.exists(destFolder)) { dir.create(destFolder)}

# For every link....
for (i in links) {
  
  # Create destination file name
  fileName <- paste0(destFolder,
                     "/",
                     sub("\\?download=1", "", basename(i)))
  
  # If file doesn't exist, then download
  if (!file.exists(fileName)) {
  
    # Jump to download link
    download <- session_jump_to(page,  i)
    
    # Extract content and save in destination folder
    writeBin(download$response$content, fileName )
  }
}

### CURRENT DAILY MINIMUM TEMPERATURE --------------------------

# NAVIGATE TO DOWNLOAD PAGE & EXTRACT DOWNLOAD LINKS

# CEH download page
linkCEDA <- paste0(dataCEDA, "/1km/tasmin/day/latest")

# Navigate to URL
page <- session_jump_to(pgsession, linkCEDA)

# Extract all links
links <- html_attr(html_nodes(page, "a"), "href")

# Subset to only download links 
links <- links[grep("download",links)]

# Extract year of data
linkYear <- str_extract(links, "(?<=tasmin_hadukgrid_uk_1km_day_).{4}") %>% 
  as.integer()

# Subset to years above cutoff year 
links <- links[linkYear > yearCutoff]

# DOWNLOAD FILES

# Set folder in which to save files
destFolder <- "../Data/Spatial_data/HAD/tasmin"

# Create destination folder if it doesn't exist
if(!file.exists(destFolder)) { dir.create(destFolder)}

# For every link....
for (i in links) {
  
  # Create destination file name
  fileName <- paste0(destFolder,
                     "/",
                     sub("\\?download=1", "", basename(i)))
  
  # If file doesn't exist, then download
  if (!file.exists(fileName)) {
    
    # Jump to download link
    download <- session_jump_to(page,  i)
    
    # Extract content and save in destination folder
    writeBin(download$response$content, fileName )
  }
}

### CURRENT DAILY MAXIMUM TEMPERATURE --------------------------

# NAVIGATE TO DOWNLOAD PAGE & EXTRACT DOWNLOAD LINKS

# CEH download page
linkCEDA <- paste0(dataCEDA, "/1km/tasmax/day/latest")

# Navigate to URL
page <- session_jump_to(pgsession, linkCEDA)

# Extract all links
links <- html_attr(html_nodes(page, "a"), "href")

# Subset to only download links 
links <- links[grep("download",links)]

# Extract year of data
linkYear <- str_extract(links, "(?<=tasmax_hadukgrid_uk_1km_day_).{4}") %>% 
  as.integer()

# Subset to years above cutoff year 
links <- links[linkYear > yearCutoff]

# DOWNLOAD FILES

# Set folder in which to save files
destFolder <- "../Data/Spatial_data/HAD/tasmax"

# Create destination folder if it doesn't exist
if(!file.exists(destFolder)) { dir.create(destFolder)}

# For every link....
for (i in links) {
  
  # Create destination file name
  fileName <- paste0(destFolder,
                     "/",
                     sub("\\?download=1", "", basename(i)))
  
  # If file doesn't exist, then download
  if (!file.exists(fileName)) {
    
    # Jump to download link
    download <- session_jump_to(page,  i)
    
    # Extract content and save in destination folder
    writeBin(download$response$content, fileName )
  }
}
