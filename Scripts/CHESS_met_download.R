# HEADER ---------------------------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: CHESS met download
#
# Script Description: Download files quickly from CEH CHESS met database 
# through webscraping approach for years 1970 to 2017
# (2018 and 2019 have to be downloaded manually separately 
# as not uploaded to the EIDC website yet)
#

### LOAD LIBRARIES & INSTALL PACKAGES ---------------------------------

# Change  library to C: (permission issues on R:)
.libPaths("R:/rsrch/cb751/lab/Charles/R/PackageLibrary")

# Load packages
library(rvest)

### ENTER URLs AND FILE LOCATIONS ------------------------------------

# Login page
CEHloginPage <- "https://catalogue.ceh.ac.uk/sso/login"

# CEH download page
linkCEH <- "https://catalogue.ceh.ac.uk/datastore/eidchub/2ab15bf0-ad08-415c-ba64-831168be7293/"
 
# Specific climate variable to download
climateVar <- "tas" # tas = surface air temperature

# Username and password (ENTER HERE)
username <- "*****"
password <- "*****"

# Set folder in which to save files
destFolder <- "../Data/Spatial_data/CHESS/CHESS_met/Current/chess_tas_1970_2017"

# Create destination folder if it doesn't exist
if(!file.exists(destFolder)) { dir.create(destFolder)}

### OPEN SESSION AND LOGIN  ----------------------------------------

# Open session using login page
pgsession <- session(CEHloginPage)

# Extract login form for session (in this case the submit is the 1st form)
pgform <- html_form(pgsession)[[1]]

# Fill in login form using username/password
filled_form <- set_values(pgform, username = username, password = password)

# Submit form
submit_form(pgsession, filled_form)

### NAVIGATE TO DOWNLOAD PAGE & EXTRACT DOWNLOAD LINKS  ----------

# Create URL containing downloadable files
URL <- paste0(linkCEH,  climateVar, "/")

# Navigate to URL
page <- session_jump_to(pgsession, URL)

# Extract all links
links <- html_attr(html_nodes(page, "a"), "href")

# Subset to only download links containing CHESS data, and then only >= 1970
links <- links[grep("chess-met",links)]
links <- links[grep("daily_196",links, invert = TRUE)]

### DOWNLOAD FILES -----------------------------------------------

# For every link....
for (i in links) {
  
  # Jump to download link
  download <- session_jump_to(page, paste0(URL, i))
  
  # Extract content and save in destination folder
  writeBin(download$response$content, paste0(destFolder, "/", i ))
  
}
