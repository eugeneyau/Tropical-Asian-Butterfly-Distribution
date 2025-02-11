### This script was annotated with the help of generative AI with access to internet ###

# Load necessary libraries
library(maptools)
library(BBmisc)
library(dplyr)
library(terra)

# Set working directory
wd <- "/lustre1/g/sbs_bonebrake/Eugene/SDMin"
setwd(wd)

# Define study area
Xmini <- 69      # Left x-axis limit
Xmaxi <- 161.6   # Right x-axis limit
Ymini <- -10     # Lower y-axis limit
Ymaxi <- 36      # Upper y-axis limit

#Import species data
data <- read.csv("10Feb2025_SEA_Occurrence_dataset_allinclusive_EY.csv")

# Standardize column names for latitude and longitude
colnames(data) <- gsub("decimalLat", "decimalLatitude", colnames(data))
colnames(data) <- gsub("decimalLon", "decimalLongitude", colnames(data))

# Convert latitude and longitude to numeric
data$decimalLatitude <- as.numeric(data$decimalLatitude)
data$decimalLongitude <- as.numeric(data$decimalLongitude)

# Convert year to numeric and filter out records before 1970
data$year <- as.numeric(data$year)
Cleaned.data <- data[!c(1:nrow(data)) %in% which(data$year < 1970),]
Cleaned.data <- Cleaned.data[,c('lohman_final_genus_species',"decimalLongitude","decimalLatitude")]

# Load shapefile of world map
Global <- terra::vect("world.shp")

# Crop the global shapefile to the defined study area
Global <- terra::crop(Global,terra::ext(Xmini,Xmaxi,Ymini,Ymaxi))
Global <- terra::project(Global, "epsg:6933")  # Project to a specific coordinate reference system
r.rast <- terra::rast(Global, resolution = c(10000, 10000))  # Create a raster with 10km x 10km resolution
r <- terra::rasterize(Global, r.rast)  # Rasterize the cropped shapefile

# Load MASS library for statistical functions
library(MASS)

Cleaned.data <- Cleaned.data[!is.na(Cleaned.data$decimalLatitude),] # Remove negative/NA data
XY <- Cleaned.data[,c("decimalLongitude","decimalLatitude")] # Draws out latitude and longitude for each records
alldata <- terra::vect(XY, geom=c("decimalLongitude","decimalLatitude"), crs="epsg:4326") # Converts to point file and give CRS, for it to be identified as GIS stuff
alldata.proj <- terra::project(alldata,"epsg:6933")
alldata.proj <- terra::crop(alldata.proj,ext(r))
alldata.proj.raster <- terra::rasterize(alldata.proj, r, vals=1) # Use base raster to make raster which has data only at occurance point
alldata.proj.raster[alldata.proj.raster>=1] <- 1 # Limit values to 1 in every 10km x 10km grid square 

# Reproject the raster back to original CRS
alldata.rast.reproj <- terra::project(alldata.proj.raster,"epsg:4326")

# Extract coordinates of occurrence points
df <- crds(alldata.rast.reproj, na.rm=TRUE) 

# Perform kernel density estimation on occurrence points
dens <- kde2d(df[,1], df[,2], n = c(nrow(alldata.rast.reproj), ncol(alldata.rast.reproj)))
dens.ras <- raster::raster(dens)  # Convert density estimate to raster
dens.ras <- terra::rast(dens.ras)
dens.ras <- terra::project(dens.ras, "epsg:6933")
dens.ras <- terra::resample(dens.ras, r, method = "average")  # Resample to match the original raster resolution

# Set output directory for results
shortcut_wd <- "/lustre1/g/sbs_bonebrake/Eugene/SDMin/shortcut"#SDM output result file destination
setwd(shortcut_wd)

# Write the density raster to a file
terra::writeRaster(dens.ras, "dens_ras.tif", overwrite=FALSE)
