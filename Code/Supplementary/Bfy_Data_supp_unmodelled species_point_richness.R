### This script was annotated with the help of generative AI with internet access ###

# Load necessary libraries
library(maptools)
library(BBmisc)
library(dplyr)
library(terra)

# Define extent of the study area
Xmini <- 69       # Left x-axis limit
Xmaxi <- 161.6    # Right x-axis limit
Ymini <- -10      # Lower y-axis limit
Ymaxi <- 36       # Upper y-axis limit

# Set parameters for output files
resolution <- 10000  # Resolution for raster output
buffer_res <- 30000  # Buffer size for spatial data points
outputname <- paste0("point_richness_unmodelled",buffer_res)
existingoutput_dirname <- "biomod_4.2.4_datapp_2025" # Directory name of existing outputs

# Create a directory for model output files
model.files.dir <- paste0("results_", outputname)
dir.create(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/",model.files.dir)) 

# List modelled species
splist <- list.files(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/clipped_maps_",existingoutput_dirname))

# Set working directory for data processing
wd <- "/lustre1/g/sbs_bonebrake/Eugene/SDMin"
setwd(wd)

# Import our occurrence dataset
data <- read.csv("Occurrence Records of Tropical Asian Butterflies-1970-2024_12Feb2025EJ.csv")

# Format data
colnames(data) <- gsub("decimalLat","decimalLatitude",colnames(data))
colnames(data) <- gsub("decimalLon","decimalLongitude",colnames(data))

data$decimalLatitude <- as.numeric(data$decimalLatitude)
data$decimalLongitude <- as.numeric(data$decimalLongitude)

# Truncate data to study area
data <- subset(data,decimalLatitude>=Ymini)
data <- subset(data,decimalLatitude<=Ymaxi)
data <- subset(data,decimalLongitude>=Xmini)
data <- subset(data,decimalLongitude<=Xmaxi)

# Filter out records with years before 1970
data$year <- as.numeric(data$year)
data <- data[! c(1:nrow(data)) %in% which(data$year<1970),]

# Format species names
data[data$lohman_final_genus_species == "Polygonia c-album","lohman_final_genus_species"] <- "Polygonia c album"
data[data$lohman_final_genus_species == "Polygonia c-aureum","lohman_final_genus_species"] <- "Polygonia c aureum"

# List all species names
Species.list <- sort(unique(data$lohman_final_genus_species))

# Format splist (modelled species)
splist <- gsub("_ensemble_TSSbin.tif", "", splist)
splist <- gsub("\\.", " ", splist)

# Identify unmodelled species, record as Species.list
bool <- Species.list %in% splist
Pending.species.list <- Species.list[!bool]
Species.list <- sort(Pending.species.list)

# Load and prepare the world shapefile for future use
Global1 <- terra::vect("world.shp")
Global <- terra::crop(Global1,ext(Xmini,Xmaxi,Ymini,Ymaxi))
Global <- terra::project(Global,"epsg:6933") 
r.rast <- terra::rast(Global, resolution=c(resolution,resolution))
r <- terra::rasterize(Global, r.rast)

# Set the working directory for saving output files
save_wd <- paste("/lustre1/g/sbs_bonebrake/Eugene/SDMout/",model.files.dir, sep="")
setwd(save_wd) 

# Function to create raster point buffer layer of occurrence points for each unmodelled species
MyBiomodSF <- function(m) {
  
  # Extract data for the current species
  SDM_data <- subset(data,lohman_final_genus_species == Species.list[m]) 
  SDM_data <- SDM_data[!is.na(SDM_data$decimalLatitude),] #Remove negative/NA data
  myRespXY <- SDM_data[,c("decimalLongitude","decimalLatitude")] #Draws out latitude and longitude for each records
  
  # Convert to spatial points with defined CRS
  Projected.co <- terra::vect(myRespXY, geom=c("decimalLongitude","decimalLatitude"), crs="epsg:4326")
  T.Projected.co <- terra::project(Projected.co,"epsg:6933")
  T.Projected.co <- buffer(T.Projected.co, width=buffer_res) #Create a buffer around points
  
  # Rasterize the buffered points onto the base raster
  T.Projected.co.raster <- terra::rasterize(T.Projected.co, r, vals=1)
  T.Projected.co.raster[T.Projected.co.raster>=1] <- 1 #Ensure binary raster map output 
  r[r>0] <- 0 #Reset raster values
  
  # Mask the raster to only include areas with knwon species occurrences (accessible regions)
  map.raster <- mask(r, T.Projected.co.raster, inverse=FALSE, maskvalues=1, updatevalue=1)
  
  # Save the resulting raster as output
  writeRaster(map.raster, paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/",model.files.dir,"/",gsub(" |_",".",Species.list[m]),".tif"), overwrite=TRUE)
  
}

# Use parallel processing to generate point density raster for all unmodelled species
library(parallel)
m <- 1:length(Species.list)
mclapply(m, MyBiomodSF, mc.cores=8) #Assuming 8 cores available on your computer

# List, extract, stack, and plot output (point density rasters)
maps.list.out <- list.files(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/",model.files.dir,"/"))
clippedmaps <- paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/",model.files.dir,"/",maps.list.out)

my_stack <- terra::rast(clippedmaps)  #Create a raster stack from output files
alpha_point <- terra::app(my_stack, fun=sum, cores=8) #Calculate point richness
plot(alpha_point) #Plot the point richness

# Save point richness raster
writeRaster(alpha_point, paste0("/lustre1/g/sbs_bonebrake/Eugene/Others_out/",outputname,"_alpha_point_density.tif"), overwrite=FALSE)

# List, extract, and stack our SDM maps
splist <- list.files(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/clipped_maps_",existingoutput_dirname))
clippedmaps <- paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/clipped_maps_",existingoutput_dirname,"/",splist)
clipped_stack <- terra::rast(clippedmaps)
clipped_stack_mean <- clipped_stack["EMmean"] #Use results of mean ensemble algorithm
SDMoutput <- terra::app(clipped_stack_mean, fun=sum, cores=8)  #Get alpha diversity of modelled species

# Stack, calculate and plot overall diversity from point richness and SDM output 
overall_diversity <- c(SDMoutput, alpha_point) #Combine point richness output and SDM alpha diversity
overall_diversity <- terra::app(overall_diversity, fun=sum, cores=8) #Sum overall diversity
plot(overall_diversity) #Plot overall diversity

# Save overall diversity raster
writeRaster(overall_diversity, paste0("/lustre1/g/sbs_bonebrake/Eugene/Others_out/",outputname,"_overall_diversity.tif"), overwrite=FALSE)
