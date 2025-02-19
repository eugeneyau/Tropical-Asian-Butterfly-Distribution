### This script was annotated with the help of generative AI with internet access ###

########## Terra Plot Alpha Diversity ##########
# Load necessary libraries
library(terra)
library(parallel)

filename <- "biomod_4.2.4_datapp_2025" # Name of the file storing your SDM outputs
cores <- 8 # Number of cores available on your PC, for speeding up computations

# Define the extent of the study area
Xmini <- 69 #left x-axis
Xmaxi <- 161.6 #right x-axis
Ymini <- -10 #lower y-axis
Ymaxi <- 36 #upper y-axis

# Create a new directory to store output files
clippedrast.dir <- paste0("clipped_maps_", filename) 
dir.create(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/",clippedrast.dir))

# Generating list of species with SDM output
basepaths <- paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filename,"/") # Path to SDM outputs
projection <- "/proj_Current/proj_Current_" # Path for navigating SDM results
species.list.out0 <- list.files(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filename)) # List of species names
poss0 <- paste0(basepaths,species.list.out0,projection,species.list.out0,"_ensemble_TSSbin.tif") # Generate file paths for SDM predictions
tf <- file.exists(poss0) # Check which species have SDM output 
species.list.out <- species.list.out0[tf] # Filter species with valid SDM outputs
unlist(species.list.out) # Flatten the list of species
specieslistout <- gsub("\\.", " ", species.list.out) # Format species names

# Read csv file which documents presence/absence of each species in each of our zones
boolspzone <- read.csv("/lustre1/g/sbs_bonebrake/Eugene/SDMin/Species_zones_lohman_bool.csv") # Pleasae download the file
landmask <- terra::vect("/lustre1/g/sbs_bonebrake/Eugene/SDMin/mask_all_landmass.shp") # Pleasae download the file
landmask <- na.omit(landmask, geom=TRUE) # Remove NA geometries
landmask <- terra::project(landmask,"epsg:6933")
landmask <- terra::sort(landmask, "layer")

# Generate raster map of world
Global1 <- terra::vect("/lustre1/g/sbs_bonebrake/Eugene/SDMin/world.shp")
terra::crs(Global1)  <- "epsg:4326"
Global <- terra::crop(Global1,terra::ext(Xmini,Xmaxi,Ymini,Ymaxi)) # Crop raster file to define study area
Global <- terra::project(Global,"epsg:6933") 
r.rast <- terra::rast(Global, resolution=c(10000,10000))
r <- terra::rasterize(Global, r.rast)

# Function to clip output rasters for each species
Clipmapbysp <- function(m) {
  
  #Generate mask of all zones where the species is known to occur
  bool <- subset(boolspzone,Species == specieslistout[m])[,3:13] == TRUE # Get the zones where the species is known to occur in
  maskselect <- c(1,2,3,4,5,6,7,8,9,10,11)[bool]
  spmask <- landmask[maskselect] # Create custom mask based on selected zones
  spmask.r <- terra::rasterize(spmask, r.rast)
  
  # Load the model prediction for the species
  modprediction <- terra::rast(paste0(basepaths,species.list.out[m],projection,species.list.out[m],"_ensemble_TSSbin.tif"))
  
  # Clip the SDM output map to preserve only zones where the species is known to occur (therefore accessible to them)
  clippedmap <- terra::mask(modprediction, spmask.r, inverse=FALSE, maskvalues=NA, updatevalue=NA)
  clippedmap[is.na(clippedmap)]<-0 # Replace NA values with 0
  clippedmap <- terra::mask(clippedmap, r, inverse=TRUE, maskvalues=1, updatevalue=NA) # Mask out oceans and seas
  
  #Save clipped map as .tif file
  writeRaster(clippedmap, paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/",clippedrast.dir,"/",species.list.out[m],"_ensemble_TSSbin.tif"), overwrite=FALSE)
}

# Apply the clipping function to all species in parallel
m <- 1:length(specieslistout)
mclapply(m, Clipmapbysp, mc.cores=cores) # Use parallel processing to clip maps


### Plot alpha diversity ###

# List the clipped maps and create a raster stack
maps.list.out <- list.files(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/",clippedrast.dir,"/"))
clippedmaps <- paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/",clippedrast.dir,"/",maps.list.out)

my_stack <- terra::rast(clippedmaps) # Create a raster stack from clipped maps
my_stack_mean <- my_stack["EMmean"] # Extract mean ensemble model results
alpha_mean <- terra::app(my_stack_mean, fun=sum, cores=cores) # Calculate alpha diversity
plot(alpha_mean) # Plot ensemble model alpha diversity results


### Save as raster ###

wd <- "/lustre1/g/sbs_bonebrake/Eugene"
setwd(wd)

writeRaster(alpha_ca, paste0(filename,"_alpha_ca.tif"), overwrite=FALSE)
writeRaster(alpha_mean, paste0(filename,"_alpha_mean.tif"), overwrite=FALSE)


