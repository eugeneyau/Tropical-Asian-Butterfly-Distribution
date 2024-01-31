########## Terra Plot Alpha Diversity ##########
library(terra)
library(parallel)

filename <- "biomod_4.2.4_mask_height_NDVImeansd" #Name of the file storing your SDM outputs
cores <- 8 #Number of cores available on your PC, for speeding up things

Xmini <- 69 #left x-axis, originally 64
Xmaxi <- 161.6 #right x-axis, originally 166
Ymini <- -10 #lower y-axis, originally -13
Ymaxi <- 36 #upper y-axis, originally 37

#Creating new file to store output
clippedrast.dir <- paste0("clipped_maps_", filename) 
dir.create(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/",clippedrast.dir))

#Generating list of species with SDM output
basepaths <- paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filename,"/") #Where I store my SDM outputs, you will have to modify this path
projection <- "/proj_Current/proj_Current_" #Helps R navigate in your file storing SDM results, need to change if your file structure is different
species.list.out0 <- list.files(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filename)) #Extract file names of all files in SDM results folder, should be a list of species names
poss0 <- paste0(basepaths,species.list.out0,projection,species.list.out0,"_ensemble_TSSbin.tif") #Generate file path of SDM-produced prediction (.tif raster file) for all species
tf <- file.exists(poss0) #Check which species has SDM output (.tif file)
species.list.out <- species.list.out0[tf] #Generate list of species with SDM output (.tif file), we only want to work on them
unlist(species.list.out)
specieslistout <- gsub("\\.", " ", species.list.out) #Format list of species names, substitute "." with " ".

#Read csv file which documents presence/absence of each species in each of our 11 zones
boolspzone <- read.csv("/lustre1/g/sbs_bonebrake/Eugene/SDMin/Species_zones_lohman_bool.csv") #Pleasae download the file "Species_zones_lohman_bool.csv"
landmask <- terra::vect("/lustre1/g/sbs_bonebrake/Eugene/SDMin/mask_all_landmass.shp")
landmask <- na.omit(landmask, geom=TRUE)
landmask <- terra::project(landmask,"epsg:6933")

#Generate raster map of world
Global1 <- terra::vect("/lustre1/g/sbs_bonebrake/Eugene/SDMin/world.shp")
terra::crs(Global1)  <- "epsg:4326"
Global <- terra::crop(Global1,ext(Xmini,Xmaxi,Ymini,Ymaxi)) #Crop raster file to define study area // extent(left x-axis, right x-axis, lower y-axis, upper y-axis)
Global <- terra::project(Global,"epsg:6933") 
r.rast <- terra::rast(Global, resolution=c(10000,10000))
r <- terra::rasterize(Global, r.rast)

#Clip output rasters for each species
Clipmapbysp <- function(m) {
  
  #Generate mask of all zones where the species is known to occur
  bool <- subset(boolspzone,Species == specieslistout[m])[,3:13] == TRUE
  maskselect <- c(1,2,3,4,5,6,7,8,9,10,11)[bool]
  spmask <- landmask[maskselect]
  spmask.r <- terra::rasterize(spmask, r.rast)
  
  modprediction <- terra::rast(paste0(basepaths,species.list.out[m],projection,species.list.out[m],"_ensemble_TSSbin.tif"))
  
  #Use the mask to clip SDM output map to preserve only zones where the species is known to occur in
  clippedmap <- terra::mask(modprediction, spmask.r, inverse=FALSE, maskvalues=NA, updatevalue=NA)
  clippedmap[is.na(clippedmap)]<-0
  clippedmap <- terra::mask(clippedmap, r, inverse=TRUE, maskvalues=1, updatevalue=NA)
  
  #Save clipped map as .tif file
  writeRaster(clippedmap, paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/",clippedrast.dir,"/",species.list.out[m],"_ensemble_TSSbin.tif"), overwrite=FALSE)
}

m <- 1:length(specieslistout)
mclapply(m, Clipmapbysp, mc.cores=cores)


### Plot alpha diversity ###

maps.list.out <- list.files(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/",clippedrast.dir,"/"))
clippedmaps <- paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/",clippedrast.dir,"/",maps.list.out)

my_stack <- terra::rast(clippedmaps)
my_stack_ca <- my_stack["EMca"] 
alpha_ca <- terra::app(my_stack_ca, fun=sum, cores=cores)
plot(alpha_ca) #Committee averaging ensemble model results

my_stack <- terra::rast(clippedmaps)
my_stack_mean <- my_stack["EMmean"] 
alpha_mean <- terra::app(my_stack_mean, fun=sum, cores=cores)
plot(alpha_mean) #Mean ensemble model results
