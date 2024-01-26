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

#Read csv file which documents presence/absence of each species in each of our 10 zones
boolspzone <- read.csv("/lustre1/g/sbs_bonebrake/Eugene/SDMin/Species_zones_lohman_bool.csv") #Pleasae download the file "Species_zones_lohman_bool.csv"
names(boolspzone) <- c("ID","Species","borneo","easia","india","indonesia","japan","malaysia","newguinea","philippines","seasia","taiwan") 

zones0 <- list("borneo","eastasia","india","indonesia","japan","malaysia","newguinea","philippines","seasia","taiwan") #List of all zones

#Generate raster map of world
Global1 <- terra::vect("/lustre1/g/sbs_bonebrake/Eugene/SDMin/world.shp")
terra::crs(Global1)  <- "epsg:4326"
Global <- terra::crop(Global1,ext(Xmini,Xmaxi,Ymini,Ymaxi)) #Crop raster file to define study area // extent(left x-axis, right x-axis, lower y-axis, upper y-axis)
Global <- terra::project(Global,"epsg:6933") 
r.rast <- terra::rast(Global, resolution=c(10000,10000))

#Clip output rasters for each species
Clipmapbysp <- function(m) {
  
  #Get all zones where the species is known to occur and merge them as one raster mask
  bool <- subset(boolspzone,Species == specieslistout[m])[,3:12] == TRUE
  zones <- zones0[bool]
  as.character(zones)
  
  try({Mask1 <- terra::vect(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMin/landmasknew/mask_",zones[1],".shp"))
  Mask1 <- na.omit(Mask1, geom=TRUE)
  Mask1 <- terra::project(Mask1,"epsg:6933")
  Mask1.r <- terra::rasterize(Mask1, r.rast)
  }, silent = TRUE)
  
  try({Mask2 <- terra::vect(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMin/landmasknew/mask_",zones[2],".shp"))
  Mask2 <- na.omit(Mask2, geom=TRUE)
  Mask2 <- terra::project(Mask2,"epsg:6933")
  Mask2.r <- terra::rasterize(Mask2, r.rast)
  }, silent = TRUE)
  
  try({Mask3 <- terra::vect(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMin/landmasknew/mask_",zones[3],".shp"))
  Mask3 <- na.omit(Mask3, geom=TRUE)
  Mask3 <- terra::project(Mask3,"epsg:6933")
  Mask3.r <- terra::rasterize(Mask3, r.rast)
  }, silent = TRUE)
  
  try({Mask4 <- terra::vect(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMin/landmasknew/mask_",zones[4],".shp"))
  Mask4 <- na.omit(Mask4, geom=TRUE)
  Mask4 <- terra::project(Mask4,"epsg:6933")
  Mask4.r <- terra::rasterize(Mask4, r.rast)
  }, silent = TRUE)
  
  try({Mask5 <- terra::vect(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMin/landmasknew/mask_",zones[5],".shp"))
  Mask5 <- na.omit(Mask5, geom=TRUE)
  Mask5 <- terra::project(Mask5,"epsg:6933")
  Mask5.r <- terra::rasterize(Mask5, r.rast)
  }, silent = TRUE)
  
  try({Mask6 <- terra::vect(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMin/landmasknew/mask_",zones[6],".shp"))
  Mask6 <- na.omit(Mask6, geom=TRUE)
  Mask6 <- terra::project(Mask6,"epsg:6933")
  Mask6.r <- terra::rasterize(Mask6, r.rast)
  }, silent = TRUE)
  
  try({Mask7 <- terra::vect(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMin/landmasknew/mask_",zones[7],".shp"))
  Mask7 <- na.omit(Mask7, geom=TRUE)
  Mask7 <- terra::project(Mask7,"epsg:6933")
  Mask7.r <- terra::rasterize(Mask7, r.rast)
  }, silent = TRUE)
  
  try({Mask8 <- terra::vect(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMin/landmasknew/mask_",zones[8],".shp"))
  Mask8 <- na.omit(Mask8, geom=TRUE)
  Mask8 <- terra::project(Mask8,"epsg:6933")
  Mask8.r <- terra::rasterize(Mask8, r.rast)
  }, silent = TRUE)
  
  try({Mask9 <- terra::vect(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMin/landmasknew/mask_",zones[9],".shp"))
  Mask9 <- na.omit(Mask9, geom=TRUE)
  Mask9 <- terra::project(Mask9,"epsg:6933")
  Mask9.r <- terra::rasterize(Mask9, r.rast)
  }, silent = TRUE)
  
  try({Mask10 <- terra::vect(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMin/landmasknew/mask_",zones[10],".shp"))
  Mask10 <- na.omit(Mask10, geom=TRUE)
  Mask10 <- terra::project(Mask10,"epsg:6933")
  Mask10.r <- terra::rasterize(Mask10, r.rast)
  }, silent = TRUE)
  
  try(Mask1.r <- c(Mask1.r,Mask2.r), silent = TRUE)
  try(Mask1.r <- c(Mask1.r,Mask3.r), silent = TRUE)
  try(Mask1.r <- c(Mask1.r,Mask4.r), silent = TRUE)
  try(Mask1.r <- c(Mask1.r,Mask5.r), silent = TRUE)
  try(Mask1.r <- c(Mask1.r,Mask6.r), silent = TRUE)
  try(Mask1.r <- c(Mask1.r,Mask7.r), silent = TRUE)
  try(Mask1.r <- c(Mask1.r,Mask8.r), silent = TRUE)
  try(Mask1.r <- c(Mask1.r,Mask9.r), silent = TRUE)
  try(Mask1.r <- c(Mask1.r,Mask10.r), silent = TRUE)
  
  Mask1.r[is.na(Mask1.r)]<-0
  
  try(Mask1.r <- terra::app(Mask1.r, fun=sum, cores=8), silent = TRUE)
  
  modprediction <- terra::rast(paste0(basepaths,species.list.out[m],projection,species.list.out[m],"_ensemble_TSSbin.tif"))
  
  #Use the mask to clip SDM output map to preserve only zones where the species is known to occur in
  clippedmap <- mask(modprediction, Mask1.r, inverse=TRUE, maskvalues=1:10, updatevalue=NA)
  clippedmap[is.na(clippedmap)]<-0
  
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
plot(alpha_ca)

my_stack <- terra::rast(clippedmaps)
my_stack_mean <- my_stack["EMmean"] 
alpha_mean <- terra::app(my_stack_mean, fun=sum, cores=cores)
plot(alpha_mean)


### Publish Graph ###


wd <- "/lustre1/g/sbs_bonebrake/Eugene/SDMin"
setwd(wd)

#select truncation
Xmini <- 69 #left x-axis, originally 64
Xmaxi <- 161.6 #right x-axis, originally 166
Ymini <- -10 #lower y-axis, originally -13
Ymaxi <- 36 #upper y-axis, originally 37

#Import species data
Global1 <- terra::vect("world_full.shp")
Global <- terra::crop(Global1,ext(Xmini,Xmaxi,Ymini,Ymaxi)) #Crop raster file to define study area // extent(left x-axis, right x-axis, lower y-axis, upper y-axis)
Global <- terra::project(Global,"epsg:6933") 

clipped_alpha_ca <- mask(alpha_ca, Global, inverse=FALSE, updatevalue=NA)
plot(clipped_alpha_ca)

clipped_alpha_mean <- mask(alpha_mean, Global, inverse=FALSE, updatevalue=NA)
plot(clipped_alpha_mean)

writeRaster(clipped_alpha_ca, paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/",filename,"_alpha_ca_cropped.tif"), overwrite=FALSE)
writeRaster(clipped_alpha_mean, paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/",filename,"_alpha_mean_cropped.tif"), overwrite=FALSE)
