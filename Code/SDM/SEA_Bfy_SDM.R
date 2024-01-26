library(biomod2)
library(maptools)
library(BBmisc)
library(dplyr)
library(terra)

###select truncation

Xmini <- 69 #left x-axis, originally 64
Xmaxi <- 161.6 #right x-axis, originally 166
Ymini <- -10 #lower y-axis, originally -13
Ymaxi <- 36 #upper y-axis, originally 37


###Creating output directory in designated location

filename <- paste0(Sys.Date(),"_SDM")

model.files.dir <- paste0("results_", filename)
dir.create(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/",model.files.dir)) 
species.list.dir <- paste0("species_list_", filename)
dir.create(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/",species.list.dir))

#Setup working directory (Data source)
wd <- "/lustre1/g/sbs_bonebrake/Eugene/SDMin"
setwd(wd)


###Import, format and extract data, produce alphabetically sorted list of species to be modeled

#Import species data
data <- read.csv("Occurrence_records_16Aug2023_lohman_EJ.csv")

#Limit data
Xmini <- 69 #left x-axis, originally 64
Xmaxi <- 161.6 #right x-axis, originally 166
Ymini <- -10 #lower y-axis, originally -13
Ymaxi <- 36 #upper y-axis, originally 37

#Truncate data
data <- subset(data,decimalLatitude>=Ymini)
data <- subset(data,decimalLatitude<=Ymaxi)
data <- subset(data,decimalLongitude>=Xmini)
data <- subset(data,decimalLongitude<=Xmaxi)

data$year <- as.numeric(data$year)
data <- subset(data,year>=1970)
data[data$lohman_final_genus_species == "Polygonia c-album","lohman_final_genus_species"] <- "Polygonia c album"
data[data$lohman_final_genus_species == "Polygonia c-aureum","lohman_final_genus_species"] <- "Polygonia c aureum"

Count <- data %>% count(lohman_final_genus_species) #Count number of records
Count <- subset(Count, n > 9) #Only species with 10 or more records are analyzed
Cleaned.data <- data[which(data$lohman_final_genus_species %in% Count$lohman_final_genus_species),] #Returns position/index of the value which satisfies given condition
Species.list <- unique(Cleaned.data$lohman_final_genus_species) #Create list of species with at least 10 records

#Cleaning Species.list
Species.list[Species.list == "Polygonia c-album"] <- "Polygonia c_album" #R cannot handle hyphen
Species.list[Species.list == "Polygonia c-aureum"] <- "Polygonia c_aureum"
Species.list <- Species.list[ !Species.list == 'not present'] # Remove the "not present" value (if any) from Species.list
Species.list <- sort(Species.list) # Alphabetically sort list


###Load and prepare study area map

Global1 <- terra::vect("world.shp")
Global <- terra::crop(Global1,ext(Xmini,Xmaxi,Ymini,Ymaxi)) #Crop raster file to define study area
Global <- terra::project(Global,"epsg:6933") 
r.rast <- terra::rast(Global, resolution=c(10000,10000)) #Creates raster with 10km resolution, might look strange if we plot it out but doesn't matter for our purpose
r <- terra::rasterize(Global, r.rast)


###Set output destination and file path of current climate files

save_wd <- paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/",model.files.dir) #SDM output result file destination
bioclimsource <- "worldclim"
current_wd <- paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMin/",bioclimsource) #Current climate data file location


###Select and prepare desired bioclim variables for modeling

file.current.list <-list.files(current_wd)
bio.seq <- grepl("01|04|05|06|12|13|14|15",file.current.list) 
file.current.list <- file.current.list[bio.seq]
current.climate <- c(paste0(current_wd,"/",file.current.list)) # Stacks all previously selected variables
current.climate <- terra::rast(current.climate)
current.climate <- terra::crop(current.climate,ext(Xmini,Xmaxi,Ymini,Ymaxi))
current.climate <- terra::project(current.climate,"epsg:6933") 
current.climate <- terra::resample(current.climate, r, method="average")


###Select and prepare other environmental variables for modeling

NDVImean <- terra::rast(paste0(current_wd,"/NDVImean.tif"))
names(NDVImean) <- "NDVImean"
NDVImean <- terra::project(NDVImean,"epsg:4326")
NDVImean <- terra::crop(NDVImean,ext(Xmini,Xmaxi,Ymini,Ymaxi))
NDVImean <- terra::project(NDVImean,"epsg:6933")
NDVImean <- terra::resample(NDVImean, r, method="average")

ph <- terra::rast(paste0(current_wd,"/Soilgrids_phh2o.tif"))
names(ph) <- "Soilgrids_soil_pH"
ph <- terra::project(ph,"epsg:4326")
ph <- terra::crop(ph,ext(Xmini,Xmaxi,Ymini,Ymaxi))
ph <- terra::project(ph,"epsg:6933")
ph <- terra::resample(ph, r, method="average")

SOC <- terra::rast(paste0(current_wd,"/Soilgrids_soc.tif"))
names(SOC) <- "Soilgrids_soil_organic_carbon"
SOC <- terra::project(SOC,"epsg:4326")
SOC <- terra::crop(SOC,ext(Xmini,Xmaxi,Ymini,Ymaxi))
SOC <- terra::project(SOC,"epsg:6933")
SOC <- terra::resample(SOC, r, method="average")

NITRO <- terra::rast(paste0(current_wd,"/Soilgrids_nitrogen.tif"))
names(NITRO) <- "Soilgrids_total_nitrogen"
NITRO <- terra::project(NITRO,"epsg:4326")
NITRO <- terra::crop(NITRO,ext(Xmini,Xmaxi,Ymini,Ymaxi))
NITRO <- terra::project(NITRO,"epsg:6933")
NITRO <- terra::resample(NITRO, r, method="average")

c.height <- terra::rast(paste0(current_wd,"/Cheight.tif"))
names(c.height) <- "Canopy_Height"
c.height <- terra::project(c.height,"epsg:4326")
c.height <- terra::crop(c.height,ext(Xmini,Xmaxi,Ymini,Ymaxi))
c.height <- terra::project(c.height,"epsg:6933")
c.height <- terra::resample(c.height, r, method="average")

#Stack all layers of environmental variables
current.climate <- c(current.climate, NDVImean, c.height, ph, SOC, NITRO)


###Create bias mask (for PA records generation)

library(MASS)

SDMdata <- Cleaned.data[!is.na(Cleaned.data$decimalLatitude),] #Remove negative/NA data
SDMdataXY <- SDMdata[,c("decimalLongitude","decimalLatitude")] #Draws out latitude and longitude for each records
alldata <- terra::vect(SDMdataXY, geom=c("decimalLongitude","decimalLatitude"), crs="epsg:4326") #Converts to point file and give CRS
alldata.proj <- terra::project(alldata,"epsg:6933")
alldata.proj.raster <- terra::rasterize(alldata.proj, r, vals=1) #Use base raster to make raster which shows "1" only at cells with occurance data, and "0" at cells with no data 
alldata.proj.raster[alldata.proj.raster>=1] <- 1 #Ensure the maximum cell value is limited to "1" in every 10km x 10km grid square 
alldata.rast <- alldata.proj.raster

alldata.rast.reproj <- terra::project(alldata.rast,"epsg:4326")

df <- crds(alldata.rast.reproj, na.rm=TRUE) 

nrow(alldata.rast.reproj)

dens <- kde2d(df[,1], df[,2], n = c(nrow(alldata.rast.reproj), ncol(alldata.rast.reproj))) #Perform two-dimensional kernel density estimation
dens.ras <- raster::raster(dens)
dens.ras <- terra::rast(dens.ras)
dens.ras <- terra::project(dens.ras,"epsg:6933")
dens.ras <- terra::resample(dens.ras, r, method="average")


###Clean environment for optimized RAM usage, redirect wd to desired output location

rm(NDVImean, c.height, ph, SOC, NITRO)
rm(SDMdata, SDMdataXY, alldata, alldata.proj, alldata.proj.raster, alldata.rast, dens, df)
rm(Count, data, t.Cleaned.data)
gc(full = TRUE)

setwd(save_wd) #Ensure saved to right location


###Function to run SDMs individually for each species, designed to support parallel computing as long as your device have multiple CPU cores and sufficient RAM

MyBiomodfun <- function(m) {
  message("Species", Species.list[[m]]) #Let me know which SDM the program is running
  error.msg <- NULL 
  SDM_data <- subset(Cleaned.data,lohman_final_genus_species == Species.list[m]) #(Subsets/draws out data of single species for analysis)
  SDM_data <- SDM_data[!is.na(SDM_data$decimalLatitude),] #Remove negative/NA data
  myRespXY <- SDM_data[,c("decimalLongitude","decimalLatitude")] #Draws out latitude and longitude for each records
  Projected.co <- terra::vect(myRespXY, geom=c("decimalLongitude","decimalLatitude"), crs="epsg:4326") #Converts to point file and give CRS
  T.Projected.co <- terra::project(Projected.co,"epsg:6933")
  T.Projected.co.raster <- terra::rasterize(T.Projected.co, r, vals=1) #Use base raster to make raster which has data only at occurrence point
  T.Projected.co.raster[T.Projected.co.raster>=1] <- 1 #Limit the max. number of data points to 1 in every 10km x 10km grid square 
  
  #Create a dataframe to record 10km x 10km coordinates of grid squares with >=1 data points
  df <- crds(T.Projected.co.raster, na.rm=TRUE) 
  
  if (nrow(df) > 9) { #Check if we have more than 10 EFFECTIVE records (10 cells with occurrence records at 10km x 10km resolution)
    myRespXY <- crds(T.Projected.co.raster, na.rm=TRUE)
    vectXY <- terra::vect(myRespXY, crs="epsg:6933")
    
    message("modeling current climate") #Starting of the SDM part
    
    #Options for SDM construction, refer to biomod2 package manual for details.
    myBiomodOption <- BIOMOD_ModelingOptions(
      
      GAM = list(algo = 'GAM_mgcv', 
                 type = 's_smoother',
                 k = -1, 
                 interaction.level = 0, 
                 myFormula = NULL,
                 family = binomial(link = 'logit'),
                 method = 'REML', 
                 optimizer = c('outer','newton'), 
                 select = TRUE, 
                 knots = NULL, 
                 paraPen = NULL), 
      
      ANN = list(NbCV = 5, 
                 size = NULL, 
                 decay = NULL, 
                 maxit = 200), 
      
      MARS = list(type = 'simple', 
                  interaction.level = 0, 
                  myFormula = NULL,
                  nk = NULL, 
                  penalty = 2, 
                  thresh = 0.001, 
                  nprune = NULL, 
                  pmethod = 'backward'), 
      
      CTA = list(control= list(minsplit = 5, 
                               minbucket = round(5/3), 
                               xval = 5, 
                               cp = 0.001, 
                               maxdepth = 25)),
      
      RF = list(do.classif = TRUE, 
                ntree = 2000, 
                mtry = 'default', 
                nodesize = 5), 
      
      XGBOOST = list(max.depth = 5, 
                     eta = 0.01, 
                     nrounds = 1000, 
                     objective = "binary:logistic", 
                     nthread = 1) 
    ) 
    
    #Partition and reserve model evaluation data
    eval.no <- floor(nrow(df)*0.10) #Calculates the number of occurrence data entries to be reserved (10% of all data available), actual number of entries rounded down to the nearest integer
    eval.cols <- sample(nrow(df), eval.no) #Randomly select occurrence data entries to be reserved
    
    #Splits all available occurrence data into training data(90%) and evaluation data(10%)
    base.vectXY <- terra::vect(df, crs="epsg:6933")
    eval.df <- df[eval.cols,]
    mod.df <- df[-eval.cols,]
    vectXY <- terra::vect(mod.df, crs="epsg:6933") #Training data
    eval.vectXY <- erase(base.vectXY, vectXY) #Evaluation data
    
    #Assign value 1 to presence
    values(eval.vectXY) <- data.frame(rec=rep(1, nrow(crds(eval.vectXY))))
    values(base.vectXY) <- data.frame(rec=rep(1, nrow(df)))
    
    #Generate pseudo absence (PA) records for model training
    PA.dens <- mask(dens.ras, T.Projected.co.raster, inverse=FALSE, maskvalues=1, updatevalue=NA) #Exclude cells with occurrence data from possible PA points
    PA.dens <- mask(dens.ras, r, inverse=TRUE, maskvalues=1, updatevalue=NA) #Exclude non-land cells from possible PA points
    PA.dens.r <- raster::raster(PA.dens)
    PAcoord <- raster::rasterToPoints(PA.dens.r, spatial = FALSE) #Need rasterToPoints() function to extract coordinates for each cell along with cell value
    
    PAsets <- 5 #Indicate number of sets of PA records we want to generate, should be equal to CV.nb.rep in the BIOMOD_Modeling() function, which will be used later to run the SDMs
    PAsamplesets <- PAsets+2
    
    ###Generate Model Training PA
    PAdf <- PAcoord[sample(seq(1:nrow(PAcoord)), size = nrow(df)*PAsamplesets, replace = FALSE, prob = PAcoord[, "layer"]), 1:2]
    myRespXY <- rbind(df,PAdf)
    
    ###Generate Evaluation PA
    Ev_PAdf_plus <- PAcoord[sample(seq(1:nrow(PAcoord)), size = eval.no+1, replace = FALSE, prob = PAcoord[, "layer"]), 1:2]
    Ev_PAdf_plus <- as.data.frame(Ev_PAdf_plus)
    Ev_PAdf <- Ev_PAdf_plus[c(1:eval.no),]
    #Assign value 0 to eval PA, because biomod2 currently only does not accept "NA" as absence input.
    eval.vectpa <- terra::vect(Ev_PAdf, geom=c("x", "y"), crs="epsg:6933")
    na <- data.frame(rec=rep(0, nrow(Ev_PAdf)))
    values(eval.vectpa) <- na
    
    #Combine presence and PA points for model evaluation
    eval.vectXYfull <- rbind(eval.vectXY, eval.vectpa)
    
    #Sample PA points and format the sampled points
    presence.rep <- seq(1,1, length.out=nrow(df))
    PA.rep <- seq(0,0, length.out=nrow(df)*PAsamplesets) # Since the coming PA lists randomly sample from available PA points, need more PA than actually needed
    myResp = c(presence.rep, PA.rep)
    myResp.PA <- ifelse(myResp == 1, 1, NA)
    myPAtable <- data.frame(PA1 = ifelse(myResp == 1, TRUE, FALSE),
                            PA2 = ifelse(myResp == 1, TRUE, FALSE),
                            PA3 = ifelse(myResp == 1, TRUE, FALSE),
                            PA4 = ifelse(myResp == 1, TRUE, FALSE),
                            PA5 = ifelse(myResp == 1, TRUE, FALSE)
    )
    myPAtable[sample(which(myPAtable[, 1] == FALSE), nrow(df), replace = FALSE), 1] = TRUE
    myPAtable[sample(which(myPAtable[, 2] == FALSE), nrow(df), replace = FALSE), 2] = TRUE
    myPAtable[sample(which(myPAtable[, 3] == FALSE), nrow(df), replace = FALSE), 3] = TRUE
    myPAtable[sample(which(myPAtable[, 4] == FALSE), nrow(df), replace = FALSE), 4] = TRUE
    myPAtable[sample(which(myPAtable[, 5] == FALSE), nrow(df), replace = FALSE), 5] = TRUE
    
    ### Further format the occurrence and PA data for biomod2 to read
    myBiomodData <- BIOMOD_FormatingData(resp.var = myResp.PA, 
                                         expl.var = current.climate, 
                                         resp.xy = myRespXY,
                                         resp.name = Species.list[m], 
                                         PA.strategy = "user.defined", 
                                         PA.user.table = myPAtable,
                                         eval.resp.var = eval.vectXYfull,
                                         eval.expl.var = current.climate,
                                         eval.resp.xy = eval.df,
                                         filter.raster = TRUE) 
    
    ### Training SDMs
    myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
                                        models = c("GLM","CTA","MARS","MAXNET","XGBOOST"), #Indicate algorithms
                                        bm.options = myBiomodOption,
                                        CV.strategy = "random", #How to partition data into training set and validation set for cross-validation
                                        CV.nb.rep = 5, #Number of times to run cross-validation
                                        CV.perc = 0.889, #Percentage of data in training set, 80% of all occurrence data available, given by 0.8/0.9=0.889
                                        prevalence = 0.5, #Default setting
                                        metric.eval = c("TSS"), #Model evaluation method(s)
                                        scale.models = TRUE, #models predictions scaled with a binomial GLM or not 
                                        CV.do.full.models = FALSE, 
                                        modeling.id = paste("FirstModeling",sep=""), #Name of model (Species name_Firstmodeling)
                                        var.import = 5, #number of permutations to be done for each variable to estimate variable importance
                                        seed.val = 66) #Can provide seed for better reproducibility 
    
    #csv file for performance of each models built
    write.csv(get_evaluations(myBiomodModelOut),
              file.path(paste(gsub(" |_",".",Species.list[m]),"/",Species.list[m],"_formal_models_evaluation.csv", sep="")))
    #csv file for importance of each variables
    write.csv(get_variables_importance(myBiomodModelOut),
              file.path(paste(gsub(" |_",".",Species.list[m]),"/",Species.list[m],"_formal_models_variables_importance.csv", sep="")))
    
    #Merging all models into one model (ensemble model)
    myBiomodEM <- try(BIOMOD_EnsembleModeling( #Used try() since it is possible to have no model with TSS > your threshold, which can result in error
      bm.mod = myBiomodModelOut,
      models.chosen = "all",
      em.by = "all", 
      metric.select = c("TSS"),
      metric.select.thresh = 0.7, #Only consider a model for averaging when its TSS>0.7
      em.algo = c('EMmean', 'EMca'),
      prob.ci.alpha = 0.05,
      prob.mean.weight.decay = "proportional")) 
    
    #Extract Boyce index and TSS for the ensemble model(s)
    myBiomodPO <- try(BIOMOD_PresenceOnly(bm.mod = myBiomodModelOut, #Used try() since it is possible to have no model with TSS > your threshold, which can result in error
                                          bm.em = myBiomodEM,
                                          bg.env = current.climate))
    try(write.csv(myBiomodPO,
                  file.path(paste(gsub(" |_",".",Species.list[m]),"/",Species.list[m],"_po_evaluation.csv", sep=""))), silent=TRUE)
    
    #Each model individually estimates habitat suitability
    myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                      new.env = current.climate, #Specify the environment to project SDM in
                                      proj.name = "Current",
                                      models.chosen = "all",
                                      metric.binary = "TSS",
                                      compress = "xz",
                                      build.clamping.mask = F,
                                      output.format = ".grd")
    
    #If we have 1 model in ensemble model, project the ensemble model as well
    if (length(myBiomodEM@em.computed) == 1) {
      myBiomodEF <- try(BIOMOD_EnsembleForecasting(
        bm.em = myBiomodEM,
        bm.proj = myBiomodProj,
        models.chosen = 'all',
        metric.binary = "TSS"))
      write.table("Finished", #For progress monitoring and gives info on whether the SDM run was successful
                  file.path(paste("/lustre1/g/sbs_bonebrake/Eugene/SDMout/",species.list.dir,"/",Species.list[m],"_Single_Model_Finished.txt", sep="")))} 
    
    #If we have multiple models merged as ensemble model(s), project the ensemble model(s) as well
    if (length(myBiomodEM@em.computed) > 1) { #Combine 15 models to do projection
      myBiomodEF <- BIOMOD_EnsembleForecasting(
        bm.em = myBiomodEM,
        bm.proj = myBiomodProj,
        models.chosen = 'all',
        metric.binary = "TSS")
      
      write.table("Finished", #For progress monitoring and gives info on whether the SDM run was successful
                  file.path(paste("/lustre1/g/sbs_bonebrake/Eugene/SDMout/",species.list.dir,"/",Species.list[m],"_Finally_Finished.txt", sep="")))
      
    } else {write.table("Failed", #For progress monitoring and gives info on whether the SDM run was successful
                        file.path(paste("/lustre1/g/sbs_bonebrake/Eugene/SDMout/",species.list.dir,"/",Species.list[m],"Model_Score_Low_EM_Fail.txt", sep="")))}
    
    write.table("Finished", #For progress monitoring and gives info on whether the SDM run was successful
                file.path(paste(gsub(" |_",".",Species.list[m]),"/",Species.list[m],"_Finished.txt", sep="")))
    
  } else {
    write.table("Failed", #For progress monitoring and gives info on whether the SDM run was successful
                file.path(paste("/lustre1/g/sbs_bonebrake/Eugene/SDMout/",species.list.dir,"/",Species.list[m],"_Insufficient_Records_Failed.txt", sep="")))
  }
}


###Now specify which species you would want to run SDM with, you can, of course, run all species with sufficient data using this:

m <- 1:1294

#OR... If you're interested in any particular species on the Species.list, you can use match() to obtain the location of your particular species on the Species.list.
#Note that match() is not exactly meant to be a search tool, so it's a good idea to copy the exact species name from the Species.list.
#For example, if you only want to model Papilio memnon:

species <- "Papilio memnon" 
m <- match(species, Species.list)

#"I can wait forever"
#  - If you have only one species to run, you can use the built-in method to run the MyBiomodfun function in the simplest manner
lapply(m, MyBiomodfun)

#Parallel computing -- UNLIMITED (computing) POWER!!!
#  - For running multiple species
#We've tried various methods to make our SDM modelling parallel, and mclapply is currently the easiest way to implement parallel computing (But Windows users beware, mclapply() function uses forked processes to parallelize things and therefore cannot be executed on Windows OS. We suggest following the guidelines of biomod2 and use the "snowfall" package, based on socket process, to implement parallel computing on Windows OS. Still, Mac and Linux users are recommended to use the mclapply() function since fork is more efficient than socket and requires less preparation work). Other methods (such as "future" package) also works.
#  - For normal (PC) users: Your PC/laptop might advertise a CPU capable of running 10+ threads, but your biggest obstacle is actually your RAM. mc.cores=4 will be a good start.
#  - For HPC users: From my experience a node with 180GB RAM can run this script with mc.cores=15 steadily.

parallel::mclapply(m, MyBiomodfun, mc.cores=15)





