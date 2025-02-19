### This script was annotated with the help of generative AI with internet access ###

########################################################
##################      Read  Me      ##################

# 1. Generate cleaned_occurrence.shp (vector) with this script

# 2. Input vector generated into QGIS

# 3. Use point sampling tool (plugin) 
#    Get values from null_mask 

# 4. Export:
#    occ rec vector -> speciesname_lohman.csv
#    point sampling tool output -> spzone_lohman.csv

########################################################
########################################################

### 1. Generate cleaned_occurrence.shp

# Load necessary libraries
library(maptools)
library(BBmisc)
library(dplyr)
library(terra)

# Set working directory
wd <- "/lustre1/g/sbs_bonebrake/Eugene/SDMin"
setwd(wd)


# Define truncation limits for the study area
Xmini <- 69 #left x-axis, originally 64
Xmaxi <- 161.6 #right x-axis, originally 166
Ymini <- -10 #lower y-axis, originally -13
Ymaxi <- 36 #upper y-axis, originally 37

#Import species data
data <- read.csv("Occurrence Records of Tropical Asian Butterflies-1970-2024_12Feb2025EJ.csv")

# Standardize column names for latitude and longitude
colnames(data) <- gsub("decimalLat","decimalLatitude",colnames(data))
colnames(data) <- gsub("decimalLon","decimalLongitude",colnames(data))

# Convert latitude and longitude to numeric
data$decimalLatitude <- as.numeric(data$decimalLatitude)
data$decimalLongitude <- as.numeric(data$decimalLongitude)

# Prepare data by filtering based on defined geographical limits
data <- subset(data,decimalLatitude>=Ymini)
data <- subset(data,decimalLatitude<=Ymaxi)
data <- subset(data,decimalLongitude>=Xmini)
data <- subset(data,decimalLongitude<=Xmaxi)

# Convert year to numeric and filter out records before 1970
data$year <- as.numeric(data$year)
data <- data[! c(1:nrow(data)) %in% which(data$year<1970),]

# Correct species names (R does not like "-"s)
data[data$lohman_final_genus_species == "Polygonia c-album","lohman_final_genus_species"] <- "Polygonia c album"
data[data$lohman_final_genus_species == "Polygonia c-aureum","lohman_final_genus_species"] <- "Polygonia c aureum"

# Do the analysis only for species with >= 10 records
Count <- data %>% count(lohman_final_genus_species) #Count number of records, Feeds data into count(), saves result in count
Count <- subset(Count, n > 9) #Threshold was 25, changed to 10 for more species

# Clean the data to include only species with sufficient records
Cleaned.data <- data[which(data$lohman_final_genus_species %in% Count$lohman_final_genus_species),]

# Create a sorted list of unique species names
Species.list <- sort(unique(Cleaned.data$lohman_final_genus_species)) #Alphabetically sort species name

# Further prepare data for analysis
SDM_data <- Cleaned.data[!is.na(Cleaned.data$decimalLatitude),] #Remove negative/NA data
myRespXY <- SDM_data[,c("lohman_final_genus_species","decimalLatitude","decimalLongitude")] #Draws out latitude and longitude for each records
Projected.co <- terra::vect(myRespXY, geom=c("decimalLongitude","decimalLatitude"),crs="epsg:4326") #Converts to point file and give CRS

writeVector(Projected.co, "/lustre1/g/sbs_bonebrake/Eugene/SDMout/cleanedoccLohman/cleaned_occurrence.shp", filetype="ESRI Shapefile", layer=NULL, insert=FALSE,
            overwrite=FALSE, options="ENCODING=UTF-8")


### [QGIS] ###
### 2. Input vector generated into QGIS
### 3. Use point sampling tool (plugin), get values from null_mask
### 4. Generate speciesname_lohman.csv and spzone_lohman.csv
### [QGIS] ###


### 5. Get speciesname_lohman.csv and spzone_lohman.csv

# Set working directory
wd <- "/lustre1/g/sbs_bonebrake/Eugene/SDMin/landmaskdata_2025"
setwd(wd)

# Read species name and zone data from CSV files
recname <- read.csv("speciesname_lohman.csv")
reczone <- read.csv("spzone_lohman.csv")

# Initialize a data frame to store occurrence zones
columns <- c("Species","Zone") 
occzones <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(occzones) <- columns

# Combine species names and zones into a single data frame
occzones <- rbind(occzones, data.frame("Species"=recname[,1], 
                                       "Zone"=reczone[,1]))

# Initialize a data frame for recording species presence in different zones 
columns <- c("Species","borneo","borneoperc","eastasia","eastasiaperc","india","indiaperc",
             "indonesia","indonesiaperc","japan","japanperc","malaysia","malaysiaperc",
             "newguinea","newguineaperc","oceania","oceaniaperc","philippines","philippinesperc",
             "seasia","seasiaperc","taiwan","taiwanperc") 
spzones <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(spzones) <- columns

# Loop through each species to record occurrence in different zones
for(m in 1:length(Species.list)) {
  sp_data <- subset(occzones, Species == Species.list[m]) # Get data for the current species
  
  total <- nrow(sp_data) # Total occurrences for the species
  
  # Calculate occurrences and percentages for each zone
  borneo <- nrow(subset(sp_data,Zone == "borneo"))
  borneoperc <- borneo/total
  
  eastasia <- nrow(subset(sp_data,Zone == "eastasia"))
  eastasiaperc <- eastasia/total
  
  india <- nrow(subset(sp_data,Zone == "india"))
  indiaperc <- india/total
  
  indonesia <- nrow(subset(sp_data,Zone == "indonesia"))
  indonesiaperc <- indonesia/total
  
  japan <- nrow(subset(sp_data,Zone == "japan"))
  japanperc <- japan/total
  
  malaysia <- nrow(subset(sp_data,Zone == "malaysia"))
  malaysiaperc <- malaysia/total
  
  oceania <- nrow(subset(sp_data,Zone == "oceania"))
  oceaniaperc <- oceania/total
  
  newguinea <- nrow(subset(sp_data,Zone == "newguinea"))
  newguineaperc <- newguinea/total
  
  philippines <- nrow(subset(sp_data,Zone == "philippines"))
  philippinesperc <- philippines/total
  
  seasia <- nrow(subset(sp_data,Zone == "seasia"))
  seasiaperc <- seasia/total
  
  taiwan <- nrow(subset(sp_data,Zone == "taiwan"))
  taiwanperc <- taiwan/total
  
  # Record data to data frame
  spzones <- rbind(spzones, data.frame("Species"=Species.list[[m]], 
                                       "borneo"=borneo,
                                       "borneoperc"=borneoperc,
                                       "eastasia"=eastasia,
                                       "eastasiaperc"=eastasiaperc,
                                       "india"=india,
                                       "indiaperc"=indiaperc,
                                       "indonesia"=indonesia,
                                       "indonesiaperc"=indonesiaperc,
                                       "japan"=japan,
                                       "japanperc"=japanperc,
                                       "malaysia"=malaysia,
                                       "malaysiaperc"=malaysiaperc,
                                       "newguinea"=newguinea,
                                       "newguineaperc"=newguineaperc,
                                       "oceania"=oceania,
                                       "oceaniaperc"=oceaniaperc,
                                       "philippines"=philippines,
                                       "philippinesperc"=philippinesperc,
                                       "seasia"=seasia,
                                       "seasiaperc"=seasiaperc,
                                       "taiwan"=taiwan,
                                       "taiwanperc"=taiwanperc))
}

# Initialize a boolean data frame for zones each single species has been observed in
columns <- c("Species","borneo","eastasia","india","indonesia","japan","malaysia",
             "newguinea","oceania","philippines","seasia","taiwan") 
boolspzone <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(boolspzone) <- columns

# Loop through each row to determine presence in each zone
for(i in 1:nrow(data)){
  
  Species <- data$Species[i]
  
  # Species considered present and established in zones containing >= 1% of all records of that species
  if(data$borneoperc[i]>=0.01){
    borneo="TRUE"
  }else{
    borneo="FALSE"
  }
  
  if(data$eastasiaperc[i]>=0.01){
    eastasia="TRUE"
  }else{
    eastasia="FALSE"
  }
  
  if(data$indiaperc[i]>=0.01){
    india="TRUE"
  }else{
    india="FALSE"
  }
  
  if(data$indonesiaperc[i]>=0.01){
    indonesia="TRUE"
  }else{
    indonesia="FALSE"
  }
  
  if(data$japanperc[i]>=0.01){
    japan="TRUE"
  }else{
    japan="FALSE"
  }
  
  if(data$malaysiaperc[i]>=0.01){
    malaysia="TRUE"
  }else{
    malaysia="FALSE"
  }
  
  if(data$newguineaperc[i]>=0.01){
    newguinea="TRUE"
  }else{
    newguinea="FALSE"
  }
  
  if(data$oceaniaperc[i]>=0.01){
    oceania="TRUE"
  }else{
    oceania="FALSE"
  }
  
  if(data$philippinesperc[i]>=0.01){
    philippines="TRUE"
  }else{
    philippines="FALSE"
  }
  
  if(data$seasiaperc[i]>=0.01){
    seasia="TRUE"
  }else{
    seasia="FALSE"
  }
  
  if(data$taiwanperc[i]>=0.01){
    taiwan="TRUE"
  }else{
    taiwan="FALSE"
  }
  
  # Record species presences in data frame
  boolspzone <- rbind(boolspzone, data.frame("Species"=Species, 
                                             "borneo"=borneo,
                                             "eastasia"=eastasia,
                                             "india"=india,
                                             "indonesia"=indonesia,
                                             "japan"=japan,
                                             "malaysia"=malaysia,
                                             "newguinea"=newguinea,
                                             "oceania"=oceania,
                                             "philippines"=philippines,
                                             "seasia"=seasia,
                                             "taiwan"=taiwan))
  
}

# Save result to a CSV file
write.csv(boolspzone, "Species_zones_lohman_bool.csv")


