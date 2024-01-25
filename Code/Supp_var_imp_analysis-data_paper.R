##############################################################################
library(dplyr)

wd <- "/lustre1/g/sbs_bonebrake/Eugene/SDMin"
setwd(wd)

#select truncation
Xmini <- 69 #left x-axis, originally 64
Xmaxi <- 161.6 #right x-axis, originally 166
Ymini <- -10 #lower y-axis, originally -13
Ymaxi <- 36 #upper y-axis, originally 37

#Import species data
data <- read.csv("Occurrence_records_16Aug2023_lohman_EJ.csv")

#Truncate data
data <- subset(data,decimalLatitude>=Ymini)
data <- subset(data,decimalLatitude<=Ymaxi)
data <- subset(data,decimalLongitude>=Xmini)
data <- subset(data,decimalLongitude<=Xmaxi)

data$year <- as.numeric(data$year)
data <- subset(data,year>=1970)
data[data$lohman_final_genus_species == "Polygonia c-album","lohman_final_genus_species"] <- "Polygonia c album"
data[data$lohman_final_genus_species == "Polygonia c-aureum","lohman_final_genus_species"] <- "Polygonia c aureum"

Count <- data %>% count(lohman_final_genus_species) #Count number of records, Feeds data into count(), saves result in count
Count <- subset(Count, n > 9) #Threshold was 25, changed to 10 for more species

Cleaned.data <- data[which(data$lohman_final_genus_species %in% Count$lohman_final_genus_species),] #Returns  position/index of the value which satisfies given condition

Species.list <- unique(Cleaned.data$lohman_final_genus_species) #Enters species name

### Cleaning Species.list
Species.list[Species.list == "Polygonia c-album"] <- "Polygonia c album" #R cannot handle hyphen, have to change to underscore 
Species.list[Species.list == "Polygonia c-aureum"] <- "Polygonia c aureum"
Species.list <- Species.list[ !Species.list == 'not present'] # Remove the "not present" value in Species.list
Species.list <- sort(Species.list) # Alphabetically sort list

#Get list of species with SDM output
filepath <- "biomod_4.2.4_varimp_realnewfull"

finished <- list.files(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/species_list_", filepath))

finishedsp <- sub("_Finally_Finished.txt", "", finished)
bool1 <- Species.list %in% finishedsp
done <- Species.list[bool1]
length(done)

g.done <- gsub(" ", ".", done)


########## Use mean of species to generate data, list out by variable ###########

#Combine all var imp data as a table (1 model per row)
varimp <- read.csv(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_", filepath, "/",gsub("_",".",g.done[[1]]),"/",done[[1]],"_formal_models_variables_importance.csv"))
for(i in 2:length(done)) {
  new <- read.csv(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_", filepath, "/",gsub("_",".",g.done[[i]]),"/",done[[i]],"_formal_models_variables_importance.csv"))
  varimp <- rbind(varimp, new)
}
species.name <- gsub("_.*","",varimp$full.name)
varimp$full.name <- species.name

#Get list of species and variables
listed.sp <- unique(species.name)
var.list <- unique(varimp$expl.var)

#Create data frame to document var imp
columns <- c("Species","Variable","SPmean") 
varimp.weighted <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(varimp.weighted) <- columns

#Get mean variable importance for each species
for(x in 1:length(listed.sp)) { sp.data <- subset(varimp, varimp$full.name==listed.sp[[x]])
    for(y in 1:length(var.list)) {
        sp.vardata <- subset(sp.data, sp.data$expl.var==var.list[[y]])
        SPmean <- mean(sp.vardata$var.imp, na.rm=TRUE)
        varimp.weighted <- rbind(varimp.weighted, data.frame("Species"=listed.sp[[x]], "Variable"=var.list[[y]], "SPmean"=SPmean))
}}

columns <- c("Variable","Importance.mean", "Importance.median", "Importance.sd") 
varimp.breakdown <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(varimp.breakdown) <- columns

#Get mean, sd, median of var importance among all species for each variables
for(i in 1:length(var.list)) {
  data <- subset(varimp.weighted, varimp.weighted$Variable==var.list[[i]])
  Imp.mean <- mean(data$SPmean, na.rm=TRUE)
  Imp.sd <- sd(data$SPmean, na.rm=TRUE)
  Importance.median <- median(data$SPmean, na.rm=TRUE)
  varimp.breakdown = rbind(varimp.breakdown, data.frame("Variable"=var.list[[i]], "Importance.mean"=Imp.mean, "Importance.median"=Importance.median, "Importance.sd"=Imp.sd))
}

varimp.breakdown

#Save var imp data as .csv file
write.csv(varimp.breakdown, file.path(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/", filepath, "_variables_importance.csv")))
