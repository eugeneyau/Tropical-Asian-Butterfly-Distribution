### This script was annotated with the help of generative AI with internet access ###

library(dplyr)

### Prepare for analysis
# Define the file path for SDM output 
filepath <- "biomod_4.2.4_datapp_2025"
basepath <- paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filepath)
species.list.out <- list.files(basepath) # List out all species from the base path

BASEprojection <- "/proj_Current/proj_Current_"
poss <- paste0(basepath,"/",species.list.out,BASEprojection,species.list.out,"_ensemble_TSSbin.tif") # Full path of SDM-projected maps
tf <- file.exists(poss) # Check if the output files exist
sp.list <- species.list.out[tf] # Filter the species list to only those with existing files
sp_list <- gsub("\\."," ",sp.list) # Format species list

### Producing variable importance table by adding rows
# Read the variable importance data for the first species
varimp <- read.csv(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_", filepath, "/",sp.list[[1]],"/",sp_list[[1]],"_formal_models_variables_importance.csv"))

# Loop through the remaining species and append their variable importance data
for(i in 2:length(sp.list)) {
  new <- read.csv(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_", filepath, "/",sp.list[[i]],"/",sp_list[[i]],"_formal_models_variables_importance.csv"))
  varimp <- rbind(varimp, new)
}

# Format full.name column
species.name <- gsub("_.*","",varimp$full.name)
varimp$full.name <- species.name

### Analyze variable importance table
# Get mean variable importance for each species
listed.sp <- unique(species.name) # Get unique species names
var.list <- unique(varimp$expl.var) # Get unique explanatory variables

# Initialize a data frame to store weighted variable importance
columns <- c("Species","Variable","SPmean") 
varimp.weighted <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(varimp.weighted) <- columns

# Loop through each species to calculate mean variable importance
for(x in 1:length(listed.sp)) { 
  sp.data <- subset(varimp, varimp$full.name==listed.sp[[x]]) # Subset data for the current species
for(y in 1:length(var.list)) {
  sp.vardata <- subset(sp.data, sp.data$expl.var==var.list[[y]]) # Subset data for the current variable
  SPmean <- mean(sp.vardata$var.imp, na.rm=TRUE) # Calculate mean variable importance
  varimp.weighted <- rbind(varimp.weighted, data.frame("Species"=listed.sp[[x]], "Variable"=var.list[[y]], "SPmean"=SPmean))
}}


# Summarise variable importance for all species
columns <- c("Variable","Importance.mean", "Importance.median", "Importance.sd") 
varimp.breakdown <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(varimp.breakdown) <- columns

# Loop through each variable to calculate summary statistics
for(i in 1:length(var.list)) {
  data <- subset(varimp.weighted, varimp.weighted$Variable==var.list[[i]])
  Imp.mean <- mean(data$SPmean, na.rm=TRUE)
  Imp.sd <- sd(data$SPmean, na.rm=TRUE)
  Importance.median <- median(data$SPmean, na.rm=TRUE)
  varimp.breakdown <- rbind(varimp.breakdown, data.frame("Variable"=var.list[[i]], "Importance.mean"=Imp.mean, "Importance.median"=Importance.median, "Importance.sd"=Imp.sd))
}

# Display the final variable importance breakdown
varimp.breakdown 

