### This script was annotated with the help of generative AI with internet access ###

# Setup working directory and load data
wdpc <- "C:/Users/Eugene/Desktop"
setwd(wdpc)

data <- read.csv("spname_unharmonized_dataset.csv") 
changelist <- read.csv("Bonebrake_SEA_species_DJLvalidated4_EY.csv") #Make sure to download this file too
changelistPH <- read.csv("PH_species_DJLvalidated_DJL.csv") #Make sure to download this file too

# Merge change lists
# Initialize new columns in changelistPH for harmonization
changelistPH$Lohman_validated_subspecies <- NA
changelistPH$accepted_genus_species <- NA
changelistPH$Genus <- NA
changelistPH$specific_epithet <- NA
changelistPH$Notes <- NA

changelistPH <- changelistPH[,-1]
colnames(changelistPH) <- c("genus_species","Authorities","Lohman_validated_binomial","Lohman_validated_subspecies",
                            "accepted_genus_species","Genus","specific_epithet","Notes")
changelistPH <- changelistPH[,colnames(changelist)] # Reorder columns to match the original dataset

changelist <- rbind(changelist,changelistPH) # Append the PH change list to the main change list

# Create new columns in the main dataset for final genus and species names
data$lohman_final_genus_species <- NA
data$lohman_final_subspecies <- NA

# For each "genus_species" name (scientific name as in GBIF data), add our validated names in the new columns
for(m in 1:nrow(changelist)) {
  oldname <- changelist[m,]$genus_species
  newname <- changelist[m,]$Lohman_validated_binomial
  newsubname <- changelist[m,]$Lohman_validated_subspecies
  
  # Update the main dataset with the new names where the old name matches
  data[data$source_genus_species == oldname,"lohman_final_genus_species"] <- newname
  data[data$source_genus_species == oldname,"lohman_final_subspecies"] <- newsubname
}

# Initialize a new column for authorities
data$authorities <- NA

for(m in 1:nrow(changelist)) {
  original.spname <- changelist[m,]$genus_species
  lohman.authorities <- changelist[m,]$Authorities
  
  # Update the main dataset with the authorities where the species name matches
  data[data$source_genus_species == original.spname,"authorities"] <- lohman.authorities
}

# Save the harmonized dataset to a new CSV file
write.csv(data, "spname_harmonized_dataset.csv")

