
# Setup working directory and load data
wdpc <- "C:/Users/Eugene/Desktop"
setwd(wdpc)

data <- read.csv("Occurrence_records_16Aug2023_EJ.csv") #Or your GBIF data
changelist <- read.csv("Bonebrake_SEA_species_DJLvalidated3.csv") #Make sure to download this file from our Github repository: Tropical-Asian-Butterfly-Distribution/Code/Files

# Create new columns
data$lohman_final_genus_species <- NA
data$lohman_final_subspecies <- NA
data$authorities <- NA

# For each "genus_species" name (scientific name as in GBIF data), add our validated names in the new columns
for(m in 1:nrow(changelist)) {
  oldname <- changelist[m,]$genus_species
  newname <- changelist[m,]$Lohman_validated_binomial
  newsubname <- changelist[m,]$Lohman_validated_subspecies
  lohman.authorities <- changelist[m,]$Authorities
  
  
  data[data$genus_species == oldname,"lohman_final_genus_species"] <- newname
  data[data$genus_species == oldname,"lohman_final_subspecies"] <- newsubname
  data[data$genus_species == original.spname,"authorities"] <- lohman.authorities
}

write.csv(data, "OccurrenceRecords_22Jan2024_EY.csv")

