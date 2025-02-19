# Setup working directory and load data
wdpc <- "C:/Users/Eugene/Desktop"
setwd(wdpc)

data <- read.csv("spname_harmonized_dataset.csv") 
changelist <- read.csv("clean_family_EJ.csv") #Make sure to download this file from our Github repository: Tropical-Asian-Butterfly-Distribution/Code/Files

# Update family names in the main dataset 
for(m in 1:nrow(changelist)) {
  oldname <- changelist[m,]$family.sp.list
  newname <- changelist[m,]$final_family_EJ

  data[data$lohman_final_genus_species == oldname & !(is.na(data$lohman_final_genus_species)),"source_family"] <- newname
}

# Save the updated dataset to a new CSV file
write.csv(data, "updated_family_spname_harmonized_dataset.csv")