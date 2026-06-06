
################## Get all GBIF names to harmonize ##################

library(dplyr)

# Setup working directory and load data
wdpc <- "C:/Users/Eugene/Desktop"
wdG15 <- "C:/Users/skywa/Desktop" #Laptop
setwd(wdG15)

changelist <- read.csv("Bonebrake_SEA_species_DJLvalidated5_YFL.csv") # Make sure to download this file from our Figshare repository (SDMsupp_files.zip)
changelistPH <- read.csv("PH_species_DJLvalidated_DJL.csv") # Make sure to download this file from our Figshare repository (SDMsupp_files.zip)
GBIFdata <- read.delim("0169126-240321170329656.csv") # Raw GBIF download 

#GBIFdata <- GBIFdata[,c("family","genus","species","coordinateUncertaintyInMeters","infraspecificEpithet","taxonRank","scientificName",
#                        "verbatimScientificName","verbatimScientificNameAuthorship","issue")]
cleanedGBIFdata <- subset(GBIFdata, coordinateUncertaintyInMeters <= 100000 | is.na(coordinateUncertaintyInMeters))
cleanedGBIFdata <- subset(cleanedGBIFdata, verbatimScientificName != "")
#subset_data <- cleanedGBIFdata[!(cleanedGBIFdata$taxonRank %in% c("GENUS","FAMILY")), ]
#subset_data <- cleanedGBIFdata[(cleanedGBIFdata$taxonRank %in% c("VARIETY","FORM")), ] # 45 records
#subset_data <- cleanedGBIFdata[(cleanedGBIFdata$taxonRank %in% c("UNRANKED")), ] # 7607 records with only GBIF-asssigned name and BOLD identifier

# Simple filtering of names to harmonize
cleanedGBIFdata$shortverbatimScientificName <- sub("^([^ ]+ [^ ]+).*", "\\1", cleanedGBIFdata$verbatimScientificName)

rawnames <- unique(cleanedGBIFdata$shortverbatimScientificName)
knownspnames <- unique(c(
  changelist$genus_species, 
  changelist$Lohman_validated_binomial, 
  changelistPH$source_genus_species, 
  changelistPH$lohman_genus_species))

present_names <- rawnames[rawnames %in% knownspnames]
missing_names <- rawnames[!rawnames %in% knownspnames]

sum(grepl("\\(|\\)", missing_names))
sum(!grepl("BOLD:", missing_names))

pendingnames <- missing_names[!grepl("BOLD:", missing_names)]

pending_harmonization <- cleanedGBIFdata[cleanedGBIFdata$shortverbatimScientificName %in% pendingnames, ]
pending_harmonization <- unique(pending_harmonization[, c("shortverbatimScientificName","verbatimScientificName","scientificName","species")])
pending_harmonization <- pending_harmonization[order(pending_harmonization$shortverbatimScientificName), ]

#write.csv(pending_harmonization, "pending_GBIF_harmonization.csv", row.names = FALSE, fileEncoding = "UTF-8")


################################################################################
################################## START HERE ##################################
################################################################################


################## Harmonize GBIF names ##################

# Setup working directory and load data
wdpc <- "C:/Users/Eugene/Desktop"
wdG15 <- "C:/Users/skywa/Desktop" # Laptop
setwd(wdG15)

changelistGBIF <- read.csv("GBIF_harmonization_YFL.csv") # Updated harmonization list, make sure to download this file from our GitHub repository (https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/tree/main/Version%201.1%20update)
GBIFdata <- read.delim("0169126-240321170329656.csv") # Raw GBIF download, make sure to download this file from our Figshare repository

GBIFdata <- subset(GBIFdata, coordinateUncertaintyInMeters <= 100000 | is.na(coordinateUncertaintyInMeters))

# Create new columns
GBIFdata$harmonized_genus_species <- ifelse(
  grepl("BOLD:", GBIFdata$verbatimScientificName),
  GBIFdata$species,
  sub("^([^ ]+ [^ ]+).*", "\\1", GBIFdata$verbatimScientificName) # Only keep the first two names (exclude subspecies name)
)
GBIFdata$harmonized_subspecies <- NA

# For each raw "verbatimScientificName" name in GBIF, add our validated names in the new columns
for(m in 1:nrow(changelistGBIF)) {
  oldname <- changelistGBIF[m,]$verbatimScientificName
  newname <- changelistGBIF[m,]$FinalNameInDatabase
  newsubname <- changelistGBIF[m,]$FinalSubspeciesInDatabase
  
  GBIFdata[GBIFdata$verbatimScientificName == oldname,"harmonized_genus_species"] <- newname
  GBIFdata[GBIFdata$verbatimScientificName == oldname,"harmonized_subspecies"] <- newsubname
}

# Reformat GBIF data to match our data set format
GBIFdata$recordIndex <- ""
GBIFdata$enteredBy <- "Eugene"
GBIFdata$reference <- "Derived dataset GBIF.org (20 May 2024). Filtered export of GBIF occurrence data. https://doi.org/10.15468/dd.nvw5wr"
GBIFdata$source_subfamily <- ""
GBIFdata$source_species <- ""
GBIFdata$COL_genus_species <- ""
GBIFdata$lohman_final_genus_species <- ""
GBIFdata$lohman_final_subspecies <- ""
GBIFdata$country <- ""
GBIFdata$locality_type <- ""
GBIFdata$how_coordinates_were_chosen <- ""
GBIFdata$enteredBy <- ""
GBIFdata$authorities <- ""

mapping <- c(
  recordIndex                   = "recordIndex",
  enteredBy                     = "enteredBy",
  reference                     = "reference",
  datasetKey                    = "datasetKey",
  source_family                 = "family",
  source_subfamily              = "source_subfamily",
  source_genus                  = "genus",
  source_species                = "source_species",
  source_genus_species          = "harmonized_genus_species",
  COL_genus_species             = "COL_genus_species",
  lohman_final_genus_species    = "lohman_final_genus_species",
  lohman_final_subspecies       = "lohman_final_subspecies",
  source_authority              = "verbatimScientificNameAuthorship",
  countryCode                   = "countryCode",
  country                       = "country",
  locality                      = "locality",
  locality_type                 = "locality_type",
  decimalLat                    = "decimalLatitude",
  decimalLon                    = "decimalLongitude",
  how_coordinates_were_chosen   = "how_coordinates_were_chosen",
  year                          = "year",
  coordinateUncertaintyInMeters = "coordinateUncertaintyInMeters",
  month                         = "month",
  basisOfRecord                 = "basisOfRecord",
  gbifID                        = "gbifID",
  licenseType                   = "license",
  licenseHolder                 = "rightsHolder",
  gbif_issue                    = "issue",
  authorities                   = "authorities"
)

GBIFdata <- GBIFdata[, mapping]
colnames(GBIFdata) <- names(mapping)
GBIFdata <- subset(GBIFdata, source_genus_species != "")

# Get list for species name harmonization
changelist <- read.csv("Bonebrake_SEA_species_DJLvalidated4mod1_EY_YFL.csv") # Make sure to download this file from our Figshare repository (SDMsupp_files.zip)
changelistPH <- read.csv("PH_species_DJLvalidated_DJL.csv") # Make sure to download this file from our Figshare repository (SDMsupp_files.zip)

changelistPH$Lohman_validated_subspecies <- ""
changelistPH$accepted_genus_species <- ""
changelistPH$Genus <- ""
changelistPH$specific_epithet <- ""
changelistPH$Notes <- ""

changelistPH <- changelistPH[,-1]
colnames(changelistPH) <- c("genus_species","Authorities","Lohman_validated_binomial","Lohman_validated_subspecies",
                            "accepted_genus_species","Genus","specific_epithet","Notes")
changelistPH <- changelistPH[,colnames(changelist)]

changelist <- rbind(changelist,changelistPH)

# Harmonize GBIF data using source_genus_species
for(m in 1:nrow(changelist)) {
  oldname <- changelist[m,]$genus_species
  newname <- changelist[m,]$Lohman_validated_binomial
  newsubname <- changelist[m,]$Lohman_validated_subspecies
  
  target_rows <- which(GBIFdata$source_genus_species == oldname | 
                         GBIFdata$source_genus_species == newname)
  
  if(length(target_rows) > 0) {
    GBIFdata[target_rows, "lohman_final_genus_species"] <- newname
    GBIFdata[target_rows, "lohman_final_subspecies"] <- newsubname}
}

# Merge GBIF data and our non-GBIF data
figsharedata <- read.csv("Occurrence Records of Tropical Asian Butterflies - 1970-2024_v03June2025EJ.csv") # Make sure to download this file from our Figshare repository
figsharedata <- subset(figsharedata, reference != "Derived dataset GBIF.org (20 May 2024). Filtered export of GBIF occurrence data. https://doi.org/10.15468/dd.nvw5wr")

bigdata<- rbind(figsharedata,GBIFdata)
rownames(bigdata) <- NULL
bigdata$recordIndex <- 1:nrow(bigdata)

write.csv(bigdata, "Occurrence Records of Tropical Asian Butterflies - 1970-2024 V2.csv", row.names = FALSE, fileEncoding = "UTF-8")





