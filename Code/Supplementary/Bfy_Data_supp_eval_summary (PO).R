### This script was annotated with the help of generative AI with internet access ###

###### Import .csv files documenting evaluation data ######

# Define the directory containing the evaluation data files
filedir <- "biomod_4.2.4_datapp_2025"
basepaths <- paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_", filedir)

# List all finished species files in the specified directory
sp.templist <- list.files(basepaths) # Finished species list

# Format species names
sptemplist <- gsub("\\.", " ", sp.templist) # Or gsub('[.]', '', species.list.out0)

# Construct paths to the evaluation CSV files for each species
poss <- paste0(basepaths,"/",sp.templist,"/",sptemplist,"_po_evaluation.csv")

# Check if the evaluation files exist
tf <- file.exists(poss)

# Filter the species lists to include only those with existing evaluation files
sp.list <- sp.templist[tf]
splist <- sptemplist[tf]

# Select the ensemble model to examine
algoselect = "EMmean" # Options: "EMca" or "EMmean"

# Read the files documenting evaluation scores and filter by the selected algorithm
EvalPO <- read.csv(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filedir,"/",sp.list[[1]],"/",splist[[1]],"_po_evaluation.csv"))
EvalPO <- subset(EvalPO, algo==algoselect)

# Loop through the remaining species to read and combine their evaluation data
for(i in 2:length(sp.list)) {
  new <- read.csv(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filedir,"/",sp.list[[i]],"/",splist[[i]],"_po_evaluation.csv"))
  new <- subset(new, algo==algoselect)
  EvalPO <- rbind(EvalPO, new)}


###### Extract TSS & Boyce index results separately ######

# Extract TSS evaluation metrics
TSS <- subset(EvalPO, metric.eval=="TSS")
TSS.mean <- mean(TSS$evaluation, na.rm=TRUE)
TSS.sd <- sd(TSS$evaluation, na.rm=TRUE)

# Extract Boyce index evaluation metrics
BOYCE <- subset(EvalPO, metric.eval=="BOYCE")
BOYCE.mean <- mean(BOYCE$evaluation, na.rm=TRUE)
BOYCE.sd <- sd(BOYCE$evaluation, na.rm=TRUE)

# Show the results
TSS.mean # Mean TSS among ensemble models (all species)
TSS.sd # Standard deviation of TSS among ensemble models (all species)
BOYCE.mean # Mean Boyce index among ensemble models (all species)
BOYCE.sd # Standard deviation of Boyce index among ensemble models (all species)

