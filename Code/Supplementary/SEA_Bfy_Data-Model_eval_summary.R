
###### Import .csv files documenting evaluation data ######

filedir <- "biomod_4.2.4_POevaltest"
sp.templist <- list.files(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filedir)) # Finished species list
sptemplist <- gsub("\\.", " ", sp.templist) # Or gsub('[.]', '', species.list.out0)
basepaths <- paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filedir)
unlist(sp.templist)

poss0 <- paste0(basepaths,"/",sp.templist,"/",sptemplist,"_po_evaluation.csv")
tf <- file.exists(poss0)
sp.list <- sp.templist[tf]
splist <- sptemplist[tf]

algoselect="EMca" # "EMca" or "EMmean": Select which ensemble model to examine
EvalPO <- read.csv(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filedir,"/",sp.list[[1]],"/",splist[[1]],"_po_evaluation.csv"))
EvalPO <- subset(EvalPO, algo==algoselect)

for(i in 2:length(sp.list)) {
  new <- read.csv(paste0("/lustre1/g/sbs_bonebrake/Eugene/SDMout/results_",filedir,"/",sp.list[[i]],"/",splist[[i]],"_po_evaluation.csv"))
  new <- subset(new, algo==algoselect)
  EvalPO <- rbind(EvalPO, new)}


###### Extract TSS & Boyce index results separately ######

TSS <- subset(EvalPO, metric.eval=="TSS")
TSS.mean <- mean(TSS$evaluation, na.rm=TRUE)
TSS.sd <- sd(TSS$evaluation, na.rm=TRUE)

BOYCE <- subset(EvalPO, metric.eval=="BOYCE")
BOYCE.mean <- mean(BOYCE$evaluation, na.rm=TRUE)
BOYCE.sd <- sd(BOYCE$evaluation, na.rm=TRUE)

TSS.mean # Mean TSS among ensemble models (all species)
TSS.sd # Standard deviation of TSS among ensemble models (all species)
BOYCE.mean # Mean Boyce index among ensemble models (all species)
BOYCE.sd # Standard deviation of Boyce index among ensemble models (all species)


