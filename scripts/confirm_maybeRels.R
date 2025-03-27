### CONFIRM IDENTIY OF MAYBE RELATIVES 
# goal: can I use size to determine identity of parents and offspring in MaybeRel pairs? 

### Load libraries  --------------------------------------------------------------------
library(sequoia)

### Define custom functions  --------------------------------------------------------------------

### Define key variables --------------------------------------------------------------------

date = "02082024"

### Load data --------------------------------------------------------------------

# "extra" metadata for analyses
load(paste0("../pedigree_exploration/ELF_expanded_metaData.Robj"), verbose = T) # ELF_exp_metaData
load(paste0("../pedigree_exploration/PCC_expanded_metaData.Robj"), verbose = T) # PCC_exp_metaData

# load pedigrees
load(file  = paste0("../pedigree_reconstruction/PCC_pedigree_results_", date, ".Robj"), verbose = TRUE)
# PCC_results <- list(PCC_gt_for_ped_drop2, PCC_meta, PCC_dupe_inds, ped.PCC.ped, error_rate, date) #  PCC.ped.conf,

# load raw metadata
load(file = "../pedigree_reconstruction/PCC_raw_metadata_2013_2020.Robj", verbose = T)
PCC.raw.metaData <- data

### Validate methods --------------------------------------------------------------------

### PCC --------------------------------------------------------------------
PCC.gt = PCC_results[[1]]
PCC.meta = PCC_results[[2]]
PCC.ped.output = PCC_results[[4]]
error_rate = PCC_results[[5]]

if(PCC_results[[6]] != date){
  print("date in this script and date from the pedigree do not match!")
}

# identify maybe relatives
PCC_maybe_rel_ped <- GetMaybeRel(GenoM = PCC.gt, LifeHistData=PCC.meta, Pedigree = PCC.ped.output$Pedigree, Module = "ped", Err = error_rate)

# When were maybe relatives captured and how old were they?


ind.1 <- PCC_maybe_rel_ped$MaybeRel[4,1]
ind.2 <- PCC_maybe_rel_ped$MaybeRel[4,2]

plot_TimeAge_MaybeRels <- function(ind.1, ind.2){
  plot(NULL, ylim = c(0, 70), xlim = c(2008, 2020), xlab = "year", ylab = "SVL")
  for(ind in c(ind.1, ind.2)){
    if(ind == ind.1){
      color = "blue" 
    }else{
      color = "red"
    }
    if(nrow(PCC_exp_metaData[[ind]]) > 0){
      for(j in 1:nrow(PCC_exp_metaData[[ind]])){
        if(PCC_exp_metaData[[ind]][j,"age.class"] == "A"){
          text(x = as.numeric(PCC_exp_metaData[[ind]][j,"year"]), y = as.numeric(PCC_exp_metaData[[ind]][j,"SVL"]), labels = "A", col = color, cex = 0.75, pch = 19)
        }
        if(PCC_exp_metaData[[ind]][j,"age.class"] == "J"){
          text(x = as.numeric(PCC_exp_metaData[[ind]][j,"year"]), y = as.numeric(PCC_exp_metaData[[ind]][j,"SVL"]), labels = "J", col = color, cex = 0.75, pch = 19)
        }
        if(PCC_exp_metaData[[ind]][j,"age.class"] == "N"){
          text(x = as.numeric(PCC_exp_metaData[[ind]][j,"year"]), y = as.numeric(PCC_exp_metaData[[ind]][j,"SVL"]), labels = "N", col = color, cex = 0.75, pch = 19)
        }
      }
    }
    }
  legend("bottomright", legend = c(ind.1, ind.2), col = c("blue", "red"), pch = 19)
  mtext(paste0("Capture time and age of ", ind.1, " and ", ind.2), line = 1)
}

for(n in 1:nrow(PCC_maybe_rel_ped$MaybeRel)){
  plot_TimeAge_MaybeRels(PCC_maybe_rel_ped$MaybeRel[n,1], PCC_maybe_rel_ped$MaybeRel[n,2])
}

# fit a growth curve for each individual? 



###
PCC.raw.metaData[which(PCC.raw.metaData$ID == ind.1),]
PCC.raw.metaData[which(PCC.raw.metaData$ID == ind.2),]




