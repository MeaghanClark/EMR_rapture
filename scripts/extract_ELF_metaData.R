## File name: extract_ELF_metaData.R
## Purpose: Pull out metadata for pedigree exploration from field recapture data
## M. I. Clark, November 2022
## Last updated: 05/25/2023

## note: most NA warnings are okay to ignore

## Loading ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
data <- read.csv("~/Desktop/EMR_rapture/pedigree_reconstruction/ELF_LifeHistData.csv")[,1:45]

load(paste0("../vcf_filtering/ELF_filtered_gt_all_snps_noZ_05262023.Robj"), verbose = T) # ELF_filt_gt_all
seq_inds <- colnames(ELF_filt_gt_all)
seq_inds <- seq_inds[-which(seq_inds == "ELF_544_P10")]
seq_inds <- str_remove(seq_inds, "_P[0-9]")
seq_inds <- as.integer(unlist(strsplit(seq_inds, split = "_"))[c(FALSE, TRUE)])

# ------------------------------------------------------------------------------------------------------------

## Wrangle metadata ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
# remove capture events with non-numeric IDs
data <- data[-which(is.na(as.integer(data$ELF.ID))),]

# isolate rows contain data for sequenced individual 

metaData <- vector(mode = "list", length = length(seq_inds))
for(i in 1:length(seq_inds)){
  # i = 363
  snake_data <- data[which(as.integer(data$ELF.ID) == seq_inds[i]),]
  
  # info to extract: year, number of recaptures (exclude recent recaps), age, age class, cords, doy, SVL, mass 
  # filter out recent recaptures
  # snake_data <- snake_data[which(!snake_data$Recapture == "RR"),]
  
  no.caps <- dim(snake_data)[1]
  years <- snake_data$Year
  ages <- snake_data$Age
  age.class <- snake_data$Age.Class
  easting <- snake_data$UTM.Easting
  northing <- snake_data$UTM.Northing
  captive <- snake_data$Held.For.Captive.Parturition.
  doy <- snake_data$DOY
  SVL <- snake_data$SVL..cm.
  mass <- snake_data$Mass..g.
  
  # note what capture events were radio telemetry and snake was still underground 
  emerged_status <- rep("emerged", no.caps)
  if(seq_inds[i] == 518 | seq_inds[i] == 533 | seq_inds[i] == 537){ # radio tracked individuals 
    emerged_status[which(grepl("underground", snake_data$Other.Notes))] <- "underground"
  }

  # note when the snake is gravid
  gravid_status <- rep(0, no.caps)
  gravid_status[which(snake_data$Multistate.Gravid == 1)] <- 1
  
  # put info into dataframe
  info <- data.frame(matrix(data = c(years, doy, ages, age.class, easting, northing, captive, SVL, mass, gravid_status, emerged_status), nrow = no.caps, ncol = 11))
  colnames(info) <- c("year", "doy", "age", "age.class", "UTM.easting", "UTM.northing", "held_captive", "SVL", "mass", "gravid_status", "emerged_status")
  
  # assign to list element
  names(metaData)[[i]] <- paste0("ELF_", seq_inds[i])
  metaData[[i]] <- info
} # end loop

# done! 

# ------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------
# Save ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
ELF_exp_metaData <- metaData
save(ELF_exp_metaData, file = "../pedigree_exploration/ELF_expanded_metaData.Robj")

# ------------------------------------------------------------------------------------------------------------

