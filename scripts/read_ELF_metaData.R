## File name: read_ELF_metaData.R
## Purpose: Pull out metadata needed for pedigree reconstruction from field recapture data
## M. I. Clark, October 2022
## Last updated: 2/10/2023

## note: most NA warnings are okay to ignore

## Loading ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
#date = 05252023
date = "05262023"
# re-run on 05252023 to account for new individuals dropped in filter_vcf.R

library(stringr)
data <- read.csv("~/Desktop/EMR_rapture/pedigree_reconstruction/ELF_LifeHistData.csv")[,1:45]

load(paste0("../vcf_filtering/ELF_filtered_gt_all_snps_noZ05252023.Robj"), verbose = T) # ELF_filt_gt_all
seq_inds <- colnames(ELF_filt_gt_all)
seq_inds <- seq_inds[-which(seq_inds == "ELF_544_P10")]
seq_inds <- str_remove(seq_inds, "_P[0-9]")
seq_inds <- as.integer(unlist(strsplit(seq_inds, split = "_"))[c(FALSE, TRUE)])

# ------------------------------------------------------------------------------------------------------------

## Wrangle metadata ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
# remove capture events with non-numeric IDs
#dim(data)
# data[which(is.na(as.integer(data$ELF.ID))),"ELF.ID"]
data <- data[-which(is.na(as.integer(data$ELF.ID))),]
#dim(data)

# isolate rows contain data for sequenced individual 

metaData <- data.frame(matrix(data = -9, nrow = length(seq_inds), ncol = 7))
colnames(metaData) <- c("ID", "Sex", "BirthYear", "BY.max", "BY.min", "Age", "snake_id")
for(i in 1:length(seq_inds)){
     # i = 511
    snake_data <- data[which(as.integer(data$ELF.ID) == seq_inds[i]),]
  
    # filter out recent recaptures
    snake_data <- snake_data[which(!snake_data$Recapture == "RR"),]
    
    no.caps <- dim(snake_data)[1]
  # snake id the same across capture events? 
  if(!sum(duplicated(snake_data$Snake.ID)) == no.caps-1){
    print(paste0("problem with snake IDs for individual row number ", i))
  }
  
  # sex
  if(sum(is.na(as.integer(snake_data$M))) == 0 & sum(is.na(as.integer(snake_data$F))) == 0){
    if(sum(as.integer(snake_data$M == 0 & as.integer(snake_data$F) == 1)) == no.caps){
      sex <- 1
    }else if(sum(as.integer(snake_data$M == 1 & as.integer(snake_data$F) == 0)) == no.caps) {
      sex <- 2
    }else(print(paste0("error with sex info for individual row number ", i)))
  }else if(sum(snake_data$M == ".") == no.caps & sum(snake_data$F == ".") == no.caps){
    sex <- -9
  }else if(grepl("1", snake_data$M)){
    sex <- 2
  }else if(grepl("1", snake_data$F)){
    sex <- 1
  }else if(grepl("?", snake_data$F) | grepl("?", snake_data$M)){
    sex = -9
    print(paste0("there is a ? in the sex information for individual row number ", i))
  }else(print(paste0("error with sex information for individual row number ", i, " ; sex ID likely contains weird characters")))
  
  # age at first capture
  min_ages <- as.integer(unlist(strsplit(snake_data$Age.at.1st.capture, split = "+")))
  if(all(is.na(min_ages))){
    print(paste0("all age information missing for individual row number ", i))
    if(snake_data$Age.Class == "A"){
      min_age <- 3
    }else(print(paste0("no age or age class information from individual row number ", i)))
  }else(min_age <- min(min_ages, na.rm=TRUE))
  
  # year of first capture
  cap_year <- snake_data[which(snake_data$Recapture == "N" | snake_data$Recapture == "N "),"Year"]
  if(sum(snake_data$Recapture == rep("Y", no.caps)) == no.caps){
    cap_year <- min(snake_data$Year)
  }
  if(length(cap_year) == 0){
    print(paste0("issue with capture year in individual row ", i))
    cap_year <- -9
  }
  if(length(cap_year > 1)){
    cap_year <- min(cap_year)
  }
  
  if(min_age == 0){
    birthYear <- cap_year
    by.min = -9
    by.max = -9
  }else(birthYear <- -9)
  
  if(min_age > 0){
    if(min_age == 1){
      birthYear <- cap_year - min_age
      by.min = -9
      by.max = -9
    }else if (min_age == 2){
      birthYear <- -9
      by.min <- cap_year - 3
      by.max <- cap_year - 1
    }else if (min_age == 3){
      birthYear <- -9
      by.max <- cap_year - min_age
      by.min <- -9
    }
  }else(by.max <- -9)
  
  age_class <- snake_data$Age.at.1st.capture[1]
  
  metaData[i,] <- c(seq_inds[i], sex, birthYear, by.max, by.min, age_class, snake_data$Snake.ID[1])
} # end loop

# done! 

# ELF_65 (i = 511) throws a warning message for snake ID because one record doesn't have "-" 
# separating the PIT tag numbers in the Snake.ID column, can be ignored. 

save(metaData, file = paste0("../pedigree_reconstruction/ELF_rawmetaData_,",date,".Robj"))


# ------------------------------------------------------------------------------------------------------------

# Format for sequoia ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
#metaData$BY.min <- rep(NA, length(metaData))
metaData <- metaData[,-6]
metaData$ID <- paste0("ELF_", metaData$ID)

test <- rbind.data.frame(metaData, metaData[rep(which(metaData$ID == "ELF_544"), 2),])
test2 <- test[-which(test$ID == "ELF_544")[1],]
test2[which(test2$ID == "ELF_544")[1],]$ID <- "ELF_544_P5"
test2[which(test2$ID == "ELF_544")[1],]$ID <- "ELF_544_P10"

ELF_LifeHistData <- test2

# ------------------------------------------------------------------------------------------------------------
# Save ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
save(ELF_LifeHistData, file = paste0("../pedigree_reconstruction/ELF_metaData_,",date,".Robj"))

# ------------------------------------------------------------------------------------------------------------

# Graveyard ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

# issues with ELF 544 and 554
data[which(as.integer(data$ELF.ID) == "544"),]


# fixed issues: 
prob_inds <- c(453, 497, 505, 520, 581, 582, 664)
for (i in prob_inds){
  snake_data <- data[which(as.integer(data$ELF.ID) == seq_inds[i]),]
  print(i)
  print(snake_data$M)
  print(snake_data$F)
}
# ------------------------------------------------------------------------------------------------------------
