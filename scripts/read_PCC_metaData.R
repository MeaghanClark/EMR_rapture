## File name: read_PCC_metaData.R
## Purpose: Pull out metadata needed for pedigree reconstruction from field recapture data
## M. I. Clark, November 2022
## Last updated: 2/10/2023

## note: most NA warnings are okay to igore

## revising 4/12/2024 because of observed weirdness in some "known" birth years 

## Loading ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
date = "05262023"
# re-run on 05252023 to account for new individuals dropped in filter_vcf.R

# load raw spreadsheets containing metadata into R
# main metadata file
data_1318 <- read.csv("~/Desktop/EMR_rapture/pedigree_reconstruction/PCC_LifeHistData_2013-18.csv")

# 2019 data
data_19 <- read.csv("~/Desktop/EMR_rapture/pedigree_reconstruction/PCC_LifeHistData_2019.csv")

# 2020 data
data_20 <- read.csv("~/Desktop/EMR_rapture/pedigree_reconstruction/PCC_LifeHistData_2020.csv")


# load list of sequenced individuals 
vcf_date = 02082024
load(paste0("../vcf_filtering/PCC_filtered_gt_all_snps_noZ_", vcf_date, ".Robj"), verbose = T) # PCC_filt_gt_all
seq_inds <- colnames(PCC_filt_gt_all)

# remove duplicated individual ids 
uni.inds <- unlist(lapply(strsplit(seq_inds[which(grepl("_P", seq_inds))], split = "_"), FUN = function(x){paste0(x[1], "_", x[2])}))
seq_inds[which(grepl("_P", seq_inds))] <- uni.inds
seq_inds <- unique(seq_inds)

# ------------------------------------------------------------------------------------------------------------

## Wrangle metadata ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

# match 2013-2018 metadata columnts and 2019 metadata columns
data_1318 <- subset(data_1318, select = -c(X, X.1, X.2))

# What are columns that are in data_1318 but not data_19
colnames(data_1318)[which(is.na(match(colnames(data_1318), colnames(data_19))))] # columns in 2013-18 data that are not in 2019 data (at least by the same name)
# [1] "ID"                            "GPS.Waypoint"                  "EPE..ft."                     
# [4] "Release.Time"                  "Other.Notes.Comments"          "Natural.Community.Description"

# delete release.time from data_1318
# delete "EPE..ft." from data_1318
# delete "Natural.Community.Description" from data_1318

data_1318 <- subset(data_1318, select = -c(Release.Time, EPE..ft., Natural.Community.Description))

# add "Other.Notes.Comments" to data_19
# add column for ID to data_19
# rename GPS.Waypoint.ID" in data_19 to "GPS.Waypoint" 
colnames(data_19)[25] <- "GPS.Waypoint"
data_19$Other.Notes.Comments <- NA
data_19$ID <- NA

dim(data_19)
dim(data_1318)

# What are columns that are in data_19 but not data_1318

colnames(data_19)[which(is.na(match(colnames(data_19), colnames(data_1318))))] # columns in 2013-18 data that are not in 2019 data (at least by the same name)

# [1] "Site.Name"                                                                             
# [2] "Cautery.Brand..location"                                                               
# [3] "Recapture.in.current.year...N..initial.cap.this.year..Y.captured.previously.this.year."
# [4] "Date.s..Previously.captured"                                                           
# [5] "If.female..gravid."                                                                    
# [6] "If.gravid.how.many."                                                                   
# [7] "Subcaudal.Count"                                                                       
# [8] "Age..if.known."                                                                        
# [9] "Rattle.Description"                                                                    
# [10] "Saddle.Description"                                                                    
# [11] "GPS.error..ft."  

# remove "Site.Name" from data_19
# remove "Cautery.Brand..location" from data_19
# remove "Recapture.in.current.year...N..initial.cap.this.year..Y.captured.previously.this.year." from data_19
# remove "Date.s..Previously.captured"from data_19
# remove "If.female..gravid." from data_19
# remove  "If.gravid.how.many." from data_19
# remove "Subcaudal.Count" from data_19
# remove "Age..if.known." (all NAs) from data_19
# remove "Rattle.Description" from data_19
# remove "Saddle.Description" from data_19
# remove "GPS.error..ft." from data_19

data_19 <- subset(data_19, select = -c(Site.Name, Cautery.Brand..location, Recapture.in.current.year...N..initial.cap.this.year..Y.captured.previously.this.year.,
                              Date.s..Previously.captured,  If.female..gravid., 
                              If.gravid.how.many., Subcaudal.Count, Age..if.known., 
                              Rattle.Description, Saddle.Description, GPS.error..ft.))

dim(data_19)
dim(data_1318)

# are the columns the same now? 
colnames(data_19)[which(is.na(match(colnames(data_19), colnames(data_1318))))] 
    # yes! 


# now add in 2020 data... 
colnames(data_19)[which(is.na(match(colnames(data_19), colnames(data_20))))] # columns in 19 data that are not in 20 data (at least by the same name)
# [1] "GPS.Waypoint" "ID"     
colnames(data_20)[21] <- "GPS.Waypoint"
data_20$ID <- NA

colnames(data_20)[which(is.na(match(colnames(data_20), colnames(data_19))))] # columns in 20 data that are not in 19 data
# [1] "Site.Name"                                  "If.female..gravid."                        
# [3] "If.gravid.how.many."                        "Subcaudal.Count"                           
# [5] "Age..if.known."                             "Rattle.Description"                        
# [7] "GPS.error..ft."                             "Release.Time"                              
# [9] "SFD.Notes"                                  "GPS.Coordinates.from.Original.Excel.File"  
# [11] "GPS.Coordinates.from.Original.Excel.File.1"
data_20 <- subset(data_20, select = -c(Site.Name, If.female..gravid., If.gravid.how.many., Subcaudal.Count,
                                       Age..if.known., Rattle.Description, GPS.error..ft., Release.Time,
                                       SFD.Notes, GPS.Coordinates.from.Original.Excel.File, GPS.Coordinates.from.Original.Excel.File.1))


### ID issues

# get DNA ids from lab notebooks for 2019 individuals 

# if any 2019 individuals do not have DNA ids, see if they have an ID from data_1318
dna_id_19 <- c("PCC_307", "PCC_308", "PCC_309", "PCC_310", "PCC_311", "PCC_312", "PCC_313", "PCC_314", "PCC_315")
pit_id_19 <- c(601268594, 601635563, 601262602, 601270888, 601113086, 601771530, 601768262, 600784009, 027078531)

data_19$ID[match(pit_id_19, data_19$PIT.Tag..)] <- dna_id_19

data_19[,c("PIT.Tag..", "ID")]

data_19_hold <- data_19
for(i in 1:length(data_19$PIT.Tag..)){
  if(is.na(data_19$ID[i])){
    if(length(which(data_1318$PIT.Tag.. == data_19$PIT.Tag..[i])) > 0){
      print(which(data_1318$PIT.Tag.. == data_19$PIT.Tag..[i]))
      data_19$ID[i] <- data_1318$ID[which(data_1318$PIT.Tag.. == data_19$PIT.Tag..[i])]
    } # warning means that the Pit tag was assigned more than one DNA id. Can be ignored and will be dealt with later
  }
}

data_1318$ID[is.na(match(data_1318$ID, seq_inds))]

# combine dataframes 
data <- rbind.data.frame(data_1318, data_19, data_20)

# fix weird "PCC_" dna id
data$ID[which(data$ID == "PCC_")] <- data$ID[which(data$PIT.Tag.. == data$PIT.Tag..[which(data$ID == "PCC_")])[2]]

# do we have metadata for all of the sequenced individuals? 

seq_inds[!seq_inds %in% data$ID] # list of sequenced individuals that are not in metadata 

# PCC_297, 301 and 306 removed in vcf_filtering.R 
# "PCC_295" 836329312, same as PCC_202, no changes

# "PCC_296" 601267361
data$ID[which(data$PIT.Tag.. == 601267361)] <- "PCC_296"

# "PCC_298" 027320018
data$ID[which(data$PIT.Tag.. == 027320018)] <- "PCC_298"

# "PCC_299" 601637895
data$ID[which(data$PIT.Tag.. == 601637895)] <- "PCC_299"

# "PCC_300" 601627787
data$ID[which(data$PIT.Tag.. == 601627787)] <- "PCC_300"

# "PCC_302" 836870072
data$ID[which(data$PIT.Tag.. == 836870072)] <- "PCC_302"

# "PCC_303" 601630073
data$ID[which(data$PIT.Tag.. == 601630073)] <- "PCC_303"

# "PCC_304" 27257271
data$ID[which(data$PIT.Tag.. == 27257271)] <- "PCC_304"

# "PCC_305" 601630517, no metadata with that pit tag
data$ID[which(data$PIT.Tag.. == 601630517)] <- "PCC_305"

# if the IDs above are duplicates and metadata is assigned to alt dna id, do nothing
# 
which(data$PIT.Tag.. == 836329312)
data[which(data$PIT.Tag.. == 836329312),]

# data is a dataframe of combine recapture events from 2013 to 2020. I have recorded DNA ids when appropriate. 
# Some DNA ids are duplicates

## Save final metadata data frame ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

save(seq_inds, data, file = "../pedigree_reconstruction/PCC_raw_metadata_2013_2020.Robj")


## Extract information for pedigree reconstruction  ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

metaData <- data.frame(matrix(data = -9, nrow = length(seq_inds), ncol = 5))
colnames(metaData) <- c("ID", "Sex", "BirthYear", "BY.min", "BY.max")
for(i in 1:length(seq_inds)){
  # i = 31
  seq_ind <- seq_inds[i]
  pit <- data$PIT.Tag..[which(data$ID == seq_ind)]
  
  if(sum(duplicated(pit)) != length(pit) -1){
    print(paste0("Multiple or missing PIT tags for individual ", seq_ind, " with index ", i))
    # get info based on id
    snake_data <- data[which(data$ID == seq_ind),]
  }else if(sum(pit == rep("", length = length(pit))) == length(pit)){
    print(paste0("Multiple or missing PIT tags for individual ", seq_ind, " with index ", i))
    # get info based on id
    snake_data <- data[which(data$ID == seq_ind),]
  }else if(length(pit) > 0){
    pit <- pit[1]
    snake_data <- data[which(data$PIT.Tag.. == pit),]
  }else if(length(pit) == 0){ 
    print(paste0("No metadata for ", seq_ind, " with index ", i))
    break
  }
  
    no.caps <- dim(snake_data)[1]
    
    # sex 
    if(sum(snake_data$Sex..M.F.U. == "M") == nrow(snake_data)){
      sex <- 2
    }else if(sum(snake_data$Sex..M.F.U. == "F") == nrow(snake_data)){
      sex <- 1 
    }else if(any(grepl("M", snake_data$Sex..M.F.U.) & !grepl("F", snake_data$Sex..M.F.U.))){
      sex <- 2
    }else if(any(grepl("F", snake_data$Sex..M.F.U.) & !grepl("M", snake_data$Sex..M.F.U.))){
      sex <- 1
    }else{
      sex <- NA 
      print(paste0("There is an error with sex information for individual row number ", i))
    }
    
    # age at first capture
    ages <- snake_data$Age.Class...Adult..Sub.adult.Juvenile..Yearling..Neonate
    if(all(ages == "")){
      print(paste0("all age information missing for individual row number ", i))
    }else if(any(ages == "N")){
      min_age = 0
      max_age = NA
    }else if(any(ages == "Y")){
      min_age = 0
      max_age = 1
    }else if(any(ages == "J")){
      min_age = 1
      max_age = 3
    }else{
      min_age = 3
      max_age = NA
    }
    
    # year of first capture
    if(any(snake_data$Capture.Type..Initial.Capture..Recapture..Recent.Recapture. == "Initial Capture") ){
      cap_year <- snake_data[which(snake_data$Capture.Type..Initial.Capture..Recapture..Recent.Recapture. == "Initial Capture"),"year"]
    }else if(any(snake_data$Capture.Type..Initial.Capture..Recapture..Recent.Recapture. == "Initial capture")){
      cap_year <- snake_data[which(snake_data$Capture.Type..Initial.Capture..Recapture..Recent.Recapture. == "Initial capture"),"year"]
      }else if(any(snake_data$Capture.Type..Initial.Capture..Recapture..Recent.Recapture. == "Recapture")){
      cap_year <- min(snake_data$year)
      }else{
      cap_year <- -9
      print(paste0("issue with capture year in individual row ", i))
    }

    if(min_age == 0){
      birthYear <- as.numeric(cap_year)
      by.max <- -9
      by.min <- -9 
    }else(birthYear <- -9)
    
    if(min_age > 0){
      if(!is.na(max_age)){
        by.max <- as.numeric(cap_year) - min_age
        by.min <- as.numeric(cap_year) - max_age
      }else if(is.na(max_age)){
        by.max <- as.numeric(cap_year) - min_age
        by.min <- -9
      }
    }else(by.max <- NA)
    
    metaData[i,] <- c(seq_inds[i], sex, birthYear, by.min, by.max)
} # end loop

# done! 


# ------------------------------------------------------------------------------------------------------------

# "Errors" thrown by loop
# [1] "There is an error with sex information for individual row number 24"
# [1] "all age information missing for individual row number 24"
# ---> no metadata for this individual at all, just a date of capture and ID

# [1] "There is an error with sex information for individual row number 31"
# --> only captured as neonate, no sex info 

# [1] "There is an error with sex information for individual row number 142"
# [1] "all age information missing for individual row number 142"
# ---> no metadata for this individual at all, just a date of capture and ID

# [1] "There is an error with sex information for individual row number 146"
# [1] "all age information missing for individual row number 146"
# ---> no metadata for this individual at all, just a date of capture and ID


# [1] "There is an error with sex information for individual row number 150"
# [1] "all age information missing for individual row number 150"
# ---> no metadata for this individual at all, just a date of capture and ID

# [1] "Multiple or missing PIT tags for individual PCC_249 with index 163"
# ---> no pit tag id for this individual because it was too small

# [1] "There is an error with sex information for individual row number 177"
# [1] "all age information missing for individual row number 177"
# ---> no metadata for this individual at all, just a date of capture and ID

# [1] "There is an error with sex information for individual row number 180"
# [1] "all age information missing for individual row number 180"
# ---> no metadata for this individual at all, just a date of capture and ID

# [1] "There is an error with sex information for individual row number 181"
# [1] "all age information missing for individual row number 181"
# ---> no metadata for this individual at all, just a date of capture and ID

# [1] "There is an error with sex information for individual row number 183"
# [1] "all age information missing for individual row number 183"
# ---> no metadata for this individual at all, just a date of capture and ID

# [1] "Multiple or missing PIT tags for individual PCC_295 with index 213"
# [1] "all age information missing for individual row number 213"
# ---> no metadata entries for this individual id

# [1] "There is an error with sex information for individual row number 264"
# [1] "all age information missing for individual row number 264"
# ---> no metadata for this individual at all, just a date of capture and ID

# [1] "Multiple or missing PIT tags for individual PCC_97 with index 286"
# [1] "There is an error with sex information for individual row number 286"
# [1] "all age information missing for individual row number 286"
# ---> no metadata for this individual at all, just a date of capture and ID



# Format for sequoia ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
#metaData$BY.min <- rep(NA, nrow(metaData))
#metaData <- metaData[,c(1,2,3,5,4)]

# add duplicates back in 

seq_inds_all <- colnames(PCC_filt_gt_all)
dupes <- seq_inds_all[which(grepl("_P", seq_inds_all))] 
# PCC_300
test <- rbind.data.frame(metaData, metaData[rep(which(metaData$ID == "PCC_300"), 2),])
test2 <- test[-which(test$ID == "PCC_300")[1],]
test2[which(test2$ID == "PCC_300")[1],]$ID <- "PCC_300_P3"
test2[which(test2$ID == "PCC_300")[1],]$ID <- "PCC_300_P8"
metaData <- test2

# PCC_302
test <- rbind.data.frame(metaData, metaData[rep(which(metaData$ID == "PCC_302"), 2),])
test2 <- test[-which(test$ID == "PCC_302")[1],]
test2[which(test2$ID == "PCC_302")[1],]$ID <- "PCC_302_P3"
test2[which(test2$ID == "PCC_302")[1],]$ID <- "PCC_302_P9"
metaData <- test2

#PCC_62
test <- rbind.data.frame(metaData, metaData[rep(which(metaData$ID == "PCC_62"), 2),])
test2 <- test[-which(test$ID == "PCC_62")[1],]
test2[which(test2$ID == "PCC_62")[1],]$ID <- "PCC_62_P10"
test2[which(test2$ID == "PCC_62")[1],]$ID <- "PCC_62_P3"
metaData <- test2

# compare to hand-done file
load(paste0("../pedigree_reconstruction/PCC_pedigree_results10182022.Robj"), verbose = T)

PCC_meta <- PCC_results[[2]]

# same number of individuals?
dim(PCC_meta)
dim(metaData)
# 3 individuals removed when re-filtering vcf (now added back), 3 lab duplicates removed

# same information? 
PCC_meta$ID
metaData$ID
for (i in 1:nrow(PCC_meta)){
  print(PCC_meta[i,] == metaData[which(metaData$ID == PCC_meta$ID[i]),])
}

PCC_LifeHistData <- metaData

# BY.min - Earliest year in which individual may have been born, if exact year is unknown. Ignored when BirthYear is non-missing.
# BY.max (optional) - Latest year in which individual may have been born
# ------------------------------------------------------------------------------------------------------------
# Save ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
save(PCC_LifeHistData, file = paste0("../pedigree_reconstruction/PCC_metaData_",date,".Robj"))

# ------------------------------------------------------------------------------------------------------------

# Graveyard ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------
