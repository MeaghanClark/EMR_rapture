## File name: organize_metadata.R
## Purpose: Load and organize mark recapture metadata for EMR project pedigree reconstruction and spatial analyses
## M. I. Clark, November 2022
## Last updated: 06/26/2024
    # incorporating Danielle's birth year estimates
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
#
# -----------------------------------------------------------------------------------------------------------------------------------------------------------

## Script set up --------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
library(stringr)

# load Danielle's by estimates 
load("../pedigree_reconstruction/DRB_byEstimates.Robj", verbose = T) # export_temp

# make data frame of asymptotic size values
asymptotic_L <- data.frame(matrix(NA, nrow = 4, ncol = 5))
colnames(asymptotic_L) <- c("site", "sex", "lower_CI", "upper_CI", "L_est")
asymptotic_L[1,] <- c("cass", "female", 57.5, 59.8, 58.7) # 58.7 cm (95% CI = 57.5 – 59.8) 
asymptotic_L[2,] <- c("cass", "male", 60.3, 63.3, 61.8) # 61.8 cm (95% CI = 60.3 – 63.3)
asymptotic_L[3,] <- c("barry", "female", 58.9, 61.6, 60.3) # 60.3 cm (95% CI = 58.9 – 61.6)
asymptotic_L[4,] <- c("barry", "male", 61.8, 65.1, 63.4) # 63.4 cm (95% CI = 61.8 – 65.1) 

asymptotic_L$lower_CI <- as.numeric(asymptotic_L$lower_CI)
asymptotic_L$upper_CI <- as.numeric(asymptotic_L$upper_CI)
asymptotic_L$L_est <- as.numeric(asymptotic_L$L_est)


infer_birth_year_wSVL <- function(ind, location, ind_SVL, ind_year, ind_sex){
  # ind = PCC_known_by_inds[i]
  # location = "barry"
  # ind_SVL = as.numeric(PCC_exp_metaData[[PCC_known_by_inds[i]]][n,"SVL"])
  # ind_year =PCC_exp_metaData[[PCC_known_by_inds[i]]][n,"year"] 
  # ind_sex = 1
  
  ind_year <- as.numeric(ind_year)
  if(ind_sex == 1){
    ind_sex = "female"
  }else if(ind_sex == 2){
    ind_sex = "male"
  }
  ind_age <- est_age_from_curve(L = ind_SVL, Loo = subset(asymptotic_L, site == location & sex == ind_sex, select = L_est))
  ind_age_upper <- est_age_from_curve(L = ind_SVL, Loo = subset(asymptotic_L, site == location & sex == ind_sex, select = upper_CI))
  ind_age_lower <- est_age_from_curve(L = ind_SVL, Loo = subset(asymptotic_L, site == location & sex == ind_sex, select = lower_CI))
  
  if(is.null(ind_age)){
    ind_age = NA # some very large number
  }
  if(is.null(ind_age_upper)){
    ind_age_upper = NA  # some very large number
  }
  if(is.null(ind_age_lower)){
    ind_age_lower = NA # some very large number
  }
  
  return(c(as.numeric(ind_year - ind_age), as.numeric(ind_year - ind_age_upper), as.numeric(ind_year - ind_age_lower)))
}

# -----------------------------------------------------------------------------------------------------------------------------------------------------------

## Barry County: Load metadata  -----------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# load raw spreadsheets containing metadata into R
# main metadata file
data_1318 <- read.csv("~/Desktop/EMR_rapture/pedigree_reconstruction/PCC_LifeHistData_2013-18.csv")

# 2019 data
data_19 <- read.csv("~/Desktop/EMR_rapture/pedigree_reconstruction/PCC_LifeHistData_2019.csv")

# 2020 data
data_20 <- read.csv("~/Desktop/EMR_rapture/pedigree_reconstruction/PCC_LifeHistData_2020.csv")
  # can ignore warning message

# blank/missing individuals 
data_missing <- read.csv("~/Desktop/EMR_rapture/pedigree_reconstruction/Data for Meaghan.csv")


# load list of sequenced individuals 
vcf_date = "02082024"
load(paste0("../vcf_filtering/PCC_filtered_gt_all_snps_noZ_", vcf_date, ".Robj"), verbose = T) # PCC_filt_gt_all
seq_inds <- colnames(PCC_filt_gt_all)

# remove technical duplicate individual ids 
uni.inds <- unlist(lapply(strsplit(seq_inds[which(grepl("_P", seq_inds))], split = "_"), FUN = function(x){paste0(x[1], "_", x[2])}))
seq_inds[which(grepl("_P", seq_inds))] <- uni.inds
seq_inds <- unique(seq_inds)

# -----------------------------------------------------------------------------------------------------------------------------------------------------------

## Barry County: Wrangle metadata ---------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# match 2013-2018 metadata columns and 2019 metadata columns so we can combine

data_1318 <- subset(data_1318, select = -c(X, X.1, X.2)) # remove empty columns

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

# fix weird "PCC_" dna id using PIT tag to get actual DNA id! 
data$ID[which(data$ID == "PCC_")] <- data$ID[which(data$PIT.Tag.. == data$PIT.Tag..[which(data$ID == "PCC_")])[2]] 

# do we have metadata for all of the sequenced individuals? 

seq_inds[!seq_inds %in% data$ID] # list of sequenced individuals that are not in metadata 
# "PCC_295"x "PCC_296"x "PCC_298"x "PCC_299"x "PCC_300"x "PCC_302"x "PCC_303"x "PCC_304"x
# "PCC_305"x

# Use PIT tags to add missing DNA ids
# "PCC_295" 836329312, same as PCC_202, no changes
which(data$PIT.Tag.. == 836329312)
data[which(data$PIT.Tag.. == 836329312),]

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

### incorporate missing metadata received from Jen 5/2/2024
colnames(data)
colnames(data_missing)

colnames(data)[which(is.na(match(colnames(data), colnames(data_missing))))] # columns in data that are not in data_missing
# [1] "ID"                                                         
# [2] "Capture.Type..Initial.Capture..Recapture..Recent.Recapture."
# [3] "GPS.Waypoint"                                               
# [4] "Decimal.Degrees.Northing"                                   
# [5] "Decimal.Degrees.Easting"                                    
# [6] "Survey.Team"   
data_missing$ID <- NA
data_missing$Capture.Type..Initial.Capture..Recapture..Recent.Recapture. <- NA
data_missing$GPS.Waypoint <- data_missing$GPS.Waypoint.ID
data_missing <- subset(data_missing, select = -GPS.Waypoint.ID)
data_missing$Decimal.Degrees.Northing <- data_missing$Lat
data_missing <- subset(data_missing, select = -Lat)
data_missing$Decimal.Degrees.Easting <- data_missing$Long
data_missing <- subset(data_missing, select = -Long)
data_missing$Survey.Team <- NA

colnames(data_missing)[which(is.na(match(colnames(data_missing), colnames(data))))] # columns in data_missing that are not in data
# [1] "Cautery.Brand..location" "male"                   
# [3] "female"                  "If.female..gravid."     
# [5] "If.gravid.how.many."     "Recent.meal."           
# [7] "Subcaudal.Count"         "Rattle.Description"  
data_missing <- subset(data_missing, select = -c(Cautery.Brand..location, male, female, If.female..gravid., If.gravid.how.many., Recent.meal., Subcaudal.Count, Rattle.Description))

# column names now match
# combine dataframes
data <- rbind.data.frame(data, data_missing)

# add DNA IDs for data_missing individuals (now NA) and remove empty entries

# PCC_121
data[which(data$PIT.Tag..==126),] 
# remove first entry
data <- data[-200,]
# add DNA ID to new, complete entry
data[which(data$PIT.Tag..==126),"ID"] <- "PCC_121"
data[which(data$PIT.Tag..==126),"Capture.Type..Initial.Capture..Recapture..Recent.Recapture."] <- "Initial Capture"


# PCC_229
data[which(data$PIT.Tag..==027265582),] 
# remove first entry
data <- data[-64,]
# add DNA ID to new, complete entry
data[which(data$PIT.Tag..==027265582),"ID"] <- "PCC_229"
data[which(data$PIT.Tag..==027265582),"Capture.Type..Initial.Capture..Recapture..Recent.Recapture."] <- "Initial Capture"


# PCC_232
data[which(data$PIT.Tag..==027082020),] 
# remove first entry
data <- data[-65,]
# add DNA ID to new, complete entry
data[which(data$PIT.Tag..==027082020),"ID"] <- "PCC_232"
data[which(data$PIT.Tag..==027082020),"Capture.Type..Initial.Capture..Recapture..Recent.Recapture."] <- "Initial Capture"


# PCC_237
data[which(data$PIT.Tag..==027085319),] 
# remove first entry
data <- data[-which(data$ID == "PCC_237"),]
# add DNA ID to new, complete entry
data[which(data$PIT.Tag..==027085319),"ID"] <- "PCC_237"
data[which(data$PIT.Tag..==027085319),"Capture.Type..Initial.Capture..Recapture..Recent.Recapture."] <- "Initial Capture"


# PCC_261
data[which(data$PIT.Tag..==842113124),] 
# remove first entry
data <- data[-which(data$ID == "PCC_261"),]
# add DNA ID to new, complete entry
data[which(data$PIT.Tag..==842113124),"ID"] <- "PCC_261"
data[which(data$PIT.Tag..==842113124),"Capture.Type..Initial.Capture..Recapture..Recent.Recapture."] <- "Initial Capture"


# PCC_267
data[which(data$PIT.Tag..==601634258),] 
# remove first entry
data <- data[-which(data$ID == "PCC_267"),]
# add DNA ID to new, complete entry
data[which(data$PIT.Tag..==601634258),"ID"] <- "PCC_267"
data[which(data$PIT.Tag..==601634258),"Capture.Type..Initial.Capture..Recapture..Recent.Recapture."] <- "Initial Capture"


# PCC_77, real PIT 011865563, 11863563 was a typo
data[which(data$PIT.Tag..==11863563),] 
data[which(data$PIT.Tag..==011865563),] 
# remove first entry
data <- data[-which(data$ID == "PCC_77"),]
# add DNA ID to new, complete entry
data[which(data$PIT.Tag..==011865563),"ID"] <- "PCC_77"
data[which(data$PIT.Tag..==011865563),"Capture.Type..Initial.Capture..Recapture..Recent.Recapture."] <- "Initial Capture"

# PCC_97, original PIT tag, 27287796, was lost (very rare), new PIT tag is 027268822
data[which(data$PIT.Tag..==27287796),] 
data[which(data$PIT.Tag..==027268822),] 
# remove first entry
data <- data[-which(data$ID == "PCC_97"),]
# add DNA ID to new, complete entry
data[which(data$PIT.Tag..==027268822),"ID"] <- "PCC_97"
data[which(data$PIT.Tag..==027268822),"Capture.Type..Initial.Capture..Recapture..Recent.Recapture."] <- "Recapture"

# still missing metadata for PCC_264 and PCC_265

# data is a dataframe of combine recapture events from 2013 to 2020. I have recorded DNA ids when appropriate. 
# Some DNA ids are duplicates

# update data with incorrect SVL for individual PCC_303 PIT 601630073
# 18.4 inches is 46.736 cm 
data[which(data$ID == "PCC_303"), "SVL..cm."] <- as.character(as.numeric(data[which(data$ID == "PCC_303"), "SVL..cm."]) * 2.54)

PCC_raw_meta <- list(seq_inds, data)

save(PCC_raw_meta, file = "../pedigree_reconstruction/PCC_raw_metadata.Robj")

# -----------------------------------------------------------------------------------------------------------------------------------------------------------

## Barry County: make metadata for pedigree reconstruction  -------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# want a dataframe with columns: "ID", "Sex", "BirthYear", "BY.min", and "BY.max"
# BY.min - Earliest year in which individual may have been born, if exact year is unknown. Ignored when BirthYear is non-missing.
# BY.max (optional) - Latest year in which individual may have been born

# ORIGINAL metadata based on age classes -----------------------------------------------------------------------------------------------------------------------------
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
    sex <- -9
    print(paste0("There is an error with sex information for individual row number ", i))
  }
  
  # age at first capture
  ages <- snake_data$Age.Class...Adult..Sub.adult.Juvenile..Yearling..Neonate
  if(all(ages == "")){
    print(paste0("all age information missing for individual row number ", i))
  }else if(any(grepl("N", ages))){
    min_age = 0
    max_age = NA
  }else if(any(grepl("Y", ages))){
    min_age = 0
    max_age = 1
  }else if(any(grepl("J", ages))){
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
  
  if(min_age == 0 & is.na(max_age)){
    birthYear <- as.numeric(cap_year)
    by.max <- -9
    by.min <- -9 
  }else(birthYear <- -9)
  
  if(cap_year != -9 & birthYear == -9){
    if(!is.na(max_age)){
      by.max <- as.numeric(cap_year) - min_age
      by.min <- as.numeric(cap_year) - max_age
    }else if(is.na(max_age)){ # max age will be na if you are a neonate (dealt with above), or adult/other
      by.max <- as.numeric(cap_year) - min_age
      by.min <- -9
    }
  }else{
    by.max <- -9
    by.min <- -9
  }
  
  metaData[i,] <- c(seq_inds[i], sex, birthYear, by.min, by.max)
} # end loop

# "Errors" thrown by loop -----------------------------------------------------------------------------------------------------------------------------------

data[which(data$ID == seq_inds[31]),]
# [1] "There is an error with sex information for individual row number 31"
# --> only captured as neonate, no sex info 

data[which(data$ID == seq_inds[35]),]
# [1] "issue with capture year in individual row 35"
# recent recapture, but this is the only capture record that we have 
data[which(data$PIT.Tag.. == "027082277"),]
metaData[35,"BY.max"] <- 2015
metaData[35,"BY.min"] <- 2014

data[which(data$ID == seq_inds[163]),]
# [1] "Multiple or missing PIT tags for individual PCC_249 with index 163"
# ---> no pit tag id for this individual because it was too small

data[which(data$ID == seq_inds[180]),]
# [1] "There is an error with sex information for individual row number 180"
# [1] "all age information missing for individual row number 180"
# ---> no metadata for this individual at all, just a date of capture and ID, PCC_264, PIT 842079064

data[which(data$ID == seq_inds[181]),]
# [1] "There is an error with sex information for individual row number 181"
# [1] "all age information missing for individual row number 181"
# ---> no metadata for this individual at all, just a date of capture and ID, PCC_265, PIT 312

data[which(data$ID == seq_inds[212]),]
# [1] "Multiple or missing PIT tags for individual PCC_295 with index 212"
# [1] "all age information missing for individual row number 212"
# ---> no metadata entries for this individual id, PCC_295, duplicated individual

data[which(data$ID == seq_inds[233]),]
# [1] "issue with capture year in individual row 233"
# recent recapture, but this is the only capture record that we have 
data[which(data$PIT.Tag.. == "075529074"),]
metaData[233,"BY.max"] <- 2010
metaData[233,"BY.min"] <- -9


# Format for sequoia ------------------------------------------------------------------------------------------------------------------------------------------

# add duplicates back in so we can verify detection during pedigree reconstruction
seq_inds_all <- colnames(PCC_filt_gt_all)
dupes <- seq_inds_all[which(grepl("_P", seq_inds_all))] 
# PCC_300
test <- rbind.data.frame(metaData, metaData[rep(which(metaData$ID == "PCC_300"), 2),])
test2 <- test[-which(test$ID == "PCC_300")[1],]
test2[which(test2$ID == "PCC_300")[1],]$ID <- "PCC_300_P3"
test2[which(test2$ID == "PCC_300")[1],]$ID <- "PCC_300_P8"
PCC_LifeHistData <- test2

# PCC_302
test <- rbind.data.frame(metaData, metaData[rep(which(metaData$ID == "PCC_302"), 2),])
test2 <- test[-which(test$ID == "PCC_302")[1],]
test2[which(test2$ID == "PCC_302")[1],]$ID <- "PCC_302_P3"
test2[which(test2$ID == "PCC_302")[1],]$ID <- "PCC_302_P9"
PCC_LifeHistData <- test2

#PCC_62
test <- rbind.data.frame(metaData, metaData[rep(which(metaData$ID == "PCC_62"), 2),])
test2 <- test[-which(test$ID == "PCC_62")[1],]
test2[which(test2$ID == "PCC_62")[1],]$ID <- "PCC_62_P10"
test2[which(test2$ID == "PCC_62")[1],]$ID <- "PCC_62_P3"
PCC_LifeHistData <- test2

# change -9s to NAs
PCC_LifeHistData[which(PCC_LifeHistData$BirthYear == -9),"BirthYear"] <- NA
PCC_LifeHistData[which(PCC_LifeHistData$BY.min == -9),"BY.min"] <- NA
PCC_LifeHistData[which(PCC_LifeHistData$BY.max == -9),"BY.max"] <- NA


# Save metadata ----------------------------------------------------------------------------------------------------------------------------------------------
save(PCC_LifeHistData, file = paste0("../pedigree_reconstruction/PCC_metaData_original_",vcf_date,".Robj"))

# ------------------------------------------------------------------------------------------------------------------------------------------------------------

# NEW method: use SVL to update metadata 

# editting metaData
# change -9s to NAs
metaData[which(metaData$BirthYear == -9),"BirthYear"] <- NA
metaData[which(metaData$BY.min == -9),"BY.min"] <- NA
metaData[which(metaData$BY.max == -9),"BY.max"] <- NA
orig_metaData <- metaData

# get stats for original dataframe
sum(is.na(orig_metaData$BY.min) & is.na(orig_metaData$BirthYear)) # 168 inds have no by.min
sum(is.na(orig_metaData$BY.max)& is.na(orig_metaData$BirthYear)) # 3 inds have no by.max


PCC_DRB_by <- merge(x = metaData, y = export_temp, by.x = "ID", by.y = "id", all.x = TRUE)
plot(PCC_DRB_by$estBirthYear~PCC_DRB_by$BirthYear.x, xlim = c(1996, 2023)) # known birth year individuals known in estBirthYears: yes
abline(a = 0, b =1 )

plot(PCC_DRB_by$BirthYear.y~PCC_DRB_by$BirthYear.x, xlim = c(1996, 2023)) # known birth year individuals known in estBirthYears: yes
abline(a = 0, b =1 )

# update metaData data frame with inferred age information from PCC_DRB_by

# orig: PCC_43   1        -9     -9   2010
# new: PCC_43   1        -9   2004   2008
for(i in 1:length(seq_inds)){
  ind <- seq_inds[i]
  # if birth year is known, do nothing
  if(is.na(metaData[which(metaData$ID == ind), "BirthYear"])){  # if birth year is not known
    if(!is.na(PCC_DRB_by[which(PCC_DRB_by$ID == ind), "estBirthYear"])){    # check if we have a birth year estimate,
      # if so, calculate BY.min and BY.max with estBirthYear -2 and estBirthYear + 2
      # get +/- 2 year range around inferred birth year
      i_by_max <- as.numeric(PCC_DRB_by[which(PCC_DRB_by$ID == ind), "estBirthYear"]) + 2
      i_by_min <- as.numeric(PCC_DRB_by[which(PCC_DRB_by$ID == ind), "estBirthYear"]) - 2
    
      if(is.na(metaData[which(metaData$ID == ind), "BY.max"])){ # if by.max is NA, replace with i_by_max
        metaData[metaData$ID == ind,"BY.max"] <- i_by_max
      }else if(i_by_max < subset(metaData, ID == ind, select = BY.max)){ # if inferred by.max is less than original (hypothetically CLOSER to by)
        metaData[metaData$ID == ind,"BY.max"] <- i_by_max
      }
      if(is.na(metaData[which(metaData$ID == ind), "BY.min"])){ # if by.min is NA, replace with i_by_min
        metaData[metaData$ID == ind,"BY.min"] <- i_by_min
      }else if(i_by_min > subset(metaData, ID == ind, select = BY.min)){ # if by.min is not NA, replace with i_by_min if i_by_min is greater than original (hypothetically CLOSER to by)
        metaData[metaData$ID == ind,"BY.min"] <- i_by_min
      }
    }
  } # if not, keep BY.max as informed by last recap 
}


# two longevity outliers... revert back to original metadata due to likely error surrounding by estimates
# 11868056 --> large SVL, 63.4, PCC_12
# 27331044 --> large SVL, 63.3, PCC_134

metaData[which(metaData$ID == "PCC_12"),] <- orig_metaData[which(orig_metaData$ID == "PCC_12"),]
metaData[which(metaData$ID == "PCC_134"),] <- orig_metaData[which(orig_metaData$ID == "PCC_134"),]


# PCCI metadata updated! 
# Meaghan's original method: update data frame with inferred age information
#--------------------------------------------------------------------------------------------------------------------
# for(i in 1:length(seq_inds)){
#   ind <- seq_inds[i]
#   
#   if(nrow(data[which(data$ID == ind),])> 0){
#     # don't want min SVL... want SVL at earliest capture year? 
#     
#     SVL <- subset(data, ID == ind & year == min(subset(data, ID==ind, select = year)$year), select = SVL..cm.)$SVL..cm.
#     
#     # get year at min SVL
#     cap_year <- subset(data, ID == ind & SVL..cm. == SVL, select = year)
#     
#     # get sex information
#     sex <- subset(data, ID == ind & SVL..cm. == SVL, select = Sex..M.F.U.)
#     
#     if(sex == "M"){
#       sex = 2
#     }else if(sex == "F"){
#       sex = 1
#     }else{
#       sex = NA
#     }
#     # If SVL is less than 50 (fairly accurate inference), use growth curve to approx age
#     
#     if(is.na(subset(metaData, ID == ind, select = BirthYear)) & SVL <= 50 & !is.na(sex) & !is.na(as.numeric(SVL))){
#       print(paste0("updating metaData for individual ", ind, " i is ", i))
#       inferred_by <- round(infer_birth_year_wSVL(ind, "barry", ind_SVL = as.numeric(SVL), ind_year = cap_year, ind_sex = sex)[1])
#       # get +/- 2 year range around inferred birth year
#       i_by_max <- inferred_by + 2
#       i_by_min <- inferred_by - 2
#       
#       # update PCC_LifeHistData with inferred information
#       if(is.na(subset(metaData, ID == ind, select = BY.max))){ # if by.max is NA, replace with i_by_max
#         metaData[metaData$ID == ind,"BY.max"] <- i_by_max
#       }else if(i_by_max < subset(metaData, ID == ind, select = BY.max)){ # if inferred by.max is less than original (hypothetically CLOSER to by)
#         metaData[metaData$ID == ind,"BY.max"] <- i_by_max
#       }
#       if(is.na(subset(metaData, ID == ind, select = BY.min))){ # if by.min is NA, replace with i_by_min
#         metaData[metaData$ID == ind,"BY.min"] <- i_by_min
#       }else if(i_by_min > subset(metaData, ID == ind, select = BY.min)){ # if by.min is not NA, replace with i_by_min if i_by_min is greater than original (hypothetically CLOSER to by)
#         metaData[metaData$ID == ind,"BY.min"] <- i_by_min
#       }
#       
#     }
#   }
# }
#--------------------------------------------------------------------------------------------------------------------

# get stats for updated dataframe 
sum(is.na(metaData$BY.min) & is.na(metaData$BirthYear)) # 37 inds have no by.min
sum(is.na(metaData$BY.max)& is.na(metaData$BirthYear)) # 3 inds have no by.max

# Format for sequoia ------------------------------------------------------------------------------------------------------------------------------------------

# add duplicates back in so we can verify detection during pedigree reconstruction
seq_inds_all <- colnames(PCC_filt_gt_all)
dupes <- seq_inds_all[which(grepl("_P", seq_inds_all))] 
# PCC_300
test <- rbind.data.frame(metaData, metaData[rep(which(metaData$ID == "PCC_300"), 2),])
test2 <- test[-which(test$ID == "PCC_300")[1],]
test2[which(test2$ID == "PCC_300")[1],]$ID <- "PCC_300_P3"
test2[which(test2$ID == "PCC_300")[1],]$ID <- "PCC_300_P8"
PCC_lifeHistDataInfer <- test2

# PCC_302
test <- rbind.data.frame(metaData, metaData[rep(which(metaData$ID == "PCC_302"), 2),])
test2 <- test[-which(test$ID == "PCC_302")[1],]
test2[which(test2$ID == "PCC_302")[1],]$ID <- "PCC_302_P3"
test2[which(test2$ID == "PCC_302")[1],]$ID <- "PCC_302_P9"
PCC_lifeHistDataInfer <- test2

#PCC_62
test <- rbind.data.frame(metaData, metaData[rep(which(metaData$ID == "PCC_62"), 2),])
test2 <- test[-which(test$ID == "PCC_62")[1],]
test2[which(test2$ID == "PCC_62")[1],]$ID <- "PCC_62_P10"
test2[which(test2$ID == "PCC_62")[1],]$ID <- "PCC_62_P3"
PCC_lifeHistDataInfer <- test2

# change -9s to NAs
PCC_lifeHistDataInfer[which(PCC_lifeHistDataInfer$BirthYear == -9),"BirthYear"] <- NA
PCC_lifeHistDataInfer[which(PCC_lifeHistDataInfer$BY.min == -9),"BY.min"] <- NA
PCC_lifeHistDataInfer[which(PCC_lifeHistDataInfer$BY.max == -9),"BY.max"] <- NA


# Save metadata ----------------------------------------------------------------------------------------------------------------------------------------------
save(PCC_lifeHistDataInfer, file = paste0("../pedigree_reconstruction/PCC_metaData_infer_DRB",vcf_date,".Robj"))


## Barry County: get expanded metadata ----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# use data and seq_inds from above

exp_metaData <- vector(mode = "list", length = length(seq_inds))
for(i in 1:length(seq_inds)){
  # i = 11
  
  # isolate data for given individual 
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
    print(paste0("No exp_metaData for ", seq_ind, " with index ", i))
    break
  }
  
  # info to extract: year, number of recaptures (exclude recent recaps), age, age class, cords
  
  no.caps <- dim(snake_data)[1]
  years <- snake_data$year
  age.class <- snake_data$Age.Class...Adult..Sub.adult.Juvenile..Yearling..Neonate
  easting <- snake_data$Decimal.Degrees.Easting
  northing <- snake_data$Decimal.Degrees.Northing
  SVL <- snake_data$SVL..cm.
  mass <- snake_data$Mass..g.
  # gravid_status <- snake_data$ # don't have gravidity status for most PCC recaptures! ping Jen about this... seems likely that this data exists somewhere! 
  
  # put info into dataframe
  info <- data.frame(matrix(data = c(years, age.class, easting, northing, SVL, mass), nrow = no.caps, ncol = 6))
  colnames(info) <- c("year", "age.class", "DD.easting", "DD.northing", "SVL", "mass")
  
  # assign to list element
  names(exp_metaData)[[i]] <- paste0(seq_inds[i])
  exp_metaData[[i]] <- info
} # end loop

# done! 

# [1] "Multiple or missing PIT tags for individual PCC_249 with index 163"
# --> individual too small for PIT tag

# [1] "Multiple or missing PIT tags for individual PCC_295 with index 212"
# --> duplicated indvidual, this is okay 

# [1] "Multiple or missing PIT tags for individual PCC_97 with index 283"
# --> PIT tag written in Encounter.date column, but no records are available for this individual

PCC_exp_metaData <- exp_metaData
save(PCC_exp_metaData, file = "../pedigree_exploration/PCC_expanded_metaData.Robj")

# ------------------------------------------------------------------------------------------------------------



# -----------------------------------------------------------------------------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# ***********************************************************************************************************************************************************
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------

## Cass County: Load metadata  --------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
data <- read.csv("~/Desktop/EMR_rapture/pedigree_reconstruction/ELF_LifeHistData.csv")[,1:45]

load(paste0("../vcf_filtering/ELF_filtered_gt_all_snps_noZ_", vcf_date, ".Robj"), verbose = T) # ELF_filt_gt_all
seq_inds <- colnames(ELF_filt_gt_all)
seq_inds <- seq_inds[-which(seq_inds == "ELF_544_P10")]
seq_inds <- str_remove(seq_inds, "_P[0-9]")
seq_inds <- as.integer(unlist(strsplit(seq_inds, split = "_"))[c(FALSE, TRUE)])

## Cass County: Wrangle metadata  --------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# remove capture events with non-numeric IDs
#dim(data)
# data[which(is.na(as.integer(data$ELF.ID))),"ELF.ID"]
data <- data[-which(is.na(as.integer(data$ELF.ID))),]
#dim(data)

# isolate rows contain data for sequenced individual 

metaData <- data.frame(matrix(data = NA, nrow = length(seq_inds), ncol = 7))
colnames(metaData) <- c("ID", "Sex", "BirthYear", "BY.max", "BY.min", "Age", "snake_id")
for(i in 1:length(seq_inds)){
  # i = 11
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
  if(length(cap_year) >= 1){
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

# [1] "all age information missing for individual row number 713"
# no "Age.at.1st.capture" information but age class was A

# [1] "all age information missing for individual row number 716"
# no "Age.at.1st.capture" information but age class was A

save(metaData, file = paste0("../pedigree_reconstruction/ELF_raw_metadata_,",vcf_date,".Robj"))


## Cass County: make metadata for pedigree reconstruction  --------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------

# ORIGINAL: metadate based on age classes -------------------------------------------------------------------------------------------------------------------
metaData <- metaData[,-6]
metaData$ID <- paste0("ELF_", metaData$ID)

test <- rbind.data.frame(metaData, metaData[rep(which(metaData$ID == "ELF_544"), 2),])
test2 <- test[-which(test$ID == "ELF_544")[1],]
test2[which(test2$ID == "ELF_544")[1],]$ID <- "ELF_544_P5"
test2[which(test2$ID == "ELF_544")[1],]$ID <- "ELF_544_P10"

# change -9s to NAs
test2[which(test2$BirthYear == -9),"BirthYear"] <- NA
test2[which(test2$BY.min == -9),"BY.min"] <- NA
test2[which(test2$BY.max == -9),"BY.max"] <- NA

ELF_LifeHistData <- test2

# save metadata----------------------------------------------------------------------------------------------------------------------------------------------
save(ELF_LifeHistData, file = paste0("../pedigree_reconstruction/ELF_metaData_original_",vcf_date,".Robj"))

# -----------------------------------------------------------------------------------------------------------------------------------------------------------

# ALT: Use SVL to estimate BY.min and BY.max ----------------------------------------------------------------------------------------------------------------
# NEW method: use SVL when SVL is less than 50 cm to update metadata 

# editting metaData
# change -9s to NAs
metaData[which(metaData$BirthYear == -9),"BirthYear"] <- NA
metaData[which(metaData$BY.min == -9),"BY.min"] <- NA
metaData[which(metaData$BY.max == -9),"BY.max"] <- NA
orig_metaData <- metaData


# get stats for original dataframe
sum(is.na(orig_metaData$BY.min) & is.na(orig_metaData$BirthYear)) # 217 inds have no by.min
sum(is.na(orig_metaData$BY.max)& is.na(orig_metaData$BirthYear)) # 0 inds have no by.max

#    **** --------------------------------------------------------------------------------------------------------------------
ELF_DRB_by <- merge(x = metaData, y = export_temp, by.x = "ID", by.y = "id", all.x = TRUE)
plot(ELF_DRB_by$estBirthYear~ELF_DRB_by$BirthYear.x, xlim = c(1996, 2023)) # known birth year individuals known in estBirthYears: yes
abline(a = 0, b =1 )

plot(ELF_DRB_by$BirthYear.y~ELF_DRB_by$BirthYear.x, xlim = c(1996, 2023)) # known birth year individuals known in estBirthYears: yes
abline(a = 0, b =1 )

# update metaData data frame with inferred age information from ELF_DRB_by

# orig: ELF_101   1      <NA>   2007   <NA> 047-015-793
# new: PCC_43   1        -9   2004   2008
for(i in 1:nrow(metaData)){
  ind <- metaData[i, "ID"]
  # if birth year is known, do nothing
  if(is.na(metaData[which(metaData$ID == ind), "BirthYear"])){  # if birth year is not known
    if(!is.na(ELF_DRB_by[which(ELF_DRB_by$ID == ind), "estBirthYear"])){    # check if we have a birth year estimate,
      # if so, calculate BY.min and BY.max with estBirthYear -2 and estBirthYear + 2
      # get +/- 2 year range around inferred birth year
      i_by_max <- as.numeric(ELF_DRB_by[which(ELF_DRB_by$ID == ind), "estBirthYear"]) + 2
      i_by_min <- as.numeric(ELF_DRB_by[which(ELF_DRB_by$ID == ind), "estBirthYear"]) - 2
      
      if(is.na(metaData[which(metaData$ID == ind), "BY.max"])){ # if by.max is NA, replace with i_by_max
        metaData[metaData$ID == ind,"BY.max"] <- i_by_max
      }else if(i_by_max < subset(metaData, ID == ind, select = BY.max)){ # if inferred by.max is less than original (hypothetically CLOSER to by)
        metaData[metaData$ID == ind,"BY.max"] <- i_by_max
      }
      if(is.na(metaData[which(metaData$ID == ind), "BY.min"])){ # if by.min is NA, replace with i_by_min
        metaData[metaData$ID == ind,"BY.min"] <- i_by_min
      }else if(i_by_min > subset(metaData, ID == ind, select = BY.min)){ # if by.min is not NA, replace with i_by_min if i_by_min is greater than original (hypothetically CLOSER to by)
        metaData[metaData$ID == ind,"BY.min"] <- i_by_min
      }
    }
  } # if not, keep BY.max as informed by last recap 
}


#    **** --------------------------------------------------------------------------------------------------------------------

# Meaghan's original method: update data frame with inferred age information
#--------------------------------------------------------------------------------------------------------------------
# update data frame with inferred age information
for(i in 1:length(seq_inds)){
  ind <- seq_inds[i]
  
  # limit data to years when SVL was measured
  
  data_sub <- subset(data, ELF.ID == ind & SVL..cm. != ".", select = c(Year, SVL..cm.))
  
  if(dim(data_sub)[1] > 0){
    # don't want min SVL... want SVL at earliest capture year
    SVL <- subset(data_sub, Year == min(Year), select = SVL..cm.)
    
    # get year at min SVL
    cap_year <- subset(data, ELF.ID == ind & SVL..cm. == SVL$SVL..cm.[1], select = Year)$Year[1]
    
    # get sex information
    sex <- as.integer(subset(data, ELF.ID == ind & SVL..cm. == SVL$SVL..cm.[1], select = M)$M[1]) + 1
    
    # if individual was captured twice in one year (SVL is >1 value)... take the average SVL
    if(length(SVL$SVL..cm.) > 1){
      SVL <- mean(as.numeric(SVL$SVL..cm.))
    }
    
    ind <- paste0("ELF_", ind)
    # If SVL is less than 50 (fairly accurate inference), use growth curve to approx age
    if(is.na(subset(metaData, ID == ind, select = BirthYear)) & SVL <= 50 & !is.na(sex) & !is.na(as.numeric(SVL))){
      print(paste0("updating metaData for individual ", ind, " i is ", i))
      inferred_by <- round(infer_birth_year_wSVL(ind, "cass", ind_SVL = as.numeric(SVL), ind_year = cap_year, ind_sex = sex)[1])
      # get +/- 2 year range around inferred birth year
      i_by_max <- inferred_by + 2
      i_by_min <- inferred_by - 2
      
      # update PCC_LifeHistData with inferred information
      if(is.na(subset(metaData, ID == ind, select = BY.max))){ # if by.max is NA, replace with i_by_max
        metaData[metaData$ID == ind,"BY.max"] <- i_by_max
      }else if(i_by_max < subset(metaData, ID == ind, select = BY.max)){ # if inferred by.max is less than original (hypothetically CLOSER to by)
        metaData[metaData$ID == ind,"BY.max"] <- i_by_max
      }
      if(is.na(subset(metaData, ID == ind, select = BY.min))){ # if by.min is NA, replace with i_by_min
        metaData[metaData$ID == ind,"BY.min"] <- i_by_min
      }else if(i_by_min > subset(metaData, ID == ind, select = BY.min)){ # if by.min is not NA, replace with i_by_min if i_by_min is greater than original (hypothetically CLOSER to by)
        metaData[metaData$ID == ind,"BY.min"] <- i_by_min
      }
    }
  }
}

#--------------------------------------------------------------------------------------------------------------------

# get stats for updated dataframe

sum(is.na(metaData$BY.min) & is.na(metaData$BirthYear)) # 20 inds have no by.min, 68 more individuals now have at least some more metaData! 
sum(is.na(metaData$BY.max)& is.na(metaData$BirthYear)) # 0 inds have no by.max

# prep for sequoia 

metaData <- metaData[,-6]

test <- rbind.data.frame(metaData, metaData[rep(which(metaData$ID == "ELF_544"), 2),])
test2 <- test[-which(test$ID == "ELF_544")[1],]
test2[which(test2$ID == "ELF_544")[1],]$ID <- "ELF_544_P5"
test2[which(test2$ID == "ELF_544")[1],]$ID <- "ELF_544_P10"

ELF_LifeHistDataInfer <- test2

# save metadata----------------------------------------------------------------------------------------------------------------------------------------------
save(ELF_LifeHistDataInfer, file = paste0("../pedigree_reconstruction/ELF_metaData_infer_DRB",vcf_date,".Robj"))

# -----------------------------------------------------------------------------------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------------------------------------------------------------------------------

## Cass County: get expanded metadata -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# use seq_inds and data as loaded above with capture events with non-numeric IDs removed

# isolate rows contain data for sequenced individual 

exp_metaData <- vector(mode = "list", length = length(seq_inds))
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
  names(exp_metaData)[[i]] <- paste0("ELF_", seq_inds[i])
  exp_metaData[[i]] <- info
} # end loop

# done! 

ELF_exp_metaData <- exp_metaData
save(ELF_exp_metaData, file = "../pedigree_exploration/ELF_expanded_metaData.Robj")

# -----------------------------------------------------------------------------------------------------------------------------------------------------------

## Cass County: known relatives ----------------------------------------------------------------------------------------------------------------------------
data[which(data$ID == seq_inds[31]),]

# -----------------------------------------------------------------------------------------------------------------------------------------------------------

