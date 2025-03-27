## File name: extract_PCC_metaData.R
## Purpose: Pull out metadata for pedigree exploration from field recapture data
## M. I. Clark, November 2022
## Last updated: 11/15/2022

## note: most NA warnings are okay to ignore

## Loading ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

# load combined metadata from read_PCC_metaData.R
load("../pedigree_reconstruction/PCC_raw_metadata_2013_2020.Robj", verbose = TRUE)
        # Loading objects:
        #   seq_inds
        #   data

# ------------------------------------------------------------------------------------------------------------

## Wrangle metadata ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

# isolate rows contain data for sequenced individual 

metaData <- vector(mode = "list", length = length(seq_inds))
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
    print(paste0("No metadata for ", seq_ind, " with index ", i))
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
  names(metaData)[[i]] <- paste0(seq_inds[i])
  metaData[[i]] <- info
} # end loop

# done! 

# ------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------
# Save ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
PCC_exp_metaData <- metaData
save(PCC_exp_metaData, file = "../pedigree_exploration/PCC_expanded_metaData.Robj")

# ------------------------------------------------------------------------------------------------------------

