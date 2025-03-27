
### Make processed data object for EMR inbreeding analyses and publication
# M. Clark 03/17/24

source("./pedigree_funcs.R") 

# [1] load capture metadata 
load(paste0("../pedigree_exploration/ELF_expanded_metaData.Robj"), verbose = T) # ELF_exp_metaData
load(paste0("../pedigree_exploration/PCC_expanded_metaData.Robj"), verbose = T) # PCC_exp_metaData

load(paste0("../pedigree_reconstruction/PCC_metaData_infer_DRB", "02082024",".Robj"), verbose = T) # PCC_lifeHistDataInfer
PCC_lifeHistDataInfer$ID <- str_remove(PCC_lifeHistDataInfer$ID, "_P[0-9]+")
PCC_lifeHistDataInfer <- PCC_lifeHistDataInfer[-286,]

load(paste0("../pedigree_reconstruction/ELF_metaData_infer_DRB", "02082024",".Robj"), verbose = T) # ELF_LifeHistDataInfer
ELF_LifeHistDataInfer$ID <- str_remove(ELF_LifeHistDataInfer$ID, "_P[0-9]+")
ELF_LifeHistDataInfer <- ELF_LifeHistDataInfer[-779,]

# [2] read in raw data to get mom info for ELF pedigree filtering
ELF_rawData <- read.csv("~/Desktop/EMR_rapture/pedigree_reconstruction/ELF_LifeHistData.csv")[,1:45]

# [3] Load birth year estimates and years contributing from DRB
# load("../pedigree_reconstruction/DRB_byEstimates.Robj", verbose = T) # export_temp, repeated above in life history data objects

# [4] load genetic data 
  # process genetic data in a separate script 
load("../pedigree_reconstruction/PCC_filt_gt_ped.Robj") #PCC_filt_gt_ped
load("../pedigree_reconstruction/ELF_filt_gt_ped.Robj") # ELF_filt_gt_ped
# use t_gt_flt from filter_vcf_popgen.R instead of these, some technical dupe individuals renamed 

# [5] organize data 

# combine genotype data 
gt_data <- rbind.data.frame(t(PCC_filt_gt_ped), t(ELF_filt_gt_ped))

# make empty data object
data <- data.frame(matrix(nrow = nrow(gt_data), ncol = (ncol(gt_data) + 6)))
colnames(data) <- c("id", "site", "sex", "birthYear", "BY.min", "BY.max", colnames(gt_data))

for(i in 1:nrow(data)){
  ind <- rownames(gt_data)[i] # get individual names
  
  # if ind is from PCC/Barry
  if(grepl("PCC", ind)){
    # replace site code with county code
  
    # find metadata entry 
    data[i, "site"] <- "BAR"
    data[i, "sex"] <- subset(PCC_lifeHistDataInfer, ID == ind)$Sex
    data[i, "birthYear"] <- subset(PCC_lifeHistDataInfer, ID == ind)$BirthYear
    data[i, "BY.min"] <- subset(PCC_lifeHistDataInfer, ID == ind)$BY.min
    data[i, "BY.max"] <- subset(PCC_lifeHistDataInfer, ID == ind)$BY.max
    
    # extract and save data 
    data[i, "id"] <- str_replace(ind, "PCC", "BAR")
  }
  
  # if ind is from ELF/Cass
  if(grepl("ELF", ind)){
    # find metadata entry 
    data[i, "site"] <- "CAS"
    data[i, "sex"] <- subset(ELF_LifeHistDataInfer, ID == ind)$Sex
    data[i, "birthYear"] <- subset(ELF_LifeHistDataInfer, ID == ind)$BirthYear
    data[i, "BY.min"] <- subset(ELF_LifeHistDataInfer, ID == ind)$BY.min
    data[i, "BY.max"] <- subset(ELF_LifeHistDataInfer, ID == ind)$BY.max
    
    # extract and save data 
    data[i, "id"] <- str_replace(ind, "ELF", "CAS")
    
  }
  
  # add genotype data from combine object
  
}

# adjust coordinates
PCC_coords <- get_coords(exp_metaData =PCC_exp_metaData, target_inds = colnames(PCC_filt_gt_ped), easting_col ="DD.easting", northing_col ="DD.northing", doy = FALSE)
PCC_coords$site <- "BAR"
PCC_coords$ID <- str_replace(PCC_coords$ID, "PCC", "BAR")
PCC_coords <- PCC_coords[, c("ID", "site", "easting", "northing")]

ELF_coords <- get_coords(exp_metaData =ELF_exp_metaData, target_inds = colnames(ELF_filt_gt_ped), easting_col ="UTM.easting", northing_col ="UTM.northing", doy = TRUE)
ELF_coords$site <- "CAS"
ELF_coords$ID <- str_replace(ELF_coords$ID, "ELF", "CAS")
ELF_coords <- ELF_coords[, c("ID", "site", "easting", "northing")]

# quality filter coordinates: 
PCC_coords_flt <- na.omit(PCC_coords)
PCC_coords_flt$easting[PCC_coords_flt$easting > 0] <- PCC_coords_flt$easting[PCC_coords_flt$easting > 0]*-1 # easting should be negative 

PCC_coords_flt <- PCC_coords_flt[-which(PCC_coords_flt$easting < -85.3583),] # PCC_268, outside study bounds
PCC_coords_flt <- PCC_coords_flt[-which(PCC_coords_flt$northing < 42.518),] # PCC_303, outside study bounds

ELF_coords_flt <- na.omit(ELF_coords)
ELF_coords_flt <- ELF_coords_flt[-which(ELF_coords_flt$easting > 583500),] # ELF_361, coord is ecolab

# Shift coordinates
set.seed(13523456)

PCC_coords_flt$easting <- PCC_coords_flt$easting * runif(1, -1, 1)
PCC_coords_flt$northing <- PCC_coords_flt$northing * runif(1, -1, 1)
colnames(PCC_coords_flt) <- c("id", "site", "easting_shifted", "northing_shifted")

ELF_coords_flt$easting <- ELF_coords_flt$easting * runif(1, -1, 1)
ELF_coords_flt$northing <- ELF_coords_flt$northing * runif(1, -1, 1)
colnames(ELF_coords_flt) <- c("id", "site", "easting_shifted", "northing_shifted")

BAR_coords <- PCC_coords_flt
CAS_coords <- ELF_coords_flt
# Output: 
  # (1) data frame with DNA ID | site | sex | birth year estimates | genotypes
  # (2) data frame with DNA ID | site | capture coordinates (long format)

# Save

data_list <- list(data, BAR_coords, CAS_coords)
save(data_list, file = "../final_data_objects.Robj")


