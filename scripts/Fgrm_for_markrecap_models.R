## File name: Fgrm_for_markrecap_models.R
## Purpose: Load genotype data and calculate Fgrm for mark recapture models
## Based on vcfR code by R.H. Toczydlowski
## M. I. Clark, October 2023
## Last updated: 02/08/2024

## loading ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
## load libraries
library(vcfR)
library(tidyverse)
library(tidyfast)
library("MetBrewer")
library(adegenet)
library(popgenstuff) # bradburd lab dir 
library(sequoia)
library(conStruct)
library(inbreedR)

load(file = "../vcf_filtering/pop_gen_stats_02082024.Robj", verbose = T)

load(file = "../pedigree_reconstruction/PCC_raw_metadata_2013_2020.Robj", verbose = T)
PCC_data <- data

load(file = "../pedigree_reconstruction/ELF_rawmetaData_,05262023.Robj", verbose = T)
ELF_data <- metaData

# ------------------------------------------------------------------------------------------------------------
pop_gen_stats$id
pop_gen_stats$Fgrm

# ------------------------------------------------------------------------------------------------------------
# match with metadata IDs

data <- data.frame(matrix(data = -9, nrow = length(pop_gen_stats$id), ncol = 4))
colnames(data) <- c("Site", "DNA_ID", "Site_ID", "Fgrm")

data$DNA_ID <- pop_gen_stats$id
data$Fgrm <- pop_gen_stats$Fgrm

data$Site <- unlist(lapply(strsplit(data$DNA_ID, "_"), FUN = function(x){x[1]}))

data[data$Site == "ELF", "Site_ID"]
data[data$Site == "ELF", "Site_ID"] <- ELF_data[match(data[data$Site == "ELF", "DNA_ID"], paste0("ELF_", ELF_data$ID)),"snake_id"] # match returns a vector of the positions of (first) matches of its first argument in its second.
data[data$Site == "PCC", "Site_ID"] <- PCC_data[match(data[data$Site == "PCC", "DNA_ID"], PCC_data$ID),"PIT.Tag.."] # match returns a vector of the positions of (first) matches of its first argument in its second.

# ------------------------------------------------------------------------------------------------------------

# Save data objects ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

save(data, file = "../inbreeding_data_02082024.RObj")

write.csv(data, "../inbreeding_data_02082024.csv", row.names = FALSE)

# ------------------------------------------------------------------------------------------------------------


