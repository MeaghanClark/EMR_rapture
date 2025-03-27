## File name: make_pedigree.R
## Purpose: Make pedigrees for EMR project
## M. I. Clark, April 2023
## Last updated: 06/26/2024
  # updating birth years to use Danielle's estimates

## Loading ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

## load libraries
library(sequoia) # 2.5.3
library(kinship2)

## date
# date = 040502023
# date = "05262023"
# date = "02082024"
# date = "04132024"
# date = "05072024"
date = "06262024"

## load data
# PCC
load(paste0("../vcf_filtering/PCC_filtered_gt_all_snps_noZ_", "02082024", ".Robj"), verbose = T) # PCC_filt_gt_all

load(paste0("../pedigree_reconstruction/PCC_metaData_infer_DRB", "02082024",".Robj"), verbose = T) # PCC_LifeHistData 
  # "original" has added metadata for missing individuals, but not SVL inference
  # "infer" version has additional metadata inferred using SVL 
  # "infer_DRB" is using the birth years inferred by Danielle

# ELF
load(paste0("../vcf_filtering/ELF_filtered_gt_all_snps_noZ_", "02082024", ".Robj"), verbose = T) # ELF_filt_gt_all

load(paste0("../pedigree_reconstruction/ELF_metaData_infer_DRB", "02082024",".Robj"), verbose = T) # ELF_LifeHistDataInfer
# "original" has added metadata for missing individuals, but not SVL inference

# duplicate individuals for PCC
PCC_dupe_inds <- read.csv("../pedigree_reconstruction/PCC_LifeHistData_dupes.csv", header=FALSE)
# only one known duplicate for ELF
ELF_inds <- colnames(ELF_filt_gt_all)
ELF_dupes <- ELF_inds[which(grepl("P", ELF_inds))]

dupe_inds <- rbind.data.frame(PCC_dupe_inds, c("ELF_544_P5", "ELF_544_P10", ""))

# "extra" metadata for analyses
load(paste0("../pedigree_exploration/ELF_expanded_metaData.Robj"), verbose = T) # ELF_exp_metaData
load(paste0("../pedigree_exploration/PCC_expanded_metaData.Robj"), verbose = T) # PCC_exp_metaData

# load inferred age information for ambiguous pairs from previous pedigree run
# generated in pedigree_size_exploration.R 
# load("../pedigree_reconstruction/PCC_inferred_by_04132024.Robj", verbose = T)
# load("../pedigree_reconstruction/ELF_inferred_by_04132024.Robj", verbose = T)

# load genetic data to maintain continuity in ids 
load(file = "../vcf_filtering/pop_gen_stats_02082024.Robj", verbose = T)

# ------------------------------------------------------------------------------------------------------------

## Define custom functions ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

# transpose data and replace NAs with -9 
prep_gt_data <- function(filt_gt){ 
  gt <- t(filt_gt)
  gt[is.na(gt)] <- -9 # replace NAs with -9 
  return(gt)
}

# prepares metadata for pedigree reconstruction
prep_metaData <- function(metaData, inds){
  meta <-  as.data.frame(matrix(data=NA, ncol= ncol(metaData), nrow=length(inds))) # create data frame
  for (i in 1:length(inds)){
    index <- match(inds[i], metaData$ID) # matches ind ids from snp data to metadata
    meta[i,] <- metaData[index,]
  }
  colnames(meta) <- colnames(metaData)
  #meta[is.na(meta)] <- NA # replace NAs with -9
  meta[which(is.na(meta$Sex)), "Sex"] <- -9
  meta[which(is.na(meta$BY.max)), "BY.max"] <- -9
  meta[which(is.na(meta$BY.min)), "BY.min"] <- -9
  return(meta)
}

# calculates genotyping error rate between duplicate individuals
calc_gt_error_rate <- function(dupe_ind_ids, genotypes, tri = TRUE){
  if(tri == TRUE){
    for (real_dupes in 1:nrow(dupe_ind_ids)){
      if(grepl("PCC_[:digit:]*", dupe_ind_ids[real_dupes,3])==TRUE){ # if this is a triplicate
        ind1 <- dupe_ind_ids[real_dupes,1]
        ind2 <- dupe_ind_ids[real_dupes,2]
        ind3 <- dupe_ind_ids[real_dupes,3]
        p1 <- c(ind1, ind2, "")
        p2 <- c(ind1, ind3, "")
        p3 <- c(ind2, ind3, "")
        dupe_ind_ids <- rbind.data.frame(dupe_ind_ids, p1, p2, p3)
      }
    }
    for (real_dupes in 1:nrow(dupe_ind_ids)){ # remove original triplicate row 
      if(grepl("PCC_[:digit:]*", dupe_ind_ids[real_dupes,3])==TRUE){
        dupe_ind_ids <- dupe_ind_ids[-real_dupes,]
      }
    }
  }
  error_rate <- vector(length = nrow(dupe_ind_ids))
  for(j in 1:nrow(dupe_ind_ids)){
    sites <- 0
    same <- 0
    dupes <- as.character(dupe_ind_ids[j,])
    for(k in 1:ncol(genotypes)){
      if(genotypes[which(rownames(genotypes) == dupes[1]),k] != -9){
        if(genotypes[which(rownames(genotypes) == dupes[2]),k] != -9){
          sites = sites + 1 
          if(genotypes[which(rownames(genotypes) == dupes[1]),k] == genotypes[which(rownames(genotypes) == dupes[2]),k]){
            same = same + 1
          }
        }
      }
    }
    error_rate[j] <- (1-(same/sites))
  }
  return(error_rate)
}

# returns a dataframe describing which confirmed duplicate pairs were detected by sequoia 
detect_dupe_pairs <- function(pairs, ped.out, tri = TRUE){
  # testing variables 
  #   pairs <- PCC_dupe_inds
  #   ped.out <-  ped.PCC.dup
  
  matches <- vector(length = nrow(pairs))
  mismatches <- vector(length = nrow(pairs))
  for (real_dupes in 1:nrow(pairs)){
    if(grepl("PCC_[:digit:]*", pairs[real_dupes,3])==TRUE){ # if this is a triplicate
      ind1 <- pairs[real_dupes,1]
      ind2 <- pairs[real_dupes,2]
      ind3 <- pairs[real_dupes,3]
      p1 <- c(ind1, ind2, "")
      p2 <- c(ind1, ind3, "")
      p3 <- c(ind2, ind3, "")
      pairs <- rbind.data.frame(pairs, p1, p2, p3)
    }
  }
  for (real_dupes in 1:nrow(pairs)){ # remove original triplicate row 
    if(grepl("PCC_[:digit:]*", pairs[real_dupes,3])==TRUE){
      pairs <- pairs[-real_dupes,]
    }
  }
  for (real_dupes in 1:nrow(pairs)){
    matches[real_dupes] <- unlist(match_pairs(pairs = pairs, real_dupes = real_dupes, ped.out = ped.out))[1]
    mismatches[real_dupes] <- unlist(match_pairs(pairs = pairs, real_dupes = real_dupes, ped.out = ped.out))[2]
  }
  return(cbind.data.frame(pairs, matches, mismatches))
}

match_pairs <- function(pairs, real_dupes, ped.out){
  # testing vars
    # pairs = pairs
    # real_dupes = 24
    # ped.out
  ids <- as.character(pairs[real_dupes,])
  ind1 <- ids[1]
  ind2 <- ids[2]
  if(nrow(ped.out$DupGenotype[which(ind1 == ped.out$DupGenotype[,3]),]) >= 1){ # first id column
    if(nrow(ped.out$DupGenotype[which(ind2 == ped.out$DupGenotype[,4]),]) >= 1){# second id should be in the second id column
      match <- TRUE
      mismatch <- ped.out$DupGenotype[which(ind1 == ped.out$DupGenotype[,3] & ind2 == ped.out$DupGenotype[,4]),5]
    }else{
      match <- FALSE
      mismatch <- NA
    }
  }else if(nrow(ped.out$DupGenotype[which(ind1 == ped.out$DupGenotype[,4]),]) >= 1){ # second id column
    if(nrow(ped.out$DupGenotype[which(ind2 == ped.out$DupGenotype[,3]),]) >= 1){# second id should be in the first id column
      match <- TRUE
      mismatch <- ped.out$DupGenotype[which(ind1 == ped.out$DupGenotype[,4] & ind2 == ped.out$DupGenotype[,3]),5]
    }else{
      match <- FALSE
      mismatch <- NA
    }
  }
  return(list(match, mismatch))
}

# ------------------------------------------------------------------------------------------------------------

## Prepare SNP data ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
#PCC
PCC_prep_gt_all <- prep_gt_data(PCC_filt_gt_all)

PCC_SnpStats <- SnpStats(PCC_prep_gt_all)
dim(PCC_SnpStats)
dim(PCC_prep_gt_all)
PCC_gt_for_ped <- PCC_prep_gt_all[,which(PCC_SnpStats[,1] >= 0.1)] # 288 1452

#ELF
ELF_prep_gt_all <- prep_gt_data(ELF_filt_gt_all)
ELF_SnpStats <- SnpStats(ELF_prep_gt_all)
ELF_gt_for_ped <- ELF_prep_gt_all[,which(ELF_SnpStats[,1] >= 0.1)] # .1: 779 1472

ELF_final_snp_stats <- SnpStats(ELF_gt_for_ped)
#CalcMaxMismatch(error_rate, MAF = ELF_final_snp_stats[,1])
# ------------------------------------------------------------------------------------------------------------

## Prepare meta data ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

# # concatenate inferred birth years
# for(i in 1:nrow(PCC_PO_by)){
#   if(!is.na(PCC_PO_by$by.min[i])){
#     PCC_LifeHistData[which(PCC_LifeHistData$ID == PCC_PO_by$ID[i]),"BY.min"] <- PCC_PO_by$by.min[i]
#     PCC_LifeHistData[which(PCC_LifeHistData$ID == PCC_PO_by$ID[i]),"BY.max"] <- PCC_PO_by$by.max[i]
#   }
# }
# 
# for(i in 1:nrow(ELF_PO_by)){
#   if(!is.na(ELF_PO_by$by.min[i])){
#     ELF_LifeHistData[which(ELF_LifeHistData$ID == ELF_PO_by$ID[i]),"BY.min"] <- ELF_PO_by$by.min[i]
#     ELF_LifeHistData[which(ELF_LifeHistData$ID == ELF_PO_by$ID[i]),"BY.max"] <- ELF_PO_by$by.max[i]
#   }
# }


#PCC
# PCC_meta <- prep_metaData(PCC_LifeHistData, rownames(PCC_prep_gt_all))
#PCC_meta <- PCC_LifeHistData # Sequoia version wants NAs not -9s
PCC_meta <- PCC_lifeHistDataInfer
#ELF
#ELF_meta <- prep_metaData(ELF_LifeHistData, rownames(ELF_prep_gt_all))
ELF_meta <- ELF_LifeHistDataInfer
# ------------------------------------------------------------------------------------------------------------

## Estimate genotyping error rate ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

# mean error rate from filter_vcf.R is 0.004136376 shouldn't I use the most SNPs? 
# Is filtering based on gt error rate even valid? I don't think it is. Assumes that SNPs with high 
# error rates in duplicated individuals have high error over all. I think it's more conservative to 
# use the overall error rate... .4% is actually not that high! 
error_rate = 0.002348661


# ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

## PEDIGREE RECONSTRUCTION

# ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

## PCC ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

# input data: 
# SNP: PCC_prep_gt_all_ger_cutoff
# Meta: PCC_meta
# Descriptions from Sequoia (https://jiscah.github.io/articles/vignette_main/book/running-sequoia.html)

### 1: 'PRE' --------------------------------------------------------------------
# check that the genotype data, life history data and input parameters are in a valid format, create 'Specs' element of output list

#
ped.PCC.pre <- sequoia(GenoM = PCC_gt_for_ped, LifeHistData=PCC_meta, Module = "pre", Err = error_rate)


### 2: 'DUP' --------------------------------------------------------------------
# Check for (nearly) identical genotypes, and for duplicated IDs in the genotype and life history data

ped.PCC.dup <- sequoia(GenoM = PCC_gt_for_ped, LifeHistData=PCC_meta, Module = "dup", Err = error_rate)

# There were 30 likely duplicate genotypes found, consider removing

#### remove duplicates
# PCC_59    PCC_217        not detected

detect_dupe_pairs(PCC_dupe_inds, ped.PCC.dup) # all known duplicate pairs are detected including PCC_66 and PCC_44, both are the same as PCC_179


#### remove known duplicates
PCC_dupe_inds

# sample individuals to throw out at random
# toss_inds <- vector()
# for(i in 1:nrow(PCC_dupe_inds)){
#   ind1 <- PCC_dupe_inds[i,1]
#   ind2 <- PCC_dupe_inds[i,2]
#   ind3 <- PCC_dupe_inds[i,3]
#   if(nchar(ind3) == 0){ # if there are two sets of the same ind
#     # check for missing metadata for first sample
#     ind1_meta <- PCC_meta[which(PCC_meta[,1] == PCC_dupe_inds[i,1]),]
#     # check for missing metadata for second sample
#     ind2_meta <- PCC_meta[which(PCC_meta[,1] == PCC_dupe_inds[i,2]),]
#     if(sum(ind1_meta == ind2_meta, na.rm = TRUE) == 4){ # if the metadata for both samples is the same... 
#       toss <- as.character(sample(PCC_dupe_inds[i,1:2], 1)) # sample an individual
#       toss_inds <- c(toss_inds, toss) # and add it to the list of individuals to remove
#     }else{ # if the metadata for both samples is not the same... 
#       if (sum(!is.na(ind1_meta)) > sum(!is.na(ind2_meta))){ # if ind 1 has more complete metadata... 
#         toss <- as.character(ind2_meta$ID) # toss ind 2
#         toss_inds <- c(toss_inds, toss)
#       }else if (sum(!is.na(ind2_meta)) > sum(!is.na(ind1_meta))){ # if ind 2 has more complete metadata... 
#         toss <- as.character(ind1_meta$ID) # toss ind 1
#         toss_inds <- c(toss_inds, toss)
#       } else if (!is.na(ind1_meta$BY.max) & !is.na(ind2_meta$BY.max ) & ind1_meta$BY.max != ind2_meta$BY.max){ # if the samples have different BY.max
#         if(ind1_meta$BY.max > ind2_meta$BY.max){ # if ind 1 has a higher BY.max than ind 2
#           toss <- as.character(ind1_meta$ID) # toss ind 1, keep the more conservative BY.max
#           toss_inds <- c(toss_inds, toss)
#         }else if (ind2_meta$BY.max > ind1_meta$BY.max){ # if ind 2 has a higher BY.max than ind 1
#           toss <- as.character(ind2_meta$ID) # toss ind 2, keep the more conservative BY.max
#           toss_inds <- c(toss_inds, toss)
#         }else{
#           print(paste0("Attention to i = ", i))
#           print(paste0(ind1_meta, " ", ind2_meta))
#         }
#       }
#     }
#   }else if(nchar(ind3) > 0){ # if there are three sets of the same ind
#     # check for missing metadata for first sample
#     ind1_meta <- PCC_meta[which(PCC_meta[,1] == PCC_dupe_inds[i,1]),]
#     # check for missing metadata for second sample
#     ind2_meta <- PCC_meta[which(PCC_meta[,1] == PCC_dupe_inds[i,2]),]
#     # check for missing metadata for third sample
#     ind3_meta <- PCC_meta[which(PCC_meta[,1] == PCC_dupe_inds[i,3]),]
#     if(sum(ind1_meta == ind2_meta, na.rm = TRUE) == 4 && sum(ind1_meta == ind3_meta) == 4){ # if the metadata for all samples is the same... 
#       toss <- as.character(sample(PCC_dupe_inds[i,1:3], 2)) # sample two individuals
#       toss_inds <- c(toss_inds, toss) # and add it to the list of individuals to remove
#     }else{ # if the metadata for both samples is not the same... 
#       print(paste0("Attention to i = ", i))
#       print(paste0(ind1_meta, " ", ind2_meta, " ", ind3_meta))
#     }
#   }
# }

# toss same individuals as previous runs 

rownames(PCC_gt_for_ped) %in% subset(pop_gen_stats, site == "PCC", select = id)$id



PCC_gt_for_ped_drop <- PCC_gt_for_ped[rownames(PCC_gt_for_ped) %in% subset(pop_gen_stats, site == "PCC", select = id)$id,]


# load(paste0("../pedigree_reconstruction/PCC_pedigree_results_", "04132024", ".Robj"), verbose = T)
# PCC_dupe_vec <- c(PCC_dupe_inds[,1], PCC_dupe_inds[,2], "PCC_179", "PCC_200")
# 
# toss_inds <- PCC_dupe_vec[!PCC_dupe_vec %in% PCC_results[[4]]$Pedigree$id]
# 
# 
# PCC_gt_for_ped_drop <- PCC_gt_for_ped[-match(toss_inds, rownames(PCC_gt_for_ped)),] # 262 1537 262 2223 dim: 263 2468 
# 
# ped.PCC.dup.2 <- sequoia(GenoM = PCC_gt_for_ped_drop, LifeHistData=PCC_meta, Module = "dup", Err = error_rate)
# 
# #### remove inferred duplicates
# 
# # what are the other pairs of dupes?
# 
# ped.PCC.dup.2$DupGenotype
# 
# toss_inferred <- vector()
# for(i in 1:nrow(ped.PCC.dup.2$DupGenotype)){
#   if(sum(ped.PCC.dup.2$DupGenotype[i,c(3,4)] %in% as.vector(unlist(PCC_dupe_inds))) == 0){
#     pair <- ped.PCC.dup.2$DupGenotype[i,c(3,4)]
#     ind1 <- ped.PCC.dup.2$DupGenotype[i,3]
#     ind2 <- ped.PCC.dup.2$DupGenotype[i,4]
#     # check for missing metadata for first sample
#     ind1_meta <- PCC_meta[which(PCC_meta[,1] == ind1),]
#     # check for missing metadata for second sample
#     ind2_meta <- PCC_meta[which(PCC_meta[,1] == ind2),]
#     if(sum(ind1_meta == ind2_meta, na.rm = TRUE) == 4){ # if the metadata for both samples is the same... 
#       toss <- as.character(sample(pair, 1)) # sample an individual
#       toss_inferred <- c(toss_inferred, toss)
#     }else{ # if the metadata for both samples is not the same... 
#       if (sum(!is.na(ind1_meta)) > sum(!is.na(ind2_meta))){ # if ind 1 has more complete metadata... 
#         toss <- as.character(ind2_meta$ID) # toss ind 2
#         toss_inferred <- c(toss_inferred, toss)
#       }else if (sum(!is.na(ind2_meta)) > sum(!is.na(ind1_meta))){ # if ind 2 has more complete metadata... 
#         toss <- as.character(ind1_meta$ID) # toss ind 1
#         toss_inferred <- c(toss_inferred, toss)
#       } else if (!is.na(ind1_meta$BY.max) & !is.na(ind2_meta$BY.max) &ind1_meta$BY.max != ind2_meta$BY.max){ # if the samples have different BY.max
#         if(ind1_meta$BY.max > ind2_meta$BY.max){ # if ind 1 has a higher BY.max than ind 2
#           toss <- as.character(ind1_meta$ID) # toss ind 1, keep the more conservative BY.max
#           toss_inferred <- c(toss_inferred, toss)
#         }else if (ind2_meta$BY.max > ind1_meta$BY.max){ # if ind 2 has a higher BY.max than ind 1
#           toss <- as.character(ind2_meta$ID) # toss ind 2, keep the more conservative BY.max
#           toss_inferred <- c(toss_inferred, toss)
#         }else{
#           print(paste0("Attention to i = ", i))
#           print(paste0(ind1_meta, " ", ind2_meta))
#         }
#       }
#     }
#   }
#   print(paste0("For individuals ", pair[1], " and ", pair[2], ", ", toss, " is added to toss list. i is ", i))
# }

# remove inferred duplicates that are actually different individuals from list
# but I want to remove the individual with the least metadata

# need to check these in metadata! 
# DIFFERENT "For individuals PCC_210 and PCC_99, PCC_99 is added to toss list. i is 11"
# 210 : 836852544, F
# 99: 027327627, F
# 1 SNP mismatch 
# SAME "For individuals PCC_202 and PCC_295, PCC_295 is added to toss list. i is 12"
# 202: 836329312 M
# 295: no PIT ID or sex, 
# 6 mismatch, hard to tell but likely same individual/no evidence of different individuals
# SAME "For individuals PCC_290 and PCC_314, PCC_314 is added to toss list. i is 33"
# 290: 600784009, F
# 314: 600784009 or 600789009 (hard to read on blood vial)
# 31 mismatches likely same pit tag, likely same individual 

# updated 01292024: dropping 295, 290
# updated 05252023: made the decision to drop PCC_283, PCC_86, and PCC_99
# 283/4: different PIT tags, 149/86: different pit tags, 210/99: different pit tags
#incorrect_dupes <- c("PCC_283", "PCC_4", "PCC_149", "PCC_86", "PCC_210", "PCC_99")
#toss_inferred <- toss_inferred[is.na(match(toss_inferred, incorrect_dupes))]

# 02062024: 
# [1] "For individuals PCC_202 and PCC_295, PCC_295 is added to toss list. i is 1"
# [1] "For individuals PCC_290 and PCC_314, PCC_290 is added to toss list. i is 2"

# 04132024: PCC_290 and PCC_314 not inferred.. dropping PCC_290 because individuals are likely the same 
# toss_inferred <- c("PCC_295", "PCC_290")
# PCC_gt_for_ped_drop2 <- PCC_gt_for_ped_drop[-match(toss_inferred, rownames(PCC_gt_for_ped_drop)),] # 260 1537
# 
# 
# # verify that ids match genetic ids 
# #individuals in pop_gen_stats that are not in the ped list
# subset(pop_gen_stats, site == "PCC", select = id)$id[!subset(pop_gen_stats, site == "PCC", select = id)$id %in% rownames(PCC_gt_for_ped_drop2)]
# # [1] "PCC_140"    "PCC_156"    "PCC_16"     "PCC_17"     "PCC_188"    "PCC_1"      "PCC_206"    "PCC_217"   
# # [9] "PCC_255"    "PCC_42"     "PCC_57"     "PCC_62_P10"
# 
# check_ids <- function(ind, genetic_data, pedigree_ids){
  # ind = "PCC_140"
  # genetic_data = pop_gen_stats
  # pedigree_ids = PCC_gt_for_ped_drop2
#   if (nrow(subset(genetic_data, id == ind)) > 0 ){
#     print(paste0("ID ", ind, " is in genetic data"))
#   }
#   if(sum(rownames(PCC_gt_for_ped_drop2)==ind) == 0){
#     print(paste0("ID ", ind, " is not in pedigree list"))
#   }
#   if(sum(PCC_dupe_inds == ind) == 1){
#     print(paste0("ID ", ind, " is a duplicate"))
#   }
# }
# 
# lapply(subset(pop_gen_stats, site == "PCC", select = id)$id[!subset(pop_gen_stats, site == "PCC", select = id)$id %in% rownames(PCC_gt_for_ped_drop2)], 
#        FUN = check_ids, genetic_data = pop_gen_stats, pedigree_ids = rownames(PCC_gt_for_ped_drop2))
# # ids don't match because of duplicates selected 
# 
# # individuals in ped list that are not in pop_gen_stats
# rownames(PCC_gt_for_ped_drop2)[!rownames(PCC_gt_for_ped_drop2) %in% subset(pop_gen_stats, site == "PCC", select = id)$id]
# subset(pop_gen_stats, id == "PCC_105")
# sum(rownames(PCC_gt_for_ped_drop2) == "PCC_105")
# sum(PCC_dupe_inds == "PCC_105")
# lapply(rownames(PCC_gt_for_ped_drop2)[!rownames(PCC_gt_for_ped_drop2) %in% subset(pop_gen_stats, site == "PCC", select = id)$id], 
#        FUN = check_ids, genetic_data = pop_gen_stats, pedigree_ids = rownames(PCC_gt_for_ped_drop2))
# 

### 3: 'PAR' --------------------------------------------------------------------
# Assign genotyped parents to genotyped individuals. Includes call to MakeAgePrior() to estimate AgePriors based on the just-assigned parents.

ped.PCC.par <- sequoia(GenoM = PCC_gt_for_ped_drop, LifeHistData=PCC_meta, Module = "par", Err = error_rate)

# 10182022 : assigned 118 dams and 88 sires to 267 individuals
# 11142022 : assigned 110 dams and 87 sires to 264 individuals
# 02032023 : assigned 113 dams and 94 sires to 264 individuals
# 02102023 : assigned 121 dams and 89 sires to 264 individuals
---
# 04052023 : assigned 94 dams and 73 sires to 264 individuals
# 05252023 : assigned 95 dams and 74 sires to 261 individuals
# 05262023 : assigned 95 dams and 74 sires to 261 individuals 
# 01292024 : assigned 95 dams and 74 sires to 260 individuals 2k snps
# 01292024 : assigned 94 dams and 73 sires to 260 individuals 1k snps
# 01292024 : assigned 95 dams and 74 sires to 260 individuals 1.5k snps
# 02062024 : assigned 95 dams and 74 sires to 260 individuals
# 04132024 : assigned 102 dams and 77 sires to 260 individuals
# 05072024 ORIGINAL : assigned 96 dams and 78 sires to 260 individuals
# 05072024 ORIGINAL, fixed metadata bugs : assigned 96 dams and 78 sires to 260 individuals
# 05082024 INFER, fixed metadata bugs : assigned 104 dams and 79 sires to 260 individuals
# 06262024 INFER DRB : assigned 117 dams and 80 sires to 260 individuals

### 4: 'PED' --------------------------------------------------------------------
# Cluster half- and full-siblings and assign each cluster a dummy-parent; assign grandparents to sibships and singletons.

ped.PCC.ped <- sequoia(GenoM = PCC_gt_for_ped_drop, LifeHistData=PCC_meta, Module = "ped", Err = error_rate)

ped.PCC.ped$Pedigree

# 10182022 : assigned 158 dams and 150 sires to 267 + 27 individuals (real + dummy)
# 11142022 : assigned 156 dams and 155 sires to 264 + 28 individuals (real + dummy)
# 02032023: assigned 158 dams and 154 sires to 264 + 31 individuals (real + dummy)
# 02102023: assigned 162 dams and 158 sires to 264 + 32 individuals (real + dummy)
# 04052023: assigned 176 dams and 180 sires to 264 + 45 individuals (real + dummy)
# 05252023: assigned 181 dams and 182 sires to 261 + 47 individuals (real + dummy)
# 05262023: assigned 176 dams and 175 sires to 261 + 45 individuals (real + dummy)
# 01292024: assigned 177 dams and 174 sires to 260 + 45 individuals (real + dummy) 2k snps
# 01292024: assigned 174 dams and 172 sires to 260 + 41 individuals (real + dummy) 1k snps
# 01292024: assigned 179 dams and 177 sires to 260 + 45 individuals (real + dummy) 1.5k snps
# 02062024: assigned 175 dams and 173 sires to 260 + 42 individuals (real + dummy)
# 02062024: assigned 171 dams and 169 sires to 261 + 46 individuals (real + dummy)
# 04132024: assigned 172 dams and 165 sires to 260 + 44 individuals (real + dummy)
# 05072024 ORIGINAL : assigned 171 dams and 168 sires to 260 + 44 individuals (real + dummy)
# 05072024 ORIGINAL updated version : assigned 168 dams and 166 sires to 260 + 45 individuals (real + dummy)
# 05072024 ORIGINAL updated version, fixed metadata bugs : assigned 168 dams and 166 sires to 260 + 45 individuals (real + dummy)
# 05082024 INFER, fixed metadata bugs : assigned 171 dams and 166 sires to 260 + 45 individuals (real + dummy)
# 06262024 INFER DRB : assigned 182 dams and 175 sires to 260 + 49 individuals (real + dummy)

### 5: Rescue pairs  --------------------------------------------------------------------

# if working from previously generated pedigree: 
# load previous run 
# if(load_previous == TRUE){
#   load(paste0("../pedigree_reconstruction/PCC_pedigree_results_", date, ".Robj"), verbose = T)
#   PCC_gt_for_ped_drop2 <- PCC_results[[1]]
#   PCC_meta <- PCC_results[[2]]
#   ped.PCC.ped <- PCC_results[[4]]
#   error_rate <- PCC_results[[5]]
# }

# 05072024 ORIGINAL updated version, fixed metadata bugs: Found 12 likely parent-offspring pairs, and 219 other non-assigned pairs of possible relatives


PCC_maybe_rel_ped <- GetMaybeRel(GenoM = PCC_gt_for_ped_drop2, LifeHistData=PCC_meta, Pedigree = ped.PCC.ped$Pedigree, Module = "ped", Err = error_rate)

PCC_maybe_po <- PCC_maybe_rel_ped$MaybeRel[PCC_maybe_rel_ped$MaybeRel$TopRel == "PO",]

# subset(PCC_meta, ID == "PCC_155")
# subset(PCC_meta, ID == "PCC_169")
# subset(data, ID == "PCC_155")
# subset(data, ID == "PCC_169")
# infer_birth_year_wSVL("PCC_155", "barry", ind_SVL = as.numeric(65.4), ind_year = 2015, ind_sex = 2)
# infer_birth_year_wSVL("PCC_169", "barry", ind_SVL = as.numeric(47.8), ind_year = 2015, ind_sex = 1)

### 6: Pedigree Visualization  --------------------------------------------------------------------
# input data: 
# SNP: PCC_prep_gt_all_drop2
# Meta: PCC_meta

pdf(file = paste0("../pedigree_reconstruction/PCC_snp_stats_", date, ".pdf"), height = 6, width = 10)
out <- SnpStats(GenoM = PCC_gt_for_ped_drop2, ErrFlavour = mean(error_rate))
dev.off()

pdf(file = paste0("../pedigree_reconstruction/PCC_age_prior_", date, ".pdf"), height = 6, width = 6)
MakeAgePrior(ped.PCC.ped$Pedigree, LifeHistData = PCC_meta)
dev.off()

# # inferred dupes w/ conflicting metadata: 
# 
# # "PCC_283" and "PCC_4", "PCC_149" and "PCC_86" , "PCC_210" and "PCC_99" )
# 
# # PCC_283 and PCC_4 (2 mismatches)
# # who are their assigned parents? 
# ped.PCC.ped$Pedigree[which(ped.PCC.par$Pedigree$id == "PCC_283"),] # same assigned parents
# ped.PCC.ped$Pedigree[which(ped.PCC.par$Pedigree$id == "PCC_4"),]
# # are they assigned offspring? 
# ped.PCC.ped$Pedigree[which(ped.PCC.par$Pedigree$dam == "PCC_283"),] # No
# ped.PCC.ped$Pedigree[which(ped.PCC.par$Pedigree$dam == "PCC_4"),] # PCC_223 is assigned as offspring
# 
# # Assigned as full sibs

# Capture time
# PCC_283: initial cap 2018, SVL 53.7, F, A
# PCC_4: initial cap 2013, SVL 44.7, F, J

# PCC_149 and PCC_86 (4 mismatches)
# who are their assigned parents? 
# ped.PCC.ped$Pedigree[which(ped.PCC.par$Pedigree$id == "PCC_149"),] # 
# ped.PCC.ped$Pedigree[which(ped.PCC.par$Pedigree$id == "PCC_86"),] #
# # are they assigned offspring? 
# ped.PCC.ped$Pedigree[which(ped.PCC.par$Pedigree$dam == "PCC_149"),] # no offpsring
# ped.PCC.ped$Pedigree[which(ped.PCC.par$Pedigree$dam == "PCC_86"),] # no offspring

# Assigned as full sibs

# # Capture time
# # PCC_149: initial cap 2015, SVL 41.2, M, J, 027298014
# # PCC_86: initial cap 2014, SVL 38.8, F, J, 027323569
# PCC_meta[which(PCC_meta$ID == "PCC_149"),] # BirthYear = 2015
# PCC_meta[which(PCC_meta$ID == "PCC_86"),] # BY.max = 2011
# # both were captured as juveniles, but have different age estimates
# 
# # PCC_210 and PCC_99 (1 mismatch)
# # who are their assigned parents? 
# ped.PCC.par$Pedigree[which(ped.PCC.par$Pedigree$id == "PCC_210"),] # PCC_76 assigned as sire
# ped.PCC.par$Pedigree[which(ped.PCC.par$Pedigree$id == "PCC_99"),] # PCC_76 assigned as sire,
# # are they assigned offspring? 
# ped.PCC.par$Pedigree[which(ped.PCC.par$Pedigree$dam == "PCC_210"),] # PCC_291 is assigned as offspring
# ped.PCC.par$Pedigree[which(ped.PCC.par$Pedigree$dam == "PCC_99"),] #
# 

# Suspected Relatives 

maybeSibs1 <- c("PCC_102", "PCC_103", "PCC_104", "PCC_111", "PCC_112", "PCC_113", "PCC_119")

maybeSibs2 <- c("PCC_101", "PCC_120", "PCC_127", "PCC_128", "PCC_129", "PCC_138", "PCC_139")

ped.PCC.ped$Pedigree[match(maybeSibs1, ped.PCC.ped$Pedigree$id),]
ped.PCC.ped$Pedigree[match(maybeSibs2, ped.PCC.ped$Pedigree$id),]

# that's cool! some are related, some are not! 

# neonates found together 
#          8922 run                 #81122 run    101822 run
#   ID      dam     sire        dam     sire      dam    sire
# PCC_101   F0012 PCC_161     PCC_140 PCC_161     PCC_36 PCC_161    ### 36 and 140 are the same ind!
# PCC_120    <NA>  PCC_81     PCC_209  PCC_81     PCC_209  PCC_81
# PCC_127  PCC_36 PCC_161     PCC_140 PCC_161     PCC_36 PCC_161
# PCC_138  PCC_36 PCC_161     PCC_140 PCC_161     PCC_36 PCC_161
# PCC_139    <NA> PCC_161     PCC_140 PCC_161     PCC_36 PCC_161
# 
#               11142022 run      #02032023 run             # 02102023
#   ID      dam     sire 
# PCC_101   PCC_140 PCC_161   PCC_101   PCC_140 PCC_161     PCC_101   PCC_140 PCC_161
# PCC_120   PCC_209  PCC_81   PCC_120   PCC_209  PCC_81     PCC_120   PCC_209  PCC_81
# PCC_127   PCC_140 PCC_161   PCC_127   PCC_140 PCC_161     PCC_127   PCC_140 PCC_161
# PCC_138   PCC_140 PCC_161   PCC_138   PCC_140 PCC_161     PCC_138   PCC_140 PCC_161 
# PCC_139   PCC_140 PCC_161   PCC_139   PCC_140 PCC_161     PCC_139   PCC_140 PCC_161
#
#     04052023                # 05252023                    # 02082024
# PCC_101 PCC_140 PCC_161     PCC_101  PCC_36 PCC_161       PCC_101 PCC_140 PCC_161
# PCC_120 PCC_209  PCC_81     PCC_120 PCC_209  PCC_81       PCC_120 PCC_209  PCC_81
# PCC_127 PCC_140 PCC_161     PCC_127 PCC_36 PCC_161        PCC_140 PCC_161
# PCC_138 PCC_140 PCC_161     PCC_138 PCC_36 PCC_161        PCC_140 PCC_161
# PCC_139 PCC_140 PCC_161     PCC_139  PCC_36 PCC_161       PCC_140 PCC_161

pdf(file = paste0("../pedigree_reconstruction/PCC_summary_seq_", date, ".pdf"), height = 6, width = 10)
SummarySeq(ped.PCC.ped)
dev.off()

pdf(file = paste0("../pedigree_reconstruction/PCC_dyad_", date, ".pdf"), height = 12, width = 12)
PlotRelPairs(GetRelM(ped.PCC.ped$Pedigree, GenBack = 2, patmat=TRUE))
dev.off()

pdf(file = paste0("../pedigree_reconstruction/PCC_dyad_simple_", date, ".pdf"), height = 12, width = 12)
PlotRelPairs(GetRelM(ped.PCC.ped$Pedigree, GenBack = 1, patmat=TRUE))
dev.off()


gen <- getGenerations(Ped = ped.PCC.par$Pedigree)
hist(gen)

PCC_snp_stats <- SnpStats(GenoM = PCC_prep_gt_all_drop2, Pedigree = ped.PCC.par$Pedigree)
hist(PCC_snp_stats$Err.hat)
abline(v = mean(PCC_snp_stats$Err.hat), col = "blue", lty = 2)
abline(v = error_rate, col = "green", lty = 2)
mean(PCC_snp_stats$Err.hat) # 0.02014067

### 7: Compare to previous pedigree  --------------------------------------------------------------------
# 01292024, compare 1K to 2K snps
#load(file = "../pedigree_reconstruction/PCC_pedigree_results_01292024_2KSNPs.Robj", verbose = T)
#PCC_results[[4]]$Pedigree

#pdf(file = paste0("../pedigree_reconstruction/PCC_PedCompare_", date, ".pdf"), height = 6, width = 10)
#PedCompare(Ped1 = PCC_results[[4]]$Pedigree, Ped2 = ped.PCC.ped$Pedigree)
#dev.off()

# 40502023
load(file = "../pedigree_reconstruction/PCC_pedigree_results_40502023.Robj", verbose = T)
#PCC_results[[4]]$Pedigree

# 05252023
load(file = "../pedigree_reconstruction/PCC_pedigree_results_05252023.Robj", verbose = T)

# 02082024
load(file = "../pedigree_reconstruction/PCC_pedigree_results_02082024.Robj", verbose = T)

# 05072024
load(file = "../pedigree_reconstruction/PCC_pedigree_results_05072024.Robj", verbose = T)

pdf(file = paste0("../pedigree_reconstruction/PCC_PedCompare_", date, ".pdf"), height = 6, width = 10)
PedCompare(Ped1 = PCC_results[[4]]$Pedigree, Ped2 = ped.PCC.ped$Pedigree)
dev.off()

### 8: Pedigree Confidence  --------------------------------------------------------------------

load(paste0("../pedigree_reconstruction/PCC_pedigree_results_", date, ".Robj"), verbose = T)

# using PCC_pedigree and PCC_meta from PCC_pedigree_results_06262024.Robj
PCC_pedigree <- PCC_results[[4]]$Pedigree
PCC_meta <- PCC_results[[2]]

PCC.ped.conf <- EstConf(Pedigree = PCC_pedigree, LifeHistData = PCC_meta)

# row, column, matrix
par(mfrow= c(1,2))
barplot(height = PCC.ped.conf[["PedErrors"]][,3,2], 
        main = "Sire mismatches", 
        ylab = "Proportion mismatch", 
        xlab = "individual + parent codes")

barplot(height = PCC.ped.conf[["PedErrors"]][,3,1], 
        main = "Dam mismatches", 
        ylab = "Proportion mismatch", 
        xlab = "individual + parent codes")


PCC.ped.conf[["ConfProb"]]

save(PCC.ped.conf, file = "../pedigree_reconstruction/pcc_conf_06262024.robj")

### 9: Traditional Pedigree  --------------------------------------------------------------------
# To get the output pedigree into kinship2 compatible format:

# if using previously run pedigree: 
ped.PCC.ped <- PCC_results[[4]]

PedP <- sequoia::PedPolish(ped.PCC.ped$Pedigree, DropNonSNPd=FALSE,
                           FillParents = TRUE)
PedP$Sex <- with(PedP, ifelse(id %in% dam, "female",  "male"))
# default to 'male' to avoid warning: "More than 25% of the gender values are
#  'unknown'"

Ped.fix <- with(PedP, kinship2::fixParents(id=id, dadid=sire, momid=dam,
                                           sex=Sex))
Ped.k <- with(Ped.fix, kinship2::pedigree(id, dadid, momid, sex, missid=0))

genotyped <- grep("PCC", Ped.k$id)
not_geno <- grep("PCC", Ped.k$id, invert=TRUE)
affected <- vector(length=length(Ped.k$id))
affected[genotyped] <- 1
affected[not_geno] <- 0

pdf(file = paste0("../pedigree_reconstruction/PCC_pedigree_", date, ".pdf"), height = 10, width = 10)
par(mar=c(5,4,7,2)+0.1)
plot(Ped.k, id = rep("", 323), symbolsize=1, 
     density = c(-1, 35, 65, 20), 
     affected = affected)
dev.off()

### 10: Save results  --------------------------------------------------------------------
PCC_results <- list(PCC_gt_for_ped_drop, PCC_meta, PCC_dupe_inds, ped.PCC.ped, error_rate, date) #  PCC.ped.conf,
save(PCC_results, file = paste0("../pedigree_reconstruction/PCC_pedigree_results_", date, ".Robj"))


## ELF ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

# input data: 
# SNP: ELF_prep_gt_all # 780 2223 
# Meta: ELF_meta
# Descriptions from Sequoia (https://jiscah.github.io/articles/vignette_main/book/running-sequoia.html)

### 1: 'PRE' --------------------------------------------------------------------
# check that the genotype data, life history data and input parameters are in a valid format, create 'Specs' element of output list

ped.ELF.pre <- sequoia(GenoM = ELF_gt_for_ped, LifeHistData=ELF_meta, Module = "pre", Err = error_rate)

# 'BY.max' must be greater than or equal to 'BY.min', or NA

# 10182022
# Warning:  There are 7 monomorphic (fixed) SNPs, these will be excluded 
# After exclusion, There are  788  individuals and  2364  SNPs.
# Ageprior: Flat 0/1, overlapping generations, MaxAgeParent = 13,13

# 11142022
# Warning:  There are 9 monomorphic (fixed) SNPs, these will be excluded 
# After exclusion, There are  788  individuals and  2362  SNPs.
# Ageprior: Flat 0/1, overlapping generations, MaxAgeParent = 13,13

# 2032023
# Warning:  There are 213 monomorphic (fixed) SNPs, these will be excluded 
# After exclusion, There are  788  individuals and  2255  SNPs.
# Ageprior: Flat 0/1, overlapping generations, MaxAgeParent = 13,13

### 2: 'DUP' --------------------------------------------------------------------
# Check for (nearly) identical genotypes, and for duplicated IDs in the genotype and life history data

ped.ELF.dup <- sequoia(GenoM = ELF_gt_for_ped, LifeHistData=ELF_meta, Module = "dup", Err = error_rate)

#There were 2 likely duplicate genotypes found, consider removing

ped.ELF.dup$DupGenotype
# ELF_131    ELF_504 *** need to check in metadata
# ELF_544_P10 ELF_544_P5
ELF_meta[which(ELF_meta$ID == "ELF_131"),] # same sex, same by.max
ELF_meta[which(ELF_meta$ID == "ELF_504"),]

# ELF_504: 0A01701972
# ELF_131: 171

# ***I will throw one out, but need to double check with Eric!!! 

# sample individuals to throw out at random

sum(ELF_gt_for_ped[which(rownames(ELF_gt_for_ped) == "ELF_544_P10"),] == -9)
sum(ELF_gt_for_ped[which(rownames(ELF_gt_for_ped) == "ELF_544_P5"),] == -9)

toss_inds <- c("ELF_544_P10", "ELF_504") # keeping ELF_131 b/c more precise age info and ELF_544_P10 because less missing data! 

ELF_gt_for_ped_drop <- ELF_gt_for_ped[-match(toss_inds, rownames(ELF_gt_for_ped)),] # dim: 778 2223

#### no other inferred duplicates! 

# make sure names match genetic data 
sum(rownames(ELF_gt_for_ped_drop) %in% subset(pop_gen_stats, site == "ELF", select = id)$id) / length(rownames(ELF_gt_for_ped_drop))


### 3: 'PAR' --------------------------------------------------------------------
# Assign genotyped parents to genotyped individuals. Includes call to MakeAgePrior() to estimate AgePriors based on the just-assigned parents.
ped.ELF.par <- sequoia(GenoM = ELF_gt_for_ped_drop, LifeHistData=ELF_meta, Module = "par", Err = error_rate)

# 10182022 : assigned 515 dams and 335 sires to 786 individuals
# 11142022 : assigned 513 dams and 337 sires to 786 individuals
# 02032023 : assigned 528 dams and 339 sires to 786 individuals
# 02102023 : assigned 548 dams and 347 sires to 786 individuals
# 04052023 : assigned 507 dams and 328 sires to 786 individuals
# 05252023 : assigned 505 dams and 326 sires to 783 individuals
# 01292024 : assigned 500 dams and 324 sires to 778 individuals 1k snps (.2)
# 01292024 : assigned 500 dams and 324 sires to 778 individuals 1.5k snps (.1)
# 02082023 : assigned 499 dams and 323 sires to 777 individuals
# 04132024 : assigned 502 dams and 325 sires to 777 individuals
# 05082024 INFER : assigned 512 dams and 329 sires to 777 individuals

### 4: 'PED' --------------------------------------------------------------------
# Cluster half- and full-siblings and assign each cluster a dummy-parent; assign grandparents to sibships and singletons.

ped.ELF.ped <- sequoia(GenoM = ELF_gt_for_ped_drop, LifeHistData=ELF_meta, Module = "ped", Err = error_rate)

ped.ELF.ped$Pedigree

# 11142022 : assigned 728 dams and 718 sires to 786 + 87 individuals (real + dummy)
# 02032023 : assigned 732 dams and 711 sires to 786 + 88 individuals (real + dummy)
# 02102023 : assigned 737 dams and 730 sires to 786 + 89 individuals (real + dummy)
# 04052023 : assigned 790 dams and 779 sires to 786 + 111 individuals (real + dummy)
# 05252023 : assigned 763 dams and 756 sires to 783 + 100 individuals (real + dummy)
# 05262023 : assigned 762 dams and 754 sires to 782 + 104 individuals (real + dummy)
# 01292024 : assigned 746 dams and 752 sires to 778 + 100 individuals (real + dummy) 1k snps
# 01292024 : assigned 755 dams and 768 sires to 778 + 108 individuals (real + dummy) 1.5k snps
# 02082024 : assigned 747 dams and 747 sires to 777 + 107 individuals (real + dummy)
# 04132024 : assigned 753 dams and 755 sires to 777 + 107 individuals (real + dummy)
# 05082024 INFER : assigned 758 dams and 758 sires to 777 + 116 individuals (real + dummy)
# 06262024 INFER DRB: assigned 745 dams and 752 sires to 777 + 113 individuals (real + dummy)

ped.ELF.ped$Pedigree[grep("ELF_544", ped.ELF.ped$Pedigree$id),"id"] <- "ELF_544"
### 5: Rescue pairs  --------------------------------------------------------------------

# if working from previously generated pedigree: 
# load previous run 
# if(load_previous == TRUE){
#   load(paste0("../pedigree_reconstruction/ELF_pedigree_results_", date, ".Robj"), verbose = T)
#   ELF_gt_for_ped_drop <- ELF_results[[1]]
#   ELF_meta <- ELF_results[[2]]
#   ped.ELF.ped <- ELF_results[[4]]
#   error_rate <- PCC_results[[5]]
# }

ELF_maybe_rel_ped <- GetMaybeRel(GenoM = ELF_gt_for_ped_drop, LifeHistData=ELF_meta, Pedigree = ped.ELF.ped$Pedigree, Module = "ped", Err = error_rate)

ELF_maybe_p0 <- ELF_maybe_rel_ped$MaybeRel[ELF_maybe_rel_ped$MaybeRel$TopRel == "PO",]

ELF_exp_metaData[["ELF_181"]] # first cap 2011 as adult
ELF_exp_metaData[["ELF_859"]] # first cap 2016 as Y
ELF_exp_metaData[["ELF_841"]] # 2016 2 yo

subset(ELF_meta, ID == "ELF_181")
subset(ELF_meta, ID == "ELF_859")
subset(ELF_meta, ID == "ELF_841")

subset(ped.ELF.ped$Pedigree, id == "ELF_181")
subset(ped.ELF.ped$Pedigree, id == "ELF_859") # already assigned dam
subset(ped.ELF.ped$Pedigree, id == "ELF_841") # already assigned dam

subset(ELF_meta, ID == "ELF_524")

test <- CalcPairLL(Pairs = data.frame("ELF_181", "ELF_859"), GenoM = ELF_gt_for_ped_drop, Pedigree = NULL, LifeHistData = ELF_meta)

ELF_exp_metaData[Maybe_rel_ped[[1]][i,1]]
ELF_meta[which(ELF_meta$ID == Maybe_rel_ped[[1]][i,1]),]

ELF_exp_metaData[Maybe_rel_ped[[1]][i,2]]
ELF_meta[which(ELF_meta$ID == Maybe_rel_ped[[1]][i,2]),]

Maybe_rel_ped[[1]][i,]

# 15 : maybe can use size? 
# 14 : maybe can use size? 
# 13 : 12 OH, seems w/in margin of error? but 607 already has a mother, so 326 can't be dam
# 12 : maybe can use size? 
# 11: 12 OH, seems w/in margin of error? but 614 already has a father, so 130 can't be sire
# 10 : 12 OH, seems w/in margin of error? 
# 9 : seems like a good pair, ELF_144 is dam, ELF_346 is offspring
# 8 : 12 OH, seems w/in margin of error? 
# 7 : 12 OH, seems w/in margin of error? 
# 6 : 12 OH, seems w/in margin of error? 
# 5 : 12 OH, seems w/in margin of error? 
# 4 : 12 OH, seems w/in margin of error? 
# 3 : 12 OH, seems w/in margin of error? 
# 2 : probably can discern based on capture times... ELF_923 is offspring, ELF_347 is sire, already has sire! 
# 1 : 12 OH, seems w/in margin of error? 


# add in maybe pairs 13, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1

Maybe_rel_ped[[1]][i,]

### 9 
# which is offspring? 144: parent  436: offspring

# find offspring in pedigree
ped.ELF.ped$Pedigree[which(ped.ELF.ped$Pedigree$id == "ELF_346"),]

# add parent
ped.ELF.ped$Pedigree[which(ped.ELF.ped$Pedigree$id == "ELF_346"),"dam"] <- "ELF_144"
ped.ELF.ped$Pedigree[which(ped.ELF.ped$Pedigree$id == "ELF_346"),"LLRdam"] <- 20.05

### 6: Pedigree Visualization  --------------------------------------------------------------------
# input data: 
# SNP: ELF_prep_gt_all_drop
# Meta: ELF_meta

pdf(file = paste0("../pedigree_reconstruction/ELF_snp_stats_", date, ".pdf"), height = 6, width = 10)
out <- SnpStats(GenoM = ELF_gt_for_ped_drop, ErrFlavour = mean(error_rate))
dev.off()

pdf(file = paste0("../pedigree_reconstruction/ELF_age_prior_", date, ".pdf"), height = 6, width = 6)
MakeAgePrior(ped.ELF.ped$Pedigree, LifeHistData = ELF_meta)
dev.off()

pdf(file = paste0("../pedigree_reconstruction/ELF_summary_seq_", date, ".pdf"), height = 6, width = 10)
SummarySeq(ped.ELF.ped)
dev.off()

pdf(file = paste0("../pedigree_reconstruction/ELF_dyad_", date, ".pdf"), height = 12, width = 12)
PlotRelPairs(GetRelM(ped.ELF.ped$Pedigree, GenBack = 2, patmat=TRUE))
dev.off()

pdf(file = paste0("../pedigree_reconstruction/ELF_dyad_simple_", date, ".pdf"), height = 12, width = 12)
PlotRelPairs(GetRelM(ped.ELF.ped$Pedigree, GenBack = 1, patmat=TRUE))
dev.off()


gen <- getGenerations(Ped = ped.ELF.ped$Pedigree)
hist(gen)

### 7: Known relatives  --------------------------------------------------------------------
# creates list where the list name is the name of the captive mom, and the contents of the list entry are the years she was held for partuition
ELF_pedigree <- ped.ELF.ped$Pedigree
captive.inds <- list()
count = 1
for(i in 1:length(ELF_exp_metaData)){
  if(any(ELF_exp_metaData[[i]]$held_captive == "Y")){
    captive.inds[[count]] <- subset(ELF_exp_metaData[[i]], held_captive == "Y", select = year)$year
    names(captive.inds)[count] <- names(ELF_exp_metaData)[[i]]
    count = count + 1
  }
}

names(captive.inds)
captive.kids <- vector(mode = "list", length = length(captive.inds))
names(captive.kids) <- names(captive.inds)
# read in raw data to get mom info! 
ELF_rawData <- read.csv("~/Desktop/EMR_rapture/pedigree_reconstruction/ELF_LifeHistData.csv")[,1:45]

for(i in 1:length(captive.inds)){
  print(ELF_rawData[grepl(strsplit(names(captive.inds)[i], "_")[[1]][2], ELF_rawData$Other.Notes),"Other.Notes"])
  ELF_rawData[grepl(strsplit(names(captive.inds)[i], "_")[[1]][2], ELF_rawData$Other.Notes) & (grepl("mom", ELF_rawData$Other.Notes, ignore.case = TRUE) | grepl("mother", ELF_rawData$Other.Notes, ignore.case = TRUE)),"Other.Notes"]  
  captive.kids[[i]] <- ELF_rawData[grepl(strsplit(names(captive.inds)[i], "_")[[1]][2], ELF_rawData$Other.Notes) & (grepl("mom", ELF_rawData$Other.Notes, ignore.case = TRUE) | grepl("mother", ELF_rawData$Other.Notes, ignore.case = TRUE)),"ELF.ID"]
}

# filter didn't quite work for ELF_200
captive.kids[["ELF_200"]] <- captive.kids[["ELF_200"]][1:8]
# ELF_513 was held in lab for partuition in 2013, but passed 9-10 slugs, no neonates
captive.kids <- captive.kids[-which(names(captive.kids) == "ELF_513")]

# remove slug record
captive.kids[["ELF_161"]] <- captive.kids[["ELF_161"]][-6]

ELF_correct_assignment <- vector(mode = "list", length = length(captive.kids))
names(ELF_correct_assignment) <- names(captive.kids)
# what % do we correctly assign in pedigree? 
for(i in 1:length(captive.kids)){
  assigned_count = 0
  in_ped_count = 0
  for(j in 1:length(captive.kids[[i]])){
    kid <- paste0("ELF_", captive.kids[[i]][j])
    if(any(grepl(kid, ELF_pedigree$id))){ # if kid made it into the pedigree
      in_ped_count <- in_ped_count + 1
      if(subset(ELF_pedigree, id == kid, select = dam) == names(captive.kids)[i]){ # was it assigned the right mom?
        assigned_count = assigned_count + 1 
      }
    }
  }
  ELF_correct_assignment[[i]] <- assigned_count/in_ped_count
}

mean(unlist(ELF_correct_assignment) *100)

# 100% correct assignment! 

### 7: Compare to previous pedigree  --------------------------------------------------------------------
load(file = "../pedigree_reconstruction/ELF_pedigree_results_05072024.Robj", verbose = T)
ELF_results[[4]]$Pedigree

pdf(file = paste0("../pedigree_reconstruction/ELF_PedCompare_", date, ".pdf"), height = 6, width = 10)
PedCompare(Ped1 = ELF_results[[4]]$Pedigree, Ped2 = ped.ELF.ped$Pedigree)
dev.off()

### 8: Pedigree Confidence  --------------------------------------------------------------------

load(paste0("../pedigree_reconstruction/ELF_pedigree_results_", date, ".Robj"), verbose = T)

# using ELF_pedigree and ELF_meta from ELF_pedigree_results_06262024.Robj
ELF_pedigree <- ELF_results[[4]]$Pedigree
ELF_meta <- ELF_results[[2]]

ELF.ped.conf <- EstConf(Pedigree = ELF_pedigree, LifeHistData = ELF_results[[2]])

save(ELF.ped.conf, file = "../pedigree_reconstruction/elf_conf_06262024.robj")

# row, column, matrix
par(mfrow= c(1,2))
barplot(height = ELF.ped.conf[["PedErrors"]][,3,2], 
        main = "Sire mismatches", 
        ylab = "Proportion mismatch", 
        xlab = "individual + parent codes")

barplot(height = ELF.ped.conf[["PedErrors"]][,3,1], 
        main = "Dam mismatches", 
        ylab = "Proportion mismatch", 
        xlab = "individual + parent codes")


ELF.ped.conf[["ConfProb"]]
ELF.ped.conf[["PedErrors"]]
### 9: Traditional Pedigree  --------------------------------------------------------------------

# if using previous pedigree run: 
ped.ELF.ped <- ELF_results[[4]]

# To get the output pedigree into kinship2 compatible format:
PedP <- sequoia::PedPolish(ped.ELF.ped$Pedigree, DropNonSNPd=FALSE,
                           FillParents = TRUE)
PedP$Sex <- with(PedP, ifelse(id %in% dam, "female",  "male"))
# default to 'male' to avoid warning: "More than 25% of the gender values are
#  'unknown'"

Ped.fix <- with(PedP, kinship2::fixParents(id=id, dadid=sire, momid=dam,
                                           sex=Sex))
Ped.k <- with(Ped.fix, kinship2::pedigree(id, dadid, momid, sex, missid=0))

genotyped <- grep("ELF", Ped.k$id)
not_geno <- grep("ELF", Ped.k$id, invert=TRUE)
affected <- vector(length=length(Ped.k$id))
affected[genotyped] <- 1
affected[not_geno] <- 0
 
pdf(file = paste0("../pedigree_reconstruction/ELF_pedigree_", date, ".pdf"), height = 10, width = 10)
plot(Ped.k, id = rep("", 926), symbolsize=1, 
     density = c(-1, 35, 65, 20), 
     affected = affected)
dev.off()

### 10: Save results  --------------------------------------------------------------------
ELF_results <- list(ELF_gt_for_ped_drop, ELF_meta, ELF_dupes, ped.ELF.ped, date) # ped.conf
save(ELF_results, file = paste0("../pedigree_reconstruction/ELF_pedigree_results_", date, ".Robj"))

# ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------



## GRAVEYARD ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

## old code
## Remove known duplicates-- under construction
# this code just removes the known duplicates used to calculate the genotyping error... 
# need to change it so that it also removes duplicates that were detected by the software! 

PCC_infer_dupes <- ped.out.all$DupGenotype[,3:4] # inferred and known duplicates

PCC_dupe_inds # known duplicates

# add missing dupes
## PCC_66 PCC_44

# remove incorrect dupes: 
## PCC_283	PCC_4
## PCC_149	PCC_86
ped.out.all$DupGenotype[which(ped.out.all$DupGenotype[,4] == "PCC_4" | ped.out.all$DupGenotype[,3] == "PCC_4"),]
ped.out.all$DupGenotype[which(ped.out.all$DupGenotype[,4] == "PCC_149" | ped.out.all$DupGenotype[,3] == "PCC_149"),]

dupe_inds
toss_inds <- vector()
for(i in 1:nrow(PCC_dupe_inds)){
  ind3 <- PCC_dupe_inds[i,3]
  if(nchar(ind3) == 0){ # if there are two sets of the same ind
    toss <- sample(PCC_dupe_inds[i,1:2], 1)
    toss_inds <- c(toss_inds, toss)
  }else if(nchar(ind3) > 0){ # if there are three sets of the same ind
    toss <- sample(PCC_dupe_inds[i,1:3], 2)
    toss_inds <- c(toss_inds, toss)
  }
}

gt <- prep_gt_all # let's run with full set of SNPs for now, dim: 297 2438

gt <- PCC_prep_gt_all[-match(unlist(toss_inds), rownames(PCC_prep_gt_all)),] # dim: 273 2438

# ## Remove inferred duplicates

dupes <- ped.out.all$DupGenotype[,3:4] # inferred and known duplicates

dupe_inds # known duplicates

toss_inds <- vector()

for(i in 1:nrow(ped.out.all$DupGenotype)){
  # was it already filtered out with dupe_inds?
  ind1 <- ped.out.all$DupGenotype[i,3]
  ind2 <- ped.out.all$DupGenotype[i,4]
  
  if(sum(dupe_inds[,1:2] == ind1) == 0){
    if(sum(dupe_inds[1:2] == ind2) == 0){
      print(paste0(ind1, " and ", ind2, " need to be removed"))
      toss <- sample(ped.out.all$DupGenotype[i,3:4], 1)
      toss_inds <- c(toss_inds, toss)
    }
  }
}

# remove incorrect dupes: 
## PCC_283	PCC_4
## PCC_149	PCC_86

inferred_dupes <- unlist(toss_inds)[-na.omit(match(c("PCC_283", "PCC_4", "PCC_149", "PCC_86"), unlist(toss_inds)))]

gt <- gt[-na.omit(match(inferred_dupes, rownames(gt))),] # dim: 267 2438


