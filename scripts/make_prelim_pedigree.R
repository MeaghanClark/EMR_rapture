## File name: make_prelim_pedigree.R
## Purpose: Make pedigrees for EMR project
## M. I. Clark, June 2022
## Last updated: 10/17/2022

# this script was built for use with sequoia version 2.3.3. Updated script for sequoia 2.5.3 will make_pedigree.R

## Loading ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

## load libraries
library(sequoia)
library(kinship2)

## date we want to analyze
# date <- 11142022
date <- "02032023"

## load data
# PCC
load(paste0("../vcf_filtering/PCC_filtered_gt_all_snps_noZ", date, ".Robj"), verbose = T) # PCC_filt_gt_all
# load(paste0("../vcf_filtering/PCC_filtered_gt_low_snps_noZ", date, ".Robj"), verbose = T) # PCC_filt_gt_low
# load(paste0("../vcf_filtering/PCC_filtered_gt_high_snps_noZ", date, ".Robj"), verbose = T) # PCC_filt_gt_high

load("../pedigree_reconstruction/PCC_metaData_02102023.Robj", verbose = T) # PCC_LifeHistData

# ELF
#load(paste0("../vcf_filtering/ELF_filtered_gt_high_snps_noZ", date, ".Robj"), verbose = T) # ELF_filt_gt_high
# load(paste0("../vcf_filtering/ELF_filtered_gt_low_snps_noZ", date, ".Robj"), verbose = T) # ELF_filt_gt_low
load(paste0("../vcf_filtering/ELF_filtered_gt_all_snps_noZ", date, ".Robj"), verbose = T) # ELF_filt_gt_all

load("../pedigree_reconstruction/ELF_metaData_02102023.Robj", verbose = T) # ELF_LifeHistData

# duplicate individuals for PCC
PCC_dupe_inds <- read.csv("../pedigree_reconstruction/PCC_LifeHistData_dupes.csv", header=FALSE)
# only one known duplicate for ELF
ELF_inds <- colnames(ELF_filt_gt_all)
ELF_dupes <- ELF_inds[which(grepl("P", ELF_inds))]

dupe_inds <- rbind.data.frame(PCC_dupe_inds, c("ELF_544_P5", "ELF_544_P10", ""))
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
  meta[is.na(meta)] <- -9 # replace NAs with -9
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
      matches[real_dupes] <- unlist(match_pairs(pairs, real_dupes, ped.out))[1]
      mismatches[real_dupes] <- unlist(match_pairs(pairs, real_dupes, ped.out))[2]
  }
  return(cbind.data.frame(pairs, matches, mismatches))
}

match_pairs <- function(pairs, real_dupes, ped.out){
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
# PCC_prep_gt_low <- prep_gt_data(PCC_filt_gt_low)
# PCC_prep_gt_high <- prep_gt_data(PCC_filt_gt_high)
#ELF
ELF_prep_gt_all <- prep_gt_data(ELF_filt_gt_all)
# ELF_prep_gt_low <- prep_gt_data(ELF_filt_gt_low)
# ELF_prep_gt_high <- prep_gt_data(ELF_filt_gt_high)

# ------------------------------------------------------------------------------------------------------------

## Prepare meta data ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

#PCC
PCC_meta <- prep_metaData(PCC_LifeHistData, rownames(PCC_prep_gt_all))

#ELF
ELF_meta <- prep_metaData(ELF_LifeHistData, rownames(ELF_prep_gt_all))

# ------------------------------------------------------------------------------------------------------------

## Estimate genotyping error rate ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

# using PCC data because most duplicated individuals come from that

PCC_error_rate <- calc_gt_error_rate(PCC_dupe_inds, PCC_prep_gt_all)
mean(PCC_error_rate) # 0.003440915

# look for SNPS with high genotyping error rate

PCC_dupe_pairs <- PCC_dupe_inds
for (real_dupes in 1:nrow(PCC_dupe_pairs)){
  if(grepl("PCC_[:digit:]*", PCC_dupe_pairs[real_dupes,3])==TRUE){ # if this is a triplicate
    ind1 <- PCC_dupe_pairs[real_dupes,1]
    ind2 <- PCC_dupe_pairs[real_dupes,2]
    ind3 <- PCC_dupe_pairs[real_dupes,3]
    p1 <- c(ind1, ind2, "")
    p2 <- c(ind1, ind3, "")
    p3 <- c(ind2, ind3, "")
    PCC_dupe_pairs <- rbind.data.frame(PCC_dupe_pairs, p1, p2, p3)
  }
}
for (real_dupes in 1:nrow(PCC_dupe_pairs)){ # remove original triplicate row 
  if(grepl("PCC_[:digit:]*", PCC_dupe_pairs[real_dupes,3])==TRUE){
    PCC_dupe_pairs <- PCC_dupe_pairs[-real_dupes,]
  }
}

no.pairs <- vector(length = ncol(PCC_prep_gt_all))
no.matches <- vector(length = ncol(PCC_prep_gt_all))
# for each pair of individuals
for(j in 1:nrow(PCC_dupe_pairs)){
  dupes <- as.character(PCC_dupe_pairs[j,])
  for(k in 1:ncol(PCC_prep_gt_all)){
    if(PCC_prep_gt_all[which(rownames(PCC_prep_gt_all) == dupes[1]),k] != -9){
      if(PCC_prep_gt_all[which(rownames(PCC_prep_gt_all) == dupes[2]),k] != -9){
        no.pairs[k] =  no.pairs[k] + 1  # track number of pairs sequenced 
        if(PCC_prep_gt_all[which(rownames(PCC_prep_gt_all) == dupes[1]),k] == PCC_prep_gt_all[which(rownames(PCC_prep_gt_all) == dupes[2]),k]){
          no.matches[k] = no.matches[k] + 1
        }
      }
    }
  }
}

geno.success <- no.matches/no.pairs
hist(geno.success)
hist(no.matches)

PCC_prep_gt_all_ger_cutoff_95 <- PCC_prep_gt_all[,geno.success >= 0.95] # filter snps 292 2306
PCC_prep_gt_all_ger_cutoff_98 <- PCC_prep_gt_all[,geno.success >= 0.98] # filter snps 292 2245
PCC_prep_gt_all_ger_cutoff_90 <- PCC_prep_gt_all[,geno.success >= 0.90] # filter snps 
PCC_prep_gt_all_ger_cutoff_88 <- PCC_prep_gt_all[,geno.success >= 0.88] # filter snps 
PCC_prep_gt_all_ger_cutoff_92 <- PCC_prep_gt_all[,geno.success >= 0.92] # filter snps 292 2262

cut_offs <- seq(0.8, 0.99, by = 0.01)
#pdf(file = paste0("../pedigree_reconstruction/error_rate_", date, ".pdf"), height = 6, width = 10)
par(mfrow=c(1,3))
plot(NULL, xlim = range(cut_offs), ylim = c(2000, 2400), main = "No. SNPs", 
     xlab = "cut-off", ylab = "number of SNPs")
for(i in cut_offs){
  gt_cutoff <- PCC_prep_gt_all[,geno.success > i] 
  points(x = i , y = dim(gt_cutoff)[2])
}
# did that improve error rate? 
plot(NULL, xlim = range(cut_offs), ylim = c(0, 0.005), main = "Error Rate", 
     xlab = "cut-off", ylab = "genotyping error rate")
abline(h = 0, col = "red", lty = 2)
for(i in cut_offs){
  gt_cutoff <- PCC_prep_gt_all[,geno.success > i] 
  points(x = i , y = mean(calc_gt_error_rate(PCC_dupe_inds, gt_cutoff)))
}
abline(v = 0.92, col = "red", lty = 2)
# do we detect all dupes in pedigree? 
plot(NULL, xlim = range(cut_offs), ylim = c(0, 4), main = "dupes detected", 
     xlab = "cut-off", ylab = "no. undetected")
for(i in cut_offs){
  gt_cutoff <- PCC_prep_gt_all[,geno.success > i] 
  error_rate <- mean(calc_gt_error_rate(PCC_dupe_inds, gt_cutoff))
  ped.PCC.dup <- sequoia(GenoM = gt_cutoff, LifeHistData=PCC_meta, Module = "dup", Err = error_rate, quiet = TRUE)
  detected <- detect_dupe_pairs(PCC_dupe_inds, ped.PCC.dup) # all known duplicate pairs are detected except PCC_66 and PCC_44, but both are the same as PCC_179
  points(x = i, y = sum(detected[,4] == 0))
}
abline(v = 0.92, col = "green", lty = 2)
abline(v = 0.93, col = "red", lty = 2)
#dev.off()
## cut off of >= 0.91 should result in (almost) all duplicates detected and lowers genotyping error rate (to 0.003149121), but is likely higher than what I was 
## using before accounting for triplicate dupes

# > 0.95
PCC_error_rate_ger_95 <- calc_gt_error_rate(PCC_dupe_inds, PCC_prep_gt_all_ger_cutoff_95)
mean(PCC_error_rate_ger_95) # 0.001127076

# > 0.92 
PCC_error_rate_ger_92 <- calc_gt_error_rate(PCC_dupe_inds, PCC_prep_gt_all_ger_cutoff_92)
mean(PCC_error_rate_ger_92) # 0.002258857
hist(PCC_error_rate_ger_92)
abline(v = mean(PCC_error_rate_ger_92))

# > 0.98
PCC_error_rate_ger_98 <- calc_gt_error_rate(PCC_dupe_inds, PCC_prep_gt_all_ger_cutoff_98)
mean(PCC_error_rate_ger_98) # 0

# > 0.90 
PCC_error_rate_ger_90 <- calc_gt_error_rate(PCC_dupe_inds, PCC_prep_gt_all_ger_cutoff_90)
mean(PCC_error_rate_ger_90) # 0.002974235

# > 0.88
PCC_error_rate_ger_88 <- calc_gt_error_rate(PCC_dupe_inds, PCC_prep_gt_all_ger_cutoff_88)
mean(PCC_error_rate_ger_88) # 0.003006724


# running with >= 0.92 cut off
PCC_error_rate_ger <- PCC_error_rate_ger_92

# visualize cut off impacts
par(mfrow=c(1,2))
hist(PCC_error_rate)
hist(PCC_error_rate_ger_92)
plot(y = PCC_error_rate, x = 1:length(PCC_error_rate), xlab = "Duplicated individuals", ylab = "Genotyping error rate")
plot(y = PCC_error_rate_ger, x = 1:length(PCC_error_rate_ger), xlab = "Duplicated individuals", ylab = "Genotyping error rate")


# filter both PCC and ELF SNPS using GER 
PCC_prep_gt_all_ger_cutoff <- PCC_prep_gt_all[,geno.success >= 0.91] # filter snps 
ELF_prep_gt_all_ger_cutoff <- ELF_prep_gt_all[,geno.success >= 0.91]

# investigate outliers
cbind.data.frame(PCC_dupe_pairs, PCC_error_rate_ger_92)
outlier1 <- PCC_dupe_inds[8,]
outlier2 <- PCC_dupe_inds[14,]

calc_gt_error_rate(outlier1, PCC_prep_gt_all) # all, 0.01666667

# pairwise, all these are high 
calc_gt_error_rate(outlier1[c(1,2)], PCC_prep_gt_all) # 0.01666667
calc_gt_error_rate(outlier1[c(1,3)], PCC_prep_gt_all) # 0.01209503
calc_gt_error_rate(outlier1[c(2,3)], PCC_prep_gt_all) # 0.008225108

calc_gt_error_rate(outlier2, PCC_prep_gt_all) # 0.01622807

# how many differences are there between non-dupe inds?
background_div <- vector(length = 100)
for(i in 1:100){
  background_div[i] <- calc_gt_error_rate(t(as.data.frame(sample(rownames(PCC_prep_gt_all_ger_cutoff_92), 2))), PCC_prep_gt_all_ger_cutoff_92, tri = FALSE) # values around 0.3-0.4
}
hist(background_div) # 0.4074497
abline(v = mean(background_div))


# ELF duplicate
ELF_error_rate <- calc_gt_error_rate(dupe_inds[22,1:2], ELF_prep_gt_all, tri = FALSE) # 0.001234568
ELF_error_rate_ger <- calc_gt_error_rate(dupe_inds[22,1:2], ELF_prep_gt_all_ger_cutoff, tri = FALSE) # 0.001238646

error_rate <- c(PCC_error_rate_ger, ELF_error_rate_ger) 

mean(error_rate) # 0.001322022 with ger cut off SNPs, change from 0.002880964 with all SNPs, 0.002305068 with 0.90 cut off SNPs, 
                # 0.002907482 with 0.90 cut off SNPS when triplicate dupes are considered, 0.001196314 with 0.92 cut off SNPs when 
                # triplicate dupes are considered

# 0.90 cut off when triplicate dupes are considers catches all known duplicates in the data, but this does not improve genotyping error 
# rate

# mean error rate from filter_vcf.R is 0.004455166... shouldn't I use the most SNPs? 
# Is filtering based on gt error rate even valid? I don't think it is. Assumes that SNPs with high 
# error rates in duplicated individuals have high error over all. I think it's more conservative to 
# use the overall error rate... .4% is actually not that high! 
error_rate = 0.004136376


# ------------------------------------------------------------------------------------------------------------

## SNP set selection ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
## no need to run this code again, going with all snps

# what set of SNPs results in the most identified duplicate individuals
# out <- SnpStats(GenoM = PCC_prep_gt_all, ErrFlavour = mean(error_rate))
# out <- SnpStats(GenoM = PCC_prep_gt_low, ErrFlavour = mean(error_rate))
# out <- SnpStats(GenoM = PCC_prep_gt_high, ErrFlavour = mean(error_rate))

ped.out.all <- sequoia(GenoM = PCC_prep_gt_all, LifeHistData=PCC_meta, Module = "dup", Err = mean(error_rate)) # 30 dupes
ped.out.low <- sequoia(GenoM = PCC_prep_gt_low, LifeHistData=PCC_meta, Module = "dup", Err = mean(error_rate)) # 30 dupes
ped.out.high <- sequoia(GenoM = PCC_prep_gt_high, LifeHistData=PCC_meta, Module = "dup", Err = mean(error_rate)) # 32 dupes

ped.out.all$DupGenotype
ped.out.low$DupGenotype
ped.out.high$DupGenotype

# duplicates
two <- PCC_dupe_inds[-c(8, 17),1:2] 

detect_all <- detect_dupe_pairs(two, ped.out.all)
detect_low <- detect_dupe_pairs(two, ped.out.low)
detect_high <- detect_dupe_pairs(two, ped.out.high)

# thruples, will finish later! 

three <- dupe_inds[c(8, 17),]

# end: count of correctly identified pairs out of the 18 duplicate pairs, and 2 duplicate thruples in the dataset 

# what about thruple with highest genotyping error? PCC_66, PCC_44, PCC_179
ped.out.all$DupGenotype[which(ped.out.all$DupGenotype[,4] == "PCC_44" | ped.out.all$DupGenotype[,3] == "PCC_44"),]
# PCC_179 is the same as 44 and 66, but 44 and 66 are not marked independently, but they are all the same according to PIT tags! 

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

ped.PCC.pre <- sequoia(GenoM = PCC_prep_gt_all, LifeHistData=PCC_meta, Module = "pre", Err = error_rate)

# Warning:  There are 15 monomorphic (fixed) SNPs, these will be excluded 
# After exclusion, There are  292  individuals and  2359  SNPs.
# Ageprior: Flat 0/1, overlapping generations, MaxAgeParent = 10,10

### 2: 'DUP' --------------------------------------------------------------------
        # Check for (nearly) identical genotypes, and for duplicated IDs in the genotype and life history data

ped.PCC.dup <- sequoia(GenoM = PCC_prep_gt_all, LifeHistData=PCC_meta, SeqList = ped.PCC.pre, Module = "dup", Err = error_rate)

#### remove duplicates
detect_dupe_pairs(PCC_dupe_inds, ped.PCC.dup) # all known duplicate pairs are detected except PCC_66 and PCC_44, but both are the same as PCC_179

#### remove known duplicates
PCC_dupe_inds

# sample individuals to throw out at random
toss_inds <- vector()
for(i in 1:nrow(PCC_dupe_inds)){
  ind1 <- PCC_dupe_inds[i,1]
  ind2 <- PCC_dupe_inds[i,2]
  ind3 <- PCC_dupe_inds[i,3]
  if(nchar(ind3) == 0){ # if there are two sets of the same ind
    # check for missing metadata for first sample
    ind1_meta <- PCC_meta[which(PCC_meta[,1] == PCC_dupe_inds[i,1]),]
    # check for missing metadata for second sample
    ind2_meta <- PCC_meta[which(PCC_meta[,1] == PCC_dupe_inds[i,2]),]
    if(sum(ind1_meta == ind2_meta) == 4){ # if the metadata for both samples is the same... 
      toss <- as.character(sample(PCC_dupe_inds[i,1:2], 1)) # sample an individual
      toss_inds <- c(toss_inds, toss) # and add it to the list of individuals to remove
    }else{ # if the metadata for both samples is not the same... 
      if (sum(ind1_meta != -9) > sum(ind2_meta != -9)){ # if ind 1 has more complete metadata... 
        toss <- as.character(ind2_meta$ID) # toss ind 2
        toss_inds <- c(toss_inds, toss)
      }else if (sum(ind2_meta != -9) > sum(ind1_meta != -9)){ # if ind 2 has more complete metadata... 
        toss <- as.character(ind1_meta$ID) # toss ind 1
        toss_inds <- c(toss_inds, toss)
      } else if (ind1_meta$BY.max != ind2_meta$BY.max){ # if the samples have different BY.max
        if(ind1_meta$BY.max > ind2_meta$BY.max){ # if ind 1 has a higher BY.max than ind 2
          toss <- as.character(ind1_meta$ID) # toss ind 1, keep the more conservative BY.max
          toss_inds <- c(toss_inds, toss)
        }else if (ind2_meta$BY.max > ind1_meta$BY.max){ # if ind 2 has a higher BY.max than ind 1
          toss <- as.character(ind2_meta$ID) # toss ind 2, keep the more conservative BY.max
          toss_inds <- c(toss_inds, toss)
        }else{
          print(paste0("Attention to i = ", i))
          print(paste0(ind1_meta, " ", ind2_meta))
        }
      }
    }
  }else if(nchar(ind3) > 0){ # if there are three sets of the same ind
    # check for missing metadata for first sample
    ind1_meta <- PCC_meta[which(PCC_meta[,1] == PCC_dupe_inds[i,1]),]
    # check for missing metadata for second sample
    ind2_meta <- PCC_meta[which(PCC_meta[,1] == PCC_dupe_inds[i,2]),]
    # check for missing metadata for third sample
    ind3_meta <- PCC_meta[which(PCC_meta[,1] == PCC_dupe_inds[i,3]),]
    if(sum(ind1_meta == ind2_meta) == 4 && sum(ind1_meta == ind3_meta) == 4){ # if the metadata for all samples is the same... 
      toss <- as.character(sample(PCC_dupe_inds[i,1:3], 2)) # sample two individuals
      toss_inds <- c(toss_inds, toss) # and add it to the list of individuals to remove
    }else{ # if the metadata for both samples is not the same... 
          print(paste0("Attention to i = ", i))
          print(paste0(ind1_meta, " ", ind2_meta, " ", ind3_meta))
        }
      }
    }

PCC_prep_gt_all_drop <- PCC_prep_gt_all[-match(toss_inds, rownames(PCC_prep_gt_all)),] # dim: 269 2438

#### remove inferred duplicates

# what are the other pairs of dupes?
toss_inferred <- vector()
for(i in 1:nrow(ped.PCC.dup$DupGenotype)){
  if(sum(ped.PCC.dup$DupGenotype[i,c(3,4)] %in% as.vector(unlist(PCC_dupe_inds))) == 0){
    pair <- ped.PCC.dup$DupGenotype[i,c(3,4)]
    ind1 <- ped.PCC.dup$DupGenotype[i,3]
    ind2 <- ped.PCC.dup$DupGenotype[i,4]
      # check for missing metadata for first sample
      ind1_meta <- PCC_meta[which(PCC_meta[,1] == ind1),]
      # check for missing metadata for second sample
      ind2_meta <- PCC_meta[which(PCC_meta[,1] == ind2),]
      if(sum(ind1_meta == ind2_meta) == 4){ # if the metadata for both samples is the same... 
        toss <- as.character(sample(pair, 1)) # sample an individual
        toss_inferred <- c(toss_inferred, toss)
      }else{ # if the metadata for both samples is not the same... 
        if (sum(ind1_meta != -9) > sum(ind2_meta != -9)){ # if ind 1 has more complete metadata... 
          toss <- as.character(ind2_meta$ID) # toss ind 2
          toss_inferred <- c(toss_inferred, toss)
        }else if (sum(ind2_meta != -9) > sum(ind1_meta != -9)){ # if ind 2 has more complete metadata... 
          toss <- as.character(ind1_meta$ID) # toss ind 1
          toss_inferred <- c(toss_inferred, toss)
        } else if (ind1_meta$BY.max != ind2_meta$BY.max){ # if the samples have different BY.max
          if(ind1_meta$BY.max > ind2_meta$BY.max){ # if ind 1 has a higher BY.max than ind 2
            toss <- as.character(ind1_meta$ID) # toss ind 1, keep the more conservative BY.max
            toss_inferred <- c(toss_inferred, toss)
          }else if (ind2_meta$BY.max > ind1_meta$BY.max){ # if ind 2 has a higher BY.max than ind 1
            toss <- as.character(ind2_meta$ID) # toss ind 2, keep the more conservative BY.max
            toss_inferred <- c(toss_inferred, toss)
          }else{
            print(paste0("Attention to i = ", i))
            print(paste0(ind1_meta, " ", ind2_meta))
          }
        }
  }
  }
  print(paste0("For individuals ", pair[1], " and ", pair[2], ", ", toss, " is added to toss list. i is ", i))
}

# remove inferred duplicates that are actually different individuals from list
    # but I want to remove the individual with the least metadata

# need to check these in metadata! 
# SAME "For individuals PCC_315 and PCC_88, PCC_315 is added to toss list. i is 3" 
    # same PIT tag in metadata, confirmed same individual 
# DIFFERENT "For individuals PCC_210 and PCC_99, PCC_99 is added to toss list. i is 11"
    # 210 : 836852544, F
    # 99: 027327627, F
    # 1 SNP mismatch 
# SAME "For individuals PCC_202 and PCC_295, PCC_295 is added to toss list. i is 12"
    # 202: 836329312 M
    # 295: no PIT ID or sex, 
    # 6 mismatch, hard to tell but likely same individual/no evidence of different individuals
# SAME "For individuals PCC_298 and PCC_80, PCC_298 is added to toss list. i is 13" SAME
    # 298: 27320018
    # 80: 027320018
    # 1 mismatch, same PIT tag in metadata, confirmed same individual 
# SAME "For individuals PCC_132 and PCC_304, PCC_132 is added to toss list. i is 23"
    # 132: 027257271
    # 304: 27257271
    # 1 mismatch, same PIT tage in metadata, confirmed same individual 
# SAME "For individuals PCC_290 and PCC_314, PCC_314 is added to toss list. i is 33"
    # 290: 600784009, F
    # 314: 600784009 or 600789009 (hard to read on blood vial)
    # 31 mismatches likely same pit tag, likely same individual 

# 283/4: different PIT tags, 149/86: different pit tags 
incorrect_dupes <- c("PCC_283", "PCC_4", "PCC_149", "PCC_86", "PCC_210", "PCC_99")
toss_inferred <- toss_inferred[is.na(match(toss_inferred, incorrect_dupes))]
PCC_prep_gt_all_drop2 <- PCC_prep_gt_all_drop[-match(toss_inferred, rownames(PCC_prep_gt_all_drop)),] # 263 2468


### 3: 'PAR' --------------------------------------------------------------------
        # Assign genotyped parents to genotyped individuals. Includes call to MakeAgePrior() to estimate AgePriors based on the just-assigned parents.

ped.PCC.par <- sequoia(GenoM = PCC_prep_gt_all_drop2, LifeHistData=PCC_meta, Module = "par", Err = mean(error_rate))

# 10182022 : assigned 118 dams and 88 sires to 267 individuals
# 11142022 : assigned 110 dams and 87 sires to 264 individuals
# 02032023 : assigned 113 dams and 94 sires to 264 individuals
# 02102023 : assigned 121 dams and 89 sires to 264 individuals

### 4: 'PED' --------------------------------------------------------------------
        # Cluster half- and full-siblings and assign each cluster a dummy-parent; assign grandparents to sibships and singletons.

ped.PCC.par <- sequoia(GenoM = PCC_prep_gt_all_drop2, LifeHistData=PCC_meta, SeqList = ped.PCC.par, Module = "ped", Err = mean(error_rate))

ped.PCC.par$Pedigree

# 10182022 : assigned 158 dams and 150 sires to 267 + 27 individuals (real + dummy)
# 11142022 : assigned 156 dams and 155 sires to 264 + 28 individuals (real + dummy)
# 02032023: assigned 158 dams and 154 sires to 264 + 31 individuals (real + dummy)
# 02102023: assigned 162 dams and 158 sires to 264 + 32 individuals (real + dummy)

### 5: Pedigree Visualization  --------------------------------------------------------------------
  # input data: 
    # SNP: PCC_prep_gt_all_drop2
    # Meta: PCC_meta

pdf(file = paste0("../pedigree_reconstruction/PCC_snp_stats_", date, ".pdf"), height = 6, width = 10)
out <- SnpStats(GenoM = PCC_prep_gt_all_drop2, ErrFlavour = mean(error_rate))
dev.off()

pdf(file = paste0("../pedigree_reconstruction/PCC_age_prior_", date, ".pdf"), height = 6, width = 6)
PlotAgePrior(ped.PCC.par$AgePriors)
dev.off()

# inferred dupes w/ conflicting metadata: 

# "PCC_283" and "PCC_4", "PCC_149" and "PCC_86" , "PCC_210" and "PCC_99" )

# PCC_283 and PCC_4 (2 mismatches)
  # who are their assigned parents? 
  ped.PCC.par$Pedigree[which(ped.PCC.par$Pedigree$id == "PCC_283"),]
  ped.PCC.par$Pedigree[which(ped.PCC.par$Pedigree$id == "PCC_4"),]
  # are they assigned offspring? 
  ped.PCC.par$Pedigree[which(ped.PCC.par$Pedigree$dam == "PCC_283"),] # PCC_4 is assigned as offspring
  ped.PCC.par$Pedigree[which(ped.PCC.par$Pedigree$dam == "PCC_4"),] # PCC_223 is assigned as offspring
  
  # Assigned as mother-offspring pair 
  
  # Capture time
  # PCC_283: initial cap 2018, SVL 53.7, F, A
  # PCC_4: initial cap 2013, SVL 44.7, F, J
  
# PCC_149 and PCC_86 (4 mismatches)
  # who are their assigned parents? 
  ped.PCC.par$Pedigree[which(ped.PCC.par$Pedigree$id == "PCC_149"),] # PCC_200 is assigned as dam
  ped.PCC.par$Pedigree[which(ped.PCC.par$Pedigree$id == "PCC_86"),] # PCC_200 is assigned as dam
  # are they assigned offspring? 
  ped.PCC.par$Pedigree[which(ped.PCC.par$Pedigree$dam == "PCC_149"),] # no offpsring
  ped.PCC.par$Pedigree[which(ped.PCC.par$Pedigree$dam == "PCC_86"),] # no offspring
  
  # Assigned as mother-offspring pair
  
  # Capture time
  # PCC_149: initial cap 2015, SVL 41.2, M, J, 027298014
  # PCC_86: initial cap 2014, SVL 38.8, F, J, 027323569
  PCC_meta[which(PCC_meta$ID == "PCC_149"),] # BirthYear = 2015
  PCC_meta[which(PCC_meta$ID == "PCC_86"),] # BY.max = 2011
    # both were captured as juveniles, but have different age estimates
  
# PCC_210 and PCC_99 (1 mismatch)
  # who are their assigned parents? 
  ped.PCC.par$Pedigree[which(ped.PCC.par$Pedigree$id == "PCC_210"),] # PCC_76 assigned as sire
  ped.PCC.par$Pedigree[which(ped.PCC.par$Pedigree$id == "PCC_99"),] # PCC_76 assigned as sire, PCC_210 as dam
  # are they assigned offspring? 
  ped.PCC.par$Pedigree[which(ped.PCC.par$Pedigree$dam == "PCC_210"),] # PCC_99 is assigned as offspring
  ped.PCC.par$Pedigree[which(ped.PCC.par$Pedigree$dam == "PCC_99"),] # PCC_291 is assigned as offspring

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


pdf(file = paste0("../pedigree_reconstruction/PCC_summary_seq_", date, ".pdf"), height = 6, width = 10)
SummarySeq(ped.PCC.par)
dev.off()

pdf(file = paste0("../pedigree_reconstruction/PCC_dyad_", date, ".pdf"), height = 12, width = 12)
PlotRelPairs(GetRelM(ped.PCC.par$Pedigree, GenBack = 2, patmat=TRUE))
dev.off()

pdf(file = paste0("../pedigree_reconstruction/PCC_dyad_simple_", date, ".pdf"), height = 12, width = 12)
PlotRelPairs(GetRelM(ped.PCC.par$Pedigree, GenBack = 1, patmat=TRUE))
dev.off()


gen <- getGenerations(Ped = ped.PCC.par$Pedigree)
hist(gen)

PCC_snp_stats <- SnpStats(GenoM = PCC_prep_gt_all_drop2, Pedigree = ped.PCC.par$Pedigree)
hist(PCC_snp_stats$Err.hat)
abline(v = mean(PCC_snp_stats$Err.hat), col = "blue", lty = 2)
abline(v = error_rate, col = "green", lty = 2)
mean(PCC_snp_stats$Err.hat) # 0.01832852

### 6: Pedigree Confidence  --------------------------------------------------------------------

PCC.ped.conf <- EstConf(Pedigree = ped.PCC.par$Pedigree, LifeHistData = PCC_meta)

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

### 7: Traditional Pedigree  --------------------------------------------------------------------
# To get the output pedigree into kinship2 compatible format:
PedP <- sequoia::PedPolish(ped.PCC.par$Pedigree, DropNonSNPd=FALSE,
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
plot(Ped.k, id = rep("", 330), symbolsize=1, 
     density = c(-1, 35, 65, 20), 
     affected = affected)
dev.off()

### 8: Save results  --------------------------------------------------------------------
PCC_results <- list(PCC_prep_gt_all_drop2, PCC_meta, PCC_dupe_inds, ped.PCC.par, PCC.ped.conf, error_rate, date)
save(PCC_results, file = paste0("../pedigree_reconstruction/PCC_pedigree_results_", date, ".Robj"))


## ELF ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

# input data: 
# SNP: ELF_prep_gt_all # 788 2371 02032023: ELF_prep_gt_all_ger_cutoff
# Meta: ELF_meta
# Descriptions from Sequoia (https://jiscah.github.io/articles/vignette_main/book/running-sequoia.html)

### 1: 'PRE' --------------------------------------------------------------------
# check that the genotype data, life history data and input parameters are in a valid format, create 'Specs' element of output list

ped.ELF.pre <- sequoia(GenoM = ELF_prep_gt_all, LifeHistData=ELF_meta, Module = "pre", Err = mean(error_rate))

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

ped.ELF.dup <- sequoia(GenoM = ELF_prep_gt_all, LifeHistData=ELF_meta, SeqList = ped.ELF.pre, Module = "dup", Err = mean(error_rate))

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

toss_inds <- c(sample(c("ELF_544_P10", "ELF_544_P5"), 1), "ELF_504") # keeping ELF_131 b/c more precise age info

ELF_prep_gt_all_drop <- ELF_prep_gt_all[-match(toss_inds, rownames(ELF_prep_gt_all)),] # dim: 786 2371

#### no other inferred duplicates! 

### 3: 'PAR' --------------------------------------------------------------------
# Assign genotyped parents to genotyped individuals. Includes call to MakeAgePrior() to estimate AgePriors based on the just-assigned parents.

ped.ELF.par <- sequoia(GenoM = ELF_prep_gt_all_drop, LifeHistData=ELF_meta, SeqList = ped.ELF.dup, Module = "par", Err = mean(error_rate))

# 10182022 : assigned 515 dams and 335 sires to 786 individuals
# 11142022 : assigned 513 dams and 337 sires to 786 individuals
# 02032023 : assigned 528 dams and 339 sires to 786 individuals
# 02102023 : assigned 548 dams and 347 sires to 786 individuals

### 4: 'PED' --------------------------------------------------------------------
# Cluster half- and full-siblings and assign each cluster a dummy-parent; assign grandparents to sibships and singletons.

ped.ELF.par <- sequoia(GenoM = ELF_prep_gt_all_drop, LifeHistData=ELF_meta, SeqList = ped.ELF.par, Module = "ped", Err = mean(error_rate))

ped.ELF.par$Pedigree

# 11142022 : assigned 728 dams and 718 sires to 786 + 87 individuals (real + dummy)
# 02032023 : assigned 732 dams and 711 sires to 786 + 88 individuals (real + dummy)
# 02102023 : assigned 737 dams and 730 sires to 786 + 89 individuals (real + dummy)

### 5: Pedigree Visualization  --------------------------------------------------------------------
# input data: 
# SNP: ELF_prep_gt_all_drop
# Meta: ELF_meta

pdf(file = paste0("../pedigree_reconstruction/ELF_snp_stats_", date, ".pdf"), height = 6, width = 10)
out <- SnpStats(GenoM = ELF_prep_gt_all_drop, ErrFlavour = mean(error_rate))
dev.off()

pdf(file = paste0("../pedigree_reconstruction/ELF_age_prior_", date, ".pdf"), height = 6, width = 6)
PlotAgePrior(ped.ELF.par$AgePriors)
dev.off()

pdf(file = paste0("../pedigree_reconstruction/ELF_summary_seq_", date, ".pdf"), height = 6, width = 10)
SummarySeq(ped.ELF.par)
dev.off()

pdf(file = paste0("../pedigree_reconstruction/ELF_dyad_", date, ".pdf"), height = 12, width = 12)
PlotRelPairs(GetRelM(ped.ELF.par$Pedigree, GenBack = 2, patmat=TRUE))
dev.off()

pdf(file = paste0("../pedigree_reconstruction/ELF_dyad_simple_", date, ".pdf"), height = 12, width = 12)
PlotRelPairs(GetRelM(ped.ELF.par$Pedigree, GenBack = 1, patmat=TRUE))
dev.off()


gen <- getGenerations(Ped = ped.ELF.par$Pedigree)
hist(gen)

### 6: Pedigree Confidence  --------------------------------------------------------------------

ped.conf <- EstConf(Pedigree = ped.ELF.par$Pedigree, LifeHistData = ELF_meta)

# row, column, matrix
par(mfrow= c(1,2))
barplot(height = ped.conf[["PedErrors"]][,3,2], 
        main = "Sire mismatches", 
        ylab = "Proportion mismatch", 
        xlab = "individual + parent codes")

barplot(height = ped.conf[["PedErrors"]][,3,1], 
        main = "Dam mismatches", 
        ylab = "Proportion mismatch", 
        xlab = "individual + parent codes")


ped.conf[["ConfProb"]]

### 7: Traditional Pedigree  --------------------------------------------------------------------
# To get the output pedigree into kinship2 compatible format:
PedP <- sequoia::PedPolish(ped.ELF.par$Pedigree, DropNonSNPd=FALSE,
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
plot(Ped.k, id = rep("", 936), symbolsize=1, 
     density = c(-1, 35, 65, 20), 
     affected = affected)
dev.off()

### 8: Save results  --------------------------------------------------------------------
ELF_results <- list(ELF_prep_gt_all_drop, ELF_meta, ELF_dupes, ped.ELF.par, ped.conf, date)
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


