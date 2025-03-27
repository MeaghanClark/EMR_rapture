## File name: filter_vcf_popgen.R
## Purpose: Load and filter VCF file from STACKS for EMR project
## Based on code by R.H. Toczydlowski
## M. I. Clark, June 2022
## Last updated: 02/16/2022

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
library(strataG)
library(hierfstat)


## define color palette 
palette = "Archambault"
met.brewer(palette, n=6, type="continuous")
MetBrewer::colorblind.friendly(palette)
colors <- met.brewer(palette, n = 6)

# date <- format(Sys.Date(), "%m%d%Y")
# date <- "02032023"
# date = 040502023
# date = "05262023"
date = "02082024"
# ------------------------------------------------------------------------------------------------------------

## Load custom functions ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
source("./vcf_funcs.R")

# ------------------------------------------------------------------------------------------------------------

## Load genetic data  from vcf_filter.R ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

load(file = paste0("../vcf_filtering/PCC_filtered_gt_all_snps_noZ_", date, ".Robj"), verbose = T) # PCC_filt_gt_all
#colnames(PCC_filt_gt_all) <- str_remove(colnames(PCC_filt_gt_all), "_P[0-9]+")
load(file = paste0("../vcf_filtering/ELF_filtered_gt_all_snps_noZ_", date, ".Robj"), verbose = T) # ELF_filt_gt_all

# ------------------------------------------------------------------------------------------------------------

## Extract genotypes ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
## Extract genotypes
# genotypes <- vcfR::extract.gt(vcf, return.alleles = F)
# 
# # create list of loci 
# positions <- vcf@fix
# 
# # process genotypes
# #gt.unflt <- process_gt(genotypes)
# gt <- process_gt(genotypes)

# # fix column names
# inds <- unlist(lapply(colnames(gt), FUN = function(x){
#   y <- unlist(strsplit(x, split = "_"))
#   if(length(y) == 12){
#     z <- paste0(y[1], "_", y[2], "_", y[3])
#   }
#   if(length(y) == 8){
#     z <- paste0(y[1], "_", y[2])
#   }
#   return(z)
# }))
# 
# colnames(gt) <- inds

# colnames(vcf@gt) <- c("FORMAT", inds)

# ------------------------------------------------------------------------------------------------------------

## Drop individuals dropped in pedigree SNP filtering, including "high het" individuals ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

# gt <- subset(gt, select = -c(PCC_297, PCC_301, PCC_306, PCC_283, PCC_86, PCC_99, ELF_274, ELF_401, ELF_138, ELF_824))
# # 2223 1080
# 
# ## Remove individuals with low mean coverage
# gt <- subset(gt, select = -c(ELF_141, ELF_242, ELF_289, ELF_324, ELF_333, ELF_342, ELF_893, ELF_906, ELF_917, PCC_285, PCC_308, PCC_73))
# 
# dim(gt) # 5607 1068

# ------------------------------------------------------------------------------------------------------------

## Filter SNPs based on pedigree ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

# need to filter out individuals dropped in pedigree filtering
# and drop duplicate individuals 

load(paste0("../pedigree_reconstruction/ELF_pedigree_results_", "05072024", ".Robj"), verbose = T)
load(paste0("../pedigree_reconstruction/PCC_pedigree_results_", "05072024", ".Robj"), verbose = T)

# get list of individuals retained in pedigree
ELF_ped_inds <- ELF_results[[4]]$Pedigree$id
ELF_ped_inds <- ELF_ped_inds[grep("ELF", ELF_ped_inds)] # 777


PCC_ped_inds <- PCC_results[[4]]$Pedigree$id
PCC_ped_inds <- PCC_ped_inds[grep("PCC", PCC_ped_inds)] # 260

all_ped_inds <- c(ELF_ped_inds, PCC_ped_inds) # length: 1044

# filter genotypes to retain just pedigree individuals -- VCF
gt_flt <- gt[,match(all_ped_inds, inds)]

vcf_LD_flt <- vcf_LD
vcf_LD_flt@gt <- vcf_LD_flt@gt[,match(c("FORMAT", all_ped_inds), c("FORMAT", inds))]

inds <- colnames(gt_flt)

# filter genotypes to retain just pedigree individuals -- genotype matrices 

colnames(ELF_filt_gt_all) # just need to drop "ELF_544_P10" and "ELF_504"

ELF_filt_gt_ped <- ELF_filt_gt_all[,-which(colnames(ELF_filt_gt_all)== "ELF_544_P10")]
ELF_filt_gt_ped <- ELF_filt_gt_ped[,-which(colnames(ELF_filt_gt_ped) == "ELF_504")] # 2223  777, matches! 
colnames(ELF_filt_gt_ped) <- str_remove(colnames(ELF_filt_gt_ped), "_P[0-9]+")

colnames(PCC_filt_gt_all) # 288
length(PCC_ped_inds) # 260

dim(PCC_filt_gt_all[,colnames(PCC_filt_gt_all) %in% PCC_ped_inds]) # 2223  260

PCC_filt_gt_ped <- PCC_filt_gt_all[,colnames(PCC_filt_gt_all) %in% PCC_ped_inds]
colnames(PCC_filt_gt_ped) <- str_remove(colnames(PCC_filt_gt_ped), "_P[0-9]+")

save(PCC_filt_gt_ped, file = "../pedigree_reconstruction/PCC_filt_gt_ped.Robj")
save(ELF_filt_gt_ped, file = "../pedigree_reconstruction/ELF_filt_gt_ped.Robj")

# USE ELF_filt_gt_ped and PCC_filt_gt_ped MOVING FORWARD
load("../pedigree_reconstruction/PCC_filt_gt_ped.Robj")
load("../pedigree_reconstruction/ELF_filt_gt_ped.Robj")

# ------------------------------------------------------------------------------------------------------------

# maf viz-----------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#maf <- calc_maf(t(t_gt_flt))
#hist(maf)

ELF_maf <- calc_maf(ELF_filt_gt_ped)
PCC_maf <- calc_maf(PCC_filt_gt_ped)

hist(ELF_maf)
hist(PCC_maf)

ELF_filt_gt_ped_maf <- ELF_filt_gt_ped[which(ELF_maf > 0.05),]
PCC_filt_gt_ped_maf <- PCC_filt_gt_ped[which(ELF_maf > 0.05),]

# ------------------------------------------------------------------------------------------------------------


## Basic popgen stats ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
# USING ELF_filt_gt_ped and PCC_filt_gt_ped, 
# need gt matrix with N rows and L columns, where N is the number of individuals and L is the number of loci
  # i.e. rows = individuals, L = loci 

# create joint matrix 

t_gt_flt <- rbind(t(ELF_filt_gt_ped), t(PCC_filt_gt_ped))

calcThetaW(t_gt_flt) # Wu and Watterson's theta at polymorphic sites
calcThetaW(matrix(t_gt_flt[1,])) # Wu and Watterson's theta at polymorphic sites

theta <- apply(t_gt_flt, 1, FUN <- function(x) calcThetaW(matrix(x)))

# PWP
#### PCC 
PCC_freqs <- t_gt_flt[grep("PCC", rownames(t_gt_flt)),]/2

PCC_coGeno <- matrix(NA,nrow(PCC_freqs),nrow(PCC_freqs))
for (i in 1:nrow(PCC_freqs)){
  for (j in 1:nrow(PCC_freqs)){
    PCC_coGeno[i,j] <- sum(!is.na(PCC_freqs[i,]) & !is.na(PCC_freqs[j,]))
    PCC_coGeno[j,i] <- sum(!is.na(PCC_freqs[i,]) & !is.na(PCC_freqs[j,]))
  }
}

PCC_PWP <- freqs2pairwisePi(freqs = PCC_freqs, coGeno = PCC_coGeno)
colnames(PCC_PWP) <- rownames(PCC_freqs)
rownames(PCC_PWP) <- rownames(PCC_freqs)

# no missing data 
PCC_freqs_nMD <- PCC_freqs[, colSums(is.na(PCC_freqs)) == 0]
PCC_PWP_nMD <- freqs2pairwisePi(freqs = PCC_freqs_nMD)

plot(PCC_PWP, PCC_PWP_nMD)

#### ELF
ELF_freqs <- t_gt_flt[grep("ELF", rownames(t_gt_flt)),]/2

ELF_coGeno <- matrix(NA,nrow(ELF_freqs),nrow(ELF_freqs))
for (i in 1:nrow(ELF_freqs)){
  for (j in 1:nrow(ELF_freqs)){
    ELF_coGeno[i,j] <- sum(!is.na(ELF_freqs[i,]) & !is.na(ELF_freqs[j,]))
    ELF_coGeno[j,i] <- sum(!is.na(ELF_freqs[i,]) & !is.na(ELF_freqs[j,]))
  }
}

ELF_PWP <- freqs2pairwisePi(freqs = ELF_freqs, coGeno = ELF_coGeno)
colnames(ELF_PWP) <- rownames(ELF_freqs)
rownames(ELF_PWP) <- rownames(ELF_freqs)

# no missing data 
ELF_freqs_nMD <- ELF_freqs[, colSums(is.na(ELF_freqs)) == 0]
ELF_PWP_nMD <- freqs2pairwisePi(freqs = ELF_freqs_nMD)

plot(ELF_PWP, ELF_PWP_nMD)

### 
PWP <- list(PCC_PWP, ELF_PWP)

########## 
sites <- c(rep("ELF", sum(grepl("ELF", rownames(t_gt_flt)))), rep("PCC", sum(grepl("PCC", rownames(t_gt_flt)))))
hist(theta[grep("ELF", names(theta))])
hist(theta[grep("PCC", names(theta))], add = T)

boxplot(theta~sites, main = "individual theta at polymorphic sites", col = c("blue", "orange"))

het <- calcHet(t_gt_flt)

pdf(file =paste0("../vcf_filtering/site_het_", date, ".pdf"), width = 6, height = 6)
boxplot(het~sites, main = "individual het at polymorphic sites", col = colors[c(2,1)])
dev.off()


## PCA # ------------------------------------------------------------------------------------------------------------
# drop ELF_269
ELF_filt_gt_ped_maf <- ELF_filt_gt_ped_maf[,-which(colnames(ELF_filt_gt_ped_maf) == "ELF_269")]

t_gt_flt_maf <- rbind(t(ELF_filt_gt_ped_maf), t(PCC_filt_gt_ped_maf))

doPCA <- function(gt, nPCs=4)
{
  sampleCov <- stats::cov(t(gt)/2, use = "pairwise.complete.obs")
  pcAxes <- eigen(sampleCov)$vectors[, 2:(nPCs + 1)]
  row.names(pcAxes) <- row.names(gt)
  return(pcAxes)
}
getPCAVarExplained <- function(gt){
  sampleCov <- stats::cov(t(gt)/2, use = "pairwise.complete.obs")
  eig <- eigen(sampleCov)
  return(eig$values/sum(eig$values))
}

bb_pca <- doPCA(t_gt_flt_maf, nPCs = 20)
var_explained <- getPCAVarExplained(t_gt_flt_maf)
var_explained[1:20] *100

# apparent discrepancies between PCAs with this code and PCAs from genlight objects 
# as in: pca <- glPca(snps_genlight, nf = 10) # you can select 10 components

all_genlight <- new("genlight", gen = t_gt_flt_maf, ind.names = rownames(t_gt_flt_maf))
all_pca <- glPca(all_genlight, nf = 10)

par(mfrow=c(1,2))
plot(bb_pca[,1], bb_pca[,2], pch = 19, col = alpha(c(rep(colors[2], 777), rep(colors[1], 260)), alpha = 0.5), 
     xlab = "PC1", ylab = "PC2")
plot(x = all_pca$scores[,1], y = all_pca$scores[, 2], pch = 19,
     col = alpha(c(rep(colors[2], 777), rep(colors[1], 260)), alpha = 0.5), 
     xlab = "PC1", ylab = "PC2")
legend("bottomright", c("ELF", "PCC"), pch = 19, col = c("blue", "orange"))
# these are slightly different, but pretty similar

# plot for supplement
pdf(file =paste0("../vcf_filtering/PCA_both_maf_", date, ".pdf"), width = 6, height = 6)
plot(bb_pca[,1], bb_pca[,2], pch = 19, col = alpha(c(rep(colors[2], 777), rep(colors[1], 260)), alpha = 0.5), 
     xlab = "PC1", ylab = "PC2")
#points(bb_pca[136,1], bb_pca[136,2], col = "red")
legend("bottomright", legend = c("Cass", "Barry"), pch = 19, col = c(colors[2], colors[1]))
#points(bb_pca[136,1], bb_pca[136,2], col ="red")
dev.off()

plot(bb_pca[,7], bb_pca[,8], pch = 19, col = alpha(c(rep(colors[2], 776), rep(colors[1], 260)), alpha = 0.5), 
     xlab = "PC1", ylab = "PC2")
#points(bb_pca[136,1], bb_pca[136,2], col = "red")
legend("bottomright", legend = c("Cass", "Barry"), pch = 19, col = c(colors[2], colors[1]))

hist(bb_pca[grepl("ELF", rownames(bb_pca)),5], col = colors[2], main = NULL, xlim = c(-0.1,0.1))
hist(bb_pca[grepl("PCC", rownames(bb_pca)),5], col = colors[1], main = NULL, add = T)

# PCs that show kinship structure: 
# 2
# 3
# 5
# 6
joint_PCA_outlier <- list(bb_pca, var_explained)
save(joint_PCA_outlier, file = "../vcf_filtering/joint_pca_w_outlier.RObj")



ELF_pca <- doPCA(t(ELF_filt_gt_ped_maf), nPCs = 4)
ELF_var_explained <- getPCAVarExplained(ELF_filt_gt_ped_maf)
ELF_var_explained[1:20] *100

# compare methods

ELF_genlight <- new("genlight", gen = t(ELF_filt_gt_ped_maf), ind.names = colnames(ELF_filt_gt_ped_maf))
ELF_genlight_PCA <- glPca(ELF_genlight, nf = 10)

plot(ELF_pca[,1], ELF_pca[,2], pch = 19, col = alpha(colors[2], alpha = 0.5), 
     xlab = "PC1", ylab = "PC2")

plot(ELF_genlight_PCA$scores[,1], ELF_genlight_PCA$scores[,2], pch = 19, col = alpha(colors[2], alpha = 0.5), 
     xlab = "PC1", ylab = "PC2")

pdf(file =paste0("../vcf_filtering/PCA_ELF_maf", date, ".pdf"), width = 6, height = 6)
plot(ELF_pca[,1], ELF_pca[,2], pch = 19, col = alpha(colors[2], alpha = 0.5), 
     xlab = "PC1", ylab = "PC2")
#points(ELF_pca[136,1], ELF_pca[136,2], col = "red")

dev.off()

PCC_pca <- doPCA(t(PCC_filt_gt_ped_maf), nPCs = 4)
PCC_var_explained <- getPCAVarExplained(PCC_filt_gt_ped_maf)
PCC_var_explained[1:20] *100

# compare methods
PCC_genlight <- new("genlight", gen = t(PCC_filt_gt_ped_maf), ind.names = colnames(PCC_filt_gt_ped_maf))
PCC_genlight_PCA <- glPca(PCC_genlight, nf = 10)


plot(PCC_pca[,1], PCC_pca[,2], pch = 19, col = alpha(colors[2], alpha = 0.5), 
     xlab = "PC1", ylab = "PC2")

plot(PCC_genlight_PCA$scores[,1], PCC_genlight_PCA$scores[,2], pch = 19, col = alpha(colors[2], alpha = 0.5), 
     xlab = "PC1", ylab = "PC2")

pca <- list(ELF_pca, PCC_pca, ELF_genlight_PCA$scores, PCC_genlight_PCA$scores)
save(pca, file = "../vcf_filtering/PCA_loadings_05212024.RObj")
# names(ELF.off[[1]]) %in% names(het)
# match(names(ELF.off[[1]]), names(het))
# het[match(names(ELF.off[[1]]), names(het))]

## prelim inbreeding modeling 

# pdf(file =paste0("../vcf_filtering/het_v_off_talk", date, ".pdf"), width = 6, height = 6)
# par(mfrow = c(2,2))
# plot(x = het[match(names(ELF.off[[1]]), names(het))], y = as.vector(ELF.off[[1]]), 
#      ylab = "number of offspring", 
#      xlab = "heterozygosity", 
#      main = "Cass dams", 
#      pch = 19, 
#      col = alpha(colors[2], alpha = 0.7))
# abline(lm(ELF.off[[1]]~het[match(names(ELF.off[[1]]), names(het))]), col = "red", lty = 2)
# 
# # summary(lm(ELF.off[[1]]~het[match(names(ELF.off[[1]]), names(het))]))
# 
# plot(x = het[match(names(ELF.off[[2]]), names(het))], y = as.vector(ELF.off[[2]]), 
#      ylab = "number of offspring", 
#      xlab = "heterozygosity", 
#      main = "Cass sires", 
#      pch = 19, 
#      col = alpha(colors[2], alpha = 0.7))
# abline(lm(ELF.off[[2]]~het[match(names(ELF.off[[2]]), names(het))]), col = "red", lty = 2)
# #summary(lm(ELF.off[[2]]~het[match(names(ELF.off[[2]]), names(het))]))
# 
# plot(x = het[match(names(PCC.off[[1]]), names(het))], y = as.vector(PCC.off[[1]]), 
#      ylab = "number of offspring", 
#      xlab = "heterozygosity", 
#      main = "Barry dams", 
#      pch = 19, 
#      col = alpha(colors[1], alpha = 0.7))
# abline(lm(PCC.off[[1]]~het[match(names(PCC.off[[1]]), names(het))]), col = "red", lty = 2)
# #summary(lm(PCC.off[[1]]~het[match(names(PCC.off[[1]]), names(het))]))

# plot(x = het[match(names(PCC.off[[2]]), names(het))], y = as.vector(PCC.off[[2]]), 
#      ylab = "number of offspring", 
#      xlab = "heterozygosity", 
#      main = "Barry sires", 
#      pch = 19, 
#      col = alpha(colors[1], alpha = 0.7))
# abline(lm(PCC.off[[2]]~het[match(names(PCC.off[[2]]), names(het))]), col = "red", lty = 2)
# #summary(lm(PCC.off[[2]]~het[match(names(PCC.off[[2]]), names(het))]))
# dev.off()
# 
# 
# # high het removed
# which(het[match(names(ELF.off[[2]]), names(het))] > 0.5) # 22
# names(ELF.off[[2]])[22]
# het[match(names(ELF.off[[2]]), names(het))][-22]
# plot(x = het[match(names(ELF.off[[2]]), names(het))][-22], y = as.vector(ELF.off[[2]])[-22], 
#      ylab = "number of offspring", 
#      xlab = "heterozygosity", 
#      main = "Cass sires", 
#      pch = 19, 
#      col = alpha(colors[2], alpha = 0.7))
# abline(lm(ELF.off[[2]][-22]~het[match(names(ELF.off[[2]]), names(het))][-22]), col = "red", lty = 2)

# ------------------------------------------------------------------------------------------------------------

# FST
all_genind_nmaf <- df2genind(X = t_gt_flt, ind.names = rownames(t_gt_flt), ncode = 1)

# add population information
all_genind_nmaf$pop <- as.factor(unlist(lapply(X = strsplit(rownames(t_gt_flt), "_"), FUN= function(x){x[1]})))

wc(all_genind_nmaf)

# Fgrm ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
# F can also be measured using the diagonal elements of a genomic relatedness matrix (Fgrm) (Bérénos et al., 2016; Huisman et al., 2016; Powell, Visscher, & Goddard, 2010; Pryce, Haile- Mariam, Goddard, & Hayes, 2014; Yang et al., 2010), which can be calculated with unmapped loc
# calculate in GCTA
# genotypes: t_gt_flt, 1050 individuals, 2307 loci 

# separate genotypes into PCC and ELF 
t_gt_PCC <- t_gt_flt[which(grepl("PCC_", rownames(t_gt_flt))),] # 214 2307
t_gt_ELF <- t_gt_flt[which(grepl("ELF_", rownames(t_gt_flt))),] # 782 2307

# calculate population-wide allele frequency of reference alleles, make vector, p, length = i 

get_pop_AF <- function(t_gt){
  # t_gt is a genotype matrix where each row is an individual and each column is a SNP. Genotypes are encoded as 0/1/2 or NA for missing data
  p <- vector(length = ncol(t_gt))
  for(i in 1:ncol(t_gt)){
    p[i] <- sum(t_gt[,i], na.rm = TRUE) / (2 * sum(!is.na(t_gt[,i])))
  }
  return(p)
}

calc_Fgrm <- function(t_gt){
  # t_gt is a genotype matrix where each row is an individual and each column is a SNP. Genotypes are encoded as 0/1/2 or NA for missing data
  
  # filter out monomorphic sites 
  t_gt_poly <- t_gt[,which(!colSums(t_gt, na.rm = TRUE)==0)]
  
  # calculate population-wide allele frequencies 
  p <- get_pop_AF(t_gt_poly)
  Fgrm <- vector(length = nrow(t_gt_poly))
  for(j in 1:nrow(t_gt_poly)){
    Fgrm_i <- vector(length = ncol(t_gt_poly))
    for(i in 1:ncol(t_gt_poly)){
      x <- t_gt_poly[j,i]
      Fgrm_i[i] <- (x^2 - (1 + 2*p[i])*x + 2*(p[i]^2)) / (2*p[i])*(1-p[i]) # sites have to be polymorphic? 
    } 
    Fgrm[j] <- sum(Fgrm_i, na.rm = TRUE) / sum(!is.na(Fgrm_i))
  }
  return(Fgrm)
}

Fgrm_ELF <- calc_Fgrm(t_gt_ELF)
Fgrm_PCC <- calc_Fgrm(t_gt_PCC)

# 
which(Fgrm_ELF > 0.8)

# t_gt is a genotype matrix where each row is an individual and each column is a SNP. Genotypes are encoded as 0/1/2 or NA for missing data

t_gt_migrant <- t_gt_ELF[which(Fgrm_ELF > 0.8),]

# filter out monomorphic sites 
t_gt_poly_ELF <- t_gt_ELF[,which(!colSums(t_gt_ELF, na.rm = TRUE)==0)]

# calculate population-wide allele frequencies 
p <- get_pop_AF(t_gt_ELF)
Fgrm <- vector(length = nrow(t_gt_migrant))

  Fgrm_i <- vector(length = length(t_gt_migrant))
  for(i in 1:length(t_gt_migrant)){
    x <- t_gt_migrant[i]
    Fgrm_i[i] <- (x^2 - (1 + 2*p[i])*x + 2*(p[i]^2)) / (2*p[i])*(1-p[i]) # sites have to be polymorphic? 
  } 

# test

plot(x = p, y = t_gt_migrant, pch = 19, col = alpha("black", alph = 0.05), xlab = "population  allele frequencies", ylab = "potential migrant genotypes")
plot(x = p, y = t_gt_ELF[40,], pch = 19, col = alpha("black", alph = 0.05), xlab = "population  allele frequencies", ylab = "genotypes")

# check to make sure it correlates with ML heterozygosity 

# het <- calcHet(t_gt_flt)

rownames(t_gt_flt)[1:10] == names(het)[1:10]

rownames(t_gt_ELF)[which(Fgrm_ELF > 0.5)]


# from pop_gen_stats
het <- pop_gen_stats$het
Fgrm <- pop_gen_stats$Fgrm
pdf(file = paste0("../vcf_filtering/Fgrm_vs_het_", date, ".pdf"), width = 10, height = 5)
# mat <- matrix(c(1,1,1,1,2,3,4,5), nrow = 4, ncol = 2, byrow = TRUE)
mat <- matrix(c(1,1,1,1,2,3,4,5), nrow = 2, ncol = 4, byrow = FALSE)

layout(mat)

plot(het~Fgrm, pch = 19, col = alpha("black", alpha = 0.25), xlab = "Fgrm", ylab = "heterozygosity")
legend("topright", c("both sites", "Barry County", "Cass County"), col = c("black", colors[1], colors[2]), pch = c(19, 15, 15))
hist(Fgrm[which(grepl("PCC", pop_gen_stats$site))], col = colors[1], xlab = "Fgrm", main = NULL)
hist(het[which(grepl("PCC", pop_gen_stats$site))], col = colors[1], xlab = "heterozygosity", main = NULL)
hist(Fgrm[which(grepl("ELF", pop_gen_stats$site))], col = colors[2], xlab = "Fgrm", main = NULL)
hist(het[which(grepl("ELF", pop_gen_stats$site))], col = colors[2], xlab = "heterozygosity", main = NULL)
dev.off()

summary(lm(het~c(Fgrm_ELF, Fgrm_PCC)))

range(Fgrm[which(grepl("PCC", pop_gen_stats$site))])
sum(Fgrm[which(grepl("PCC", pop_gen_stats$site))] > 0.1) / length(Fgrm[which(grepl("PCC", pop_gen_stats$site))])
range(Fgrm[which(grepl("ELF", pop_gen_stats$site))])
sum(Fgrm[which(grepl("ELF", pop_gen_stats$site))] > 0.1) / length(Fgrm[which(grepl("ELF", pop_gen_stats$site))])
hist(Fgrm[which(grepl("ELF", pop_gen_stats$site))])
# ------------------------------------------------------------------------------------------------------------
# get Ne estimates

# install strataG for LDNe in R capabilities 
# devtools::install_github('ericarcher/sprex')
# 
# options(repos = c(
#   zkamvar = 'https://zkamvar.r-universe.dev',
#   CRAN = 'https://cloud.r-project.org'))
# 
# install.packages('strataG')

# load data with estimated birth years 
load("../inbreeding_models/data_for_analyses_07132024.Robj")

table(as.numeric(subset(data, site == "PCC")$estBirthYear)) # max: 2015, 59 inds
table(as.numeric(subset(data, site == "ELF")$estBirthYear)) # max: 2013

PCC_cohort <- subset(data, site == "PCC" & estBirthYear == "2015")$id
t_gt_PCC_cohort <- t_gt_PCC[which(rownames(t_gt_PCC) %in% PCC_cohort),]

# wrangle data into gtype object
PCC_genlight <- new("genlight", gen = t_gt_PCC_cohort, ind.names = rownames(t_gt_PCC_cohort))
PCC_gtypes <- genlight2gtypes(PCC_genlight)

# use: ldNe 
PCC_ldNe <- ldNe(PCC_gtypes, drop.missing = TRUE, num.cores = 3)
# PCC_ldNe_maf_0.02 <- ldNe(PCC_gtypes, maf.threshold = 0.02, drop.missing = TRUE, num.cores = 3)


# wrangle data into gtype object
ELF_cohort <- subset(data, site == "ELF" & estBirthYear == "2013")$id
t_gt_ELF_cohort <- t_gt_ELF[which(rownames(t_gt_ELF) %in% ELF_cohort),]

ELF_genlight <- new("genlight", gen = t_gt_ELF_cohort, ind.names = rownames(t_gt_ELF_cohort))
ELF_gtypes <- genlight2gtypes(ELF_genlight)

# use: ldNe 
ELF_ldNe <- ldNe(ELF_gtypes, drop.missing = TRUE, num.cores = 3)
ELF_ldNe_maf_0.02 <- ldNe(ELF_gtypes, maf.threshold = 0.02, drop.missing = TRUE, num.cores = 3)


# implement jackknife method to incorporate sampling uncertainty into 95% confidence intervals 


jackknife <- function(gtypes, n_samples){
  Ne_vec <- vector(length = n_samples)
  for(i in 1:n_samples){
    Ne_vec[i] <- ldNe(gtypes[-i], drop.missing = TRUE, num.cores = 3)$Ne
  }
  return(Ne_vec)
}

calcCI_LDNe <- function(LDNe_obj, vec, ci = 0.95){
  theta <- LDNe_obj$Ne
  n <- LDNe_obj$S
  partials <- vec
  pseudos <- (n*theta) - (n-1)*partials
  jack.est <- mean(pseudos)
  jack.se <- sqrt(var(pseudos)/n)
  alpha = 1-ci
  CI <- qt(alpha/2,n-1,lower.tail=FALSE)*jack.se
  jack.ci <- c(jack.est - CI, jack.est + CI)

  return(jack.ci)
}


PCC_Ne_jacknife <- jackknife(PCC_gtypes, n_samples = 59)
save(PCC_Ne_jacknife, file = "../vcf_filtering/PCC_cohort_Ne_jackknife.Robj")

calcCI(PCC_Ne_jacknife)



hist(PCC_Ne_jacknife, xlim = c(20,70))
abline(v = calcCI_LDNe(PCC_ldNe, PCC_Ne_jacknife), col = "red", lty = 2)
abline(v = t.test(PCC_Ne_jacknife)$conf.int, col = "orange", lty = 3)
abline(v = mean(PCC_Ne_jacknife), col = "blue")

ELF_Ne_jacknife <- jackknife(ELF_gtypes, n_samples = 125)
save(ELF_Ne_jacknife, file = "../vcf_filtering/ELF_cohort_Ne_jackknife.Robj")

calcCI_LDNe(ELF_ldNe, ELF_Ne_jacknife)
hist(ELF_Ne_jacknife, xlim = c(22, 35))
abline(v = calcCI(ELF_Ne_jacknife), col = "red", lty = 2)
abline(v = mean(ELF_Ne_jacknife), col = "blue")

t.test(PCC_Ne_jacknife)$conf.int
# # Calculate jackknife estimate and it's variance
# jk_median_length = np.mean(median_lengths)
# jk_var = (n-1)*np.var(median_lengths)
# 
# # Assuming normality, calculate lower and upper 95% confidence intervals
# jk_lower_ci = jk_median_length - 1.96*np.sqrt(jk_var)
# jk_upper_ci = jk_median_length + 1.96*np.sqrt(jk_var)


## Collate all popgen-----------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

pop_gen_stats <- cbind.data.frame(rownames(t_gt_flt), sites, theta, het, c(Fgrm_ELF, Fgrm_PCC))
colnames(pop_gen_stats) <- c("id", "site", "theta", "het", "Fgrm")
save(pop_gen_stats, file = paste0("../vcf_filtering/pop_gen_stats_", date, ".Robj"))

save(PWP, file = paste0("../vcf_filtering/pwp_", date, ".Robj"))
# ------------------------------------------------------------------------------------------------------------


## Compare inbreeding w/ and w/o maf-----------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------


Fgrm_ELF_maf <- calc_Fgrm(t(ELF_filt_gt_ped_maf))
Fgrm_PCC_maf <- calc_Fgrm(t(PCC_filt_gt_ped_maf))
Fgrm_maf <- c(Fgrm_ELF_maf, Fgrm_PCC_maf)

plot(Fgrm_ELF, Fgrm_ELF_maf)

het_maf <- calcHet(t_gt_flt_maf)

pdf(file = paste0("../vcf_filtering/Fgrm_vs_het_maf", date, ".pdf"), width = 10, height = 5)
# mat <- matrix(c(1,1,1,1,2,3,4,5), nrow = 4, ncol = 2, byrow = TRUE)
mat <- matrix(c(1,1,1,1,2,3,4,5), nrow = 2, ncol = 4, byrow = FALSE)

layout(mat)

plot(het_maf~Fgrm_maf, pch = 19, col = alpha("black", alpha = 0.25), xlab = "Fgrm", ylab = "heterozygosity")
legend("topright", c("both sites", "Barry County", "Cass County"), col = c("black", colors[1], colors[2]), pch = c(19, 15, 15))
hist(Fgrm_maf[which(grepl("PCC", pop_gen_stats$site))], col = colors[1], xlab = "Fgrm", main = NULL)
hist(het_maf[which(grepl("PCC", pop_gen_stats$site))], col = colors[1], xlab = "heterozygosity", main = NULL)
hist(Fgrm_maf[which(grepl("ELF", pop_gen_stats$site))], col = colors[2], xlab = "Fgrm", main = NULL)
hist(het_maf[which(grepl("ELF", pop_gen_stats$site))], col = colors[2], xlab = "heterozygosity", main = NULL)
dev.off()

plot(het_maf~het, pch = 19, col = alpha("black", alpha = 0.25), xlab = "het", ylab = "het_maf")
legend("topright", c("both sites", "Barry County", "Cass County"), col = c("black", colors[1], colors[2]), pch = c(19, 15, 15))

pdf(file = paste0("../vcf_filtering/Fgrm_het_maf_vs_no", date, ".pdf"), width = 10, height = 10)
par(mfrow=c(2,2))
plot(het_maf[which(grepl("PCC", pop_gen_stats$site))]~ het[which(grepl("PCC", pop_gen_stats$site))], col = alpha(colors[1], alpha = 0.5), pch = 19, xlab = "het, no maf", ylab = "het, maf")
lines(x = c(0, 1), y = c(0, 1), col = "gray", lty = 2)
plot(het_maf[which(grepl("ELF", pop_gen_stats$site))]~ het[which(grepl("ELF", pop_gen_stats$site))], col = alpha(colors[2], alpha = 0.5), pch = 19, xlab = "het, no maf", ylab = "het, maf")
lines(x = c(0, 1), y = c(0, 1), col = "gray", lty = 2)

plot(Fgrm_maf[which(grepl("PCC", pop_gen_stats$site))]~ Fgrm[which(grepl("PCC", pop_gen_stats$site))],col = alpha(colors[1], alpha = 0.5), pch = 19, xlab = "Fgrm, no maf", ylab = "Fgrm, maf")
lines(x = c(0, 1), y = c(0, 1), col = "gray", lty = 2)
plot(Fgrm_maf[which(grepl("ELF", pop_gen_stats$site))]~ Fgrm[which(grepl("ELF", pop_gen_stats$site))],col = alpha(colors[1], alpha = 0.5), pch = 19, xlab = "Fgrm, no maf", ylab = "Fgrm, maf")
lines(x = c(0, 1), y = c(0, 1), col = "gray", lty = 2)
legend("topright", c("Barry County", "Cass County"), col = c(colors[1], colors[2]), pch = 19)
dev.off()

# ELF_269

ELF_results[[4]]$Pedigree[which(ELF_results[[4]]$Pedigree$sire == "ELF_269"),]
# offspring: ELF_473, ELF_534, ELF_666, ELF_774, ELF_828, but all assigned with negative LLRsire... should check in updated pedigree

Fgrm[which(colnames(ELF_filt_gt_ped_maf) == "ELF_473")] # all have low inbreeding
Fgrm[which(colnames(ELF_filt_gt_ped_maf) == "ELF_534")]
Fgrm[which(colnames(ELF_filt_gt_ped_maf) == "ELF_666")]
Fgrm[which(colnames(ELF_filt_gt_ped_maf) == "ELF_774")]
Fgrm[which(colnames(ELF_filt_gt_ped_maf) == "ELF_828")]

het[which(colnames(ELF_filt_gt_ped_maf) == "ELF_473")] # all have decently high het
het[which(colnames(ELF_filt_gt_ped_maf) == "ELF_534")]
het[which(colnames(ELF_filt_gt_ped_maf) == "ELF_666")]
het[which(colnames(ELF_filt_gt_ped_maf) == "ELF_774")]
het[which(colnames(ELF_filt_gt_ped_maf) == "ELF_828")]

# ------------------------------------------------------------------------------------------------------------

# load("../vcf_filtering/pop_gen_stats.Robj", verbose = T)

# make hist of inbreeding
hist(pop_gen_stats$Fgrm[grepl("PCC", pop_gen_stats$id)])
hist(pop_gen_stats$Fgrm[grepl("ELF", pop_gen_stats$id)])

pdf(file = paste0("../pedigree_exploration/Fgrm_site_poster_", date, ".pdf"), height = 8, width = 6)
par(mfrow=c(2,1))
hist(pop_gen_stats$Fgrm[grepl("PCC", pop_gen_stats$id)], 
     xlab = "Fgrm", 
     main = paste0("Barry County"), col = colors[1], cex.lab = 2, cex.axis = 2, cex.main = 2)
hist(pop_gen_stats$Fgrm[grepl("ELF", pop_gen_stats$id)], 
     xlab = "Fgrm", 
     main = paste0("Cass County"), col = colors[2], cex.lab = 2, cex.axis = 2, cex.main = 2)
dev.off()


## Nonspatial conStruct? ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
# get coords
load("../pedigree_exploration/coords40502023.Robj", verbose = TRUE)

# rename individuals with technical duplicate names 
which(grepl("_P", rownames(t_gt_flt)))
rownames(t_gt_flt)[396] <- "ELF_544"
rownames(t_gt_flt)[989] <- "PCC_300"
rownames(t_gt_flt)[990] <- "PCC_302"
rownames(t_gt_flt)[1017] <- "PCC_62"

# get PCC coordinates
PCC_coords <- coords[[1]][match(unique(coords[[1]]$ID), coords[[1]]$ID),]

PCC_gt_flt <- t_gt_flt[na.omit(match(PCC_coords$ID, rownames(t_gt_flt))),] 
PCC_AF = PCC_gt_flt/2
PCC_conStruct <- conStruct(spatial = FALSE, K = 2, freqs = PCC_AF, geoDist = NULL, coords = as.matrix(PCC_coords[,2:3]), prefix = "PCC", n.chains = 1, n.iter = 1, make.figs= TRUE, save.files = TRUE)


##### 

# troubleshooting conStruct AF code 
PCC_AF_drop <- conStruct:::drop.invars(PCC_AF) # same dims
PCC_AF_drop <- conStruct:::drop.missing(PCC_AF) # same dims 
n.loci <- ncol(PCC_AF)
obsCov <- conStruct:::calc.covariance(PCC_AF)
is.na(obsCov)

# ------------------------------------------------------------------------------------------------------------

# GRAVEYARD
## LD  ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
LD <- read.table("../vcf_filtering/filtered_LD_allpairs_nolimit.ld", header = TRUE)

# calculate pairwise distance between SNPs
LD$dist <- abs(LD$BP_A - LD$BP_B)

# remove Z chromosome
LD_noZ <- LD[which(LD$CHR_A != "Scate-Z"),]
LD_noZ <- LD_noZ[which(LD_noZ$CHR_B != "Scate-Z"),]

hist(LD_noZ$R2, main = "Histogram of r^2", col = "tomato3")
hist(LD_noZ$dist, main = "Histogram of distance", col = "tomato3")

hist(LD_noZ$dist[which(LD_noZ$dist < 1e5)])

# randomly sample snps to reduce size
dim(LD_noZ) # 31422628        8
keep <- sample(x = 1:31422628, size = 1e5, replace = FALSE)

LD_samp <- LD_noZ[keep,] # 
plot(x = LD_samp$dist, y = LD_samp$R2, col = alpha("tomato3", alpha = 0.1), pch = 19, 
     xlab = "genomic distance (bp)", ylab = "r^2", 
     xlim = c(0, 1e5))

lw1 <- loess(LD_samp$R2[which(LD_samp$dist < 1e5)] ~ LD_samp$dist[which(LD_samp$dist < 1e5)])
j <- order(LD_samp$dist[which(LD_samp$dist < 1e5)])
lines(LD_samp$dist[which(LD_samp$dist < 1e5)][j],lw1$fitted[j],col="darkolivegreen3",lwd=3)


plot(x = LD_noZ$dist, y = LD_noZ$R2, col = alpha("tomato3", alpha = 0.1), pch = 19, 
     xlab = "genomic distance (bp)", ylab = "r^2")
abline(h = mean(LD_noZ[which(LD_noZ$dist > 1e4),"R2"]), col = "darkolivegreen3", lwd = 4, lty = 2)

plot(x = LD_noZ$dist, y = LD_noZ$R2, col = alpha("tomato3", alpha = 0.2), pch = 19, 
     xlab = "genomic distance (bp)", ylab = "r^2", 
     xlim = c(0, 100000))
plot(x = LD_noZ$dist, y = LD_noZ$R2, col = alpha("tomato3", alpha = 0.2), pch = 19, 
     xlab = "genomic distance (bp)", ylab = "r^2", 
     xlim = c(0, 1000))

plot(x = LD_noZ$dist, y = LD_noZ$R2, col = alpha("tomato3", alpha = 0.05), pch = 19, 
     xlab = "genomic distance (bp)", ylab = "r^2", 
     xlim = c(0, 1e4))
abline(h = mean(LD_noZ[which(LD_noZ$dist > 1e4),"R2"]), col = "darkolivegreen3", lwd = 4, lty = 2)
lw1 <- loess(LD_noZ$R2 ~ LD_noZ$dist)
j <- order(LD_noZ$dist)
lines(LD_noZ$dist[j],lw1$fitted[j],col="thistle",lwd=3)


# 1285 snps 
# 3565 comparions 
LD_ma1 <- LD_noZ[which(LD_noZ$CHR_A == "Scate-ma1"),]
LD_ma1 <- LD_ma1[which(LD_ma1$CHR_B == "Scate-ma1"),]


hist(LD_ma1$R2, main = "Histogram of r^2", col = "tomato3")
hist(LD_ma1$dist, main = "Histogram of distance", col = "tomato3")

subset <- LD_noZ[which(LD_noZ$dist < 1000),]
plot(x = subset$dist, y = subset$R2, col = alpha("tomato3", alpha = 0.05), pch = 19, 
     xlab = "genomic distance (bp)", ylab = "r^2")
abline(v = 875, col = "thistle4", lwd = 4, lty = 2)

abline(h = mean(LD_noZ[,"R2"]), col = "thistle4", lwd = 4, lty = 2)
abline(h = median(LD_noZ[,"R2"]), col = "thistle1", lwd = 4, lty = 3)
legend("topright", c("mean", "median"), col = c("thistle4", "thistle1"), lwd = 4, lty = c(2, 3))


abline(h = mean(LD_noZ[which(LD_noZ$dist > 1e4),"R2"]), col = "darkolivegreen3", lwd = 4, lty = 2)
abline(h = mean(LD_noZ[which(LD_noZ$dist > 1e4),"R2"]) +sd(LD_noZ[which(LD_noZ$dist > 1e4),"R2"]) , col = "darkolivegreen4", lwd = 4, lty = 3)
abline(h = median(LD_noZ[which(LD_noZ$dist > 1e4),"R2"]), col = "darkolivegreen1", lwd = 4, lty = 3)


legend("topright", c("mean+sd, dist > 1e4", "mean, dist > 1e4", "median > 1e4", "mean", "median"), col = c("darkolivegreen4", "darkolivegreen3", "darkolivegreen1", "thistle4", "thistle1"), lwd = 4, lty = c(4, 2, 3, 2, 3))

mean(LD_noZ[,"R2"]) # 0.01344609
median(LD_noZ[,"R2"]) # 0.00464271

length(LD_noZ[which(LD_noZ$R2 < 0.01344609),"dist"])
# ------------------------------------------------------------------------------------------------------------

## Run PCA ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
# need to redo to work with genotype matrix! 
# code adapted from https://wurmlab.com/genomicscourse/2019/practicals/population_genetics/popgen

# convert vcfR format to genind
snps_genlight <- vcfR2genlight(x = vcf)

sites <- as.factor(c(rep("ELF", 793), rep("PCC", 297)))

snps_genlight@pop <- sites

# cool visualization
glPlot(snps_genlight)

# run PCA
pca <- glPca(snps_genlight, nf = 10) # you can select 10 components

# Quick plot
scatter(pca, posi = "none")

## plot eigenvalues (scree plot)
barplot(pca$eig[1:100], main="eigenvalues", col=heat.colors(length(pca$eig[1:100])))

# Plot of the first few axes coloured by population

pca$eig[1]/sum(pca$eig) *100
pca$eig[2]/sum(pca$eig) *100

pdf(file =paste0("../vcf_filtering/PCA_all_", date, ".pdf"), width = 6, height = 6)
plot(x = pca$scores[, 1], y = pca$scores[, 2],
     col = alpha(c(rep("blue", 793), rep("orange", 297)), 0.25), 
     cex = 1, pch = 19, 
     xlab = "PC1, 10.3% variation", ylab = "PC2, 2.4% variation")
legend("bottomright", c("ELF", "PCC"), pch = 19, col = c("blue", "orange"))
dev.off()


plot(pca$scores[,1], pca$scores[,3],
     col = c(rep("blue", 793), rep("orange", 297)), cex = 2)


# Separate pop PCA
sep_pop <- seppop(snps_genlight)

sep_pop$ELF

# run PCA for ELF
pca_ELF <- glPca(sep_pop$ELF, nf = 10) # you can select 10 components

# Quick plot
scatter(pca_ELF, posi = "none")

## plot eigenvalues (scree plot)
barplot(pca_ELF$eig[1:100], main="eigenvalues", col=heat.colors(length(pca_ELF$eig[1:100])))

# Plot of the first few axes coloured by population

pca_ELF$eig[1]/sum(pca_ELF$eig) *100
pca_ELF$eig[2]/sum(pca_ELF$eig) *100

pdf(file =paste0("../vcf_filtering/PCA_ELF_", date, ".pdf"), width = 6, height = 6)
plot(x = pca_ELF$scores[, 1], y = pca_ELF$scores[, 2],
     col = alpha(rep("blue", 793), 0.25), 
     cex = 1, pch = 19, 
     xlab = "PC1, 3.6% variation", ylab = "PC2, 3.1% variation")
dev.off()
plot(pca_ELF$scores[,1], pca_ELF$scores[,3],
     col = rep("blue", 793), cex = 2)

# run PCA for PCC
pca_PCC <- glPca(sep_pop$PCC, nf = 10) # you can select 10 components

# Quick plot
scatter(pca_PCC, posi = "none")

## plot eigenvalues (scree plot)
barplot(pca_PCC$eig[1:100], main="eigenvalues", col=heat.colors(length(pca_ELF$eig[1:100])))

# Plot of the first few axes coloured by population

pca_PCC$eig[1]/sum(pca_PCC$eig) *100
pca_PCC$eig[2]/sum(pca_PCC$eig) *100

pdf(file =paste0("../vcf_filtering/PCA_PCC_", date, ".pdf"), width = 6, height = 6)
plot(x = pca_PCC$scores[, 1], y = pca_PCC$scores[, 2],
     col = alpha(rep("orange", 793), 0.25), 
     cex = 1, pch = 19, 
     xlab = "PC1, 3.2% variation", ylab = "PC2, 3.0% variation")
dev.off()

plot(pca_PCC$scores[,1], pca_PCC$scores[,3],
     col = rep("blue", 793), cex = 2)
# ------------------------------------------------------------------------------------------------------------

## Run MAF PCA ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

# convert vcfR format to genind
snps_maf_genlight <- vcfR2genlight(x = vcf_maf)

sites <- as.factor(c(rep("ELF", 793), rep("PCC", 297)))

snps_maf_genlight@pop <- sites

# run PCA
pca_maf <- glPca(snps_maf_genlight, nf = 10) # you can select 10 components

# Quick plot
scatter(pca_maf, posi = "none")

## plot eigenvalues (scree plot)
barplot(pca_maf$eig[1:100], main="eigenvalues", col=heat.colors(length(pca_maf$eig[1:100])))

# Plot of the first few axes coloured by population

pca_maf$eig[1]/sum(pca_maf$eig) *100
pca_maf$eig[2]/sum(pca_maf$eig) *100

pdf(file =paste0("../vcf_filtering/PCA_maf_all_", date, ".pdf"), width = 6, height = 6)
plot(x = pca_maf$scores[, 1], y = pca_maf$scores[, 2],
     col = alpha(c(rep("blue", 793), rep("orange", 297)), 0.25), 
     cex = 1, pch = 19, 
     xlab = "PC1, 10.4% variation", ylab = "PC2, 2.4% variation")
legend("bottomright", c("ELF", "PCC"), pch = 19, col = c("blue", "orange"))
dev.off()
# ------------------------------------------------------------------------------------------------------------

## Run LD PCA ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

# convert vcfR format to genind
snps_ld_genlight <- vcfR2genlight(x = vcf_LD_flt)

snps_ld_genlight@ind.names <- all_ped_inds

#snps_ld_genlight <- snps_ld_genlight[snps_ld_genlight@ind.names %in% all_ped_inds]

sites <- as.factor(c(rep("ELF", sum(grepl("ELF", all_ped_inds))), rep("PCC", sum(grepl("PCC", all_ped_inds)))))
snps_ld_genlight@pop <- sites

# run PCA
pca_ld <- glPca(snps_ld_genlight, nf = 10) # you can select 10 components

# Quick plot
scatter(pca_ld, posi = "none")

## plot eigenvalues (scree plot)
barplot(pca_ld$eig[1:100], main="eigenvalues", col=heat.colors(length(pca_ld$eig[1:100])))

# Plot of the first few axes coloured by population

pca_ld$eig[1]/sum(pca_ld$eig) *100
pca_ld$eig[2]/sum(pca_ld$eig) *100

pdf(file =paste0("../vcf_filtering/PCA_maf_LD_flt_all_", date, ".pdf"), width = 6, height = 6)
plot(x = pca_ld$scores[, 1], y = pca_ld$scores[, 2],
     col = alpha(c(rep(colors[2], 786), rep(colors[1], 264)), 0.5), 
     cex = 1, pch = 19, cex.lab = 1.5,
     xlab = "PC1, 7.92% variation", ylab = "PC2, 2.4% variation")
legend("topright", c("Cass", "Barry"), pch = 19, col = c(colors[2], colors[1]), cex = 1.5)
dev.off()

# Separate pop PCA
sep_pop <- seppop(snps_ld_genlight)

sep_pop$ELF

# run PCA for ELF
pca_ELF <- glPca(sep_pop$ELF, nf = 10) # you can select 10 components

## plot eigenvalues (scree plot)
barplot(pca_ELF$eig[1:100], main="eigenvalues", col=heat.colors(length(pca_ELF$eig[1:100])))

# Plot of the first few axes coloured by population

pca_ELF$eig[1]/sum(pca_ELF$eig) *100
pca_ELF$eig[2]/sum(pca_ELF$eig) *100

pdf(file =paste0("../vcf_filtering/ELF_PCA_maf_LD_flt_all_", date, ".pdf"), width = 6, height = 6)
plot(x = pca_ELF$scores[, 1], y = pca_ELF$scores[, 2],
     col = alpha(rep(colors[2], 793), 0.5), 
     cex = 1, pch = 19, cex.lab = 1.5,
     xlab = "PC1, 3.5% variation", ylab = "PC2, 3.0% variation")
dev.off()

# run PCA for PCC
pca_PCC <- glPca(sep_pop$PCC, nf = 10) # you can select 10 components

## plot eigenvalues (scree plot)
barplot(pca_PCC$eig[1:100], main="eigenvalues", col=heat.colors(length(pca_ELF$eig[1:100])))

# Plot of the first few axes coloured by population

pca_PCC$eig[1]/sum(pca_PCC$eig) *100
pca_PCC$eig[2]/sum(pca_PCC$eig) *100

pdf(file =paste0("../vcf_filtering/PCC_PCA_maf_LD_flt_all_", date, ".pdf"), width = 6, height = 6)
plot(x = pca_PCC$scores[, 1], y = pca_PCC$scores[, 2],
     col = alpha(rep(colors[1], 793), 0.5), 
     cex = 1, pch = 19, cex.lab = 1.5,
     xlab = "PC1, 3.2% variation", ylab = "PC2, 3.0% variation")
dev.off()

PCA_scores <- list(pca_PCC$scores, pca_ELF$scores)
save(PCA_scores, file = paste0("../pedigree_exploration/PCA_scores_", date, ".Robj"))

# ------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------

## Compare PCAs ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

pdf(file =paste0("../vcf_filtering/PCA_comp_", date, ".pdf"), width = 6, height = 6)
par(mfrow = c(1,3))
# all 
plot(x = pca$scores[, 1], y = pca$scores[, 2],
     col = alpha(c(rep("blue", 793), rep("orange", 297)), 0.25), 
     cex = 1, pch = 19, 
     xlab = "PC1, 10.3% variation", ylab = "PC2, 2.4% variation", 
     main = "all SNPs, 10096")
legend("bottomright", c("ELF", "PCC"), pch = 19, col = c("blue", "orange"))

#MAF
plot(x = pca_maf$scores[, 1], y = pca_maf$scores[, 2],
     col = alpha(c(rep("blue", 793), rep("orange", 297)), 0.25), 
     cex = 1, pch = 19, 
     xlab = "PC1, 10.4% variation", ylab = "PC2, 2.4% variation", 
     main = "maf cut off, 8269")
legend("bottomright", c("ELF", "PCC"), pch = 19, col = c("blue", "orange"))


# LD
plot(x = pca_ld$scores[, 1], y = pca_ld$scores[, 2],
     col = alpha(c(rep("blue", 793), rep("orange", 297)), 0.25), 
     cex = 1, pch = 19, 
     xlab = "PC1, 8.26% variation", ylab = "PC2, 2.3% variation", 
     main = "maf and LD pruned, 2307")
legend("bottomright", c("ELF", "PCC"), pch = 19, col = c("blue", "orange"))

dev.off()

# ------------------------------------------------------------------------------------------------------------

## PCAs colored by birth year --------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
# isolate individual with known birth years
ELF_by <- as.numeric(ELF_meta$BirthYear)
PCC_by <- as.numeric(PCC_meta$BirthYear)

ELF_inds_w_by <- ELF_meta$ID[which(!ELF_by == -9)]
PCC_inds_w_by <- PCC_meta$ID[which(!PCC_by == -9)]

ELF_inds_w_by <- cbind.data.frame(ELF_inds_w_by, ELF_meta$BirthYear[which(!ELF_by == -9)])
PCC_inds_w_by <- cbind.data.frame(PCC_inds_w_by, PCC_meta$BirthYear[which(!PCC_by == -9)])

colnames(ELF_inds_w_by) <- c("id","by")
colnames(PCC_inds_w_by) <- c("id", "by")

ELF_inds_w_by$by <- as.numeric(ELF_inds_w_by$by)
PCC_inds_w_by$by <- as.numeric(PCC_inds_w_by$by)

pdf(file =paste0("../vcf_filtering/ELF_by_hist_", date, ".pdf"), width = 6, height = 6)
hist(ELF_inds_w_by$by, xlab = "birth year", col = colors[2], main = NULL, cex.lab = 2, cex.axis = 1.5)
dev.off()

pdf(file =paste0("../vcf_filtering/PCC_by_hist", date, ".pdf"), width = 6, height = 6)
hist(PCC_inds_w_by$by, xlab = "birth year", col = colors[1], main = NULL, cex.lab = 2, cex.axis = 1.5)
dev.off()

# PCC
PCC_pca_inds <- sep_pop$PCC@ind.names
PCC_inds_w_by <- PCC_inds_w_by[-which(PCC_inds_w_by$id == "PCC_300_P3"),]
match(PCC_inds_w_by$id, PCC_pca_inds) # PCC_pca_inds has PCC_300_P8
color_vec <- rep("gray", length(PCC_pca_inds))
color_vec[match(PCC_inds_w_by$id, PCC_pca_inds)] <- met.brewer("Greek", n = 70)

pdf(file =paste0("../vcf_filtering/PCC_PCA_maf_LD_flt_all_by_", date, ".pdf"), width = 6, height = 4)
mat <- matrix(c(1,1,1,1,1,1,1,1,1,2,2,2), ncol = 4)  
layout(mat)

plot(x = pca_PCC$scores[, 1], y = pca_PCC$scores[, 2],
     col = alpha(color_vec, 0.7), 
     cex = 1, pch = 19, cex.lab = 1.5,
     xlab = "PC1, 3.2% variation", ylab = "PC2, 3.1% variation")
hist(PCC_inds_w_by$by, xlab = "birth year", col = colors[1], main = NULL)

dev.off()

#ELF
ELF_pca_inds <- sep_pop$ELF@ind.names
ELF_inds_w_by[which(ELF_inds_w_by$id == "ELF_544_P10"),"id"] <- "ELF_544_P5"
color_vec <- rep("gray", length(ELF_pca_inds))
color_vec[match(ELF_inds_w_by$id, ELF_pca_inds)] <- met.brewer("Greek", n = 491)

pdf(file =paste0("../vcf_filtering/ELF_PCA_maf_LD_flt_all_by_", date, ".pdf"), width = 6, height = 4)
mat <- matrix(c(1,1,1,1,1,1,1,1,1,2,2,2), ncol = 4)  
layout(mat)
plot(x = pca_ELF$scores[, 1], y = pca_ELF$scores[, 2],
     col = alpha(color_vec, 0.7), 
     cex = 1, pch = 19, cex.lab = 1.5,
     xlab = "PC1, 3.5% variation", ylab = "PC2, 3.0% variation")

hist(ELF_inds_w_by$by, xlab = "birth year", col = colors[2], main = NULL)

dev.off()
# ------------------------------------------------------------------------------------------------------------

## Calculate Inbreeding ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
# multiple- locus heterozygosity (MLH), which is calculated as the proportion of genotyped loci that are heterozygous (Szulkin, Bierne, & David, 2010
# F can also be measured using the diagonal elements of a genomic relatedness matrix (Fgrm) (Bérénos et al., 2016; Huisman et al., 2016; Powell, Visscher, & Goddard, 2010; Pryce, Haile- Mariam, Goddard, & Hayes, 2014; Yang et al., 2010), which can be calculated with unmapped loc
# Kardos et al. 2016

# # MLH
t_gt_flt_mlh <- t_gt_flt
t_gt_flt_mlh[which(t_gt_flt_mlh == 2)] <- 0
# 
# mlh <- MLH(t_gt_flt_mlh)

# get snps and fitness for each population
gt_PCC_g2 <- t_gt_flt_mlh[which(grepl("PCC_", rownames(t_gt_flt_mlh))),]
gt_ELF_g2 <- t_gt_flt_mlh[which(grepl("ELF_", rownames(t_gt_flt_mlh))),]

PCC_off_g2 <- unlist(PCC.off)[which(grepl("PCC_", names(unlist(PCC.off))))]
gt_PCC_g2_sub <- gt_PCC_g2[match(names(PCC_off_g2), rownames(gt_PCC_g2)),]

ELF_off_g2 <- unlist(ELF.off)[which(grepl("ELF_", names(unlist(ELF.off))))]
gt_ELF_g2_sub <- gt_ELF_g2[match(names(ELF_off_g2), rownames(gt_ELF_g2)),]

# calculate g2 and boostraps to estimate CI

g2_PCC <- g2_snps(gt_PCC_g2_sub, nboot = 1e3)
g2_ELF <- g2_snps(gt_ELF_g2_sub, nboot = 1e3)

# variance in het
het_PCC <- sMLH(gt_PCC_g2_sub)
het_ELF <- sMLH(gt_ELF_g2_sub)

het_var_PCC <- var(het_PCC)
het_var_ELF <- var(het_ELF)

# linear model
mod_PCC <- lm(PCC_off_g2 ~ het_PCC)
mod_ELF <- lm(ELF_off_g2 ~ het_ELF)

coef(mod_PCC)[2] # beta

# r^2 between fitness and het
Wh_PCC <- cor(PCC_off_g2, predict(mod_PCC))^2
Wh_ELF <- cor(ELF_off_g2, predict(mod_ELF))^2

# r^2 between inbreeding and sMLH 
hf_PCC <- r2_hf(genotypes = gt_PCC_g2_sub, type = "snps", nboot = 1e3)
hf_ELF <- r2_hf(genotypes = gt_ELF_g2_sub, type = "snps", nboot = 1e3)

Wf_PCC <- r2_Wf(genotypes=gt_PCC_g2_sub, trait=PCC_off_g2, family="gaussian", type = "snps", nboot = 1e3)
Wf_ELF <- r2_Wf(genotypes=gt_ELF_g2_sub, trait=ELF_off_g2, family="gaussian", type = "snps", nboot = 1e3)

# 
# ------------------------------------------------------------------------------------------------------------
