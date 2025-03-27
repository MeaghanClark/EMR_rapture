## File name: explore_pedigree.R
## Purpose: explore pedigrees for EMR project
## M. I. Clark, October 2022
## Last updated: 10/21/2022

## Loading ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
library(MetBrewer)

## define color palette 
palette = "Archambault"
met.brewer(palette, n=6, type="continuous")
MetBrewer::colorblind.friendly(palette)
colors <- met.brewer(palette, n = 6)


## date we want to analyze
# date <- 10182022
# date <- 11142022
# date <- "02102023"
# date <- 40502023
date <- "05252023"

## load data from make_prelim_pedigree.R
# contain: 
# 1. genotype information
# 2. metadata
# 3. duplicate individuals
# 4. results from full pedigree run
# 5. vcf filtering date

load(paste0("../pedigree_reconstruction/ELF_pedigree_results_", date, ".Robj"), verbose = T)
load(paste0("../pedigree_reconstruction/PCC_pedigree_results_", date, ".Robj"), verbose = T)

ELF_ConfProb <- ELF_results[[5]]$ConfProb
PCC_ConfProb <- PCC_results[[5]]$ConfProb


ELF_results[[5]][["PedErrors"]]
PCC_results[[5]][["PedErrors"]]

pdf(file = paste0("../pedigree_reconstruction/ped_conf_mismatches_", date, ".pdf"), height = 8, width = 8)
par(mfrow= c(2,2))
barplot(height = PCC_results[[5]][["PedErrors"]][,3,2], 
        main = "PCC Sire mismatches", 
        ylab = "Proportion mismatch", 
        xlab = "individual + parent codes",         
        col = colors[1])

barplot(height = PCC_results[[5]][["PedErrors"]][,3,1], 
        main = "PCC Dam mismatches", 
        ylab = "Proportion mismatch", 
        xlab = "individual + parent codes", 
        col = colors[1])

barplot(height = ELF_results[[5]][["PedErrors"]][,3,2], 
        main = "ELF Sire mismatches", 
        ylab = "Proportion mismatch", 
        xlab = "individual + parent codes", 
        col = colors[2])

barplot(height = ELF_results[[5]][["PedErrors"]][,3,1], 
        main = "ELF Dam mismatches", 
        ylab = "Proportion mismatch", 
        xlab = "individual + parent codes", 
        col = colors[2])
dev.off()


pdf(file = paste0("../pedigree_reconstruction/ped_conf_false_neg_", date, ".pdf"), height = 8, width = 8)
par(mfrow= c(2,2))
barplot(height = PCC_results[[5]][["PedErrors"]][,1,2], 
        main = "PCC Sire", 
        ylab = "False Negative", 
        xlab = "individual + parent codes",         
        col = colors[1])
barplot(height = PCC_results[[5]][["PedErrors"]][,1,1], 
        main = "PCC Dam", 
        ylab = "False Negative", 
        xlab = "individual + parent codes", 
        col = colors[1])
barplot(height = ELF_results[[5]][["PedErrors"]][,1,2], 
        main = "ELF Sire", 
        ylab = "False Negative", 
        xlab = "individual + parent codes", 
        col = colors[2])
barplot(height = ELF_results[[5]][["PedErrors"]][,1,1], 
        main = "ELF Dam", 
        ylab = "False Negative", 
        xlab = "individual + parent codes", 
        col = colors[2])
dev.off()

pdf(file = paste0("../pedigree_reconstruction/ped_conf_false_pos_", date, ".pdf"), height = 8, width = 8)
par(mfrow= c(2,2))
barplot(height = PCC_results[[5]][["PedErrors"]][,2,2], 
        main = "PCC Sire", 
        ylab = "False Positive", 
        xlab = "individual + parent codes",         
        col = colors[1])
barplot(height = PCC_results[[5]][["PedErrors"]][,2,1], 
        main = "PCC Dam", 
        ylab = "False Positive", 
        xlab = "individual + parent codes", 
        col = colors[1])
barplot(height = ELF_results[[5]][["PedErrors"]][,2,2], 
        main = "ELF Sire", 
        ylab = "False Positive", 
        xlab = "individual + parent codes", 
        col = colors[2])
barplot(height = ELF_results[[5]][["PedErrors"]][,2,1], 
        main = "ELF Dam", 
        ylab = "False Positive", 
        xlab = "individual + parent codes", 
        col = colors[2])
dev.off()
