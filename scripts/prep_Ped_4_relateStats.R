library(dplyr)
library(nadiv)

# prepping pedigrees for RelateStat

## load pedigrees ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
load("../pedigree_reconstruction/ELF_pedigree_results_05262023.Robj", verbose = T)
ELF_pedigree <- ELF_results[[4]]$Pedigree

load("../pedigree_reconstruction/PCC_pedigree_results_05262023.Robj", verbose = T)
PCC_pedigree <- PCC_results[[4]]$Pedigree

# ------------------------------------------------------------------------------------------------------------

## change NA to * ----------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
ELF_pedigree[which(is.na(ELF_pedigree$dam)),"dam"] <- "*"
ELF_pedigree[which(is.na(ELF_pedigree$sire)),"sire"] <- "*"

PCC_pedigree[which(is.na(PCC_pedigree$dam)),"dam"] <- "*"
PCC_pedigree[which(is.na(PCC_pedigree$sire)),"sire"] <- "*"

ELF_ped_out <- ELF_pedigree[,1:3]
PCC_ped_out <- PCC_pedigree[,1:3]

# ------------------------------------------------------------------------------------------------------------

## add sex and year metadata ----------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

ELF_ped_out <- merge(ELF_ped_out, ELF_results[[4]]$LifeHistPar, by= "id", all.x = TRUE, all.y = TRUE)[,c("id", "dam", "sire", "Sex", "BirthYear", "BY.est")]
ELF_ped_out[which(ELF_ped_out$BirthYear == -999),"BirthYear"] <- "*"
ELF_ped_out[which(is.na(ELF_ped_out$BY.est)),"BY.est"] <- "*"
ELF_ped_out[which(ELF_ped_out$Sex == 3),"Sex"] <- "*"
ELF_ped_out[which(is.na(ELF_ped_out$Sex)),"Sex"] <- "*"

ELF_merged_df <- merge(x = ELF_ped_out, 
                   y = ELF_results[[4]]$DummyIDs, 
                   by = "id", all.x = TRUE, all.y = TRUE)

ELF_merged_df$sex[!is.na(ELF_merged_df$Sex.x)] <- ELF_merged_df$Sex.x[!is.na(ELF_merged_df$Sex.x)]
ELF_merged_df$sex[!is.na(ELF_merged_df$Sex.y)] <- ELF_merged_df$Sex.y[!is.na(ELF_merged_df$Sex.y)]

ELF_final_df <- ELF_merged_df[,c("id", "dam.x", "sire.x", "Sex.x", "BY.est.x")]
colnames(ELF_final_df) <- c("id", "dam", "sire", "sex", "cohort")


PCC_ped_out <- merge(PCC_ped_out, PCC_results[[4]]$LifeHistPar, by= "id", all.x = TRUE, all.y = TRUE)[,c("id", "dam", "sire", "Sex", "BirthYear", "BY.est")]
PCC_ped_out[which(PCC_ped_out$BirthYear == -999),"BirthYear"] <- "*"
PCC_ped_out[which(is.na(PCC_ped_out$BY.est)),"BY.est"] <- "*"
PCC_ped_out[which(PCC_ped_out$Sex == 3),"Sex"] <- "*"
PCC_ped_out[which(is.na(PCC_ped_out$Sex)),"Sex"] <- "*"

PCC_merged_df <- merge(x = PCC_ped_out, 
                       y = PCC_results[[4]]$DummyIDs, 
                       by = "id", all.x = TRUE, all.y = TRUE)

PCC_merged_df$sex[!is.na(PCC_merged_df$Sex.x)] <- PCC_merged_df$Sex.x[!is.na(PCC_merged_df$Sex.x)]
PCC_merged_df$sex[!is.na(PCC_merged_df$Sex.y)] <- PCC_merged_df$Sex.y[!is.na(PCC_merged_df$Sex.y)]

PCC_final_df <- PCC_merged_df[,c("id", "dam.x", "sire.x", "Sex.x", "BY.est.x")]
colnames(PCC_final_df) <- c("id", "dam", "sire", "sex", "cohort")
PCC_final_df[which(is.na(PCC_final_df$sex)), "sex"] <- "*"

# convert 1s and 2s to Fs and Ms 

PCC_final_df[which(PCC_final_df$sex == "1"),"sex"] <- "F"
PCC_final_df[which(PCC_final_df$sex == "2"),"sex"] <- "M"
ELF_final_df[which(ELF_final_df$sex == "1"),"sex"] <- "F"
ELF_final_df[which(ELF_final_df$sex == "2"),"sex"] <- "M"

# ------------------------------------------------------------------------------------------------------------

## define "ancestors" and "focal cohort" ----------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
# ancestors
hist(as.numeric(ELF_final_df$cohort))
ELF_anc <- ELF_final_df[which(ELF_final_df$cohort < 2010),]
ELF_anc <- ELF_anc %>% filter(cohort != "*") %>% select(id)

hist(as.numeric(PCC_final_df$cohort))
PCC_anc <- PCC_final_df[which(PCC_final_df$cohort < 2010),]
PCC_anc <- PCC_anc %>% filter(cohort != "*") %>% select(id)


# focal cohort

hist(as.numeric(ELF_final_df$cohort))
ELF_focal <- ELF_final_df[which(ELF_final_df$cohort > 2015),]
ELF_focal <- ELF_focal %>% filter(cohort != "*") %>% select(id)

hist(as.numeric(PCC_final_df$cohort))
PCC_focal <- PCC_final_df[which(PCC_final_df$cohort < 2015),]
PCC_focal <- PCC_focal %>% filter(cohort != "*") %>% select(id)

## make pedigree relatedness matrix ----------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
PCC_nadiv_ped <- prepPed(PCC_final_df, gender = "sex")
PCC_relMat <- makeA(PCC_nadiv_ped[1:3])

ELF_nadiv_ped <- prepPed(ELF_final_df, gender = "sex")
ELF_relMat <- makeA(ELF_nadiv_ped[1:3])
as.matrix(ELF_relMat)
# ------------------------------------------------------------------------------------------------------------

## export output ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
write.table(ELF_final_df, file = "../pedigree_skew/ELF_ped_df.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(PCC_final_df, file = "../pedigree_skew/PCC_ped_df.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(ELF_anc, file = "../pedigree_skew/ELF_anc.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(PCC_anc, file = "../pedigree_skew/PCC_anc.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(as.matrix(PCC_relMat), file = "../pedigree_skew/PCC_relMat.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(as.matrix(ELF_relMat), file = "../pedigree_skew/ELF_relMat.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(ELF_focal, file = "../pedigree_skew/ELF_focal.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(PCC_focal, file = "../pedigree_skew/PCC_focal.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)



# ------------------------------------------------------------------------------------------------------------
