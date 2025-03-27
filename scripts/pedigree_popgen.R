## File name: pedigree_popgen.R
## Purpose: calculate basic popgen stats from RAPTURE data and relate to the pedigree (including pedigree relatedness)
## M. I. Clark, November 2022
## Last updated: 11/16/2022

# note: not independent, still requires some data objects from explore_pedigrees.R

## Loading ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
library(devtools)
library(popgenstuff)
library(data.table)
library(hexbin)
library(scales)

date <- 11142022

load(file = paste0("../vcf_filtering/PCC_filtered_gt_all_snps_noZ", date, ".Robj"), verbose = T) # PCC_filt_gt_all
load(file = paste0("../vcf_filtering/ELF_filtered_gt_all_snps_noZ", date, ".Robj"), verbose = T) # ELF_filt_gt_all

load(file = paste0("../vcf_filtering/unfiltered_gt", date, ".Robj"), verbose = T)

PCC.off <- offspring_per_ind(PCC_pedigree)
ELF.off <- offspring_per_ind(ELF_pedigree_flt)

## Heterozygosity  ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

PCC.het <- calcHet(t(gt_PCC_samp))

ELF.het <- calcHet(t(gt_ELF_samp))

pdf(file = "../pedigree_exploration/het_at_ped_snps.pdf", width = 10, height = 6)
par(mfrow=c(1,2))
plot(PCC.het, ylim = c(0.1, 0.55), pch = 19, col = alpha("black", 0.25), 
     xlab = "individuals", ylab = "heterozygosity at pedigree SNPs", 
     main = "PCC")
abline(h = mean(PCC_het), col = "red", lty = 2)
plot(ELF.het, ylim = c(0.1, 0.55), pch = 19, col = alpha("black", 0.25), 
     xlab = "individuals", ylab = "heterozygosity at pedigree SNPs", 
     main = "ELF")
abline(h = mean(ELF_het), col = "red", lty = 2)
dev.off()

# PCC
PCC.off.gt <- unlist(PCC.off)[grepl("PCC", names(unlist(PCC.off)))]
PCC.het[match(names(PCC.off.gt), names(PCC.het))]
PCC_dat <- cbind.data.frame(PCC.off.gt, PCC.het[match(names(PCC.off.gt), names(PCC.het))])
colnames(PCC_dat) <- c("offspring", "het")
par(mfrow=c(1,1))
plot(PCC_dat$offspring~PCC_dat$het, pch = 19, col = alpha("black", 0.25), main = "PCC")
abline(lm(offspring~het, data = PCC_dat), col = "red", lty = 2)


# ELF
ELF.off.gt <- unlist(ELF.off)[grepl("ELF", names(unlist(ELF.off)))]
ELF.het[match(names(ELF.off.gt), names(ELF.het))]
ELF_dat <- cbind.data.frame(ELF.off.gt, ELF.het[match(names(ELF.off.gt), names(ELF.het))])
colnames(ELF_dat) <- c("offspring", "het")
par(mfrow=c(1,1))
plot(ELF_dat$offspring~ELF_dat$het, pch = 19, col = alpha("black", 0.25), main = "ELF")
abline(lm(offspring~het, data = ELF_dat), col = "red", lty = 2)

pdf(file = "../pedigree_exploration/het_v_off_at_ped_snps.pdf", width = 10, height = 6)
par(mfrow=c(1,2))
plot(PCC_dat, pch = 19, col = alpha("black", 0.25), main = "PCC")
abline(lm(het~offspring, data = PCC_dat), col = "red", lty = 2)
plot(ELF_dat, pch = 19, col = alpha("black", 0.25), main = "ELF")
abline(lm(het~offspring, data = ELF_dat), col = "red", lty = 2)
dev.off()

## Theta  ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

calcThetaW(t(gt_PCC_samp))
calcThetaW(t(gt_ELF_samp))

## Calculate pedigree relatedness  ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

ELF.Rped <- CalcRped(ELF_pedigree, OUT="DF")
hist(ELF.Rped$R.ped)
range(ELF.Rped$R.ped)

PCC.Rped <- CalcRped(PCC_pedigree, OUT="DF")
hist(PCC.Rped$R.ped)
range(PCC.Rped$R.ped)

## Pairwise pi  ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

PCC.PWP <- freqs2pairwisePi(t(gt_PCC_samp))

ELF.PWP <- freqs2pairwisePi(t(gt_ELF_samp))

colnames(PCC.PWP) <- colnames(gt_PCC_samp)
rownames(PCC.PWP) <- colnames(gt_PCC_samp)
PCC.PWP_df <- as.data.frame(as.table(PCC.PWP))
colnames(PCC.PWP_df) <- c("IID1", "IID2", "PWP")


colnames(ELF.PWP) <- colnames(gt_ELF_samp)
rownames(ELF.PWP) <- colnames(gt_ELF_samp)
ELF.PWP_df <- as.data.frame(as.table(ELF.PWP))
colnames(ELF.PWP_df) <- c("IID1", "IID2", "PWP")


# plot combine with pedigree relatedness
PCC.Rel.both <- merge(data.table(PCC.PWP_df, key=c("IID1", "IID2")),
                  data.table(PCC.Rped, key=c("IID1", "IID2")), all.x=TRUE)
PCC.Rel.both <- as.data.frame(PCC.Rel.both)  # turn back into regular dataframe

round(cor(PCC.Rel.both[, c("PWP","R.ped")], 
          use="pairwise.complete"), 2)

plot(PCC.Rel.both$PWP ~ PCC.Rel.both$R.ped, pch = 19, col = alpha("black", 0.05), 
     xlab = "pedigree relatedness", 
     ylab = "pairwise pi", 
     main = "PCC")


ELF.Rel.both <- merge(data.table(ELF.PWP_df, key=c("IID1", "IID2")),
                      data.table(ELF.Rped, key=c("IID1", "IID2")), all.x=TRUE)
ELF.Rel.both <- as.data.frame(ELF.Rel.both)  # turn back into regular dataframe

round(cor(ELF.Rel.both[, c("PWP","R.ped")], 
          use="pairwise.complete"), 2)

plot(ELF.Rel.both$PWP ~ ELF.Rel.both$R.ped, pch = 19, col = alpha("black", 0.05), 
     xlab = "pedigree relatedness", 
     ylab = "pairwise pi", 
     main = "ELF")


# save plot
pdf(file = "../pedigree_exploration/pwp_v_ped_rel.pdf", width = 10, height = 6)
par(mfrow=c(1,2))
plot(PCC.Rel.both$PWP ~ PCC.Rel.both$R.ped, pch = 19, col = alpha("black", 0.05), 
     xlab = "pedigree relatedness", 
     ylab = "pairwise pi", 
     main = "PCC")
plot(ELF.Rel.both$PWP ~ ELF.Rel.both$R.ped, pch = 19, col = alpha("black", 0.05), 
     xlab = "pedigree relatedness", 
     ylab = "pairwise pi", 
     main = "ELF")
dev.off()

# hex plots 
pdf(file = "../pedigree_exploration/pwp_v_ped_rel_hex.pdf", width = 10, height = 6)
par(mfrow=c(1,2))

hexbin::hexbinplot(ELF.Rel.both$PWP ~ ELF.Rel.both$R.ped, 
                   xbins=100, aspect=1,
                   xlim=c(-.05,.9), ylim=c(-.2, .9),
                   xlab="Pedigree relatedness", ylab="pairwise pi",
                   trans=log10, inv=function(x) 10^x, 
                   colorcut=seq(0,1,length=14), maxcnt=10^6.5,
                   colramp = function(n) {grDevices::hcl.colors(n, palette='Berlin')})

hexbin::hexbinplot(PCC.Rel.both$PWP ~ PCC.Rel.both$R.ped, 
                   xbins=100, aspect=1,
                   xlim=c(-.05,.9), ylim=c(-.2, .9),
                   xlab="Pedigree relatedness", ylab="pairwise pi",
                   trans=log10, inv=function(x) 10^x, 
                   colorcut=seq(0,1,length=14), maxcnt=10^6.5,
                   colramp = function(n) {grDevices::hcl.colors(n, palette='Berlin')})
dev.off()


## Save data objects  ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
save(PCC.Rped, ELF.Rped, PCC.het, ELF.het, PCC.PWP, ELF.PWP, file = "../pedigree_exploration/popgen_stats.Robj")

