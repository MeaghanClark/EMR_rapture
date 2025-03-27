## File name: filter_vcf.R
## Purpose: Load and filter VCF file from STACKS for EMR project
## Based on code by R.H. Toczydlowski
## M. I. Clark, June 2022
## Last updated: 04/22/2024

# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# This script loads a .vcf file and performs filtering and other analyses
#
# -----------------------------------------------------------------------------------------------------------------------------------------------------------

## Script set up --------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------

## define color palette 
palette = "Archambault"
met.brewer(palette, n=6, type="continuous")
MetBrewer::colorblind.friendly(palette)
colors <- met.brewer(palette, n = 6)

date <- format(Sys.Date(), "%m%d%Y")
# ----------------------------------------------------------------------------------------------------------------------------------------------------------

## Load custom functions and required libraries  -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------

source("./vcf_funcs.R")

# ----------------------------------------------------------------------------------------------------------------------------------------------------------

## Read in files -------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------
## read vcf
vcf <- read.vcfR("../rapture_ref_map_output/emr_rapture_filtered_15.vcf", verbose = T)

## read filtered bed file
targets <- read.csv("../vcf_filtering/Spaced_baits.fas_mrg_buf.bed", sep = "\t", header=FALSE)
colnames(targets) <- c("chrom", "start", "end")

## load in depth files from vcftools

# generated using: vcftools --vcf ./rapture_ref_map_output/RE_RUN_725.populations.snps.vcf --site-mean-depth --out ./vcf_filtering/depth_vcftools
depth <- read.csv("/Users/meaghan/Desktop/EMR_rapture/vcf_filtering/depth_vcftools_12924.ldepth.mean", sep = "\t")

# generated using: vcftools --vcf ./rapture_ref_map_output/RE_RUN_725.populations.snps.vcf --depth --out ./vcf_filtering/depth_vcftools
ind_depth <- read.csv("/Users/meaghan/Desktop/EMR_rapture/vcf_filtering/depth_vcftools_12924.idepth", sep = "\t")

site_depth <- read.csv("/Users/meaghan/Desktop/EMR_rapture/vcf_filtering/site_depth.ldepth", sep = "\t")

# ------------------------------------------------------------------------------------------------------------

## Depth Exploration ---------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

# site  depth
hist(site_depth$SUM_DEPTH)

# Mean individual depth

# consider removing these individuals with mean depth less than 15
low_depth_inds <- ind_depth[which(ind_depth$MEAN_DEPTH < 15),]

# should drop these individuals before filtering SNPs any further

## depth from VCFtools, per SNP
mean(depth$MEAN_DEPTH) # mean site depth
range(depth$MEAN_DEPTH)
range(depth$VAR_DEPTH)

hist(depth$MEAN_DEPTH, main = NULL, xlab = "average depth per SNP")
hist(depth$VAR_DEPTH)

chromosomes <- unique(targets[,1])
pdf(file = "../vcf_filtering/filtered_vcf/raw_depth_012924.pdf", width = 8, height = 8)
for(chrom in chromosomes){
  data <- depth[which(depth$CHROM == chrom),]
  if(dim(data)[1] > 0){
    par(mfrow=c(1,1))
    plot(x = data$POS, y = data$MEAN_DEPTH, 
         xlab = "position", 
         ylab = "mean depth", 
         ylim = c(0,90),
         xlim = c(135901900, 135902520),
         main = chrom, 
         pch = 19, 
         col = adjustcolor("black", 0.2))
    abline(h = 10, col = "red", lwd = 3)
  }
}
dev.off()

chromosomes <- c("Scate-ma1", "Scate-Z", "Scate-un351")
pdf(file = "../vcf_filtering/filtered_vcf/raw_depth_example_012924.pdf", width = 12, height = 5)
par(mfrow=c(1,3))
for(chrom in chromosomes){
  data <- depth[which(depth$CHROM == chrom),]
  if(dim(data)[1] > 0){
    plot(x = data$POS, y = data$MEAN_DEPTH, 
         xlab = "position", 
         ylab = "mean depth", 
         ylim = c(0,90),
         main = chrom, 
         pch = 19, 
         col = adjustcolor("black", 0.2))
    abline(h = 10, col = "red", lwd = 3)
  }
}

dev.off()


# depth per ind 
pdf(file = "../vcf_filtering/filtered_vcf/depth_per_ind_12924.pdf", width = 12, height = 5)
par(mfrow=c(1,2))
hist(ind_depth$MEAN_DEPTH, main = NULL, xlab = "average depth per individual")
plot(x = ind_depth$MEAN_DEPTH, y = ind_depth$N_SITES, main = NULL, 
     xlab = "average depth per individual", 
     ylab = "Number of site with genotype calls")
abline(v = 15, col = "red")
dev.off()

# ----------------------------------------------------------------------------------------------------------------------------------------------------------

## Verify that all SNPs are biallelic ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------

## Minor allele frequency, pre-filtering
maf <- maf(vcf)
hist(maf[,4], xlab = "MAF from unfiltered SNPs", main = paste0("average maf = ", mean(maf[,4])))

#check if there are more than 2 alleles at any loci
gt.temp <- vcfR::extract.gt(vcf, return.alleles = T)
# NA handling (NA here means SNP was not called)
gt.temp[gt.temp == "."] <- NA
gt.temp <- as.data.frame(t(gt.temp), stringsAsFactors = FALSE)
gt.temp$sampid = rownames(gt.temp) 
gt.temp <- gt.temp %>% dplyr::select(sampid,everything())
gt.temp.long <- gt.temp %>% tidyfast::dt_pivot_longer(.,names_to="SNPid",values_to="genotype",cols=2:ncol(gt.temp), factor_key = T) %>% na.omit()
gt.temp.long %>% distinct(SNPid,genotype) %>% group_by(SNPid) %>% 
  mutate(all_alleles = paste0(genotype, collapse = "")) %>% mutate(all_alleles = gsub("/","",all_alleles)) %>% 
  mutate(number_unique_alleles_atSNP = length(unique(strsplit(all_alleles, "")[[1]]))) %>% 
  group_by(number_unique_alleles_atSNP) %>% summarise(number_of_SNPs=n())

# all SNPs are biallelic! 

# ----------------------------------------------------------------------------------------------------------------------------------------------------------

## Filter VCF ----------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------

## Extract genotypes
genotypes <- vcfR::extract.gt(vcf, return.alleles = F)

## Create different data subsets based on depth

# all: genotypes, vcf@fix
positions <- vcf@fix

## Transpose and make into data.frame 
gt.unflt <- process_gt(genotypes)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------

## Test for Excess Het -------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------

# convert to genotype count format? 
geno.counts <- matrix(nrow = dim(gt.unflt)[1], ncol = 3)
colnames(geno.counts) <- c("AA", "AB", "BB")
for(i in 1:nrow(gt.unflt)){
  geno.counts[i,"AA"] <- sum(gt.unflt[i,] == 0, na.rm = TRUE)
  geno.counts[i,"AB"] <- sum(gt.unflt[i,] == 1, na.rm = TRUE)
  geno.counts[i,"BB"] <- sum(gt.unflt[i,] == 2, na.rm = TRUE)
}

# Warning, this can take a while to run! 
ExHet_pvals <- HWExactStats(X = geno.counts, x.linked = FALSE, alternative = "greater", plinkcode=FALSE, verbose = FALSE)

sum(ExHet_pvals < 0.05)

gt.hetflt <- gt.unflt[-which(ExHet_pvals < 0.05),]

# ----------------------------------------------------------------------------------------------------------------------------------------------------------

## Generate SFS# -------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------

pdf(file ="../vcf_filtering/filtered_vcf/SFS_ExHetflt.pdf", width = 10, height = 6)
par(mfrow = c(1,2))
plot_SFS(gt.unflt, col = colors[2], 
         main = "all SNPs") # all
plot_SFS(gt.hetflt, col = colors[2], 
         main = "ExHet Filter") # all

dev.off()

# population specific
raw_names <- unlist(strsplit(colnames(gt.hetflt), "_"))
sites <- raw_names[which(raw_names == "ELF" | raw_names == "PCC")]
gt_ELF_hetflt <- gt.hetflt[,which(sites == "ELF")]
gt_ELF_hetflt <- gt.hetflt[,which(sites == "PCC")]

pdf(file ="../vcf_filtering/filtered_vcf/SFS_pops_12924.pdf", width = 10, height = 6)
par(mfrow = c(1,2))
plot_SFS(gt_ELF_hetflt, col = colors[3], 
         main = "ELF, all SNPs", 
         #xlim = c(0, 1000), 
         #ylim = c(0, 65), 
         rm.zero = TRUE) # all
plot_SFS(gt_PCC_hetflt, col = colors[5], 
         main = "PCC, all SNPs", 
         #xlim = c(0, 1000), 
         #ylim = c(0,250), 
         rm.zero = TRUE) # lower
dev.off()

# ----------------------------------------------------------------------------------------------------------------------------------------------------------

## Calculate genotyping error rate
# ----------------------------------------------------------------------------------------------------------------------------------------------------------
# duplicate individuals for PCC
PCC_dupe_inds <- read.csv("../pedigree_reconstruction/PCC_LifeHistData_dupes.csv", header=FALSE)
dupe_inds <- rbind.data.frame(PCC_dupe_inds, c("ELF_544_P5", "ELF_544_P10", ""))

error_rate <- calc_gt_error_rate_flt(dupe_inds, gt.hetflt, tri = TRUE)
# hist(error_rate)
# abline(v= mean(error_rate), col = "green") # 0.002348661

# ----------------------------------------------------------------------------------------------------------------------------------------------------------

## Select a single SNP per targeted region # 
# ----------------------------------------------------------------------------------------------------------------------------------------------------------
## Choose a single SNP randomly per targeted region (sample_snps)

gt_samp <- sample_snps(genotypes = gt.hetflt, target_regions = targets, loci_info = vcf@fix, n_inds = dim(gt.hetflt)[2], z.chrom = TRUE) # "67.2007255139057 percent of targeted regions have a sequenced SNP"

#67.2007255139057 percent of targeted regions have a sequenced SNP

# ----------------------------------------------------------------------------------------------------------------------------------------------------------

## Binomial tests for contaminated samples  ----------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------

pdf(file = "../vcf_filtering/binom_test_all_inds.pdf")
for(i in 2:length(colnames(vcf@gt))){
  binom_test_for_allele_counts(vcf@gt[,which(colnames(vcf@gt) == colnames(vcf@gt)[i])], colnames(vcf@gt)[i]) # weirdo
}
dev.off()

# ELF_639
binom_test_for_allele_counts(vcf@gt[,which(colnames(vcf@gt) == "ELF_639")], "ELF_639") # weirdo
# look into details
ind = "ELF_328"
ind_gt <- vcf@gt[,which(colnames(vcf@gt) == "ELF_328")]
p_val <- vector(length = length(ind_gt))

for(l in 1:length(ind_gt)){
  gt <- strsplit(ind_gt[l], split = ":")[[1]][1]
  if(gt == "0/1" | gt == "1/0"){
    ac <- strsplit(ind_gt[l], split = ":")[[1]][3]
    ac_ref <- strsplit(ac, split = ",")[[1]][1]
    ac_alt <- strsplit(ac, split = ",")[[1]][2]
    if(!is.na(as.numeric(ac_ref))){
      test <- binom.test(x = c(as.numeric(ac_ref), as.numeric(ac_alt)), p = 0.5, alternative = "two.sided")
      p_val[l] <- test$p.value
    }else{
      p_val[l] <- NA
    }
  }else{p_val[l] <- NA}
}

which(p_val < 0.05)
length(na.omit(p_val)) * 0.05
sum(na.omit(p_val) < 0.05) / length(na.omit(p_val))

load(file="../vcf_filtering/pop_gen_stats.Robj", verbose = T)
pop_gen_stats[which(pop_gen_stats$id == "ELF_139"),] # high Fgrm = ELF_269


par(mfrow=c(2,3))

binom_test_for_allele_counts(vcf@gt[,which(colnames(vcf@gt) == "ELF_274")], "ELF_274") # weirdo
binom_test_for_allele_counts(vcf@gt[,which(colnames(vcf@gt) == "ELF_824")], "ELF_824") # wow these look different 
binom_test_for_allele_counts(vcf@gt[,which(colnames(vcf@gt) == "ELF_165")], "ELF_165") # normal
binom_test_for_allele_counts(vcf@gt[,which(colnames(vcf@gt) == "ELF_401")], "ELF_401") # weirdo
binom_test_for_allele_counts(vcf@gt[,which(colnames(vcf@gt) == "ELF_729")], "ELF_729") # normal
binom_test_for_allele_counts(vcf@gt[,which(colnames(vcf@gt) == "PCC_138")], "PCC_138") # weird
# unexpected read proportions for an uncontaminated individual 

par(mfrow=c(2,3))
binom_test_for_allele_counts(vcf@gt[,which(colnames(vcf@gt) == "ELF_135")], "ELF_135")
binom_test_for_allele_counts(vcf@gt[,which(colnames(vcf@gt) == "ELF_358")], "ELF_358")
binom_test_for_allele_counts(vcf@gt[,which(colnames(vcf@gt) == "PCC_100")], "PCC_100")
binom_test_for_allele_counts(vcf@gt[,which(colnames(vcf@gt) == "PCC_277")], "PCC_277")
binom_test_for_allele_counts(vcf@gt[,which(colnames(vcf@gt) == "PCC_220")], "PCC_220")
binom_test_for_allele_counts(vcf@gt[,which(colnames(vcf@gt) == "ELF_508")], "ELF_508")

# excess at 1 --> lower coverage sites with less power to detect heterozygosity? 

# what about very very similar individuals? 
genotypes_big <- vcfR::extract.gt(vcf, return.alleles = F)
gt_big <- process_gt(genotypes_big)

sum(gt_big[,which(colnames(gt_big)=="PCC_4")] == gt_big[,which(colnames(gt_big)=="PCC_283")], na.rm = T) # 9904 / 9917

sum(unlist(lapply(X = 1:nrow(gt_big), FUN = function(x) {
  !is.na(gt_big[,which(colnames(gt_big)=="PCC_4")][x]) && !is.na(gt_big[,which(colnames(gt_big)=="PCC_283")][x])
})))

sum(gt_big[,which(colnames(gt_big)=="PCC_86")] == gt_big[,which(colnames(gt_big)=="PCC_149")], na.rm = T) # 9933 / 9955
sum(unlist(lapply(X = 1:nrow(gt_big), FUN = function(x) {
  !is.na(gt_big[,which(colnames(gt_big)=="PCC_86")][x]) && !is.na(gt_big[,which(colnames(gt_big)=="PCC_149")][x])
})))

sum(gt_big[,which(colnames(gt_big)=="PCC_99")] == gt_big[,which(colnames(gt_big)=="PCC_210")], na.rm = T) # 9885 / 9898
sum(unlist(lapply(X = 1:nrow(gt_big), FUN = function(x) {
  !is.na(gt_big[,which(colnames(gt_big)=="PCC_99")][x]) && !is.na(gt_big[,which(colnames(gt_big)=="PCC_210")][x])
})))

sum(gt_big[,which(colnames(gt_big)=="PCC_138")] == gt_big[,which(colnames(gt_big)=="PCC_210")], na.rm = T) # 6351 /9117
sum(unlist(lapply(X = 1:nrow(gt_big), FUN = function(x) {
  !is.na(gt_big[,which(colnames(gt_big)=="PCC_138")][x]) && !is.na(gt_big[,which(colnames(gt_big)=="PCC_210")][x])
})))

# ----------------------------------------------------------------------------------------------------------------------------------------------------------

## Remove erroneous individuals from gt_samp
# ----------------------------------------------------------------------------------------------------------------------------------------------------------
# remove non-PCC individuals and individuals without metadata from PCC
# "PCC_297" no pit tag, not assigned in preliminary pedigree runs, remove
# "PCC_301" 601263359 remove, not from PCC
# "PCC_306" 604078853 not from PCC? 
# added 5/25/2023: 
# "PCC_283", likely same individual as PCC_4
# "PCC_86", likely same individual as 149
# "PCC_99", likely same individual as 210
# "ELF_274", likely contaminated, see binomial tests
# "ELF_401", likely contaminated, see binomial tests
# "ELF_138", likely contaminated, see binomial tests
# "ELF_824", likely contaminated, see binomial tests
# "ELF_639", likely contaminated, see binomial tests

gt_samp <- subset(gt_samp, select = -c(PCC_297, PCC_301, PCC_306, PCC_283, PCC_86, PCC_99, ELF_274, ELF_401, ELF_138, ELF_824, ELF_639))
# 2223 1080

## Remove individuals with low mean coverage
gt_samp <- subset(gt_samp, select = -c(ELF_141, ELF_242, ELF_289, ELF_324, ELF_333, ELF_342, ELF_893, ELF_906, ELF_917, PCC_285, PCC_308, PCC_73))

save(gt_samp, 
     file = paste0("../vcf_filtering/unfiltered_gt_", date, ".Robj"))

# ----------------------------------------------------------------------------------------------------------------------------------------------------------

## FILTERING
# ----------------------------------------------------------------------------------------------------------------------------------------------------------

# Plot missingness 
pdf(file =paste0("../vcf_filtering/filtered_vcf/pre_filter_missingness_noZ_12924", date, ".pdf"), width = 10, height = 6)
par(mfrow = c(1, 2))
plot_snp_missingness(gt_samp, 0.1, main = "targets, all SNPs", xlim = c(0,0.2))
plot_ind_missingness(gt_samp, 0.2, main = "targets, all SNPs", xlim = c(0,0.8))
dev.off()

# filter SNPS then individuals
dim(gt_samp) # 2223 1067
gt_samp_f1 <- filter_snp_missingness(gt_samp, 0.1, plot = TRUE) # drop loci with > 10% missing data
dim(gt_samp_f1) # 2223 1067
gt_samp_f2 <- filter_ind_missingness(gt_samp_f1, 0.2, plot = TRUE) # drop individuals with > 30% missing data
dim(gt_samp_f2) # 2223 1067

# ----------------------------------------------------------------------------------------------------------------------------------------------------------

## Split into populations
# ----------------------------------------------------------------------------------------------------------------------------------------------------------
raw_names <- unlist(strsplit(colnames(gt_samp_f2), "_"))
sites <- raw_names[which(raw_names == "ELF" | raw_names == "PCC")]
gt_flt_ELF <- gt_samp_f2[,which(sites == "ELF")] # 2468  784 # 2223  780
gt_flt_PCC <- gt_samp_f2[,which(sites == "PCC")] # 2468  289 # 2223  288
# ----------------------------------------------------------------------------------------------------------------------------------------------------------

## Calculate MAF
# ----------------------------------------------------------------------------------------------------------------------------------------------------------
maf <- calc_maf(gt_flt_PCC)
hist(maf)
sum(maf > 0.1)/length(maf) # proportion of SNPs with maf > 0.1

# maf ELF
maf <- calc_maf(gt_flt_ELF)
hist(maf)
sum(maf > 0.1)/length(maf) # proportion of SNPs with maf > 0.1
# ----------------------------------------------------------------------------------------------------------------------------------------------------------

## Save filtered data --------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------
# PCC
PCC_filt_gt_all <- gt_flt_PCC
save(PCC_filt_gt_all, file = paste0("../vcf_filtering/PCC_filtered_gt_all_snps_noZ_", date, ".Robj"))

# ELF
ELF_filt_gt_all <- gt_flt_ELF
save(ELF_filt_gt_all, file = paste0("../vcf_filtering/ELF_filtered_gt_all_snps_noZ_", date, ".Robj"))

# ----------------------------------------------------------------------------------------------------------------------------------------------------------

### END SCRIPT 