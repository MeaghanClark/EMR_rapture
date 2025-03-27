## File name: filter_vcf.R
## Purpose: Load and filter VCF file from STACKS for EMR project
## Based on code by R.H. Toczydlowski
## M. I. Clark, June 2022
## Last updated: 04/22/2024


## loading ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
## load libraries
library(vcfR)
library(dplyr)
library(tidyr)
library(tidyfast)
library("MetBrewer")
library(HardyWeinberg)

## define color palette 
palette = "Archambault"
met.brewer(palette, n=6, type="continuous")
MetBrewer::colorblind.friendly(palette)
colors <- met.brewer(palette, n = 6)

date <- format(Sys.Date(), "%m%d%Y")
# ------------------------------------------------------------------------------------------------------------

## Define custom functions ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

sample_snps <- function(genotypes, target_regions, loci_info, n_inds, z.chrom = TRUE){
  #testing
  # gt_ELF_samp <- sample_snps(genotypes = gt_ELF, target_regions = targets, loci_info = vcf@fix, n_inds = dim(gt_ELF)[2], z.chrom = TRUE)
   # genotypes <- gt.hetflt
   # target_regions <- targets
   # n_inds = dim(gt_PCC)[2]
   # loci_info <- vcf@fix
   # z.chrom <- TRUE
  # 
  #get rid of monomorphic loci
  bi_allelic <- vector()
  for (j in 1:nrow(genotypes)){
    if(length(na.omit(unique(genotypes[j,]))) > 1){
      bi_allelic <- c(bi_allelic, j)
    }else(print(paste0("found a monomorphic locus, j is ", j)))
  }
  
  genotypes <- genotypes[bi_allelic,]
  loci_info <- loci_info[bi_allelic,]
  
  #sample one snp per targeted region 
  gt.target <- as.data.frame(matrix(data=-1, nrow = 0, ncol = n_inds))
  colnames(gt.target) <- colnames(genotypes)
  no.data <- data.frame()
  if(z.chrom == TRUE){ # removes Z chromosome
    chromosomes <- unique(target_regions[,1])
    chromosomes <- chromosomes[-which(chromosomes == "Scate-Z")]
  }else{
    chromosomes <- unique(target_regions[,1])
  }
  
  for (chrom in chromosomes){ # need to edit code to go chromosome by chromosome because numbers will overlap
    chunk <- target_regions[which(target_regions[,1] == chrom),]
    gt.chunk <- data.frame()
    ids <- vector()
    for (i in 1:nrow(chunk)){
      start <- chunk[i,2]
      end <- chunk[i,3] 
      snps <- as.numeric(loci_info[which(loci_info[,1] == chrom & as.numeric(loci_info[,2]) > start & as.numeric(loci_info[,2]) < end),2])
      #print(length(snps))
      if(length(snps) == 0 ){
        print(paste0("no SNPs between ", start, " and ", end, " on chromosome ", chrom, ", i is ", i))
          no.data <- rbind(no.data, chunk[i,])
          # next
      }
      if(length(snps) >= 1){ # if there is a snp in the targeted region
        if(length(snps) > 1){ # subsample one snp
          target.snp <- sample(snps,1)
          keep.data <- genotypes[which(as.numeric(loci_info[,2]) == target.snp & loci_info[,1] == chrom),] # isolate genotypes for snp
          
        }
        if(length(snps) == 1 ){ # or select only snp
          target.snp <- snps
          keep.data <- genotypes[which(as.numeric(loci_info[,2]) == target.snp & loci_info[,1] == chrom),] # isolate genotypes for snp
        }
        
        snp_id <- rownames(genotypes)[which(as.numeric(loci_info[,2]) == target.snp & loci_info[,1] == chrom)]
        gt.chunk <- rbind.data.frame(gt.chunk, keep.data)
        
        # text for troubleshooting
        #print(paste0(i, " ", chrom, " ", target.snp)) 
        #print(paste0(snp_id, " | ", length(ids)))
        #print(dim(gt.chunk))
        
        ids <- c(ids, snp_id)
        colnames(gt.chunk) <- names(keep.data)
        rownames(gt.chunk) <- ids
      } # end if loop
    } # end loop through chunks
    gt.target <- rbind.data.frame(gt.target, gt.chunk) 
  } # end loop through chromosomes
  print(paste0(100*(nrow(gt.target)/nrow(target_regions)), " percent of targeted regions have a sequenced SNP"))
  return(gt.target)
} # end function

process_gt <- function(gt_matrix){
  # take matrix of genotype calls from vcfR::extract.gt and convert it into a dataframe of numeric values
  # 0 and 2 are homozygotes, 1 is heterozygote 
  
  # NA here means SNP not called
  gt <- as.data.frame(gt_matrix, stringsAsFactors = FALSE) %>% as.matrix()
  
  #check which syntax is used to denote genotypes
  genos <- if(any(grepl("/",gt))){
    c("0/0","0/1","1/0","1/1")
  } else if(any(grepl("|",gt))){
    c("0\\|0","0\\|1","1\\|0","1\\|1")
  }
  gt <- gsub(genos[1],0,gt)	
  gt <- gsub(genos[2],1,gt)
  gt <- gsub(genos[3],1,gt)
  gt <- gsub(genos[4],2,gt)
  class(gt) <- "numeric"
  return(gt)
}

filter_snp_missingness <- function(gt_df, miss_threshold, plot = TRUE){ # what fraction of missing data is okay?
  # testing vars
  #  gt_df = gt_samp
  #  miss_threshold = 0.1
  # 
  keep <- vector()
  miss <- vector()
  for(i in 1:nrow(gt_df)){ # for each SNP
    perc_miss <- sum(is.na(gt_df[i,]))/ncol(gt_df) # number of missing genotype calls / number of individuals 
    miss <- c(miss, perc_miss) # add to vector of missingness %s (one per SNP)
    if(perc_miss <= miss_threshold){
      keep <- c(keep, i) 
    } # end if loop
  } # end for loop
  if(plot == TRUE){
    hist(miss, main = "Missingness per SNP")
    abline(v = miss_threshold, col = "red")
  }
  return(gt_df[keep,])
}

filter_ind_missingness <- function(gt_df, miss_threshold, plot = TRUE){
  # testing vars
  #gt_df = gt_PCC_samp_f1
  #miss_threshold = 0.4
  # 
  keep <- vector()
  miss <- vector()
  for(i in 1:ncol(gt_df)){ # for each SNP
    perc_miss <- sum(is.na(gt_df[,i]))/nrow(gt_df) # number of missing genotype calls / number of individuals 
    miss <- c(miss, perc_miss)
    if(perc_miss <= miss_threshold){
      keep <- c(keep, i) 
    } else{
      print(paste0("tossing individual ", colnames(gt_df)[i]))
    }  # end if loop
  } # end for loop
  if(plot == TRUE){
    hist(miss, main = "Missingness per individual")
    abline(v = miss_threshold, col = "red")
  }
  return(gt_df[,keep])
}

plot_ind_missingness <- function(gt_df, miss_threshold, ...){
  # testing vars
  #  gt_df = gt_PCC_samp
  #  miss_threshold = 0.8
  # 
  miss <- vector()
  for(i in 1:ncol(gt_df)){ # for each SNP
    perc_miss <- sum(is.na(gt_df[,i]))/nrow(gt_df) # number of missing genotype calls / number of individuals 
    miss <- c(miss, perc_miss)
  } # end for loop
  hist(miss, xlab = "missingness per individual", ...)
  abline(v = miss_threshold, col = "red")
  mtext(paste0("no. inds: ", dim(gt_df)[2]))
}

plot_snp_missingness <- function(gt_df, miss_threshold, ...){
  # testing vars
  #  gt_df = gt_PCC_samp
  #  miss_threshold = 0.8
  # 
  miss <- vector()
  for(i in 1:nrow(gt_df)){ # for each SNP
    perc_miss <- sum(is.na(gt_df[i,]))/ncol(gt_df) # number of missing genotype calls / number of individuals 
    miss <- c(miss, perc_miss)
  } # end for loop
  hist(miss, xlab = "missingness per SNP", ...)
  abline(v = miss_threshold, col = "red")
  mtext(paste0("no. snps: ", dim(gt_df)[1]))
}

make_SFS <- function(genotypes){
  mac <- vector(length = nrow(genotypes))
  ac <- vector(length = nrow(genotypes))
  
  for(i in 1:nrow(genotypes)){
    # y-axis: number of alleles
    ac[i] <- 2*length(na.omit(genotypes[i,]))
    
    # x-axis: count of minor alleles  
    # count of alt alleles: 
    alt <- sum(genotypes[i,], na.rm = TRUE)
    null <- sum(na.omit(genotypes[i,]) == 0)*2 + sum(na.omit(genotypes[i,]) == 1)
    
    if (alt >= null){ # if there are more "alt" (1) alleles than "null" (0) alleles
      mac[i] <- null
    }
    if (null > alt){ # if there are more "null" (0) alleles than "alt" (1) alleles
      mac[i] <- alt
    }
  }
  SFS_data <- cbind.data.frame(ac, mac)
  return(SFS_data) # dataframe with two columns (count of minor alleles and number of alleles) and one row per SNP. 
}

plot_SFS <- function(genotypes, color = colors[1], rm.zero = FALSE, ylim = c(0, 250), ...){
  sfs <- make_SFS(genotypes)
  if(rm.zero == TRUE){
    sfs <- sfs[which(!sfs$mac == 0), ]
    
    af <- round(sfs$mac/sfs$ac, digits = 2)
    table.sfs <- table(af)
    barplot(table.sfs,
            ylim = ylim,
            ylab = "number of sites", 
            xlab = "minor allele count", 
            col = color, 
            border = NA, 
            density = NULL,
            space = NULL, 
            ...)
  }else{
    af <- round(sfs$mac/sfs$ac, digits = 2)
    table.sfs <- table(af)
    barplot(table.sfs, 
            ylim = ylim,
            ylab = "number of sites", 
            xlab = "minor allele count", 
            col = color, 
            border = NA, 
            density = NULL,
            space = NULL, 
            ...)
  }
}

calc_maf <- function(gt_df){
  maf_vec <- vector(length=ncol(gt_df))
  for(i in 1:ncol(gt_df)){
    # testing
    # gt_df <- gt_flt_ELF
    # count total alleles at a locus
    total <- 2* sum(!is.na(gt_df[i,]))
    
    # count ref allele 
    alt <- sum(gt_df[i,], na.rm=TRUE)
    
    # count alt allele 
    ref <- total - alt
    
    # calculate minor allele frequency 
    
    if(alt > ref){
      maf <- ref / total
    }
    if(ref > alt){
      maf <- alt / total
    }
    maf_vec[i] <- maf
  }
  return(maf_vec)
}

calc_gt_error_rate_flt <- function(dupe_ind_ids, genotypes, tri = TRUE){
  #dupe_ind_ids = dupe_inds
  #genotypes = gt.unflt
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
  error_rate <- vector(length = nrow(genotypes))
  all_same <- vector(length = nrow(genotypes))
  all_pairs <- vector(length = nrow(genotypes))
  for(k in 1:nrow(genotypes)){
    pairs <- 0
    same <- 0
    for(j in 1:nrow(dupe_ind_ids)){
      dupes <- as.character(dupe_ind_ids[j,])
      if(!is.na(genotypes[k,which(colnames(genotypes) == dupes[1])])){
        if(!is.na(genotypes[k,which(colnames(genotypes) == dupes[2])])){ # if there is a gt call for both individuals... 
          pairs = pairs + 1
          if(genotypes[k,which(colnames(genotypes) == dupes[1])] == genotypes[k,which(colnames(genotypes) == dupes[2])]){
            same = same + 1
          }
        }
      }
    }
    error_rate[k] <- (1-(same/pairs))
    all_same[k] <- same
    all_pairs[k] <- pairs
  }
  return(error_rate)
}

target_coverage <- function(loci_info, target_loci = targets){ # need to separate targets and positions.alt (loci_info) bah
  on_target <- vector(length=nrow(target_loci))
  for(chunk in 1:nrow(target_loci)){
    chrom <- loci_info[chunk,1]
    start <- loci_info[chunk,2]
    end <- loci_info[chunk,3]
    snps <- as.numeric(loci_info[which(loci_info[,1] == chrom & as.numeric(loci_info[,2]) > start & as.numeric(loci_info[,2]) < end),2])
    on_target[chunk] <- length(snps)
  }
  return(sum(on_target != 0)/length(on_target))
}

# ------------------------------------------------------------------------------------------------------------


## Read in files ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
## read vcf
vcf <- read.vcfR("../rapture_ref_map_output/emr_rapture_filtered_15.vcf", verbose = T)

## read filtered bed file
targets <- read.csv("../vcf_filtering/Spaced_baits.fas_mrg_buf.bed", sep = "\t", header=FALSE)
colnames(targets) <- c("chrom", "start", "end")

## depth file from vcftools
    # generated using: vcftools --vcf ./rapture_ref_map_output/RE_RUN_725.populations.snps.vcf --site-mean-depth --out ./vcf_filtering/depth_vcftools
depth <- read.csv("/Users/meaghan/Desktop/EMR_rapture/vcf_filtering/depth_vcftools_12924.ldepth.mean", sep = "\t")

# generated using: vcftools --vcf ./rapture_ref_map_output/RE_RUN_725.populations.snps.vcf --depth --out ./vcf_filtering/depth_vcftools
ind_depth <- read.csv("/Users/meaghan/Desktop/EMR_rapture/vcf_filtering/depth_vcftools_12924.idepth", sep = "\t")

site_depth <- read.csv("/Users/meaghan/Desktop/EMR_rapture/vcf_filtering/site_depth.ldepth", sep = "\t")
## Depth Exploration ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

# site  depth------------------------------------------------------------------------------------------------------------
hist(site_depth$SUM_DEPTH)

# Mean individual depth------------------------------------------------------------------------------------------------------------

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

# ------------------------------------------------------------------------------------------------------------

## Verify that all SNPs are biallelic ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

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

# ------------------------------------------------------------------------------------------------------------

## Filter VCF ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

## Extract genotypes
genotypes <- vcfR::extract.gt(vcf, return.alleles = F)

## Create different data subsets based on depth

# all: genotypes, vcf@fix
positions <- vcf@fix

## Transpose and make into data.frame 
gt.unflt <- process_gt(genotypes)

# ------------------------------------------------------------------------------------------------------------

## Test for Excess Het

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

## Generate SFS# ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

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

# target_coverage(loci_info = positions, target_loci = targets) # 0.4697703 there is at least one snp in 52% of targeted regions
# ------------------------------------------------------------------------------------------------------------

## Calculate genotyping error rate
# ------------------------------------------------------------------------------------------------------------
# duplicate individuals for PCC
PCC_dupe_inds <- read.csv("../pedigree_reconstruction/PCC_LifeHistData_dupes.csv", header=FALSE)
dupe_inds <- rbind.data.frame(PCC_dupe_inds, c("ELF_544_P5", "ELF_544_P10", ""))

error_rate <- calc_gt_error_rate_flt(dupe_inds, gt.hetflt, tri = TRUE)
# hist(error_rate)
# abline(v= mean(error_rate), col = "green") # 0.002348661

# ------------------------------------------------------------------------------------------------------------

## Select a single SNP per targeted region # 
# ------------------------------------------------------------------------------------------------------------
## Choose a single SNP randomly per targeted region (sample_snps)

gt_samp <- sample_snps(genotypes = gt.hetflt, target_regions = targets, loci_info = vcf@fix, n_inds = dim(gt.hetflt)[2], z.chrom = TRUE) # "67.2007255139057 percent of targeted regions have a sequenced SNP"

#67.2007255139057 percent of targeted regions have a sequenced SNP

# ------------------------------------------------------------------------------------------------------------

## Remove erroneous individuals from gt_samp
# ------------------------------------------------------------------------------------------------------------
# remove non-PCC individuals and individuals without metadata from PCC
# "PCC_297" no pit tag, not assigned in preliminary pedigree runs, remove
# "PCC_301" 601263359 remove, not from PCC
# "PCC_306" 604078853 not from PCC? 
# added 5/25/2023: 
# "PCC_283", likely same individual as PCC_4
# "PCC_86", likely same individual as 149
# "PCC_99", likely same individual as 210
# "ELF_274", likely contaminated, see binomial tests in filter_vcf_popgen.R
# "ELF_401", likely contaminated, see binomial tests in filter_vcf_popgen.R
# "ELF_138", likely contaminated, see binomial tests in filter_vcf_popgen.R
# "ELF_824", likely contaminated, see binomial tests in filter_vcf_popgen.R
# "ELF_639", likely contaminated, see binomial tests in filter_vcf_popgen.R

gt_samp <- subset(gt_samp, select = -c(PCC_297, PCC_301, PCC_306, PCC_283, PCC_86, PCC_99, ELF_274, ELF_401, ELF_138, ELF_824, ELF_639))
# 2223 1080


## Remove individuals with low mean coverage
gt_samp <- subset(gt_samp, select = -c(ELF_141, ELF_242, ELF_289, ELF_324, ELF_333, ELF_342, ELF_893, ELF_906, ELF_917, PCC_285, PCC_308, PCC_73))

save(gt_samp, 
     file = paste0("../vcf_filtering/unfiltered_gt_", date, ".Robj"))

# ------------------------------------------------------------------------------------------------------------

## FILTER FOR PEDIGREE
# ------------------------------------------------------------------------------------------------------------

## Plot missingness 
# ------------------------------------------------------------------------------------------------------------
pdf(file =paste0("../vcf_filtering/filtered_vcf/pre_filter_missingness_noZ_12924", date, ".pdf"), width = 10, height = 6)
par(mfrow = c(1, 2))
plot_snp_missingness(gt_samp, 0.1, main = "targets, all SNPs", xlim = c(0,0.2))
plot_ind_missingness(gt_samp, 0.2, main = "targets, all SNPs", xlim = c(0,0.8))
dev.off()

# ------------------------------------------------------------------------------------------------------------

# FILTERING AND VIZ
### ------------------------------------------------------------------------------------------------

# filter SNPS then individuals
dim(gt_samp) # 2223 1067
gt_samp_f1 <- filter_snp_missingness(gt_samp, 0.1, plot = TRUE) # drop loci with > 10% missing data
dim(gt_samp_f1) # 2223 1067
gt_samp_f2 <- filter_ind_missingness(gt_samp_f1, 0.2, plot = TRUE) # drop individuals with > 30% missing data
dim(gt_samp_f2) # 2223 1067

# ------------------------------------------------------------------------------------------------------------
## Split into populations
# ------------------------------------------------------------------------------------------------------------
raw_names <- unlist(strsplit(colnames(gt_samp_f2), "_"))
sites <- raw_names[which(raw_names == "ELF" | raw_names == "PCC")]
gt_flt_ELF <- gt_samp_f2[,which(sites == "ELF")] # 2468  784 # 2223  780
gt_flt_PCC <- gt_samp_f2[,which(sites == "PCC")] # 2468  289 # 2223  288
# ------------------------------------------------------------------------------------------------------------

## Calculate MAF
# ------------------------------------------------------------------------------------------------------------
maf <- calc_maf(gt_flt_PCC)
hist(maf)
sum(maf > 0.1)/length(maf) #  0.6597222

# maf ELF
maf <- calc_maf(gt_flt_ELF)
hist(maf)
sum(maf > 0.1)/length(maf) # 0.6923077

## Save filtered data ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
# PCC
PCC_filt_gt_all <- gt_flt_PCC
save(PCC_filt_gt_all, file = paste0("../vcf_filtering/PCC_filtered_gt_all_snps_noZ_", date, ".Robj"))

# ELF
ELF_filt_gt_all <- gt_flt_ELF
save(ELF_filt_gt_all, file = paste0("../vcf_filtering/ELF_filtered_gt_all_snps_noZ_", date, ".Robj"))

# ------------------------------------------------------------------------------------------------------------


## Graveyard ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
# this code was put into the graveyard when I switched from separating sites before selecting a single snp per targeted
# region to selecting the same set of snps for each site so that I can filter based on genotyping error rate in make_prelim_pedigree.R
# ------------------------------------------------------------------------------------------------------------
## split file into ELF and PCC --------------------------------------------------------------------------------
raw_names <- unlist(strsplit(colnames(gt.unflt), "_"))
sites <- raw_names[which(raw_names == "ELF" | raw_names == "PCC")]

### do not split yet! 
# ELF
gt_ELF <- gt.unflt[,which(sites == "ELF")]
gt_ELF_low <- gt[,which(sites == "ELF")]
gt_ELF_high <- gt.alt[,which(sites == "ELF")]
# PCC
gt_PCC <- gt.unflt[,which(sites == "PCC")]
gt_PCC_low <- gt[,which(sites == "PCC")]
gt_PCC_high <- gt.alt[,which(sites == "PCC")]


save(gt_ELF_samp, gt_ELF_samp_low, gt_ELF_samp_high, 
              gt_PCC_samp, gt_PCC_samp_low, gt_PCC_samp_high, 
              file = paste0("../vcf_filtering/unfiltered_gt", date, ".Robj"))

## plot missingness 
#PCC ------------------------------------------------
pdf(file =paste0("../vcf_filtering/filtered_vcf/pre_filter_missingness_PCC_noZ", date, ".pdf"), width = 10, height = 6)
par(mfrow = c(3, 4))
# all snps
plot_snp_missingness(gt_PCC_samp, 0.1, main = "targets, all SNPs", xlim = c(0,0.8))
plot_ind_missingness(gt_PCC_samp, 0.5, main = "targets, all SNPs", xlim = c(0,0.8))

plot_snp_missingness(gt_PCC, 0.1, main = "unfiltered, all SNPs", xlim = c(0,0.8))
plot_ind_missingness(gt_PCC, 0.5, main = "unfiltered, all SNPs", xlim = c(0,0.8))

# high
plot_snp_missingness(gt_PCC_samp_high, 0.1, main = "targets, high", xlim = c(0,0.8))
plot_ind_missingness(gt_PCC_samp_high, 0.5, main = "targets, high", xlim = c(0,0.8))

plot_snp_missingness(gt_PCC_high, 0.1, main = "unfiltered, high", xlim = c(0,0.8))
plot_ind_missingness(gt_PCC_high, 0.5, main = "unfiltered, high", xlim = c(0,0.8))

# low
plot_snp_missingness(gt_PCC_samp_low, 0.1, main = "targets, low", xlim = c(0,0.8))
plot_ind_missingness(gt_PCC_samp_low, 0.5, main = "targets, low", xlim = c(0,0.8))

plot_snp_missingness(gt_PCC_low, 0.1, main = "unfiltered, low", xlim = c(0,0.8))
plot_ind_missingness(gt_PCC_low, 0.5, main = "unfiltered, low", xlim = c(0,0.8))

dev.off()
### ------------------------------------------------------------------------------------------------
#ELF ------------------------------------------------
pdf(file =paste0("../vcf_filtering/filtered_vcf/pre_filter_missingness_ELF_noZ", date, ".pdf"), width = 10, height = 6)
par(mfrow = c(3, 4))
# all snps
plot_snp_missingness(gt_ELF_samp, 0.1, main = "targets, all SNPs", xlim = c(0,0.8))
plot_ind_missingness(gt_ELF_samp, 0.5, main = "targets, all SNPs", xlim = c(0,0.8))

plot_snp_missingness(gt_ELF, 0.1, main = "unfiltered, all SNPs", xlim = c(0,0.8))
plot_ind_missingness(gt_ELF, 0.5, main = "unfiltered, all SNPs", xlim = c(0,0.8))

# high
plot_snp_missingness(gt_ELF_samp_high, 0.1, main = "targets, high", xlim = c(0,0.8))
plot_ind_missingness(gt_ELF_samp_high, 0.5, main = "targets, high", xlim = c(0,0.8))

plot_snp_missingness(gt_ELF_high, 0.1, main = "unfiltered, high", xlim = c(0,0.8))
plot_ind_missingness(gt_ELF_high, 0.5, main = "unfiltered, high", xlim = c(0,0.8))

# low
plot_snp_missingness(gt_ELF_samp_low, 0.1, main = "targets, low", xlim = c(0,0.8))
plot_ind_missingness(gt_ELF_samp_low, 0.5, main = "targets, low", xlim = c(0,0.8))

plot_snp_missingness(gt_ELF_low, 0.1, main = "unfiltered, low", xlim = c(0,0.8))
plot_ind_missingness(gt_ELF_low, 0.5, main = "unfiltered, low", xlim = c(0,0.8))

dev.off()


# Recommendation call rate filtering from Sequioa 
par(mfrow=c(2,2))
plot_snp_missingness(gt_ELF_samp, 0.5, main = "targets, all SNPs", xlim = c(0,0.8))
plot_ind_missingness(gt_ELF_samp, 0.2, main = "targets, all SNPs", xlim = c(0,0.8))

plot_snp_missingness(gt_PCC_samp, 0.5, main = "targets, all SNPs", xlim = c(0,0.8))
plot_ind_missingness(gt_PCC_samp, 0.2, main = "targets, all SNPs", xlim = c(0,0.8))

### ------------------------------------------------------------------------------------------------

## Filter for missingness 

# FILTERING AND VIZ
### ------------------------------------------------------------------------------------------------
# prelim pass: used filter_snp_missingness with a 0.1 threshol
# changing to filtering snps with 0.1 threshold and filter inds with 0.2 threshold
#PCC

gt_PCC_samp_f1 <- filter_ind_missingness(gt_PCC_samp, 0.2, plot = TRUE) # 
gt_PCC_samp_low_f1 <- filter_ind_missingness(gt_PCC_samp_low, 0.2, plot = TRUE) # 
gt_PCC_samp_high_f1 <- filter_ind_missingness(gt_PCC_samp_high, 0.2, plot = TRUE) # 

dim(gt_PCC_samp_f1)
dim(gt_PCC_samp_low_f1)
dim(gt_PCC_samp_high_f1)

gt_PCC_samp_f2 <- filter_snp_missingness(gt_PCC_samp_f1, 0.1, plot = TRUE) # drop loci with > 10% missing data
gt_PCC_samp_low_f2 <- filter_snp_missingness(gt_PCC_samp_low_f1, 0.1, plot = TRUE) # drop loci with > 10% missing data
gt_PCC_samp_high_f2 <- filter_snp_missingness(gt_PCC_samp_high_f1, 0.1, plot = TRUE) # drop loci with > 10% missing data

dim(gt_PCC_samp_f2)
dim(gt_PCC_samp_low_f2)
dim(gt_PCC_samp_high_f2)

#ELF
gt_ELF_samp_f1 <- filter_ind_missingness(gt_ELF_samp, 0.2, plot = TRUE) # drop loci with > 10% missing data
gt_ELF_samp_low_f1 <- filter_ind_missingness(gt_ELF_samp_low, 0.2, plot = TRUE) # drop loci with > 10% missing data
gt_ELF_samp_high_f1 <- filter_ind_missingness(gt_ELF_samp_high, 0.2, plot = TRUE) # drop loci with > 10% missing data

dim(gt_ELF_samp_f1)
dim(gt_ELF_samp_low_f1)
dim(gt_ELF_samp_high_f1)

gt_ELF_samp_f2 <- filter_snp_missingness(gt_ELF_samp_f1, 0.1, plot = TRUE) # drop loci with > 10% missing data
gt_ELF_samp_low_f2 <- filter_snp_missingness(gt_ELF_samp_low_f1, 0.1, plot = TRUE) # drop loci with > 10% missing data
gt_ELF_samp_high_f2 <- filter_snp_missingness(gt_ELF_samp_high_f1, 0.1, plot = TRUE) # drop loci with > 10% missing data

dim(gt_ELF_samp_f2) # 2368  788
dim(gt_ELF_samp_low_f2)
dim(gt_ELF_samp_high_f2)

# pdf("../vcf_filtering/stepwise_filters_PCC.pdf", width = 6, height = 12)
# par(mfrow=c(5,2))
# # pre=filtering
# plot_snp_missingness(gt_PCC_samp, 1.0, main = "pre-filtering", xlim = c(0,0.8))
# plot_ind_missingness(gt_PCC_samp, 1.0, main = "pre-filtering", xlim = c(0,0.8))
# 
# # drop loci with > 80% missing data
# plot_snp_missingness(gt_PCC_samp_f1, 0.8, main = "snps >= 80% missing dropped", xlim = c(0,0.8))
# plot_ind_missingness(gt_PCC_samp_f1, 0.8, main = "snps >= 80% missing dropped", xlim = c(0,0.8))
# 
# # drop individuals with > 50% missing data
# plot_snp_missingness(gt_PCC_samp_f2, 0.5, main = "inds >= 50% missing dropped", xlim = c(0,0.8))
# plot_ind_missingness(gt_PCC_samp_f2, 0.5, main = "inds >= 50% missing dropped", xlim = c(0,0.8))
# 
# # drop loci with > 50% missing data
# plot_snp_missingness(gt_PCC_samp_f3, 0.5, main = "snps >= 50% missing dropped", xlim = c(0,0.8))
# plot_ind_missingness(gt_PCC_samp_f3, 0.5, main = "snps >= 50% missing dropped", xlim = c(0,0.8))
# 
# # drop loci with > 10% missing data
# plot_snp_missingness(gt_PCC_samp_f4, 0.1, main = "inds >= 10% missing dropped", xlim = c(0,0.8))
# plot_ind_missingness(gt_PCC_samp_f4, 0.1, main = "snps >= 10% missing dropped", xlim = c(0,0.8))
# dev.off()
# 
# 
# dim(gt_PCC_samp_f4)
# 
# # retain 297 individuals and 2566 loci from PCC 
# 
#   # less dropped individuals when filter for depth FIRST 
# 
# 
### ------------------------------------------------------------------------------------------------

# maf PCC
maf <- calc_maf(gt_PCC_samp_f1)
hist(maf)
sum(maf > 0.1)/length(maf)

# maf ELF
maf <- calc_maf(gt_ELF_samp_f1)
hist(maf)
sum(maf > 0.1)/length(maf)

## Save filtered data ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
# PCC
PCC_filt_gt_all <- gt_PCC_samp_f2
save(PCC_filt_gt_all, file = paste0("../vcf_filtering/PCC_filtered_gt_all_snps_noZ", date, ".Robj"))

PCC_filt_gt_low <- gt_PCC_samp_low_f2
save(PCC_filt_gt_low, file = paste0("../vcf_filtering/PCC_filtered_gt_low_snps_noZ", date, ".Robj"))

PCC_filt_gt_high <- gt_PCC_samp_high_f2
save(PCC_filt_gt_high, file = paste0("../vcf_filtering/PCC_filtered_gt_high_snps_noZ", date, ".Robj"))

# ELF
ELF_filt_gt_all <- gt_ELF_samp_f2
save(ELF_filt_gt_all, file = paste0("../vcf_filtering/ELF_filtered_gt_all_snps_noZ", date, ".Robj"))

ELF_filt_gt_low <- gt_ELF_samp_low_f2
save(ELF_filt_gt_low, file = paste0("../vcf_filtering/ELF_filtered_gt_low_snps_noZ", date, ".Robj"))

ELF_filt_gt_high <- gt_ELF_samp_high_f2
save(ELF_filt_gt_high, file = paste0("../vcf_filtering/ELF_filtered_gt_high_snps_noZ", date, ".Robj"))


# ------------------------------------------------------------------------------------------------------------

gt.wide <- as.data.frame(gt)
gt.wide$sampid <- rownames(gt.wide)

pops <- unlist(strsplit(gt.wide$sampid, "_"))[c(TRUE,FALSE)]
gt.PCC <- gt.wide[which(pops == "PCC"),]
gt.ELF <- gt.wide[which(pops == "ELF"),]

# compare duplicate individuals
dupes <- gt.wide[which(grepl("_P", gt.wide$sampid)),]
rownames(dupes)
inds <- vector(length = length(rownames(dupes)))
for (i in 1:length(rownames(dupes))){
  
  inds[i] <- paste0(strsplit(rownames(dupes)[i], "_")[[1]][1], "_", strsplit(rownames(dupes)[i], "_")[[1]][2])
}
u.inds <- unique(inds)
dupes$sampid <- inds

not.equal <- vector(length = length(u.inds))
equal <- vector(length = length(u.inds))
nas <- vector(length = length(u.inds))
for (i in 1:length(u.inds)){
  ind <- u.inds[i]
  comp.rows <- dupes[which(dupes$sampid==ind),]
  #comp.rows <- rbind.data.frame(c(1:10), c(1, 2, 3, 4, 10, 9, 8, 7, 6, 5))
  count.false = 0
  count.true = 0
  count.na = 0 
  for (j in 1:ncol(comp.rows)) {
    if (!is.na(comp.rows[1,j]) & !is.na(comp.rows[2,j])){
      if (comp.rows[1,j] != comp.rows[2,j]){
        count.false = count.false + 1
      }else(count.true = count.true + 1)
    }else(count.na = count.na + 1)
  }
  not.equal[i] <- count.false
  equal[i] <- count.true
  nas[i] <- count.na
}

dupes_percent <- cbind.data.frame(u.inds, 100*(not.equal/(not.equal + equal+ nas)))
colnames(dupes_percent) <- c("individual", "percent different")

# ------------------------------------------------------------------------------------------------------------
## not sure I actually want to work with long format data
#convert genotype matrix from wide to long format
gt.long <- gt %>% as.data.frame()
gt.long$sampid = rownames(gt.long) 
gt.long <- gt.long %>% dplyr::select(sampid,everything())
gt.long <- gt.long %>% tidyfast::dt_pivot_longer(.,names_to="SNPid",values_to="genotype",cols=2:ncol(gt.long), factor_key = T)

#merge on population key

#get pop names
gt.long <- gt.long %>% tidyfast::dt_separate(., sampid, into = c("pop","garbage", "dup"), sep = "_", remove = F)

#get a dataframe of genotype freqs per SNP within each population (removing NAs here removes individuals that didn't have this SNP scored)
freq.perSNP <- gt.long %>% group_by(pop,SNPid,genotype) %>% summarise(N=n()) %>% ungroup() %>% na.omit()
freq.perSNP <- freq.perSNP %>% tidyfast::dt_pivot_wider(., names_from = genotype, values_from = N) %>% as.data.frame()
#replace NA's with 0's so math works later (NA here means SNP called but 0 of that genotype present)
freq.perSNP[is.na(freq.perSNP)==T] <- 0
#make column N - how many indivs is SNP scored/sequenced in?
freq.perSNP <- freq.perSNP %>% mutate(N=`0`+`1`+`2`)
#calculate allele counts from genotype frequencies
freq.perSNP <- freq.perSNP %>% mutate(count_A = (`1` + `0`*2),
                                      count_a = (`1` + `2`*2),
                                      #count_minor = pmin(count_A, count_a),
                                      propindivshet = `1`/N,
                                      freq_A = count_A/(2*N),
                                      freq_a = count_a/(2*N),
                                      He = 2*freq_A*freq_a,
                                      FIS = ((He-propindivshet)/He))




#find private alleles (and those only scored in 1 pop)
private <- freq.perSNP %>% dplyr::select(pop,SNPid,N,count_A,count_a) %>% 
  #count how many pops each SNP is scored in (so we can ID SNPs only scored in 1 pop)
  group_by(SNPid) %>% add_tally(name = "total_pops_scored_in") %>% 
  #to be a private allele (across all pops) only 1 pop can have an allele count > 0
  #so count how many pops have allele count > 0
  mutate(N_pops_with_A = sum(count_A > 0)) %>% 
  mutate(N_pops_with_a = sum(count_a > 0)) %>% 
  #mark which pop private allele is in
  mutate(SNPstatus = ifelse(count_A > 0 & N_pops_with_A == 1, "private", NA)) %>% 
  mutate(SNPstatus = ifelse(count_a > 0 & N_pops_with_a == 1, "private", SNPstatus)) %>% 
  #denote SNPs only scored in 1 pop
  mutate(SNPstatus = ifelse(N_pops_with_A == 1 & N_pops_with_a == 1 & total_pops_scored_in == 1, "only_scored_in_this_pop", SNPstatus))

#count how many are in each pop
pri.all <- private %>% group_by(pop,SNPstatus) %>% summarise(n=n()) %>% filter(is.na(SNPstatus)==F) %>% 
  pivot_wider(., names_from = "SNPstatus", values_from = "n") %>% merge(., freq.perSNP %>% distinct(pop), by = "pop", all.x = T, all.y = T)

# note - NA in above output means 0 private alleles

pri.all$pop <- as.factor(pri.all$pop)
barplot(pri.all$private~pri.all$pop)


