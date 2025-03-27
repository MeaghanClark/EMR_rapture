## load libraries
library(vcfR)
library(dplyr)
library(tidyr)
library(tidyfast)
library("MetBrewer")
library(HardyWeinberg)

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

binom_test_for_allele_counts <- function(ind_gt, ind){
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
  hist(p_val, main = ind)
  if(exists("pop_gen_stats")){
    mtext(paste0("het = ", round(pop_gen_stats[which(pop_gen_stats$id == ind),"het"], digits = 3)), col = "red")
  }
  
  #return(p_val)
}

