# get Ne estimates

# install strataG for LDNe in R capabilities 
# I remember this install being somewhat of a pain, but it eventually worked! 

load(strataG)
load(adegenet) # used to convert from allele counts to strataG input object

# load data with estimated birth years
load("../inbreeding_models/data_for_analyses_07132024.Robj")

# genotype data is called "t_gt_PCC" and "t_gt_ELF". These are dataframes of allele counts where
# each row is an individual and each column is a SNP, and the value is the count of the minor 
# allle (i.e. all 0s, 1s, and 2s)

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
# looks like I tried doing this w/ and w/0 a maf threshold, and I don't think it made a difference

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
