## File name: vcf_depth.R
## Purpose: Load VCF file from STACKS and calculate depth per locus per individual for EMR project
## M. I. Clark, July 2022
## Last updated: 07/28/2022

## Set up ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
Sys.getenv("R_LIBS_USER")

# capture command line variables
args <- commandArgs(trailingOnly = TRUE)
print(c("My args are ", args))

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

INDIR <- args[1]
VCF <- args[2]
OUTDIR <- args[3]

## loading ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
## load libraries
library(vcfR)

# ------------------------------------------------------------------------------------------------------------

## read in files ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
## read vcf
vcf <- read.vcfR(paste0(INDIR, "/", VCF), verbose = T)

# ------------------------------------------------------------------------------------------------------------

## Depth------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

data <- vcf@gt
depth <- as.data.frame(matrix(data=NA, ncol= ncol(data), nrow=nrow(data)))
for (i in 2:ncol(data)){
  for (j in 1:nrow(data)){
    depth[j,i] <-  as.numeric(strsplit(data[j,i], ":")[[1]][2])
  }
}

save(depth, file = paste0(OUTDIR, "depth.Robj"))
