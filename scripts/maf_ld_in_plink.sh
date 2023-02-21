#!/bin/bash -login

# script written by M. Clark, 2/21/2023
# filters vcf for maf cut off and calculate pairwise linkage disequilibrium 
# run from directory with vcf file to be filtered. Outputs will be put in the same directory

# usage: 
# maf_ld_in_plink INPUT.vcf MAF_OUTPUT LD_OUTPUT

# load plink
module load PLINK/1.9b_4.1-x86_64

INPUT=$1
OUTPUT_NAME=$2


## load .vcf, do maf cut off, output vcf

plink2 --vcf $INPUT --maf 0.05 --recode vcf --out $MAF_OUTPUT

## load maf filtered .vcf and calculate pairwise ld 

# Load the MAF filtered VCF file and turn into plink binary output (bed file)
plink2 --vcf ${MAF_OUTPUT}.vcf --make-bed --out $MAF_OUTPUT

# Calculate pairwise linkage disequilibrium
plink2 --bfile $MAF_OUTPUT --r2 --out $LD_OUTPUT # maybe not matrix?


#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)

echo ----------------------------------------------------------------------------------------
