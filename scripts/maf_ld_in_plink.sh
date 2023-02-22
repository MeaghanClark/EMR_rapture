#!/bin/bash -login

# script written by M. Clark, 2/21/2023
# filters vcf for maf cut off and calculate pairwise linkage disequilibrium 
# run from directory with vcf file to be filtered. Outputs will be put in the same directory

# usage: 
# maf_ld_in_plink INPUT.vcf MAF_OUTPUT LD_OUTPUT

# load plink
module load PLINK/1.9b_4.1-x86_64

INPUT=$1
MAF_OUTPUT=$2
LD_OUTPUT=$3


## load .vcf, do maf cut off, output vcf

plink --vcf $INPUT --allow-extra-chr --double-id --maf 0.05 --recode vcf --out $MAF_OUTPUT

echo ----------------------------------------------------------------------------------------

## load maf filtered .vcf and calculate pairwise ld 

# Load the MAF filtered VCF file and turn into plink binary output (bed file)
plink --vcf ${MAF_OUTPUT}.vcf --allow-extra-chr --double-id --make-bed --out $MAF_OUTPUT

echo ----------------------------------------------------------------------------------------

# Calculate pairwise linkage disequilibrium
plink --bfile $MAF_OUTPUT --allow-extra-chr --double-id --r2 --out $LD_OUTPUT # maybe not matrix?


echo ----------------------------------------------------------------------------------------
# calculate pairwise LD between SNPs in windows of 50 kb with a step size of 5 Kb
plink --bfile $MAF_OUTPUT --indep-pairwise 50 5

# Calculate physical distance between neighboring SNPs
plink --bfile $MAF_OUTPUT --extract plink.prune.in --make-bed --out ${LD_OUTPUT}_pruned
awk '{print $1,$4,$4+$3-$2}' ${LD_OUTPUT}_pruned.bim > ${LD_OUTPUT}_pruned.dist


echo ----------------------------------------------------------------------------------------
echo ----------------------------------------------------------------------------------------
