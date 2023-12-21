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

echo
echo input file is $INPUT, maf output is $MAF_OUTPUT, and ld output is $LD_OUTPUT
echo

echo ----------------------------------------------------------------------------------------
echo load .vcf, do maf cut off, output vcf

# plink --vcf $INPUT --allow-extra-chr --double-id --maf 0.05 --recode vcf --out $MAF_OUTPUT

echo 

## load maf filtered .vcf and calculate pairwise ld 

# Load the MAF filtered VCF file and turn into plink binary output (bed file)
plink --vcf ${MAF_OUTPUT}.vcf --allow-extra-chr --double-id --make-bed --out $MAF_OUTPUT

echo
echo ----------------------------------------------------------------------------------------
echo Calculate pairwise linkage disequilibrium

# Calculate pairwise linkage disequilibrium
# plink --bfile $MAF_OUTPUT --allow-extra-chr --double-id --r2  --inter-chr --ld-window-r2 0 --out $LD_OUTPUT # the inter-chr option specifies to calculate r2 for all snp pairs. This could take a while! 
# from these results, I want to use a window size of 900 bp 
echo
echo ----------------------------------------------------------------------------------------
echo Calculate pairwise LD between SNPs in windows and then calculate physical distance between neighboring SNPs

# calculate pairwise LD between SNPs in windows of 50 kb with a step size of 5 Kb
plink --bfile $MAF_OUTPUT --allow-extra-chr --double-id --indep-pairwise 900 100 0.2 --out ${LD_OUTPUT} # unclear if this is the right setting?

# filter VCF based on SNP list from above command
plink --bfile $MAF_OUTPUT --allow-extra-chr --double-id --extract ${LD_OUTPUT}.prune.in --recode vcf --out ${LD_OUTPUT}_pruned

echo
echo ----------------------------------------------------------------------------------------
echo ----------------------------------------------------------------------------------------
