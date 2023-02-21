#!/bin/bash -login

# script written by M. Clark, 8/2/2022
# filters vcf output by STACKs for downstream pedigree analysis 

# usage: 
# filter_vcf.sh INPUT.vcf OUTPUT.vcf 

#load programs we want to use
module purge
module load GCC/6.4.0-2.28
module load OpenMPI/2.1.2
module load bcftools/1.9.64
module list

INPUT=$1
OUTPUT=$2

# filter to retain SNPs where 90% of individuals have > 7 reads
# filter to retain SNPs where 90% of individuals have a GQ > 19

bcftools view -i 'F_PASS(FMT/DP>7) > 0.9' $INPUT |
bcftools view -i 'F_PASS(FMT/GQ>19) > 0.9' -o $OUTPUT




#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)

echo ----------------------------------------------------------------------------------------
