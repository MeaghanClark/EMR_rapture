#!/bin/bash


#load programs we want to use
# CHANGE
module purge
module load GCC/8.3.0
module load OpenMPI/3.1.4
module load  R/4.0.2


module list

export R_LIBS_USER=~/R_Lib/4.0.2

Rscript ${RSCRIPT} --args ${INDIR} ${VCR} ${OUTDIR}

echo done running rscript

#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)

echo ----------------------------------------------------------------------------------------

