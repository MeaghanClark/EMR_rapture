#!/bin/bash

#--------------- EXECUTABLE ---------------

# for RAPTURE data.
#this script calcuates coverage at targeted loci specified in a .bed file for an alignment file
		# Last updated 03/23/2022 by MI Clark, originally written by MI Clark and R Toczydlowski based on script by B Reid

#load programs we want to use
module purge
module load GCC/10.2.0
module load BEDTools/2.30.0
module load SAMtools
module load powertools
module list


# text for troubleshooting
echo I am importing: 
echo bed as $BED
echo cpus as $CPUS
echo alignment file as $alignment
echo run name as $RUN_NAME
echo sample name as $SAMPLE_NAME
echo storage node as $STORAGENODE
echo indir as $INDIR
echo outdir as $OUTDIR
echo logfile dir as $LOGFILESDIR


# calculate coverage using bedtools
bedtools coverage -a $BED -b ${INDIR}/${SAMPLE_NAME}.sort.flt1.bam > ${OUTDIR}/${SAMPLE_NAME}.cov


#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)

echo ----------------------------------------------------------------------------------------
seff ${SLURM_JOBID}
