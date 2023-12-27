#!/bin/bash

# This executable will run cutadapt on paired end sequencing data
# Intended to be called by wrapper-run_cutadapt.sh
# M. Clark, last updated 12/27/2023
module purge
module load GCC/7.3.0-2.30  OpenMPI/3.1.1
module load cutadapt/2.1-Python-3.6.6
module list

--export=FORWARDREAD=$forwardread,REVERSEREAD=$reverseread,CPUS=$cpus,RUN_NAME=$run_name,SAMPLE_NAME=$sample_name,STORAGENODE=$storagenode,SCRATCHNODE=$scratchnode,INDIR=$indir,OUTDIR=$outdir,LOGFILESDIR=$logfilesdir \

echo I am importing: $FORWARDREAD, $REVERSEREAD, $CPUS, $RUN_NAME, $SAMPLE_NAME, $STORAGENODE, $SCRATCHNODE, $INDIR, $OUTDIR, $LOGFILESDIR \

/mnt/research/Fitz_Lab/software/cutadapt/bin/cutadapt \
	-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    --report=full    \
    -o trimmed.${FORWARDREAD} -p trimmed.${REVERSEREAD} \
	$FORWARDREAD $REVERSEREAD
	
wait

echo I am done

#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)

echo ----------------------------------------------------------------------------------------	
