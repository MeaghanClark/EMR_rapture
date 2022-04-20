#!/bin/bash -login

#load programs we want to use
module purge
module load GCC/9.3.0  
module load OpenMPI/4.0.3
module load iccifort/2020.1.217
module load impi/2019.7.217
module load Stacks/2.59

echo I am importing: $JOBNAME, $CPUS, $RUN_NAME, $STORAGENODE, $LOGFILESDIR, $INDIR, $OUTDIR, $POPMAP, $M, $N \

ref_map.pl --samples $INDIR --popmap $POPMAP -o $OUTDIR -T $CPUS --rm-pcr-duplicates -X "populations: -p 2 -r 0.75 --hwe --fstats --vcf --structure --phylip"

echo I have finished running ref_map.pl


#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)

echo ----------------------------------------------------------------------------------------
