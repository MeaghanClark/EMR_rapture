#!/bin/bash -login

#load programs we want to use
module purge
module load GCC/9.3.0  
module load OpenMPI/4.0.3
module load iccifort/2020.1.217
module load impi/2019.7.217
module load Stacks/2.59

echo I am importing: $JOBNAME, $CPUS, $RUN_NAME, $STORAGENODE, $SCRATCHNODE, $LOGFILESDIR, $OUTDIR, $OUTDIR_CLONE, $FORWARDREAD, $REVERSEREAD, $BARCODES \

process_radtags -1 ${INDIR}/${FORWARDREAD} -2 ${INDIR}/${REVERSEREAD} --paired -o ${OUTDIR} -i gzfastq -y gzfastq -b $BARCODES -q -r --bestrad --barcode_dist_1 1 -c -D -e sbfI --inline_null

wait 

echo I am done demultiplexing

#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)

echo ----------------------------------------------------------------------------------------
