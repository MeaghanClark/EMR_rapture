#!/bin/bash

#--------------- EXECUTABLE ---------------

# for RAD data with a reference genome.
#this script aligns pairs of forward and reverse sequence files listed in a text file to the given reference
		# Last updated 03/16/2022 by MI Clark, originally written by MI Clark and R Toczydlowski based on script by B Reid
# input: Reference genome path, forward and reverse reads from wrapper-align_to_genome.sh 
#
# output: ${SAMPLE_NAME}.sort.flt1.bam file, forward and reverse reads aligned to reference genome
#			final .bam files are stored in directory at <storagenode>/<run_name>/alignments 
#			example: /mnt/home/clarkm89/massasauga/alignments/PCC_321.rmdup.bam 

#load programs we want to use
module purge
module load GCC/6.4.0-2.28
module load OpenMPI/2.1.1
module load SAMtools
module load BWA/0.7.17
module load powertools
module list 
 
# text for troubleshooting
echo I am importing: 
echo reference as $REFERENCE
echo forwardread as $FORWARDREAD
echo reverseread as $REVERSEREAD
echo cpus as $CPUS
echo run name as $RUN_NAME
echo sample name as $SAMPLE_NAME
echo storage node as $STORAGENODE
echo scratch node as $SCRATCHNODE
echo indir as $INDIR
echo outdir as $OUTDIR
echo logfile dir as $LOGFILESDIR


echo using $REFERENCE 
echo trying to locate files ${INDIR}/$FORWARDREAD and ${INDIR}/$REVERSEREAD



#run bwa 
bwa mem -t $CPUS $REFERENCE ${INDIR}/$FORWARDREAD ${INDIR}/$REVERSEREAD > ${SCRATCHNODE}/${SAMPLE_NAME}.aln-pe.sam 
 
echo I have created a .sam file for $SAMPLE_NAME and put it in $SCRATCHNODE

# run samtools
samtools view -Sb -o ${SCRATCHNODE}/${SAMPLE_NAME}.aln-pe.bam ${SCRATCHNODE}/${SAMPLE_NAME}.aln-pe.sam #change sam to bam
echo I have converted the .sam file to a .bam file and saved it in ${SCRATCHNODE}

samtools sort -o ${SCRATCHNODE}/${SAMPLE_NAME}.sort.bam ${SCRATCHNODE}/${SAMPLE_NAME}.aln-pe.bam #sort by readname
echo I have sorted the .bam file and saved it as ${SCRATCHNODE}/${SAMPLE_NAME}.sort.bam

samtools view -f 0x2 -b ${SCRATCHNODE}/${SAMPLE_NAME}.sort.bam > ${OUTDIR}/${SAMPLE_NAME}.sort.flt1.bam # filter for proper pairs of reads
echo I have filtered the .bam file and saved it as ${SCRATCHNODE}/${SAMPLE_NAME}.sort.flt1.bam

echo I am all done with using bwa and samtools to align $FORWARDREAD and $REVERSEREAD to $REFERENCE

#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)

echo ----------------------------------------------------------------------------------------
seff ${SLURM_JOBID}
