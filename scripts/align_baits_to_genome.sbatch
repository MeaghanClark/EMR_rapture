#!/bin/bash

#--------------- EXECUTABLE ---------------

# This script is designed to align fasta files containing bait information to a reference genome and generate .bed files corresponding to the targeted loci
#this script aligns sequence files listed in a text file to the given reference
		# Last updated 03/21/2022 by MI Clark, originally written by MI Clark and R Toczydlowski based on script by B Reid
# input: Reference genome path, fasta files, in and out directories from wrapper-align_to_genome.sh 
# output:	${OUTDIR}/${BAITS}.aln.sam		# alignment
#		${OUTDIR}/${BAITS}.aln.bam		# converted to .bam
#		${OUTDIR}/${BAITS}.sort.bam		# sorted
#		${OUTDIR}/${BAITS}.bed			# convert .bam to .bed
#		${OUTDIR}/${BAITS}_merged.bed		# merge overlapping/adjacent baits
#		${OUTDIR}/${BAITS}_mrg_buf.bed		# add 500 bp buffer



#load modules for alignment, sorting, and indexing
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
echo baits as $BAITS
echo run name as $RUN_NAME
echo storage node as $STORAGENODE
echo indir as $INDIR
echo outdir as $OUTDIR
echo logfile dir as $LOGFILESDIR

echo using $REFERENCE 

#run bwa to align baits to reference genome
bwa mem $REFERENCE ${INDIR}/${BAITS} > ${OUTDIR}/${BAITS}.aln.sam 
 
echo I have created a .sam file for $BAITS and put it in $OUTDIR

# run samtools to sort and index alignment
samtools view -Sb -o ${OUTDIR}/${BAITS}.aln.bam ${OUTDIR}/${BAITS}.aln.sam #change sam to bam
echo I have converted the .sam file to a .bam file and saved it in ${OUTDIR}

samtools sort -o ${OUTDIR}/${BAITS}.sort.bam ${OUTDIR}/${BAITS}.aln.bam #sort by readname
echo I have sorted the .bam file and saved it as ${OUTDIR}/${BAITS}.sort.bam

samtools index ${OUTDIR}/${BAITS}.sort.bam
echo I have indexed the bam file. 

# load modules for making .bed files
module purge
module load GCC/10.2.0
module load BEDTools/2.30.0
module load SAMtools
module list

# create .bed file and merge all overlapping/adjacent baits

bedtools bamtobed -i ${OUTDIR}/${BAITS}.sort.bam > ${OUTDIR}/${BAITS}.bed
bedtools merge -i ${OUTDIR}/${BAITS}.bed > ${OUTDIR}/${BAITS}_merged.bed

# add a +500 bp buffer
samtools faidx $REFERENCE
awk -v OFS='\t' {'print $1,$2'} ${REFERENCE}.fai > ${OUTDIR}/genomeFile.txt
slopBed -i ${OUTDIR}/${BAITS}_merged.bed -g ${OUTDIR}/genomeFile.txt -r 500 -l 0 > ${OUTDIR}/${BAITS}_mrg_buf.bed

#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)

echo ----------------------------------------------------------------------------------------
seff ${SLURM_JOBID}
