#!/bin/bash
		
# Last updated 03/22/2022 by MI Clark, originally written by R Toczydlowski 

# run from project directory (where you want output directory to be created)
# Submits jobs (one per alignment file) that filters for loci specified in a provided .bed file 

# usage: ./scripts/wrapper-align_to_genome.sh [site]

#define variables:
storagenode=/mnt/home/clarkm89 #path to top level of dir where input/output files live

logfilesdir=logfiles_iso_rap #name of directory to create and then write log files to
executable=./scripts/isolate_rapture_loci.sbatch #script to run 

cpus=1 #number of CPUs to request/use per dataset
<<<<<<< HEAD
ram_per_cpu=4G #amount of RAM to request/use per CPU CHANGE

run_name=EMR_RAPTURE #label to use on output files
indir=$storagenode/$run_name/filtered_alignments
outdir=$storagenode/$run_name/rapture_alignments_keep_pairs
=======
ram_per_cpu=24G #amount of RAM to request/use per CPU CHANGE

run_name=EMR_RAPTURE #label to use on output files
indir=$storagenode/$run_name/filtered_alignments
outdir=$storagenode/$run_name/rapture_alignments
>>>>>>> e9790f51bef6627bc0b7de6a2badc90b976bd1a8

### change me!
site=$1 
bed=$storagenode/$run_name/baits/alignment/Spaced_baits.fas_mrg_buf.bed #filepath of .bed file detailing genomic location of targeted loci
list_of_alignments=$storagenode/$run_name/${site}_filterListRapture # list of alignment files from wrapper-align_to_genome.sh
#---------------------------------------------------------

#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d ./$logfilesdir ]; then mkdir ./$logfilesdir; fi
if [ ! ./$list_of_alignments ]; then echo WARNING list of inds does not exist, go investigate

else
 
#for each line of the list_of_seqfiles file: 
#	execute job that XXXXXX (whatever it does)  
while read alignment
do 

	#define sample name and forward and reverse reads
	sample_name=$(echo $alignment | cut -d "." -f 1) 	
	
	jobname=${sample_name}_iso_rap

	#if input directory doesn't contain at least 1 .bam file; print warning, otherwise proceed with files that are there
	n_inputfiles=($(ls $indir/*.bam | wc -l))
	if [ $n_inputfiles = 0 ]
		then echo WARNING - there are no .bam files in $indir, go investigate
	else	

	#submit job to cluster
	sbatch --job-name=$jobname \
			--export=BED=$bed,CPUS=$cpus,ALIGNMENT=$alignment,RUN_NAME=$run_name,SAMPLE_NAME=$sample_name,STORAGENODE=$storagenode,INDIR=$indir,OUTDIR=$outdir,LOGFILESDIR=$logfilesdir \
			--cpus-per-task=$cpus \
			--mem-per-cpu=$ram_per_cpu \
			--output=./$logfilesdir/${jobname}_%A.out \
			--error=./$logfilesdir/${jobname}_%A.err \
<<<<<<< HEAD
			--time=4:00:00 \
=======
			--time=12:00:00 \
>>>>>>> e9790f51bef6627bc0b7de6a2badc90b976bd1a8
			$executable
			
	echo submitted a job to filter $alignment, reading targeted loci from $bed
	fi		
done < $list_of_alignments

fi
