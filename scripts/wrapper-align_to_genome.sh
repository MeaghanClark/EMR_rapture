#!/bin/bash
		
# Last updated 03/16/2022 by MI Clark, originally written by R Toczydlowski 

#  run from project directory (where you want output directory to be created)

# usage: ./scripts/wrapper-align_to_genome.sh P1 # library number
#define variables:
storagenode=/mnt/home/clarkm89 #path to top level of dir where input/output files live
scratchnode=/mnt/scratch/clarkm89/EMR_RAPTURE_alignments

logfilesdir=logfiles_aligntogenome #name of directory to create and then write log files to
executable=./scripts/align_to_genome.sbatch #script to run 

cpus=4 #number of CPUs to request/use per dataset
ram_per_cpu=2G #amount of RAM to request/use per CPU CHANGE

run_name=EMR_RAPTURE #label to use on output files
reference=$storagenode/massasauga/reference/Scatenatus_HiC_v1.1.fasta #filepath of reference file

lib_ID=$1 # libraryID, e.g. "P1" 
list_of_inds=./barcodes/${lib_ID}_inds.txt # list of individuals in lib_ID

#---------------------------------------------------------

#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d ./$logfilesdir ]; then mkdir ./$logfilesdir; fi
if [ ! ./$list_of_inds ]; then echo WARNING list of inds does not exist, go investigate

else
 
#for each line of the list_of_seqfiles file: 
#	execute job that XXXXXX (whatever it does)  
while read seqfiles
do 

	#define sample name and forward and reverse reads
	sample_name=$seqfiles	
	forwardread=${sample_name}.1.fq.gz
	reverseread=${sample_name}.2.fq.gz
	
	jobname=${lib_ID}_align

	#define path to where input sequence data live for this dataset
	#and where clean/processed reads should be written to
	indir=/mnt/scratch/clarkm89/$run_name/demultiplexed/$lib_ID/
	outdir=$storagenode/$run_name/alignments
	
	#if input directory doesn't contain at least 1 .gz file; print warning, otherwise proceed with files that are there
	n_inputfiles=($(ls $indir/*.gz | wc -l))
	if [ $n_inputfiles = 0 ]
		then echo WARNING - there are no .gz files in $indir, go investigate
	else	

	#submit job to cluster
	sbatch --job-name=$jobname \
			--export=REFERENCE=$reference,FORWARDREAD=$forwardread,REVERSEREAD=$reverseread,CPUS=$cpus,RUN_NAME=$run_name,SAMPLE_NAME=$sample_name,STORAGENODE=$storagenode,SCRATCHNODE=$scratchnode,INDIR=$indir,OUTDIR=$outdir,LOGFILESDIR=$logfilesdir \
			--cpus-per-task=$cpus \
			--mem-per-cpu=$ram_per_cpu \
			--output=./$logfilesdir/${jobname}_${sample_name}_%A.out \
			--error=./$logfilesdir/${jobname}_${sample_name}_%A.err \
			--time=22:00:00 \
			$executable
			
	echo submitted a job to align forward read $forwardread and reverse read $reverseread from individual $sample_name to $reference
	fi		
done < $list_of_inds

fi
