#!/bin/bash
		
# Last updated 03/02/2022 by MI Clark, script format by R Toczydlowski 

# This script is the first in series of scripts going through the STACKS pipeline.
# This script runs the process_radtags functionality of STACKS. 
		
# run from project directory (where you want output directory to be created)

#define variables:
storagenode=/mnt/home/clarkm89 #path to top level of dir where input/output files live
scratchnode=/mnt/scratch/clarkm89/EMR_RAPTURE
#jobname=process_radtags #label for SLURM book-keeping 
run_name=EMR_RAPTURE #label to use on output files
logfilesdir=logfiles_radtags #name of directory to create and then write log files to

executable=$storagenode/$run_name/scripts/run_Process_radtags.sbatch #script to run

# define where the data is
indir=Rawdata
#outdir=demultiplexed
lib_ids=lib_id.txt

cpus=1 #number of CPUs to request/use per dataset 
ram_per_cpu=25G #25 amount of RAM to request/use per CPU 


#---------------------------------------------------------
#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d ./$logfilesdir ]; then mkdir ./$logfilesdir; fi
#if [ ! -d ./$outdir ]; then mkdir ./$outdir; fi

#submit job to cluster

while read lib_ids
do
	forwardread=EMR_${lib_ids}_*R1_001.fastq.gz
	reverseread=EMR_${lib_ids}_*R2_001.fastq.gz
	jobname=${lib_ids}_process_radtags
	barcodes=$storagenode/$run_name/barcodes/EMR_barcodes_${lib_ids}.txt	
	outdir=demultiplexed/${lib_ids}
	if [ ! -d ./$outdir ]; then mkdir ./$outdir; fi

	sbatch --job-name=$jobname \
	--export=JOBNAME=$jobname,CPUS=$cpus,RUN_NAME=$run_name,STORAGENODE=$storagenode,SCRATCHNODE=$scratchnode,LOGFILESDIR=$logfilesdir,OUTDIR=$outdir,INDIR=$indir,BARCODES=$barcodes,FORWARDREAD=$forwardread,REVERSEREAD=$reverseread \
	--cpus-per-task=$cpus \
	--mem-per-cpu=$ram_per_cpu \
	--output=./$logfilesdir/${jobname}_%A.out \
	--error=./$logfilesdir/${jobname}_%A.err \
	--time=12:00:00 \
	$executable
	

echo ----------------------------------------------------------------------------------------
echo I am explorting: $jobname, $cpus, $run_name, $storagenode, $logfilesdir, $outdir, $outdir_clone, $forwardread, $reverseread, $barcodes 
echo My executable is $executable		

done < $lib_ids
