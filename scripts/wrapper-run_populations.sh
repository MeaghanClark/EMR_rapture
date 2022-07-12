#!/bin/bash
		
# Last updated 07/12/2022 by MI Clark, script format by R Toczydlowski 

# This script is in a series of scripts going through the STACKS pipeline.
		
# run from project directory (where you want output directory to be created)

#define variables:
storagenode=/mnt/home/clarkm89 #path to top level of dir where input/output files live
jobname=populations #label for SLURM book-keeping 
run_name=EMR_RAPTURE #label to use on output files
logfilesdir=logfiles_re_pop #name of directory to create and then write log files to

executable=$storagenode/$run_name/scripts/run_populations.sbatch #script to run
popmap=popmap.txt

# define where the data is
indir=rapture_ref_map_keep_pairs_output
outdir=rapture_ref_map_keep_pairs_output_rerun_71222

cpus=20 #20 #number of CPUs to request/use per dataset 
ram_per_cpu=15G #15 #amount of RAM to request/use per CPU 

#---------------------------------------------------------
#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d ./$logfilesdir ]; then mkdir ./$logfilesdir; fi
if [ ! -d ./$outdir ]; then mkdir ./$outdir; fi

#submit job to cluster

sbatch --job-name=$jobname \
--export=JOBNAME=$jobname,CPUS=$cpus,RUN_NAME=$run_name,STORAGENODE=$storagenode,LOGFILESDIR=$logfilesdir,OUTDIR=$outdir,INDIR=$indir,POPMAP=$popmap,M=$M,N=$n \
--cpus-per-task=$cpus \
--mem-per-cpu=$ram_per_cpu \
--output=./$logfilesdir/${jobname}_%A.out \
--error=./$logfilesdir/${jobname}_%A.err \
--time=168:00:00 \
$executable


echo ----------------------------------------------------------------------------------------
echo I am explorting: $jobname, $cpus, $run_name, $storagenode, $logfilesdir, $outdir, $indir, $popmap, $M, $n \
echo My executable is $executable		
