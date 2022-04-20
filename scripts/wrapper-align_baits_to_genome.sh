#!/bin/bash
		
# Last updated 03/16/2022 by MI Clark, originally written by R Toczydlowski 

#  run from project directory (where you want output directory to be created)

# usage: ./scripts/wrapper-align_baits_to_genome.sh baits.fas # library number

#define variables:
storagenode=/mnt/home/clarkm89 #path to top level of dir where input/output files live

logfilesdir=logfiles_alignbaits #name of directory to create and then write log files to
executable=./scripts/align_baits_to_genome.sbatch #script to run 

cpus=1
ram_per_cpu=15G #amount of RAM to request/use per CPU CHANGE

run_name=EMR_RAPTURE #label to use on output files
reference=$storagenode/massasauga/reference/Scatenatus_HiC_v1.1.fasta #filepath of reference file

indir=$storagenode/$run_name/baits
outdir=$storagenode/$run_name/baits/alignment

#### things to change with each run
baits=Spaced_baits.fas
jobname=baits_spaced

#---------------------------------------------------------

#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d ./$logfilesdir ]; then mkdir ./$logfilesdir; fi
if [ ! ${indir}/${baits} ]; then echo WARNING baits file does not exist, go investigate

else
	#submit job to cluster
	sbatch --job-name=$jobname \
			--export=REFERENCE=$reference,BAITS=$baits,RUN_NAME=$run_name,STORAGENODE=$storagenode,INDIR=$indir,OUTDIR=$outdir,LOGFILESDIR=$logfilesdir \
			--cpus-per-task=$cpus \
			--mem-per-cpu=$ram_per_cpu \
			--output=./$logfilesdir/${jobname}_${baits}_%A.out \
			--error=./$logfilesdir/${jobname}_${baits}_%A.err \
			--time=22:00:00 \
			$executable
			
	echo submitted a job to align $baits located at $indir to $reference

fi
