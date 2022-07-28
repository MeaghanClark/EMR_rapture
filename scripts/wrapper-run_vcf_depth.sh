#!/bin/bash

# Last updated 09/23/2020 by MI Clark, script format by R Toczydlowski 

#  run from project directory (where you want output directory to be created)

#define variables:

jobname=vcf_depth #label for SLURM book-keeping 
run_name=EMR_RAPTURE #label to use on output files

# locations

storagenode=/mnt/home/clarkm89 #path to top level of dir where input/output files live
logfilesdir=logfiles_vcfdepth #name of directory to create and then write log files to
indir=$storagenode/$run_name/rapture_ref_map_keep_pairs_output_rerun_72522
outdir=$indir/vcf_depth

# files
executable=$storagenode/$run_name/scripts/run_vcf_depth.sbatch #script to run 
rscript=$storagenode/$run_name/scripts/vcf_depth.R #R script called by executable
vcf=$indir/populations.snps.vcf

cpus=1 #number of CPUs to request/use per dataset 
ram_per_cpu=50G #amount of RAM to request/use per CPU;


#---------------------------------------------------------
#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d ./$logfilesdir ]; then mkdir ./$logfilesdir; fi

#submit job to cluster

sbatch --job-name=$jobname \
	--export=JOBNAME=$jobname,CPUS=$cpus,RUN_NAME=$run_name,STORAGENODE=$storagenode,OUTDIR=$outdir,INDIR=$indir,VCF=$vcf,LOGFILESDIR=$logfilesdir,RSCRIPT=$rscript \
    --cpus-per-task=$cpus \
    --mem-per-cpu=$ram_per_cpu \
    --output=./$logfilesdir/${jobname}_%A.out \
    --error=./$logfilesdir/${jobname}_%A.err \
    --time=168:00:00 \
    $executable