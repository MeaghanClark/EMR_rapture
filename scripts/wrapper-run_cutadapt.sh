#!/bin/bash

# Last updated 12/27/23 by MI Clark, script format by R Toczydlowski 

# This script runs the run_cutadapt.sbatch file on demultiplexed reads

# run from project directory 

#define variables:
storagenode=/mnt/home/clarkm89 #path to top level of dir where input/output files live
scratchnode=/mnt/scratch/clarkm89/EMR_RAPTURE
#jobname=process_radtags #label for SLURM book-keeping 
run_name=EMR_RAPTURE #label to use on output files
logfilesdir=logfiles_cutadapt #name of directory to create and then write log files to

executable=$storagenode/$run_name/scripts/run_cutadapt.sbatch #script to run

# define where the data is
indir=${scratchnode}/demultiplexed

# define run 
lib_ID=$1 # libraryID, e.g. "P1" 
list_of_inds=./barcodes/${lib_ID}_inds.txt # list of individuals in lib_ID

#---------------------------------------------------------

#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d ./$logfilesdir ]; then mkdir ./$logfilesdir; fi
if [ ! ./$list_of_inds ]; then echo WARNING list of inds does not exist, go investigate

else

#for each line of the list_of_seqfiles file: 
#       execute job that XXXXXX (whatever it does)  
while read seqfiles
do

        #define sample name and forward and reverse reads
        sample_name=$seqfiles
        forwardread=${sample_name}.1.fq.gz
        reverseread=${sample_name}.2.fq.gz

        jobname=${lib_ID}_align

        #define path to where input sequence data live for this dataset
        #and where clean/processed reads should be written to
        indir=$scratchnode/demultiplexed/$lib_ID/
        outdir=$scratchnode/demult_cut

        #if input directory doesn't contain at least 1 .gz file; print warning, otherwise proceed with files that are there
        n_inputfiles=($(ls $indir/*.gz | wc -l))
        if [ $n_inputfiles = 0 ]
                then echo WARNING - there are no .gz files in $indir, go investigate
        else

        #submit job to cluster
        sbatch --job-name=$jobname \
                        --export=FORWARDREAD=$forwardread,REVERSEREAD=$reverseread,CPUS=$cpus,RUN_NAME=$run_name,SAMPLE_NAME=$sample_name,STORAGENODE=$storagenode,SCRATCHNODE=$scratchnode,INDIR=$indir,OUTDIR=$outdir,LOGFILESDIR=$logfilesdir \
                        --cpus-per-task=$cpus \
                        --mem-per-cpu=$ram_per_cpu \
                        --output=./$logfilesdir/${jobname}_${sample_name}_%A.out \
                        --error=./$logfilesdir/${jobname}_${sample_name}_%A.err \
                        --time=22:00:00 \
                        $executable

        echo submitted a job to cut adapters from forward read $forwardread and reverse read $reverseread from individual $sample_name
        fi
done < $list_of_inds

fi

#---------------------------------------------------------











