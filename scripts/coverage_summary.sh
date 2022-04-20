#!/bin/bash

#--------------- EXECUTABLE ---------------

# for RAPTURE data. 
#this script calculates coverage summaries from the output of wrapper-compute_coverage.sh

# create summary of coverage
list=filterList
storagenode=/mnt/home/clarkm89/EMR_RAPTURE/coverage_post_alignment
echo "individual,5x,10x,20x,50x,mapped2rapt" >> ${storagenode}/summary/indivcovsum.csv
seqno=$(wc -l ${list} | awk '{print $1}')
x=1
while [ $x -le $seqno ] 
do
	string="sed -n ${x}p ${list}" 
	seqname=$($string)  
	seqname=$(echo $seqname | cut -d "." -f 1)
 	five=$(awk '$4>4' ${storagenode}/${seqname}.cov | wc -l)
	ten=$(awk '$4>9' ${storagenode}/${seqname}.cov | wc -l)
	twenty=$(awk '$4>19' ${storagenode}/${seqname}.cov | wc -l)
	fifty=$(awk '$4>49' ${storagenode}/${seqname}.cov | wc -l)
	total=$(awk 'NR > 3 { print $4 }' ${storagenode}/${seqname}.cov | paste -sd+ - | bc)
#	fq=$(cat ${alignment_location}/${seqname}.fq | wc -l)
#	let filt=$fq/4
	echo $seqname $five $ten $twenty $fifty $total | tr ' ' , >> ${storagenode}/summary/indivcovsum.csv
	x=$(( $x + 1 ))
done



