#!/bin/bash

## create sample name file with unique sample names (without lane or R1/R2)


## loop through sample_names.txt and concatenate L1 and L2
nsamples=$(tail -1 sample_names.txt | cut -f 1)
i=0
for filename in ` awk 'NR>1 {print $2}' sample_names.txt`; do
	i=$(($i+1))
	echo "file $i of $nsamples"
	cat ${filename}*L00*_R1*.fastq.gz > ${filename}_R1.fastq.gz
	cat ${filename}*L00*_R2*.fastq.gz > ${filename}_R2.fastq.gz
done

