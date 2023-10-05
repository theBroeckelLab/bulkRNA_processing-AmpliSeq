#!/bin/bash
#SBATCH --job-name=fc24-alignment
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --mem-per-cpu=4gb
#SBATCH --cpus-per-task=8
#SBATCH --time=00:15:00
#SBATCH --array=1-160%6
#SBATCH --account=ubroecke
#SBATCH --output=job_%j.out

#### these settings ####
#Time: 1h37min
#Cores per node: 48
#CPU Utilized: 00:31:49
#CPU Efficiency: 17.44% of 03:02:24 core-walltime
#Job Wall-clock time: 00:03:48
#Memory Utilized: 1.37 GB
#Memory Efficiency: 0.72% of 192.00 GB#	
########################


# Any Slurm directives that you specify for resources at the top of your script e.g. --cpus-per-task or --nodes or --mem etc will be applied to each array task i.e. you only need to specify the resources required for a single sample/task, not the entire array.
# HPC has 60 nodes and 48 cores (360G RAM) per node


## load modules
module load bwa
module load samtools

## define genome reference
REF="./BWAIndex/genome.fa"

# Specify input directory
#

# Specify the path to the config file
config="./sample_names.txt"

# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

## Run alignment, sorting, and QC
bwa mem -M -t 8 -R "@RG\tID:${sample}\tPL:ILLUMINA\tSM:${sample}" -L 5 $REF ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz | samtools view -bS -@ 8 - > ${sample}.bam
samtools sort -m 4G -@ 8 ${sample}.bam -o ${sample}.sorted.bam
samtools index -@ 8 ${sample}.sorted.bam
samtools flagstat ${sample}.sorted.bam > ${sample}_alignstats.txt

## unload modules
module unload bwa
module unload samtools
