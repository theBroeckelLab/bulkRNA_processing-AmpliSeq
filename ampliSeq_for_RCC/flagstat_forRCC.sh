#!/bin/bash
#SBATCH --job-name=fc24-flagstat
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --mem-per-cpu=2gb
#SBATCH --cpus-per-task=4
#SBATCH --time=00:01:00
#SBATCH --array=1-160%6
#SBATCH --account=ubroecke
#SBATCH --output=job_%j.out

## load modules
module load samtools

# Specify the path to the config file
config="./sample_names.txt"

# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

## Run alignment, sorting, and QC
samtools flagstat ${sample}.sorted.bam > ${sample}_alignstats.txt

## unload modules
module unload bwa
module unload samtools
