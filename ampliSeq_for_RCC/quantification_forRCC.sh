#!/bin/bash
#SBATCH --job-name=test-quantification
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --mem-per-cpu=4gb
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --array=1-160%6
#SBATCH --account=ubroecke
#SBATCH --output=job_%j.out

# HPC has 60 nodes and 48 cores (360G RAM) per node

#### Job Output w/ current settings ########
#Time: 8.5hours
#Job ID: 1723657
#Array Job ID: 1723657_160
#Cluster: hpc2020
#User/Group: jblamer/sg-ubroecke
#State: COMPLETED (exit code 0)
#Nodes: 1
#Cores per node: 6
#CPU Utilized: 00:20:08
#CPU Efficiency: 16.63% of 02:01:06 core-walltime
#Job Wall-clock time: 00:20:11
#Memory Utilized: 1.17 GB
#Memory Efficiency: 4.88% of 24.00 GB
#############################################

## load modules
module load htseq
module load samtools

## define genome annotation file
ANNOTATION="Transcriptome.rna_manifest.20190313.geneExpression.targets.gtf"

# Specify input directory
#

# Specify the path to the config file
config="./sample_names.txt"

# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

## Run alignment, sorting, and QC
htseq-count -f bam -r pos -s no ${sample}.sorted.bam $ANNOTATION > ${sample}.counts

## unload modules
module unload htseq
module unload samtools


##may need to run sacct -j <jobID> to confirm all array jobs completed successfully
