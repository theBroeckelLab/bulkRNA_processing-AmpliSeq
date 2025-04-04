Author: Ryan Gallagher (04-04-2025) 

## NWGC AmpliSeq Quantification Documentation


This document discusses each step in our [BulkRNA AmpliSeq Pipeline](https://github.com/theBroeckelLab/bulkRNA_processing-AmpliSeq/blob/main/ampliSeq_for_RCC/ampliSeq_master.txt). 

## Preliminary Steps

#### Data Transfer

Our RNA data is sequenced by NWGC. We get that data via transfer over Globus connect. See the [Globus Connect for Broeckel Lab]() tutorial for setting up Globus for file transfers.

#### Extract .tar file

Extracting from the archive will produce 4 files per sample. R1 & R2 are paired-end reads each with 2 lanes:

* R1
  * Lane 1
  * Lane 2
* R2
  * Lane 1
  * Lane 2

The command to extract:

`tar -xvf BROECKEL_RUN_TAR.tar`

#### Parse Sample Names & Create Job Array Task File

To properly create SLURM jobs, we will create a job array task file that includes trimmed sample names from our datasets.

This task array will have two columns with the first being an index (1, 2, ..., n) and the second column being the trimmed sample names. The lanes, replicate numbers, and extensions will be trimmed - i.e. ***1088\_CER\_1x\_S35\_L001\_R1\_001.fastq.gz*** will become ***1088\_CER\_1x\_S35***. The following command accomplishes this:

```
ls -1 | 
sed 's/_L00._R._001.fastq.gz//' | uniq |
awk '{printf("%1d\t%s\n", NR, $0)}' | 
awk 'BEGIN {OFS="\t"; print "ArrayTaskID","SampleName"} {print $0, ""}' > sample_names.txt
```

#### Copy Scripts into Scratch Directory

Data processing is run in the `/scratch` directory. So, we need to copy processing scripts from group (as they are written to reference items in the working directory). 

```
cp /group/ubroecke/work/MAPLE/scripts/ampliSeq_scripts/ampliSeq_for_RCC/combineLanes.sh ./
cp /group/ubroecke/work/MAPLE/scripts/ampliSeq_scripts/ampliSeq_for_RCC/alignment_forRCC.sh ./
cp /group/ubroecke/work/MAPLE/scripts/ampliSeq_scripts/ampliSeq_for_RCC/quantification_forRCC.sh ./
cp /group/ubroecke/work/MAPLE/scripts/ampliSeq_scripts/ampliSeq_for_RCC/qc_generate.R ./
cp /group/ubroecke/work/MAPLE/scripts/ampliSeq_scripts/ampliSeq_for_RCC/combine_ampliSeq.R ./

```

#### Copy References into Scratch Directory

Similarly, we will need the hg ampliseq reference and annotation files in our working directory.

```
cp -R /group/ubroecke/work/MAPLE/references/illumina_ampliseq_ref/Genomes/Illumina/hg19/genome/Sequence/BWAIndex ./
cp /group/ubroecke/work/MAPLE/references/illumina_ampliseq_ref/Transcriptome.rna_manifest.20190313.geneExpression.targets.gtf ./
cp /group/ubroecke/work/MAPLE/references/illumina_ampliseq_ref/Transcriptome.rna_gene_targets2.txt ./

```

## Main Steps

### 1. Concatenate Lanes 1 & 2

Lanes for each sample replicate must be concatenated. For example, ***1088\_CER\_1x\_S35\_L001\_R1\_001.fastq.gz*** and ***1088\_CER\_1x\_S35\_L002\_R1\_001.fastq.gz*** should be concatenated into one large fastq.gz file. The `combineLanes.sh` script will do this.

#### 1.1 `combinedLanes.sh` 

This file contains a simple for loop to extract the sample names from the list of files that were extracted from the .tar. This creates a sample names file with unique samples that do not include the lane or R1/R2 tag and creates a `sample_names.txt` file.

**NOTE: Sometimes the output includes the `sample_names.txt` file! This will create an empty column in your final quantification table names called "sample_names". Make sure it's not in this file, and we should fix this!**



### 2. Alignment

This is the first of two computationally intensive steps in our pipeline. We will create a job script to submit to the RCC scheduler using `$ sbatch alignment_forRCC.sh`. Estimated time to complete: **1.5 - 2 hours**.

#### 2.1 `alignment_forRCC.sh`

The header of this script designates the SLURM job directives:

```
#SBATCH --job-name=fc24-alignment # Names job for monitoring
#SBATCH --nodes=1 # Requests 1 Node
#SBATCH --ntasks=6 # Requests 6 tasks in parallel
#SBATCH --mem-per-cpu=4gb # Allocate 4gb of CPU memory
#SBATCH --cpus-per-task=8 # Each task gets 8 cores
#SBATCH --time=00:15:00 # Set a task time limit 
#SBATCH --array=1-160%6 # Submits 160 tasks but allow 6 to run at a time
#SBATCH --account=ubroecke # Charges specified account 
#SBATCH --output=job_%j.out # Logs output file
```

The SLURM directives specified will be applied to each array tasks (i.e. you only need to specify the resources required for a single sample/task, not the entire array). 

The HPC has 60 nodes and 48 cores (360G RAM) per node. 



The body of this file sets the BWAIndex/genome.fa as the reference, then specfies a sample name from the `sample_names.txt` file. The sample selected gets aligned to the reference genome to build a `.fastq.gz` file, which is then converted into a `.bam` file via a `samtools view` command. This `.bam` is then sorted, indexed, and output as a `.sorted.bam` file. Alignment statistics are also produced. 

This process is run across all items in the `sample_names.txt` file as structured SLURM job submissions. 



### 3. Expression Quantification

We use ***htseq-count*** for quantification. This is our other most computationally intensive step. We will build another job script that will be submitted to the RCC. This script is `quantification_forRCC.sh`. 



#### 3.1 `quantification_forRCC.sh`

For this script, our SLURM options are set the same as the previous file. This file loads the `htseq` and `samtools` modules and defines the genome annotation file (imported in the preliminary steps). 

This file uses `sample_names.txt` as a config to run across all of our samples. 

The quantification exists in a single line command:

`htseq-count -f bam -r pos -s no ${sample}.sorted.bam $ANNOTATION > ${sample}.counts`

This command outputs a `.bam` sorted by `position` while identifying our data as "non-strang specific".

Documentation is found [here](https://htseq.readthedocs.io/en/latest/htseqcount.html)



### 4. Generate QC Tables

In the **Alignment** step, we created files named `*_alignstats.txt`. These files yield QC metrics. This step will combine these stats into one file which will be output as `qc_out.xlsx`. The file used to do this is `qc_generate.R`. 



#### 4.1 `qc_generate.R`

This R script reads in all files that are labeled `*_alignstats.txt`, and feeds them into an R dataframe. This dataframe has quality control information about our reads like "% Properly Paired" or "% Mapped to Target". The resulting dataframe has one row per alignment and is saved as `qc_out.xlsx` and `qc_out.txt`. 



### 5. Combine Counts Files into Counts Matrix

The final step is to merge all of our counts files into one large combined counts file. The resulting table will have genes as rows and samples as columns. The file used to do this is `combine_ampliSeq.R`. 



#### 5.1 `combine_ampliSeq.R`

This R script takes the names of our samples from the files in our working directory that end with `*.counts` and takes the names of the targeted genes from our imported file `Transcriptome.rna_gene_targets2.txt`. This script initializes a matrix with the samples as column and genes as rows, then proceeds to fill based on the counts data.



The resulting, final file is `combined_counts.txt`. 



