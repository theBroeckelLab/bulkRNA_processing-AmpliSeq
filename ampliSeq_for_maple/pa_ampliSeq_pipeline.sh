## Pipeline for AmpliSeq data

##**## Run Prelim Steps One at a Time Manually ##**##
## Prelim step 1 data transfers
#-transfer data from NWGC to oak.phys.mcw.edu via Globus
#-then transfer data from oak.phys.mcw.edu to maple.phys.mcw.edu via CoreFTP

## Prelim step 2: extract the tar archive
#tar -xvf broeckel_run2.tar

##Prelim step 3: Parse Sample Names 
##not a great way to do this, need to trim sample name after lane (L001_)
#ls -1 *fastq* | cut -f 1,2,3,4 -d "_" | uniq > sample_names.txt
#ls -1 |  sed 's/_L00._R._001.fastq.gz//' > sample_names.txt


##**## Remaining steps can be run at once ##**##
##Step 1 : Merge fastq.gz files for R1 and R2 from Lane 1 and Lane 2
cat ./sample_names.txt | xargs -L 1 | parallel -j 4 'cat ./{}*L00*_R1*.fastq.gz > {}_R1.fastq.gz'
cat ./sample_names.txt | xargs -L 1 | parallel -j 4 'cat ./{}*L00*_R2*.fastq.gz > {}_R2.fastq.gz'

##Step 2 : Align to reference genome
cat ./sample_names.txt | xargs -L 1 | parallel -j 3 'bwa mem -M -t 8 -R "@RG\tID:{}\tPL:ILLUMINA\tSM:{}" -L 5 /data/ref/illumina_ampliseq_ref/Genomes/Illumina/hg19/genome/Sequence/BWAIndex/genome.fa {}_R1.fastq.gz {}_R2.fastq.gz | samtools view -bS -@ 8 - > {}.bam'

##Step 3 : sort bam
cat ./sample_names.txt | xargs -L 1 | parallel -j 4 'samtools sort -m 4G -@ 8 {}.bam -o {}.sorted.bam'



##Step 4 : calculate alignment, qc, etc
cat ./sample_names.txt | xargs -L 1 | parallel -j 24 'samtools flagstat {}.sorted.bam > {}_alignstats.txt'
grep "read1\|read2\|mapped\|properly paired" *alignstats.txt > all_samples_flagstat.txt


##Step 5 : Quantification with htSeq
cat ./sample_names.txt | xargs -L 1 | parallel -j 24 'htseq-count -f bam -r pos -s no {}.sorted.bam /data/ref/illumina_ampliseq_ref/Transcriptome.rna_manifest.20190313.geneExpression.targets.gtf > {}.counts'


##Step 6 : Parsing mis-target counts
##pull EOF values for non-target reads and sum counts per sample
grep "^__" *.counts > all_samples_nontarget_reads.txt
my_awk='{total +=$2} END{print total}'
ls -1 *.counts | xargs -L 1 | parallel -P 1 "grep -v "^__" {} | awk '$my_awk'" > summed_reads_perSample.txt

##Step 6 : generate QC tables
Rscript /data/scripts/ampliSeq_scripts/qc_generate.R

##Step 7 : combine counts files
Rscript /data/scripts/ampliSeq_scripts/combine_ampliSeq.R


