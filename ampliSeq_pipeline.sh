## PA Pipeline for AmpliSeq data
##**##test run21: xGB 10.25 hours

##Step 0 : Parse Sample Names 
##**## do this manually ##**##
#ls -1 *fastq* | cut -f 1,2,3,4 -d "_" | uniq > sample_names.txt

##Step 1 : Merge fastq.gz files for R1 and R2 from Lane 1 and Lane 2
cat ./sample_names.txt | xargs -L 1 | parallel -P 4 'cat ./{}*L00*_R1*.fastq.gz > {}_R1.fastq.gz'
cat ./sample_names.txt | xargs -L 1 | parallel -P 4 'cat ./{}*L00*_R2*.fastq.gz > {}_R2.fastq.gz'

##Step 2 : Align to reference genome
cat ./sample_names.txt | xargs -L 1 | parallel -P 3 'bwa mem -M -t 8 -R "@RG\tID:{}\tPL:ILLUMINA\tSM:{}" -L 5 /data/ref/illumina_ampliseq_ref/Genomes/Illumina/hg19/genome/Sequence/BWAIndex/genome.fa {}_R1.fastq.gz {}_R2.fastq.gz | samtools view -bS -@ 8 - > {}.bam'

##Step 3 : sort bam
cat ./sample_names.txt | xargs -L 1 | parallel -P 4 'samtools sort -m 4G -@ 1 {}.bam -o {}.sorted.bam'
##[paggarwal@maple:/data/Illumina/misc/ampliseq_test_pipeline]$samtools index 110_reprep.sorted



##Step 4 : calculate alignment, qc, etc
cat ./sample_names.txt | xargs -L 1 | parallel -P 4 'samtools flagstat {}.sorted.bam > {}_alignstats.txt'
Grep "read1\|read2\|mapped\|properly paired" *alignstats.txt > all_samples_flagstat.txt
#ls -1 *sorted.bam | xargs -L 1 | parallel -P 4 'htseq-qa -t bam {} -m 60 -o {.}.pdf' 


##Step 5 : Quantification with htSeq
cat ./sample_names.txt | xargs -L 1 | parallel -P 10 'htseq-count -f bam -r pos -s no {}.sorted.bam /data/ref/illumina_ampliseq_ref/Transcriptome.rna_manifest.20190313.geneExpression.targets.gtf > {}.counts'


##Step 6 : Parsing mis-target counts
#pull EOF values for non-target reads per sample
grep "^__" c*.counts > all_samples_nontarget_reads.txt
#summed counts per sample
my_awk='{total +=$2} END{print total}'
ls -1 *.counts | xargs -L 1 | parallel -P 1 "grep -v "^__" {} | awk '$my_awk'" > summed_reads_perSample.txt

##Step 6 : generate QC tables
Rscript /data/scripts/ampliSeq_scripts/qc_generate.R

##Step 7 : combine counts files
Rscript /data/scripts/ampliSeq_scripts/combine_ampliSeq.R

#sed 's/[.].*$//g' raw_counts.counts | cut -f 1 > run4_names.txt
#paste run4_names.txt raw_counts.counts | cut -f 1,3-162 > raw_counts_nwgc_run4_genenames.counts
#paste raw_counts_nwgc_run1_genenames.counts raw_counts_nwgc_run2_genenames.counts raw_counts_nwgc_run3_genenames.counts raw_counts_nwgc_run4_genenames.counts | cut -f 1-161,163-322,324-483,485-644 > raw_counts_nwgc_run1_2_3_4_combined_genenames.counts
#cat ../../../../../run1_good_111419/samplenames_run1.txt run2_samplenames.txt > samplenames_run1_run2.txt 

