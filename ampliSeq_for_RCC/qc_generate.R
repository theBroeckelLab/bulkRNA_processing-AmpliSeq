library(openxlsx)

##read-in sample names
#samplenames=as.character(read.csv("./sample_names.txt", sep="\n", row.names = NULL, header=F)$V1)
#setwd("/scratch/g/ubroecke/fc24_processing/")
samplenames <- list.files("./", pattern="*alignstats.txt$", full.names=FALSE)
samplenames=gsub("_alignstats.txt","",samplenames)

##read-in alignstats
files.alignstats <- list.files("./", pattern="*alignstats.txt$", full.names=TRUE)
files.alignstats=lapply(files.alignstats, read.csv, sep="\t", header=F, stringsAsFactors=F)
names(files.alignstats)=samplenames

##read-in counts
files.counts <- list.files("./", pattern="*.counts$", full.names=TRUE)
files.counts=lapply(files.counts, read.csv, sep="\t", header=F, stringsAsFactors=F)
names(files.counts)=samplenames
## remove mis-aligned reads from ends of files
for (i in 1:length(files.counts)) {files.counts[[i]]=files.counts[[i]][-grep("^__", files.counts[[i]]$V1),]}

## initiate dataframe
qc.cols=c("sample","project","total_reads","mapped","perc_mapped","read1","read2","properly_paired", "perc_properly_paired","mapped_to_target","pairs_mapped_to_target","perc_mapped_to_target",	"perc_on_target",	"genes_gt10")
qc.df=data.frame(matrix(NA, ncol=length(qc.cols), nrow=length(files.alignstats)))
colnames(qc.df)=qc.cols
qc.df$sample=samplenames

for (i in 1:nrow(qc.df)) {
  split.aligntxt=strsplit(files.alignstats[[i]]$V1, split=" ")
  qc.df$total_reads[i]=as.numeric(split.aligntxt[[1]][1])
  qc.df$mapped[i]=as.numeric(split.aligntxt[[7]][1])
  qc.df$read1[i]=as.numeric(split.aligntxt[[10]][1])
  qc.df$read2[i]=as.numeric(split.aligntxt[[11]][1])
  qc.df$properly_paired[i]=as.numeric(split.aligntxt[[12]][1])
  qc.df$mapped_to_target[i]=sum(files.counts[[i]]$V2)
  qc.df$genes_gt10[i]=length(which(files.counts[[i]]$V2>=10))
}
qc.df$perc_mapped=qc.df$mapped/qc.df$total_reads
qc.df$perc_properly_paired=qc.df$properly_paired/qc.df$mapped
qc.df$pairs_mapped_to_target=qc.df$mapped_to_target*2
qc.df$perc_mapped_to_target=qc.df$pairs_mapped_to_target/qc.df$mapped
qc.df$perc_on_target=qc.df$pairs_mapped_to_target/qc.df$properly_paired


## write table
write.xlsx(qc.df, "./qc_out.xlsx")
write.table(qc.df, "./qc_out.txt", col.names=T, row.names = F, quote=F, sep="\t")


