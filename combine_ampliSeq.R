
### combine counts files / generate table ###

## read-in gene and sample names]
snames=as.character(read.csv("./sample_names.txt", header=F)[,1])
gnames=as.character(read.csv("/data/ref/illumina_ampliseq_ref/Transcriptome.rna_gene_targets2.txt", header=F)[,1])
## read-in counts files
filenames.out=list.files("./", pattern="*.counts$", full.names=TRUE)
count.list=lapply(filenames.out, read.csv, row.names=NULL, sep="\t", header=F, stringsAsFactors = F)
##initiate and fill counts table
comb=matrix(NA, nrow=length(gnames), ncol=length(snames))
rownames(comb)=gnames; colnames(comb)=snames
for (i in 1:length(count.list)) {comb[,i]=count.list[[i]][,2]}
##remove QC values at EOFs
comb=comb[-grep("^__", rownames(comb)),]
write.table(comb, "./combined.counts", row.names=T, col.names=T, quote=F, sep="\t")
##check reliabality of counts
#comb2=read.csv("./combined.counts", row.names=1, header=T, sep="\t")
