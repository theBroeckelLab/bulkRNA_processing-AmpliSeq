###### ANALYSIS #########
library(ggplot2)
library(openxlsx)
qc.df=read.xlsx("Z:/Projects/Project Management/Cytotoxicity/Data/TKI_UWGS/run1_to_23_combined_and_SynthegoKO/qc_out.xlsx")

##Total Reads
ggdf=qc.df
ggdf$sample=factor(ggdf$sample, levels=ggdf$sample[order(qc.df$read1, decreasing=F)])
ggplot(ggdf, aes(x=sample, y=read1, color=project))+
  geom_point(size=4)+
  ylab("read pairs")+xlab("")+
  theme_bw()+
  geom_hline(yintercept=quantile(ggdf$read1, 0.50), col="black", alpha=1, lty=2, lwd=1) +
  geom_hline(yintercept=quantile(ggdf$read1, 0.75), col="grey", alpha=0.80, lty=2, lwd=1) +
  geom_hline(yintercept=quantile(ggdf$read1, 0.25), col="grey", alpha=0.80, lty=2, lwd=1) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=6))

##Percent mapped to ref, properly paired, and mapped to ampliseq target gene
ggdf=data.frame(sample=rep(qc.df$sample, 3), 
                qc=rep(c("Perc_mapped","Perc_paired","Perc_on_target"), each=nrow(qc.df)),
                value=c(qc.df$perc_mapped, qc.df$perc_properly_paired, qc.df$perc_on_target))
ggdf$sample=factor(ggdf$sample, levels=qc.df$sample[order(qc.df$perc_mapped, decreasing=F)])
ggdf$qc=factor(ggdf$qc, levels=unique(ggdf$qc))
ggplot(ggdf, aes(x=sample, y=value, color=qc, group=qc))+
  geom_point()+
  geom_line()+
  #ylim(0,1)+
  ylab("percent")+xlab("")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=6))

## total reads all FlowCells
qc.comb=read.xlsx("Z:/Projects/Project Management/Cytotoxicity/Data/TKI_UWGS/QC/All_runs_combined JB.xlsx")
ggdf=qc.comb
ggdf$SeqRun=factor(ggdf$SeqRun, levels=unique(ggdf$SeqRun))
ggplot(ggdf, aes(x=SeqRun, y=read1, group=SeqRun)) +
  geom_boxplot()+
  ylab("read pairs")+
  xlab("FlowCell Run")+
  theme_bw()

##Percent mapped to ref, properly paired, and mapped to ampliseq target gene
ggdf=data.frame(fc=rep(qc.comb$SeqRun, 3), 
                qc=rep(c("Perc_mapped","Perc_paired","Perc_on_target"), each=nrow(qc.comb)),
                value=c(qc.comb$`%.mapped`, qc.comb$`%.properly.paired`, qc.comb$`%.on.target`))
ggdf$fc=factor(ggdf$fc, levels=unique(ggdf$fc))
ggdf$qc=factor(ggdf$qc, levels=unique(ggdf$qc))
ggdf$value[which(ggdf$fc%in%c("21", "22","23"))]=(100*(ggdf$value[which(ggdf$fc%in%c("21", "22","23"))]))
ggplot(ggdf, aes(x=fc, y=value, fill=qc))+
  geom_boxplot()+
  ylab("percent")+xlab("FlowCell Run")+
  theme_bw()
View(qc.df)
