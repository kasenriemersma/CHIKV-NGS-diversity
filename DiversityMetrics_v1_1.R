## upload package libraries
library("diverse")
library("dplyr")

#### Enter the sample name and output directory for the sample
args<-commandArgs()
outputdir<-args[6]
samplename<-args[7]
outputsample<-paste(outputdir,samplename,sep="/")

## import Perbase text file
table<-paste("ntcounts_",samplename,sep="")
table<-paste(table,".csv",sep="")
ntcounts<-read.csv(table)

## rename column headers
columns<-c("pos","ref","reads_all","mismatches","deletions","insertions","A","C","T","G")
ntcounts<-ntcounts[columns]
colnames(ntcounts)<-c("position","reference","coverage","mismatches","deletions","insertions","A","C","T","G")

## generate absolute SNV columns
ntcounts$A.SNV<-as.numeric(ifelse(ntcounts$reference=="A",0,
                                  ifelse(ntcounts$reference=="C",ntcounts$A,
                                         ifelse(ntcounts$reference=="G",ntcounts$A,
                                                ifelse(ntcounts$reference=="T",ntcounts$A,
                                                       ifelse(ntcounts$reference=="",0,"na"))))))

ntcounts$C.SNV<-as.numeric(ifelse(ntcounts$reference=="C",0,
                                  ifelse(ntcounts$reference=="A",ntcounts$C,
                                         ifelse(ntcounts$reference=="G",ntcounts$C,
                                                ifelse(ntcounts$reference=="T",ntcounts$C,
                                                       ifelse(ntcounts$reference=="",0,"na"))))))

ntcounts$G.SNV<-as.numeric(ifelse(ntcounts$reference=="G",0,
                                  ifelse(ntcounts$reference=="A",ntcounts$G,
                                         ifelse(ntcounts$reference=="C",ntcounts$G,
                                                ifelse(ntcounts$reference=="T",ntcounts$G,
                                                       ifelse(ntcounts$reference=="",0,"na"))))))

ntcounts$T.SNV<-as.numeric(ifelse(ntcounts$reference=="T",0,
                                  ifelse(ntcounts$reference=="A",ntcounts$T,
                                         ifelse(ntcounts$reference=="C",ntcounts$T,
                                                ifelse(ntcounts$reference=="G",ntcounts$T,
                                                       ifelse(ntcounts$reference=="",0,"na"))))))

## generate SNV frequency column
ntcounts$snv.freq<-(ntcounts$mismatches/ntcounts$coverage)

## generate Squared Deviation column
ntcounts$sq.dev<-((ntcounts$snv.freq)^2)

## generate unique SNV column
ntcounts$unique.A.SNV<-ifelse(ntcounts$A.SNV>0,1,0)
ntcounts$unique.C.SNV<-ifelse(ntcounts$C.SNV>0,1,0)
ntcounts$unique.G.SNV<-ifelse(ntcounts$G.SNV>0,1,0)
ntcounts$unique.T.SNV<-ifelse(ntcounts$T.SNV>0,1,0)

unique.A.SNV<-sum(ntcounts$unique.A.SNV)
unique.C.SNV<-sum(ntcounts$unique.C.SNV)
unique.G.SNV<-sum(ntcounts$unique.G.SNV)
unique.T.SNV<-sum(ntcounts$unique.T.SNV)
unique.SNV<-sum(unique.A.SNV,unique.C.SNV,unique.G.SNV,unique.T.SNV)
unique.SNV

## subset data for coverage >= 300
variant<-subset(ntcounts,ntcounts$coverage>299)
variant$snv.freq<-((variant$A.SNV + variant$C.SNV + variant$G.SNV + variant$T.SNV)/variant$coverage)
variant$sq.dev<-((variant$snv.freq)^2)

## calculate root mean squared deviation
rmsd<-sqrt(mean(variant$sq.dev))
rmsd

## calculate shannon entropy
nucvar<-subset(variant,select=c("A","C","G","T"))
nucvar<-as.matrix(nucvar)

entropy.nucvar<-diversity(nucvar, type="entropy")
entropy.mean<-mean(entropy.nucvar$entropy)
entropy.sd<-sd(entropy.nucvar$entropy)

## calculate Gini-Simpson index

ginisimpson<-diversity(nucvar, type="gini-simpson")
gs.mean<-mean(ginisimpson$gini.simpson.C)
ginisimpson$mean.gs.C<-gs.mean
gs.sd<-sd(ginisimpson$gini.simpson.C)

## calculate mutation frequency (SNV/10,000 nt sequenced)
mut.freq<-(sum(variant$A.SNV,variant$C.SNV,variant$G.SNV,variant$T.SNV)/sum(variant$coverage))*10000
mut.freq
unique.mut.freq<-(unique.SNV/sum(variant$coverage))*10000
unique.mut.freq

## calculate nucleotide diversity (pi)
variant$dij<-((variant$snv.freq*variant$coverage)/(((variant$coverage^2)-variant$coverage)/2))
pi<-sum(variant$dij)

## Calculate coverage metrics
percentcovered<-(nrow(variant)/11811*100)
meandepth<-mean(variant$coverage)
percentcovered
meandepth

## merge diversity indices and write file
diversity<-c(percentcovered,meandepth,mut.freq,unique.mut.freq,rmsd,entropy.mean,entropy.sd)
diversity<-t(diversity)
colnames(diversity)<-c("%Covered","MeanDepth","Mut.Freq.per.10K","Unique.Mut.Freq","RMSD","Mean.Shannon.Entropy", "SD.Shannon.Entropy")
diversity<-as.data.frame(diversity)
outputdiversity<-paste(outputsample,"_coverage_and_diversity_metrics.csv",sep="")
write.csv(diversity, file = outputdiversity, row.names = FALSE)