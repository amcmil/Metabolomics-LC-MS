#Amy McMillan 08/04/2016

#-------------------------------------------
#summary statistics for untargeted LC-MS metabolomics
#-------------------------------------------

#Note: assumes data is on log2 scale. Metabolite table with no m/z or rt info is "met_t", #original output from xcms with rt, mz etc is "peaks"
#metadata table is "mdata"
#identify metabolites whos abundance differs significantly between two or more groups
#for 2 groups use "wilcox.test", for >2 use "kruskal.test"

#read in original peaktable with m/z and rt values
peaks <- read.table("peaklist_raw.txt",  header=T, check.names=F, row.names=1, sep="\t")

#read in zero replaced and logged peaktable 
met <- read.table("peaklist_log2.txt",  header=T, check.names=F, row.names=1, sep="\t")

#transpose
met_t<-t(met)

#create variable "order" with sample order from metabolite table
order<-rownames(met_t)

#import metadata with disease vs control, fungal species etc, order same as peaklist
mdata <- read.table("metadata.txt",  header=T, check.names=F, row.names=1, sep="\t") 
mdata<-as.data.frame(mdata[order,])

#wilcox test
out<-matrix(data=NA, nrow =1, ncol=ncol(met_t))
colnames(out)<-colnames(met_t)

for(i in 1:ncol(met_t)){
t<-wilcox.test(met_t[,i] ~ mdata$"group")
out[,i]<-t$p.value 
}

#fdr correction on pvalues to account for multiple testing
td<-as.data.frame(out)
fdr=p.adjust(td,method="fdr")
pval<-t(rbind(out,fdr))

#calculate average peak area within conditions
av_cond1<-matrix(data=NA, nrow =1, ncol=ncol(met_t))
colnames(av_cond1)<-colnames(met_t)

for(i in 1:ncol(met_t)){
av<-(mean(met_t[,i] [which (mdata$"group"=="1" )]))
av_cond1[,i]<-2^av 
}

av_cond2<-matrix(data=NA, nrow =1, ncol=ncol(met_t))
colnames(av_cond2)<-colnames(met_t)

for(i in 1:ncol(met_t)){
av<-(mean(met_t[,i] [which (mdata$"group"=="2" )]))
av_cond2[,i]<-2^av 
}

#calculate average fold change of metabolite using geometric mean  
out_fc_1<-matrix(data=NA, nrow =1, ncol=ncol(met_t))
colnames(out_fc_1)<-colnames(met_t)

for(i in 1:ncol(met_t)){
fc<-(mean(met_t[,i] [which (mdata$"group"=="1" )]))-((mean(met_t[,i] [which (mdata$"group"=="2" )])))
2^fc
out_fc_1[,i]<-2^fc 
}

out_fc_2<-matrix(data=NA, nrow =1, ncol=ncol(met_t))
colnames(out_fc_2)<-colnames(met_t)

for(i in 1:ncol(met_t)){
fc<-(mean(met_t[,i] [which (mdata$"group"=="2" )]))-((mean(met_t[,i] [which (mdata$"group"=="1" )])))
2^fc
out_fc_2[,i]<-2^fc 
}

#transpost pvalue table 
pval_t<-t(pval)

#Contatinate and add mz and rt in minutes from "peaks"
sig<-t(rbind(peaks[,1],(peaks[,4])/60,pval_t,av_cond1,av_cond2, out_fc_1,out_fc_2))

#write table
write.table(sig,"stats_summary.txt",sep="\t",col.names=NA)

