#Amy McMillan 16/09/2015

#-------------------------------------------
#summary statistics for untargeted LC-MS metabolomics
#-------------------------------------------

#Note: assumes data is on log2 scale. Metabolite table with no m/z or rt info is "met_t", #original output from xcms with rt, mz etc is "peaks"
#metadata table is "mdata"
#identify metabolites whos abundance differs significantly between two or more groups
#for 2 groups use "wilcox.test", for >2 use "kruskal.test"

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
#repeat for number of conditions

#calculate prevalence of metabolite withing each conditions  
#percent of samples within each condition with E5 or greater for example. 16.60964
#16.60964 = log2(1E5)
prev_cond1<-matrix(data=NA, nrow =1, ncol=ncol(met_t))
colnames(prev_cond1)<-colnames(met_t)

for(i in 1:ncol(met_t)){
prev<-(length(which(mdata$"group"==1 & met_t[,i]>16.59))/length(mdata$"group"==1))
prev_cond1[,i]<-prev 
}

prev_cond2<-matrix(data=NA, nrow =1, ncol=ncol(met_t))
colnames(prev_cond2)<-colnames(met_t)

for(i in 1:ncol(met_t)){
prev<-(length(which(mdata$"group"==2 & met_t[,i]>16.59))/length(mdata$"group"==2))
prev_cond2[,i]<-prev 
}

#transpost pvalue table 
pval_t<-t(pval)

#Contatinate and add mz and rt in minutes from "peaks"
sig<-t(rbind(peaks[,1],(peaks[,4])/60,pval_t,av_cond1,av_cond2,prev_cond1,prev_cond2))

#write table
write.table(sig,"stats_summary.txt",sep="\t",col.names=NA)

