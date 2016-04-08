#Amy McMillan 16/09/2015
#-----------------------------
# PCA
#-----------------------------

#requires two packages. FactoMineR is for PCA. MetabolAnalyze does pareto scaling.
library(FactoMineR)
library(MetabolAnalyze)

#read in metabolite table. Make sure it is log scale with no mz or rt info 
met <- read.table("peaklist_log2_pos_working_Study_IDs_58_final.txt",  header=T, check.names=F, row.names=1, sep="\t")

#PCA function requires Rows are sampleIDs, columns are metabolites
met_t<-t(met)

#create variable "order" with sample order from metabolite table
order<-rownames(met_t)

#import metadata with disease vs control, fungal species etc, order same as peaklist
mdata <- read.table("../../../../metadata_recalculated_wfh_final_58.txt",  header=T, check.names=F, row.names=1, sep="\t") 
mdata<-as.data.frame(mdata[order,])

#pareto scale data
met_s<-scaling(met_t,type="pareto")

#do PCA
pca<-PCA(met_s, graph = FALSE, scale.unit=FALSE)

#scoreplot (samples)
plot(pca,label=c("none"))

#loadings plot (metabolites))
plot(pca,choix="var",pch=0.8)

#color scoreplot by metadata columns 
pdf("scoreplot_comp3_4_SAM_cont.pdf",height=5.5,width=5)
plot(pca,col.ind=mdata$"group",label=c("none"),cex=1.8,axes=c(3,4),palette=palette(c("red2","black")))
legend("topright",legend=unique(mdata$"group"),cex=1.2,pt.cex=1.8,col=unique(mdata$"group"),pch=20)
dev.off()

pdf("scoreplot_comp3_4_SAM_cont_with_labels.pdf",height=5.5,width=5)
plot(pca,col.ind=mdata$"group",cex=0.6,axes=c(3,4),palette=palette(c("red2","black")))
legend("topright",legend=unique(mdata$"group"),cex=1.2,pt.cex=1.8,col=unique(mdata$"group"),pch=20)
dev.off()

#for heat colors use palette=palette(heat.colors(10))
#to plot other components use axes=c(2,3) to plot comp 2 vs 3 for example

