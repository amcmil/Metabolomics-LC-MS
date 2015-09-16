#Amy McMillan 16/09/2015
#-----------------------------
# PCA
#-----------------------------

#requires two packages. FactoMineR is for PCA. MetabolAnalyze does pareto scaling.
library(FactoMineR)
library(MetabolAnalyze)

#read in metabolite table. Make sure it is log scale with no mz or rt info 
met <- read.table("peaktable.txt",  header=T, check.names=F, row.names=1, sep="\t")

#PCA function requires Rows are sampleIDs, columns are metabolites
met_t<-t(met)

#create variable "order" with sample order from metabolite table
order<-rownames(met_t)

#import metadata with disease vs control, fungal species etc, order same as peaklist
mdata <- read.table("metadata.txt",  header=T, check.names=F, row.names=1, sep="\t")
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
plot(pca, col.ind=mdata$"columnID",cex=0.8,palette=palette(c("seagreen3","black")))
legend("topright",legend=unique(mdata$group),col=unique(mdata$group),pch=20)

#for heat colors use palette=palette(heat.colors(10))
#to plot other components use axes=c(2,3) to plot comp 2 vs 3 for example

