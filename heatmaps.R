#Amy McMillan 16/09/2015

#----------------
#heatmaps
#---------------
#for extensive list of adjustable parameters see ?heatmap.2
library(gplots)

met<- read.table("peaklist.txt",  header=T, check.names=F, row.names=1, sep="\t")
dm<-as.matrix(met)
 
pdf("heatmap.pdf",height=20, width=15)
heatmap.2(dm,hclustfun = function(x) hclust(x,method = "average"),distfun = function(x) dist(x,method = "manhattan"),dendrogram=("row"),trace="none", col=rev(heat.colors(1000)), scale="row", margins=c(8,20),Colv="NA",keysize=0.8,cexRow=0.8, cexCol=0.8,key.title=NULL, key.xlab=NULL,key.ylab=NULL,)
dev.off()