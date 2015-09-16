#Amy McMillan 16/09/2015
#stripcharts to compare individual metabolite abundances between groups

#read in metabolite and metadata table and create variable order with rownames or met 
met <- read.table("peaklist.txt", sep="\t", comment.char="", check.names=FALSE, row.names=1, skip=1)

order<-rownames(met_t)

#import metadata with disease vs control, fungal species etc, order same as peaklist
mdata <- read.table("metadata.txt",  header=T, check.names=F, row.names=1, sep="\t")
mdata<-as.data.frame(mdata[order,])

#transpose
td<-as.data.frame(t(met))

#plotting 
stripchart(as.numeric(as.character(td$metabolite ID))~mdata$"group", data=td, vertical=TRUE, method="jitter", jitter=0.25, pch=20, col="red",main="metaboliteID")

#to loop through metabolites and create stripchart for each 
for(i in 1:ncol(td)){
 	nm <- colnames(td)[i]
	pdf(paste(file=nm,".pdf", sep=""))
	stripchart(as.numeric(as.character(td[,i]))~td$"group", data=td, vertical=TRUE, 	 method="jitter", jitter=0.25, pch=20, col="red",ylim=c(-15,1))
	dev.off()
}
