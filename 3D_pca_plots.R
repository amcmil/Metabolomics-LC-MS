#Amy McMillan 16/09/2015

#--------------------
#3D pca plots
#--------------------

library(rgl)

#do pca as normal here

#make variable "score" containing only the score plot coordinates
scores<-as.data.frame(pca$ind$coor)

#plot 3d pca. 
plot3d(scores$Dim.1,scores$Dim.2,scores$Dim.3,size=5,col=as.integer(mdata$group))
rgl.postscript( "testplot.pdf", fmt = "pdf", drawText = TRUE )