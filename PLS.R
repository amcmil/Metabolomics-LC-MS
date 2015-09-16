#for good summary of pls functions
#http://beyondvalence.blogspot.ca/2014/02/r-partial-least-squares-regression.html

library(pls)

#construct pls model
pls <- plsr(mdata$"group"~
as.matrix(met_t),ncomp=10,validation="CV",segments=10,segment.type=c("random"),jackknife=TRUE,trace=TRUE)

#to get percent variation explained
summary(pls)

#plot scores
plot(pls$scores[,1],pls$scores[,2],col=mdata$"group",pch=20,main="PLS-DA")
summary(pls)

#plot loadings
plot(pls$loadings[,1],pls$loadings[,2],col=mdata$"group",pch=20,main="PLS-DA")

#cross validation
#select number of component where RMSEP is lowest
pls.RMSEP<-RMSEP(pls,estimate="CV")
min<-which.min(pls.RMSEP$val)-1
plot(RMSEP(pls),legendpos="topright",main="RMSEP curve")
points(min,min(pls.RMSEP$val),pch=20,col="red",cex=1.2)

#CV plot with selected number of components
plot(pls, ncomp = 3, asp = 1, line = TRUE,main="CV")

#Q2
R2(pls) 

#R2
R2(pls, estimate = "train") 

#variation explained by all variables in metadata in loop 

for(i in 1:nrow(mdata)){
	s <- rownames(mdata)[i]
	pls_out <- plsr(mdata[i,] ~ as.matrix(met_t),
	ncomp=10,validation="CV",segments=10,segment.type=c("random"), jackknife=TRUE)
	s <- append(s, capture.output(summary(pls_out)), after=length(s))
	write(s, file="pls_percent_var_explained_mdata.txt", append=TRUE) 
}



