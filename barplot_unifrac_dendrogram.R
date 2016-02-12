# Barplot from OTU table and UniFrac distance matrix
# Mar 2013, JM
#----------------------------------------------------------------------------------------

#Expects QIIME input formatted table:
d <- read.table("td_OTU_lacto_sep_rare_rem_131.txt", header=T, sep="\t", row.names=1, skip=0, comment.char="", check.names=FALSE)
##header here
##OTU_ID	015B 	018B 	019B	taxonomy

#Sum OTU counts for each column and get fraction values
#For QIIME input formatted table:
#y <- apply( d[,1:ncol(d)-1] , 2, function(x) { x / sum(x) } )
#excludes the taxonomy column

#FOR THE summarize_taxa.py OUTPUT FROM QIIME USE THIS ONE INSTEAD:
y <- apply( d[,1:ncol(d)] , 2, function(x) { x / sum(x) } )


#---------------------------------------------------------------------------------------------------
# The following chunk of code will group OTUs together that are in low abundance (<0.01 hard-coded)
# Find and replace all ^#
#---------------------------------------------------------------------------------------------------
##"rem" row has to exist
newrowname<-"rem"
#
## create a one-row matrix the same length as data containing "0" for each value (rep.int)
temprow <- matrix(c(rep.int(0,ncol(y))),nrow=1,ncol=ncol(y))
#
## make it a data.frame and give cols the same names as data
newrow <- data.frame(temprow)
colnames(newrow) <- colnames(y)
rownames(newrow) <-newrowname
#
## rbind the empty row to data
#
y <- rbind(y,newrow)
#
##Get the taxa names
##Since rownames have to be unique, paste together OTU number and taxonomic lineage
#
bugnames<-append(as.vector(paste(rownames(d), d$taxonomy, sep="_")), newrowname)
##length(bugnames), class(bugnames)
#
rownames(y)<-bugnames
##duplicate rownames not allowed
#
##### move all bugs present at less than 1% into the "rem" fraction
y <- apply( y , 2, function(x) {
	small <- ( x < 0.01 ) #& ( bugnames != "rem" )
	rare <- sum( x[small] )
	x[small] <- 0 # *** NA
	x["rem"] <- x["rem"] + rare
	return(x)
##	return(rare)
} )
### Check that everyting sums to 1
total <- apply( y, 2, function(x) { sum(x,na.rm=TRUE) } ) # should be all equal to one
if ( any( abs(total-1) > (.Machine$double.eps*16) ) ) stop("bugs were erroneously created or destroyed")
#---------------------------------------------------------------------------------------------------
# End chunk
#---------------------------------------------------------------------------------------------------

#colorscheme
colours <- c("steelblue3", "skyblue1","#cfebff","cornflowerblue","cyan","indianred1", "mediumpurple1","#FFED6F","ivory2","tan1","olivedrab3","royalblue4","mediumorchid3", "darkcyan","pink", "deeppink3", "#C0C0C0", "#999933", "#CCCC99","#9999CC", "#006666", "#3399FF", "#993300", "#CCCC99", "#666666", "#FFCC66", "#6699CC", "#663366", "#9999CC", "#CCCCCC", "#669999", "#CCCC66", "#CC6600", "#9999FF", "#0066CC", "#99CCCC", "#999999", "#FFCC00", "#009999", "#FF9900", "#999966", "#66CCCC", "#C0C0C0", "purple","yellow","green","orange","red","white","purple","yellow","black")


#-----------------------------
# Dendogram from UniFrac dm
#-----------------------------
#get as a matrix
dm<-as.matrix(read.table("weighted_unifrac_dm_no_met_rem.txt", sep="\t", header=T, row.names=1, comment.char="", check.names=FALSE))

#convert to dist structure
dmf<-as.dist(dm)
#"average" is most similar to UPGMA, apparently
dendo <- hclust(dmf, method="average")

# PLOT TOGETHER!!!!
pdf('barplot_dendo_weighted_average_131_rare_otu_rem.pdf', height=8, width=15)

#GG legacy code. Fix size, margins, position
par(mfrow=c(2,1), mar=c(1, 3, 2, 1) + 0.1)
plot(dendo, axes=F, ylab=NULL, ann=F)
#order the barplot 
barplot(y[,dendo$order], space=0, col=colours, las=2)
dev.off()

#-----------
# Optional legend
#---------------
#plot the legend
pdf('legend_barplot_nugent_ord.pdf', height=60, width=16)
barplot(y, space=0, col=colours, ylim=c(-20,0), legend=T)
dev.off()

z<-y[,dendo$order]
write.table(z,"dendo_order.txt",sep="\t",col.names=NA)

