#Amy McMillan 16/09/2015
#--------------------------------------
#Workflow for peaklist generation from .raw files using xc-ms in R
#--------------------------------------

#use proteowizard (http://proteowizard.sourceforge.net) to convert files to centroid mode (.mzml or .mzxml files) by selecting "peak picking" and level =1 option
#description of centroid vs profile ms data acquisiton
#http://acdlabs.typepad.com/elucidation/2008/03/what-is-the-dif.html
#note xcms nature protocol paper Table 1 suggests HPLC/orbitrap parameters: ppm=2.5, peak width=c(10,60), bw=5, mzwid-0.015,Prefilter=c(3,5000)
#setting ppm at 2.5 gives may give warning to lower ppm. seting ppm=1.0 avoids this

#load library
library(xcms)

#select files 
files<-list.files("/path to files",recursive=TRUE,full.names=TRUE)      

#peak detection
xset<-xcmsSet(files,method="centWave", polarity="positive",prefilter=c(3,5000),ppm=2.5, snthresh=5,peakwidth=c(5,20),noise=100000)

#rt correction
xset1 <- retcor(xset, method="obiwarp", plottype = c("deviation"))

#group samples
xset2 <- group(xset1, bw = 5, minfrac = 0.25, mzwid = 0.015)

#fill peaks--> integrates masses below sn threshold 
xset3 <- fillPeaks(xset2)

#write to peak table
#note RT is reported in seconds. Need to convert to minutes
peaks<-peakTable(xset3,filebase="/path to files")

#write peaklist to file
write.table(peaks,"peaklist_pos_all.txt",sep="\t",col.names=NA)

#to replace zeros with 2/3 min value on per metabolite basis for human data(required before data can be converted to log scale)
#[,9:89] idicates range of sample columns. xc-ms output samples start at column 9.
which.m <-  apply(d[,9:89],1, function(x){min(x[which(x>0)])*0.6667})
for(i in 1:nrow(d)){
x <- which(d[i,] == 0)
d[i,x]  <-  which.m[i]
}

#write log peaklist to file
write.table(peaks_log,"peaklist_pos_log2_all.txt",sep="\t",col.names=NA)

#--------------------------------------------
# Details of xc-ms parameters (from xcms online)
#-----------------------------------------------

#ppm		maximal tolerated m/z deviation in consecutive scans, in ppm (parts per million)
#minimum peak width		minimum chromatographic peak width in seconds note: must be less than max peak width. See also here.
#maximum peak width		maximum chromatographic peak width in seconds note: must be greater than min peak width. See also here.
#mzdiff		minimum difference in m/z for peaks with overlapping retention times, can be negative to allow overlap
#Signal/Noise threshold		Signal/Noise threshold
#Integration method		Integration method. If =1 peak limits are found through descent on the mexican hat filtered data, if =2 the descent is done on the real data. Method 2 is very accurate but prone to noise, while method 1 is more robust to noise but less exact.
#prefilter peaks		Prefilter step for the first phase. Mass traces are only retained if they contain at least [prefilter peaks] peaks with intensity >= [prefilter intensity]
#prefilter intensity		Prefilter step for the first phase. Mass traces are only retained if they contain at least [prefilter peaks] peaks with intensity >= [prefilter intensity]
#Noise Filter		optional argument which is useful for data that was centroided without any intensity threshold, centroids with intensity < noise are omitted from ROI detection

#paramaters groups
#bw		Allowable retention time deviations, in seconds. In more detail: bandwidth (standard deviation or half width at half maximum) of gaussian smoothing kernel to apply to the peak density chromatogram
#mzwid		width of overlapping m/z slices to use for creating peak density chromatograms and grouping peaks across samples
#minfrac		minimum fraction of samples necessary in at least one of the sample groups for it to be a valid group
#View Advanced Options
#max		maximum number of groups to identify in a single m/z slice
#minsamp		minimum number of samples necessary in at least one of the sample groups for it to be a valid group

#fillpeaks
#For each sample, identify peak groups where that sample is not represented. For each of those peak groups, integrate the signal in the region of that peak group and create a new peak.