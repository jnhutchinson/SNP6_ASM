#script takes 3 arguments
# 1 - rawUNCUTdata, model data in the last column
# 2 - name of the file with the quantile normalized data
# 3 - name of the file with the adjustment factors

library(affy)	#load the library that does the quantile norm calcs
args<-commandArgs() #load command line arguments into variable

#read in the data, grab the header separately and stow in a variable for now
header<-read.table(args[3], nrows=1)
raw<-read.table(args[3], header=TRUE, row.names=1)

len=length(raw)-1
zeros=rep(0, each=len)
weightvector=c(zeros, 1)
weightvector
#convert data to matrix and do all calcs
rawmatrix<-data.matrix(raw)
quantmatrix<-normalize.quantiles.robust(rawmatrix, copy=TRUE, weights=weightvector)
adjustmatrix<-quantmatrix/rawmatrix

#round up results and convert back to data frames
quant<-data.frame(round(quantmatrix, digits=1))
adjust<-data.frame(round(adjustmatrix, digits=4))

#add back rownames
row.names(quant)<-row.names(raw)
row.names(adjust)<-row.names(raw)

#print to files
write.table(header, args[4], col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
write.table(quant, args[4], sep="\t", append=TRUE, col.names=FALSE, row.names=TRUE, quote=FALSE)
write.table(header, args[5], col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
write.table(adjust, args[5], sep="\t", append=TRUE, col.names=FALSE, row.names=TRUE, quote=FALSE)
