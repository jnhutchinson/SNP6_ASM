## ----functions-----------------------------------------------------------
# tests if library installed: if yes, loads library; if no, installs library
pkgTest <- function(x)  {
    if (!require(x,character.only = TRUE))
    {
      install.packages(x,dep=TRUE, repo="http://lib.stat.cmu.edu/R/CRAN/" )
      
        if(!require(x,character.only = TRUE)) stop("Package not found")
    }
  }


## ----libraries-----------------------------------------------------------
pkgTest("reshape")
pkgTest("plyr")
pkgTest("getopt")
require(reshape)
require(plyr)
require(getopt)


## ----variables-----------------------------------------------------------
# the directory with the raw non-median summed intensities, passed from python script
args <- commandArgs(trailingOnly = TRUE)
dataDir=args[1]
if(file.exists(dataDir)){
  print(dataDir)
} else {
  print("data directory not found")
}



## ----process-------------------------------------------------------------
# get list of raw intensity files
rawintFiles <- list.files(file.path(dataDir), pattern="rawints" )
# load rawintensity files and median summarize one at a time
medints.l <- llply(rawintFiles, function(rawintFile) {
  print(rawintFile)
  rawints <- read.table(file.path(dataDir, rawintFile), header=T)
  # get the index of the column that contains the intensities
  intCol <- grep(sub(".rawints.txt","", rawintFile), names(rawints))
  # subset to genotype probesets only
  rawints <- subset(rawints, probeset_type=="GenoType")
  # aggregate by probeset_id and allele and take median value for intensities 
  rawints <- aggregate(rawints[,intCol], by=list(probeset_id=rawints$probeset_id, alleleID=rawints$block), median)
  return(rawints)
})
print("data imported")
names(medints.l) <- sub(".rawints.txt", "", rawintFiles)

medints.m <- melt(medints.l, id.vars=c("probeset_id", "alleleID"), measure.vars="x")
print("data melted")
medints.d <- cast(medints.m, probeset_id + alleleID ~ L1)
print("data cast")
head(medints.d)

medints.d <- medints.d[grep("SNP_A", medints.d$probeset_id),]
medints.d$alleleID[medints.d$alleleID==0] <- "A"
medints.d$alleleID[medints.d$alleleID==1] <- "B"
rownames(medints.d) <- paste(medints.d$probeset_id, medints.d$alleleID, sep="_")
medints.d$probeset_id <- NULL
medints.d$alleleID <- NULL
print("writing data")
write.table(medints.d, file=file.path(dataDir, "median.summarized.intensities.tab"), quote=F, sep="\t", row.names=T, col.names=NA)



