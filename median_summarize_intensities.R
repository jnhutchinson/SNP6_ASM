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
  # intCol <- grep(sub(".rawints.txt","", rawintFile), names(rawints)) # hardcoded to column 8 
  # subset to genotype probesets only
  rawints <- subset(rawints, probeset_type=="GenoType")
  # aggregate by probeset_id and allele and take me dian value for intensities 
  medints <- aggregate(rawints[,8], by=list(probeset_id=rawints$probeset_id, alleleID=rawints$block), median)
  medints <- medints[grep("SNP_A", medints$probeset_id),]
  medints$alleleID[medints$alleleID==0] <- "A"
  medints$alleleID[medints$alleleID==1] <- "B"
  rownames(medints) <- paste(medints$probeset_id, medints$alleleID, sep="_")
  medints$probeset_id <- NULL
  medints$alleleID <- NULL
  names(medints) <- sub("rawints.txt", "", rawintFile)
  write.table(medints, file=file.path(dataDir, sub("rawints.txt", "medsumints.txt", rawintFile)), quote=F, sep="\t", row.names=T, col.names=NA)
  print(paste("wrote median summarized data for ", rawintFile, sep=""))
  return(medints)
})

medints.all <- do.call(cbind, medints.l)
write.table(medints.all, file=file.path(dataDir, "median.summarized.intensities.tab"), quote=F, sep="\t", row.names=T, col.names=NA)



