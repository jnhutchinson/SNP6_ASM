
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
pkgTest("aroma.affymetrix")
require(reshape)
require(plyr)
require(aroma.affymetrix)


## ----variables-----------------------------------------------------------
# the base of the aroma.affymetrix directory structure
baseDir <- "/Volumes/2TB/backups/Broad_backup/pipeline" # change this to match your own basedirectory
if(!file.exists(file.path(baseDir, "rawData", "CELfiles"))){
  print("CELfile directory structure incorrect")
  q()
  }
if(!file.exists(file.path(baseDir, "annotationData", "chipTypes", "GenomeWideSNP_6.Full", "GenomeWideSNP_6.Full.cdf" ))){
  print("No Affy6 library (cdf) file or in incorrect directory")
  q()                                                                                                                       }
setwd(baseDir)


## ----annotations---------------------------------------------------------
# load in probe annotations, detailing allele identities of the "A" and "B" probesetIDs
A.B.probe.annots <- read.table(file.path(baseDir, "annotationData", "simplified.annot.A.B.hg18.tab"), header=T, sep="\t")
# reformat annotations and label
A.B.probe.annots.m <- melt(A.B.probe.annots, id.vars="SNP_ID")
names(A.B.probe.annots.m) <- c("probeset", "alleleID", "allele")

# load array library (locations of probes on chip etc.)
cdf <- AffymetrixCdfFile$byChipType("GenomeWideSNP_6.Full") 
# extract the indices of all SNP probes
SNPprobe.indices <- grep("SNP_A", getUnitNames(cdf))


## ----process-------------------------------------------------------------
# list CEL files
celfiles <- list.files(file.path(baseDir, "rawData", "CELfiles"), pattern="CEL|cel" )

# load CEL files and process one at a time
sapply(celfiles, function(celfile) {
  cf <- AffymetrixCelFile(file.path(baseDir, "rawData", "CELfiles", celfile), cdf=cdf)
  cs <- AffymetrixCelSet(files=list(cf))
  # get individual probesetID intensities
  y <- getUnitIntensities(cs, units=SNPprobe.indices[1:2], verbose=T)
  # reformat list of lists into a informatively labelled dataframe 
  intensities <- melt(y)
  names(intensities)[c(1:2, 5,6)] <- c("probe", "sample", "allele", "probeset")
  # take median value for all probes within a probesetID
  intensities.median <- aggregate(intensities$value, by=list(sample=intensities$sample, allele=intensities$allele, probeset=intensities$probeset), median)
  # merge in the AlleleIDs
  intensities.median <- merge(intensities.median, A.B.probe.annots.m, by=c("allele", "probeset"))
  intensities.median$probesetID <- paste(intensities.median$probeset, intensities.median$alleleID, sep="_") 
  
  data <- cast(intensities.median, probesetID ~ sample, value="x")
  write.table(data, file=file.path(baseDir, "rawData", "rawintensities", paste("raw", celfile, "tab", sep=".")))
  })


