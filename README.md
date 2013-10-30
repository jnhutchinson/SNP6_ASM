#SNP6_ASM

Scripts for calling Allele specific methylation from Affymetrix SNP6.0 arrays.

##Requirements
Affymetrix SNP6.0 library (cdf) file - download from [Netaffx](http://www.affymetrix.com/analysis/index.affx) or from [here](http://dl.dropboxusercontent.com/u/4253254/resources/SNP6/GenomeWideSNP_6.Full.cdf)

##Usage

###Extract all probes raw intensities and median summarize
python extract.probe.intensities.py -h for help

optional arguments:
  -h, --help            show this help message and exit

required arguments:
  -c CELDIR, --CELdir CELDIR  
                        directory containing CEL files  
  -a APTDIR, --aptdir APTDIR  
                        directory containing apt-cel-extract binary  
  -d CDFFILE, --CDFfile CDFFILE  
                        full path to CDF file for AffySNP6 chip  

Make sure all your CEL file extensions are in UPPERCASE

Should return a tab delimited file with the median summarized intensities of all samples (one row per probeid, samples in columns)

##Notes

This is a work in progress. 

The initial R methods using aroma.affymetrix are too slow. Replaced them with the affy power tools to extract initial intensities and then median summarize in R.
 
Next steps are to properly integrate the quantile normalization steps.





Currently very slow, with the main slowdown happening at the reformat of the data structure after extracting the raw intensities. 
