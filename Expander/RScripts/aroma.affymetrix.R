# aroma.affymetrix

# Path: rawData/experiment_name/HuEx-1_0-st-v2/
#(downloaded HuEx-1_0-st-v2,coreR2,A20070914,EP.cdf - the CDF defining "core" units.)
pkgTest <- function(x)
  {
    if (!require(x,character.only = TRUE))
    {
	source("http://bioconductor.org/biocLite.R")
	biocLite(x)
    }
  }
cmd_args = commandArgs();
if (cmd_args[3]=="false"){
	if(!require("affy",character.only = TRUE)||!require("aroma.affymetrix",character.only = TRUE)){
		stop("Missing affy and aroma.affymetrix package")
	}
}	
pkgTest("affy")

aromaTest<- function(x)
  {
    if (!require(x,character.only = TRUE))
    {
	source("http://www.aroma-project.org/hbLite.R")
	hbInstall("aroma.affymetrix")
    }
  }
		
aromaTest("aroma.affymetrix")

library(affy)
library(aroma.affymetrix)
verbose <- Arguments$getVerbose(-8, timestamp=TRUE)	# Coerces to Verbose object. (Verbose: Class to writing verbose messages to a connection or file)

#Getting annotation data files
chipType <- "HuEx-1_0-st-v2" 
cdf <- AffymetrixCdfFile$byChipType(chipType, tags="coreR2,A20070914,EP")

# Defining CEL set
cs <- AffymetrixCelSet$byName("experiment_name", cdf=cdf) 

# Background Adjustment and Normalization
bc <- RmaBackgroundCorrection(cs, tag="coreR2")
csBC <- process(bc,verbose=verbose)
qn <- QuantileNormalization(csBC, typesToUpdate="pm")
csN <- process(qn, verbose=verbose) # This will take approx 30-60s per array.  

# Summarization
# To fit a summary of the entire transcript (i.e. estimate the overall expression for the transcript, not exon-by-exon):
plmTr <- ExonRmaPlm(csN, mergeGroups=TRUE)
fit(plmTr, verbose=verbose)		# This will roughly take a few minutes per array if you are using the core probesets only.

# To extract the estimates (transcript or probeset) use either extractMatrix() or extractDataFrame() on the ChipEffectSet that corresponds to the plm object:
cesTr <- getChipEffectSet(plmTr)
trFit <- extractDataFrame(cesTr, units=NULL, addNames=TRUE) # this is for probeset

write.table(trFit, file = "aromaRes4.txt", quote = FALSE, sep = '\t', row.names = FALSE)








# possible parameters:
# 'experiment_name'
# name of cdf file that user downloaded >> we can parse it to get the tags...



# source("G:\\cseagull\\R_info\\aroma.affymetrix\\aroma.affymetrix.R")






