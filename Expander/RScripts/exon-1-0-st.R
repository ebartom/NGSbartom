# aroma.affymetrix

# Path: rawData/experiment_name/HuEx-1_0-st-v2/
#(downloaded HuEx-1_0-st-v2,coreR2,A20070914,EP.cdf - the CDF defining "core" units.)

aromaTest<- function(x)
  {
    if (!require(x,character.only = TRUE))
    {
	source("http://www.aroma-project.org/hbLite.R")
	hbInstall("aroma.affymetrix")
    }
  }
cmd_args = commandArgs();
if (cmd_args[3]=="false"){
	if(!require("aroma.affymetrix",character.only = TRUE)){
		stop("Missing aroma.affymetrix package")
	}
}					
aromaTest("aroma.affymetrix")

params <- as.matrix(read.table("aromaParameters.txt")) # adds a double quote to each cell read

cdfPath <- params[1,2]
cdfName <- params[2,2]
myExperimentName <- params[3,2]
resultFilePath <- params[4,2]


# library(affy)
library(aroma.affymetrix)
verbose <- Arguments$getVerbose(-8, timestamp=TRUE)	# Coerces to Verbose object. (Verbose: Class to writing verbose messages to a connection or file)




#Getting annotation data files
#cdf <- AffymetrixCdfFile$byChipType(chipType, tags=myTags)
cdf <- AffymetrixCdfFile$fromFile(cdfName, path=cdfPath)

# Defining CEL set
cs <- AffymetrixCelSet$byName(myExperimentName, cdf=cdf) 

# Background Adjustment and Normalization
# bc <- RmaBackgroundCorrection(cs, tags=myTags)
bc <- RmaBackgroundCorrection(cs, tag=cdfName)		# or should it be tags and not tag? since this name seems like a few tags
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

# Remove columns 2, 3, 4, 5 (corresponding to groupName, unit, group, cell)
trFit <- trFit[,c(1, 6:ncol(trFit))]

write.table(trFit, file = resultFilePath, quote = FALSE, sep = '\t', row.names = FALSE)








# source("G:\\cseagull\\R_info\\aroma.affymetrix\\exon-1-0-st.R")






