# Gene 1.0 ST array analysis
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

library(aroma.affymetrix)
verbose <- Arguments$getVerbose(-8, timestamp=TRUE)

# Getting annotation data files
# cdf <- AffymetrixCdfFile$byChipType(chipType, tags="r3")
cdf <- AffymetrixCdfFile$fromFile(cdfName, path=cdfPath)

# Defining CEL set
cs <- AffymetrixCelSet$byName(myExperimentName, cdf=cdf)

# Background Adjustment and Normalization
bc <- RmaBackgroundCorrection(cs, tag=cdfName)
csBC <- process(bc,verbose=verbose)
qn <- QuantileNormalization(csBC, typesToUpdate="pm")
csN <- process(qn, verbose=verbose)

# Summarization
plm <- RmaPlm(csN)
fit(plm, verbose=verbose) 

# To extract the estimates, you can use 'extractDataFrame' on the ChipEffectSet object that corresponds to the plm object:
ces <- getChipEffectSet(plm)
gExprs <- extractDataFrame(ces, units=NULL, addNames=TRUE)

# Remove columns 2, 3, 4, 5 (corresponding to groupName, unit, group, cell)
gExprs <- gExprs[,c(1, 6:ncol(gExprs))]


write.table(gExprs, file = resultFilePath, quote = FALSE, sep = '\t', row.names = FALSE)


# source("G:\\cseagull\\R_info\\aroma.affymetrix\\gene-1-0-st.R")








