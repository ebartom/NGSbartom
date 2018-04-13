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
	if(!require("affy",character.only = TRUE)||!require("gcrma",character.only = TRUE)){
		stop("Missing affy and gcrma package")
	}
}
#print(cmd_args)
pkgTest("affy")
pkgTest("gcrma")
writeLines("Reading arguments...", con = stdout(), sep = "\n")
#args = commandArgs(TRUE)

print(length(cmd_args))


if( length(cmd_args) > 3 ) {
	celDirPath = cmd_args[4]
	print(celDirPath)
	setwd(celDirPath)
}

params <- as.matrix(read.table("celPreprocParameters.txt")) # adds a double quote to each cell read
preprocMethod <- params[1,2]			## rma or gcrma
IsSpecificCdf <- (params[2,2] == "T") 	## bool
CdfPackageName <- ""

if(IsSpecificCdf) {
	CdfPackageName <- params[3,2]
}

print(preprocMethod)
print(IsSpecificCdf)
print(CdfPackageName)

writeLines("Opening affy package...", con = stdout(), sep = "\n")
library(affy)

if(preprocMethod == "gcrma"){
	library(gcrma)
}

if(!IsSpecificCdf){
	Data<-ReadAffy()
} else {	# make sure to put } and else on the same line! 
	Data<-ReadAffy(cdfname = CdfPackageName)
	print("User specified a cdfname\n")
}	

writeLines("Performing RMA...", con = stdout(), sep = "\n")

if(preprocMethod == "rma"){
	print("Using RMA as preprocessing method\n")
	eset<-rma(Data)
} else if(preprocMethod == "gcrma"){
	print("Using GCRMA as preprocessing method\n")
	eset<-gcrma(Data)
} else {
	print("** ERROR unrecognized preprocessing method ** \n")
}


writeLines("Saving matrix...", con = stdout(), sep = "\n")
write.exprs(eset, file="unlikelyNameForRes.txt")


writeLines("Done!", con = stdout(), sep = "\n")
