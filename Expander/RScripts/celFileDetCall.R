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
	if(!require("affy",character.only = TRUE)){
		stop("Missing affy package")
	}
}	
pkgTest("affy")
library(affy)
Data<-ReadAffy()
eset<-rma(Data)
Calls <- mas5calls(Data)
exprsMat <- exprs(eset)
callMat <- exprs(Calls)

combMat <- matrix(rbind(exprsMat,callMat),nrow(callMat))
colnames(combMat) <- rbind(colnames(exprsMat),colnames(callMat))
rownames(combMat) <- rownames(exprsMat)

write.table(combMat, file="unlikelyNameForRes.txt", quote = FALSE, sep = "\t", col.names = NA)
