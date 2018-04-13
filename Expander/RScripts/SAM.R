pkgTest <- function(x)
  {
    if (!require(x,character.only = TRUE))
    {
	install.packages("samr", repos = "http://cran.us.r-project.org")
	source("http://bioconductor.org/biocLite.R")
	biocLite("impute")
    }
  }

cmd_args = commandArgs();
if (cmd_args[3]=="false"){
	if(!require("samr",character.only = TRUE)){
		stop("Missing samr package")
	}
}		
pkgTest("samr")
print(Sys.time())
library(samr)
params <- as.matrix(read.table("samParameters.txt"))

### Get parameters: to remove the double quote, use 'as.numeric'
log2 = (params[1] == "T")        # boolean: is data given in log2? F or T
FDR = as.numeric(params[2])         # FDR percent. Will be multiplied by 2 if only one-tailed is asked (since here we do two-tailed)
minFC = as.numeric(params[3])       # The minimum fold change desired; should be >1; default is zero, meaning no fold change criterion is applied
groups = params[4]      # A string containing the number of conditions in each group
regulation = as.numeric(params[5])  # give up/down genes or both up and down

print("log2: "); print(log2);
print("FDR: "); print(FDR);
print("Regulation Type [0=up, 1=down, 2=diff]: "); print(regulation);
print(groups);


### Parse groups

groupsList = unlist(strsplit(groups, ","))
isFirst = 1
for (i in 1:length(groupsList)) { 
	if(isFirst) {
		y <- rep(i, groupsList[i])
		isFirst = 0
	}
	else
		y <- c(y, rep(i, groupsList[i]))
}

### Regulation type:
UP = 1
DOWN = 2
DIFF = 0

if(length(groupsList) == 2) { type = "Two class unpaired" } else { type = "Multiclass"; regulation = DIFF; }

print(type);
 


### We do here 2-tailed, so adjust the FDR cutoff accordingly:
if ((regulation == UP) || (regulation == DOWN)) {
	FDR = FDR*2
	print("FDR for 2-tailed test: "); print(FDR);
}

### Read and manipulate input file
	
data <- as.matrix(read.table("samInputFile.txt", header=T, ,row.names = 1, sep = '\t'))
x <- data

DataList <- list(x = x, y = y, geneid=as.character(1:nrow(x)), genenames=row.names(data),  logged2=log2)
samr.obj<-samr(DataList,  resp.type = type, nperms=100)

delta.table <- samr.compute.delta.table(samr.obj, min.foldchange=minFC)
tableFDR_rowNum <- (which(delta.table[,5] <= FDR))[1]
del <- delta.table[tableFDR_rowNum, 1]
siggenes.table<- samr.compute.siggenes.table(samr.obj, del, DataList, delta.table)

### Write output file
if(regulation == UP) {
	bigMat <- rbind(siggenes.table$genes.up)
	qvalColNum = ncol(bigMat );
	write.table(siggenes.table$genes.up[,c(2,qvalColNum,qvalColNum-1)], file = "samOutputFile.txt", quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE);
	#write.table(siggenes.table$genes.up[,2], file = "samOutputFile.txt", quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE);
	print("Outputing 'up' diff genes");
}
if(regulation == DOWN) {
	bigMat <- rbind(siggenes.table$genes.lo)
	qvalColNum = ncol(bigMat );
	write.table(siggenes.table$genes.lo[,c(2,qvalColNum,qvalColNum-1)], file = "samOutputFile.txt", quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE);
	#write.table(siggenes.table$genes.lo[,2], file = "samOutputFile.txt", quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE);
	print("Outputing 'down' diff genes");
}
if(regulation == DIFF) {
	bigMat <- rbind(siggenes.table$genes.up, siggenes.table$genes.lo) # combine the up regulated genes with the down regulated
	qvalColNum = ncol(bigMat)  # the q-value is at the last column of the table (the number of columns changes according to num of groups)
	sortedMat <- bigMat[order(bigMat[,qvalColNum]), ] # sort according to the 'q-value' column
	#write.table(sortedMat[,2], file = "samOutputFile.txt", quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
	write.table(sortedMat[,c(2,qvalColNum,qvalColNum-1)], file = "samOutputFile.txt", quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
	#write.table(bigMat, file = "bigSamOutputFile.txt", quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
	print("Outputing diff genes (up and down)");
}

print(Sys.time())


### Assumptions:
# groups parameter contains at least two groups.
# input file, parameters file and output file have fixed names.
