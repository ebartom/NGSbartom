pkgTest <- function(x)
  {
    if (!require(x,character.only = TRUE))
    {
	source("http://bioconductor.org/biocLite.R")
	biocLite("eisa")
    }
  }
cmd_args = commandArgs();
if (cmd_args[3]=="false"){
	if(!require("isa2",character.only = TRUE)){
		stop("Missing eisa package")
	}
}	
pkgTest("isa2")
library("isa2")
getwd()
data <- as.matrix(read.table("isaInputFile.txt", header=T ,row.names = 1, sep = '\t'))
isa.result=isa(data)
rows = isa.result$rows
cols = isa.result$columns
scores = isa.result$seeddata$rob

write.table(scores, file = "isa_scores.txt", quote = FALSE, sep = '\t', row.names = TRUE, col.names = FALSE)
write.table(rows, file = "isa_rows.txt", quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
write.table(cols, file = "isa_cols.txt", quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
