# This script computes differential expressed genes from RNA-Seq count matrix 
# using Deseq2 package which is available under LGPL license

pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x, repos = "http://cran.us.r-project.org")
    source("http://bioconductor.org/biocLite.R")
    biocLite(x)
  }
}

cmd_args = commandArgs();
if (cmd_args[3]=="false"){
  if(!require("DESeq2",character.only = TRUE)){
    stop("Missing edgeR package")
  }
}


#debug

#setwd("D:/acgt_lab/Expander/RScripts")

#install required packages
pkgTest("DESeq2")

print(Sys.time())
library(DESeq2)

# Top tags for tagwise analysis
options( digits = 3 ) # print only 3 digits

#read parameters
params <- as.matrix(read.table("rnaSeqParameters.txt"))
# Get parameters: to remove the double quote, use 'as.numeric'
inAnnotFile = params[1]
inFDR = as.numeric(params[2]) # FDR percent. Will be multiplied by 2 if only one-tailed is asked (since here we do two-tailed)
varDiff = params[3]      # A string containing the variable name of comparison
inFirstLabel = params[4]      # A string containing the first label of comparison
inSecondLabel = params[5]      # A string containing the second label of comparison
adjustFDR = params[6]=="1"         # a boolean flag to afjust by FDR?
forbVars <- unlist(strsplit(params[7],split=" ")) # forbidden variables

print(inAnnotFile)
print(inFDR)
print(varDiff)
print(inFirstLabel)
print(inSecondLabel)
print(adjustFDR)
print(forbVars)

# Read and manipulate input file
raw.data <- read.table( file = "rnaSeqInputFile.txt" , header = TRUE ,sep="\t",row.names=1)
head( raw.data )
counts <- raw.data

#prints some information on the data
print(paste("counts:",dim( counts )))
print(paste("Library size:",colSums( counts ))) # Library Sizes
print(paste("Library sizes in millions of reads:",colSums( counts ) / 1e06)) # Library Sizes in millions of reads


#create Deseq2 object dds
countdata = as.matrix(counts)
annotData = read.table( file = inAnnotFile , header = TRUE ,sep="\t",row.names=1)
annotData <- annotData[, !names(annotData) %in% forbVars,drop = FALSE]
head(annotData)
#coldata = data.frame(group = as.factor(y), row.names = colnames(counts))
vars <- paste0(colnames(annotData))
fmla <- as.formula(paste("~ ", paste(vars, collapse= "+")))
ddsMat <- DESeqDataSetFromMatrix(countData = countdata, colData = annotData, design = fmla)
#ddsMat$group <- relevel(ddsMat$group, 2)

# remove rows with all zeros or low counts
dds <- ddsMat[ rowSums(counts(ddsMat)) > 1, ]

dds <- DESeq(dds)

#FDR
res.fdr <- results(dds, contrast = c(varDiff,inFirstLabel,inSecondLabel))
# remove NAs
res.fdr.noNA <- res.fdr[complete.cases(res.fdr),]
res.fdr.noNA <- res.fdr.noNA[ order(res.fdr.noNA$log2FoldChange, decreasing=TRUE), ]

#FDR filtering and sorting - strongest down regulation
resSig <- res.fdr.noNA
if(adjustFDR)
  resSig <- subset(res.fdr.noNA, padj < inFDR)


#head(resSig[ order(resSig$log2FoldChange), ])
#strongest up-regulation
#head(resSig[ order(resSig$log2FoldChange, decreasing=TRUE), ])
resSig <- resSig[ order(resSig$log2FoldChange, decreasing=TRUE), ]
#output the results
head(resSig)
write.table(resSig[,c(2,5,6)],file = "rnaSeqOutputFile.txt",quote = FALSE, sep="\t",row.names = TRUE)

print(Sys.time())















