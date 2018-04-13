# This script gets rsem.genes.results files, combine all the expected_counts columns
# in each file into one matrix, round the values to the nearest integer value
# and runs differential expression on two groups
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
	if(!require("edgeR",character.only = TRUE)){
		stop("Missing edgeR package")
	}
}
if (cmd_args[3]=="false"){
	if(!require("limma",character.only = TRUE)){
		stop("Missing limma package")
	}
}
#install required packages
pkgTest("limma")		
pkgTest("edgeR")

print(Sys.time())
library(edgeR)
library(limma)

# Top tags for tagwise analysis
options( digits = 3 ) # print only 3 digits

#read parameters
params <- as.matrix(read.table("rnaSeqParameters.txt"))

# Get parameters: to remove the double quote, use 'as.numeric'
inFDR = as.numeric(params[1])         # FDR percent. Will be multiplied by 2 if only one-tailed is asked (since here we do two-tailed)
groups = params[2]      # A string containing the number of conditions in each group
dispersion = as.numeric(params[3])  # Dispersion type ( Tagwise = 0, Common = 1, Poisson = 2
priorN = as.numeric(params[4])		# Prior n (only relevant when using Tagwise dispersion)

print(inFDR)

# Parse groups
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
print(y)#debug

# Dispersion type:
TagWise = 0
Common = 1
Poisson = 2

if(length(groupsList) != 2) { stop("Number of groups is not 2") }

# Read and manipulate input file
raw.data <- read.table( file = "rnaSeqInputFile.txt" , header = TRUE ,sep="\t",row.names=1)
head( raw.data )


counts <- raw.data
#prints some information on the data
print(paste("counts:",dim( counts )))
print(paste("Library size:",colSums( counts ))) # Library Sizes
print(paste("Library sizes in millions of reads:",colSums( counts ) / 1e06)) # Library Sizes in millions of reads
print("Number of genes with low counts")
print(table( rowSums( counts ) )[ 1:30 ]) # Number of genes with low counts

#create edgeR object cds
groups = as.matrix(y)
dim(groups)
print(groups)#debug
counts = as.matrix(counts)

cds <- DGEList( counts , group = groups )


cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) > 1) >= 3, ]
cds <- calcNormFactors( cds )
# effective library sizes
print(paste("Effective libaray size: ",cds$samples$lib.size * cds$samples$norm.factors))

# estimate dispersion and run differential expression test
if(dispersion == Common){
cds <- estimateCommonDisp( cds )
#print(paste("Common dispersion: ",cds$common.dispersion))
diffexp <- exactTest( cds , pair = c( "2" , "1" ) )
# Store full topTags results table
resultsTbl <- topTags( diffexp , n = nrow( diffexp$table ) )$table
# Names/IDs of DE genes
de.genes <- rownames( resultsTbl )[ resultsTbl$FDR <= inFDR ]
resultsTbl<- subset(resultsTbl, FDR <= inFDR)
# Amount significant
print(paste("Amount of significant: ",length( de.genes )))
# Percentage of total genes
print(paste("Percentage of total genes: ",length( de.genes ) / nrow( resultsTbl) * 100))
}


if(dispersion == TagWise){
cds <- estimateTagwiseDisp( cds  )
#print(paste("TagWise dispersion: " ,cds$tagwise.dispersion ))
diffexp <- exactTest( cds  , pair = c( "2" , "1" ) )
# Store full topTags results table
resultsTbl <- topTags( diffexp , n = nrow( diffexp$table ) )$table
# Names/IDs of DE genes
de.genes <- rownames( resultsTbl )[ resultsTbl$FDR <= inFDR ]
resultsTbl<- subset(resultsTbl, FDR <= inFDR)
# Amount significant
print(paste("Amount of significant: ",length( de.genes )))
# Percentage of total genes
print(paste("Percentage of total genes: ",length( de.genes ) / nrow( resultsTbl) * 100))
}


if(dispersion == Poisson){
diffexp <- exactTest( cds , dispersion = 1e-06 , pair = c( "2" , "1" ) )
# Store full topTags results table
resultsTbl <- topTags( diffexp , n = nrow( diffexp$table ) )$table
# Names/IDs of DE genes
de.genes <- rownames( resultsTbl )[ resultsTbl$FDR <= inFDR ]
resultsTbl<- subset(resultsTbl, FDR <= inFDR)
# Amount significant
print(paste("Amount of significant: ",length( de.genes )))
# Percentage of total genes
print(paste("Percent-age of total genes: ",length( de.genes ) / nrow( resultsTbl) * 100))
}

# Output the results
# Change column names to be specific to the analysis, logConc and logFC are the same in both.
colnames( resultsTbl ) <- c( "logFC" , "logCPM" , "pValue" , "FDR" )
# Below provides the info to re-order the count matrix to be in line with the order of the results.
wh.rows <- match( de.genes , rownames( cds$counts ) )
print(wh.rows)
# Tagwise Results
if(dispersion == TagWise){
combResults <- cbind( resultsTbl,
"Disp" = cds$tagwise.dispersion[ wh.rows ] ,
"UpDown" = decideTestsDGE( diffexp , p.value = inFDR )[ wh.rows ]  )
}
if(dispersion == Common){
combResults <- cbind( resultsTbl,
"Disp" = cds$common.dispersion[ wh.rows ] ,
"UpDown" = decideTestsDGE( diffexp , p.value = inFDR )[ wh.rows ] )
}
if(dispersion == Poisson){
combResults <- cbind( resultsTbl,
"UpDown" = decideTestsDGE( diffexp , p.value = inFDR )[ wh.rows ] )
}

write.table(combResults,file = "rnaSeqOutputFile.txt", sep="\t",row.names = TRUE)

print(Sys.time())
