args <- commandArgs()

assembly <-sub('--assembly=', '', args[grep('--assembly=', args)])
outputDirectory <-sub('--outputDirectory=', '',args[grep('--outputDirectory=',args)])
gnModelFile <-sub('--gnModelFile=', '', args[grep('--gnModelFile=', args)])
countFile <-sub('--countFile=', '', args[grep('--countFile=', args)])
#sample <- sub('--sample=', '', args[grep('--sample=', args)])
#sample <- sub('$', '.bam', sample)
numCores <- as.integer(sub('--numCores=', '', args[grep('--numCores=',args)]))
comparisonFile <-sub('--comparisonFile=', '', args[grep('--comparisonFile=',args)])
runMDS <-sub('--runMDS=', '',args[grep('--runMDS=',args)])

# If runMDS isn't specified, default is to run it.
if (identical(runMDS,character(0))){
   runMDS<-1
}

# If comparison file isn't specified, default is to skip it.
runComp<-1
if (identical(comparisonFile,character(0))){
   runComp<-0
}

# If output directory isn't specified, default is current directory
if (identical(outputDirectory,character(0))){
   outputDirectory<-""
} else { 
  outputDirectory<-paste(outputDirectory,"/",sep="")
}

library(edgeR)
library(GenomicRanges)
library(biomaRt)

##----------load counts --------#
##load count data
print("Loading count data")
print(countFile)

if(grepl('rda',countFile)){
   counts <- get(load(file=countFile))
}
if(grepl('txt',countFile)){
    counts <- read.delim(file=countFile,header=TRUE,sep="\t")
}
head(counts)
dim(counts)

print("Loading gene model data")
gnModel <- get(load(gnModelFile))
print("Peek at gnModel")
head(gnModel)

hostMart <- ""

if (assembly == "hg19") {
   organismStr <- "hsapiens"
   hostMart <- "feb2014.archive.ensembl.org"
}
if (assembly == "mm9") {
   organismStr <- "mmusculus"
   hostMart <- "feb2014.archive.ensembl.org"
}
if (assembly == "mm10") {
   organismStr <- "mmusculus"
   hostMart <- "feb2014.archive.ensembl.org"
}
if (assembly == "sacCer3") {
   organismStr <- "scerevisiae"
   hostMart <- "feb2014.archive.ensembl.org"
}
if (assembly == "dm3") {
   organismStr <- "dmelanogaster"
   hostMart <- "feb2014.archive.ensembl.org"
}
organismStr

listMarts(host=hostMart)

bm <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host=hostMart)
ds <- listDatasets(bm)
dataset <- ds[grep(organismStr,ds$description),]$dataset
dataset

print("Loading biomart")
bm <- useDataset(paste(organismStr,"_gene_ensembl",sep=""), mart=bm)

anno <- getBM(mart=bm, attributes=c('ensembl_gene_id','gene_biotype','external_gene_id','description'))

##----------MDS plot--------------------#

allData<-data.frame(counts[,6:length(colnames(counts))])	

print (runMDS)
if (runMDS==1){
   print("Generating MDS plot.")
   newnames <- counts[,6:length(colnames(counts))]
   newnames <- gsub("\\.\\d+$","",newnames)
   dge <- DGEList(allData)      #filter low counts & calc lib.sizes
   dge <- calcNormFactors(dge)  #calcs normalization factors based on lib.size (RNA composition effect)
   dge <- estimateCommonDisp(dge,verbose=T)
   dge <- estimateTagwiseDisp(dge)
   print("printing Normalized counts")
   dge$pseudo.counts
   write.table(dge$pseudo.counts, file=paste(outputDirectory,"allData.normCounts.txt",sep=""), row.names = TRUE, col.names = TRUE, sep = "\t" )
   pdfname <-paste(outputDirectory,"allData.MDS.pdf",sep="")
   pdf(pdfname)
   plotMDS(dge, method="bcv")
   dev.off()
}

##-----------gene expression--------------#

runEdgeR <- function(data,comparison){
#	 print(head(data))
#	 print(comparison)
	 group1 <- c()
	 group2 <- c()
	 for(i in 1:length(comparison)){
	       value <- as.integer(comparison[i])
	       if (value == 1){
	       	  group1<-c(group1,(colnames(data)[i]))
		  }						
		  if (value == -1){				
	       	  group2<-c(group2,(colnames(data)[i]))
		  }	  		  
	 }
	 group1 <- data.frame(subset(data,select=group1))
#	 print(head(group1))
	 group2 <- data.frame(subset(data,select=group2))
#	 print(head(group2))
	 print (colnames(group1))
	 print (colnames(group2))
	 subdata <-cbind(group1,group2)
	 setup<- c(rep("group1",each=length(group1)),rep("group2",each=length(group2)))
	 print(setup)
	 grp <- factor(setup)
	 grp <- relevel(grp,ref="group1")
	 design <- model.matrix(~grp)
#	 print(design)
	 colnames(design)<-levels(grp)
#	 print ("Colnames of design")
#	 print(colnames(design))
	 print(head(subdata))
	 dge <- DGEList(subdata[rowSums(cpm(subdata)>1)>=2,], group=grp) #filter low counts & calc lib.sizes	
    	 dge <- calcNormFactors(dge)  #calcs normalization factors based on lib.size (RNA composition effect)
	 dge <- estimateGLMCommonDisp(dge,design)  #Estimates a common negative binomial dispersion
    	 dge <- estimateGLMTrendedDisp(dge,design) #Estimates the abundance-dispersion trend by Cox-Reid approximate profile likelihood
    	 dge <- estimateGLMTagwiseDisp(dge,design) #Compute an empirical Bayes estimate of the negative binomial dispersion parameter for each transcript, with expre
    	 fit <- glmFit(dge,design)    #Fit a negative binomial generalized log-linear model to counts
    	 lrt <- glmLRT(fit,coef="group2") 
#	 print("making annotated dataframe\n")
 	 foo <- rownames(lrt$table)                #List of differentially expressed genes
    	 stopifnot(rownames(lrt$table)==rownames(gnModel[foo,]))
    	 stopifnot(rownames(dge$counts)==rownames(gnModel[foo,]))
    	 RPKM <- rpkm(subdata[foo,],gene.length=gnModel[foo,]$width)
    	 stopifnot(rownames(lrt$table)==rownames(RPKM))
    	 df <- cbind(as.data.frame(gnModel[foo,]),RPKM,lrt$table,adj.p=p.adjust(lrt$table$PValue,method="BH"))
    	 iv <- match(rownames(df),anno$ensembl_gene_id)
    	 df$gene <- anno[iv,'external_gene_id']
    	 up <- df$adj.p < 0.01 & df$logFC > 0
    	 cat("adding flags\n")
    	 flag <- rep(0, nrow(df))
    	 flag[up] <- 1
    	 df$up <- flag
    	 dn <- df$adj.p < 0.01 & df$logFC < 0
    	 flag1 <- rep(0, nrow(df))    
    	 flag1[dn] <- 1
    	 df$dn <- flag1
    	 df$biotype <- anno[iv,'gene_biotype']
    	 df$description <- anno[iv,'description']
    	 df
}

if (runComp == 1) {
   print(comparisonFile)
   comparisons <- read.table(file=comparisonFile,sep=",",as.is=T,header=T,row.names=1)
   print(comparisons)
   print(dim(comparisons))
   print("Comparisons file:")
#   print(comparisons[0,])
#   print(colnames(comparisons))
   colnames(comparisons) <- gsub("\\_","\\.",colnames(comparisons))
   colnames(allData) <- gsub("\\_","\\.",colnames(allData))
   allDataSampleOrder <- colnames(allData)
   print("SampleOrder:")
   print(allDataSampleOrder)
#   rownames(comparisons) <- comparisons$Comparisons
#   print(rownames(comparisons))
   print ("ComparisonSampleOrder:")
   print(colnames(comparisons))
   print ("Reorder comparisons based on counts column order")
   comparisons <- comparisons[,allDataSampleOrder]
   print(allDataSampleOrder)
   print(colnames(comparisons))
   print(rownames(comparisons))
   for(c in 1:length(rownames(comparisons))){
     	comparison <- (paste(c(comparisons[c,1:length(colnames(comparisons))])))
    	print("Running EdgeR comparison")
    	print(paste(c(comparisons[c,1:length(colnames(comparisons))])))
    	print(rownames(comparisons)[c])
    	comp <-runEdgeR(data=allData,comparison)
   	label1<- paste(outputDirectory,paste(rownames(comparisons)[c],"edgeR.df.rda",sep="."),sep="")
	label2<- paste(outputDirectory,paste(rownames(comparisons)[c],"edgeR.txt",sep="."),sep="")
	save(comp,file=label1)
   	write.table(comp,file=label2,sep="\t",col.names=NA)
   }
}

