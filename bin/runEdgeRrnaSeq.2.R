
args <- commandArgs()

assembly <-sub('--assembly=', '', args[grep('--assembly=', args)])
outputDirectory <-sub('--outputDirectory=', '',args[grep('--outputDirectory=',args)])
countFile <-sub('--countFile=', '', args[grep('--countFile=', args)])
numCores <- as.integer(sub('--numCores=', '', args[grep('--numCores=',args)]))
comparisonFile <-sub('--comparisonFile=', '', args[grep('--comparisonFile=',args)])
runMDS <-sub('--runMDS=', '',args[grep('--runMDS=',args)])
filterOff <-sub('--filterOff=', '',args[grep('--filterOff=',args)])

# If runMDS isn't specified, default is to run it.
if (identical(runMDS,character(0))){
   runMDS<-1
}

# If filterOff isn't specified, default is to filter genes by counts.
if (identical(filterOff,character(0))){
    filterOff<-0
} else {
    print("Genes not being filtered for minimum number of counts.");
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
if (rownames(counts)[1] == "1"){
    print("Problem with rownames, re-adjusting.")
    newcounts<-data.frame(counts[,2:length(colnames(counts))])
    rownames(newcounts) <- counts[,1]
    print(head(newcounts))
    counts <- newcounts
}
dim(counts)

hostMart <- ""

if (assembly == "hg19") {
   organismStr <- "hsapiens"
   hostMart <- "feb2014.archive.ensembl.org"
   attributes <- c('ensembl_gene_id','gene_biotype','external_gene_id','description','ensembl_transcript_id')
}
if (assembly == "hg38") {
   organismStr <- "hsapiens"
   hostMart <- "feb2014.archive.ensembl.org"
   attributes <- c('ensembl_gene_id','gene_biotype','external_gene_id','description','ensembl_transcript_id')
}
if (assembly == "hg38.mp") {
   organismStr <- "hsapiens"
   hostMart <- "feb2014.archive.ensembl.org"
   attributes <- c('ensembl_gene_id','gene_biotype','external_gene_id','description','ensembl_transcript_id')
}
if (assembly == "mm9") {
    organismStr <- "mmusculus"
    hostMart <- "feb2014.archive.ensembl.org"
    attributes <- c('ensembl_gene_id','gene_biotype','external_gene_id','description','ensembl_transcript_id')
}
if (assembly == "mm10") {
    organismStr <- "mmusculus"
    hostMart <- "feb2014.archive.ensembl.org"
    attributes <- c('ensembl_gene_id','gene_biotype','external_gene_id','description','ensembl_transcript_id')
}
if (assembly == "rn6") {
    organismStr <- "rnorvegicus"
    hostMart <- "mar2016.archive.ensembl.org"
    attributes <- c('ensembl_gene_id','gene_biotype','external_gene_name','description','ensembl_transcript_id')
}
if (assembly == "grcm38") {
   organismStr <- "mmusculus"
   hostMart <- "dec2016.archive.ensembl.org"
   attributes <- c('ensembl_gene_id','gene_biotype','external_gene_name','description','ensembl_transcript_id')
}
if (assembly == "sacCer3") {
   organismStr <- "scerevisiae"
   hostMart <- "feb2014.archive.ensembl.org"
   attributes <- c('ensembl_gene_id','gene_biotype','external_gene_id','description','ensembl_transcript_id')
}
if (assembly == "dm3") {
   organismStr <- "dmelanogaster"
   hostMart <- "feb2014.archive.ensembl.org"
   attributes <- c('ensembl_gene_id','gene_biotype','external_gene_id','description','ensembl_transcript_id')
}
organismStr

listMarts(host=hostMart)

bm <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host=hostMart)
ds <- listDatasets(bm)
dataset <- ds[grep(organismStr,ds$description),]$dataset
dataset

print("Loading biomart")
bm <- useDataset(paste(organismStr,"_gene_ensembl",sep=""), mart=bm)
bm
#listAttributes(bm)
#anno <- getBM(mart=bm, attributes=c('ensembl_gene_id','gene_biotype','external_gene_id','description','ensembl_transcript_id'))
anno <- getBM(mart=bm, attributes=attributes)

##----------MDS plot--------------------#

geneData <- data.frame(counts[,1:5])
print ("Gene data from counts table")
colnames(geneData)[4] <- "width"
print(head(geneData))
countData<-data.frame(counts[,6:length(colnames(counts))])
print ("Counts data from counts table")
print(head(countData))

method <- gsub(".all.counts.txt","",countFile)
method <- gsub("^.*\\/","",method)

print (runMDS)
if (runMDS==1){
   print("Generating MDS plot.")
   newnames <- counts[,6:length(colnames(counts))]
   newnames <- gsub("\\.\\d+$","",newnames)
   if(filterOff==0){
#       countData <- na.omit(countData)
       dge <- DGEList(countData[rowSums(cpm(countData)>1)>=2,]) #filter low counts & calc lib.sizes
   } else {
       print ("Not filtering out low counts.")
       dge <- DGEList(countData)      #calc lib.sizes
   }
   dge <- calcNormFactors(dge)  #calcs normalization factors based on lib.size (RNA composition effect)
   dge <- estimateCommonDisp(dge,verbose=T)
   dge <- estimateTagwiseDisp(dge)
   print("printing Normalized counts")
   dge$pseudo.counts
   method <- gsub(".all.counts.txt","",countFile)
   method <- gsub("^.*\\/","",method)
   write.table(dge$pseudo.counts, file=paste(outputDirectory,method,".normCounts.txt",sep=""), row.names = TRUE, col.names = TRUE, sep = "\t" )
   print ("Calculating and printing FPKM table")
   RPKM <- rpkm(countData[,],gene.length=geneData$width)
   head(RPKM)
   write.table(RPKM, file=paste(outputDirectory,method,".fpkms.txt",sep=""), row.names = TRUE, col.names = TRUE, sep = "\t" )
   pdfname <-paste(outputDirectory,method,".MDS.pdf",sep="")
   print(pdfname)
   pdf(pdfname)
#   plotMDS(dge, method="bcv")
   plotMDS(dge)
   dev.off()
}

##-----------gene expression--------------#

runEdgeR <- function(data,comparison){
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
    group2 <- data.frame(subset(data,select=group2))
    print (colnames(group1))
    print (colnames(group2))
    subdata <-cbind(group1,group2)
#    print(head(subdata))
    setup<- c(rep("group1",each=length(group1)),rep("group2",each=length(group2)))
    print(setup)
    grp <- factor(setup)
    grp <- relevel(grp,ref="group1")
    design <- model.matrix(~grp)
    colnames(design)<-levels(grp)
    print(head(subdata))
    if (filterOff == 1){
        dge <- DGEList(subdata, group=grp) # calc lib.sizes
    }else{
        dge <- DGEList(subdata[rowSums(cpm(subdata)>1)>=2,], group=grp) #filter low counts & calc lib.sizes
    }
                                        #    dge <- DGEList(subdata, group=grp) #calc lib.sizes
    print("Effective library size for each sample in comparison")
    print(dge$samples)
    dge <- calcNormFactors(dge)  #calcs normalization factors based on lib.size (RNA composition effect)
    print("Normalization factors for each sample in comparison, estimated with TMM method.")
    print(dge$samples)
    dge <- estimateGLMCommonDisp(dge,design)  #Estimates a common negative binomial dispersion
    dge <- estimateGLMTrendedDisp(dge,design) #Estimates the abundance-dispersion trend by Cox-Reid approximate profile likelihood
    dge <- estimateGLMTagwiseDisp(dge,design) #Compute an empirical Bayes estimate of the negative binomial dispersion parameter for each transcript, with expression
    fit <- glmFit(dge,design)    #Fit a negative binomial generalized log-linear model to counts
    lrt <- glmLRT(fit,coef="group2") 
    foo <- rownames(lrt$table)                #List of differentially expressed genes
    stopifnot(rownames(dge$counts)==rownames(geneData[foo,]))
#    RPKM <- rpkm(subdata[foo,],gene.length=geneData[foo,]$width)
    CPM <- cpm(subdata)[foo,]
    stopifnot(rownames(lrt$table)==rownames(CPM))
#    df <- cbind(as.data.frame(geneData[foo,]),RPKM,lrt$table,adj.p=p.adjust(lrt$table$PValue,method="BH"))
#    print(topTags(lrt))
    # Changed RPKM to CPM on July 13, 2016
    df <- cbind(as.data.frame(geneData[foo,]),CPM,lrt$table,adj.p=p.adjust(lrt$table$PValue,method="BH"))
    if (identical(method,"rsem")){
        iv <- match(rownames(df),anno$ensembl_transcript_id)
            df$gene <- anno[iv,'external_gene_id']
    } else {
        iv <- match(rownames(df),anno$ensembl_gene_id)
            df$gene <- anno[iv,'external_gene_id']
    }
    if ((identical(assembly,"hg38.mp")) || (identical(assembly,"rn6"))) {
        head(anno)				
        iv <- match(rownames(df),anno$external_gene_id)
            df$gene <- anno[iv,'ensembl_gene_id']
    }
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
   print("Comparisons file:")
   print(comparisonFile)
   comparisons <- read.table(file=comparisonFile,sep=",",as.is=T,header=T,row.names=1)
   print(comparisons)
   print(dim(comparisons))
   colnames(comparisons) <- gsub("\\_","\\.",colnames(comparisons))
   colnames(countData) <- gsub("\\_","\\.",colnames(countData))
   countDataSampleOrder <- colnames(countData)
   print("SampleOrder:")
   print(countDataSampleOrder)
   print ("ComparisonSampleOrder:")
   print(colnames(comparisons))
   print ("Reorder comparisons based on counts column order")
   comparisons <- comparisons[,countDataSampleOrder]
   print(countDataSampleOrder)
   print(colnames(comparisons))
   print(rownames(comparisons))
   for(c in 1:length(rownames(comparisons))){
     	comparison <- (paste(c(comparisons[c,1:length(colnames(comparisons))])))
    	print("Running EdgeR comparison")
    	print(paste(c(comparisons[c,1:length(colnames(comparisons))])))
    	print(rownames(comparisons)[c])
    	comp <-runEdgeR(data=countData,comparison)
        method <- gsub(".all.counts.txt","",countFile)
        method <- gsub("^.*\\/","",method)
   	label1<- paste(outputDirectory,paste(rownames(comparisons)[c],method,"edgeR.df.rda",sep="."),sep="")
	label2<- paste(outputDirectory,paste(rownames(comparisons)[c],method,"edgeR.txt",sep="."),sep="")
	save(comp,file=label1)
   	write.table(comp,file=label2,sep="\t",col.names=NA)
   }
}

