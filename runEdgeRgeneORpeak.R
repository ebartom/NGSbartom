args <- commandArgs()

assembly <-sub('--assembly=', '', args[grep('--assembly=', args)])
outputDirectory <-sub('--outputDirectory=', '',args[grep('--outputDirectory=',args)])
txdbfile <-sub('--txdbfile=','',args[grep('--txdbfile=',args)])
countFile <-sub('--countFile=', '', args[grep('--countFile=', args)])
numCores <- as.integer(sub('--numCores=', '', args[grep('--numCores=',args)]))
comparisonFile <-sub('--comparisonFile=', '', args[grep('--comparisonFile=',args)])
runMDS <-sub('--runMDS=', '',args[grep('--runMDS=',args)])
type <-sub('--type=','',args[grep('--type=',args)])
hostMart <-sub('--hostMart=','',args[grep('--hostMart=',args)])

# If hostMart isn't specified, use the current default at ensembl.org
if (identical(hostMart,character(0))){
   hostMart<-"ensembl.org"
}

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

# If type isn't specified, default is RNA
if (identical(type,character(0))){
    type<- "RNA"
}

library(edgeR)
library(GenomicRanges)
library(biomaRt)
library(GenomicFeatures)

##----------load counts --------#
##load count data
print("Loading count data")
print(countFile)

if(grepl('rda',countFile)){
   counts <- get(load(file=countFile))
}
if(grepl('txt',countFile)){
    counts <- read.delim(file=countFile,header=TRUE,sep="\t")
    if (type == "chipseq"){
        rownames(counts) <- counts$name
    }
}
head(counts)
dim(counts)

if (assembly == "hg19") {
   organismStr <- "hsapiens"
}
if (assembly == "mm9") {
   organismStr <- "mmusculus"
}
if (assembly == "mm10") {
   organismStr <- "mmusculus"
}
if (assembly == "sacCer3") {
   organismStr <- "scerevisiae"
}
if (assembly == "dm3") {
   organismStr <- "dmelanogaster"
}
organismStr

print(assembly)
if (assembly == "hg19") { organismStr <- "Hsapiens" }
if ((assembly == "mm9") || (assembly == "mm10")) { organismStr <- "Mmusculus" }
if (assembly == "sacCer3") { organismStr <- "Scerevisiae"}
if (assembly == "dm3") { organismStr <- "Dmelanogaster"}

assemblyLibrary <- paste("BSgenome.", organismStr, ".UCSC.", assembly, sep="")
print(assemblyLibrary)
library(assemblyLibrary,character.only=TRUE)

if (assembly == "hg19") { organism <- Hsapiens }
if ((assembly == "mm9") || (assembly == "mm10")) { organism <- Mmusculus }
if (assembly == "sacCer3") { organism <- Scerevisiae}
if (assembly == "dm3") { organism <- Dmelanogaster}

if (type == "RNA"){
    txdb <- loadDb(txdbfile)
    txdb
    
    seqlevels(txdb,force=TRUE) <- seqlevels(txdb)[grep("_|\\d+.1$",seqlevels(txdb), invert=TRUE)]
    seqlevels(txdb) <- sub("^","chr", seqlevels(txdb))
    seqlevels(txdb) <- sub("MT","M", seqlevels(txdb))
    seqlevels(txdb) <- sub("Mito","M", seqlevels(txdb))
    
    gnModel <- transcriptsBy(txdb, 'gene')
    gnModel <- unlist(range(gnModel))## Gets the genomic region convered by transcripts
    seqinfo(gnModel) <- seqinfo(organism)[seqlevels(gnModel)]
    gnModel<-as.data.frame(gnModel)
    print("Peek at gnModel")
    head(gnModel)

    listMarts(host=hostMart)
    
    print("Loading biomart")
    bm <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host=hostMart)
    ds <- listDatasets(bm)
    dataset = paste(tolower(organismStr),"gene_ensembl",sep="_")
    bm <- useDataset(dataset=dataset, mart=bm)
    
    anno <- getBM(mart=bm, attributes=c('ensembl_gene_id','gene_biotype','external_gene_name','description'))
}

##----------MDS plot--------------------#

if (type == "RNA"){
    allData<-data.frame(counts[,6:length(colnames(counts))])
} else if (type == "chipseq"){
    allData<-data.frame(counts[,5:length(colnames(counts))])
}
head(allData)
print(type)

print (runMDS)
if (runMDS==1){
    print("Generating MDS plot.")
    if (type == "RNA"){
        print ("Skipping first five columns of metadata")
        newnames <- counts[,6:length(colnames(counts))]
    } else if (type == "chipseq"){
        print ("Skipping first four columns of metadata")
        newnames <- counts[,5:length(colnames(counts))]
    }
    newnames <- gsub("\\.\\d+$","",newnames)
    head(allData)
    head(newnames)
   dge <- DGEList(allData)      #filter low counts & calc lib.sizes
    dge <- calcNormFactors(dge)  #calcs normalization factors based on lib.size (RNA composition effect)
    prefix <- print(basename(countFile))
    prefix <- sub(".counts.txt","",prefix)
   pdfname <-paste(outputDirectory,prefix,".MDS.pdf",sep="")
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
    	 dge <- calcNormFactors(dge) #calcs normalization factors based on lib.size (RNA composition effect)
	 dge <- estimateGLMCommonDisp(dge,design)  #Estimates a common negative binomial dispersion
    	 dge <- estimateGLMTrendedDisp(dge,design) #Estimates the abundance-dispersion trend by Cox-Reid approximate profile likelihood
    	 dge <- estimateGLMTagwiseDisp(dge,design) #Compute an empirical Bayes estimate of the negative binomial dispersion parameter for each transcript, with expre
    	 fit <- glmFit(dge,design)    #Fit a negative binomial generalized log-linear model to counts
    	 lrt <- glmLRT(fit,coef="group2") 
 	 foo <- rownames(lrt$table)                #List of differentially expressed genes/peaks
         if (type == "RNA"){
             print("making annotated dataframe\n")
             stopifnot(rownames(lrt$table)==rownames(gnModel[foo,]))
             stopifnot(rownames(dge$counts)==rownames(gnModel[foo,]))
             RPKM <- rpkm(subdata[foo,],gene.length=gnModel[foo,]$width)
             stopifnot(rownames(lrt$table)==rownames(RPKM))
             df <- cbind(as.data.frame(gnModel[foo,]),RPKM,lrt$table,adj.p=p.adjust(lrt$table$PValue,method="BH"))
             iv <- match(rownames(df),anno$ensembl_gene_id)
             df$gene <- anno[iv,'external_gene_id']
             df$biotype <- anno[iv,'gene_biotype']
             df$description <- anno[iv,'description']
         } else if (type == "chipseq"){
             print("making annotated dataframe\n")
             df <- cbind(counts[rownames(lrt$table),1:3],subdata[rownames(lrt$table),],lrt$table,adj.p=p.adjust(lrt$table$PValue,method="BH"))
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
    	 df
}

if (runComp == 1) {
   print(comparisonFile)
#   comparisons <- read.table(file=comparisonFile,sep=",",as.is=T,header=T,row.names=1)
   comparisons <- read.table(file=comparisonFile,sep=",",header=T,row.names=1)
   print(comparisons)
   colnames(comparisons) <- gsub("\\_","\\.",colnames(comparisons))
   colnames(allData) <- gsub("\\_","\\.",colnames(allData))
#   print(colnames(comparisons))
#   print(colnames(allData))
   allDataSampleOrder <- colnames(allData)
   print(allDataSampleOrder)
   print("here")
#   comparisons <- comparisons[,allDataSampleOrder]
   print(allDataSampleOrder)
   print(colnames(comparisons))
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

