args <- commandArgs()

assembly <-sub('--assembly=', '', args[grep('--assembly=', args)])
outputDirectory <-sub('--outputDirectory=', '',args[grep('--outputDirectory=',args)])
countFile <-sub('--countFile=', '', args[grep('--countFile=', args)])
numCores <- as.integer(sub('--numCores=', '', args[grep('--numCores=',args)]))
comparisonFile <-sub('--comparisonFile=', '', args[grep('--comparisonFile=',args)])
groupFile <-sub('--groupFile=', '', args[grep('--groupFile=',args)])
runMDS <-sub('--runMDS=', '',args[grep('--runMDS=',args)])
filterOff <-sub('--filterOff=', '',args[grep('--filterOff=',args)])

# Rscript /projects/p20742//tools/runEdgeRwithGroups.R --assembly=mm9 --countFile=/projects/b1025/etb/deqingMLL2/newResults/RNAseqMLL2wtVmut//RNAseqMLL2wtVmut/bam//htseq.all.counts.corr.txt --comparisonFile=/projects/b1025/etb/deqingMLL2/newResults/RNAseqMLL2wtVmut/MLL2repeat.comparisons.csv --groupFile=/projects/b1025/etb/deqingMLL2/newResults/RNAseqMLL2wtVmut/Xin.design.csv --numCores=10 --outputDirectory=/projects/b1025/etb/deqingMLL2/newResults/RNAseqMLL2wtVmut//RNAseqMLL2wtVmut/analysis --runMDS=0

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
dim(counts)

hostMart <- ""

if (assembly == "hg19") {
   organismStr <- "hsapiens"
   hostMart <- "feb2014.archive.ensembl.org"
}
if (assembly == "hg38") {
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

anno <- getBM(mart=bm, attributes=c('ensembl_gene_id','gene_biotype','external_gene_id','description','ensembl_transcript_id'))

#geneData <- data.frame(counts[,1:5])
#colnames(geneData)[4] <- "width"
#print ("Gene data from counts table")
#print(head(geneData))
countData<-data.frame(counts[,2:length(colnames(counts))])
rownames(countData) <- counts$Gene
print ("Counts data from counts table")
print(head(countData))

method <- gsub(".all.counts.txt","",countFile)
method <- gsub("^.*\\/","",method)


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
    dge <- DGEList(subdata[rowSums(cpm(subdata)>1)>=2,], group=grp) #filter low counts & calc lib.sizes
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
    if(filterOff==0){
        dge <- DGEList(countData[rowSums(cpm(countData)>1)>=2,]) #filter low counts & calc lib.sizes
    } else {
        dge <- DGEList(countData)      #calc lib.sizes
    }
    print(dim(countData))
    print("Groups file:")
    print(groupFile)
    groups <- read.table(file=groupFile,sep=",",as.is=T,header=T,row.names=1)
    print(groups)
    genotype1 <- factor(groups$Genotype)
    celltype1 <- factor(groups$Celltype)
    print("genotype")
    print(genotype1)
    print("celltype")
    print(celltype1)

 #   print(groups)
 #   print(dim(groups))
 #   print (dge$samples)
    countDataSampleOrder <- colnames(countData)
    print("SampleOrder:")
    print(countDataSampleOrder)
                                        #   groups <- groups[countDataSampleOrder,]
 #   genotype <- factor(c('MLL2KO','MLL2KO','MLL2KO','MLL2KO','MLL2KO','MLL2KO','MLL2KO','MLL2KO','Y2604A','Y2604A','mCXXC','mCXXC','WT','WT','WT','WT','WT','WT','WT','WT'))
#    genotype <- factor(c('Y2604A','Y2604A','mCXXC','mCXXC','MLL2KO','MLL2KO','MLL2KO','MLL2KO','MLL2KO','MLL2KO','MLL2KO','MLL2KO','Y2604A','Y2604A','mCXXC','mCXXC','WT','WT','WT','WT','WT','WT','WT','WT'))
 #   celltype <- factor(c('ES','ES','ES','ES','EB3d','EB6d','EB6d','EB9d','ES','EB6d','ES','EB6d','ES','ES','ES','ES','EB3d','EB3d','EB6d','EB9d'))
    #celltype <- factor(c('ES','ES','ES','ES','ES','ES','ES','ES','EB3d','EB3d','EB6d','EB6d','ES','EB3d','ES','EB3d','ES','ES','ES','ES','EB3d','EB3d','EB6d','EB6d'))
#    dim(genotype)
#    dim(celltype)
#    design2 <- model.matrix(~0+genotype+celltype, data=dge$samples)
    design <- model.matrix(~0+genotype1+celltype1, data=dge$samples)
#    print("Manual design")
#    print(design2)
    print("Table design")
    print(design)
#    dge$samples$group <- groups
#    print (dge$samples)
#    print("Comparisons file:")
#    print(comparisonFile)
#    comparisons <- read.table(file=comparisonFile,sep=",",as.is=T,header=T,row.names=1)
#    print(comparisons)
#    print(dim(comparisons))
#    colnames(comparisons) <- gsub("\\_","\\.",colnames(comparisons))
#    colnames(countData) <- gsub("\\_","\\.",colnames(countData))
#    print ("ComparisonSampleOrder:")
#    print(colnames(comparisons))
#    print("GroupOrder:")
#    print(colnames(groups))
#    design <- model.matrix(0+dge$samples$group)
#    colnames(design) <- levels(dge$samples$group)
#    design
    dge <- calcNormFactors(dge)  #calcs normalization factors based on lib.size (RNA composition effect)
    print("Normalization factors for each sample in comparison, estimated with TMM method.")
    print(dge$samples)
    dge <- estimateGLMCommonDisp(dge,design)  #Estimates a common negative binomial dispersion
    dge <- estimateGLMTrendedDisp(dge,design) #Estimates the abundance-dispersion trend by Cox-Reid approximate profile likelihood
    dge <- estimateGLMTagwiseDisp(dge,design) #Compute an empirical Bayes estimate of the negative binomial dispersion parameter for each transcript, with expression
    fit <- glmFit(dge,design)    #Fit a negative binomial generalized log-linear model to counts
    WTvMLL2KO.ES <- makeContrasts(genotype1MLL2KO - genotype1WT, levels=design)
    lrt <- glmLRT(fit,contrast=WTvMLL2KO.ES)
    foo <- rownames(lrt$table)                #List of differentially expressed genes
#    stopifnot(rownames(dge$counts)==rownames(geneData[foo,]))
    label2<- paste(outputDirectory,paste("WTvMLL2KO.ES",method,"edgeR.txt",sep="."),sep="")
    df <- cbind(lrt$table,adj.p=p.adjust(lrt$table$PValue,method="BH"))
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
    write.table(df,file=label2,sep="\t",col.names=NA)
#    foo <- rownames(lrt$table)                #List of differentially expressed genes   
#    for(c in 1:length(rownames(comparisons))){
#     	comparison <- (paste(c(comparisons[c,1:length(colnames(comparisons))])))
#    	print("Running EdgeR comparison")
#    	print(paste(c(comparisons[c,1:length(colnames(comparisons))])))
#    	print(rownames(comparisons)[c])
#    	comp <-runEdgeR(data=countData,comparison)
#        method <- gsub(".all.counts.txt","",countFile)
#        method <- gsub("^.*\\/","",method)
#   	label1<- paste(outputDirectory,paste(rownames(comparisons)[c],method,"edgeR.df.rda",sep="."),sep="")
#	label2<- paste(outputDirectory,paste(rownames(comparisons)[c],method,"edgeR.txt",sep="."),sep="")
#	save(comp,file=label1)
#   	write.table(comp,file=label2,sep="\t",col.names=NA)
#   }
}

print (runMDS)
if (runMDS==1){
   print("Generating MDS plot.")
   newnames <- counts[,6:length(colnames(counts))]
   newnames <- gsub("\\.\\d+$","",newnames)
   if(filterOff==0){
       dge <- DGEList(countData[rowSums(cpm(countData)>1)>=2,]) #filter low counts & calc lib.sizes
   } else {
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
