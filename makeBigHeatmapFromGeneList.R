args <- commandArgs()

geneListFile <-sub('--geneListFile=', '', args[grep('--geneListFile=', args)])
countFile<-sub('--countFile=','',args[grep('--countFile=',args)])
assembly<-sub('--assembly=','',args[grep('--assembly=',args)])
geneOrder<-sub('--geneOrder=','',args[grep('--geneOrder=',args)])
clusterGenes<-sub('--clusterGenes=','',args[grep('--clusterGenes=',args)])
outputCDT <- sub('--outputCDT=','',args[grep('--outputCDT=',args)])

# This Rscript takes a list of interesting genes, and a matrix of gene
# expression values.  It will pull out the set of interesting genes from
# the matrix and make a heatmap from their expression.

# NB: The items in the gene list need to be row names in the count file.

# Assembly and geneOrder are optional parameters.
# Assembly is used to convert from Ensembl IDs to gene names for better figures.
# geneOrder should be set to a particular gene name in the gene list, if
# desired.  If set, then the samples will be sorted by expression of the named
# gene in the final heatmap.

if (identical(assembly,character(0))){
    assembly<-"na"
    organismStr<- "na"
}
if (identical(geneOrder,character(0))){
    geneOrder<-"na"
}
if (identical(clusterGenes,character(0))){
    clusterGenes <- 1
}
if (identical(outputCDT,character(0))){
    outputCDT <- 0
}

library(gplots)
library(edgeR)
library(biomaRt)

# Set hostMart based on 
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

##----------load differentially expressed genes --------#
print("Loading gene list")
print(geneListFile)

if(grepl('rda',countFile)){
   allCounts <- get(load(file=countFile))
}
if(grepl('txt',countFile)){
    allCounts <- read.table(file=countFile,header=TRUE,sep="\t",row.names=1)
#    head(allCounts)
}
if (countFile == "genomicMatrix"){
    allCounts <- read.table(file=countFile,header=TRUE,sep="\t",row.names=1)
    if(colnames(allCounts)[1]=="probe"){
        rownames(allCounts) = allCounts$probe
    }
    if(colnames(allCounts)[1]=="sample"){
        rownames(allCounts) = allCounts$sample
    }
}
gl <- read.delim(file=geneListFile,header=FALSE,sep="\t")
gl

dim(gl)
dim(allCounts)

if((grepl('normCounts',countFile)) || (grepl('genomicMatrix',countFile)) ||
   (grepl('fpkms',countFile))){
    print("Keep all columns of data")
    sampleNum <- dim(allCounts)[2]
    counts<-allCounts
}else {
    print("Skip first 5 columns of metadata")
    counts<-allCounts[,6:(5+sampleNum)]
    sampleNum <- dim(allCounts)[2]-5
}
#counts<-counts[,order(names(counts))]

dim(counts)
#head(counts)

row.names(gl)<-gl$V1
print ("Gene list")
head(gl)
gl.counts <-counts[row.names(gl),]
print ("Gene list counts")
head(gl.counts)
dim(gl.counts)
print("removing genes with little to no expression")
dim(gl.counts)
gl.counts <-gl.counts[rowSums(abs(gl.counts))>=2,]
print("removing genes with NAs")
gl.counts<-na.omit(gl.counts)
colnames(gl.counts)<- gsub("\\.\\d+$","",colnames(gl.counts))
dim(gl.counts)

if (assembly != "na"){
    print("Replacing Ensembl IDs with common names")
#    hostMart<-"ensembl.org"
    hostMart <- "feb2014.archive.ensembl.org"
    listMarts(host=hostMart)
#    dataset = paste(tolower(organismStr),"_gene_ensembl",sep="")
    bm <- useMart("ENSEMBL_MART_ENSEMBL",host=hostMart)
#    bm <- useMart("ensembl")
    ds <- listDatasets(bm)
    ds <- ds[grep(organismStr,ds$dataset),]$dataset
    print("Loading biomart")
    bm <- useDataset(paste(ds), mart=bm)  
    dataset = paste(tolower(organismStr),"gene_ensembl",sep="_")
    anno <- getBM(mart=bm, attributes=c('ensembl_gene_id','gene_biotype','external_gene_id','description'))
    print(row.names(gl.counts))
    iv <- match(row.names(gl.counts),anno$ensembl_gene_id)
    print("Changing ensembl ids to external gene ids.")
    for (i in 1:length(row.names(gl.counts))){
        if (identical(anno[iv[i],'external_gene_id'],NA_character_)){
#            print(paste("External gene name is missing for ",row.names(gl.counts[i,]),sep=""))
        } else {
            row.names(gl.counts)[i] <- anno[iv[i],'external_gene_id']
        }
    }
}
#print ("Gene list counts")
#head(gl.counts)

if (geneOrder != "na"){
    print(paste("Re-ordering samples according to expression of gene ",geneOrder,sep=""))
    gl.counts <- gl.counts[,order(gl.counts[geneOrder,])]
}
#head(gl.counts)

listName <- gsub("\\.geneList.txt$","",geneListFile)
listName <- gsub("\\.txt$","",listName)
pdfFile <- paste(listName,"exp.heatmap.pdf",sep=".")
if (geneOrder != "na"){
    pdfFile <- paste(listName,geneOrder,"heatmap.pdf",sep=".")
}
print(pdfFile)
print(head(rownames(gl.counts)))
if(grepl('_at$',rownames(gl.counts))){
    rownames(gl.counts) <- gl[,2]
    rownames(gl.counts) <- gsub(",",".",rownames(gl.counts))
}

listLabel <- gsub("^.*/","",listName)
pdf(pdfFile)
#heatmap.2(as.matrix(gl.counts),trace="none",dendrogram="col",scale="row",Rowv=FALSE,Colv=TRUE,labRow=row.names(gl.counts),cexRow=1,colsep=0,rowsep=0,cexCol=1,margins=c(12,9),main=listLabel,col=colorRampPalette(c("blue","white","red"))(75))
if (clusterGenes == 1){
    gl.cor<-(cor(t(gl.counts),use="pairwise.complete.obs",method="pearson"))
    gl.cor.dist<-as.dist(1-gl.cor)
    gl.tree<-hclust(gl.cor.dist,method='average')
    heatmap.2(as.matrix(gl.counts),trace="none",dendrogram="row",scale="row",Rowv=as.dendrogram(gl.tree),Colv=NA,labRow=row.names(gl.counts),cexRow=1,colsep=0,rowsep=0,cexCol=1,margins=c(12,9),main=listLabel,col=colorRampPalette(c("blue","white","red"))(75))
} else {
    heatmap.2(as.matrix(gl.counts),trace="none",dendrogram="none",scale="row",Rowv=NA,Colv=NA,labRow=row.names(gl.counts),cexRow=1,colsep=0,rowsep=0,cexCol=1,margins=c(12,9),main=listLabel,col=colorRampPalette(c("blue","white","red"))(75))
}
dev.off()
if (clusterGenes == 1){
    gl.sorted <-gl.counts[rev(gl.tree$order),]
} else {
    gl.sorted <- gl.counts
}
write.table(gl.sorted,paste(listName,"clustered.txt",sep="."),sep="\t")
if(outputCDT == 1){
#    gl.sorted <- log(gl.sorted[,1:length(colnames(gl.sorted))])
    gl.sorted <- cbind(row.names(gl.sorted),gl.sorted)
    write.table(gl.sorted,paste(listName,"expression.cdt",sep="."),sep="\t")
}


