args <- commandArgs()

geneListFile <-sub('--geneListFile=', '', args[grep('--geneListFile=', args)])
countFile<-sub('--countFile=','',args[grep('--countFile=',args)])
assembly<-sub('--assembly=','',args[grep('--assembly=',args)])

if (identical(assembly,character(0))){
    assembly<-"na"
    organismStr<- "na"
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
    head(allCounts)
}
gl <- read.delim(file=geneListFile,header=FALSE,sep="\t")

dim(gl)
dim(allCounts)

if(grepl('normCounts',countFile)){
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
head(gl)
gl.counts <-counts[row.names(gl),]
head(gl.counts)
dim(gl.counts)
gl.counts <-gl.counts[rowSums(gl.counts)>=2,]
gl.counts<-na.omit(gl.counts)
colnames(gl.counts)<- gsub("\\.\\d+$","",colnames(gl.counts))
dim(gl.counts)
head(gl.counts)

if (assembly != "na"){
    hostMart<-"ensembl.org"
    listMarts(host=hostMart)
    bm <- useMart("ensembl")
    ds <- listDatasets(bm)
    ds <- ds[grep(organismStr,ds$dataset),]$dataset
    print("Loading biomart")
    bm <- useDataset(paste(ds), mart=bm)  
    dataset = paste(tolower(organismStr),"gene_ensembl",sep="_")
    anno <- getBM(mart=bm, attributes=c('ensembl_gene_id','gene_biotype','external_gene_name','description'))
    print(row.names(gl.counts))
    iv <- match(row.names(gl.counts),anno$ensembl_gene_id)
    print("Changing ensembl ids to external gene names.")
    for (i in 1:length(row.names(gl.counts))){
        if (identical(anno[iv[i],'external_gene_name'],NA_character_)){
#            print(paste("External gene name is missing for ",row.names(gl.counts[i,]),sep=""))
        } else {
            row.names(gl.counts)[i] <- anno[iv[i],'external_gene_name']
        }
    }
}
head(gl.counts)

listName <- gsub("\\.geneList.txt$","",geneListFile)
pdfFile <- paste(listName,"heatmap.pdf",sep=".")
print(pdfFile)

listLabel <- gsub("^.*/","",listName)
gl.cor<-(cor(t(gl.counts),use="pairwise.complete.obs",method="pearson"))
gl.cor.dist<-as.dist(1-gl.cor)
gl.tree<-hclust(gl.cor.dist,method='average')
pdf(pdfFile)
heatmap.2(as.matrix(gl.counts),trace="none",dendrogram="row",scale="row",Rowv=as.dendrogram(gl.tree),Colv=NA,labRow=row.names(gl.counts),cexRow=1,colsep=0,rowsep=0,cexCol=1,margins=c(12,9),main=listLabel,col=colorRampPalette(c("blue","white","red"))(75))
#heatmap.2(as.matrix(gl.counts),trace="none",dendrogram="row",scale="row",Rowv=as.dendrogram(gl.tree),Colv=NA,labRow=row.names(gl.counts),cexRow=1,colsep=0,rowsep=0,cexCol=1,margins=c(12,9),main=listLabel,col=redblue)
dev.off()
gl.sorted <-gl.counts[rev(gl.tree$order),]
write.table(gl.sorted,paste(listName,"clustered.txt",sep="."),sep="\t")

