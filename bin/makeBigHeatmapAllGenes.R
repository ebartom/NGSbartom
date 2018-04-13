args <- commandArgs()

countFile<-sub('--countFile=','',args[grep('--countFile=',args)])
filterLowCounts <- sub('--filterLowCounts=','',args[grep('--filterLowCounts=',args)])
clusterSamples<-sub('--clusterSamples=','',args[grep('--clusterSamples=',args)])
outputCDT <- sub('--outputCDT=','',args[grep('--outputCDT=',args)])

if (identical(filterLowCounts,character(0))){
    filterLowCounts <- 0
}
if (identical(clusterSamples,character(0))){
    clusterSamples <- 0
}
if (identical(outputCDT,character(0))){
    outputCDT <- 0
}

library(gplots)
library(edgeR)

print("Loading counts table")
print(countFile)

if(grepl('rda',countFile)){
   allCounts <- get(load(file=countFile))
}
if(grepl('txt',countFile)){
    allCounts <- read.delim(file=countFile,header=TRUE,sep="\t")
#    head(allCounts)
    if (colnames(allCounts)[1] == "X"){
        rownames(allCounts) <- allCounts$X
        allCounts <- allCounts[,2:dim(allCounts)[2]]
    }
    print(colnames(allCounts))
    head(allCounts)
}

if (grepl('counts.txt$',countFile)){
    print("Skip first 5 columns of metadata")
    sampleNum <- dim(allCounts)[2]-5
    counts<-allCounts[,6:(5+sampleNum)]
} else if (grepl('fractions.txt$',countFile)){
    print("Skip first two columns of metadata")
    sampleNum <- dim(allCounts)[2]-2
    print(paste("SampleNum is",sampleNum))
    counts<-allCounts[,3:(2+sampleNum)]
    rownames(counts)<-allCounts[,2]
} else if (grepl('coverages.txt$',countFile)){
    print("Skip first two columns of metadata")
    sampleNum <- dim(allCounts)[2]-2
    print(paste("SampleNum is",sampleNum))	
    counts<-allCounts[,3:(2+sampleNum)]
    rownames(counts)<-allCounts[,2]
} else {
    print("Keep all columns of data")
    sampleNum <- dim(allCounts)[2]
    counts<-allCounts
}    
#counts<-counts[,order(names(counts))]

dim(counts)
head(counts)

if (filterLowCounts == 1) {
    print("Filtering out genes with low counts")
    nonZero <- rowSums(counts) >= 1
    notLowCount <- rowSums(counts) > 2
    print(head(notLowCount))
    print(sum(notLowCount))
#    head(counts[lowCount,])
    counts<-counts[notLowCount,]
    dim(counts)
}
dim(counts)

countFile <- gsub("\\.txt$","",countFile)
pdfFile <- paste(countFile,"heatmap.pdf",sep=".")
print(pdfFile)

#combined.cor<-(cor(t(counts[,1:sampleNum]),use="pairwise.complete.obs",method="pearson"))
combined.cor<-(cor(t(counts),use="pairwise.complete.obs",method="pearson"))
combined.cor.dist<-as.dist(1-combined.cor)
combined.tree<-hclust(combined.cor.dist,method='average')
pdf(pdfFile)
if (clusterSamples==1){
    print("Clustering genes and samples")
     cor2<-(cor(counts,use="pairwise.complete.obs",method="pearson"))
    cor2.dist<-as.dist(1-cor2)
    tree2<-hclust(cor2.dist,method='average')
    heatmap.2(as.matrix(counts[1:sampleNum]),trace="none",dendrogram="both",scale="row",Rowv=as.dendrogram(combined.tree),Colv=as.dendrogram(tree2),labRow=rownames(counts),cexRow=1,colsep=0,rowsep=0,cexCol=1,margins=c(12,9),main=countFile,col=colorRampPalette(c("blue","white","red"))(75))
} else {
    heatmap.2(as.matrix(counts[,1:sampleNum]),trace="none",dendrogram="row",scale="row",Rowv=as.dendrogram(combined.tree),Colv=NA,labRow=rownames(counts),cexRow=1,colsep=0,rowsep=0,cexCol=1,margins=c(12,9),main=countFile,col=colorRampPalette(c("blue","white","red"))(75))
}
dev.off()
combined.sorted <-counts[rev(combined.tree$order),]
write.table(combined.sorted,paste(countFile,"clustered.txt",sep="."),sep="\t")

if(outputCDT == 1){
    gl.sorted <- cbind(row.names(combined.sorted),combined.sorted)
    write.table(combined.sorted,paste(countFile,"cdt",sep="."),sep="\t")
}
