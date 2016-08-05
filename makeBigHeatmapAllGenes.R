args <- commandArgs()

countFile<-sub('--countFile=','',args[grep('--countFile=',args)])

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
    counts<-allCounts[,6:(5+sampleNum)]
    sampleNum <- dim(allCounts)[2]-5
} else {
    print("Keep all columns of data")
    sampleNum <- dim(allCounts)[2]
    counts<-allCounts
}    
#counts<-counts[,order(names(counts))]

dim(counts)
head(counts)

#nonZero <- rowSums(counts) >= 1
#lowCount <- rowSums(counts) < 2
#sum(nonZero)
#head(counts[lowCount,])
#counts<-counts[nonZero,]
#dim(counts)

countFile <- gsub("\\.txt$","",countFile)
pdfFile <- paste(countFile,"heatmap.pdf",sep=".")
print(pdfFile)

combined.cor<-(cor(t(counts[,1:sampleNum]),use="pairwise.complete.obs",method="pearson"))
combined.cor.dist<-as.dist(1-combined.cor)
combined.tree<-hclust(combined.cor.dist,method='average')
pdf(pdfFile)
heatmap.2(as.matrix(counts[,1:sampleNum]),trace="none",dendrogram="row",scale="row",Rowv=as.dendrogram(combined.tree),Colv=NA,labRow=rownames(counts),cexRow=1,colsep=0,rowsep=0,cexCol=1,margins=c(12,9),main=countFile,col=colorRampPalette(c("blue","white","red"))(75))
dev.off()
combined.sorted <-counts[rev(combined.tree$order),]
write.table(combined.sorted,paste(countFile,"clustered.txt",sep="."),sep="\t")

