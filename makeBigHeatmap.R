args <- commandArgs()

degFile <-sub('--degFile=', '', args[grep('--degFile=', args)])
adjp <-as.numeric(sub('--adjp=', '', args[grep('--adjp=', args)]))
countFile<-sub('--countFile=','',args[grep('--countFile=',args)])

if (identical(adjp,character(0))){
   adjp<-0.01
}

library(gplots)
library(edgeR)

##----------load differentially expressed genes --------#
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
print("Loading differential expressed gene table")
print(degFile)
if(grepl('rda',degFile)){
   deg <- get(load(file=degFile))
}
if(grepl('txt',degFile)){
    deg <- read.delim(file=degFile,header=TRUE,sep="\t")
    rownames(deg)<-deg$X
}

print("DEGfile dimensions")
dim(deg)
print("countfile dimensions")
dim(allCounts)

print("number of genes below p-value")
combined <- deg$adj.p < adjp
print(sum(combined))
head(combined)

combined.genes <- na.omit(deg[combined,"gene"])
                                        #combined.genes
head(combined.genes)
combined.ids <- rownames(deg[combined,])
print(paste("IDs of genes with p-values < ",adjp,sep=""))
head(combined.ids)
#length(combined.ids)

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
head(counts)

combined.counts <-counts[combined.ids,]
#colnames(combined.counts)<- gsub("\\.\\d+$","",colnames(combined.counts))
dim(combined.counts)
                                        #head(combined.counts)
if (length(unique(combined.genes)) == length(unique(combined.ids))){
    rownames(combined.counts)<-combined.genes
}

#nonZero <- rowSums(combined.counts) >= 1
#lowCount <- rowSums(combined.counts) < 2
#sum(nonZero)
#head(combined.counts[lowCount,])
#combined.counts<-combined.counts[nonZero,]
#dim(combined.counts)

adjplabel <- gsub("^0\\.","",adjp)
comparison <- gsub("\\.edgeR.txt$","",degFile)
pdfFile <- paste(comparison,adjplabel,"heatmap.pdf",sep=".")
print(pdfFile)

comparisonLabel <- gsub("^.*/","",comparison)
combined.cor<-(cor(t(combined.counts[,1:sampleNum]),use="pairwise.complete.obs",method="pearson"))
combined.cor.dist<-as.dist(1-combined.cor)
combined.tree<-hclust(combined.cor.dist,method='average')
pdf(pdfFile)
heatmap.2(as.matrix(combined.counts[,1:sampleNum]),trace="none",dendrogram="row",scale="row",Rowv=as.dendrogram(combined.tree),Colv=NA,labRow=rownames(combined.counts),cexRow=1,colsep=0,rowsep=0,cexCol=1,margins=c(12,9),main=paste("DEG (<",adjp,") ",comparisonLabel,sep=""),col=colorRampPalette(c("blue","white","red"))(75))
dev.off()
combined.sorted <-combined.counts[rev(combined.tree$order),]
write.table(combined.sorted,paste(comparison,adjplabel,"combined.clustered.txt",sep="."),sep="\t")

