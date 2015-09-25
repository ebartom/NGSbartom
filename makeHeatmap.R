args <- commandArgs()

degFile <-sub('--degFile=', '', args[grep('--degFile=', args)])
adjp <-sub('--adjp=', '', args[grep('--adjp=', args)])

if (identical(adjp,character(0))){
   adjp<-0.01
}

library(gplots)
library(edgeR)

##----------load differentially expressed genes --------#
print("Loading differential expressed gene table")
print(degFile)

if(grepl('rda',degFile)){
   deg <- get(load(file=degFile))
}
if(grepl('txt',degFile)){
   deg <- read.delim(file=degFile,header=TRUE,sep="\t")
}
#head(deg)
dim(deg)
sampleNum <- dim(deg)[2] - 16
print(sampleNum)
counts<-deg[,7:(6+sampleNum)]
rownames(counts)<-deg$X
head(counts)


up <- deg$adj.p < adjp & deg$logFC > 0
sum(up)

down <- deg$adj.p < adjp & deg$logFC < 0
sum(down)

up.counts <- counts[up,]
dim(up.counts)
down.counts<-counts[down,]
dim(down.counts)

print("Check uniqueness.")
print(length(unique(deg[up,12+sampleNum])))
print(length(unique(rownames(deg[up,]))))
rownames(up.counts)<-deg[up,12+sampleNum]
head(up.counts)
rownames(down.counts)<-deg[down,12+sampleNum]
head(down.counts)

adjplabel <- gsub("^0\\.","",adjp)
comparison <- gsub("\\.edgeR.txt$","",degFile)
UpPdfFile <- paste(comparison,adjplabel,"up.heatmap.pdf",sep=".")
DnPdfFile <- paste(comparison,adjplabel,"dn.heatmap.pdf",sep=".")
print(UpPdfFile)
print(DnPdfFile)

comparisonLabel <- gsub("^.*/","",comparison)
up.cor<-(cor(t(up.counts[,1:sampleNum]),use="pairwise.complete.obs",method="pearson"))
up.cor.dist<-as.dist(1-up.cor)
up.tree<-hclust(up.cor.dist,method='average')
pdf(UpPdfFile)
heatmap.2(as.matrix(up.counts[,1:sampleNum]),trace="none",dendrogram="row",scale="row",Rowv=as.dendrogram(up.tree),Colv=NA,labRow=rownames(up.counts),cexRow=1,colsep=0,rowsep=0,cexCol=0.7,main=paste("Up (<",adjp,") ",comparisonLabel,sep=""))
dev.off()
up.sorted <-up.counts[rev(up.tree$order),]
write.table(up.sorted,paste(comparison,adjplabel,"up.clustered.txt",sep="."),sep="\t")

down.cor<-(cor(t(down.counts[,1:sampleNum]),use="pairwise.complete.obs",method="pearson"))
down.cor.dist<-as.dist(1-down.cor)
down.tree<-hclust(down.cor.dist,method='average')
pdf(DnPdfFile)
heatmap.2(as.matrix(down.counts[,1:sampleNum]),trace="none",dendrogram="row",scale="row",Rowv=as.dendrogram(down.tree),Colv=NA,labRow=rownames(down.counts),cexRow=1,cexCol=0.7,colsep=0,rowsep=0,main=paste("Down (<",adjp,") ",comparisonLabel,sep=""))
dev.off()
down.sorted <-down.counts[rev(down.tree$order),]
write.table(down.sorted,paste(comparison,adjplabel,"dn.clustered.txt",sep="."),sep="\t")
