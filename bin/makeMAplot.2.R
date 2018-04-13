
args <- commandArgs()

degFile <-sub('--degFile=', '', args[grep('--degFile=', args)])
adjp <-as.numeric(sub('--adjp=', '', args[grep('--adjp=', args)]))
labelTop <-sub('--labelTop=', '', args[grep('--labelTop=', args)])
labelList <-sub('--labelList=', '', args[grep('--labelList=', args)])

if (identical(adjp,character(0))){
   adjp<-0.01
}
if (identical(labelTop,character(0))){
   labelTop<-0
}
if (identical(labelList,character(0))){
   labelList<-""
}

if(grepl('txt',labelList)){
    labels <- read.delim(file=labelList,header=FALSE,sep="\t")
    print("Loading list of genes to label")
    print(labels)
    dim(labels)
}

library(edgeR)

##----------load differentially expressed genes --------#
print("Loading differential expressed gene table")
print(degFile)

if(grepl('rda',degFile)){
    deg <- get(load(file=degFile))
        deg <- deg[order(deg$adj.p),]
}
if(grepl('txt',degFile)){
    deg <- read.delim(file=degFile,header=TRUE,sep="\t")
    rownames(deg)<-deg$X
    deg <- deg[order(deg$adj.p),]
}
head(deg)
dim(deg)
sampleNum <- dim(deg)[2] - 16

up <- deg$adj.p < adjp & deg$logFC > 0
sum(up)

down <- deg$adj.p < adjp & deg$logFC < 0
sum(down)

if (labelTop == 1){
    topOnePercent <- deg$adj.p < adjp/100
    print(sum(topOnePercent))
    geneLabels<- as.graphicsAnnot(deg[topOnePercent,"gene"])
}
if (labelList >0){
    myGenes <- deg$X %in% labels[,1]
    print(sum(myGenes))
    geneLabels <- as.graphicsAnnot(deg[myGenes,"X"])
}

adjplabel <- gsub("^0\\.","",adjp)
comparison <- gsub("\\.edgeR.txt$","",degFile)
#pngFile <- paste(comparison,adjplabel,"MAplot.png",sep=".")
#print(pngFile)
pdfFile <- paste(comparison,adjplabel,"MAplot.pdf",sep=".")
print(pdfFile)
comparison <- gsub("^.*/","",comparison)

#png(pngFile,height=500,width=500)
pdf(pdfFile)
par(mar=c(5.1,4.6,4.1,2.1))
plot(deg$logCPM,deg$logFC
      ,col="grey87"
      ,pch=19
      ,cex=1
 #     ,ylim=c(-4,2) 
      ,main=comparison,xlab="log2(CPM)",ylab="log2FC"
      ,cex.lab=1.6, cex.axis=1.6, cex.main=1.6, cex.sub=1.6)
  points(deg[up,"logCPM"],deg[up,"logFC"], col="purple", pch=19,cex=1)
  points(deg[down,"logCPM"],deg[down,"logFC"], col="green", pch=19,cex=1)
legend("bottomright", c(paste(sum(up),"genes p.adj-val < ",adjp),paste(sum(down),"genes p.adj-val < ",adjp)),
  pch=c(19,19,1), col=c("purple","green","black"),cex=0.9,)
abline(h=0, lty="dashed", col="grey")

print(labelTop)
print(labelList)
if (labelTop == 1){
    if (sum(topOnePercent) > 0){
        text(deg[topOnePercent,"logCPM"],deg[topOnePercent,"logFC"],labels=deg[topOnePercent,"gene"],pos=3)
    }
}
if (labelList > 0){
    text(deg[myGenes,"logCPM"],deg[myGenes,"logFC"],labels=deg[myGenes,"X"],pos=3)
    points(deg[myGenes,"logCPM"],deg[myGenes,"logFC"], col="black", pch=1,cex=1)
}
dev.off()
