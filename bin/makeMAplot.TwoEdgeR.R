
args <- commandArgs()

degFile1 <-sub('--degFile1=', '', args[grep('--degFile1=', args)])
degFile2 <-sub('--degFile2=', '', args[grep('--degFile2=', args)])
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
print("Loading differential expressed gene tables")

print(degFile1)

if(grepl('rda',degFile1)){
    deg1 <- get(load(file=degFile1))
    deg1 <- deg1[order(deg1$adj.p),]
}
if(grepl('txt',degFile1)){
    deg1 <- read.delim(file=degFile1,header=TRUE,sep="\t")
    rownames(deg1)<-deg1$X
    deg1 <- deg1[order(deg1$adj.p),]
}
print(degFile1)
head(deg1)
dim(deg1)

print(degFile2)

if(grepl('rda',degFile2)){
    deg2 <- get(load(file=degFile2))
    deg2 <- deg2[order(deg2$adj.p),]
}
if(grepl('txt',degFile2)){
    deg2 <- read.delim(file=degFile2,header=TRUE,sep="\t")
    rownames(deg2)<-deg2$X
    deg2 <- deg2[order(deg2$adj.p),]
}
print(degFile2)
head(deg2)
dim(deg2)

sampleNum1 <- dim(deg1)[2] - 16
sampleNum2 <- dim(deg2)[2] - 16
sampleNum1
sampleNum2
deg1indeg2 <- na.omit(deg1[row.names(deg2),])
deg2indeg1 <- na.omit(deg2[row.names(deg1),])
deg1indeg2 <- deg1indeg2[order(row.names(deg1indeg2)),]
deg2indeg1 <- deg2indeg1[order(row.names(deg2indeg1)),]

print("deg1indeg2")
print(head(deg1indeg2))
dim(deg1indeg2)

print("deg2indeg1")
print(head(deg2indeg1))
dim(deg2indeg1)


                                        #gl.counts <-counts[row.names(gl),]


sigDeg1 <- deg1indeg2$adj.p < adjp
print("sigDeg1")
head(sigDeg1)
print(sum(sigDeg1))
genes1 <- na.omit(row.names(deg1indeg2[sigDeg1,]))
head(genes1)
#print(deg1indeg2[genes1,"gene"])
#print(deg1indeg2[genes1,"adj.p"])

sigDeg2 <- deg2indeg1$adj.p < adjp
print("sigDeg2")
head(sigDeg2)
print(sum(sigDeg2))
genes2 <- na.omit(row.names(deg2indeg1[sigDeg2,]))
head(genes2)
#print(deg2indeg1[genes2,"gene"])
#print(deg2indeg1[genes2,"adj.p"])

if (sum(sigDeg1) == 0){
    print("No significant genes in set1!")
}

if (sum(sigDeg2) == 0){
    print("No significant genes in set2!")
}

joint <- na.omit(intersect(genes1,genes2))
print("joint")
print(deg1indeg2[joint,"gene"])
print(length(joint))
print(adjp)

if (labelTop == 1){
    geneLabels<- as.graphicsAnnot(deg1[joint,"gene"])
}
#if (labelList >0){
#    myGenes <- deg1indeg2$X %in% labels[,1]
#    myGenes <- deg1indeg2[myGenes,"X"]
#    print(length(myGenes))
#    print(sum(myGenes))
    ## firstGene <- deg1[myGenes,"X"][1]
    ## print(firstGene)
    ## if (substring(firstGene,1,3) == "ENS"){
    ##     print("This is an ensembl ID")
    ##     geneLabels <- as.graphicsAnnot(deg1[myGenes,"gene"])
    ## } else{
    ##     geneLabels <- as.graphicsAnnot(deg1[myGenes,"X"])
    ##    }
#}

adjplabel <- gsub("^0\\.","",adjp)
comparison1 <- basename(gsub("\\.edgeR.txt$","",degFile1))
comparison2 <- basename(gsub("\\.edgeR.txt$","",degFile2))
#pngFile <- paste(comparison,adjplabel,"MAplot.png",sep=".")
#print(pngFile)
pdfFile <- paste(comparison1,comparison2,adjplabel,"MAplot.pdf",sep=".")
print(pdfFile)
comparison1 <- gsub("^.*/","",comparison1)
comparison2 <- gsub("^.*/","",comparison2)

#png(pngFile,height=500,width=500)
pdf(pdfFile)
par(mar=c(5.1,4.6,4.1,2.1))
plot(deg1indeg2$logFC,deg2indeg1$logFC
      ,col="grey87"
      ,pch=19
      ,cex=1
 #     ,ylim=c(-4,2) 
      ,main="Relative log2 Fold Change",xlab=paste("log2FC",comparison1),ylab=paste("log2FC",comparison2)
      ,cex.lab=1.6, cex.axis=1.6, cex.main=1.6, cex.sub=1.6)
if (length(joint) > 0){
    points(deg1indeg2[joint,"logFC"],deg2indeg1[joint,"logFC"], col="red", pch=19,cex=1)
}
legend("bottomright", c(paste(length(joint),"genes both sets p.adj <",adjp),paste(dim(deg1indeg2)[1],"genes in both datasets")),pch=c(19,19), col=c("red","grey"),cex=0.9)
abline(h=0, lty="dashed", col="grey")
abline(v=0, lty="dashed", col="grey")
print(labelTop)
print(labelList)
if (labelTop == 1){
    if (length(joint) > 0){
        text(deg1indeg2[joint,"logFC"],deg2indeg1[joint,"logFC"],labels=deg1indeg2[joint,"gene"],pos=3)
    }
}
if (labelList > 0){
    print("findMatches")
    myGenes <- deg1indeg2$X %in% labels[,1]
    print("pull out matching")
    myGenes <- deg1indeg2[myGenes,"X"]
    print(myGenes)
#    print(length(myGenes))
                                        #    print(sum(myGenes))
    print("First, one way")
    print(myGenes[1])
    if (substring(myGenes[1],1,3) == "ENS"){
        print("This is an ensembl ID")
        deg1indeg2rows <- deg1indeg2$X %in% myGenes
        deg2indeg1rows <- deg2indeg1$X %in% myGenes
        print(deg1indeg2[deg1indeg2rows,])
        print(deg2indeg1[deg2indeg1rows,])
        text(deg1indeg2[deg1indeg2rows,"logFC"],deg2indeg1[deg2indeg1rows,"logFC"],labels=deg1indeg2[deg1indeg2rows,"gene"],pos=3)
        ## print ("deg1indeg2")
        ## print(deg1indeg2[firstGene,])
        ## print ("deg1")
        ## print(deg1[firstGene,])
        ## print ("deg2indeg1")
        ## print(deg2indeg1[firstGene,])
        ## print ("deg2")
        ## print(deg2[firstGene,])
    } else{
        deg1indeg2rows <- deg1indeg2$gene %in% myGenes
        deg2indeg1rows <- deg2indeg1$gene %in% myGenes
        text(deg1indeg2[deg1indeg2rows,"logFC"],deg2indeg1[deg2indeg1rows,"logFC"],labels=deg1indeg2[deg1indeg2rows,"X"],pos=3)
    }
    points(deg1indeg2[deg1indeg2rows,"logFC"],deg2indeg1[deg2indeg1rows,"logFC"], col="black", pch=1,cex=1)
}
dev.off()
