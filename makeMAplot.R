
args <- commandArgs()

degFile <-sub('--degFile=', '', args[grep('--degFile=', args)])
adjp <-as.numeric(sub('--adjp=', '', args[grep('--adjp=', args)]))
labelTop <-sub('--labelTop=', '', args[grep('--labelTop=', args)])

if (identical(adjp,character(0))){
   adjp<-0.01
}
if (identical(adjp,character(0))){
   labelTop<-1
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

topOnePercent <- deg$adj.p < adjp/100

geneLabels<- as.graphicsAnnot(deg[topOnePercent,"gene"])

adjplabel <- gsub("^0\\.","",adjp)
comparison <- gsub("\\.edgeR.txt$","",degFile)
pngFile <- paste(comparison,adjplabel,"MAplot.png",sep=".")
print(pngFile)
comparison <- gsub("^.*/","",comparison)

png(pngFile,height=500,width=500)
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

print(sum(topOnePercent))
print(labelTop)
if ((sum(topOnePercent) > 0) && (labelTop == 1)){
    text(deg[topOnePercent,"logCPM"],deg[topOnePercent,"logFC"],labels=deg[topOnePercent,"gene"],pos=3)
}
  dev.off()
