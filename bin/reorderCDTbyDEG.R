library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(Rsamtools)
library(rtracklayer)
library(biomaRt)
#library(ChIPpeakAnno)


args <- commandArgs()

cdtFile <-sub('--cdtFile=', '', args[grep('--cdtFile=', args)])
degFile <- sub('--degFile=','',args[grep('--degFile=',args)])
order <- sub('--order=','',args[grep('--order=',args)])

deg <- as.data.frame(read.delim(file=degFile,header=TRUE,sep="\t"))
rownames(deg) <- deg$X
cdt <- as.data.frame(read.delim(file=cdtFile,header=TRUE,sep="\t"))
rownames(cdt) <- cdt$UID

littleDeg <- deg[rownames(cdt),]
if ((order == "adjp")|| (order == "adj.p")){
    littleDeg <- littleDeg[order(littleDeg$adj.p),]
}
if (order == "logFC"){
    littleDeg <- littleDeg[rev(order(littleDeg$logFC)),]
}
cdt <- cdt[order(rownames(littleDeg)),]
#head(littleDeg)

cdtFile <- gsub("cdt$",paste(order,"cdt",sep="."),cdtFile)
write.table(cdt, file=cdtFile, sep="\t",row.names=FALSE)
