library(gplots)
library(ggplot2)
#library(ggfortify)

args<-commandArgs()
fpkmTable <- sub('--fpkmTable=','',args[grep('--fpkmTable=',args)])

fpkms <- read.table(fpkmTable,header=TRUE,row.names=1)

numSamples<-dim(fpkms)[2]
numSamples
samples<-colnames(fpkms)

emptyCorrelations <- rep(0,numSamples*numSamples)
corrTable<- matrix(emptyCorrelations,nrow=numSamples,ncol=numSamples)
i<-0
for (sample1 in samples){
    i<-i+1
    i
    j<-0
    for (sample2 in samples){
      	j<- j+1
	j
	correlation <- cor(fpkms[,i],fpkms[,j])
	print(correlation)
	corrTable[i,j]<-correlation
	}
}
rownames(corrTable)<-samples
colnames(corrTable)<-samples

label<-gsub("\\.txt$","",fpkmTable)
pdfFile<-paste(label,"fpkmCorr.pdf",sep=".")
tableFile<-paste(label,"fpkmCorr.txt",sep=".")
pdfFile
tableFile

head(corrTable)
write.table(corrTable,file=tableFile)
pdf(pdfFile)
heatmap.2(corrTable,trace="none")
dev.off()

pcaFile<-paste(label,"pcaPlot.pdf",sep=".")
pdf(pcaFile)
pca <-prcomp(fpkms)
summary(pca)
#plot(pca)plot(pca$x[,1:kd
pl
dev.off()