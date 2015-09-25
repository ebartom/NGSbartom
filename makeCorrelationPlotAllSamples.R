args <- commandArgs()

countFile <-sub('--countFile=', '', args[grep('--countFile=', args)])
outputDirectory <-sub('--outputDirectory=', '',args[grep('--outputDirectory=',args)])
library(edgeR)

##----------load counts table-----------##
print("Loading counts table")
print(countFile)

if(grepl('rda',countFile)){
   counts <- get(load(file=countFile))
}
if(grepl('txt',countFile)){
   counts <- read.delim(file=countFile,header=TRUE,sep="\t")
}

sampleNum <- dim(counts)[2] - 5
sampleNum
samples<-colnames(counts)[6:(5+sampleNum)]
samples
counts<-counts[,6:(5+sampleNum)]
dim(counts)
counts <- counts[rowSums(cpm(counts)>1)>=3,]
head(counts)
dim(counts)
#counts.corr<-cor(counts,method="spearman")
#counts.corr

#sampleNum <-3
png(paste(outputDirectory,"allSamples.correlationPlot.png",sep="/"),height=(sampleNum*200),width=(sampleNum*200))
    par(mfrow=c(sampleNum,sampleNum))
    for (i in 1:sampleNum) {
        for (j in 1:sampleNum){
            plot(log10(counts[,i]), log10(counts[,j]),
                 xlab=paste("log10(",paste(samples[i]),")",sep=""),
                 ylab=paste("log10(",paste(samples[j]),")",sep=""),
                 #cex=1.5, cex.lab=1.5,cex.axis=1.8, cex.main=1.8,
                 col="#1C0DFF15",pch=20)
            legend("topleft",legend= round(cor(counts[,i],counts[,j],method="pearson"), 4), cex=1.4)
        }
    }
    dev.off()
