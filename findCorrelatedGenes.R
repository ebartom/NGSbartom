args <- commandArgs()

countFile <-sub('--countFile=', '', args[grep('--countFile=', args)])
geneName <- sub('--gene=', '', args[grep('--gene=', args)])
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
if(grepl('genomicMatrix',countFile)){
   counts <- read.delim(file=countFile,header=TRUE,sep="\t")
}

if (countFile != 'genomicMatrix'){
    sampleNum <- dim(counts)[2] - 5
    sampleNum
    samples<-colnames(counts)[6:(5+sampleNum)]
    samples
    counts<-counts[,6:(5+sampleNum)]
    dim(counts)
}
if (countFile == 'genomicMatrix'){
    print ("Input file is genomicMatrix")
    dim(counts)
#    head(counts)[1:5]
    rownames(counts)=counts$sample
#    print(rownames(counts)[1:5])
    sampleNum <- dim(counts)[2]-1
    counts <- counts[,2:sampleNum+1]
#    dim(counts)
}
#head(counts)[1:5]
counts <- counts[rowSums(counts>1)>=3,]
dim(counts)
geneName
#print(rownames(counts)[1:5])
geneCounts <- as.numeric(counts[(rownames(counts))==geneName,])
geneCounts

correlations = data.frame(gene1=rep(0,dim(counts)[1]),gene2=rep(0,dim(counts)[1]),corrValue=rep(0,dim(counts)[1]))
for (i in 1:dim(counts)[1]){
#for (i in 1:dim(counts)[1]){
#    print(i)
#    print(rownames(counts)[i])
    newCounts <- as.numeric(counts[i,])
#    print(geneCounts[1:5])
#    print(newCounts[1:5])
    correlation <- cor(geneCounts,newCounts,method="pearson")
    correlations[i,] <- c(geneName,rownames(counts)[i],correlation)
}
#head(correlations)
                                        #print(correlations$corrValue)
correlations <- na.omit(correlations)
correlations$corrValue <- as.numeric(correlations$corrValue)
correlations <- correlations[rev(order(correlations$corrValue)),]
write.table(correlations,file=paste(geneName,"corrTable.txt",sep="."),sep="\t")
