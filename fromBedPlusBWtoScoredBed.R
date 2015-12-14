library(GenomicRanges)
library(rtracklayer)

args <- commandArgs()

bwfile <- sub('--bwfile=','',args[grep('--bwfile=',args)])
bedfile <- sub('--bedfile=','',args[grep('--bedfile=',args)])
scoreMultiplier <- sub('--factor=','',args[grep('--factor=',args)])

if (identical(scoreMultiplier,character(0))){
    scoreMultiplier <- 100
}

bwfile
bedfile
scoreMultiplier <- as.numeric(scoreMultiplier)
scoreMultiplier
peaks <- read.delim(file=bedfile,header=FALSE,sep="\t")
head(peaks)
dim(peaks)
if (dim(peaks)[2] == 3){
    colnames(peaks)<-c("chr","start","end")
    peaks$name <- paste(basename(bedfile),rownames(peaks),sep="_")
} else if (dim(peaks)[2] == 4){
    colnames(peaks)<-c("chr","start","end","score")
    peaks$name <- paste(basename(bedfile),rownames(peaks),sep="_")
} else if (dim(peaks)[2] == 5){
    colnames(peaks)<-c("chr","start","end","name","score")
} else if (dim(peaks)[2] == 6){
    colnames(peaks)<-c("chr","start","end","name","score","strand")
}
rownames(peaks) <- peaks$name
print("Converting peaks to granges object")
gpeaks <- as(peaks,"GRanges")
head(gpeaks)

cat("importing:", bwfile, sep="\n")
bw.peak <- import.bw(bwfile,RangedData=FALSE,selection = BigWigSelection(gpeaks))
cat("calc coverage\n")
bw.peak.cov <- coverage(bw.peak,weight='score')

cat("get summit for peak region\n")
max.cov <- with(as.data.frame(gpeaks),{       
    mcmapply(function(seqname,start,end){
        max(bw.peak.cov[[seqname]][start:end])
    }
            ,mc.cores=8
            ,as.character(seqnames),start,end)
})
max.cov <- data.frame(max.cov)
head(max.cov)
gpeaks$score <- max.cov[,1]*scoreMultiplier

#gpeaks$length <- width(gpeaks)

head(gpeaks)

df <- data.frame(seqnames=seqnames(gpeaks),
                 starts=start(gpeaks),
                 ends=end(gpeaks),
                 names=names(gpeaks),
                 scores=score(gpeaks),
                 strands=strand(gpeaks)
     )

bedfile <- gsub(".bed$","",bedfile)
bwLabel <- gsub(".bw$","",basename(bwfile))
bedfile <- paste(bedfile,".score",bwLabel,".bed",sep="")
print(bedfile)

write.table(df,file=bedfile,quote=F,sep="\t",row.names=F,col.names=F)
