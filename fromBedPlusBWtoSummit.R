library(GenomicRanges)
library(rtracklayer)

args <- commandArgs()

bwfile <- sub('--bwfile=','',args[grep('--bwfile=',args)])
bedfile <- sub('--bedfile=','',args[grep('--bedfile=',args)])

peaks <- read.delim(file=bedfile,header=FALSE,sep="\t")
head(peaks)
dim(peaks)
if (dim(peaks)[2] == 3){
# Ignore any columns beyond the first three.
    colnames(peaks)<-c("chr","start","end")
    peaks$name <- paste(basename(bedfile),rownames(peaks),sep="_")
#    peaks$chr <- sub(" ","",peaks$chr)
} else if (dim(peaks)[2] == 4){
    colnames(peaks)<-c("chr","start","end","score")
    peaks$name <- paste(basename(bedfile),rownames(peaks),sep="_")
} else if (dim(peaks)[2] == 5){
    colnames(peaks)<-c("chr","start","end","name","score")
}
rownames(peaks) <- peaks$name
gpeaks <- as(peaks,"GRanges")
head(gpeaks)

cat("Remove chromosomes containing underscores from bedfile\n")
seqlevels(gpeaks,force=TRUE) <- seqlevels(gpeaks)[grep("_",seqlevels(gpeaks), invert=TRUE)]

cat("importing:", bwfile, sep="\n")
bw.peak <- import.bw(bwfile,RangedData=FALSE,selection = BigWigSelection(gpeaks))
cat("calc coverage\n")
bw.peak.cov <- coverage(bw.peak,weight='score')

#cat("get coverage for peak region\n")
#mean.cov <- with(as.data.frame(gpeaks),{       
#    mcmapply(function(seqname,start,end){
#        sum(bw.peak.cov[[seqname]][start:end])
#    }
#            ,mc.cores=8
#            ,as.character(seqnames),start,end)
#})
#mean.cov <- data.frame(mean.cov)
#gpeaks$mean.cov <- mean.cov[,1]/width(gpeaks)
#head(mean.cov)

cat("get max coverage for peak region\n")
max.cov <- with(as.data.frame(gpeaks),{       
    mcmapply(function(seqname,start,end){
        max(bw.peak.cov[[seqname]][start:end])
    }
            ,mc.cores=8
            ,as.character(seqnames),start,end)
})
max.cov <- data.frame(max.cov)
head(max.cov)
gpeaks$max.cov <- max.cov[,1]

cat("get summit for peak region\n")
summit.cov <- with(as.data.frame(gpeaks),{       
    mcmapply(function(seqname,start,end){
        which.max(bw.peak.cov[[seqname]][start:end])
    }
            ,mc.cores=8
            ,as.character(seqnames),start,end)
})
cat("creating data frame of summits")
summit.cov <- data.frame(summit.cov)
head(summit.cov)
gpeaks$summit <- summit.cov[,1]+start(gpeaks)

#gpeaks$length <- width(gpeaks)

head(gpeaks)

df <- data.frame(seqnames=seqnames(gpeaks),
                 starts=gpeaks$summit-1,
                 ends=gpeaks$summit,
                 names=names(gpeaks),
                 scores=score(gpeaks),
                 strands=c(rep(".", length(gpeaks))))

filename <- sub(".bed$","",bedfile)
filename <- paste(filename,".summits.bed",sep="")
filename
head(df)
write.table(df,file=filename,quote=F,sep="\t",row.names=F,col.names=F)
