library(Rsamtools)
library(GenomicAlignments)
library(parallel)
library(GenomicFeatures)
library(rtracklayer)
library(biomaRt)
library(GenomicRanges)
library(gplots)

args <- commandArgs()

bedFile <-sub('--bedFile=', '', args[grep('--bedFile=', args)])

print("bedFile")
bedFile

peaks <- read.delim(file=bedFile,header=FALSE,sep="\t")
head(peaks)
# Ignore any columns beyond the first five.
peaks <- peaks[,1:3]
colnames(peaks)<-c("chr","start","end")
peaks$chr <- sub(" ","",peaks$chr)
#head(peaks)
gpeaks <- as(peaks,"GRanges")

gpeaks$lengths <- width(gpeaks)
head(gpeaks)

pdf(paste(bedFile,"lengths.pdf",sep="."))
hist(gpeaks$lengths)
dev.off()
