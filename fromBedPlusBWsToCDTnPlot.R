library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(Rsamtools)
library(rtracklayer)
library(ChIPpeakAnno)
library(gplots)

args <- commandArgs()

bedFile <-sub('--bedFile=', '', args[grep('--bedFile=', args)])
bwlist <- sub('--bwlist=','',args[grep('--bwlist=',args)])
numBins <- 25
boundaryNum <- 10

print (paste("Read in bed file",bedFile,sep=" "))
peaks <- import.bed(bedFile)
peaks$ensembl_gene_id = gsub("\\.ENST\\d+$","",peaks$name)
head(peaks)

bws <- read.table(bwlist,header=FALSE,sep="\t")
print(bws)
sampleNum <- dim(bws)[1]
if (dim(bws)[2]>2){
    colors <- bws[,2]
} else {
    n <- sampleNum
    colors <-rainbow(n, s = 1, v = 1, start = 0, end = max(1, n - 1)/n, alpha = 1)
}
bws <- as.character(bws[,1])
print(bws)
print(colors)

report <- data.frame(sample.name=sub(".bw","",basename(bws))
                    ,bw=as.character(bws))
report$order <- sub('.*\\-|.*\\.',"",dirname(as.character(report$bw)))
report$order <- sub("/bam","",report$order)
#report$rename <- gsub(".bw$",".bam",report$bw)
report$rename <- gsub(".bw$","",report$bw)
report$rename <- basename(report$rename)
#report$rename <- paste(gsub("\\-","_",report$sample),report$order,sep="_")
print(report)

matBin <-function(peakdf,model){
### read in and filter peaks
    sapply(unique(peakdf$rename),function(x,peakdf=report){
        NCores=20
#        fname <- paste(Dir,x,".rda", sep="")
#        fname2 <- paste(Dir,x,".txt", sep="")
        cat("importing:", x, sep="\n")
        bw <- as.character(peakdf$bw[peakdf$rename %in% x])
        print(bw)
#        bw <- as.character(peakdf$bw)
        bw.peak <- import.bw(bw,RangedData=FALSE,selection = BigWigSelection(model))
        cat("calc coverage\n")
        bw.peak.cov <- coverage(bw.peak,weight='score')
        cat("get coverage for peak region\n")
        cov <- with(as.data.frame(model),{
            mcmapply(function(seqname,start,end,strand){
                r <- bw.peak.cov[[seqname]][start:end]
                if(strand == '-'){r <- rev(r)}
                return(r)
            }
                    ,mc.cores=NCores
                    ,as.character(seqnames),start,end,as.character(strand))
        })
        cat("convert list to matrix\n")
        mat <- do.call(rbind, mclapply(cov, as.numeric,mc.cores=NCores))
        cov <- data.frame(mat)
        cat("bin the matrix\n")
        window.cov <- function(row){
              window <- as.integer(ncol(cov)/numBins)
              window.coverage <- lapply(0:(window-1), function(jump)
                  rowMeans(row[(jump*(numBins-1)+1):(jump*(numBins-1)+1)+(numBins-1)])
                                        )
              t(as.matrix(unlist(window.coverage)))
          }
        win <- mclapply(1:nrow(cov), function(i)
            window.cov(cov[i,]),mc.cores=NCores) 
        bin.mat <- do.call(rbind, mclapply(win, as.numeric, mc.cores=NCores))
        df <- data.frame(bin.mat)
        rownames(df) <- model$ensembl_gene_id
        head(df)
        df
#        save(df,file=fname)
#       write.table(df,file=fname2,sep="\t")
        assign(x, df,envir=.GlobalEnv)
    })
}

names(peaks) <- NULL
print ("Running matBin for peaks")
matBin(peakdf=report,model=peaks)

print("Create CDT file")
cdt <- cbind(UID=peaks$ensembl_gene_id
            ,NAME=peaks$ensembl_gene_id
             )
for(i in 1:length(report$rename)){
    sample <- report$rename[i]
    x <- get(sample)
    print(sample)
    cdt <- cbind(cdt,sample=x,boundary=c(rep(boundaryNum, dim(cdt)[1])))
    colnames(cdt) <- gsub("^sample",sample,colnames(cdt))
}
cdt <- na.omit(cdt)
#head(cdt)
bedFile <- gsub(".bed$","",bedFile)
write.table(cdt, file=paste(bedFile,"cdt",sep="."), sep="\t",row.names=FALSE)

cdtData <- cdt[,3:dim(cdt)[2]]
#head(cdtData)

binLength <- dim(cdtData)[2]/sampleNum
dim(cdtData)
sampleNum
binLength
numBins
for (i in 1:dim(cdtData)[2]){
    if (i %% binLength == binLength-2){
        colnames(cdtData)[i] <- gsub(".X\\d+$","",colnames(cdtData)[i])
        colnames(cdtData)[i-(binLength/2)] <- colnames(cdtData)[i]
        colnames(cdtData)[i] <- ""
    } else {
        colnames(cdtData)[i] <- ""
    }
}

lmat = rbind(c(0,0),c(0,1))
lhei = c(0.1,4)
lwid = c(0.1,4)
pdfFile <- paste(bedFile,"heatmap.pdf",sep=".")
print(pdfFile)
pdf(pdfFile)
heatmap.2(as.matrix(cdtData),trace="none",dendrogram="none",scale="none",Rowv=NA,Colv=NA,labRow=rownames(cdtData),cexRow=1,colsep=0,rowsep=0,cexCol=1,margins=c(12,9),main=bedFile,key=FALSE,lmat=lmat,lwid=lwid,lhei=lhei,col=colorRampPalette(c("white","darkorchid4"))(75),add.expr=abline(v=0:sampleNum*binLength))
dev.off()

summary <- colMeans(cdtData)

justData <- summary[1:binLength-1]
for (i in 1:sampleNum-1){
    for (j in 1:binLength){
        if (j%%binLength != 0){
#            print(i*binLength+j)
 #           print(summary[i*binLength+j])
            justData <- c(justData,summary[i*binLength+j])
        }
    }
}
binLength <- binLength-1
ylabel <- "BinNumber"
bwlist <- sub('--bwlist=','',args[grep('--bwlist=',args)])
#flanks <- basename(bedFile[grep('.{1}\\d+{1,5}.{1}\\d+{1,5}.{1}tss',basename(bedFile))])
#grep(".\\d+.\\d+.tss", basename(bedFile), ignore.case = FALSE, perl = FALSE,value=TRUE)
#flanks

pdfFile <- paste(bedFile,"mean.pdf",sep=".")
print(pdfFile)
pdf(pdfFile)
plot(justData[0:binLength],type="l",col=colors[0],ylim=c(0,1.5*max(justData)),ylab="ChIP RPM",xlab="BinNumber",main=paste(basename(bedFile)," (",length(peaks),")",sep=""))
for (i in 1:sampleNum){
    dataStart <- i*binLength+1
    dataEnd <- (i+1)*binLength
 #   print(i)
 #   print(dataStart)
 #   print(dataEnd)
    lines(justData[dataStart:dataEnd],col=colors[i])
}
legend("topleft",report$rename[1:sampleNum],col=colors,lty=1)
dev.off()   

