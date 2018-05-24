library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(Rsamtools)
library(rtracklayer)
#library(ChIPpeakAnno)


args <- commandArgs()

bedFile <-sub('--bedFile=', '', args[grep('--bedFile=', args)])
bwdir <- sub('--bwdir=','',args[grep('--bwdir=',args)])
numBins <- 25

print (paste("Read in bed file",bedFile,sep=" "))
peaks <- import.bed(bedFile)
peaks$ensembl_gene_id = gsub("\\.ENST\\d+$","",peaks$name)
head(peaks)

bws <- list.files(bwdir,pattern=".bw$")
print(bws)


report <- data.frame(sample.name=sub(".bw","",basename(bws))
                    ,bw=as.character(bws))
report$order <- sub('.*\\-|.*\\.',"",dirname(as.character(report$bw)))
report$order <- sub("/bam","",report$order)
#report$rename <- gsub(".bw$",".bam",report$bw)
report$rename <- gsub(".bw$","",report$bw)
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
        bw.peak <- import.bw(paste(bwdir,bw,sep="/"),RangedData=FALSE,selection = BigWigSelection(model))
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
    cdt <- cbind(cdt,sample=x,boundary=c(rep(10, dim(cdt)[1])))
    colnames(cdt) <- gsub("^sample",sample,colnames(cdt))
}
cdt <- na.omit(cdt)
head(cdt)

bedFile <- gsub(".bed$","",bedFile)
write.table(cdt, file=paste(bedFile,"cdt",sep="."), sep="\t",row.names=FALSE)
