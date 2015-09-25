args <- commandArgs()

topHatDir <-sub('--topHatDir=', '', args[grep('--topHatDir=', args)])

nClus <-sub('--nClus=', '', args[grep('--nClus=', args)])

print(topHatDir)

topHatReport <- function(topHatDir,nClus){
    library(GenomicAlignments)
    library(Rsamtools)
    library(parallel)
    topHat.mapped <- list.files(path=topHatDir,pattern="^accepted_hits.bam$",recursive=TRUE,full.names=TRUE)
    topHat.unmapped <- file.path(dirname(topHat.mapped),'unmapped.bam')
    stopifnot(all(file.exists(topHat.unmapped)))  ## Breaks if there is nounmapped bam associated with a bam file
    getNumber <- function(mapped.file){
      unmapped.file <- file.path(dirname(mapped.file),'unmapped.bam')
      if(!file.exists(unmapped.file)){ return(NULL) }
      param <- ScanBamParam(flag = scanBamFlag(isUnmappedQuery=NA)
                            ,tag='NH')
      t <- readGAlignments(mapped.file,
                                param=param)
      tt <- table(values(t)$NH)
      tt <- tt/as.numeric(names(tt))
      counts <- c(unique=tt[as.numeric(names(tt))==1]/1e6,
                  multi.match=sum(tt[as.numeric(names(tt))>1])/1e6,
                  unmapped=countBam(unmapped.file)$records/1e6)
c(sprintf("%0.1fM",c(sum(counts),counts)),sprintf("%0.1f%%",counts/sum(counts)*100))
    }
    counts <-t(data.frame(mclapply(topHat.mapped,getNumber,mc.cores=nClus,mc.preschedule=FALSE)))
    rownames(counts) <- sub(".+/","",dirname(topHat.mapped))
    colnames(counts) <- c("total reads",
                          "singly mapped",
                          "multiply mapped",
                          "unmapped",
                          "singly mapped (%)",
                          "multiply mapped (%)",
                          "unmapped reads (%)")
    write.table(counts, file=paste(topHatDir,"aln_report.txt",sep="/"), sep="\t", col.names=NA,quote=F)
    counts
  }


topHatReport(topHatDir,nClus)
