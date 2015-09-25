args <- commandArgs()

assembly <-sub('--assembly=', '', args[grep('--assembly=', args)])
bamDir <-sub('--bamDir=', '', args[grep('--bamDir=', args)])
sample <- sub('--sample=', '', args[grep('--sample=', args)])
sample <- sub('$', '.bam', sample)

library(GenomicAlignments)
library(Rsamtools)
library(rtracklayer)
library(GenomicRanges)

print(assembly)
if (assembly == "hg19") organismStr <- "Hsapiens"
if (assembly == "mm9") organismStr <- "Mmusculus"
if (assembly == "mm10") organismStr <- "Mmusculus"
if (assembly == "sacCer3") organismStr <- "Scerevisiae"
if (assembly == "dm3") organismStr <- "Dmelanogaster"
print(organismStr)

assemblyLibrary <- paste("BSgenome.", organismStr, ".UCSC.", assembly, sep="")
print(assemblyLibrary)

library(assemblyLibrary,character.only=TRUE)

if (assembly == "hg19") organism <- Hsapiens
if (assembly == "mm9") organism <- Mmusculus
if (assembly == "mm10") organism <- Mmusculus
if (assembly == "sacCer3") organism <- Scerevisiae
if (assembly == "dm3") organism <- Dmelanogaster

#setwd to bamfile dir 
setwd(bamDir)

#----------------------------------load bamfiles-----------------------------------#
BFL <- list.files(pattern = '.bam$')
BFL <- grep(pattern=sample,BFL,value=TRUE)
BFL

param <- ScanBamParam(what='mapq',tag='NH')

bam2bw <- function(BF,organism){
    cat("opening:", BF, sep="\n")
    bd <- readGAlignments(BF,param=param)
    bd <- bd[!values(bd)$NH > 1] 
    seqlevels(bd,force=TRUE) <- seqlevels(bd)[grep("_",seqlevels(bd), invert=TRUE)]
    cat("change the strand:", BF, sep="\n")
    strand(bd) <- ifelse(strand(bd) == '+', '-', '+')
    ## plus coverage
    cat("getting coverage for plus strand\n")
    plus <- coverage(bd[strand(bd)=="+",])
    plus.rpm <- plus*(1e6/length(bd)) ## normalizing to the library size
    seqlengths(plus.rpm) <- seqlengths(organism)[names(plus.rpm)]
    ## minus coverage
    cat("getting coverage for plus minus\n")
    minus <- coverage(bd[strand(bd)=="-",])
    minus.rpm <- minus*(-1e6/length(bd)) ## normalizing to the library size
    seqlengths(minus.rpm) <- seqlengths(organism)[names(minus.rpm)]
    ## export rpm to bigWig
    plus.outfile <- sub("bam", "plus", BF)
    minus.outfile <- sub("bam", "minus", BF)
    plus.outfile <- paste(plus.outfile, ".bw", sep="")
    minus.outfile <- paste(minus.outfile, ".bw", sep="")
    cat(paste("exporting to plus bigwig", plus.outfile, "\n", sep="\t"))
    export.bw(plus.rpm, plus.outfile)
    cat(paste("exporting to minus bigwig", minus.outfile, "\n", sep="\t"))
    export.bw(minus.rpm, minus.outfile)
    cat("export complete:", BF, sep="\n")
}

# for each element of our vector, call the bam2bw function
mclapply(BFL,organism=organism,bam2bw,mc.cores=1,mc.preschedule=FALSE)
