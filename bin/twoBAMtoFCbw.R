args <- commandArgs()

assembly <-sub('--assembly=', '', args[grep('--assembly=', args)])
extLen <-sub('--extLen=', '', args[grep('--extLen=', args)])
bam1 <-sub('--bam1=', '', args[grep('--bam1=', args)])
bam2 <-sub('--bam2=', '', args[grep('--bam2=', args)])

if (identical(extLen,character(0))){
  extLen <- 150
} else {
    extLen <- as.integer(extLen)
}

library(GenomicAlignments)
library(Rsamtools)
library(rtracklayer)
library(GenomicRanges)

print(assembly)
if ((assembly == "hg19") || (assembly == "hg38")) { organismStr <- "Hsapiens" }
if ((assembly == "mm9") || (assembly == "mm10")) { organismStr <- "Mmusculus" }
if (assembly == "sacCer3") organismStr <- "Scerevisiae"
if (assembly == "dm3") organismStr <- "Dmelanogaster"
print(organismStr)

assemblyLibrary <- paste("BSgenome.", organismStr, ".UCSC.", assembly, sep="")
print(assemblyLibrary)

library(assemblyLibrary,character.only=TRUE)

if ((assembly == "hg19") || (assembly == "hg38")) { organism <- Hsapiens }
if ((assembly == "mm9") || (assembly == "mm10")) { organism <- Mmusculus }
if (assembly == "sacCer3") organism <- Scerevisiae
if (assembly == "dm3") organism <-Dmelanogaster

twobam2bw <- function(BF1,BF2,organism){
    cat("opening:", BF1, sep="\n")
    bd1 <- readGAlignments(BF1)
    print (length(bd1))
    seqlevels(bd1,force=TRUE) <- seqlevels(bd1)[grep("_",seqlevels(bd1), invert=TRUE)]
    seqlevels(bd1,force=TRUE) <- seqlevels(bd1)[grep("EBV",seqlevels(bd1), invert=TRUE)]
    seqlevels(bd1,force=TRUE) <- seqlevels(bd1)[grep("chrM",seqlevels(bd1), invert=TRUE)]
#    seqlevels(bd1,force=TRUE) <- seqlevels(bd1)[grep("chr17",seqlevels(bd1), invert=TRUE)]
    cat("convert to GRanges\n")
    mygr1 <- trim(as(bd1,"GRanges"))
    if (extLen > 0){
        print(paste("extending reads to ",extLen,sep=""))
        mygr1 <- resize(mygr1, extLen)
    }
    cat("getting coverage\n")
    cov1 <- coverage(mygr1)
    rpm1 <- cov1*(1e6/length(bd1))
    seqlengths(rpm1) <- seqlengths(organism)[names(rpm1)]
    
    cat("opening:", BF2, sep="\n")
    bd2 <- readGAlignments(BF2)
    print (length(bd2))
    seqlevels(bd2,force=TRUE) <- seqlevels(bd2)[grep("_",seqlevels(bd2), invert=TRUE)]
    seqlevels(bd2,force=TRUE) <- seqlevels(bd2)[grep("EBV",seqlevels(bd2), invert=TRUE)]
    seqlevels(bd2,force=TRUE) <- seqlevels(bd2)[grep("chrM",seqlevels(bd2), invert=TRUE)]
#    seqlevels(bd2,force=TRUE) <- seqlevels(bd2)[grep("chr17",seqlevels(bd2), invert=TRUE)]
    cat("convert to GRanges\n")
    mygr2 <- trim(as(bd2,"GRanges"))
    if (extLen > 0){
        print(paste("extending reads to ",extLen,sep=""))
        mygr2 <- resize(mygr2, extLen)
    }
    cat("getting coverage\n")
    cov2 <- coverage(mygr2)
    rpm2 <- cov2*(1e6/length(bd2))
    seqlengths(rpm2) <- seqlengths(organism)[names(rpm2)]


    cat("Calculating fold change\n")
    fc <- rpm1/rpm2
    print(fc)
    fc.gr <- as(log(fc), "GRanges")
    #print(fc.gr)
    logfc <- fc.gr[-which(is.infinite(fc.gr$score)),]
    #print(head(logfc))
    logfc <- logfc[!is.na(logfc$score)]
    #print(head(logfc))

     ## export rpm to bigWig
     sample1 <- basename(gsub(".bam", "", BF1))
     sample2 <- basename(gsub(".bam", "", BF2))
     #logfc <- log(fc)
     #print(logfc)
    outfile <- paste(sample1,sample2,"logFC.bw", sep=".")
     cat(paste("exporting to bigwig", outfile, "\n", sep="\t"))
     export.bw(logfc, outfile)
     cat("export complete:", outfile, sep="\n")
}


# for each element of our vector, call the bam2bw function
#mclapply(BF1=bam1,BF2=bam2,organism=organism,twobam2bw,mc.cores=1,mc.preschedule=FALSE)
#mclapply(bam2,organism=organism,bam2bw,mc.cores=1,mc.preschedule=FALSE)
twobam2bw(bam1,bam2,organism)
