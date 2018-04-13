args <- commandArgs()

assembly <-sub('--assembly=', '', args[grep('--assembly=', args)])
bamDir <-sub('--bamDir=', '', args[grep('--bamDir=', args)])
extLen <-sub('--extLen=', '', args[grep('--extLen=', args)])
sample <- sub('--sample=', '', args[grep('--sample=', args)])
print(sample)
sample <- sub('$', '.bam', sample)
sample <- sub('.bam.bam', '.bam', sample)
print(sample)

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

#setwd to bamfile dir 
setwd(bamDir)

#----------------------------------load bamfiles-----------------------------------#
BFL <- list.files(pattern = '.bam$')
BFL <- grep(pattern=sample,BFL,value=TRUE)
BFL


bam2bw <- function(BF,organism){
      cat("opening:", BF, sep="\n")
      bd <- readGAlignments(BF)
      seqlevels(bd,force=TRUE) <- seqlevels(bd)[grep("_",seqlevels(bd), invert=TRUE)]
      cat("convert to GRanges\n")
      mygr <- as(bd,"GRanges")
      cat("extending reads\n")
      if (extLen > 0){
          mygr <- resize(mygr, extLen)
      }
      cat("getting coverage\n")
      # get coverage                                                             
      cov <- coverage(mygr)
      rpm <- cov*(1e6/length(bd))
      seqlengths(rpm) <- seqlengths(organism)[names(rpm)]
      ## export rpm to bigWig
      outfile <- sub(".bam", "", BF)
      outfile <- paste(outfile, ".bw", sep="")
      cat(paste("exporting to bigwig", outfile, "\n", sep="\t"))
      export.bw(rpm, outfile)
      cat("export complete:", BF, sep="\n")
  }

# for each element of our vector, call the bam2bw function
mclapply(BFL,organism=organism,bam2bw,mc.cores=1,mc.preschedule=FALSE)
