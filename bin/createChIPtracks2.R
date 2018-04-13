args <- commandArgs()

assembly <-sub('--assembly=', '', args[grep('--assembly=', args)])
bamDir <-sub('--bamDir=', '', args[grep('--bamDir=', args)])
extLen <-sub('--extLen=', '', args[grep('--extLen=', args)])
sample <- sub('--sample=', '', args[grep('--sample=', args)])
sample <- sub('$', '.bam', sample)
multiMap <- sub('--multiMap=', '', args[grep('--multiMap=',args)])

if (identical(multiMap,character(0))){
    multiMap <- 0
}
multiMap <- as.integer(multiMap)

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
if (assembly == "rn6") organismStr <- "Rnorvegicus"
print(organismStr)

assemblyLibrary <- paste("BSgenome.", organismStr, ".UCSC.", assembly, sep="")
print(assemblyLibrary)

library(assemblyLibrary,character.only=TRUE)

if ((assembly == "hg19") || (assembly == "hg38")) { organism <- Hsapiens }
if ((assembly == "mm9") || (assembly == "mm10")) { organism <- Mmusculus }
if (assembly == "sacCer3") organism <- Scerevisiae
if (assembly == "dm3") organism <-Dmelanogaster
if (assembly == "rn6") organism <- Rnorvegicus

#setwd to bamfile dir 
setwd(bamDir)

#----------------------------------load bamfiles-----------------------------------#
BFL <- list.files(pattern = '.bam$')
BFL <- grep(pattern=sample,BFL,value=TRUE)
BFL

param <- ScanBamParam(what='mapq',tag='NH')
multiString <- ""
if (multiMap == 1){
    multiString <- ".multi"
    print("Keeping multiply mapped reads.")
} else if (multiMap == 0){
    print("Discarding multiply mapped reads.")
}

bam2bw <- function(BF,organism){
      cat("opening:", BF, sep="\n")
      bd <- readGAlignments(BF)
      print (length(bd))
      if (multiMap == 0){
          bd <- bd[!values(bd)$NH > 1]
          print("Multi-mapped hits are removed")
          print (length(bd))
      }      
      seqlevels(bd,force=TRUE) <- seqlevels(bd)[grep("_",seqlevels(bd), invert=TRUE)]
      seqlevels(bd,force=TRUE) <- seqlevels(bd)[grep("EBV",seqlevels(bd), invert=TRUE)]
      cat("convert to GRanges\n")
      mygr <- as(bd,"GRanges")
      if (extLen > 0){
          print(paste("extending reads to",extLen,sep=""))
          mygr <- resize(mygr, extLen)
      }
      cat("getting coverage\n")
      # get coverage                                                            
      cov <- coverage(mygr)
      rpm <- cov*(1e6/length(bd))
      seqlengths(rpm) <- seqlengths(organism)[names(rpm)]
      ## export rpm to bigWig
      outfile <- sub(".bam", "", BF)
      outfile <- paste(outfile,multiString, ".bw", sep="")
      cat(paste("exporting to bigwig", outfile, "\n", sep="\t"))
      export.bw(rpm, outfile)
      cat("export complete:", BF, sep="\n")
  }

# for each element of our vector, call the bam2bw function
mclapply(BFL,organism=organism,bam2bw,mc.cores=1,mc.preschedule=FALSE)
