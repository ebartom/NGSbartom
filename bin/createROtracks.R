args <- commandArgs()

help <- function(){
    cat("createROtracks.R :
- Create tracks for run on sequencing such as GROseq and PROseq or TT-seq\n")
    cat("Usage: \n")
    cat("--bamFile    : Sample bam file (must be a full path)                           [required]\n")    
    cat("--outName    : Prefix for your output file                                     [default = sample]\n")    
    cat("--assembly   : Genome build hg19, mm9, mm10, dm3 ect.                          [default = hg19]\n")
    cat("--sepStrands : If the strands should be separated or kept in one track (0/1)   [default = 1]\n")
    cat("--flipStrand : If the strand should be flipped. Depends where the primer is.   [default = 1]\n")
    cat("               (typically GROseq no ;TTseq, RNAseq, PROseq yes) (0/1)                       \n")
    cat("--extLen     : number of bases to extend the reads to                          [default = 0]\n")
    cat("--noStrand   : replace strand with * (ChIP-seq) if set ignore flipStrand (0/1) [default = 0]\n")
    cat("--threePrime : Only report 1 position at the 3' end of the read (PRO-seq)(0/1) [default = 0]\n")
    cat("--fivePrime  : Only report 1 position at the 5' end of the read (0/1)          [default = 0]\n")
    cat("--NHflag     : If bamFile has NH flag (tophat, bowtie2 or  star aligner (0/1)  [default = 0]\n")
    cat("--multiMap   : Number of reported hits to keep from NH flag (0/1)              [default = 0]
                         (only use if NHflag = 1 ) can not be used with bowtie1 as no
                         NH flag exists. In tophat if you use -g 10, it reports up to
                         10 mappings, but could actually be more not reported. You need
                         to report 1 more mapping than you want to keep!                             \n")    
    cat("--SpikeIn    : Fraction to normalize (multiply) reads by.                      [default = 1]
                         This fraction is control spike-in total reads over experiment.
                         If your control spikein total reads are 3549426 and experiment
                         are 4812918, your fraction is 0.737479. 1 means that is no
                         normalization.                                                              \n")    
    cat("\n")
    q()
}

## Save values of each argument
if( !is.na(charmatch("-h",args)) || !is.na(charmatch("-help",args))){
    help()
} else {
    bamFile     <- sub('--bamFile=', '',    args[grep('--bamFile=', args)])
    outName     <- sub('--outName=', '',    args[grep('--outName=', args)])
    assembly    <- sub('--assembly=', '',   args[grep('--assembly=', args)])
    sepStrands  <- sub('--sepStrands=', '', args[grep('--sepStrands=', args)])
    flipStrand  <- sub('--flipStrand=', '', args[grep('--flipStrand=', args)])
    extLen      <- sub('--extLen=', '',     args[grep('--extLen=', args)])
    noStrand    <- sub('--noStrand=', '',   args[grep('--noStrand=', args)])
    threePrime  <- sub('--threePrime=', '', args[grep('--threePrime=', args)])
    fivePrime   <- sub('--fivePrime=', '',  args[grep('--fivePrime=', args)])
    NHflag      <- sub('--NHflag=', '',     args[grep('--NHflag=', args)])
    multiMap    <- sub('--multiMap=', '',   args[grep('--multiMap=',args)])
    SpikeIn     <- sub('--SpikeIn=', '',    args[grep('--SpikeIn=', args)])
}

if (identical(outName,character(0))){
   outName <- sub(".bam$", "", bamFile)
}

if (identical(assembly,character(0))){
   assembly <- "hg19"
}

if (identical(sepStrands,character(0))){
   sepStrands <- 1
}else{
   sepStrands <- as.numeric(sepStrands)
}

if (identical(flipStrand,character(0))){
    flipStrand <- 1
}else{
    flipStrand <- as.numeric(flipStrand)
}

if (identical(extLen,character(0))){
   extLen <- 0
}else{
   extLen <- as.numeric(extLen)
}

if (identical(noStrand,character(0))){
   noStrand <- 0
}else{
   noStrand <- as.numeric(noStrand)
}

if (identical(threePrime,character(0))){
   threePrime <- 0
}else{
   threePrime <- as.numeric(threePrime)
}

if (identical(fivePrime,character(0))){
   fivePrime <- 0
}else{
   fivePrime <- as.numeric(fivePrime)
}

if (identical(NHflag,character(0))){
   NHflag <- 0
}else{
   NHflag <- as.numeric(NHflag)
}

if (identical(multiMap,character(0))){
    multiMap <- 0
}else{
    multiMap <- as.numeric(multiMap)
}

if (identical(SpikeIn,character(0))){
   SpikeIn <- 1
}else{
    SpikeIn <- as.numeric(SpikeIn)
}

cat("Running script with the following options:", sep="\n")
cat("bamFile:",    bamFile,    sep="\n")
cat("outName:",    outName,    sep="\n")
cat("assembly:",   assembly,   sep="\n")
cat("flipStrand:", flipStrand, sep="\n")
cat("extLen:",     extLen,     sep="\n")
cat("noStrand:",   noStrand,   sep="\n")
cat("threePrime:", threePrime, sep="\n")
cat("fivePrime:",  fivePrime,  sep="\n")
cat("NHflag:",     NHflag,     sep="\n")
cat("multiMap:",   multiMap,   sep="\n")
cat("SpikeIn:",    SpikeIn,    sep="\n")

library(GenomicAlignments)
library(Rsamtools)
library(rtracklayer)
library(GenomicRanges)

cat(assembly, sep="\n")
if ((assembly == "hg19") || (assembly == "hg38")) { organismStr <- "Hsapiens" }
if ((assembly == "mm9") || (assembly == "mm10")) { organismStr <- "Mmusculus" }
if (assembly == "sacCer3") organismStr <- "Scerevisiae"
if (assembly == "dm3") organismStr <- "Dmelanogaster"
cat(organismStr, sep="\n")

assemblyLibrary <- paste("BSgenome.", organismStr, ".UCSC.", assembly, sep="")
cat(assemblyLibrary, sep="\n")

library(assemblyLibrary,character.only=TRUE)

if ((assembly == "hg19") || (assembly == "hg38")) { organism <- Hsapiens }
if ((assembly == "mm9") || (assembly == "mm10")) { organism <- Mmusculus }
if (assembly == "sacCer3") organism <- Scerevisiae
if (assembly == "dm3") organism <-Dmelanogaster



bam2bw <- function(BF,organism){
    cat("opening:", BF, sep="\n")
    if( NHflag == 1 ){
        param                <- ScanBamParam(what='mapq',tag='NH')
        bd                   <- readGAlignments(BF,param=param)
        if (multiMap == 0){
            print("filter out multimapped reads") 
            bd               <- bd[!values(bd)$NH > 1]
        }else{
            bd               <- bd[!values(bd)$NH > multiMap]
            outName          <- sub("$", ".multi", outName)
        }
    }else{
        bd                   <- readGAlignments(BF)
    }
    cat("remove unassembled and viral chromosomes", sep="\n") 
    seqlevels(bd,force=TRUE) <- seqlevels(bd)[grep("_|EBV",seqlevels(bd), invert=TRUE)]
    print(seqlevels(bd))
    if( flipStrand > 0 & noStrand == 0){
        cat("change the strand", sep="\n")
        strand(bd)           <- ifelse(strand(bd) == '+', '-', '+')
    }
    if( noStrand > 0 ){
        strand(bd)           <- '*'
    }
    cat("convert to GRanges\n")
    mygr                    <- as(bd,"GRanges")
    if (extLen > 0){
        cat("extend reads", sep="\n")
        mygr                 <- resize(mygr, extLen, fix='start')
    }
    if (threePrime > 0){
        cat("take three prime position of the reads", sep="\n")
        mygr                 <- resize(mygr, width=1, fix='end')
        outName              <- sub("$", ".3prime", outName)
    }
    if ( fivePrime > 0 ){
        cat("take five prime position of the reads", sep="\n")
        mygr                 <- resize(mygr, width=1, fix='start')
        outName              <- sub("$", ".5prime", outName)
    }
    if (sepStrands > 0){        
        cat("getting coverage for separate strands\n")
        ## get plus coverage                                                             
        plus                  <- coverage(mygr[strand(mygr) == "+"])
        plus.rpm              <- plus*(1e6/length(bd))
        seqlengths(plus.rpm)  <- seqlengths(organism)[names(plus.rpm)]
        ## get minus coverage
        minus                 <- coverage(mygr[strand(mygr) == "-"])
        minus.rpm             <- minus*(-1e6/length(bd))
        seqlengths(minus.rpm) <- seqlengths(organism)[names(minus.rpm)]
        ## set the outfile name
        plus.outfile          <- sub("$", ".rpm.plus.bw", outName)
        minus.outfile         <- sub("$", ".rpm.minus.bw", outName)
        ## export rpm to bigWig
        cat(paste("exporting r.p.m. plus bigwig:", plus.outfile, sep="\n"), sep="\n")
        export.bw( plus.rpm, plus.outfile )
        cat(paste("exporting r.p.m. minus bigwig:", minus.outfile, sep="\n"), sep="\n")
        export.bw( minus.rpm, minus.outfile )      
        if( SpikeIn != 1 ){
            print("making SpikeIn normalized tracks")
            ## get plus coverage
            plus.si               <- plus*SpikeIn
            seqlengths(plus.si)   <- seqlengths(organism)[names(plus.si)]
            ## get minus coverage
            minus.si             <- minus*(-SpikeIn)
            seqlengths(minus.si) <- seqlengths(organism)[names(minus.si)]
            ## set the outfile name
            plus.si.outfile      <- sub("$", ".spikeNorm.plus.bw", outName)
            minus.si.outfile     <- sub("$", ".spikeNorm.minus.bw", outName)
            ## export rpm to bigWig
            cat(paste("exporting spikeNorm plus bigwig:", plus.si.outfile, sep="\n"), sep="\n")
            export.bw( plus.si, plus.si.outfile )
            cat(paste("exporting spikeNorm minus bigwig:", minus.si.outfile, sep="\n"), sep="\n")
            export.bw( minus.si, minus.si.outfile ) 
        }
        cat("export complete", sep="\n")
    }else{       
        cat("getting coverage for both strands\n")
        cov                  <- coverage(mygr)
        rpm                  <- cov*(1e6/length(bd))
        seqlengths(rpm)      <- seqlengths(organism)[names(rpm)]
        ## export rpm to bigWig
        outfile              <- sub("$", ".bw", outName)
        cat( paste("exporting r.p.m. bigwig:", outfile, sep="\n"), sep="\n")
        export.bw( rpm, outfile )
        if( SpikeIn != 1 ){
            cat("making SpikeIn normalized tracks", sep="\n")
            cov.si             <- cov*SpikeIn
            seqlengths(cov.si) <- seqlengths(organism)[names(cov.si)]
            ## export rpm to bigWig
            outfile.si         <- sub("$", ".spikeNorm.bw", outName)
            cat(paste("exporting spikeNorm bigwig:", outfile.si, sep="\n"), sep="\n")
            export.bw( cov.si, outfile.si )
        }
        cat("export complete", sep="\n")          
    }    
}

# for each element of our vector, call the bam2bw function
mclapply(bamFile,organism=organism,bam2bw,mc.cores=1,mc.preschedule=FALSE)
