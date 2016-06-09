args <- commandArgs()

assembly <-sub('--assembly=', '', args[grep('--assembly=', args)])
bamDir <-sub('--bamDir=', '', args[grep('--bamDir=', args)])
txdbfile <-sub('--txdbfile=', '',args[grep('--txdbfile=',args)])
numCores <- as.integer(sub('--numCores=', '', args[grep('--numCores=',args)]))
multiMap <- sub('--multiMap=', '', args[grep('--multiMap=',args)])

#assembly = "dm3"
#txdbfile = "/projects/p20742//anno/Txdb/dmelanogaster_gene_ensembl_Ens74.txdb"
#numCores = 4
#multiMap = 0
#bamDir = "/projects/b1025/etb/ryan//TANGO-046/bam/"

if (identical(multiMap,character(0))){
    multiMap <- 0
}
multiMap <- as.integer(multiMap)

library(Rsamtools)
library(GenomicAlignments)
library(parallel)
library(GenomicFeatures)
library(rtracklayer)
library(biomaRt)
library(GenomicRanges)

print(assembly)
if ((assembly == "hg19") || (assembly == "hg38")) { organismStr <- "Hsapiens" }
if ((assembly == "mm9") || (assembly == "mm10")) { organismStr <- "Mmusculus" }
if (assembly == "sacCer3") organismStr <- "Scerevisiae"
if (assembly == "dm3") organismStr <- "Dmelanogaster"

assemblyLibrary <- paste("BSgenome.", organismStr, ".UCSC.", assembly, sep="")
print(assemblyLibrary)

library(assemblyLibrary,character.only=TRUE)

if ((assembly == "hg19") || (assembly == "hg38")) { organism <- Hsapiens }
if ((assembly == "mm9") || (assembly == "mm10")) { organism <- Mmusculus }
if (assembly == "sacCer3") organism <- Scerevisiae
if (assembly == "dm3") organism <-Dmelanogaster

txdb <- loadDb(txdbfile)
txdb
seqlevels(txdb,force=TRUE) <- seqlevels(txdb)[grep("_|\\d+.1$",seqlevels(txdb), invert=TRUE)]
seqlevels(txdb) <- sub("","chr", seqlevels(txdb))
seqlevels(txdb) <- sub("MT","M", seqlevels(txdb))
seqlevels(txdb) <- sub("Mito","M", seqlevels(txdb))

#parameters that must be set before the counter or it will only do the first BF
mapq=10 #map qs of 10
gnModel <- transcriptsBy(txdb, 'gene')
gnModel <- unlist(range(gnModel))## Gets the genomic region convered by transcripts
seqinfo(gnModel) <- seqinfo(organism)[seqlevels(gnModel)]

#move to directory containing bam files and gather the names of the bam files
setwd(bamDir)
BFL<-dir(pattern='.*.bam$',recursive=FALSE)
BFL

param <- ScanBamParam(what='mapq',tag='NH')

counter <-function(BF,txdb,mapq=10){
    aln <- readGAlignments(BF,param=param)
    strand(aln) <- ifelse(strand(aln) == '+', '-', '+') #stranded library
    ovlp <- countOverlaps(aln, gnModel)          #Count how many genes each alignment has
    aln <- aln[ovlp==1 ]                         # Remove reads overlapping more than one gene
    if (multiMap == 0){
        aln <- aln[!values(aln)$NH > 1]              # Removing multimapped reads
    }
    aln <- aln[!values(aln)$mapq <= mapq]        # Removing low map quality reads
    counts <- countOverlaps(gnModel, aln)  # Counting how many reads overlap each gene
    names(counts) <- names(gnModel)        # Making sure we have the gene names
    counts                                 # Return counts for each bam
}


counts <- do.call(cbind,mclapply(BFL,counter,txdb=gnModel,mc.cores=numCores,mc.preschedule=FALSE))

colnames(counts) <-  sub(".bam","",basename(BFL))

gnModel.df <- as.data.frame(gnModel)
head(gnModel.df)
save(gnModel.df, file=file.path("gnModel.rda"))
write.table(gnModel.df, file=file.path("gnModel.txt"),sep="\t")

counts <- data.frame(counts)
counts <- cbind(gnModel.df,counts)

countFile<-"counts"
if (multiMap == 1){ countFile<-"counts.multi" }
    
## save the data to disk
save(counts, file=paste(bamDir,"/",countFile,".rda",sep=""))
write.table(counts, file=paste(bamDir,"/",countFile,".txt",sep=""),sep="\t")
