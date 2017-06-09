library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(Rsamtools)
library(rtracklayer)
library(biomaRt)
#library(BSgenome.Hsapiens.UCSC.hg19)
library(ChIPpeakAnno)

args <- commandArgs()

txdbfile <- sub('--txdbfile=','',args[grep('--txdbfile=', args)])
assembly <-sub('--assembly=', '', args[grep('--assembly=', args)])
up <- as.numeric(sub('--up=','',args[grep('--up=',args)]))
down <- as.numeric(sub('--down=','',args[grep('--down=',args)]))
geneList <-sub('--geneList=', '', args[grep('--geneList=', args)])
bwfile <- sub('--bwfile=','',args[grep('--bwfile=',args)])
keepOrder <- sub('--keepOrder=','',args[grep('--keepOrder=',args)])
up
down

if (identical(keepOrder,character(0))){
    keepOrder <- 0
}

print(assembly)

if ((assembly == "hg19") || (assembly == "hg38")) { organismStr <- "Hsapiens"
                      species <- "Homo sapiens"}
if ((assembly == "mm9") || (assembly == "mm10")) { organismStr <- "Mmusculus"
                                               species <- "Mus musculus"}
if (assembly == "sacCer3") { organismStr <- "Scerevisiae"
                             species <- "Saccharomyces cerevisiae"}
if (assembly == "dm3") { organismStr <- "Dmelanogaster"
                         species <- "Drosophila melanogaster"}
if (assembly == "rn6") { organismStr <- "Rnorvegicus"
                         species <- "Rattus norvegicus"}

assemblyLibrary <- paste("BSgenome.", organismStr, ".UCSC.", assembly, sep="")
print(assemblyLibrary)
library(assemblyLibrary,character.only=TRUE)

if ((assembly == "hg19") || (assembly == "hg38")) { organism <- Hsapiens }
if ((assembly == "mm9") || (assembly == "mm10")) { organism <- Mmusculus }
if (assembly == "sacCer3") { organism <- Scerevisiae}
if (assembly == "dm3") { organism <- Dmelanogaster}
if (assembly == "rn6") { organism <- Rnorvegicus}

txdb <- loadDb(txdbfile)
txdb

seqlevels(txdb,force=TRUE) <- seqlevels(txdb)[grep("_|\\d+.1$",seqlevels(txdb), invert=TRUE)]
seqlevels(txdb) <- sub("^","chr", seqlevels(txdb))
seqlevels(txdb) <- sub("MT","M", seqlevels(txdb))
seqlevels(txdb) <- sub("Mito","M", seqlevels(txdb))

print("Prepare gene models")
gnModel <- transcriptsBy(txdb, 'gene')
print("Retrieve seqinfo for organism")
seqinfo(gnModel) <- seqinfo(organism)[seqlevels(gnModel)]
model <- gnModel
tss <- promoters(model,upstream=up,downstream=down)

print ("Read in geneList")
geneListFile <- gsub('.txt$','',geneList)
geneList <- read.delim(file=geneList,header=FALSE)
head(tss)
head(names(tss))
                                        #geneTss <- tss[geneList[,"V1"] %in% names(tss),]
print("Matching genes in gene list to their annotation in txdbfile")
geneTss <- tss[geneList[,"V1"],]
print ("GeneTSS before unlist")
head(geneTss)
print ("GeneTSS after unlist")
geneTss <- unlist(geneTss)
geneTss$ensembl_gene_id <- names(geneTss)
head(geneTss)

### get the most highly occupied tss
Bin <- function(bw,model){
    cat("importing:", bw, sep="\n")
    bw.peak <- import.bw(bw,RangedData=FALSE,selection = BigWigSelection(model))
    cat("calc coverage\n")
    bw.peak.cov <- coverage(bw.peak,weight='score')
    cat("get coverage for peak region\n")
    mean.cov <- with(as.data.frame(model),{       
        mcmapply(function(seqname,start,end){
#            print(paste(seqname,start,end,tx_id,tx_name,sep=","))
            sum(bw.peak.cov[[seqname]][start:end])
        }
                ,mc.cores=8
                ,as.character(seqnames),start,end)
    })
    mean.cov <- data.frame(mean.cov)
    rownames(mean.cov) <- model$tx_name
    mean.cov
}

print ("Pick out the single transcript for each gene with the most signal in bwfile")
names(geneTss) <- NULL
tssCopy <- geneTss
print("Number of total transcripts")
ranges(geneTss)
peakCov <- do.call(cbind,mclapply(bwfile,model=tssCopy,Bin,mc.cores=1))
print("Pick most occupied TSS's and remove duplicates for top100.DUIV.up.tss")
## pick most occupied and remove duplicates
names(geneTss) <- geneTss$ensembl_gene_id
geneTss$tssTotalCovBWfile <- peakCov$mean.cov
or <- order(elementMetadata(geneTss)$tssTotalCovBWfile,decreasing=TRUE)
geneTss.or <- geneTss[or]
geneTss.or
reps <- duplicated(geneTss.or$ensembl_gene_id)
geneTss <- geneTss.or[!reps]
head(geneTss)
print("Number of transcripts with only one per gene")
ranges(geneTss)

bed.df <- data.frame(seqnames=seqnames(geneTss),
                     starts=start(geneTss)-1,
                     ends=end(geneTss),
                     names=paste(geneTss$ensembl_gene_id,geneTss$tx_name,sep="."),
                     scores=c(rep(".", length(geneTss))),
                     strands=strand(geneTss),
                     ids=geneTss$ensembl_gene_id)
if (keepOrder == 0){
    bed.df <- bed.df[order(bed.df$seqnames,bed.df$start),]
} else {
    geneList.clean <- geneList[setequal(geneList$V1,bed.df$ids),]
#    print(geneTss)
    print("geneList before clean")
    length(geneList)
    print(geneList)
    print("geneList after clean")
    print(geneList.clean)
    length(geneList.clean)
    print(length(bed.df$ids))
    print("bed before order")
    print(bed.df$ids)
    bed.df <- bed.df[order(match(bed.df$ids,geneList.clean)),]
    print("bed after order")
    print(bed.df$ids)
}

write.table(bed.df[,1:6], file=paste(geneListFile,up,down,"tss.bed",sep="."), quote=F, sep="\t", row.names=F, col.names=F)


