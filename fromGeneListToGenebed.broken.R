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
geneList <-sub('--geneList=', '', args[grep('--geneList=', args)])

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
seqlevels(txdb,force=TRUE) <- seqlevels(txdb)[grep("^GL",seqlevels(txdb), invert=TRUE)]
seqlevels(txdb) <- sub("^","chr", seqlevels(txdb))
seqlevels(txdb) <- sub("MT","M", seqlevels(txdb))
seqlevels(txdb) <- sub("Mito","M", seqlevels(txdb))

print("Prepare gene models")
gnModel <- transcriptsBy(txdb, 'gene')
print("gnModel")
head(gnModel)
print("Retrieve seqinfo for organism")
seqinfo(gnModel) <- seqinfo(organism)[seqlevels(gnModel)]
model <- unlist(gnModel)
model$length <- width(model)
model$ensembl_gene_id<- names(model)
print("model")
head(model)

or <- order(elementMetadata(model)$length,decreasing=TRUE)
model.or <- model[or]
model.or
reps <- duplicated(model.or$ensembl_gene_id)
model <- model.or[!reps]
head(model)
print("Number of Genes, after removing all but the longest transcript per gene")
ranges(model)

print ("Read in geneList")
geneListFile <- gsub('.txt$','',geneList)
geneList <- as.vector(read.delim(file=geneList,header=FALSE))
print("Matching genes in gene list to their annotation in txdbfile")
myKnownGenes <- as.vector(names(model))
myFoundGenes <- geneList %in% myKnownGenes
print("my list")
head(geneList)
print(typeof(geneList))
print("my known genes")
head(myKnownGenes)
print(typeof(myKnownGenes))
print("intersection")
head(myFoundGenes)
myGenes <- model[myFoundGenes,]
print ("myGenes before unlist")
head(myGenes)
print ("myGenes after unlist")
myGenes <- unlist(myGenes)
myGenes$ensembl_gene_id <- names(myGenes)
head(myGenes)

bed.df <- data.frame(seqnames=seqnames(myGenes),
                     starts=start(myGenes)-1,
                     ends=end(myGenes),
                     names=paste(myGenes$ensembl_gene_id,myGenes$tx_name,sep="."),
                     scores=c(rep(".", length(myGenes))),
                     strands=strand(myGenes),
                     ids=myGenes$ensembl_gene_id)
geneList.clean <- geneList[setequal(geneList$V1,bed.df$ids),]
#    print(myGenes)
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


write.table(bed.df[,1:6], file=paste(geneListFile,"geneBodies.bed",sep="."), quote=F, sep="\t", row.names=F, col.names=F)


