library(Rsamtools)
library(GenomicAlignments)
library(parallel)
library(GenomicFeatures)
library(rtracklayer)
library(biomaRt)
library(GenomicRanges)

args <- commandArgs()

assembly <-sub('--assembly=', '', args[grep('--assembly=', args)])
peakFile <-sub('--peakFile=', '', args[grep('--peakFile=', args)])
outDir <-sub('--outputDirectory=','',args[grep('--outputDirectory=',args)])
txdbfile <-sub('--txdbfile=','',args[grep('--txdbfile=',args)])

print("assembly")
assembly
print("peakFile")
peakFile
print("outDir")
outDir
print("txdbfile")
txdbfile
#peakFile<-"/projects/b1025/etb/andreaPRC/DUIVchip/peaks/DUIV-K27ac-071415.macsPeaks.bed"
#outDir<-"/projects/b1025/etb/andreaPRC/DUIVchip/peaks/"
#assembly<-"hg19"
#txdbfile<- "/projects/b1025/anno/Txdb/hsapiens_gene_ensembl_Ens75.txdb"

print(assembly)
if (assembly == "hg19") { organismStr <- "Hsapiens" }
if ((assembly == "mm9") || (assembly == "mm10")) { organismStr <- "Mmusculus" }
if (assembly == "sacCer3") { organismStr <- "Scerevisiae"}
if (assembly == "dm3") { organismStr <- "Dmelanogaster"}

assemblyLibrary <- paste("BSgenome.", organismStr, ".UCSC.", assembly, sep="")
print(assemblyLibrary)
library(assemblyLibrary,character.only=TRUE)

if (assembly == "hg19") { organism <- Hsapiens }
if ((assembly == "mm9") || (assembly == "mm10")) { organism <- Mmusculus }
if (assembly == "sacCer3") { organism <- Scerevisiae}
if (assembly == "dm3") { organism <- Dmelanogaster}

txdb <- loadDb(txdbfile)
txdb

seqlevels(txdb,force=TRUE) <- seqlevels(txdb)[grep("_|\\d+.1$",seqlevels(txdb), invert=TRUE)]
seqlevels(txdb) <- sub("^","chr", seqlevels(txdb))
seqlevels(txdb) <- sub("MT","M", seqlevels(txdb))
seqlevels(txdb) <- sub("Mito","M", seqlevels(txdb))

gnModel <- transcriptsBy(txdb, 'gene')
gnModel <- unlist(gnModel)
seqinfo(gnModel) <- seqinfo(organism)[seqlevels(gnModel)]

counts <- read.delim(file=peakFile,header=FALSE,sep="\t")
colnames(counts)<-c("chr","start","end","name","score")
counts$name <- sub(" ","",counts$name)
counts$chr <- sub(" ","",counts$chr)
gcounts <- as(counts,"GRanges")
genome(gcounts) <- assembly
seqlengths(gcounts) <- seqlengths(gnModel)[seqlevels(gcounts)]
head(gcounts)

counts$nearestGene<-names(gnModel[nearest(gcounts,gnModel)])
dist<-as.data.frame(distanceToNearest(gcounts,gnModel))
counts$distToNearestGene<-dist$distance
tss <-promoters(gnModel,upstream=0,downstream=1)
tss <-resize(gnModel,fix="start",width=1)
counts$nearestTSS<-names(tss[nearest(gcounts,tss)])
dist<-as.data.frame(distanceToNearest(gcounts,tss))
counts$distToNearestTSS<-dist$distance
head(counts)

dataset = paste(tolower(organismStr),"gene_ensembl",sep="_")
dataset
hostMart <- "ensembl.org"
#hostMart <- "dec2014.archive.ensembl.org"
bm <- useMart("ENSEMBL_MART_ENSEMBL",host=hostMart,dataset=dataset)
anno <- getBM(mart=bm, attributes=c('ensembl_gene_id','gene_biotype','external_gene_name'))
#anno <- getBM(mart=bm, attributes=c('ensembl_gene_id','gene_biotype','external_gene_name','description'))

counts.anno <- merge(counts,anno,by.x="nearestTSS",by.y="ensembl_gene_id")
names(counts.anno)[names(counts.anno)=="external_gene_name"]<- "tssGeneName"
names(counts.anno)[names(counts.anno)=="gene_biotype"]<- "tssGeneBiotype"
#names(counts.anno)[names(counts.anno)=="description"]<- "tssGeneDescription"
counts.anno <- merge(counts.anno,anno,by.x="nearestGene",by.y="ensembl_gene_id")
names(counts.anno)[names(counts.anno)=="external_gene_name"]<- "nearestGeneName"
names(counts.anno)[names(counts.anno)=="gene_biotype"]<- "nearestGeneBiotype"
counts.anno<-counts.anno[,c("chr","start","end","name","score","distToNearestGene","nearestGene","nearestGeneName","nearestGeneBiotype","distToNearestTSS","nearestTSS","tssGeneName","tssGeneBiotype")]
head(counts.anno)

peakPrefix <- sub(".bed$","",peakFile)
write.table(counts.anno, file=paste(peakPrefix,"anno.txt",sep="."),sep="\t",row.names=FALSE)
