args <- commandArgs()

geneListFile <-sub('--glFile=', '', args[grep('--glFile=', args)])
assembly <-sub('--assembly=', '', args[grep('--assembly=', args)])

library(biomaRt)

##----------load differentially expressed genes --------#
print("Loading gene list file")
print(geneListFile)

geneList <- read.delim(file=geneListFile,header=TRUE,sep="\t")
geneList <- unique(geneList[,1])
length(geneList)

##---------load correct Biomart------------------------#

hostMart <- ""

if (assembly == "hg19") {
   organismStr <- "hsapiens"
   hostMart <- "feb2014.archive.ensembl.org"
}
if (assembly == "hg38") {
   organismStr <- "hsapiens"
   hostMart <- "feb2014.archive.ensembl.org"
}
if (assembly == "mm9") {
   organismStr <- "mmusculus"
   hostMart <- "feb2014.archive.ensembl.org"
}
if (assembly == "mm10") {
   organismStr <- "mmusculus"
   hostMart <- "feb2014.archive.ensembl.org"
}
if (assembly == "sacCer3") {
   organismStr <- "scerevisiae"
   hostMart <- "feb2014.archive.ensembl.org"
}
if (assembly == "dm3") {
   organismStr <- "dmelanogaster"
   hostMart <- "feb2014.archive.ensembl.org"
}
organismStr

print("setting up ensembl")
dataset = paste(organismStr,"_gene_ensembl",sep="")
bm <- useMart("ENSEMBL_MART_ENSEMBL",host=hostMart,dataset=dataset)

ens2ext <- getBM(mart=bm, attributes=c('ensembl_gene_id','external_gene_id'))
head(ens2ext)

## Remove blank entries
print("Remove blank entries")
ens2ext <- ens2ext[ens2ext$external_gene_id!= '',]

print("pairing external IDs with Ensembl IDs")
geneIDconv <- by(ens2ext$ensembl_gene_id,ens2ext$external_gene_id,
                function(x) as.character(x))
print("geneIDconv List")
print(head(geneIDconv))
print("geneList")
print(head(geneList))

ensList <- geneIDconv[geneList]
print("ensList")
print(head(ensList))
                                        #write.table(ensList,file=paste(geneListFile,"ens.txt",sep="."),quote=FALSE,sep="\n")
filename <- gsub(".txt$","",geneListFile)
filename <- paste(filename,"ens.txt",sep=".")
cat(sapply(ensList, toString),file=filename,sep="\n")
