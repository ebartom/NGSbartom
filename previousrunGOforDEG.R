args <- commandArgs()

degFile <-sub('--degFile=', '', args[grep('--degFile=', args)])
adjp <-sub('--adjp=', '', args[grep('--adjp=', args)])
assembly <-sub('--assembly=', '', args[grep('--assembly=', args)])

if (identical(adjp,character(0))){
   adjp<-0.01
}

library(topGO)
library(ggplot2)
library(RColorBrewer)
library(biomaRt)

##----------load differentially expressed genes --------#
print("Loading differential expressed gene table")
print(degFile)

if(grepl('rda',degFile)){
   deg <- get(load(file=degFile))
}
if(grepl('txt',degFile)){
    deg <- read.delim(file=degFile,header=TRUE,sep="\t")
    rownames(deg)<-deg$X
}

adjplabel <- gsub("^0\\.","",adjp)
comparison <- gsub("\\.edgeR.txt$|\\.edgeR.rda","",degFile)
all <-unique(rownames(deg))

##---------load correct Biomart------------------------#

hostMart <- ""

if (assembly == "hg19") {
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

listMarts(host=hostMart)

bm <- useMart("ensembl")
ds <- listDatasets(bm)
ds <- ds[grep(organismStr,ds$dataset),]$dataset

print("Loading biomart")
bm <- useDataset(paste(ds), mart=bm)  

EG2GO <- getBM(mart=bm, attributes=c('ensembl_gene_id','external_gene_name','go_id'))
head(EG2GO)

## Remove blank entries
print("Remove blank entries")
EG2GO <- EG2GO[EG2GO$go_id != '',]

print("pairing GO IDs with Ensembl IDs")
geneID2GO <- by(EG2GO$go_id,EG2GO$ensembl_gene_id,
                function(x) as.character(x))

##-----------------------------------Functions--------------------------------------#
##write a function to run topGO
runGO <- function(otype,setName,setLength,geneList){
    print(otype)
    print(setName)
    print(setLength)
    fname <- paste(setName, setLength, otype, "GO.txt", sep="_") ## file name of go table
    GOData <- new("topGOdata", ontology=otype, allGenes=geneList, nodeSize=20,## run topGO
                    annot=annFUN.gene2GO,
                    gene2GO=geneID2GO)
    resultFisher <- runTest(GOData, algorithm = "classic", statistic = "fisher")## statistical test for topGO
    x <- GenTable(GOData, classicFisher = resultFisher, topNodes=40, orderBy="classicFisher") ## make go table
    x$enrich <- x$Significant/x$Expected ## calculate enrichment based on what you expect by chance
    x$classicFisher <- as.numeric(x$classicFisher) ## put pvalue in the table
    x$p.adj <- p.adjust(x$classicFisher,method="BH") ## calculate the adjusted pvalue with Benjamini & Hochberg correction aka "fdr"
    x$log.p.adj <- -log10(x$p.adj) ## convert adjusted pvalue to -log10 scale
    x <- x[!x$p.adj > 0.05,]   ## throw out unsignificant pvalues
    write.table(x, file=fname, sep="\t", col.names=NA)## save the go table
    printGraph(GOData,## make the tree for the go data
               resultFisher,
               firstSigNodes = 5,
               fn.prefix = paste(setName, setLength, otype, sep="_"),
               useInfo = "up",
               pdfSW = TRUE
               )
    return(x)  
}

## function to make barplot of -log10 adjusted pvalues colored by enrichment
drawBarplot <- function(go, ontology, setName, setSize){
    print(go)
#    print(ontology)
    print(setName)
    print(setSize)
    go$Term <-factor(go$Term, levels=go[order(go$log.p.adj), "Term"]) ## sort table by adjusted p-value
    ptitle <- paste(ontology, setName, setSize) ## plot title
    ptitle <- gsub("^.*/","",ptitle)
    pfname <- paste(setName,ontology,"png",sep=".")## name of png file
    png(filename=paste(pfname),height=600,width=700)
    print({
    p <- ggplot(go, aes(y=log.p.adj, x=Term, fill=enrich)) + ## ggplot barplot function
      geom_bar(stat="identity") +
      ggtitle(ptitle) +
      xlab("") + ylab("-log10(adjusted p-value)") +
      scale_fill_gradientn(colours =my_palette, limits = c(0, max(go$enrich))) +
      coord_flip()+
      theme(panel.grid.major = element_line(colour = "grey"), panel.grid.minor = element_blank(), 
                panel.background = element_blank(), axis.line = element_line(colour = "black"))+
      theme(text = element_text(size=16),
          axis.text.x = element_text(vjust=1,color="black",size=14),
          axis.text.y = element_text(color="black",size=14),
            plot.title=element_text(size=16))    
})
    dev.off()
}
   

print("get up genes and make geneList")
up <- deg$adj.p < adjp & deg$logFC > 0
up <- unique(rownames(deg[up,]))
all <-unique(deg$X)
up.geneList <-  factor(as.integer(all %in% up))
names(up.geneList) <- all
#head(up.geneList)

up.setsize <- sum(as.numeric(levels(up.geneList))[up.geneList])
print("setsize for significant genes") 
up.setsize

adjplabel <- gsub("^0\\.","",adjp)
comparison <- gsub("\\.edgeR.txt$","",degFile)

print("make GO table for the up genes") 
go.UP.BP <- runGO("BP",paste(comparison,adjplabel,"up",sep="."),up.setsize,up.geneList)
go.UP.MF <- runGO("MF",paste(comparison,adjplabel,"up",sep="."),up.setsize,up.geneList)

my_palette <- colorRampPalette(c("blue","#FF5000"))(n = 10)

print("make the png for the up genes")
drawBarplot(go=go.UP.BP,ontology="BP",setName=paste(comparison,adjplabel,"up",sep="."),setSize=up.setsize)
drawBarplot(go=go.UP.MF,ontology="MF",setName=paste(comparison,adjplabel,"up",sep="."),setSize=up.setsize)
    
print("get down genes and make geneList")
dn <- deg$adj.p < adjp & deg$logFC < 0
dn <- unique(rownames(deg[dn,]))
all <-unique(rownames(deg))
dn.geneList <-  factor(as.integer(all %in% dn))
names(dn.geneList) <- all

dn.setsize <- sum(as.numeric(levels(dn.geneList))[dn.geneList])
print("setsize for significant genes") 
dn.setsize

print("make GO table for down genes")
go.DN.BP <- runGO("BP",paste(comparison,adjplabel,"down",sep="."),dn.setsize,dn.geneList)
go.DN.MF <- runGO("MF",paste(comparison,adjplabel,"down",sep="."),dn.setsize,dn.geneList)

print("make barplot for down genes")
drawBarplot(go=go.DN.BP,ontology="BP",setName=paste(comparison,adjplabel,"down",sep="."),setSize=dn.setsize)
drawBarplot(go=go.DN.MF,ontology="MF",setName=paste(comparison,adjplabel,"down",sep="."),setSize=dn.setsize)
