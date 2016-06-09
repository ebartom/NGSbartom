args <- commandArgs()

degFile <-sub('--degFile=', '', args[grep('--degFile=', args)])
adjp <-as.numeric(sub('--adjp=', '', args[grep('--adjp=', args)]))
assembly <-sub('--assembly=', '', args[grep('--assembly=', args)])

if (identical(adjp,character(0))){
   adjp<-0.01
}

library(GO.db)
library(topGO)
library(ggplot2)
library(RColorBrewer)
library(biomaRt)

##----------load differentially expressed genes --------#
print("Loading differential expressed gene table")
print(degFile)
method <- gsub(".edgeR.txt","",degFile)
method <- gsub("^.*\\/.*\\.","",method)
print(method)

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
print("Retrieving GO information")
EG2GO <- getBM(mart=bm, attributes=c('ensembl_gene_id','external_gene_id','go_id','ensembl_transcript_id'))
head(EG2GO)

#listMarts(host=hostMart)
#bm <- useMart("ensembl")
#ds <- listDatasets(bm)
#ds <- ds[grep(organismStr,ds$dataset),]$dataset
#print("Loading biomart")
#bm <- useDataset(paste(ds), mart=bm)  
#EG2GO <- getBM(mart=bm, attributes=c('ensembl_gene_id','external_gene_name','go_id'))
#head(EG2GO)

## Remove blank entries
print("Remove blank entries")
EG2GO <- EG2GO[EG2GO$go_id != '',]

print("pairing GO IDs with Ensembl IDs")
if (identical(method,"rsem")){
    geneID2GO <- by(EG2GO$go_id,EG2GO$ensembl_transcript_id,
                    function(x) as.character(x))
}else{
    geneID2GO <- by(EG2GO$go_id,EG2GO$ensembl_gene_id,
                    function(x) as.character(x))
}

xx <- as.list(GOTERM)

##-----------------------------------Functions--------------------------------------#
##write a function to run topGO
runGO <- function(geneList,xx=xx,otype,setName){
#    geneList        <- factor(as.integer(GOterms %in% list.terms))
#    names(geneList) <- GOterms   
    setLength       <- sum(as.numeric(levels(geneList))[geneList]) 
    fname           <- paste(setName, setLength, otype, "GO.txt", sep="_")## file name of go table
    GOData          <- new("topGOdata", ontology=otype, allGenes=geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)##run topGO
    resultFisher    <- runTest(GOData, algorithm = "classic", statistic = "fisher")## statistical test for topGO
    x               <- GenTable(GOData, classicFisher = resultFisher, topNodes=length(names(resultFisher@score)), orderBy="classicFisher")## make go table
    x               <- data.frame(x)
    pVal            <- data.frame(pval=signif(resultFisher@score, 6)) ## get unrounded pvalue
    x$enrich        <- x$Significant/x$Expected ## calculate enrichment based on what you expect by chance
    x$classicFisher <- as.numeric(x$classicFisher)## put pvalue in the table
    x$p.unround     <- pVal[x$GO.ID,"pval"]## put unrounded pvalue in the table
    x$p.adj         <- signif(p.adjust(x$classicFisher, method="BH"), 6)## calculate the adjusted pvalue with Benjamini & Hochberg correction
    x$log.p     <- -log10(x$p.unround) ## convert adjusted p value to -log10 for plot magnitude
    x               <- x[x$p.unround<=0.05,] ## calculate the adjusted pvalue with Benjamini & Hochberg correction aka "fdr"
    x$Term.full     <- sapply(x$GO.ID, FUN=function(n){Term(xx[[n]])}) ## get the full term name
    write.table(x, file=fname, sep="\t", col.names=NA) ## save the table
    printGraph(GOData,## make the tree for the go data
               resultFisher,
               firstSigNodes = 5,
               fn.prefix = paste(setName, setLength, otype, sep="_"),
               useInfo = "all",
               pdfSW = TRUE
               )    
    return(x)  
}

## function to make barplot of -log10 adjusted pvalues colored by enrichment
drawBarplot <- function(go, ontology, setName, setSize){
    #print(go)
#    print(ontology)
    print(setName)
    print(setSize)
    go$Term.full <-factor(go$Term.full, levels=go[order(go$log.p), "Term.full"]) ## sort table by adjusted p-value
    ptitle <- paste(ontology, setName, setSize) ## plot title
    ptitle <- gsub("^.*/","",ptitle)
    pfname <- paste(setName,ontology,"png",sep=".")## name of png file
    if(nrow(go) < 20 ){
        toprange <- 1:nrow(go)
    }else{
        toprange <- 1:20
    }
    png(filename=paste(pfname),height=600,width=700)
    print({
    p <- ggplot(go[toprange,], aes(y=log.p, x=Term.full)) + ## ggplot barplot function
      geom_bar(stat="identity",fill="green") +
      ggtitle(ptitle) +
      xlab("") + ylab("-log10(p-value)") +
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
all <-unique(rownames(deg))
up.geneList <-  factor(as.integer(all %in% up))
names(up.geneList) <- all
print(head(up.geneList))

up.setsize <- sum(as.numeric(levels(up.geneList))[up.geneList])
print("setsize for significant genes") 
up.setsize

adjplabel <- gsub("^0\\.","",adjp)
comparison <- gsub("\\.edgeR.txt$|\\.edgeR.df.rda$","",degFile)

print("make GO table for the up genes")

go.UP.BP <- runGO(geneList=up.geneList,xx=xx,otype="BP",setName=paste(comparison,adjplabel,"up",sep="."))
go.UP.MF <- runGO(geneList=up.geneList,xx=xx,otype="MF",setName=paste(comparison,adjplabel,"up",sep="."))

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
go.DN.BP <- runGO(geneList=dn.geneList,xx=xx,otype="BP",setName=paste(comparison,adjplabel,"down",sep="."))
go.DN.MF <- runGO(geneList=dn.geneList,xx=xx,otype="MF",setName=paste(comparison,adjplabel,"down",sep="."))

print("make barplot for down genes")
drawBarplot(go=go.DN.BP,ontology="BP",setName=paste(comparison,adjplabel,"down",sep="."),setSize=dn.setsize)
drawBarplot(go=go.DN.MF,ontology="MF",setName=paste(comparison,adjplabel,"down",sep="."),setSize=dn.setsize)
