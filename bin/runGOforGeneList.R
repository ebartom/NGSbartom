args <- commandArgs()

geneListFile <-sub('--glFile=', '', args[grep('--glFile=', args)])
assembly <-sub('--assembly=', '', args[grep('--assembly=', args)])

library(GO.db)
library(topGO)
library(ggplot2)
library(RColorBrewer)
library(biomaRt)

# I think the problem is that the geneList needs to be a vector of true and false indicating whether each gene is in the group or not.  The total length of the list needs to be the universe of possible genes.

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

#print("Loading biomart")
#bm <- useDataset(paste(ds), mart=bm)  
print("Retrieving GO information")
#listAttributes(bm)
EG2GO <- getBM(mart=bm, attributes=c('ensembl_gene_id','external_gene_id','go_id'))
head(EG2GO)

## Remove blank entries
print("Remove blank entries")
EG2GO <- EG2GO[EG2GO$go_id != '',]

print("pairing GO IDs with Ensembl IDs")
geneID2GO <- by(EG2GO$go_id,EG2GO$ensembl_gene_id,
                function(x) as.character(x))

xx <- as.list(GOTERM)
#head(xx)

##-----------------------------------Functions--------------------------------------#
##write a function to run topGO
runGO <- function(geneList=geneList,xx=xx,otype,setName){
    print("starting runGO")
    print(head(geneList))
    print(otype)
    print(setName)
                                        #    geneList        <- factor(as.integer(GOterms %in% list.terms))
                                        #    names(geneList) <- GOterms   
    setLength       <- sum(as.numeric(levels(geneList))[geneList])
    #setLength <- length(geneList)
    print(setLength)
    fname           <- paste(setName, setLength, otype, "GO.txt", sep="_")## file name of go table
    print(fname)
     GOData          <- new("topGOdata", ontology=otype, allGenes= geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)##run topGO
    resultFisher    <- runTest(GOData, algorithm = "classic", statistic = "fisher")## statistical test for topGO
    x               <- GenTable(GOData, classicFisher = resultFisher, topNodes=length(names(resultFisher@score)), orderBy="classicFisher")## make go table
    x               <- data.frame(x)
    pVal            <- data.frame(pval=signif(resultFisher@score, 6)) ## get unrounded pvalue
    pVal
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


print ("Total number of unique ensembl ids")
print(length(unique(EG2GO$ensembl_gene_id)))
all <- unique(EG2GO$ensembl_gene_id)
geneList <- factor(as.integer(all %in% geneList))
names(geneList) <- all
print(head(geneList))
geneList.setsize <- sum(as.numeric(levels(geneList))[geneList])
print(length(geneList))
print("setsize for significant genes") 
geneList.setsize

print("make GO table for the genes")

go.BP <- runGO(geneList=geneList,xx=xx,otype="BP",setName=basename(geneListFile))
go.MF <- runGO(geneList=geneList,xx=xx,otype="MF",setName=basename(geneListFile))

print("make the png for the up genes")
drawBarplot(go=go.BP,ontology="BP",setName=basename(geneListFile),setSize=geneList.setsize)
drawBarplot(go=go.MF,ontology="MF",setName=basename(geneListFile),setSize=geneList.setsize)
    
