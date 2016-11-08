library(limma)
library(VennDiagram)

args <- commandArgs()

geneListA <-sub('--geneListA=', '', args[grep('--geneListA=', args)])
geneListB <-sub('--geneListB=', '', args[grep('--geneListB=', args)])
geneListC <-sub('--geneListC=', '', args[grep('--geneListC=', args)])

print ("Read in geneList A")
geneListFileA <- gsub('.geneList.txt$','',geneListA)
geneListA <- read.delim(file=geneListA,header=FALSE)
print ("geneListA")
#geneListA

print ("Read in geneList B")
geneListFileB <- gsub('.geneList.txt$','',geneListB)
geneListB <- read.delim(file=geneListB,header=FALSE)
print ("geneListB")
#geneListB

if (identical(geneListC,character(0))){
    geneListFileC <- "NA"
} else {
    print ("Read in geneList C")
    geneListFileC <- gsub('.geneList.txt$','',geneListC)
    geneListC <- read.delim(file=geneListC,header=FALSE)
    print("geneListC")
}
ab <- length(unique(intersect(geneListA[,1],geneListB[,1])))
a <- length(unique(geneListA[,1]))
a
b <- length(unique(geneListB[,1]))
b
ab
geneListFileA
geneListFileB
pdfname <- paste(geneListFileA,geneListFileB,"vennDiagram.pdf",sep=".")

geneListFileC
if (identical(geneListFileC, "NA")){
    pdfname
    pdf(pdfname)
    draw.pairwise.venn(a,b,cross.area = ab,
                       category = c(geneListFileA,geneListFileB),
                       fill = c("red","yellow"),
                       cex = 2,
                       cat.cex = 2
                       #                       lty = "blank",
   #                    cat.pos=c(-40,20),
    #                   alpha=c(0.8,0.8)
                       )
    dev.off()
} else {
    c <- length(unique(geneListC[,1]))
    c
    ac <- length(unique(intersect(geneListA[,1],geneListC[,1])))
    ac
    bc <- length(unique(intersect(geneListB[,1],geneListC[,1])))
    bc
    abc <- length(unique(intersect(intersect(geneListA[,1],geneListB[,1]),geneListC[,1])))
    abc
    pdfname <- paste(geneListFileA,geneListFileB,geneListFileC,"vennDiagram.pdf",sep=".")
    pdfname
    pdf(pdfname)
    draw.triple.venn(a,b,c,ab,bc,ac,abc,
                     category = c(geneListFileA,geneListFileB,geneListFileC),
                     fill = c("red","yellow","blue"),
                     cex = 2,
                     cat.cex = 2
                     )
    dev.off()
}
