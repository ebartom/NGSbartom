args <- commandArgs()

edgeRfile <-sub('--edgeRfile=', '', args[grep('--edgeRfile=', args)])
adjp <-as.numeric(sub('--adjp=', '', args[grep('--adjp=', args)]))

method <- gsub(".edgeR.txt","",edgeRfile)
method <- gsub("^.*\\/.*\\.","",method)
print(method)
adjplabel <- gsub("^0\\.","",adjp)
print(adjplabel)
comparison <- gsub("\\.edgeR.txt$|\\.edgeR.rda","",edgeRfile)
print(comparison)
upListFile <- paste(comparison,adjplabel,"up","geneList.txt",sep=".")
dnListFile <- paste(comparison,adjplabel,"dn","geneList.txt",sep=".")

print("Loading edgeR data")
print(edgeRfile)

if(grepl('rda',edgeRfile)){
   edger <- get(load(file=edgeRfile))
}
if(grepl('txt',edgeRfile)){
    edger <- read.delim(file=edgeRfile,header=TRUE,sep="\t")
    rownames(edger) <- edger$X
}
#edger <- data.frame(edger)
head(edger)
dim(edger)


up <- edger[edger$adj.p < adjp & edger$logFC > 0,]
dim(up)

upList <- as.data.frame(up$X)

dn <- edger[edger$adj.p < adjp & edger$logFC < 0,]
dim(dn)

dnList <- as.data.frame(dn$X)

write.table(dnList,file=dnListFile,quote=F,sep="\t",row.names=F,col.names=F)
write.table(upList,file=upListFile,quote=F,sep="\t",row.names=F,col.names=F)
