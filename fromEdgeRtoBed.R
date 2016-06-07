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
upBedFile <- paste(comparison,adjplabel,"up","bed",sep=".")
dnBedFile <- paste(comparison,adjplabel,"dn","bed",sep=".")

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

upBed <- cbind(as.data.frame(up$Chr), up$Start, up$End, as.data.frame(up$X), score=up$adj.p, up$Strand)
upBed$score <- as.integer(log(as.numeric(upBed$score),2)*-20)
highScore <- upBed$score > 1000
upBed$score[highScore]= 1000
upBed <- upBed[rev(order(upBed$score)),]

dn <- edger[edger$adj.p < adjp & edger$logFC < 0,]
dim(dn)

dnBed <- cbind(as.data.frame(dn$Chr), dn$Start, dn$End, as.data.frame(dn$X), score=dn$adj.p, dn$Strand)
dnBed$score <- as.integer(log(as.numeric(dnBed$score),2)*-20)
highScore <- dnBed$score > 1000
dnBed$score[highScore]= 1000
dnBed <- dnBed[rev(order(dnBed$score)),]

write.table(dnBed,file=dnBedFile,quote=F,sep="\t",row.names=F,col.names=F)
write.table(upBed,file=upBedFile,quote=F,sep="\t",row.names=F,col.names=F)
