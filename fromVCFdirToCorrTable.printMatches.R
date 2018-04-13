
library("VariantAnnotation")
library(gplots)

library(RColorBrewer)
monochrome <- brewer.pal(9, "Greys")

args<-commandArgs()
vcfDir <- sub('--vcfDir=','',args[grep('--vcfDir=',args)])
assembly <- sub('--assembly=','',args[grep('--assembly=',args)])
filePattern <- sub('--pattern=','',args[grep('--pattern=',args)])
vcfDir
assembly
#filePattern<-paste(filePattern,"$",sep="")
filePattern


fileSet <- list.files(path = vcfDir, pattern = filePattern)
fileSet
length(fileSet)
emptyCorrelations <- rep(0,length(fileSet)*length(fileSet))
corrTable<- matrix(emptyCorrelations,nrow=length(fileSet),ncol=length(fileSet))
corrTable[5,5]
i<-0
for (file1 in fileSet){
    i<-i+1
    i
    print(paste("VCF1 is", file1))
    vcfFile1 <- paste(vcfDir,file1,sep="/")
    vcf1 <- readVcf(file = vcfFile1, genome = assembly)
    alt1 = alt(vcf1)
    ref1 = ref(vcf1)
    a1=rowRanges(vcf1)
#    	print(head(a1))
    snam1= paste(as.character(a1@seqnames),a1@ranges@start,as.character(ref1),as.character(unlist(alt1)),sep='_')
 #   	print(head(snam1))
    x1 = unlist(geno(vcf1)$AF)
    names(x1)<-snam1
 #  	print(head(x1))
 #   	dim(x1)
    j<-0
    for (file2 in fileSet){
      	j<- j+1
	j
    	print(paste("VCF2 is", file2))
	vcfFile2 <- paste(vcfDir,file2,sep="/")	
    	vcf2 <- readVcf(file = vcfFile2, genome = assembly)
    	alt2 = alt(vcf2)
    	ref2 = ref(vcf2)
    	a2=rowRanges(vcf2)
    	#print(head(a2))
    	snam2= paste(as.character(a2@seqnames),a2@ranges@start,as.character(ref2),as.character(unlist(alt2)),sep='_')
 #   	print(head(snam2))
    	x2 = unlist(geno(vcf2)$AF)
    	names(x2)<-snam2
#    	print(head(x2))
#    	dim(x2)

    	matchingSNPs1<-snam1[snam1 %in% snam2]
    	matchingSNPs2<-snam2[snam2 %in% snam1]
    	print("Num matchingSNPs")
	print(length(matchingSNPs1))
	print(i)
#   	print(head(matchingSNPs1))
#    	print("matchingSNPs2")
	print(j)
#    	print(head(matchingSNPs2))
#    	correlation = cor(x1[matchingSNPs1],x2[matchingSNPs2],use="pairwise.complete.obs")
	#print(head(matchingSNPs1))
	#print(head(x1[matchingSNPs1],n=10))
	#print(head(matchingSNPs2))
	#print(head(x2[matchingSNPs2],n=10))
	pairedGenotypes<-cbind(x1[matchingSNPs1],x2[matchingSNPs2])
	#write.table(pairedGenotypes,file=paste("match",i,j,"txt",sep="."))
	correlation = cor(x1[matchingSNPs1],x2[matchingSNPs2])
    	print("Pairwise correlation")
	file1
	label1 <-gsub("\\.raw.*vcf$","",file1)
	file2
	label2 <-gsub("\\.raw.*vcf$","",file2)
	colnames(pairedGenotypes) <- c(label1,label2)
	write.table(pairedGenotypes,file=paste("match",label1,label2,correlation,"txt",sep="."))
	if (label1 != label2){
		pdf(file=paste("match",label1,label2,correlation,"pdf",sep="."),width=3,height=8)
		heatmap.2(pairedGenotypes,scale="none",sepwidth=0,Rowv=NULL,Colv=NULL,trace="none",col=monochrome,xlab=NULL)
		dev.off()
		}
	print(correlation)
	corrTable[i,j]<-correlation
	}
}
rownames(corrTable)<-fileSet
colnames(corrTable)<-fileSet

label<-gsub("\\*","",filePattern)
filePrefix<-paste("genoCorr",label,sep=".")
filePrefix <- gsub ("\\.\\.","\\.",filePrefix)
pdfFile<-paste(filePrefix,"pdf",sep=".")
tableFile<-paste(filePrefix,"txt",sep=".")
pdfFile
tableFile

head(corrTable)
write.table(corrTable,file=tableFile)
pdf(pdfFile)
heatmap.2(corrTable,trace="none")
dev.off()
