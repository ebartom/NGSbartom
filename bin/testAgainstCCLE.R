library(ggplot2)
library(psych)
#library(ggcorrplot)
#install.packages("corrplot")
library(corrplot)

args <- commandArgs()

query.path <-sub('--queryDir=', '', args[grep('--queryDir=', args)])
querystring <-basename(query.path)

#ccle.path<-"/projects/b1042/BartomLab/CCLE/CCLE.NCM/"
ccle.path<-"/projects/p20742/CCLE.RNA/CCLE.NCM"
#query.path<-"/projects/p20742/celID.testData/lesniakQueries"
#query.path<-"/projects/p20742/celID.testData/Roger.NCM"
print(query.path)
print(ccle.path)
setwd(ccle.path)

#Read in the list of all ncm files at the above path.
ccle.ncm <- dir(pattern = "*.ncm$")

# Define a list to hold the data.
ldf<-list()

# First read in all of the NCM files
for (k in 1:length(ccle.ncm)){
    ldf[[k]] <- read.table(ccle.ncm[k],header=TRUE,row.names = 1)
    }

# Create a dataframe to hold all of the VAF values for CCLE
ncm <- data.frame(matrix(ncol = 1, nrow = length(ldf[[1]]$vaf)))
rownames(ncm) <- rownames(ldf[[1]])
dim(ncm)

# Then pull in the rest of vaf data.
for (k in 1:length(ccle.ncm)){
   ncm[ccle.ncm[k]]<-ldf[[k]]$vaf
}

# Remove the first empty row.
ncm<-ncm[-1]
print(dim(ncm))

setwd(query.path)

#Read in the list of all ncm files at the query path.
query.ncm <- dir(pattern = "*.ncm$")

# Define a list to hold the data.
qdf<-list()

# First read in all of the NCM files
for (k in 1:length(query.ncm)){
    qdf[[k]] <- read.table(query.ncm[k],header=TRUE,row.names = 1)
    }

# Create a dataframe to hold all of the VAF values for the query samples
queryVAF <- data.frame(matrix(ncol = 1, nrow = length(qdf[[1]]$vaf)))
rownames(queryVAF) <- rownames(qdf[[1]])
dim(queryVAF)

# Then pull in the rest of vaf data.
for (k in 1:length(query.ncm)){
   queryVAF[query.ncm[k]]<-qdf[[k]]$vaf
}

# Remove the first empty row.
queryVAF<-queryVAF[-1]

# Query VAF now has the ncm files for all of the query samples.
miniresults<-data.frame(matrix(ncol=length(query.ncm),nrow=length(query.ncm)))
rownames(miniresults)<-sub(".ncm","",query.ncm)
colnames(miniresults)<-sub(".ncm","",query.ncm)
for (i in 1:length(query.ncm)){
    for (j in 1:length(query.ncm)){
        miniresults[i,j]<-(cor(queryVAF[,i],queryVAF[,j],method="pearson",use="pairwise.complete.obs"))
	}
}
summary(miniresults)

# Within query correlation.
print(miniresults)
# Frustratingly, I cannot seem to alter the margins on this plot.
#par(mar=c(3,3,1,1))	
filename<-paste(c(querystring,"withinStudyCorrelations.pdf"),collapse=".")
pdf(filename,width=12,height=12)
corPlot(miniresults,diag=TRUE,las=2,main="Correlation Between Sample Genotypes",cex=1.25,
        cex.axis=1)
dev.off()

# Set up a table with the results of the correlations between queries and CCLE and populate it.
results <- data.frame(matrix(ncol = length(query.ncm), nrow = length(ccle.ncm)))
rownames(results)<-ccle.ncm
colnames(results)<-query.ncm

for (i in 1:length(ccle.ncm)){
    for (j in 1:length(query.ncm)){
	results[i,j]<-(cor(ncm[,i],queryVAF[,j],method="pearson",use="pairwise.complete.obs"))
	if (results[i,j] <0) { results[i,j]<- 0}
	}	
}

# It turns out we made a list; let's make a dataframe.
df <- data.frame(matrix(unlist(results),ncol=length(query.ncm)))
colnames(df) <-query.ncm
rownames(df) <-ccle.ncm

# Print a summary of the statistics for each sample (min, max correlations, etc).
summary(df)

for (j in 1:length(query.ncm)){
    plotHist <- paste(c(query.ncm[j],"pdf"),collapse=".")
    print(plotHist)
    maxMatch<-max(df[,j])
   for (i in 1:length(ccle.ncm)){
    	if (df[i,j] == maxMatch){
	   cellLine <- ccle.ncm[i]
	}
    }
    print(maxMatch)
    print(cellLine)
    sample<-query.ncm[j]
    sampleLabel<-sub(".ncm","",sample)
    print(sample)
    h1<-hist(df[,sample],
    breaks=c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1))
    h1counts<-h1$counts
        h1counts[h1counts == 0] = ""	
    #print(h1counts)
    pdf(plotHist) 
    hist(df[,sample],
	breaks=c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1),
	xlim=c(0,1),labels=h1counts,
	xlab="Pearson Correlation Between Sample and Database",main=sampleLabel,col="#4E2A84")
#    abline(v=0.8,col="red",lty=2)
     if (maxMatch > 0.5){
     	cellLabel <- sub(".ncm","",cellLine)
	cellLabel <- gsub(".SRR\\d+","",cellLabel)
	print(cellLabel)
     	text(x=0.4, y=200, adj=0, labels=expression(bold("Top Hit:")),pos=3,offset=1)
	text(x=0.4, y=200, adj=0, labels=cellLabel)
	text(x=0.4, y=200, adj=0, labels=sprintf("%.2f", maxMatch),pos=1,offset=1)
	}
    dev.off()
}

write.table(results,"allCor.txt")
