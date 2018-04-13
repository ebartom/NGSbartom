args <- commandArgs()

clinicalData <-sub('--clinicalData=', '', args[grep('--clinicalData=', args)])
probeList<-sub('--probeList=','',args[grep('--probeList=',args)])
genomicMatrix<-sub('--genomicMatrix=','',args[grep('--genomicMatrix=',args)])

library(survival)
library(gplots)

##----------load differentially expressed genes --------#
#print("Loading gene expression matrix")
#print(genomicMatrix)

allExpression <- read.table(file=genomicMatrix,header=TRUE,sep="\t",row.names=1)
#head(allExpression)
pl <- read.delim(file=probeList,header=FALSE,sep="\t")
clinical <- read.delim(file=clinicalData,header=TRUE,sep="\t",row.names=1)
#dim(allExpression)
#dim(pl)
#dim(clinical)
survival<-clinical[,'YEARS_TO_EVENT']
#print(rownames(clinical))
rownames(clinical) <- gsub("^","X",rownames(clinical))
#colnames(survival) <- rownames(clinical)
#print(head(clinical))

row.names(pl)<-pl$V1
#head(pl)[1]
pl.expression <-allExpression[row.names(pl),]
pl.expression<-na.omit(pl.expression)
#dim(pl.expression)
#head(pl.expression)
#head(clinical)
#print(clinical$YEARS_TO_EVENT)

probeNames <- row.names(pl.expression)
print(probeNames)
for(i in 1:length(probeNames)){
    probe <- probeNames[i]
    print ("=====================================================")
    print(probe)
    print(pl[probe,2])
    pl.expression<-pl.expression[,order(pl.expression[probe,])]
#    survival<-clinical[colnames(pl.expression),'YEARS_TO_EVENT']
#    head(survival)
    print("creating probe-specific data frame")
    df <- cbind(t(pl.expression[probe,]),clinical[colnames(pl.expression),'YEARS_TO_EVENT'],clinical[colnames(pl.expression),'X_EVENT'],clinical[colnames(pl.expression),'MLL_STATUS'])
    print(head(df))
    colnames(df) <- c("probe_expr","years_to_event","event","mll_status")
    print(head(df))
    write.table(df,paste(probe,".data.txt",sep=""),sep="\t")
    pdf(paste(probe,".heatmap.pdf",sep=""))
    heatmap.2(df,trace="none",col=colorRampPalette(c("blue","white","red"))(75),main=pl[probe,2],dendrogram="none",Colv="NA",Rowv="NA",cexCol=0.8)
    dev.off()
    print("Up (log2 >= 0.58)")
    up <- pl.expression[i,] >= 0.58
    if ((sum(up)) > 0){
#        print(head(up))
        print(colnames(pl.expression[,up]))
        print(clinical[colnames(pl.expression[,up]),'YEARS_TO_EVENT'])
        print(paste("mean:",mean(clinical[colnames(pl.expression[,up]),'YEARS_TO_EVENT'])),sep=" ")
        print(paste("stdev:",sd(clinical[colnames(pl.expression[,up]),'YEARS_TO_EVENT'])),sep=" ")
    } else { print("No up genes.")}
    print("Down (log2 <= -0.41)")
    down <- pl.expression[i,] <= -0.41
    if (sum(down) > 0){
    #print(sum(down))
        print(colnames(pl.expression[,down]))
        print(clinical[colnames(pl.expression[,down]),'YEARS_TO_EVENT'])
        print(paste("mean:",mean(clinical[colnames(pl.expression[,down]),'YEARS_TO_EVENT'])),sep=" ")
        print(paste("stdev:",sd(clinical[colnames(pl.expression[,down]),'YEARS_TO_EVENT'])),sep=" ")
    } else { print("No down genes.")}
    print ("Playing with survival")
    survival<-as.data.frame(df)
    for(j in 1:length(rownames(survival))){
        if (survival[j,'probe_expr'] >= 0.58) {
            survival[j,'probe_expr']<-1
        }else if(survival[j,'probe_expr']<= -0.41){
            survival[j,'probe_expr'] <--1
        }else{
            survival[j,'probe_expr'] <- 0
        }
    }
    print(head(survival))
    pdf(paste(probe,".kmplot.pdf",sep=""))
    survival$SurvObj <- with(survival, Surv(years_to_event,event))
    print(head(survival))
    km.by.expr <- survfit(SurvObj ~ probe_expr,data=survival)
    plot(km.by.expr,col=c("blue","black","red"),main=pl[probe,2],xlab="Years", ylab="Survival Percentage")
    legend("bottomleft", c(paste("Reduced Expr ",sum(down),sep=""),paste("Unchanged ",length(colnames(pl.expression))-(sum(down)+sum(up)),sep=""),paste("Increased Expr ",sum(up),sep="")), col=c("blue","black","red"),lty=c(1,1,1))
    dev.off()
}

negative <- survival[,'mll_status'] == 1
positive <- survival[,'mll_status'] == 2
pdf("MLLstatus.kmplot.pdf")
survival$SurvObj <- with(survival, Surv(years_to_event,event))
print(head(survival))
km.by.expr <- survfit(SurvObj ~ mll_status,data=survival)
plot(km.by.expr,col=c("blue","red"),main="MLL Status",xlab="Years", ylab="Survival Percentage")
legend("bottomleft", c(paste("MLL Negative ",sum(negative),sep=""),
                       paste("MLL Positive ",sum(positive),sep="")),
       col=c("blue","red"),lty=c(1,1,1))
dev.off()
