### R code from vignette source 'EBSeq_Vignette.Rnw'

###################################################
### code chunk number 1: EBSeq_Vignette.Rnw:172-173
###################################################
library(EBSeq)


###################################################
### code chunk number 2: EBSeq_Vignette.Rnw:198-200
###################################################
data(GeneMat)
str(GeneMat)


###################################################
### code chunk number 3: EBSeq_Vignette.Rnw:208-209
###################################################
Sizes=MedianNorm(GeneMat)


###################################################
### code chunk number 4: EBSeq_Vignette.Rnw:235-237
###################################################
EBOut=EBTest(Data=GeneMat, 
Conditions=as.factor(rep(c("C1","C2"),each=5)),sizeFactors=Sizes, maxround=5)


###################################################
### code chunk number 5: EBSeq_Vignette.Rnw:240-244
###################################################
EBDERes=GetDEResults(EBOut, FDR=0.05)
str(EBDERes$DEfound)
head(EBDERes$PPMat)
str(EBDERes$Status)


###################################################
### code chunk number 6: EBSeq_Vignette.Rnw:289-295
###################################################
data(IsoList)
str(IsoList)
IsoMat=IsoList$IsoMat
str(IsoMat)
IsoNames=IsoList$IsoNames
IsosGeneNames=IsoList$IsosGeneNames


###################################################
### code chunk number 7: EBSeq_Vignette.Rnw:302-303
###################################################
IsoSizes=MedianNorm(IsoMat)


###################################################
### code chunk number 8: EBSeq_Vignette.Rnw:324-327
###################################################
NgList=GetNg(IsoNames, IsosGeneNames)
IsoNgTrun=NgList$IsoformNgTrun
IsoNgTrun[c(1:3,201:203,601:603)]


###################################################
### code chunk number 9: EBSeq_Vignette.Rnw:339-345
###################################################
IsoEBOut=EBTest(Data=IsoMat, NgVector=IsoNgTrun, 
Conditions=as.factor(rep(c("C1","C2"),each=5)),sizeFactors=IsoSizes, maxround=5)
IsoEBDERes=GetDEResults(IsoEBOut, FDR=0.05)
str(IsoEBDERes$DEfound)
head(IsoEBDERes$PPMat)
str(IsoEBDERes$Status)


###################################################
### code chunk number 10: EBSeq_Vignette.Rnw:368-370
###################################################
data(MultiGeneMat)
str(MultiGeneMat)


###################################################
### code chunk number 11: EBSeq_Vignette.Rnw:378-381
###################################################
Conditions=c("C1","C1","C2","C2","C3","C3")
PosParti=GetPatterns(Conditions)
PosParti


###################################################
### code chunk number 12: EBSeq_Vignette.Rnw:389-391
###################################################
Parti=PosParti[-3,]
Parti


###################################################
### code chunk number 13: EBSeq_Vignette.Rnw:396-399
###################################################
MultiSize=MedianNorm(MultiGeneMat)
MultiOut=EBMultiTest(MultiGeneMat,NgVector=NULL,Conditions=Conditions,
AllParti=Parti, sizeFactors=MultiSize, maxround=5)


###################################################
### code chunk number 14: EBSeq_Vignette.Rnw:403-408
###################################################
MultiPP=GetMultiPP(MultiOut)
names(MultiPP)
MultiPP$PP[1:10,]
MultiPP$MAP[1:10]
MultiPP$Patterns


###################################################
### code chunk number 15: EBSeq_Vignette.Rnw:427-435
###################################################
data(IsoMultiList)
IsoMultiMat=IsoMultiList[[1]]
IsoNames.Multi=IsoMultiList$IsoNames
IsosGeneNames.Multi=IsoMultiList$IsosGeneNames
IsoMultiSize=MedianNorm(IsoMultiMat)
NgList.Multi=GetNg(IsoNames.Multi, IsosGeneNames.Multi)
IsoNgTrun.Multi=NgList.Multi$IsoformNgTrun
Conditions=c("C1","C1","C2","C2","C3","C3","C4","C4")


###################################################
### code chunk number 16: EBSeq_Vignette.Rnw:441-443
###################################################
PosParti.4Cond=GetPatterns(Conditions)
PosParti.4Cond


###################################################
### code chunk number 17: EBSeq_Vignette.Rnw:448-450
###################################################
Parti.4Cond=PosParti.4Cond[c(1,2,3,8,15),]
Parti.4Cond


###################################################
### code chunk number 18: EBSeq_Vignette.Rnw:455-459
###################################################
IsoMultiOut=EBMultiTest(IsoMultiMat,
NgVector=IsoNgTrun.Multi,Conditions=Conditions,
AllParti=Parti.4Cond, sizeFactors=IsoMultiSize, 
maxround=5)


###################################################
### code chunk number 19: EBSeq_Vignette.Rnw:463-468
###################################################
IsoMultiPP=GetMultiPP(IsoMultiOut)
names(MultiPP)
IsoMultiPP$PP[1:10,]
IsoMultiPP$MAP[1:10]
IsoMultiPP$Patterns


###################################################
### code chunk number 20: EBSeq_Vignette.Rnw:485-490 (eval = FALSE)
###################################################
## data(GeneMat)
## Sizes=MedianNorm(GeneMat)
## EBOut=EBTest(Data=GeneMat, 
## Conditions=as.factor(rep(c("C1","C2"),each=5)),sizeFactors=Sizes, maxround=5)
## EBDERes=GetDEResults(EBOut, FDR=0.05)


###################################################
### code chunk number 21: EBSeq_Vignette.Rnw:492-496
###################################################
EBDERes=GetDEResults(EBOut, FDR=0.05)
str(EBDERes$DEfound)
head(EBDERes$PPMat)
str(EBDERes$Status)


###################################################
### code chunk number 22: EBSeq_Vignette.Rnw:506-509
###################################################
GeneFC=PostFC(EBOut)
str(GeneFC)
PlotPostVsRawFC(EBOut,GeneFC)


###################################################
### code chunk number 23: EBSeq_Vignette.Rnw:530-533
###################################################
EBOut$Alpha
EBOut$Beta
EBOut$P


###################################################
### code chunk number 24: EBSeq_Vignette.Rnw:552-554
###################################################
par(mfrow=c(1,2))
QQP(EBOut)


###################################################
### code chunk number 25: EBSeq_Vignette.Rnw:570-572
###################################################
par(mfrow=c(1,2))
DenNHist(EBOut)


###################################################
### code chunk number 26: EBSeq_Vignette.Rnw:593-598 (eval = FALSE)
###################################################
## data(IsoList)
## IsoMat=IsoList$IsoMat
## IsoNames=IsoList$IsoNames
## IsosGeneNames=IsoList$IsosGeneNames
## NgList=GetNg(IsoNames, IsosGeneNames, TrunThre=3)


###################################################
### code chunk number 27: EBSeq_Vignette.Rnw:600-603
###################################################
names(NgList)
IsoNgTrun=NgList$IsoformNgTrun
IsoNgTrun[c(1:3,201:203,601:603)]


###################################################
### code chunk number 28: EBSeq_Vignette.Rnw:634-635 (eval = FALSE)
###################################################
## IsoNgTrun = scan(file="output_name.ngvec", what=0, sep="\n")


###################################################
### code chunk number 29: EBSeq_Vignette.Rnw:648-652 (eval = FALSE)
###################################################
## IsoSizes=MedianNorm(IsoMat)
## IsoEBOut=EBTest(Data=IsoMat, NgVector=IsoNgTrun, 
## Conditions=as.factor(rep(c("C1","C2"),each=5)),sizeFactors=IsoSizes, maxround=5)
## IsoEBDERes=GetDEResults(IsoEBOut, FDR=0.05)


###################################################
### code chunk number 30: EBSeq_Vignette.Rnw:654-655
###################################################
str(IsoEBDERes)


###################################################
### code chunk number 31: EBSeq_Vignette.Rnw:660-662
###################################################
IsoFC=PostFC(IsoEBOut)
str(IsoFC)


###################################################
### code chunk number 32: EBSeq_Vignette.Rnw:673-676
###################################################
IsoEBOut$Alpha
IsoEBOut$Beta
IsoEBOut$P


###################################################
### code chunk number 33: EBSeq_Vignette.Rnw:695-700
###################################################
par(mfrow=c(2,2))
PolyFitValue=vector("list",3)
for(i in 1:3)
    PolyFitValue[[i]]=PolyFitPlot(IsoEBOut$C1Mean[[i]], 
    IsoEBOut$C1EstVar[[i]],5)


###################################################
### code chunk number 34: EBSeq_Vignette.Rnw:713-722
###################################################
PolyAll=PolyFitPlot(unlist(IsoEBOut$C1Mean), unlist(IsoEBOut$C1EstVar),5)
lines(log10(IsoEBOut$C1Mean[[1]][PolyFitValue[[1]]$sort]), 
PolyFitValue[[1]]$fit[PolyFitValue[[1]]$sort],col="yellow",lwd=2)
lines(log10(IsoEBOut$C1Mean[[2]][PolyFitValue[[2]]$sort]), 
PolyFitValue[[2]]$fit[PolyFitValue[[2]]$sort],col="pink",lwd=2)
lines(log10(IsoEBOut$C1Mean[[3]][PolyFitValue[[3]]$sort]), 
PolyFitValue[[3]]$fit[PolyFitValue[[3]]$sort],col="green",lwd=2)
legend("topleft",c("All Isoforms","Ng = 1","Ng = 2","Ng = 3"),
col=c("red","yellow","pink","green"),lty=1,lwd=3,box.lwd=2)


###################################################
### code chunk number 35: EBSeq_Vignette.Rnw:735-737
###################################################
par(mfrow=c(2,3))
QQP(IsoEBOut)


###################################################
### code chunk number 36: EBSeq_Vignette.Rnw:749-751
###################################################
par(mfrow=c(2,3))
DenNHist(IsoEBOut)


###################################################
### code chunk number 37: EBSeq_Vignette.Rnw:768-772
###################################################
Conditions=c("C1","C1","C2","C2","C3","C3")
PosParti=GetPatterns(Conditions)
PosParti
PlotPattern(PosParti)


###################################################
### code chunk number 38: EBSeq_Vignette.Rnw:779-781
###################################################
Parti=PosParti[-3,]
Parti


###################################################
### code chunk number 39: EBSeq_Vignette.Rnw:787-793 (eval = FALSE)
###################################################
## data(MultiGeneMat)
## MultiSize=MedianNorm(MultiGeneMat)
## MultiOut=EBMultiTest(MultiGeneMat,
## NgVector=NULL,Conditions=Conditions,
## AllParti=Parti, sizeFactors=MultiSize, 
## maxround=5)


###################################################
### code chunk number 40: EBSeq_Vignette.Rnw:797-802
###################################################
MultiPP=GetMultiPP(MultiOut)
names(MultiPP)
MultiPP$PP[1:10,]
MultiPP$MAP[1:10]
MultiPP$Patterns


###################################################
### code chunk number 41: EBSeq_Vignette.Rnw:809-811
###################################################
MultiFC=GetMultiFC(MultiOut)
str(MultiFC)


###################################################
### code chunk number 42: EBSeq_Vignette.Rnw:820-822
###################################################
par(mfrow=c(2,2))
QQP(MultiOut)


###################################################
### code chunk number 43: EBSeq_Vignette.Rnw:830-832
###################################################
par(mfrow=c(2,2))
DenNHist(MultiOut)


###################################################
### code chunk number 44: EBSeq_Vignette.Rnw:847-850
###################################################
Conditions=c("C1","C1","C2","C2","C3","C3","C4","C4")
PosParti.4Cond=GetPatterns(Conditions)
PosParti.4Cond


###################################################
### code chunk number 45: EBSeq_Vignette.Rnw:855-858
###################################################
PlotPattern(PosParti.4Cond)
Parti.4Cond=PosParti.4Cond[c(1,2,3,8,15),]
Parti.4Cond


###################################################
### code chunk number 46: EBSeq_Vignette.Rnw:865-876 (eval = FALSE)
###################################################
## data(IsoMultiList)
## IsoMultiMat=IsoMultiList[[1]]
## IsoNames.Multi=IsoMultiList$IsoNames
## IsosGeneNames.Multi=IsoMultiList$IsosGeneNames
## IsoMultiSize=MedianNorm(IsoMultiMat)
## NgList.Multi=GetNg(IsoNames.Multi, IsosGeneNames.Multi)
## IsoNgTrun.Multi=NgList.Multi$IsoformNgTrun
## IsoMultiOut=EBMultiTest(IsoMultiMat,NgVector=IsoNgTrun.Multi,Conditions=Conditions,
## AllParti=Parti.4Cond, 
## sizeFactors=IsoMultiSize, maxround=5)
## IsoMultiPP=GetMultiPP(IsoMultiOut)


###################################################
### code chunk number 47: EBSeq_Vignette.Rnw:878-883
###################################################
names(MultiPP)
IsoMultiPP$PP[1:10,]
IsoMultiPP$MAP[1:10]
IsoMultiPP$Patterns
IsoMultiFC=GetMultiFC(IsoMultiOut)


###################################################
### code chunk number 48: EBSeq_Vignette.Rnw:894-897
###################################################
par(mfrow=c(3,4))
QQP(IsoMultiOut)



###################################################
### code chunk number 49: EBSeq_Vignette.Rnw:907-909
###################################################
par(mfrow=c(3,4))
DenNHist(IsoMultiOut)


###################################################
### code chunk number 50: EBSeq_Vignette.Rnw:941-949
###################################################
data(GeneMat)
GeneMat.norep=GeneMat[,c(1,6)]
Sizes.norep=MedianNorm(GeneMat.norep)
EBOut.norep=EBTest(Data=GeneMat.norep,
Conditions=as.factor(rep(c("C1","C2"))),
sizeFactors=Sizes.norep, maxround=5)
EBDERes.norep=GetDEResults(EBOut.norep)
GeneFC.norep=PostFC(EBOut.norep)


###################################################
### code chunk number 51: EBSeq_Vignette.Rnw:959-972
###################################################
data(IsoList)
IsoMat=IsoList$IsoMat
IsoNames=IsoList$IsoNames
IsosGeneNames=IsoList$IsosGeneNames
NgList=GetNg(IsoNames, IsosGeneNames)
IsoNgTrun=NgList$IsoformNgTrun
IsoMat.norep=IsoMat[,c(1,6)]
IsoSizes.norep=MedianNorm(IsoMat.norep)
IsoEBOut.norep=EBTest(Data=IsoMat.norep, NgVector=IsoNgTrun,
Conditions=as.factor(c("C1","C2")),
sizeFactors=IsoSizes.norep, maxround=5)
IsoEBDERes.norep=GetDEResults(IsoEBOut.norep)
IsoFC.norep=PostFC(IsoEBOut.norep)


###################################################
### code chunk number 52: EBSeq_Vignette.Rnw:981-993
###################################################
data(MultiGeneMat)
MultiGeneMat.norep=MultiGeneMat[,c(1,3,5)]
Conditions=c("C1","C2","C3")
PosParti=GetPatterns(Conditions)
Parti=PosParti[-3,]
MultiSize.norep=MedianNorm(MultiGeneMat.norep)
MultiOut.norep=EBMultiTest(MultiGeneMat.norep,
NgVector=NULL,Conditions=Conditions,
AllParti=Parti, sizeFactors=MultiSize.norep, 
maxround=5)
MultiPP.norep=GetMultiPP(MultiOut.norep)
MultiFC.norep=GetMultiFC(MultiOut.norep)


###################################################
### code chunk number 53: EBSeq_Vignette.Rnw:1005-1024
###################################################
data(IsoMultiList)
IsoMultiMat=IsoMultiList[[1]]
IsoNames.Multi=IsoMultiList$IsoNames
IsosGeneNames.Multi=IsoMultiList$IsosGeneNames
IsoMultiMat.norep=IsoMultiMat[,c(1,3,5,7)]
IsoMultiSize.norep=MedianNorm(IsoMultiMat.norep)
NgList.Multi=GetNg(IsoNames.Multi, IsosGeneNames.Multi)
IsoNgTrun.Multi=NgList.Multi$IsoformNgTrun
Conditions=c("C1","C2","C3","C4")
PosParti.4Cond=GetPatterns(Conditions)
PosParti.4Cond
Parti.4Cond=PosParti.4Cond[c(1,2,3,8,15),]
Parti.4Cond
IsoMultiOut.norep=EBMultiTest(IsoMultiMat.norep,
NgVector=IsoNgTrun.Multi,Conditions=Conditions,
AllParti=Parti.4Cond, sizeFactors=IsoMultiSize.norep, 
maxround=5)
IsoMultiPP.norep=GetMultiPP(IsoMultiOut.norep)
IsoMultiFC.norep=GetMultiFC(IsoMultiOut.norep)


