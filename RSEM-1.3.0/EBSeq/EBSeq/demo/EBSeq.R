library(EBSeq)
# 3.1
data(GeneMat)
str(GeneMat)
Sizes=MedianNorm(GeneMat)
EBOut=EBTest(Data=GeneMat,
			 Conditions=as.factor(rep(c("C1","C2"),each=5)),sizeFactors=Sizes, maxround=5)
DEOut=GetDEResults(EBOut)
str(DEOut)
#3.2
data(IsoList)
str(IsoList)
IsoMat=IsoList$IsoMat
str(IsoMat)
IsoNames=IsoList$IsoNames
IsosGeneNames=IsoList$IsosGeneNames
IsoSizes=MedianNorm(IsoMat)
NgList=GetNg(IsoNames, IsosGeneNames)
IsoNgTrun=NgList$IsoformNgTrun
IsoNgTrun[c(1:3,201:203,601:603)]
IsoEBOut=EBTest(Data=IsoMat, NgVector=IsoNgTrun,
				Conditions=as.factor(rep(c("C1","C2"),each=5)),sizeFactors=IsoSizes, maxround=5)
IsoDE=GetDEResults(IsoEBOut)
str(IsoDE)
#3.3
data(MultiGeneMat)
str(MultiGeneMat)
Conditions=c("C1","C1","C2","C2","C3","C3")
PosParti=GetPatterns(Conditions)
PosParti
Parti=PosParti[-3,]
Parti
MultiSize=MedianNorm(MultiGeneMat)
MultiOut=EBMultiTest(MultiGeneMat,NgVector=NULL,Conditions=Conditions,
					            AllParti=Parti, sizeFactors=MultiSize, maxround=5)
MultiPP=GetMultiPP(MultiOut)
names(MultiPP)
MultiPP$PP[1:10,]
MultiPP$MAP[1:10]
MultiPP$Patterns

#3.4
data(IsoMultiList)
IsoMultiMat=IsoMultiList[[1]]
IsoNames.Multi=IsoMultiList$IsoNames
IsosGeneNames.Multi=IsoMultiList$IsosGeneNames
IsoMultiSize=MedianNorm(IsoMultiMat)
NgList.Multi=GetNg(IsoNames.Multi, IsosGeneNames.Multi)
IsoNgTrun.Multi=NgList.Multi$IsoformNgTrun
Conditions=c("C1","C1","C2","C2","C3","C3","C4","C4")
PosParti.4Cond=GetPatterns(Conditions)
PosParti.4Cond
Parti.4Cond=PosParti.4Cond[c(1,2,3,8,15),]
Parti.4Cond
IsoMultiOut=EBMultiTest(IsoMultiMat,NgVector=IsoNgTrun.Multi,Conditions=Conditions,
					            AllParti=Parti.4Cond, sizeFactors=IsoMultiSize, maxround=5)
IsoMultiPP=GetMultiPP(IsoMultiOut)
names(MultiPP)
IsoMultiPP$PP[1:10,]
IsoMultiPP$MAP[1:10]
IsoMultiPP$Patterns


#4.1
data(GeneMat)
str(GeneMat)
Sizes=MedianNorm(GeneMat)
EBOut=EBTest(Data=GeneMat,
			 Conditions=as.factor(rep(c("C1","C2"),each=5)),sizeFactors=Sizes, maxround=5)
DEOut=GetDEResults(EBOut)
EBOut$Alpha
EBOut$Beta
EBOut$P
GeneFC=PostFC(EBOut)
str(GeneFC)
par(mfrow=c(2,2))
QQP(EBOut)
par(mfrow=c(2,2))
DenNHist(EBOut)
PlotPostVsRawFC(EBOut,GeneFC)

#4.2
data(IsoList)
str(IsoList)
IsoMat=IsoList$IsoMat
str(IsoMat)
IsoNames=IsoList$IsoNames
IsosGeneNames=IsoList$IsosGeneNames
IsoSizes=MedianNorm(IsoMat)
NgList=GetNg(IsoNames, IsosGeneNames)
IsoNgTrun=NgList$IsoformNgTrun
IsoNgTrun[c(1:3,201:203,601:603)]
IsoEBOut=EBTest(Data=IsoMat, NgVector=IsoNgTrun,
				Conditions=as.factor(rep(c("C1","C2"),each=5)),sizeFactors=IsoSizes, maxround=5)
IsoDE=GetDEResults(IsoEBOut)
str(IsoDE)
IsoEBOut$Alpha
IsoEBOut$Beta
IsoEBOut$P
IsoFC=PostFC(IsoEBOut)
str(IsoFC)
PlotPostVsRawFC(IsoEBOut,IsoFC)

par(mfrow=c(2,2))
PolyFitValue=vector("list",3)
for(i in 1:3)
	          PolyFitValue[[i]]=PolyFitPlot(IsoEBOut$C1Mean[[i]],
											                IsoEBOut$C1EstVar[[i]],5)
PolyAll=PolyFitPlot(unlist(IsoEBOut$C1Mean), unlist(IsoEBOut$C1EstVar),5)
lines(log10(IsoEBOut$C1Mean[[1]][PolyFitValue[[1]]$sort]),
	  PolyFitValue[[1]]$fit[PolyFitValue[[1]]$sort],col="yellow",lwd=2)
lines(log10(IsoEBOut$C1Mean[[2]][PolyFitValue[[2]]$sort]),
	  PolyFitValue[[2]]$fit[PolyFitValue[[2]]$sort],col="pink",lwd=2)
lines(log10(IsoEBOut$C1Mean[[3]][PolyFitValue[[3]]$sort]),
	  PolyFitValue[[3]]$fit[PolyFitValue[[3]]$sort],col="green",lwd=2)
legend("topleft",c("All Isoforms","Ig = 1","Ig = 2","Ig = 3"),
	   col=c("red","yellow","pink","green"),lty=1,lwd=3,box.lwd=2)
par(mfrow=c(2,3))
QQP(IsoEBOut)
par(mfrow=c(2,3))
DenNHist(IsoEBOut)


#4.3
data(MultiGeneMat)
str(MultiGeneMat)
Conditions=c("C1","C1","C2","C2","C3","C3")
PosParti=GetPatterns(Conditions)
PosParti
PlotPattern(PosParti)
Parti=PosParti[-3,]
Parti
MultiSize=MedianNorm(MultiGeneMat)
MultiOut=EBMultiTest(MultiGeneMat,NgVector=NULL,Conditions=Conditions,
					            AllParti=Parti, sizeFactors=MultiSize, maxround=5)
MultiPP=GetMultiPP(MultiOut)
names(MultiPP)
MultiPP$PP[1:10,]
MultiPP$MAP[1:10]
MultiPP$Patterns
MultiFC=GetMultiFC(MultiOut)
str(MultiFC)
par(mfrow=c(2,2))
DenNHist(MultiOut)
par(mfrow=c(2,2))
QQP(MultiOut)

#4.4
data(IsoMultiList)
IsoMultiMat=IsoMultiList[[1]]
IsoNames.Multi=IsoMultiList$IsoNames
IsosGeneNames.Multi=IsoMultiList$IsosGeneNames
IsoMultiSize=MedianNorm(IsoMultiMat)
NgList.Multi=GetNg(IsoNames.Multi, IsosGeneNames.Multi)
IsoNgTrun.Multi=NgList.Multi$IsoformNgTrun
Conditions=c("C1","C1","C2","C2","C3","C3","C4","C4")
PosParti.4Cond=GetPatterns(Conditions)
PosParti.4Cond
PlotPattern(PosParti.4Cond)
Parti.4Cond=PosParti.4Cond[c(1,2,3,8,15),]
Parti.4Cond
IsoMultiOut=EBMultiTest(IsoMultiMat,NgVector=IsoNgTrun.Multi,Conditions=Conditions,
					            AllParti=Parti.4Cond, sizeFactors=IsoMultiSize, maxround=5)
IsoMultiPP=GetMultiPP(IsoMultiOut)
names(MultiPP)
IsoMultiPP$PP[1:10,]
IsoMultiPP$MAP[1:10]
IsoMultiPP$Patterns
IsoMultiFC=GetMultiFC(IsoMultiOut)
str(IsoMultiFC)
par(mfrow=c(3,4))
DenNHist(IsoMultiOut)
par(mfrow=c(3,4))
QQP(IsoMultiOut)
IsoMultiFC=GetMultiFC(IsoMultiOut)



#4.5
data(GeneMat)
GeneMat.norep=GeneMat[,c(1,6)]
Sizes.norep=MedianNorm(GeneMat.norep)
EBOut.norep=EBTest(Data=GeneMat.norep,
			 Conditions=as.factor(rep(c("C1","C2"))),sizeFactors=Sizes.norep, maxround=5)
DE.norep=GetDEResults(EBOut.norep)
GeneFC.norep=PostFC(EBOut.norep)


#4.6
data(IsoList)
IsoMat=IsoList$IsoMat
IsoNames=IsoList$IsoNames
IsosGeneNames=IsoList$IsosGeneNames
NgList=GetNg(IsoNames, IsosGeneNames)
IsoNgTrun=NgList$IsoformNgTrun
IsoMat.norep=IsoMat[,c(1,6)]
IsoSizes.norep=MedianNorm(IsoMat.norep)
IsoEBOut.norep=EBTest(Data=IsoMat.norep, NgVector=IsoNgTrun,
				Conditions=as.factor(c("C1","C2")),sizeFactors=IsoSizes.norep, maxround=5)
IsoDE.norep=GetDEResults(IsoEBOut.norep)
IsoFC.norep=PostFC(IsoEBOut.norep)


#4.7
data(MultiGeneMat)
MultiGeneMat.norep=MultiGeneMat[,c(1,3,5)]
Conditions=c("C1","C2","C3")
PosParti=GetPatterns(Conditions)
Parti=PosParti[-3,]
MultiSize.norep=MedianNorm(MultiGeneMat.norep)
MultiOut.norep=EBMultiTest(MultiGeneMat.norep,NgVector=NULL,Conditions=Conditions,
					            AllParti=Parti, sizeFactors=MultiSize.norep, maxround=5)
MultiPP.norep=GetMultiPP(MultiOut.norep)
MultiFC.norep=GetMultiFC(MultiOut.norep)

#4.8
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
IsoMultiOut.norep=EBMultiTest(IsoMultiMat.norep,NgVector=IsoNgTrun.Multi,Conditions=Conditions,
					            AllParti=Parti.4Cond, sizeFactors=IsoMultiSize.norep, maxround=5)
IsoMultiPP.norep=GetMultiPP(IsoMultiOut.norep)
IsoMultiFC.norep=GetMultiFC(IsoMultiOut.norep)


# EOF
