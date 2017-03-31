# peak picking
library(playwith)

# Plotting ====================================Denoising part ====================================================

playwith({
			plot((1:length(curVecInt))*params$ScanInterval+params$delaytime,curVecInt,lwd=2,type="l",main=paste(smoothingWindow,"-Point Moving Average Filter",sep=""))
			points((1:length(curVecInt))*params$ScanInterval+params$delaytime,sn.ma,lwd=2,type="l",col="red")
		})



playwith(
		{ 
			
			plot((1:length(curVecInt))*params$ScanInterval+params$delaytime,curVecInt, lwd=2,type = "l", col = "blue")
			points((1:length(curVecInt))*params$ScanInterval+params$delaytime,bgs,lwd=2,type="l")
			points((1:length(curVecInt))*params$ScanInterval+params$delaytime,bsln,lwd=2,type="l",col = "red")
		})


# Plotting ====================================Peak picking ===========================================================================
delaytime <- params$delaytime
ScanInterval <- params$ScanInterval

playwith({
			
			plot(((1:nPoints)-1)*ScanInterval+delaytime,curVecInt,type="l",lwd=2,main=paste(fileName,"_TIC_", "_peak picking results",sep=""),xlab="ET",ylab="Intensity")
			points((peakInd-1)*ScanInterval+delaytime,curVecInt[peakInd],col="red",lwd=3)
			#points((valleyInd-1)*ScanInterval+delaytime,curVecInt[valleyInd],col="blue",lwd=3)
		},new=TRUE)

playwith({
			
			plot(((1:nPoints)-1)*ScanInterval+delaytime,curVecInt,type="l",lwd=2,main=paste(fileName,"_mass_",mz, "_peak picking results",sep=""),xlab="ET",ylab="Intensity")
			#points((peakInd-1)*ScanInterval+delaytime,curVecInt[peakInd],col="red",lwd=3)
			#points((valleyInd-1)*ScanInterval+delaytime,curVecInt[valleyInd],col="blue",lwd=3)
			#text((curPeakList$pkInd-1)*ScanInterval+delaytime, curPeakList$Intensity,as.character(curPeakList$pkInd))
		},new=TRUE)


# Plotting =================================================================================================================================

playwith({
			
			plot(((1:nPoints)-1)*ScanInterval+delaytime,curVecInt,type="l",lwd=2,main=paste(fileName,"_mass_",mz, "_peak picking results",sep=""),xlab="ET",ylab="Intensity")
			points((peakInd-1)*ScanInterval+delaytime,curVecInt[peakInd],col="red",lwd=3)
			points((valleyInd-1)*ScanInterval+delaytime,curVecInt[valleyInd],col="blue",lwd=3)
			#abline(v=(Update_Peak_list$Lbound-1)*ScanInterval+delaytime,col="purple")
			#abline(v=(Update_Peak_list$Rbound-1)*ScanInterval+delaytime,col="green")
			#text((curPeakList$pkInd-1)*ScanInterval+delaytime, curPeakList$Intensity,as.character(curPeakList$pkInd))
			#abline(v=(subset(curPeakList,isApex==1,select=c("lboundInd"))$lboundInd-1)*ScanInterval+delaytime,col="blue")
			#abline(v=(subset(curPeakList,isApex==1,select=c("rboundInd"))$rboundInd-1)*ScanInterval+delaytime,col="blue")
			
		},new=TRUE)

# Plotting =================================================================================================================================

x11()
hist(sort(Initial_Peak_list$Width),main="peak width distribution")
if (nrow(Initial_Peak_list)>=500) {
	legend("topright",paste("The 500th:", as.character(sort(Initial_Peak_list$Width,decreasing=TRUE)[500]),"& mean:",as.character(mean(sort(Initial_Peak_list$Width,decreasing=TRUE)[1:500])),sep=" "))
}


# Plotting =================================================================================================================================	
# TIC
curPeakList$ET<-(curPeakList$pkInd-1)*ScanInterval+delaytime
playwith({
			
			plot(((1:nPoints)-1)*ScanInterval+delaytime,curVecInt,type="l",lwd=2,main=paste(fileName,"_TIC_", "_peak picking results",sep=""),xlab="ET",ylab="Intensity")
			points(subset(curPeakList,isApex==1&isShared==1,select=c("ET","Intensity")),col="red",lwd=3)
			points(subset(curPeakList,isApex==1&isShared==0,select=c("ET","Intensity")),col="blue",lwd=3)
			points(subset(curPeakList,isApex==-1,select=c("ET","Intensity")),col="orange",lwd=3)
			points(subset(curPeakList,isApex==-2,select=c("ET","Intensity")),col="green",lwd=3)					
			abline(v=(curPeakList$Lbound-1)*ScanInterval+delaytime,col="yellow")
			abline(v=(curPeakList$Rbound-1)*ScanInterval+delaytime,col="green")
			abline(v=(subset(curPeakList,isApex==1,select=c("lboundInd"))$lboundInd-1)*ScanInterval+delaytime,col="blue")
			abline(v=(subset(curPeakList,isApex==1,select=c("rboundInd"))$rboundInd-1)*ScanInterval+delaytime,col="blue")
		},new=TRUE)

# EIC
curPeakList$ET<-(curPeakList$pkInd-1)*ScanInterval+delaytime
playwith({
			
			plot(((1:nPoints)-1)*ScanInterval+delaytime,curVecInt,type="l",lwd=2,main=paste(fileName,"_mass_",mz, "_peak picking results",sep=""),xlab="ET",ylab="Intensity")
			points(subset(curPeakList,isApex==1&isShared==1,select=c("ET","Intensity")),col="red",lwd=3)
			points(subset(curPeakList,isApex==1&isShared==0,select=c("ET","Intensity")),col="blue",lwd=3)
			points(subset(curPeakList,isApex==-1,select=c("ET","Intensity")),col="orange",lwd=3)
			points(subset(curPeakList,isApex==-2,select=c("ET","Intensity")),col="green",lwd=3)					
			abline(v=(curPeakList$Lbound-1)*ScanInterval+delaytime,col="yellow")
			abline(v=(curPeakList$Rbound-1)*ScanInterval+delaytime,col="green")
			abline(v=(subset(curPeakList,isApex==1,select=c("lboundInd"))$lboundInd-1)*ScanInterval+delaytime,col="blue")
			abline(v=(subset(curPeakList,isApex==1,select=c("rboundInd"))$rboundInd-1)*ScanInterval+delaytime,col="blue")
		},new=TRUE)
	
# Plotting =================================================================================================================================	
curPeakList$ET<-(curPeakList$pkInd-1)*ScanInterval+delaytime
playwith({
			
			plot(((1:nPoints)-1)*ScanInterval+delaytime,curVecInt,type="l")
			points(subset(curPeakList,isApex==1&isShared==0,select=c("ET","Intensity")),col="red")
			points(subset(curPeakList,isApex==1&isShared==1,select=c("ET","Intensity")),col="blue")			
			points(subset(curPeakList,isApex==-1,select=c("ET","Intensity")),col="yellow")
			points(subset(curPeakList,isApex==-2,select=c("ET","Intensity")),col="green")
			abline(v=(RTtime-1)*ScanInterval+delaytime,col="blue")
			abline(v=(curPeakList$Lbound-1)*ScanInterval+delaytime,col="yellow")
			abline(v=(curPeakList$Rbound-1)*ScanInterval+delaytime,col="green")
			abline(v=(curPeakList$lboundInd-1)*ScanInterval+delaytime,col="grey")
			abline(v=(curPeakList$rboundInd-1)*ScanInterval+delaytime,col="grey")
			#points((usedtime-1)*ScanInterval+delaytime, c(rep(1000,length(usedtime))),col="purple")
			text((RTtime-1)*ScanInterval+delaytime, c(rep(1000,length(RTtime))),as.character(sigToNoise))	
			text((curPeakList$pkInd-1)*ScanInterval+delaytime, curPeakList$Intensity,as.character(curPeakList$StN))
		},new=TRUE)	

# Plotting =================================================================================================================================	

x11()
hist(curPeakList$gss)

# Plotting =========================================== Deconvolution ==============================================================			

# ploting decon window
playwith({
			
			plot((1:length(TIC)-1)*ScanInterval+delaytime,TIC,type="l",main=fileName,xlab="ET(Mins)")
			points(TICpeaklist$RT,TIC[TICpeaklist$pkInd],col="red")
			text(TICpeaklist$RT,TIC[TICpeaklist$pkInd],c(1:nrow(TICpeaklist)))
			abline(v=(TICpeaklist$lboundInd-1)*ScanInterval+delaytime,col="blue")
			abline(v=(TICpeaklist$rboundInd-1)*ScanInterval+delaytime,col="blue")
			
		})	


playwith({		
			plot((1:length(TIC)-1)*ScanInterval+delaytime,TIC,type="l",main=fileName,xlab="ET(Mins)",lwd=6,cex=1)		
		})	
# Plotting =================================================================================================================================
# to plot both TIC and local EIC
minLbound<-min(localEICPeakList$lboundInd)
maxRbound<-max(localEICPeakList$rboundInd)	
maxInt<-max(TIC[curTICPk$lboundInd:curTICPk$rboundInd])
colVec<-rainbow(nrow(localEICPeakList))
masses <-localEICPeakList$mz

playwith({
			plot(x=0,type = "n",xlim=((c(minLbound,maxRbound)-1)*ScanInterval+delaytime),ylim=c(0,maxInt),main=paste(fileName,"TIC and local EICs"),xlab="ET",ylab="Intensity")
			
			for(i in c(1:nrow(localEICPeakList)))	
			{
				
				curID<-rownames(localEICPeakList[i,])
				points((allprofiles[[curID]]$ET-1)*ScanInterval+delaytime,allprofiles[[curID]]$int, lwd=4,type = "l",col=colVec[i])
				if (localEICPeakList[i,]$isApex==1) {
					points((localEICPeakList[i,]$pkInd-1)*ScanInterval+delaytime,localEICPeakList[i,]$Intensity,col="red",lwd=3)
				}else {
					points((localEICPeakList[i,]$pkInd-1)*ScanInterval+delaytime,localEICPeakList[i,]$Intensity,col="black",lwd=3)
				}
				
			}
			points((c(curTICPk$lboundInd:curTICPk$rboundInd)-1)*ScanInterval+delaytime,TIC[curTICPk$lboundInd:curTICPk$rboundInd],col="grey",lwd=4,type="l")
			points((localTICApexList$pkInd-1)*ScanInterval+delaytime,TIC[localTICApexList$pkInd],col="blue")
			
			legend("topright",as.character(masses),fill=colVec)
		},new=TRUE)

# to plot TIC and EIC retention time distribution
localEICPeakList$ET<-(localEICPeakList$pkInd-1)*ScanInterval+delaytime
localTICApexList$ET<-(localTICApexList$pkInd-1)*ScanInterval+delaytime
localTICApexList$mz=600
minRT <- min(localEICPeakList$ET,localTICApexList$ET)
maxRT <- max(localEICPeakList$ET,localTICApexList$ET)

playwith({
			
			plot(subset(localEICPeakList,isShared==0,select=c("ET","mz")),,xlim=c(minRT,maxRT),ylim=c(0,700),col="red")
			points(subset(localEICPeakList,isShared==1,select=c("ET","mz")))
			points(subset(localTICApexList,isShared==1,select=c("ET","mz")),col="blue")
})

# plot unique peak profiles (1 cluster - 1 color)
minLbound<-min(uniqueEICPeakList$lboundInd)
maxRbound<-max(uniqueEICPeakList$rboundInd)	
maxInt<-max(TIC[curTICPk$lboundInd:curTICPk$rboundInd])
maxInt<-max(uniqueEICPeakList$Intensity)
colVec<-ClustResult 
masses <-uniqueEICPeakList$mz

playwith({
			plot(x=0,type = "n",xlim=((c(minLbound,maxRbound)-1)*ScanInterval+delaytime),ylim=c(0,maxInt),main=paste(fileName,"TIC and local EICs"),xlab="ET",ylab="Intensity")
			for(i in c(1:nrow(uniqueEICPeakList)))	
			{
				
				curID<-rownames(uniqueEICPeakList[i,])
				points((allprofiles[[curID]]$ET-1)*ScanInterval+delaytime,allprofiles[[curID]]$int, lwd=4,type = "l",col=colVec[i])
				if (uniqueEICPeakList[i,]$isApex==1) {
					points((uniqueEICPeakList[i,]$pkInd-1)*ScanInterval+delaytime,uniqueEICPeakList[i,]$Intensity,col="red",lwd=3)
				}else {
					points((uniqueEICPeakList[i,]$pkInd-1)*ScanInterval+delaytime,uniqueEICPeakList[i,]$Intensity,col="black",lwd=3)
				}	
			}				
			legend("topright",as.character(masses),fill=colVec)
		},new=TRUE)




uniqueEICPeakList <- subset(localEICPeakList,isShared==0)
minLbound<-min(uniqueEICPeakList$lboundInd)
maxRbound<-max(uniqueEICPeakList$rboundInd)	
maxInt<-max(TIC[curTICPk$lboundInd:curTICPk$rboundInd])
maxInt<-max(uniqueEICPeakList$Intensity)
colVec<-rainbow(nrow(uniqueEICPeakList))
masses <-uniqueEICPeakList$mz

playwith({
			plot(x=0,type = "n",xlim=((c(minLbound,maxRbound)-1)*ScanInterval+delaytime),ylim=c(0,maxInt),main=paste(fileName,"TIC and local EICs"),xlab="ET",ylab="Intensity")
			
			for(i in c(1:nrow(uniqueEICPeakList)))	
			{
				
				curID<-rownames(uniqueEICPeakList[i,])
				points((allprofiles[[curID]]$ET-1)*ScanInterval+delaytime,allprofiles[[curID]]$int, lwd=4,type = "l",col=colVec[i])
				if (uniqueEICPeakList[i,]$isApex==1) {
					points((uniqueEICPeakList[i,]$pkInd-1)*ScanInterval+delaytime,uniqueEICPeakList[i,]$Intensity,col="red",lwd=3)
				}else {
					points((uniqueEICPeakList[i,]$pkInd-1)*ScanInterval+delaytime,uniqueEICPeakList[i,]$Intensity,col="black",lwd=3)
				}
				
			}
#			points((c(curTICPk$lboundInd:curTICPk$rboundInd)-1)*ScanInterval+delaytime,TIC[curTICPk$lboundInd:curTICPk$rboundInd],col="grey",lwd=4,type="l")
#			points((localTICApexList$pkInd-1)*ScanInterval+delaytime,TIC[localTICApexList$pkInd],col="blue")
			
			legend("topright",as.character(masses),fill=colVec)
		},new=TRUE)



# to plot local EIC only
minLbound<-min(localEICPeakList$lboundInd)
maxRbound<-max(localEICPeakList$rboundInd)	
maxInt2<-max(localEICPeakList$Intensity)
colVec<-rainbow(nrow(localEICPeakList))
masses <-localEICPeakList$mz

playwith({
			#plot(x=0,type = "n",xlim=((c(minLbound,maxRbound)-1)*ScanInterval+delaytime),ylim=c(0,max(maxInt2)),main="good Peak selection",xlab="ET",ylab="Intensity")
			plot(x=0,type = "n",xlim=((c(minLbound,maxRbound)-1)*ScanInterval+delaytime),ylim=c(0,max(TIC[curTICPk$lboundInd:curTICPk$rboundInd])),main="good Peak selection",xlab="ET",ylab="Intensity")
			
			for(i in c(1:nrow(localEICPeakList)))	
			{
				
				curID<-rownames(localEICPeakList[i,])
				points((allprofiles[[curID]]$ET-1)*ScanInterval+delaytime,allprofiles[[curID]]$int, lwd=4,type = "l",col=colVec[i])
				if (localEICPeakList[i,]$isApex==1) {
					points((localEICPeakList[i,]$pkInd-1)*ScanInterval+delaytime,localEICPeakList[i,]$Intensity,col="red",lwd=3)
				}else {
					points((localEICPeakList[i,]$pkInd-1)*ScanInterval+delaytime,localEICPeakList[i,]$Intensity,col="black",lwd=3)
				}
				
			}
			points(c(curTICPk$lboundInd:curTICPk$rboundInd)*ScanInterval+delaytime,TIC[curTICPk$lboundInd:curTICPk$rboundInd],col="grey",lwd=4,type="l")
			legend("topright",as.character(masses),fill=colVec)
		})



# Plotting =================================================================================================================================
minLbound<-min(CandidatePeaklist$lboundInd)
maxRbound<-max(CandidatePeaklist$rboundInd)	
maxInt2<-max(CandidatePeaklist$Intensity)
colVec<-rainbow(nrow(CandidatePeaklist))
masses <-CandidatePeaklist$mz

playwith({
			plot(x=0,type = "n",xlim=((c(minLbound,maxRbound)-1)*ScanInterval+delaytime),ylim=c(0,max(maxInt2)),
					main=paste(fileName,"candidate peaks"),xlab="ET",ylab="Intensity")
			for(i in c(1:nrow(CandidatePeaklist)))	
			{
				
				curID<-rownames(CandidatePeaklist[i,])
				points((allprofiles[[curID]]$ET-1)*ScanInterval+delaytime,allprofiles[[curID]]$int, lwd=4,type = "l",col=colVec[i])
				if (CandidatePeaklist[i,]$isApex==1) {
					points((CandidatePeaklist[i,]$pkInd-1)*ScanInterval+delaytime,CandidatePeaklist[i,]$Intensity,col="red",lwd=3)
				}else {
					points((CandidatePeaklist[i,]$pkInd-1)*ScanInterval+delaytime,CandidatePeaklist[i,]$Intensity,col="black",lwd=3)
				}
				
			}
			legend("topright",as.character(masses),fill=colVec)
		})	

# Plotting =================================================================================================================================
maxInt2 <- max(goodShapePeaklist$Intensity)+500
minLbound<-min(goodShapePeaklist$lboundInd)
maxRbound<-max(goodShapePeaklist$rboundInd)	
#masses <-goodShapePeaklist$mz
colVec<-rainbow(nrow(goodShapePeaklist))

labelVec<-paste("pkID:",rownames(goodShapePeaklist),"mz:",goodShapePeaklist$mz," ,gss:",format(goodShapePeaklist$gss,digits=2)," ,shrp:",format(goodShapePeaklist$shrp,digits=2))
dataPoints<-subset(goodShapePeaklist,select=c("pkInd","Intensity"))
dataPoints$pkInd<-(dataPoints$pkInd-1)*ScanInterval+delaytime
playwith({
			plot(x=0,type = "n",xlim=((c(minLbound,maxRbound)-1)*ScanInterval+delaytime),ylim=c(0,max(maxInt2)),main=paste(fileName,"good peaks"),xlab="ET",ylab="Intensity")
			for(i in 1:nrow(goodShapePeaklist))
			{
				curID<-rownames(goodShapePeaklist[i,])
				points((profiles[[curID]]$ET-1)*ScanInterval+delaytime,profiles[[curID]]$int, lwd=3,type = "l",col=colVec[i])
				#points((goodShapePeaklist[i,]$pkInd-1)*ScanInterval+delaytime,goodShapePeaklist[i,]$Intensity,col="green",lwd=3)
			}
			legend("topright",as.character(goodShapePeaklist$mz),fill=colVec)
		},data.points=dataPoints,labels=labelVec)
# Plotting =================================================================================================================================

x11(width=15,height=10)
par(mfrow=c(2,4))

# plot 1: local EIC peak profiles
minLbound<-min(localEICPeakList$lboundInd)
maxRbound<-max(localEICPeakList$rboundInd)	
maxInt2<-max(localEICPeakList$Intensity)
colVec<-rainbow(nrow(localEICPeakList))
masses <-localEICPeakList$mz

plot(x=0,type = "n",xlim=((c(minLbound,maxRbound)-1)*ScanInterval+delaytime),ylim=c(0,max(maxInt2)),
		main=paste(fileName,"local EIC peaks"),xlab="ET",ylab="Intensity")
for(i in c(1:nrow(localEICPeakList)))	
{
	
	curID<-rownames(localEICPeakList[i,])
	points((allprofiles[[curID]]$ET-1)*ScanInterval+delaytime,allprofiles[[curID]]$int, lwd=4,type = "l",col=colVec[i])
	points((localEICPeakList[i,]$pkInd-1)*ScanInterval+delaytime,localEICPeakList[i,]$Intensity,col="black",lwd=3)
}

# plot 2-4
plot(sort(localEICPeakList$gss),ylab="Gss")
abline(h=2.5,col="red")
legend("topleft",paste("Number of EIC Peaks (Gss <=2.5):", as.character(length(which(localEICPeakList$gss<=2.5))),sep=" "))

plot(sort(localEICPeakList$StN),ylab="StN")
abline(h=500,col="red")
legend("topleft",paste("Number of EIC Peaks (StN >= 500):", as.character(length(which(localEICPeakList$StN>=500))),sep=" "))

plot(sort(localEICPeakList$shrp),ylab="Sharpness")
abline(h=5,col="red")
legend("topleft",paste("Number of EIC Peaks (shrp >= 5):", as.character(length(which(localEICPeakList$StN>=5))),sep=" "))

# plot 5: candidate peak profiles
minLbound<-min(CandidatePeaklist$lboundInd)
maxRbound<-max(CandidatePeaklist$rboundInd)	
maxInt2<-max(CandidatePeaklist$Intensity)
colVec<-rainbow(nrow(CandidatePeaklist))
masses <-CandidatePeaklist$mz

plot(x=0,type = "n",xlim=((c(minLbound,maxRbound)-1)*ScanInterval+delaytime),ylim=c(0,max(maxInt2)),
		main="candidate Peaks",xlab="ET",ylab="Intensity")
for(i in c(1:nrow(CandidatePeaklist)))	
{
	
	curID<-rownames(CandidatePeaklist[i,])
	points((allprofiles[[curID]]$ET-1)*ScanInterval+delaytime,allprofiles[[curID]]$int, lwd=4,type = "l",col=colVec[i])
	points((CandidatePeaklist[i,]$pkInd-1)*ScanInterval+delaytime,CandidatePeaklist[i,]$Intensity,col="black",lwd=3)
}
legend("topright",as.character(masses),fill=colVec)

# plot 6: good peak profiles
if (nPks>0) {
	minLbound1<-min(goodShapePeaklist$lboundInd)
	maxRbound1<-max(goodShapePeaklist$rboundInd)	
	maxInt1<-max(goodShapePeaklist$Intensity)
	colVec<-rainbow(nrow(goodShapePeaklist))
	masses <-goodShapePeaklist$mz
	
	plot(x=0,type = "n",xlim=((c(minLbound1,maxRbound1)-1)*ScanInterval+delaytime),ylim=c(0,max(maxInt1)),
			main="good Peaks",xlab="ET",ylab="Intensity")
	for(i in c(1:nrow(goodShapePeaklist)))	
	{
		
		curID<-rownames(goodShapePeaklist[i,])
		points((allprofiles[[curID]]$ET-1)*ScanInterval+delaytime,allprofiles[[curID]]$int, lwd=4,type = "l",col=colVec[i])
		points((goodShapePeaklist[i,]$pkInd-1)*ScanInterval+delaytime,goodShapePeaklist[i,]$Intensity,col="black",lwd=3)
	}
	legend("topright",as.character(masses),fill=colVec)
}


# plot 7: good peak profiles
#x11()
minLbound2<-min(goodShapePeaklist$lboundInd)
maxRbound2<-max(goodShapePeaklist$rboundInd)	
maxInt2 <- max(goodShapePeaklist$Intensity)+500

colVec<-rainbow(nrow(goodShapePeaklist))
masses <-goodShapePeaklist$mz

plot(x=0,type = "n",xlim=((c(minLbound2,maxRbound2)-1)*ScanInterval+delaytime),ylim=c(0,max(maxInt2)),main="mirrored good peaks",xlab="ET",ylab="Intensity")
for(i in 1:nrow(goodShapePeaklist))
{
	curID<-rownames(goodShapePeaklist[i,])
	points((profiles[[curID]]$ET-1)*ScanInterval+delaytime,profiles[[curID]]$int, lwd=3,type = "l",col=colVec[i])
}
legend("topright",as.character(masses),fill=colVec)

# Plotting =================================================================================================================================
# Ploting model peak features
playwith({
			
			plot((c(lbound:rbound)-1)*ScanInterval+delaytime,S[,1],type = "l",lwd=4)
			for (i in 2:ncol(S)) {
				points((c(lbound:rbound)-1)*ScanInterval+delaytime,S[,i],type = "l",col=i+1,lwd=4)
			}
		})

#plot the profiles of model peaks
minLbound<-min(modelPkList$lboundInd)
maxRbound<-max(modelPkList$rboundInd)
maxInt2 <- max(modelPkList$Intensity)+500
colVec<-rainbow(nrow(modelPkList))

labelVec<-paste("pkid:",rownames(modelPkList),"mz:",modelPkList$mz," ,gss:",format(modelPkList$gss,digits=2)," ,shrp:",format(modelPkList$shrp,digits=2))
#					labelVec<-paste("gss=",format(modelPkList$gss,digits=2),"\n shrp=",format(modelPkList$shrp,digits=2))
dataPoints<-subset(modelPkList,select=c("pkInd","Intensity"))
dataPoints$pkInd<-(dataPoints$pkInd-1)*ScanInterval+delaytime

playwith({
			plot(x=0,type = "n",xlim=(c(minLbound,maxRbound)-1)*ScanInterval+delaytime,ylim=c(0,max(maxInt2)),main="model peaks of four clusters",xlab="ET",ylab="Intensity")
			for(i in 1:nrow(modelPkList))
			{
				curPkid<-rownames(modelPkList)[i]			
				points((profiles[[curPkid]]$ET-1)*ScanInterval+delaytime,profiles[[curPkid]]$int,lwd=6, type = "l",col=colVec[i])
			}
			legend("topright",as.character(modelPkList$mz),fill=colVec)
		},data.points=dataPoints,labels=labelVec)

# Clustering
pdf(paste("/Users/yni1/Desktop/","clustering.pdf",sep=""))	
tempProfiles <- profiles
names(tempProfiles) <- goodShapePeaklist$mz	
if(nPks>=2)
{
	#save the EIC profies in global variable for broadcasting to slave nodes for parallel computing
	assign("profiles", value=tempProfiles, envir = .GlobalEnv)
	r<-parDistCal4(cl,isUnion=F)##parallel verison
	distance <- as.dist(r)
	clustResut<-hclust(distance)
	if (length(clustResut$order)>2)
	{
		x11()
		plot(clustResut,lwd=2,cex=1.5,font.axis=2,ylim=c(0,90))
		axis(2,lwd=2)
		abline(h=100,lty=2,col="red",lwd=2)
		#title("Hierachincal clustering",xlab="mass",ylab="distance")
	}
	
}
dev.off()


x11()
par(mfrow=c(2,1))
for (i in c(1:2))
{
	plot(specList[[i]],type="h",col=i,main=paste(names(specList)[i],"(mins, mass)",sep= " "),lwd=2.5,xlim=c(min(specList[[i]]$mz),400),ylim=c(0,0.3+max(specList[[i]]$int)),ylab="relative intensity",xlab="mz")
	text(specList[[i]]$mz[order(specList[[i]]$int,decreasing=TRUE)[1:10]],specList[[i]]$int[order(specList[[i]]$int,decreasing=TRUE)[1:10]],specList[[i]]$mz[order(specList[[i]]$int,decreasing=TRUE)[1:10]],pos=3,cex=1.2,col="blue")
	axis(2,lwd=2)
	axis(1,lwd=2)
}

# plot shared peak features
x11()
plot(X[,10],type="l",ylim=c(1,10000))
for (i in c(2:70)) {
	points(X[,i],col=i,type="l")
}

