# Plotting =========================================== PK picking ==============================================================	
playwith ({
			plot(c(1:length(curVecInt)),curVecInt,lwd=2,type = "l",main=paste("mass",mz))
			points(peakInd,curVecInt[peakInd],col="yellow")
			points(curPeakList$pkInd,curVecInt[curPeakList$pkInd],col=curPeakList$isShared+1)

			points(valleyInd,curVecInt[valleyInd],col="green")
			#points(valleyInd_LM,curVecInt[valleyInd_LM],col="pink")
			text(curPeakList$pkInd,curVecInt[curPeakList$pkInd],as.character(curPeakList$StN))
		},new=TRUE)
#

playwith ({
			plot(c(1:length(curVecInt)),curVecInt,lwd=2,type = "l",main="TIC")
			points(peakInd,curVecInt[peakInd],col="yellow")
			points(curPeakList$pkInd,curVecInt[curPeakList$pkInd],col=curPeakList$isShared+1)
			
			points(valleyInd,curVecInt[valleyInd],col="green")
			points(valleyInd_LM,curVecInt[valleyInd_LM],col="pink")
			text(curPeakList$pkInd,curVecInt[curPeakList$pkInd],as.character(curPeakList$StN))
		},new=TRUE)


playwith ({
			plot((c(1:length(curVecInt))-1)*ScanInterval+delaytime,curVecInt,lwd=2,type = "l")
			points((curPeakList$pkInd-1)*ScanInterval+delaytime,curVecInt[curPeakList$pkInd],col="red")
			points((peakInd-1)*ScanInterval+delaytime,curVecInt[peakInd],col="blue")
			points((curPeakList$pkInd-1)*ScanInterval+delaytime,curVecInt[curPeakList$pkInd],col="red")
			points((valleyInd-1)*ScanInterval+delaytime,curVecInt[valleyInd],col="green")
			text((curPeakList$pkInd-1)*ScanInterval+delaytime,curVecInt[curPeakList$pkInd],as.character(curPeakList$StN))
		},new=TRUE)



playwith ({
			plot(c(1:length(curVecInt)),curVecInt,lwd=2,type = "l")
			points(curPeakList$pkInd,curVecInt[curPeakList$pkInd],col="red")
			points(curPeakList$Lbound,curVecInt[curPeakList$Lbound],col="blue")
			points(curPeakList$Rbound,curVecInt[curPeakList$Rbound],col="green")
			
			abline(v=curPeakList$lboundInd,col="grey")
			abline(v=curPeakList$rboundInd,col="grey")
			#text(curPeakList$pkInd,curVecInt[curPeakList$pkInd],as.character(curPeakList$shrp_Mean),col=curPeakList$isShared+1)
			
			
			#points(WT_result_min$x,curVecInt[WT_result_min$x],lwd=2,col="green")
			text(curPeakList$pkInd,curVecInt[curPeakList$pkInd],as.character(curPeakList$StN))
			#text(curPeakList$pkInd,curVecInt[curPeakList$pkInd],c(1:nrow(curPeakList)))
			
			#text(curPeakList$pkInd,curVecInt[curPeakList$pkInd],paste(curPeakList$shrp_old,curPeakList$shrp_Mean,curPeakList$shrp_Med,curPeakList$shrp_Max))
#					text(curPeakList$pkInd,curVecInt[curPeakList$pkInd],as.character(round(curPeakList$shrp_old)))
#					text(curPeakList$pkInd,curVecInt[curPeakList$pkInd]+10000,as.character(round(curPeakList$shrp_Mean)))
#					text(curPeakList$pkInd,curVecInt[curPeakList$pkInd]+20000,as.character(round(curPeakList$shrp_Med)))
#					text(curPeakList$pkInd,curVecInt[curPeakList$pkInd]+30000,as.character(round(curPeakList$shrp_Max)))
			
			#text(curPeakList$pkInd,curVecInt[curPeakList$pkInd],paste(curPeakList$shrp_Mean1,curPeakList$shrp_Mean2),col=curPeakList$isShared+1)
			
			#text(curPeakList$pkInd,curVecInt[curPeakList$pkInd],as.character(curPeakList$zigNum))
			
		},new=TRUE)

# Plotting =========================================== Deconvolution ==============================================================	
# ploting decon window
playwith({
			
			plot((1:length(TIC)-1)*ScanInterval+delaytime,TIC,type="l",main=fileName,xlab="ET(Mins)")
			points(TICpeaklist$RT,TIC[TICpeaklist$pkInd],col="red")
			text(TICpeaklist$RT,TIC[TICpeaklist$pkInd],TICpeaklist$window)
			abline(v=(TICpeaklist$lboundInd-1)*ScanInterval+delaytime,col="blue")
			abline(v=(TICpeaklist$rboundInd-1)*ScanInterval+delaytime,col="blue")
			
		},new=TRUE)	

# ==================================================	
# Plot clustering results - noisy peaks checking
x11()		
pkInd <- localEICPeakList$pkInd
names(pkInd) <- localEICPeakList$mz
#clustResult1 <- hclust(dist(pkInd))
plot(clustResult1,lwd=2,cex=1.5,font.axis=2,ylim=c(0,90))
axis(2,lwd=2)
abline(h=60,lty=2,col="red",lwd=2)

# ==================================================	

# Plot all EIC and TIC profiles
allprofiles<-vector("list",nrow(localEICPeakList))
for(i in 1:nrow(localEICPeakList))
{
	##get intensity vector of current mz
	##from long intensity vector
	curEICpk<-localEICPeakList[i,]			
	#startInd<-curEICpk$offset+curEICpk$lboundInd
	#endInd<-curEICpk$offset+curEICpk$rboundInd
	
	startInd<-curEICpk$offset+curEICpk$Lbound
	endInd<-curEICpk$offset+curEICpk$Rbound
	
	##get EIC peak profile
	ET<-curEICpk$Lbound:curEICpk$Rbound
	curProfile<-data.frame(ET=ET,int=vecInt[startInd:endInd]) 
	allprofiles[[i]]<-curProfile
	names(allprofiles)[i]<-rownames(curEICpk)	
	
}	

# ==================================================	
localEICPeakList <- subset(localEICPeakList,Intensity>=500)
minLbound<-min(localEICPeakList$Lbound)
maxRbound<-max(localEICPeakList$Rbound)	
maxInt<-max(TIC[curTICPk$lboundInd:curTICPk$rboundInd])
colVec<-rainbow(nrow(localEICPeakList))
colVec<-localEICPeakList$clust#rainbow(nrow(localEICPeakList))#localEICPeakList$isShared+1##

masses <-localEICPeakList$mz
localTICApexList<-subset(TICpeaklist,pkInd>=curTICPk$lboundInd&pkInd<=curTICPk$rboundInd)
playwith({
			plot(x=0,type = "n",xlim=((c(minLbound,maxRbound)-1)*ScanInterval+delaytime),ylim=c(0,maxInt),main=paste(fileName,"TIC and local EICs"),xlab="ET",ylab="Intensity")
			
			for(i in c(1:nrow(localEICPeakList)))	
			{
				curID<-rownames(localEICPeakList[i,])
				points((allprofiles[[curID]]$ET-1)*ScanInterval+delaytime,allprofiles[[curID]]$int, lwd=4,type = "l",col=colVec[i])		
			}
			points((c(curTICPk$lboundInd:curTICPk$rboundInd)-1)*ScanInterval+delaytime,TIC[curTICPk$lboundInd:curTICPk$rboundInd],col="grey",lwd=4,type="l")
			points((localTICApexList$pkInd-1)*ScanInterval+delaytime,TIC[localTICApexList$pkInd],col="blue")	
			#legend("topright",as.character(masses),fill=colVec)
		},new=TRUE)


# ==================================================	
# Plot clustering results 
#x11()
clustResut$labels <- uniqueEICPeakList$mz
playwith ({
			plot(clustResut,lwd=2,cex=1.5,font.axis=2,ylim=c(0,90))
		},new=TRUE)
plot(clustResut,lwd=2,cex=1.5,font.axis=2,ylim=c(0,90))
#axis(2,lwd=2)
abline(h=15,lty=2,col="red",lwd=2)

# ==================================================	
# Plot EIC of interest
EICPeaks <- localEICPeakList

EICPeaks <- uniqueEICPeakList#localEICPeakList#[curClusterIDs,]#noisyEICPeaks#modelPkList#
colVec<-EICPeaks$clust

EICPeaks <- uniqueEICPeakList
colVec<-EICPeaks$CompoundID
EICPeaks <- sharedEICPeakList
colVec<-rainbow(nrow(EICPeaks))
EICPeaks <-curUniqueEICPeakList

EICPeaks <- uniqueEICPeakList
colVec<-rainbow(nrow(EICPeaks))

minLbound<-min(EICPeaks$lboundInd)
maxRbound<-max(EICPeaks$rboundInd)	
maxInt<-max(EICPeaks$Intensity)
masses <-EICPeaks$mz

playwith({
			plot(x=0,type = "n",xlim=((c(minLbound,maxRbound)-1)*ScanInterval+delaytime),ylim=c(0,maxInt),main=paste(fileName,"EICs of interest"),xlab="ET",ylab="Intensity")
			
			#IDs <- which(curClusterPeakList$clust==4)
			for(i in c(1:nrow(EICPeaks)))	
			#for(i in IDs)	
			{				
				curID<-rownames(EICPeaks[i,])
				curMz <- EICPeaks[i,]$mz
				#points((allprofiles[[i]]$ET-1)*ScanInterval+delaytime,allprofiles[[i]]$int, lwd=4,type = "l",col=colVec[i])
				points((allprofiles[[curID]]$ET-1)*ScanInterval+delaytime,allprofiles[[curID]]$int, lwd=4,type = "l",col=colVec[i])		
				points((EICPeaks[i,]$pkInd-1)*ScanInterval+delaytime,EICPeaks[i,]$Intensity,col="red",lwd=3)
				#text((EICPeaks[i,]$pkInd-1)*ScanInterval+delaytime,EICPeaks[i,]$Intensity,as.character(EICPeaks[i,]$shrp))	
				#text((EICPeaks[i,]$pkInd-1)*ScanInterval+delaytime,EICPeaks[i,]$Intensity,as.character(EICPeaks[i,]$mz))	
				#legend("topright",as.character(masses),fill=colVec)
				abline(v=(curCompPkInd-1)*ScanInterval+delaytime,lwd=3,lty=2,col="blue")
			}	
		},new=TRUE)

# ==================================================	

#								playwith({
#											plot(curEICProfile$int,type="l")
#											points(10000*curEICProfile$int/max(curEICProfile$int),,type="l",col="red")
#										})
#								
#								x <- curEICProfile$int
#								y <- 1000*curEICProfile$int/max(curEICProfile$int)
#								
#								(x%*%y)/((x%*%x)^0.5*(y%*%y)^0.5)

EICPeaks <- modelPkList#[curClusterIDs,]#noisyEICPeaks#modelPkList#
minLbound<-min(EICPeaks$Lbound)
maxRbound<-max(EICPeaks$Rbound)	
maxInt<-max(EICPeaks$Intensity)
colVec<-rainbow(nrow(EICPeaks))
masses <-EICPeaks$mz

playwith({
			#plot(x=0,type = "n",xlim=((c(minLbound,maxRbound)-1)*ScanInterval+delaytime),ylim=c(0,maxInt),main="Selected model peaks",xlab="ET",ylab="Intensity")
			plot(x=0,type = "n",xlim=c(21.8,22.1),ylim=c(0,maxInt),main="Selected model peaks",xlab="ET",ylab="Intensity")
			
			#IDs <- which(curClusterPeakList$clust==4)
			for(i in c(1:nrow(EICPeaks)))	
			{				
				curEICpk <- EICPeaks[i,]
				curMz <- curEICpk$mz
				
				startInd<-curEICpk$offset+curEICpk$Lbound
				endInd<-curEICpk$offset+curEICpk$Rbound
				ET<-curEICpk$Lbound:curEICpk$Rbound
				curProfile<-data.frame(ET=ET,int=vecInt[startInd:endInd]) 
				
				
				points((curProfile$ET-1)*ScanInterval+delaytime,curProfile$int, lwd=4,type = "l",col=colVec[i])
				points((curEICpk$pkInd-1)*ScanInterval+delaytime,curEICpk$Intensity,col="red",lwd=3)
				#text((EICPeaks[i,]$pkInd-1)*ScanInterval+delaytime,EICPeaks[i,]$Intensity,as.character(EICPeaks[i,]$shrp))	
				text((curEICpk$pkInd-1)*ScanInterval+delaytime,curEICpk$Intensity,as.character(curEICpk$mz))	
			}	
		},new=TRUE)

# ==================================================	
# plotting EICs and scaled EICs for comparison

#Profiles <- allprofiles[which(rownames(localEICPeakList)%in%rownames(curUniqueEICPeakList))]	
EICs <- c(117,129)
Profiles <- allprofiles[which(rownames(localEICPeakList)%in%rownames(uniqueEICPeakList[which(uniqueEICPeakList$mz%in%EICs),]))]									
minLbound <- (curTICPk$lboundInd-1)*ScanInterval+delaytime
maxRbound <- (curTICPk$rboundInd-1)*ScanInterval+delaytime
maxInt <- max(uniqueEICPeakList$Intensity)
colVec <- rainbow(length(Profiles))

playwith ({
			plot(x=0,type = "n",xlim=c(minLbound,maxRbound),ylim=c(0,maxInt),main=paste(fileName,"unique peak features"),xlab="ET",ylab="Intensity")
			
			for(i in c(1:length(Profiles)))	
			{
				curID <- names(Profiles)[i]
				points((Profiles[[i]]$ET-1)*ScanInterval+delaytime, Profiles[[i]]$int, lwd=4,type = "l",col=colVec[i])		
				
				text((uniqueEICPeakList[curID,]$pkInd-1)*ScanInterval+delaytime,uniqueEICPeakList[curID,]$Intensity,
						paste(uniqueEICPeakList[curID,]$mz,round(1000*uniqueEICPeakList[curID,]$shrp/uniqueEICPeakList[curID,]$Intensity),sep="_"))
			}
			
			#legend("topright",as.character(paste(curUniqueEICPeakList$mz,round(1000*curUniqueEICPeakList$shrp/curUniqueEICPeakList$Intensity),sep="_")),fill=colVec)
		},new=TRUE)

playwith ({
			plot(x=0,type = "n",xlim=c(minLbound,maxRbound),ylim=c(0,10000),main=paste(fileName,"unique peak features"),xlab="ET",ylab="Intensity")
			
			for(i in c(1:length(Profiles)))	
			{
				curID <- names(Profiles)[i]	
				points((Profiles[[i]]$ET-1)*ScanInterval+delaytime, 10000*Profiles[[i]]$int/uniqueEICPeakList[curID,]$Intensity, lwd=4,type = "l",col=colVec[i])									
			}
			#legend("topright",as.character(paste(curUniqueEICPeakList$mz,round(1000*curUniqueEICPeakList$shrp/curUniqueEICPeakList$Intensity),sep="_")),fill=colVec)
		},new=TRUE)

# ==================================================
# within the clustering loop

Profiles <- allprofiles[which(rownames(localEICPeakList)%in%rownames(curUniqueEICPeakList))]	
#Profiles <- allprofiles[which(rownames(localEICPeakList)%in%rownames(curUniqueEICPeakList[which(curUniqueEICPeakList$mz%in%c(73,60)),]))]									
minLbound <- (curTICPk$lboundInd-1)*ScanInterval+delaytime
maxRbound <- (curTICPk$rboundInd-1)*ScanInterval+delaytime
maxInt <- max(curUniqueEICPeakList$Intensity)
colVec <- rainbow(length(Profiles))

playwith ({
			plot(x=0,type = "n",xlim=c(minLbound,maxRbound),ylim=c(0,maxInt),main=paste(fileName,"unique peak features"),xlab="ET",ylab="Intensity")
			
			for(i in c(1:length(Profiles)))	
			{
				curID <- names(Profiles)[i]
				points((Profiles[[i]]$ET-1)*ScanInterval+delaytime, Profiles[[i]]$int, lwd=4,type = "l",col=colVec[i])		
				
				text((curUniqueEICPeakList[curID,]$pkInd-1)*ScanInterval+delaytime,curUniqueEICPeakList[curID,]$Intensity,
						paste(curUniqueEICPeakList[curID,]$mz,round(1000*curUniqueEICPeakList[curID,]$shrp/curUniqueEICPeakList[curID,]$Intensity),sep="_"))
			}
			
			#legend("topright",as.character(paste(curUniqueEICPeakList$mz,round(1000*curUniqueEICPeakList$shrp/curUniqueEICPeakList$Intensity),sep="_")),fill=colVec)
		},new=TRUE)

playwith ({
			plot(x=0,type = "n",xlim=c(minLbound,maxRbound),ylim=c(0,10000),main=paste(fileName,"unique peak features"),xlab="ET",ylab="Intensity")
			
			for(i in c(1:length(Profiles)))	
			{
				curID <- names(Profiles)[i]	
				points((Profiles[[i]]$ET-1)*ScanInterval+delaytime, 10000*Profiles[[i]]$int/curUniqueEICPeakList[curID,]$Intensity, lwd=4,type = "l",col=colVec[i])									
			}
			#legend("topright",as.character(paste(curUniqueEICPeakList$mz,round(1000*curUniqueEICPeakList$shrp/curUniqueEICPeakList$Intensity),sep="_")),fill=colVec)
		},new=TRUE)
# ==================================================

# Plot unique peak mass feature
EICPeaks <- qualify_uniqueEICPeakList
#EICPeaks <- sharedEICPeakList
#EICPeaks <- subset(uniqueEICPeakList,Intensity>=500)
minLbound<-min(EICPeaks$lboundInd)
maxRbound<-max(EICPeaks$rboundInd)	
maxInt<-max(EICPeaks$Intensity)
#colVec<-rainbow(nrow(EICPeaks)) #colVec<-HC_ClusterResult
colVec<- FinalClustResut
masses <-EICPeaks$mz
playwith({
			plot(x=0,type = "n",xlim=((c(minLbound,maxRbound)-1)*ScanInterval+delaytime),ylim=c(0,maxInt),main=paste(fileName,"EICs of interest"),xlab="ET",ylab="Intensity")
			
			for(i in c(1:nrow(EICPeaks)))	
			{
				curID<-rownames(EICPeaks[i,])
				curMz <- EICPeaks[i,]$mz
				points((allprofiles[[curID]]$ET-1)*ScanInterval+delaytime,allprofiles[[curID]]$int, lwd=4,type = "l",col=colVec[i])
				points((EICPeaks[i,]$pkInd-1)*ScanInterval+delaytime,EICPeaks[i,]$Intensity,col="red",lwd=3)
				text((EICPeaks[i,]$pkInd-1)*ScanInterval+delaytime,EICPeaks[i,]$Intensity,as.character(curMz))	
				#text((EICPeaks[i,]$pkInd-1)*ScanInterval+delaytime,EICPeaks[i,]$Intensity,as.character(EICPeaks[i,]$StN))	
			}	
			#legend("topright",as.character(masses),fill=colVec)
		},new=TRUE)

# ==================================================	
# plot normalized unique mass peak
maxEIC_pos <- which.max(EICPeaks$Intensity)
maxEIC_pkInd <- EICPeaks[maxEIC_pos,]$pkInd
playwith({
			plot(x=0,type = "n",xlim=((c(minLbound,maxRbound)-1)*ScanInterval+delaytime),ylim=c(0,1200),main=paste(fileName,"EICs of interest"),xlab="ET",ylab="Intensity")
			
			for(i in c(1:nrow(EICPeaks)))	
			{
				curID<-rownames(EICPeaks[i,])
				curMz <- EICPeaks[i,]$mz
				curPkInd <-  EICPeaks[i,]$pkInd
				curPkInt <- EICPeaks[i,]$Intensity
				pkInd_diff <- maxEIC_pkInd-curPkInd
				points(((allprofiles[[curID]]$ET+pkInd_diff)-1)*ScanInterval+delaytime,(allprofiles[[curID]]$int/curPkInt)*1000, lwd=1,type = "l",col=colVec[i])			
				#points((allprofiles[[curID]]$ET-1)*ScanInterval+delaytime,(allprofiles[[curID]]$int/curPkInt)*1000, lwd=1,type = "l",col=colVec[i])
				#points((EICPeaks[i,]$pkInd-1)*ScanInterval+delaytime,(EICPeaks[i,]$Intensity/curPkInt)*1000,col="blue")
				#points((EICPeaks[i,]$pkInd-1)*ScanInterval+delaytime,EICPeaks[i,]$Intensity,col="red",lwd=3)
				#text((EICPeaks[i,]$pkInd-1)*ScanInterval+delaytime,(EICPeaks[i,]$Intensity/curPkInt)*1000,as.character(curMz))
				#legend("topright",as.character(masses),fill=colVec)
			}
	},new=TRUE)


lowShrpEIC <- subset(EICPeaks,shrp<=10)
playwith({
			plot(x=0,type = "n",xlim=((c(minLbound,maxRbound)-1)*ScanInterval+delaytime),ylim=c(0,1200),main=paste(fileName,"EICs of interest"),xlab="ET",ylab="Intensity")
			
			for(i in c(1:nrow(EICPeaks)))	
			{
				curID<-rownames(EICPeaks[i,])
				curMz <- EICPeaks[i,]$mz
				curPkInd <-  EICPeaks[i,]$pkInd
				curPkInt <- EICPeaks[i,]$Intensity
				pkInd_diff <- maxEIC_pkInd-curPkInd
				
				if (curID%in%rownames(lowShrpEIC)==FALSE) {
					points((allprofiles[[curID]]$ET-1)*ScanInterval+delaytime,(allprofiles[[curID]]$int/curPkInt)*1000, lwd=1,type = "l",col="black")
					points((EICPeaks[i,]$pkInd-1)*ScanInterval+delaytime,(EICPeaks[i,]$Intensity/curPkInt)*1000,col="blue")
				}else {
					points((allprofiles[[curID]]$ET-1)*ScanInterval+delaytime,(allprofiles[[curID]]$int/curPkInt)*1000, lwd=1,type = "l",col=colVec[i])
					points((EICPeaks[i,]$pkInd-1)*ScanInterval+delaytime,(EICPeaks[i,]$Intensity/curPkInt)*1000,col="blue")
				}
				
				#points(((allprofiles[[curID]]$ET+pkInd_diff)-1)*ScanInterval+delaytime,(allprofiles[[curID]]$int/curPkInt)*1000, lwd=1,type = "l",col=colVec[i])			

				#points((EICPeaks[i,]$pkInd-1)*ScanInterval+delaytime,EICPeaks[i,]$Intensity,col="red",lwd=3)
				#text((EICPeaks[i,]$pkInd-1)*ScanInterval+delaytime,(EICPeaks[i,]$Intensity/curPkInt)*1000,as.character(curMz))
				#legend("topright",as.character(masses),fill=colVec)
			}
		},new=TRUE)
# ==================================================	

for (i in 1:length(specList)) {
	specList[[i]]$int<- 999*specList[[i]]$int/max(specList[[i]]$int)
}
x11()
par(mfrow=c(length(specList),1))
#par(mfrow=c(2,1))
for (i in c(1:length(specList)))
{
	plot(specList[[i]],type="h",col=i,main=names(specList)[i],lwd=2.5,xlim=c(min(specList[[i]]$mz),max(specList[[i]]$mz)),ylim=c(0,1200),ylab="relative intensity",xlab="mz")
	text(specList[[i]]$mz[order(specList[[i]]$int,decreasing=TRUE)[1:10]],specList[[i]]$int[order(specList[[i]]$int,decreasing=TRUE)[1:10]],specList[[i]]$mz[order(specList[[i]]$int,decreasing=TRUE)[1:10]],pos=3,cex=1.2,col="blue")
	axis(2,lwd=2)
	axis(1,lwd=2)
}

x11()
par(mfrow=c(4,2))
curProfiles <- allprofiles[which(rownames(localEICPeakList)%in%rownames(uniqueEICPeakList))]			
for (i in 1:length(curProfiles)) {
	x<- curProfiles[[i]]$int
	y <- filter(x, c(0.2,0.2,0.2,0.2,0.2))

	plot(x,type="l")
	lines(y,col="Red")
	
	x[!is.na(y)] <- y[!is.na(y)]
	curProfiles[[i]]$int <- x
}


curProfiles <- allprofiles[which(rownames(localEICPeakList)%in%rownames(uniqueEICPeakList))]			
for (i in 1:length(curProfiles)) {
	y<- curProfiles[[i]]$int
	x <- c(1:length(y))
	fitresult <- bezierCurve(x,y,length(y))
	curProfiles[[i]]$int <- fitresult$y
	
#	plot(x, y,main=paste("beziercurve function-mass",curmz))
#	points(bezierCurve(x,y,length(y)), type="l", col="red",lwd=2)
#	
#	plot(x,type="l")
#	lines(y,col="Red")
#	
#	x[!is.na(y)] <- y[!is.na(y)]
#	curProfiles[[i]]$int <- x
}


library(nls)
y<- curProfiles[[i]]$int
x <- c(1:length(y))
nlmod <- nls(y~Const+A*x*x)

plot(x,y, main = "nls(*), data, true function and fit, n=100")
lines(x, predict(nlmod), col = 2)



for (i in 1:length(curProfiles)) {
	
	x11()
	par(mfrow=c(2,2))
	curmz <- uniqueEICPeakList$mz[i]
	y<- curProfiles[[i]]$int
	x <- c(1:length(y))
	
	lo <- loess(y~x)
	plot(x, y, main = paste("loess function-mass ",curmz))
	lines(x, predict(lo), col = 2,lwd=2)
	
	smoothingSpline = smooth.spline(x, y, spar=0.35)
	plot(x, y, main=paste("smoothspine function-mass ",curmz))
	lines(smoothingSpline,col="red",lwd=2)
	
	plot(x, y,main=paste("beziercurve function-mass",curmz))
	points(bezierCurve(x,y,length(y)), type="l", col="red",lwd=2)
	
	x<- curProfiles[[i]]$int
	y <- filter(x, c(0.2,0.2,0.2,0.2,0.2))
	plot(x,main=paste("filter function-mass",curmz))
	lines(y,col="Red",lwd=2)
}

# x, y: the x and y coordinates of the hull points
# n: the number of points in the curve.
bezierCurve <- function(x, y, n=10)
{
	outx <- NULL
	outy <- NULL
	
	i <- 1
	for (t in seq(0, 1, length.out=n))
	{
		b <- bez(x, y, t)
		outx[i] <- b$x
		outy[i] <- b$y
		
		i <- i+1
	}
	
	return (list(x=outx, y=outy))
}

bez <- function(x, y, t)
{
	outx <- 0
	outy <- 0
	n <- length(x)-1
	for (i in 0:n)
	{
		outx <- outx + choose(n, i)*((1-t)^(n-i))*t^i*x[i+1]
		outy <- outy + choose(n, i)*((1-t)^(n-i))*t^i*y[i+1]
	}
	
	return (list(x=outx, y=outy))
}

# Example usage
plot(x, y, "o", pch=20)
points(bezierCurve(x,y,20), type="l", col="red")

TIC <- NULL
filenames <- NULL
for(fileindex in c(1:length(DataFilelist)))
{
	
	inFilePath <- DataFilelist[[fileindex]]
	fileName<-parseFileName(inFilePath)
	filenames <- c(filenames,fileName)
	TICfile<-paste(WorkDir,"output/TIC/denoised_",fileName,"_TIC.cdf",sep="") 
	TICfile<-ifelse(file.exists(TICfile),TICfile,inFilePath)
	ncid <- nc_open(TICfile)
	TIC[[fileindex]] <- ncvar_get(ncid, varid="total_intensity" )
	nc_close(ncid)
	remove(ncid)
}

colVec <- rainbow(7)
playwith({
			
			plot((1:length(TIC[[1]])-1)*ScanInterval+delaytime,TIC[[1]],type="l",,col=colVec[1],xlab="ET(Mins)",lwd=4)
			for (i in 2:7) {
				points((1:length(TIC[[i]])-1)*ScanInterval+delaytime,TIC[[i]],col=colVec[i],type="l",lwd=4)
			}
			#legend("topleft",as.character(filenames),fill=colVec)
		},new=TRUE)	

points(TICpeaklist$RT,TIC[TICpeaklist$pkInd],col="red")
text(TICpeaklist$RT,TIC[TICpeaklist$pkInd],TICpeaklist$window)
abline(v=(TICpeaklist$lboundInd-1)*ScanInterval+delaytime,col="blue")
abline(v=(TICpeaklist$rboundInd-1)*ScanInterval+delaytime,col="blue")
