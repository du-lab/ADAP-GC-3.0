#====================COMMANDS========================
# for old laptop
library(playwith)
codeDir <<-"/Users/yanni/Desktop/yni1/work/SourceCode/ADAP3"
source(paste(codeDir,"pipeline.r",sep="/"))

params<-readParaFromCmd("/Users/yanni/DataTest/JiaLab/stds/Stds.db")
#params<-readParaFromCmd("/Users/yanni/DataTest/JiaLab/U50Stds2013/U50Stds2013.db")
params$codeDir<-codeDir
fileInfo<-read.csv(paste(params$WorkDir,"output/",params$JobName,"_RT.csv",sep=""))
params$delaytime<-fileInfo$firstRT[1]
params$ScanInterval<-fileInfo$ScanInterval[1]

# remote aether commands
codeDir <<-"code"
source(paste(codeDir,"pipeline.r",sep="/"))
#params<-readParaFromCmd("/home/yni/DataTest/U50Stds2013/U50Stds2013.db")
params<-readParaFromCmd("/home/yni/DataTest/LI102912/LI102912.db")
params$codeDir<-codeDir
parTICPeakpicking(params,denoised=F)
parEICpeakpicking(params)

codeDir<-"code"
source(paste(codeDir,"pipeline.r",sep="/"))
params<-readParaFromCmd("/home/yni/DataTest/U50Stds2013/U50Stds2013.db")
params$codeDir<-codeDir
fileInfo<-read.csv(paste(params$WorkDir,"output/",params$JobName,"_RT.csv",sep=""))
params$delaytime<-fileInfo$firstRT[1]
params$ScanInterval<-fileInfo$ScanInterval[1]
deconvolution(params,isDistParallel=F,clustingType="h",isDataParalell=T)

#====================COMMANDS========================
# for new lap
library(playwith)
codeDir <<-"/Users/yni/Desktop/yni1/work/SourceCode/ADAP3"
source(paste(codeDir,"pipeline.r",sep="/"))

params<-readParaFromCmd("/Users/yni/DataTest/stds/Stds.db")
#params<-readParaFromCmd("/Users/yni/DataTest/U50Stds2013/U50Stds2013.db")
params$codeDir<-codeDir
fileInfo<-read.csv(paste(params$WorkDir,"output/",params$JobName,"_RT.csv",sep=""))
params$delaytime<-fileInfo$firstRT[1]
params$ScanInterval<-fileInfo$ScanInterval[1]

#====================COMMANDS========================
# peak picking detection
fileindex = 1
mz=154

fileindex=7
inFilePath <- DataFilelist[[fileindex]]
fileName<-parseFileName(inFilePath)

#check if it exist denoised EIC file
rawEICFile <-paste(WorkDir,"output/EIC/",fileName,"EIC.cdf",sep="")
denoiseEICFile<-paste(WorkDir,"output/EIC/denoised_",fileName,"EIC.cdf",sep="")
EICFile<-ifelse(file.exists(denoiseEICFile),denoiseEICFile,rawEICFile)

if(file.exists(inFilePath)==FALSE)next
ncid <- nc_open(EICFile)
vecInt <- ncvar_get(ncid, varid="intVec")
mzVec<-ncvar_get(ncid, varid="mzVec")

nc_close(ncid)
remove(ncid)

cat(fileName,"start peak picking...\n")
mzGroup <- mzVec
time1<-Sys.time()
PeakList <- getPeaksGroup(mzGroup,vecInt,mzVec,params,fileName)
PeakList<-do.call(rbind,PeakList)
write.csv(PeakList,file=paste(WorkDir,"/output/peakpicking/",fileName,"_EIC_PeakList.csv",sep=""),row.names=F)	
time2<-Sys.time()
time2-time1
cat(fileName,"finish peak picking!\n")

#====================
# TIC peak picking
time1<-Sys.time()
DataFilelist <- params$DataFiles
for(fileindex in 1:length(DataFilelist)) {
	inFilePath <- DataFilelist[[fileindex]]
	TICPeakpicking_ncdf4(inFilePath,params,denoised=FALSE)
}
time2<-Sys.time()
time2-time1

# EIC peak picking

#====================
# decomposition
DataFilelist<-params$DataFiles
fileindex = 2
inFilePath <- DataFilelist[[fileindex]]
windowID=96
#====================
time1<-Sys.time()
inFilePath <- DataFilelist[[fileindex]]
decomposition_ncdf4_test2(inFilePath,params,cl,isDistParallel=T,clustingType="h")
time2<-Sys.time()
time2-time1

time1<-Sys.time()
for(fileindex in c(1:length(DataFilelist)))
{
	inFilePath <- DataFilelist[[fileindex]]
	decomposition_ncdf4_test2(inFilePath,params,cl,isDistParallel=T,clustingType="h")
}
time2<-Sys.time()
time2-time1
#====================
# library searching

lib<-readMSP2Spec(filename="/Users/yni/StdsLibrary/KQC.txt",withRT=F)
specResults <- NULL
specResults<-c(specResults,specList)
libMatching_Top10(lib,specResults,600)

index1= 5
index2 = 6
specDistCal(specList[[index1]],lib[[1]],isWeight=T,isNist=T)
specDistCal_check2(specList[[index1]],lib[[1]],isWeight=T,isNist=T)
specDistCal_check(specList[[index1]],lib[[1]],isWeight=T,isNist=T)
#====================

#						for(curCluster in Clusters)
#						{
#							clust <- clust + 1
#							uniqueEICPeakList[which(FinalClustResult==curCluster),]$clust <- clust
#							uniquePeakList[which(rownames(uniquePeakList)%in%rownames(uniqueEICPeakList[which(FinalClustResult==curCluster),])),]$clust <- clust
#							mdlPkCandidates<-uniqueEICPeakList[which(FinalClustResult==curCluster),]
#							modelPkList <- rbind(modelPkList,mdlPkCandidates[which.max(mdlPkCandidates$shrp),])	
#						}	

#curCluster <- 1

uniquePeakList$mz
write.csv(uniqueEICPeakList$CompoundID,"/Users/yni/Desktop/test2.csv")
write.csv(uniquePeakList$mz,"/Users/yni/Desktop/test1.csv")

curEICProfiles <- curProfiles[which(rownames(uniqueEICPeakList)%in%rownames(curUniqueEICPeakList))]	

shrpList <- NULL
for (index in 1:length(curEICProfiles)) {
	curEICProfile <- curEICProfiles[[index]]
	curEICProfile$int <- 10000*curEICProfile$int/max(curEICProfile$int)
	#curEICProfile$int <- curEICProfile$int/sqrt(max(curEICProfile$int))							
	curProfile <- curEICProfile$int
	curET <- curEICProfile$ET
	
	pkPos <- which.max(curProfile)
	base <- c(curProfile[1],curProfile[length(curProfile)])
	time <- c(1,length(curProfile))
	linear.model <- lm(base~time)
	coeff <- as.numeric(linear.model$coefficients)
	pkPredict <- coeff[2]*pkPos+coeff[1]
	P25Height <- round((curProfile[pkPos]-pkPredict)*0.25)
	scan_diff <- abs(curET-curET[pkPos])
	
	LSide_pair <- subset(data.frame(Int=curProfile[1:(pkPos-1)],scanDist=sort(c(1:(pkPos-1)),decreasing=TRUE)),Int>=P25Height)
	RSide_pair <- subset(data.frame(Int=curProfile[(pkPos+1):length(curProfile)],scanDist=c(1:(length(curProfile)-pkPos))),Int>=P25Height)
	
	curShrp <- mean(c(median((curProfile[pkPos]-LSide_pair$Int)/LSide_pair$scanDist),median((curProfile[pkPos]-RSide_pair$Int)/RSide_pair$scanDist)))
	shrpList <- c(shrpList,curShrp)
	print(shrpList)
}


#IntScan_pair <- subset(data.frame(Int=curProfile,scanDist=scan_diff),Int>=P25Height&scan_diff>0)
#curPeakList[i,]$shrp_Mean1 <-round(median(abs(IntScan_pair$Int-curProfile[pkPos])/IntScan_pair$scanDist))
#curPeakList[i,]$shrp_Mean2 <-round(mean(abs(IntScan_pair$Int-curProfile[pkPos])/IntScan_pair$scanDist)/sqrt(curProfile[pkPos]))
#LSide_pair <- data.frame(Int=abs(curProfile[1:(pkPos-1)]-curPk$Intensity),scanDist=sort(c(1:(pkPos-1)),decreasing=TRUE))
#RSide_pair <- data.frame(Int=abs(curProfile[(pkPos+1):length(curProfile)]-curPk$Intensity),scanDist=sort(c(1:(length(curProfile)-pkPos)),decreasing=TRUE))
#scan_diff <- abs(c(startInd:endInd)-curPk$pkInd)
#					IntScan_pair <- subset(data.frame(Int=abs(curProfile-curPk$Intensity),scanDist=scan_diff))
#					pkPos <- which(c(startInd:endInd)==curPk$pkInd)
#					LSide_pair <- subset(IntScan_pair[c(1:(pkPos-1)),],scanDist>0)
#					RSide_pair <- subset(IntScan_pair[c((pkPos+1):length(curProfile)),],scanDist>0)
#					if (nrow(LSide_pair)>0&nrow(RSide_pair)>0) {
#						curPeakList[i,]$shrp_Mean2 <- round(mean(c(mean(LSide_pair$Int/LSide_pair$scanDist),mean(RSide_pair$Int/RSide_pair$scanDist))))
#					curPeakList[i,]$shrp_Med <- round(mean(c(median(LSide_pair$Int/LSide_pair$scanDist),median(RSide_pair$Int/RSide_pair$scanDist))))
#					curPeakList[i,]$shrp_Max <- round(mean(c(max(LSide_pair$Int/LSide_pair$scanDist),max(RSide_pair$Int/RSide_pair$scanDist))))
#					}

#					curPeakList[i,]$shrp_old <-round(Sharpness(yy=curProfile))	


curPeakList$zigNum <- 0
for(i in 1:nrow(curPeakList)) {
	curPk<-curPeakList[i,]
	startInd<-curPk$Lbound
	endInd<-curPk$Rbound
	LProfile <- curVecInt[startInd:curPk$pkInd]
	RProfile <- curVecInt[curPk$pkInd:endInd]
	
	zigzag_points <- length(which(LProfile[-1]- LProfile[-length(LProfile)]<0)) +
			length(which(RProfile[-length(RProfile)]- RProfile[-1]<0))
	
	#curPeakList[i,]$zigNum <- round(zigzag_points/(endInd-startInd)*100)
	curPeakList[i,]$zigNum <- zigzag_points
	
	
}



write.csv(curPeakList,"/Users/yanni/Desktop/peaks.csv")
mzGroup <- mzVec


plot(subset(curPeakList,isShared==0)$shrp_Mean)


plot(curPeakList$shrp_Mean,col=curPeakList$isShared+1)

PeakList <- do.call("rbind",groupResult)

PeakList<-unlist(groupResult,recursive=F)

time1<-Sys.time()
time2<-Sys.time()
time2-time1

# cwt testing
scalerange <- round(c(2,200)/2)
scales <-  seq(from=scalerange[1], to=scalerange[2], by=2)
#cwt_yan <- wavCWT_ModifyByYan(curVecInt,scale=scales)
library("MassSpecWavelet")
# apply cwt function from MassSpecWavelet
cwt_yan <- MSWtoCWT(curVecInt,scales)
# identify peak apex using funcs from wtmsa package
tree <- wavCWTTree_ModifyByYan(cwt_yan,type = "maxima") 
WT_result <- wavCWTPeaks_Modify(tree)
# identify peak boundaries using funcs from wtmsa package
tree_min <- wavCWTTree_ModifyByYan_Min(cwt_yan,type = "minima") 
WT_result_min <- wavCWTPeaks_Modify(tree_min)

length(which(WT_result_min$x%in%valleyInd==FALSE))


# boundary candidates detection using local minimal algorithm
intCutoff<-0 
valleySpan<-as.integer(params$Valley_span)#5
isMin <- peaks(-curVecInt, valleySpan)

## those low-int points as valleys too since sometime the consecutive low points may
## not be detected by local minima if they have the same values
isZeroThrsh<-curVecInt<=intCutoff # after baseline correction, negative intensity would be produced.
isValley <- isMin|isZeroThrsh
valleyInd<-which(isValley==TRUE)

curPeakList$Lbound <- 0
curPeakList$Rbound <- 0
curPeakList$Intensity <- curVecInt[curPeakList$pkInd]

for (i in 1:nrow(curPeakList)) {
	curPeak <- curPeakList[i,]$pkInd 
	curBoundary_range <- findInterval(curPeak,valleyInd)
	if (curBoundary_range>0) {
		curPeakList[i,]$Lbound  <- valleyInd[curBoundary_range]
		
		if (curBoundary_range==length(valleyInd)) {
			curPeakList[i,]$Rbound <- nPoints
		}else {
			curPeakList[i,]$Rbound <- valleyInd[curBoundary_range+1]
		}
	}
}

curPeakList$Lbound <- 0
curPeakList$Rbound <- 0
curPeakList$Intensity <- curVecInt[curPeakList$pkInd]

for (i in 1:(nrow(curPeakList)-1)) {
	curPeak <- curPeakList[i,]$pkInd 
	curBoundary_range <- findInterval(curPeak,valleyInd)
	
	if (curBoundary_range>0) {
		curPeakList[i,]$Lbound  <- valleyInd[curBoundary_range]
		
		if (curBoundary_range==length(valleyInd)) {
			curPeakList[i,]$Rbound <- nPoints # no candidates further
		}else {
			curPeak_next <- curPeakList[i+1,]$pkInd 
			Rbound <- valleyInd[curBoundary_range+1]
			
			if (Rbound>curPeak&Rbound<curPeak_next) {
				curPeakList[i,]$Rbound <- Rbound	
			}else {
				print(i)
				EIC_part <- curVecInt[curPeak:curPeak_next]
				add_BoundInd <- curPeak+which.min(EIC_part)-1
				curPeakList[i,]$Rbound <- add_BoundInd
				
				valleyInd <- sort(c(valleyInd,add_BoundInd))
				
			}		
		}
	}
}

for (i in 1:(nrow(curPeakList)-1)) {
	curPeak <- curPeakList[i,]$pkInd 
	curBoundary_range <- findInterval(curPeak,valleyInd)
	
	if (curBoundary_range>0) {
		curPeakList[i,]$Lbound  <- valleyInd[curBoundary_range]
		
		if (curBoundary_range==length(valleyInd)) {
			curPeakList[i,]$Rbound <- nPoints # no candidates further
		}else {
			curPeak_next <- curPeakList[i+1,]$pkInd 
			Rbound <- valleyInd[curBoundary_range+1]
			
			if (Rbound>curPeak&Rbound<curPeak_next) {
				curPeakList[i,]$Rbound <- Rbound	
			}else {
				curPeakList[i,]$Rbound <- NA
			}		
		}
	}
}


# the last peak

curPeak <- curPeakList[nrow(curPeakList),]$pkInd 
curBoundary_range <- findInterval(curPeak,valleyInd)
if (curBoundary_range>0) {
	curPeakList[nrow(curPeakList),]$Lbound  <- valleyInd[curBoundary_range]
	
	if (curBoundary_range==length(valleyInd)) {
		curPeakList[nrow(curPeakList),]$Rbound <- nPoints
	}else {
		curPeakList[nrow(curPeakList),]$Rbound <- valleyInd[curBoundary_range+1]
	}
}


# get valley index list that could properly divide EIC into CPFs
valleyInd <- c(curPeakList$Lbound,curPeakList$Rbound[nrow(curPeakList)])
curPeakList<- Update_Peak_list

# output
write.csv(as.matrix(cwt_yan),"/Users/yanni/Desktop/cwt_yan.csv")
playwith({plot(tree)},new=TRUE)

result_cwt <- data.frame(PeakInd=WT_result$x,SNR=as.numeric(attr(WT_result,"snr")))
result_cwt <- subset(result_cwt,SNR>=20)

write.csv(result_cwt,"/Users/yanni/Desktop/peakindex_cwt.csv")
playwith ({
			plot(c(1:length(curVecInt)),curVecInt,lwd=2,type = "l")
			points(WT_result$x,curVecInt[WT_result$x],col="red")
		},new=TRUE)

playwith ({
			plot(c(1:length(curVecInt)),curVecInt,lwd=2,type = "l")
			points(result_cwt$PeakInd,curVecInt[result_cwt$PeakInd],col="red")
			
			points(valleyInd,curVecInt[valleyInd],lwd=2,col="yellow")
			points(WT_result_min$x,curVecInt[WT_result_min$x]+1000,lwd=2,col="green")
		},new=TRUE)



# Parameters
x <- curVecInt
scale.range = deltat(x) * c(1, length(x)) # delta: time series/frequency
n.scale = 100
wavelet = "gaussian2"
shift = 5
variance = 1

x <- CWTCoeff
n.octave.min = 1
tolerance = 0
type = "maxima"

wtmm <- as.matrix(x)
span.min = 5
gap.max = 3
skip = NULL
sampling.interval = 1

x<- CWTTree
snr.min = 3
scale.range = NULL
length.min = 10
noise.span = NULL
noise.fun = "quantile"
noise.min = NULL

n.scale = 100
scale.range = deltat(x) * c(1, length(x)) # delta: time series/frequency
octave <- logb(scale.range, 2) # log2
scale <- ifelse1(n.scale > 1, 2^c(octave[1] + seq(0, n.scale - 2) * diff(octave)/(floor(n.scale) - 1), octave[2]), scale.range[1])
scale <- unique(round(scale/sampling.interval) * sampling.interval)
scales <- scale[1:65]


# Testing massSpecWavelet
library("MassSpecWavelet")
scalerange <- round(c(2,200)/2)
scales <-  seq(from=scalerange[1], to=scalerange[2], by=2)
scales <- seq(1,1000,10)
SNR.Th <- 1
nearbyPeak <- TRUE

time1<-Sys.time()
wCoefs <- cwt(curVecInt,scales=scales,wavelet="mexh")


time2<-Sys.time()
time2-time1

time1<-Sys.time()
colnames(wCoefs) <- scales
localMax <- getLocalMaximumCWT(wCoefs)
ridgeList <- getRidge(localMax)
time2<-Sys.time()
majorPeakInfo <- identifyMajorPeaks(curVecInt, ridgeList, wCoefs, scales,SNR.Th = SNR.Th, nearbyPeak=nearbyPeak)## Plot the identified peaks
time2<-Sys.time()
time2-time1

# output
write.csv(as.matrix(wCoefs),"/Users/yanni/Desktop/cwt_MassSpecWavelet.csv")
playwith ({
			image(c(1:length(curVecInt)),scales,wCoefs[c(1:length(curVecInt)),],col=terrain.colors(100))
			box()
		})
plotLocalMax(localMax,wCoefs,range=c(1,length(curVecInt)))
plotRidgeList(ridgeList, wCoefs, range=c(1,length(curVecInt)))
result <- data.frame(peakInd = majorPeakInfo$allPeakIndex,peakSNR = round(majorPeakInfo$peakSNR))
write.csv(subset(result,peakSNR>=1),"/Users/yanni/Desktop/peakindex_massSpec2.csv")
peakIndex = majorPeakInfo$allPeakIndex
peakIndex = majorPeakInfo$peakIndex
playwith ({
			plot(c(1:length(curVecInt)),curVecInt,lwd=2,type = "l")
			points(peakIndex,curVecInt[peakIndex],col="red")
			#text(peakIndex,curVecInt[peakIndex],c(1:length(peakIndex)))
			#plotPeak(curVecInt, peakIndex, range=c(1,length(curVecInt)))
		},new=TRUE)

a[c(1067:1070)]

n.scale.min = 3
good <- which(unlist(lapply(ridgeList, function(x, n.scale.min) length(x[[1]]) > n.scale.min, n.scale.min = n.scale.min)))



endtime <- as.vector(unlist(lapply(ridgeList, function(x) x[1])))
ridgeList <- ridgeList[order(endtime)]
length(ridgeList[[1]])
order(c(1,3,9))
# Testing XCMS-MassSpecWavelet
source("/Users/yanni/Documents/workspace/XCMS_R/cwTools.R")

time1<-Sys.time()
wCoefs <- MSW.cwt(curVecInt,scales=scales,wavelet="mexh")
time2<-Sys.time()
time2-time1

time1<-Sys.time()
colnames(wCoefs) <- scales
localMax <- MSW.getLocalMaximumCWT(wCoefs)
ridgeList <- MSW.getRidge(localMax)
majorPeakInfo1 <- identifyMajorPeaks(curVecInt, ridgeList, wCoefs, scales,SNR.Th = SNR.Th, nearbyPeak=nearbyPeak)## Plot the identified peaks
time2<-Sys.time()
time2-time1

peakIndex = majorPeakInfo1$peakIndex
playwith ({
			plot(c(1:length(curVecInt)),curVecInt,lwd=2,type = "l")
			points(peakIndex,curVecInt[peakIndex],col="red")
			#text(peakIndex,curVecInt[peakIndex],c(1:length(peakIndex)))
			#plotPeak(curVecInt, peakIndex, range=c(1,length(curVecInt)))
		},new=TRUE)





if (nrow(modelPkList)>=1) { 
	decom = 1
	while (decom==1) {
		specList<-NULL
		componentList<-NULL
		mdlPkIDs<-rownames(modelPkList)							
		lbound<-min(sharedEICPeakList$Lbound)
		rbound<-max(sharedEICPeakList$Rbound)
		
		# get all EIC profiles within the window, X matrix format:column-EIC,row-scan
		X<-NULL
		mzVec<-unique(sharedEICPeakList$mz)
		nMz<-length(mzVec)			
		for(i in 1:nMz)
		{
			mz<-mzVec[i]
			curMz_list <- sharedEICPeakList[which(sharedEICPeakList$mz==mz),]									
			mzInd<-which(vectorMz==mz)
			startInd<-(mzInd-1)*totalscan+1
			endInd<-startInd+totalscan-1
			curVecInt <- vecInt[startInd:endInd]
			# get the shared peak feature left and right boundary
			curLbound <- curMz_list[1,]$lboundInd
			curRbound <- curMz_list[1,]$rboundInd
			curProfile <- merge(data.frame(ET=c(lbound:rbound)),data.frame(ET=c(curLbound:curRbound),int=curVecInt[curLbound:curRbound]),by.y= "ET",all.x=T)
			curProfile$int[is.na(curProfile$int)]<-0 
			X<-cbind(X,curProfile$int)
		}
		colnames(X)<-mzVec	
		nScans<-nrow(X)
		
		# Get model peak profiles as S matrix, each column is the EIC for one model peak	
		# add unique mass information to spec and component list
		S<-NULL
		for(i in 1:nrow(modelPkList))
		{
			curProfile<-merge(data.frame(ET=as.integer(lbound:rbound)),allprofiles[[rownames(modelPkList[i,])]],by.y="ET",all.x=T)
			curProfile$int[is.na(curProfile$int)]<-0 # expand the profile by filling zero values to those uncovered area
			S<-cbind(S,curProfile$int)
			#modelPkList[i,]$area<-sum(curProfile$int)		
			
			# init the component & spec list	
			specList[[i]]<-data.frame(mz=as.integer(),int=as.integer())
			names(specList)[i]<-paste(modelPkList$pkInd[i],windowID,i,modelPkList$mz[i],compNum+i-1)
			componentList[[i]]<-data.frame(mz=as.integer(),pkInd=as.integer(),Intensity=as.integer(),lboundInd=as.integer(),rboundInd=as.integer(),
					area=as.integer(),windowID=as.integer(),compoundID=as.integer(),isModel=as.integer(),isUnique=as.integer(),factor=as.numeric(),compID=as.integer())
			
			# add unique mass information to spec list and component list
			curClusterID <- modelPkList[i,]$clust
			curCluster_uniqueEICPeakList <- subset(uniquePeakList,clust==curClusterID)
			curCluster_uniqueIons <- data.frame(curCluster_uniqueEICPeakList[,c("mz","Intensity")])
			colnames(curCluster_uniqueIons) <- c("mz","int")
			specList[[i]] <- rbind(specList[[i]],curCluster_uniqueIons)
			
			curCluster_uniqueEICPeakList$windowID <- windowID
			curCluster_uniqueEICPeakList$compoundID <- i
			curCluster_uniqueEICPeakList$isModel <-	0	
			curCluster_uniqueEICPeakList$isUnique <- 1	
			curCluster_uniqueEICPeakList$factor <- 0
			curCluster_uniqueEICPeakList$compID <-compNum+i-1	
			curCluster_uniqueEICPeakList[rownames(modelPkList[i,]),]$isModel = 1				
			componentList[[i]]<-rbind(componentList[[i]],subset(curCluster_uniqueEICPeakList,select=c("mz","pkInd","Intensity","lboundInd","rboundInd","area",
									"windowID","compoundID","isModel","isUnique","factor","compID")))
		}		
		colnames(S)<-mdlPkIDs
		
		# Residual minimization
		for(i in 1:ncol(X))
		{
			curMz<-colnames(X)[i]						
			M<-X[,i]
			
			srcIDs<-mdlPkIDs
			A<-optim(par=rep(0,ncol(S)),fn=fr, gr = NULL,M,S,lower = 0, method="L-BFGS-B")$par
			A<-as.matrix(A)
			
			# restore each ion signals and save them into component list ans speclist
			for(m in 1:length(srcIDs))
			{
				curmdlPkID<-srcIDs[m]
				j<-which(mdlPkIDs==curmdlPkID)
				
				curPkInd <- modelPkList[curmdlPkID,]$pkInd
				curSharedLbound <- modelPkList[curmdlPkID,]$lboundInd
				curSharedRbound <- modelPkList[curmdlPkID,]$rboundInd
				
				if (A[m]==0) {
					componentList[[j]]<-rbind(componentList[[j]],data.frame(mz=as.integer(curMz),pkInd=curPkInd,lboundInd=curSharedLbound,rboundInd=curSharedRbound,Intensity=10,
									area=10,windowID=windowID,compoundID=j,isModel=0,isUnique=0,factor=0,compID=compNum+m-1))
				}
				
				if(A[m]>0)
				{
					curIntensity <- as.integer(modelPkList[curmdlPkID,]$Intensity*A[m,])
					curArea <- as.integer(modelPkList[curmdlPkID,]$area*A[m,])
					
					specList[[j]]<-rbind(specList[[j]],data.frame(mz=as.integer(curMz),int=curIntensity))
					componentList[[j]]<-rbind(componentList[[j]],data.frame(mz=as.integer(curMz),pkInd=curPkInd,lboundInd=curSharedLbound,rboundInd=curSharedRbound,
									Intensity=curIntensity,area=curArea,windowID=windowID,compoundID=j,isModel=0,isUnique=0,factor=A[m],compID=compNum+m-1))						
				}							
			}
		}
		
		# Component splitting correction: checking pairwise resovlved spectra similarity
		Num <- nrow(modelPkList)									
		if (Num>1) {
			RepeatModelPk <- NULL									
			indexPairs<-combinations(Num,2)
			for(i in 1:nrow(indexPairs))
			{
				index <- indexPairs[i,]
				index1 <- index[1]
				index2 <- index[2]
				
				##calculate time difference between two model peaks
				modelPKET <- modelPkList$pkInd[c(index1,index2)]	
				
				#fragments with the maximal time difference 10 scans are considered a compound
				if (abs(modelPKET[1]-modelPKET[2]) <=10&nrow(specList[[index1]])!=0&nrow(specList[[index2]])!=0)
				{
					score <- specDistCal(specList[[index1]],specList[[index2]],isWeight=T,isNist=T)
					if (0.9*score+(1-(abs(modelPKET[1]-modelPKET[2])/10))*100 >=SpecSimilarity_Th) {
						selectedModelPkID <- ifelse(modelPkList$Intensity[index1]>=modelPkList$Intensity[index2],mdlPkIDs[index2],mdlPkIDs[index1])
						RepeatModelPk <- c(RepeatModelPk,selectedModelPkID)
					}
				}
			}						
			# select unique filtered model peaks 
			RepeatModelPk <- unique(RepeatModelPk)
			modelPkList <- modelPkList[!rownames(modelPkList)%in%RepeatModelPk,]
		}
		
		if (nrow(modelPkList) == Num) {
			decom = 0
			compNum <- compNum+Num	
			## output model pk profiles
			for (model in c(1:nrow(modelPkList))) {
				curModel <- modelPkList[model,]
				curModelInt <- S[,model]
				strOutput <- paste(strOutput,getModelPkProfile(curModel,curModelInt,windowID,lbound,rbound,model,compNum),sep="\n\n")
			}
			
			# flag analyzed EIC peaks 
			EICpeaklist[rownames(localEICPeakList),]$flag<-1
			componentResults<-c(componentResults,componentList)
			specResults<-c(specResults,specList)
		}
	}# decom while statement
	
}else { 
	if (nrow(sharedEICPeakList)>=5) {
		specList<-NULL
		componentList<-NULL
		sharedPeaks_density <- densityMclust(sharedEICPeakList$pkInd)
		CompNum <- attr(sharedPeaks_density,"parameters")$variance$G
		
		if (CompNum>=1&&nrow(sharedEICPeakList)>=CompNum) {
			ClustResult <- kmeans(sharedEICPeakList$pkInd,CompNum)$cluster
			Clusters <- as.vector(which(table(ClustResult)>=5))
			#Clusters <- unique(ClustResult)
		}else {
			Clusters <- NULL
		}
		
		if (length(Clusters)>0) {
			sharedEICPeakList$clust <- ClustResult
			
			for (i in 1:length(Clusters)) {
				cluster <- Clusters[i]
				
				curCluster <- subset(sharedEICPeakList,clust==cluster)
				# dummy model peak used for component peak index, left/right boundary annotation only
				dummy_modelPeak <- curCluster[which.max(curCluster$Intensity),]
				curMz<-as.character(dummy_modelPeak$mz)
				maxInd<-dummy_modelPeak$pkInd
				curLbound<-dummy_modelPeak$Lbound
				curRbound<-dummy_modelPeak$Rbound
				
				curCluster_uniqueIons <- data.frame(curCluster[,c("mz","Intensity")])
				colnames(curCluster_uniqueIons) <- c("mz","int")
				
				specList[[i]]<-data.frame(mz=as.integer(),int=as.integer())
				names(specList)[i]<-paste(maxInd,windowID,cluster,curMz,compNum+i-1)
				specList[[i]] <- rbind(specList[[i]],curCluster_uniqueIons)
				
				curCluster$windowID <- windowID
				curCluster$compoundID <- i
				curCluster$isModel <- 0	
				curCluster$isUnique <- 0	
				curCluster$factor <- 0
				curCluster$compID <-compNum+i-1	
				curCluster[rownames(dummy_modelPeak[i,]),]$isModel = 1	
				
				componentList[[i]]<-data.frame(mz=as.integer(),pkInd=as.integer(),Intensity=as.integer(),lboundInd=as.integer(),rboundInd=as.integer(),
						area=as.integer(),windowID=as.integer(),compoundID=as.integer(),isModel=as.integer(),isUnique=as.integer(),factor=as.numeric(),compID=as.integer())
				componentList[[i]]<-rbind(componentList[[i]],subset(curCluster,select=c("mz","pkInd","Intensity","lboundInd","rboundInd","area",
										"windowID","compoundID","isModel","isUnique","factor","compID")))
			}
			# update compNum
			compNum <- compNum + length(Clusters)
			# flag analyzed EIC peaks 
			EICpeaklist[rownames(localEICPeakList),]$flag<-1
			componentResults<-c(componentResults,componentList)
			specResults<-c(specResults,specList)
		}
	}
}# Check localEICPeakList 2	