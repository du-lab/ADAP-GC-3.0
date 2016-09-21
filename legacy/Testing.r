# ADAP3 development

# commands @ laptop 
codeDir <<-"/Users/yanni/Desktop/yni1/work/SourceCode/ADAP3"
source(paste(codeDir,"pipeline.r",sep="/"))

params<-readParaFromCmd("/Users/yanni/DataTest/JiaLab/stds/Stds.db")
params$codeDir<-codeDir
fileInfo<-read.csv(paste(params$WorkDir,"output/",params$JobName,"_RT.csv",sep=""))
params$delaytime<-fileInfo$firstRT[1]
params$ScanInterval<-fileInfo$ScanInterval[1]

# commands @ local computer
codeDir<-"/Users/yni1/Desktop/work/SourceCode/ADAP3"
source(paste(codeDir,"pipeline.r",sep="/"))

params<-readParaFromCmd("/Users/yni1/DataTest/JiaLab/stds/Stds.db")
params$codeDir<-codeDir
fileInfo<-read.csv(paste(params$WorkDir,"output/",params$JobName,"_RT.csv",sep=""))
params$delaytime<-fileInfo$firstRT[1]
params$ScanInterval<-fileInfo$ScanInterval[1]

params<-readParaFromCmd("/Users/yni1/WingsTest/run7/run7.db")
codeDir<-"/Users/yni1/Desktop/work/SourceCode/ADAP3"
source(paste(codeDir,"pipeline.r",sep="/"))

params$codeDir<-codeDir
fileInfo<-read.csv(paste(params$WorkDir,"output/",params$JobName,"_RT.csv",sep=""))
params$delaytime<-fileInfo$firstRT[1]
params$ScanInterval<-fileInfo$ScanInterval[1]

fileindex=1
inFilePath<- DataFilelist[[fileindex]]
windowID=50#35#149#247#214

# commands @ aether server
setwd("/home/yni/workspace/SourceCode/")
codeDir<-"ADAP3"
source(paste(codeDir,"pipeline.r",sep="/"))
params<-readParaFromCmd("/home/yni/DataTest/LinearStds2/LinearStds2.db")
params$codeDir<-codeDir

# pk picking
parTICPeakpicking(params,denoised=F)
parEICpeakpicking(params)

# deconvolution
fileInfo<-read.csv(paste(params$WorkDir,"output/",params$JobName,"_RT.csv",sep=""))
params$delaytime<-fileInfo$firstRT[1]
params$ScanInterval<-fileInfo$ScanInterval[1]
deconvolution(params,isDistParallel=F,clustingType="h",isDataParalell=T)







mz=73

write.csv(curPeakList,"/Users/yni1/Desktop/curPeakList_modify2.csv")
write.csv(curPeakList,"/Users/yni1/Desktop/curPeakList.csv")

dest <- densityMclust(uniqueEICPeakList$pkInd)
attr(dest,"parameters")$variance$G


str(dest)
x11()
plot(dest,data=uniqueEICPeakList$pkInd)

x2 <- test

test2 <- test
test2$ET <- test$ET+10
y2 <- test2

data1<-merge(x2,y2,by.x="ET",by.y="ET",all=isUnion)
data1[is.na(data1)]<-0

x11()
plot(data1[,2],type="l",main=paste("10 scan diff-arcosine angle",round(dotProd)))
points(data1[,3],type="l",col="red")

dotProd=dotProduct(x=data1[,2], y=data1[,3])


LW <- uniqueEICPeakList$Rbound-uniqueEICPeakList$pkInd
RW <- uniqueEICPeakList$pkInd-uniqueEICPeakList$Lbound

abs(LW-RW)

x11()
plot(uniqueEICPeakList$mz,abs(LW-RW))
text(uniqueEICPeakList$mz,abs(LW-RW), as.character(uniqueEICPeakList$mz,col="blue"))


x11()
plot(uniqueEICPeakList$mz,peakWidth)
text(uniqueEICPeakList$mz,peakWidth, as.character(uniqueEICPeakList$mz,col="blue"))


x11()
plot(clustResult1,lwd=2,cex=1.5,font.axis=2,ylim=c(0,90))
axis(2,lwd=2)
abline(h=80,lty=2,col="red",lwd=2)

clustResult2 <- cutree(clustResut,h=100)


which(clustResult2==3)
which(clustResult2==4)
uniqueEICPeakList <- subset(localEICPeakList,isShared==0)
updatePkInd <- localEICPeakList$pkInd[-c(which(clustResult2==3),which(clustResult2==4))]


and 
which(ClustResult==1)

rownames(uniqueEICPeakList)

ClustResult

x11()
plot(dest,data=uniqueEICPeakList$pkInd)



clustResult1 <- hclust(dist(uniqueEICPeakList$pkInd))
clustResult2 <- cutree(clustResult1,h=80)
FinalClustResut<-unique(clustResult2)
noisyPeaks <- which(clustResult2==as.integer(which(table(clustResult2)<=2)))

uniqueEICPeakList <- uniqueEICPeakList[-noisyPeaks,]




qualify_uniqueEICPeakList <- subset(uniqueEICPeakList,Intensity>500)
uniqueProfiles<-vector("list",nrow(qualify_uniqueEICPeakList))
for(i in 1:nrow(qualify_uniqueEICPeakList))
{
	curEICpk<-qualify_uniqueEICPeakList[i,]			
	startInd<-curEICpk$offset+curEICpk$lboundInd
	endInd<-curEICpk$offset+curEICpk$rboundInd
	ET<-curEICpk$lboundInd:curEICpk$rboundInd
	curProfile<-data.frame(ET=ET,int=vecInt[startInd:endInd]) 
	uniqueProfiles[[i]]<-curProfile
	names(uniqueProfiles)[i]<-rownames(curEICpk)	
}


peakMclust <- Mclust(localEICPeakList$pkInd)
peakMclust <- mclustBIC(localEICPeakList$pkInd)
str(peakMclust,parameters = TRUE)


peakMclust <- Mclust(a)


# test
a <- c(rnorm(50,10,5))
dest <- densityMclust(a)
summary(dest,parameters=TRUE)
x11()
plot(dest)

# ADAP 2 evaluation
codeDir <<-"/Users/yni1/Desktop/work/SourceCode/ADAP2"
source(paste(codeDir,"pipeline.r",sep="/"))

params<-readParaFromCmd("/Users/yni1/DataTest/JiaLab/U50Stds2013/U50Stds2013.db")
params$codeDir<-codeDir

fileInfo<-read.csv(paste(params$WorkDir,"output/",params$JobName,"_RT.csv",sep=""))
params$delaytime<-fileInfo$firstRT[1]
params$ScanInterval<-fileInfo$ScanInterval[1]
delaytime <- params$delaytime
ScanInterval <- params$ScanInterval


fileindex=8
inFilePath<- DataFilelist[[fileindex]]


windowID=158
mz=73


libFile <- "/Users/yni1/StdsLibrary/KQC_old.txt"
filename <- "/Users/yni1/StdsLibrary/KQC.txt"
filename <- params$libFile_decon

curFragment<-trim(SpecPoints[6])#get current pair
curFragment<-strsplit(curFragment," ")[[1]]#split the current pair


delaytime<-params$delaytime
ScanInterval<-params$ScanInterval



playwith({
			x11(width=15,height=10)
			par(mfrow=c(5,5))
			
			for (j in c(6070:6099)) {
				plot(vecInt[totalscan*c(0:(length(vectorMz)-1))+j],type="h",lwd=5)
			}
			
		})


playwith({
		
			
				plot(vecInt[totalscan*c(0:560)+1872],type="h",lwd=4,col="blue")
			

			
		})

plot(vecInt[(totalscan*33):(totalscan*34)])
which(vectorMz==73)
which(vecInt>=5000)totalscan
which.max(vecInt)/length(vectorMz)
19.8
1200*(19.8-4.2)
length(vectorMz)


a <- hclust(dist(localEICPeakList$pkInd))
x11()
plot(a)

dist(c(5,2,1))

filename <- "/Users/yni1/DataTest/JiaLab/stds/output/decon_Feb202013_newlib/TargetLibMatch/S0.2_1SpectraList.txt"

refSpec <- lib

inputSpec <- reference
minSpecSimilarity <- 700

# draft keeping
##the main routine to call all the steps
#parGcProprocess <-function(params)
#{	
#	#get the first scan and scan interval from *_RT.csv file which is produced by SPlitEIC in C code
#	fileInfo<-read.csv(paste(params$WorkDir,"/output/",params$JobName,"_RT.csv",sep=""))
#	params$delaytime<-fileInfo$firstRT[1]
#	params$ScanInterval<-fileInfo$ScanInterval[1]
#	#denoise TIC before TIC peakpicking
##	parTICDenoising(params) 
#	#TIC peak picking and boundary detection in order to set the preliminary decon window	
#	parTICPeakpicking(params,denoised=F)
#	#denoising EIC before EIC peak picking
#	parDenoising(params,codeDir="code")
#	
#	#EIC peak picking and boundary detection in order to get EIC profile for each fragment ion
#	parEICpeakpicking(params)
#	#set the exluded ion list
#	params$NonUmassVec<-c(1:50,73,147,221)
#	#deconvolute the unresolved the peaks
#	deconvolution(params,denoised=F)
#	
#}



#								## >= one model peaks
#								if (nrow(modelPkList)>1) {
#									rmmodelPks<-data.frame(row=as.integer(),col=as.integer())						
#									modelPkDistance <- as.matrix(modelPkDist(modelPkList,profiles))
#									
#									for (i in c(1:nrow(modelPkDistance))) {
#										subDistance <- modelPkDistance[i,]
#										par_name <- rownames(modelPkDistance)[i]
#										par_index <- modelPkList[par_name,]$pkInd
#										
#										for (j in c(1:length(subDistance))) {
#											dau_name <- names(subDistance)[j]
#											dau_index <- modelPkList[dau_name,]$pkInd
#											
#											comb_factor <- abs(par_index - dau_index)*modelPkDistance[i,j]
#											
#											if (comb_factor<= 12||(abs(par_index - dau_index)<=6&modelPkDistance[i,j]<=10)) {
#												rmmodelPks <- rbind(rmmodelPks,data.frame(row=i,col=j))
#											}	
#										}
#									}
#									
#									rmmodelPks <- subset(rmmodelPks,row!=col)
#									#######################################
#									#Step 2: grouping pairs
#									#######################################
#									if (nrow(rmmodelPks)>=1) {
#										m=c(1:nrow(rmmodelPks))
#										l = 1
#										modelGroup <- list()
#										
#										while (length(m)>0) {
#											temp <- NULL
#											for (i in m[2:length(m)]) {
#												if (length(rmmodelPks[i,][(rmmodelPks[i,]%in%rmmodelPks[m[1],])])>=1) {
#													temp <- c(temp, i) 
#												}
#											}
#											temp <- c(temp, m[1]) 
#											m <- m[!m%in%temp]
#											modelGroup[[l]] <- temp
#											l <- l+1
#										}
#										
#										###########################################
#										#Step 3: get representative for each group
#										###########################################
#										
#										modelPkList_temp <- NULL
#										rmmodelPkList_temp <- NULL
#										for (i in c(1:length(modelGroup))) {
#											rmmodelPkList <- modelPkList[unique(unlist(rmmodelPks[modelGroup[[i]],])),]
#											modelPkList_temp <- rbind(modelPkList_temp,rmmodelPkList[which.max(rmmodelPkList$score),])
#											rmmodelPkList_temp <- rbind(rmmodelPkList_temp,modelPkList[rownames(modelPkList)%in%rownames(rmmodelPkList),])
#											##Only keep the one with highest score
#											#modelPkList <- rbind(modelPkList,modelPkList_temp)	
#										}
#										
#										## get removed peaks and update model Peak list
#										rmmodelPkList_temp <- rmmodelPkList_temp[!rownames(rmmodelPkList_temp)%in%rownames(modelPkList_temp),]
#										modelPkList <- modelPkList[!rownames(modelPkList)%in%rownames(rmmodelPkList_temp),]
#										
#									}
#								}


x11()
hist(localEICPeakList$pkInd,10)

assign("Global.curProfs", value=uniqueProfiles, envir = .GlobalEnv)
r<-DistCal4(isUnion=F)# non-parallel version
distance <- as.dist(r)
maxIntraDist<- ClusterDist_Th 
clustResut<-hclust(distance)
FinalClustResut<-cutree(clustResut,h=maxIntraDist)
Clusters<-unique(FinalClustResut)
if(length(Clusters)>1)isMultGroups<-T

clustResult1 <- hclust(dist(localEICPeakList$pkInd))
clustResult1 <- hclust(dist(updatePkInd))

clustResult2 <- cutree(clustResut,h=100)
FinalClustResut<-unique(clustResult2)

which(clustResult2==3)
which(clustResult2==4)

updatePkInd <- localEICPeakList$pkInd[-c(which(clustResult2==3),which(clustResult2==4))]

x11()
plot(clustResult1,lwd=2,cex=1.5,font.axis=2,ylim=c(0,90))
axis(2,lwd=2)
abline(h=100,lty=2,col="red",lwd=2)


if (nrow(localEICPeakList)>1) {
	
	allprofiles<-vector("list",nrow(localEICPeakList))
	for(i in 1:nrow(localEICPeakList))
	{
		# get EIC peak profile
		curEICpk<-localEICPeakList[i,]			
		startInd<-curEICpk$offset+curEICpk$lboundInd
		endInd<-curEICpk$offset+curEICpk$rboundInd
		ET<-curEICpk$lboundInd:curEICpk$rboundInd
		curProfile<-data.frame(ET=ET,int=vecInt[startInd:endInd]) 
		allprofiles[[i]]<-curProfile
		names(allprofiles)[i]<-rownames(curEICpk)	
	}	
	
	CandidatePeaklist<-subset(localEICPeakList,StN>=StN_Th&shrp>=Shrp_Th&gss<=Gss_Th)
	
	if (nrow(CandidatePeaklist)>0) {
		# linear combination of 4 scores: mass score,gassian similarity score,Signal to Noise, Intensity
		CandidatePeaklist$f1<-scale(CandidatePeaklist$mz,min(CandidatePeaklist$mz),diff(range(CandidatePeaklist$mz)))	
		CandidatePeaklist$f2<-scale(cos(CandidatePeaklist$gss*pi/180),min(cos(CandidatePeaklist$gss*pi/180)),diff(range(cos(CandidatePeaklist$gss*pi/180))))
		CandidatePeaklist$f3<-scale(log(CandidatePeaklist$StN),min(log(CandidatePeaklist$StN)),diff(range(log(CandidatePeaklist$StN))))
		CandidatePeaklist$f4<-scale(log(CandidatePeaklist$Intensity),min(log(CandidatePeaklist$Intensity)),diff(range(log(CandidatePeaklist$Intensity))))	
		CandidatePeaklist$f1[is.na(CandidatePeaklist$f1)]<-0
		CandidatePeaklist$f2[is.na(CandidatePeaklist$f2)]<-0
		CandidatePeaklist$f3[is.na(CandidatePeaklist$f3)]<-0
		CandidatePeaklist$f4[is.na(CandidatePeaklist$f4)]<-0
		
		scores<-(10/8)*(0.1*CandidatePeaklist$f1+0.3*CandidatePeaklist$f2+0.2*CandidatePeaklist$f3+0.2*CandidatePeaklist$f4)
		CandidatePeaklist$score<-scale(scores,F,0.001)
		
		goodShapePeaklist<-subset(CandidatePeaklist,score>=min(CandidatePeaklist$score)+diff(range(CandidatePeaklist$score))*TotalScore_factor)
		goodShapePeaklist<-subset(goodShapePeaklist,!mz%in%CurNonUmassVec)
		nPks<-nrow(goodShapePeaklist)		
		modelPkList<-NULL
		isMultGroups<-FALSE
		mdlPkIDs<-NULL
		
		if(nPks>0)
		{
			# get good peak mirror image profiles
			profiles<-vector("list",nrow(goodShapePeaklist))
			for(i in 1:nrow(goodShapePeaklist))
			{
				# get intensity vector of current mz from long intensity vector
				curEICpk<-goodShapePeaklist[i,]					
				startInd<-curEICpk$offset+curEICpk$lboundInd
				endInd<-curEICpk$offset+curEICpk$rboundInd
				ET<-curEICpk$lboundInd:curEICpk$rboundInd
				curProfile<-data.frame(ET=ET,int=vecInt[startInd:endInd])
				profiles[[i]]<-curProfile
				names(profiles)[i]<-rownames(curEICpk)	
			}
			
			# update each profile by creating mirror images choose the one with better GSS
			for(i in 1:length(profiles))
			{
				curPkId<-names(profiles[i])
				curPk<-goodShapePeaklist[curPkId,]
				pos<-curPk$gssPos
				profiles[[i]]<-copyShape(p=profiles[[i]],curPk$pkInd,from=pos)
				goodShapePeaklist[curPkId,]$lboundInd<-profiles[[i]]$ET[1]
				goodShapePeaklist[curPkId,]$rboundInd<-profiles[[i]]$ET[nrow(profiles[[i]])]
			}
			
			# Hierachical clustering: decide the number of groups within current decon window
			if(nPks>=2)
			{
				# save EIC profies in global variable for broadcasting to slave nodes for parallel computing
				assign("Global.curProfs", value=profiles, envir = .GlobalEnv)
				
				# calculate pairwise distance among the EICs				
				if(isDistParallel)
				{
					r<-parDistCal4(cl,isUnion=F)# parallel verison
				}else
				{
					r<-DistCal4(isUnion=F)# non-parallel version
				}	
				
				# convert distance matrix to triangle matrix
				distance <- as.dist(r)
				maxIntraDist<- ClusterDist_Th 
				
				if(clustingType=="h")
				{
					clustResut<-hclust(distance)
					FinalClustResut<-cutree(clustResut,h=maxIntraDist)
					Clusters<-unique(FinalClustResut)
					if(length(Clusters)>1)isMultGroups<-T								
				}	
			}
			
			
			# Model peak selection 
			modelPkList<-NULL
			
			# single cluster
			if (!isMultGroups && nPks >=1)
			{
				if (nPks==1) {
					modelPkList<-goodShapePeaklist
				} else 
				{
					modelPkList<- goodShapePeaklist[which.max(goodShapePeaklist$score),]
				}
			}
			
			# multiple clusters
			if(isMultGroups)
			{
				if(nPks==2)
				{
					# only two candidates for two clusters, they are simply assigned as model peaks						
					modelPkList<-goodShapePeaklist
				}else
				{	
					# select model peak from each group[!isTooBigIntraDist]
					for(curCluster in Clusters)
					{
						# get masses belong to current cluster
						curGroup<-names(FinalClustResut[FinalClustResut==curCluster])
						mdlPkCandidates<-goodShapePeaklist[curGroup,]	
						if(nrow(mdlPkCandidates)>0)
						{
							# get model peak with the max score
							currentMdlPk<-mdlPkCandidates[which.max(mdlPkCandidates$score),]
							modelPkList<-rbind(modelPkList,currentMdlPk)
						}
					} # for loop	
				}			
			}						
		}
		


