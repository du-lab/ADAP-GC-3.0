
##############################################
# write each model pk prfile into standard format
# ready for outputing text file
# Nov. 22, 2011
##############################################
fr <- function(x,M,S) {
	m<-0
	for(i in 1:ncol(S))
	{
		m<-m+x[i]*S[,i]
	}
	sum((M-m)^2)
}

getModelPkProfile<-function(curModel,curModelInt,windowID,lbound,rbound,model,compNum)
{
	
	strBlock<-sprintf("peakInd:%s",curModel$pkInd)
	strBlock<-paste(strBlock,sprintf("model mass:%d",curModel$mz),sep="\n")
	strBlock<-paste(strBlock,sprintf("model gss:%.2f",curModel$gss),sep="\n")
	strBlock<-paste(strBlock,sprintf("windowID:%d",windowID),sep="\n")
	strBlock<-paste(strBlock,sprintf("compoundID:%d",model),sep="\n")
	strBlock<-paste(strBlock,sprintf("compID:%d",compNum+model-1),sep="\n")
	strBlock<-paste(strBlock,sprintf("window lboundInd:%d",lbound),sep="\n")
	strBlock<-paste(strBlock,sprintf("window rboundInd:%d",rbound),sep="\n")
	
	ModelIntBlock<-NULL
	for(j in 1:length(curModelInt))
	{
		strPair<-paste(curModelInt[j]," ",sep="")
		ModelIntBlock<-paste(ModelIntBlock,strPair)
	}
	strBlock<-paste(strBlock,ModelIntBlock,sep="\n")
	strBlock
}

modelPkDist <- function(modelPkList,profiles) {
	if (!is.null(nrow(modelPkList))) {
		modelPkprofiles<-profiles[as.character(rownames(modelPkList))]
		assign("Global.curProfs", value=modelPkprofiles, envir = .GlobalEnv)
		r<-DistCal4(isUnion=F)
		#distance <- as.dist(r)
		#distance
		r
	}
}

matrix.index <- function(a, value) {
	idx <- which(data.frame(a)==value)
	col.num <- ceiling(idx/nrow(a))
	row.num <- idx - (col.num-1) * nrow(a)
	return(c(row.num, col.num))
}


deconvolution<-function(params,isDistParallel=T,clustingType="h",isDataParalell)
{

	library("snow")
	DataFilelist<-params$DataFiles
	
	# Decide the optimal number of nodes for usage
	params$nNode <- min(as.integer(params$nNode),length(DataFilelist))
	cl<-makeCluster(params$nNode,type=params$clustType)
	
	if(isDataParalell==T)
	{
		time1<-Sys.time()
		clusterExport(cl,"decomposition")#assign the current profiles to all nodes
		clusterExport(cl,"parseFileName")
    
		clusterApply(cl,DataFilelist,decomposition,params,cl,isDistParallel=F,clustingType)#parallel version
		time2<-Sys.time()
		time2-time1
	}else
	{
		time1<-Sys.time()
		for(fileindex in 1:length(DataFilelist))
		#for(fileindex in 6:7)
		{
			inFilePath <- DataFilelist[[fileindex]]
			decomposition(inFilePath,params,cl,isDistParallel=T,clustingType="h")
		}
		time2<-Sys.time()
		time2-time1
	}	
	stopCluster(cl)
}

decomposition<-function(inFilePath,params,cl,isDistParallel,clustingType)
{
	library(gdata)
	library(gtools)
	library(cluster)	
	library("ncdf")

	
	DataFilelist<-params$DataFiles
	delaytime<-params$delaytime
	ScanInterval<-params$ScanInterval
	WorkDir<-params$WorkDir
	libFile <- params$libFile_decon
#	StN_Th <- as.integer(params$StN_Th2)
#	Shrp_Th <- as.numeric(params$Sharpness_Th)
#	Gss_Th <- as.numeric(params$Gss_Th)	
#	TotalScore_factor <- as.numeric(params$TotalScore_factor)
	ClusterDist_Th <- as.numeric(params$ClusteringDist_Th)
	SpecSimilarity_Th <- as.numeric(params$IntroSpecSimilarity_Th)
	libMatch_Th <- as.numeric(params$libMatching_Th1)
	
	fileName<-parseFileName(inFilePath)
	
	#read TIC data
	TICfile<-paste(WorkDir,"output/TIC/denoised_",fileName,"_TIC.cdf",sep="") 
	TICfile<-ifelse(file.exists(TICfile),TICfile,inFilePath)
	ncid <- open.ncdf(TICfile)
	TIC <- get.var.ncdf(ncid, varid="total_intensity")
	close.ncdf(ncid)
	remove(ncid)
	
	# read EIC data
	rawEICFile<-paste(WorkDir,"output/EIC/",fileName,"EIC.cdf",sep="")
	denoiseEICFile<-paste(WorkDir,"output/EIC/denoised_",fileName,"EIC.cdf",sep="")
	EICFile<-ifelse(file.exists(denoiseEICFile),denoiseEICFile,rawEICFile)
	if(file.exists(EICFile)==FALSE)next
	ncid <- open.ncdf(EICFile)
	vecInt <- get.var.ncdf(ncid, varid="intVec")
	vectorMz<-get.var.ncdf(ncid, varid="mzVec")
	close.ncdf(ncid)
	remove(ncid)	
	
	# get decon window from TIC peak picking results
	cat(fileName,"reading TIC peak data...\n")
	TICPeakFile<-paste(params$WorkDir,"output/peakpicking/",fileName,"_TIC_PeakList.csv",sep="")
	TICpeaklist <- read.csv(TICPeakFile)
	TICpeaklist$RT<-(TICpeaklist$pkInd-1)*ScanInterval+delaytime
	TICpeaklist$window <- 0
	TICpeakIDs <- paste(TICpeaklist$lboundInd,TICpeaklist$rboundInd)
	uniqueID <- unique(TICpeakIDs)
	for (i in 1:length(uniqueID)) {
		curPos <- which(TICpeakIDs%in%uniqueID[i])
		TICpeaklist[curPos,]$window <- i
	}
	winIDs<-1:length(uniqueID)
	
	# read EIC peak picking result
	cat(fileName,"reading EIC peak data...\n")
	EICPeakFile<-paste(WorkDir,"/output/peakpicking/",fileName,"_EIC_PeakList.csv",sep="")
	if(file.exists(EICPeakFile)==FALSE)next
	EICpeaklist <- read.csv(EICPeakFile)
	#colnames(EICpeaklist)[which(colnames(EICpeaklist)=="curMz")]<-"mz"
	
	
	cat(fileName,"start deconvolution...\n")
	
	componentResults<-NULL
	specResults<-NULL
	mdlPkIDVec<-NULL
	#EICpeaklist$flag=0
	totalscan<-length(vecInt)/length(vectorMz)
	strOutput <- NULL
	compNum <- 0
	#outputPath<-paste(paste(params$WorkDir,"output/decon/",sep=""),paste(fileName,"_ModelPkProfiles.txt",sep=""),sep="")
	
	#for(windowID in c(1:29))
	for(windowID in winIDs)
	{
		cat("window:",windowID,"\n")
		curTICPk<-subset(TICpeaklist,window==windowID)[1,]
		minET=(curTICPk$lboundInd-1)*params$ScanInterval+params$delaytime
		maxET=(curTICPk$rboundInd-1)*params$ScanInterval+params$delaytime
		
		if(minET>=8)
		{
			CurNonUmassVec<-c(1:100,147,221)
		}else
		{
			CurNonUmassVec<-c(1:73,147,221)
		}
		
		# get all EIC peaks within current decon window
		localEICPeakList<-subset(EICpeaklist,pkInd>=curTICPk$lboundInd&pkInd<=curTICPk$rboundInd&lboundInd>0)
		#localEICPeakList <- subset(localEICPeakList,Intensity>=200&StN>=20) # done in peak picking
		#localEICPeakList<-subset(EICpeaklist,pkInd>=curTICPk$lboundInd&pkInd<=curTICPk$rboundInd)
	
		if (nrow(localEICPeakList)>=5) {
			
			# step 1: to remove noisy peaks within decon window
			clustResult1 <- hclust(dist(localEICPeakList$pkInd))
			#clustResult1 <- hclust(dist(localEICPeakList$pkInd),method="average")
			clustResult2 <- cutree(clustResult1,h=60)
			localEICPeakList$clust <- clustResult2
			
			# remove noisy peaks clustering less than two peaks
			noisyPeaks <- which(clustResult2%in%which(table(clustResult2)<=2)==TRUE)
			# remove noisy peaks: small cluster with no significant mass
			for (i in unique(clustResult2)) {
				curCluster <- subset(localEICPeakList,clust==i&Intensity>=500)
				if (nrow(curCluster)==0) {
					noisyPeaks <- c(noisyPeaks,which(localEICPeakList$clust==i))
				}
			}
			
			noisyPeaks <- unique(noisyPeaks)
			# clean local EIC peaks
			if (length(noisyPeaks)>0) {
				localEICPeakList <- localEICPeakList[-noisyPeaks,]
			}
			
			
			if (nrow(localEICPeakList)>=5) {
				
				uniqueClust <- unique(localEICPeakList$clust)
				
				# step 2: get all EIC profiles
				allprofiles<-vector("list",nrow(localEICPeakList))
				for(i in 1:nrow(localEICPeakList))
				{
					curEICpk<-localEICPeakList[i,]			
					startInd<-curEICpk$offset+curEICpk$Lbound
					endInd<-curEICpk$offset+curEICpk$Rbound
					ET<-curEICpk$Lbound:curEICpk$Rbound
					curProfile<-data.frame(ET=ET,int=vecInt[startInd:endInd]) 
					allprofiles[[i]]<-curProfile
					names(allprofiles)[i]<-rownames(curEICpk)
				}	
				
				# step 3: get qualified unique mass feature for 2nd clustering and model peak selection
				uniquePeakIDs <- which(localEICPeakList$isShared==0&localEICPeakList$StN>=100&localEICPeakList$shrp>=40)
				
				if (length(uniquePeakIDs)>0) {
					modelPkList <- NULL				
					clust=0					
					uniquePeakList <- localEICPeakList[uniquePeakIDs,]
					
					for (i in uniqueClust) {
						# get qualified unique mass for clustering and model peak selection
						curClusterIDs <- which(uniquePeakList$clust==i)	
						#print(length(curClusterIDs))
						if (length(curClusterIDs)>0) {
							uniqueEICPeakList <- uniquePeakList[curClusterIDs,]
							uniqueEICPeakList$CompoundID <- 0
							
							if (length(curClusterIDs)==1) {
								clust <- clust + 1
								uniqueEICPeakList$CompoundID <- clust
								
								if (length(uniqueClust)==1) {
									# no limitations to the model peak if there's only one compound 
									modelPkList <- rbind(modelPkList,uniqueEICPeakList)
								}else {
									if(!uniqueEICPeakList$mz%in%CurNonUmassVec) {
										modelPkList <- rbind(modelPkList,uniqueEICPeakList)
									}
								}
								
							}else{
								curProfiles <- allprofiles[which(rownames(localEICPeakList)%in%rownames(uniqueEICPeakList))]			
								assign("Global.curProfs", value=curProfiles, envir = .GlobalEnv)			
								r<-DistCal4(isUnion=F)
								distance <- as.dist(r)
								#clustResut<-hclust(distance,method="average")
								clustResut<-hclust(distance)
								FinalClustResult<-cutree(clustResut,h=17)
								Clusters<-unique(FinalClustResult)	
								uniqueEICPeakList$CompoundID <- clust+FinalClustResult		
								clust <- clust + length(Clusters)
								
								for(curCluster in Clusters)
								{
									curUniqueEICPeakList <- uniqueEICPeakList[which(FinalClustResult==curCluster),]
									if (length(uniqueClust)==1&length(Clusters)==1) {
										scaledShrp <- curUniqueEICPeakList$shrp
										modelPkList <- rbind(modelPkList,curUniqueEICPeakList[which.max(scaledShrp),])
									}else {
										FilteredCurUniqueEICPeakList <- curUniqueEICPeakList[!curUniqueEICPeakList$mz%in%CurNonUmassVec,]
										if (nrow(FilteredCurUniqueEICPeakList)>0) {
											scaledShrp <- FilteredCurUniqueEICPeakList$shrp
											modelPkList <- rbind(modelPkList,FilteredCurUniqueEICPeakList[which.max(scaledShrp),])
										}
									}
									
								}
							}
						}	
					} # for loop

					
					uniqueClust <- uniqueClust[-which(uniqueClust%in%unique(uniquePeakList$clust)==TRUE)]
					
					# Step 4: decomposition on local EIC mass 
					if (!is.null(modelPkList)) { 
						decom = 1
						while (decom==1) {
							specList<-NULL
							componentList<-NULL
							mdlPkIDs<-rownames(modelPkList)							
							
							localEICPeakList.noModels <- localEICPeakList[-which(rownames(localEICPeakList)%in%rownames(modelPkList)),]
							
							# step 4.1: get all local EIC profiles within the window, X matrix format:column-EIC,row-scan
							lbound<-min(localEICPeakList.noModels$Lbound)
							rbound<-max(localEICPeakList.noModels$Rbound)
							mzVec<-unique(localEICPeakList.noModels$mz)
							
#							lbound<-min(localEICPeakList)
#							rbound<-max(localEICPeakList)
#							mzVec<-unique(localEICPeakList$mz)
							nMz<-length(mzVec)
							X<-NULL
							
							for(i in 1:nMz)
							{
								mz<-mzVec[i]
								mzInd<-which(vectorMz==mz)
								startInd<-(mzInd-1)*totalscan+1
								endInd<-startInd+totalscan-1
								curVecInt <- vecInt[startInd:endInd]
								X<-cbind(X,curVecInt[lbound:rbound])
							}
							colnames(X)<-mzVec	
							nScans<-nrow(X)
							
#							# testing *********************
#							mzVec1<-unique(localEICPeakList$mz)
#							nMz1<-length(mzVec1)
#							XX<-NULL
#							extracted_specList <- NULL
#							
#							for(i in 1:nMz1)
#							{
#								mz<-mzVec1[i]
#								mzInd<-which(vectorMz==mz)
#								startInd<-(mzInd-1)*totalscan+1
#								endInd<-startInd+totalscan-1
#								curVecInt <- vecInt[startInd:endInd]
#								XX<-cbind(XX,curVecInt[lbound:rbound])
#							}
#							
#							colnames(XX)<-mzVec1	
#							nScans<-nrow(XX)
#							
#							if (nrow(modelPkList)>1) {
#								for (i in 1:nrow(modelPkList)) {
#									M1 <- XX[which(c(lbound:rbound)==modelPkList[i,]$pkInd),]
#									extracted_specList[[i]] <- data.frame(mz=as.numeric(names(M1)),int=as.vector(M1))
#									#specList[[3]] <- rbind(specList[[3]],data.frame(mz=modelPkList[1,]$mz,int=modelPkList[1,]$Intensity))	
#								}
#							}
#							
#
#						
#							# testing *********************

							# step 4.2: Get model peak profiles as S matrix, each column is the EIC for one model peak	
							# add model mass information to spec and component list
							S<-NULL
							for(i in 1:nrow(modelPkList))
							{
								# model peak features are unique ones, directly get profiles from allprofiles
								curProfile<-merge(data.frame(ET=as.integer(lbound:rbound)),allprofiles[[rownames(modelPkList[i,])]],by.y="ET",all.x=T)
								curProfile$int[is.na(curProfile$int)]<-0 # expand the profile by filling zero values to those uncovered area
								S<-cbind(S,curProfile$int)	
								
								# init the component & spec list	
								specList[[i]]<-data.frame(mz=as.integer(),int=as.integer())
								names(specList)[i]<-paste(modelPkList$pkInd[i],windowID,i,modelPkList$mz[i],compNum+i)
								componentList[[i]]<-data.frame(mz=as.integer(),pkInd=as.integer(),Intensity=as.integer(),lboundInd=as.integer(),rboundInd=as.integer(),
										area=as.integer(),windowID=as.integer(),compoundID=as.integer(),isModel=as.integer(),isUnique=as.integer(),factor=as.numeric(),compID=as.integer())
								
								# add model peak mass information to spec list and component list
								curModelPeak <- modelPkList[i,]
								specList[[i]] <- rbind(specList[[i]],data.frame(mz=curModelPeak$mz,int=curModelPeak$Intensity))
								
								curModelPeak$windowID <- windowID
								curModelPeak$compoundID <- i
								curModelPeak$isModel <-	1	
								curModelPeak$isUnique <- 1 # keep for normalization format only, will be removed
								curModelPeak$factor <- 0
								curModelPeak$compID <-compNum+i				
								componentList[[i]]<-rbind(componentList[[i]],subset(curModelPeak,select=c("mz","pkInd","Intensity","lboundInd","rboundInd","area",
														"windowID","compoundID","isModel","isUnique","factor","compID")))
							}		
							colnames(S)<-mdlPkIDs
							
							# step 4.3: Residual minimization
							for(i in 1:ncol(X))
							{
								curMz<-colnames(X)[i]
								# to indicate unique or not, will be removed later
								isShared <- localEICPeakList[which(localEICPeakList$mz==as.integer(curMz))[1],]$isShared
								isUnique <- ifelse(isShared==0,1,0)
								
								M<-X[,i]
								A<-optim(par=rep(0,ncol(S)),fn=fr, gr = NULL,M,S,lower = 0, method="L-BFGS-B")$par
								A<-as.matrix(A)
								
								# restore each ion signals and save them into component list ans speclist
								
								for(m in 1:length(mdlPkIDs))
								{					
									curPkInd <- modelPkList[m,]$pkInd
									curSharedLbound <- modelPkList[m,]$lboundInd
									curSharedRbound <- modelPkList[m,]$rboundInd
									
									if(A[m]>0)
									{
										curIntensity <- as.integer(modelPkList[m,]$Intensity*A[m,])
										curArea <- as.integer(modelPkList[m,]$area*A[m,])
										
										specList[[m]]<-rbind(specList[[m]],data.frame(mz=as.integer(curMz),int=curIntensity))
										componentList[[m]]<-rbind(componentList[[m]],data.frame(mz=as.integer(curMz),pkInd=curPkInd,lboundInd=curSharedLbound,rboundInd=curSharedRbound,
														Intensity=curIntensity,area=curArea,windowID=windowID,compoundID=m,isModel=0,isUnique=isUnique,factor=A[m],compID=compNum+m))						
									}							
								}
							}
							
							# step 4.4: Component splitting correction: checking pairwise resovlved spectra similarity
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
										score <- specDistCal_check(specList[[index1]],specList[[index2]],isWeight=T,isNist=T) 
										#score2 <-specDistCal(extracted_specList[[index1]],extracted_specList[[index2]],isWeight=T,isNist=T)
										print(paste(index1,index2,score,0.9*score+(1-(abs(modelPKET[1]-modelPKET[2])/10))*100))
										if (0.9*score+(1-(abs(modelPKET[1]-modelPKET[2])/20))*100 >=750) {#SpecSimilarity_Th
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
								print(paste(compNum,Num,sep="  "))
								
#						## output model pk profiles
#						for (model in c(1:nrow(modelPkList))) {
#							curModel <- modelPkList[model,]
#							curModelInt <- S[,model]
#							strOutput <- paste(strOutput,getModelPkProfile(curModel,curModelInt,windowID,lbound,rbound,model,compNum),sep="\n\n")
#						}
							}
						}# decom while statement
					} # check model peak 
				} # check unique peaks for decon
				
				if (length(uniqueClust)>0) {
					
					if (is.null(modelPkList)) { # initialization if no component detected before
						specList<-NULL
						componentList<-NULL
					}
					
					CompoundNum <- ifelse(is.null(modelPkList),0,nrow(modelPkList))
					
					for (i in uniqueClust) {
						curCompPkInd <- round(median(subset(localEICPeakList,clust==i)$pkInd))
						curLeftbound <- round(median(subset(localEICPeakList,clust==i)$Lbound))
						curRightbound <- round(median(subset(localEICPeakList,clust==i)$Rbound))
						
						if (curLeftbound < curCompPkInd & curCompPkInd < curRightbound) {
							
							# initialization
							curSpecList<-data.frame(mz=as.integer(),int=as.integer())
							curComponentList<-data.frame(mz=as.integer(),pkInd=as.integer(),Intensity=as.integer(),lboundInd=as.integer(),rboundInd=as.integer(),
									area=as.integer(),windowID=as.integer(),compoundID=as.integer(),isModel=as.integer(),isUnique=as.integer(),factor=as.numeric(),compID=as.integer())
							
							# extract mass spectrum
							for (m in 1:length(allprofiles)) {
								curProfile <- allprofiles[[m]]
								
								curArea <- sum(subset(curProfile,ET>=curLeftbound&ET<=curRightbound)$int)	
								curPk <- subset(curProfile,ET==curCompPkInd)
								if (nrow(curPk)==1) {
									curIntensity <- as.integer(curPk$int)
									curMz <- as.integer(localEICPeakList$mz[m])
									
									curSpecList<-rbind(curSpecList,data.frame(mz=as.integer(curMz),int=curIntensity))
									curComponentList<-rbind(curComponentList,data.frame(mz=as.integer(curMz),pkInd=curPkInd,lboundInd=curLeftbound,rboundInd=curRightbound,
													Intensity=curIntensity,area=curArea,windowID=windowID,compoundID=CompoundNum+1,isModel=0,isUnique=0,factor=0,compID=compNum+1))	
								}
								
							}
							
							# check quality and import spectrum and component into list
							if (nrow(curSpecList) >5) {
								mostIntenseMassPos <- which.max(curComponentList$Intensity)
								curComponentList[mostIntenseMassPos,]$isModel <- 1
								curComponentList[mostIntenseMassPos,]$isUnique <- 1
								curModelmass <- curComponentList[mostIntenseMassPos,]$mz
								specList[[CompoundNum+1]] <- curSpecList
								names(specList)[CompoundNum+1]<-paste(curCompPkInd,windowID,CompoundNum+1,curModelmass,compNum+1)
								componentList[[CompoundNum+1]]<-curComponentList
								CompoundNum <- CompoundNum+1
								compNum <- compNum +1
								print(paste(compNum,"X",sep="  "))
							}
						} # check qualify
					}	# for loop
				} # end check missing compounds
				
				componentResults<-c(componentResults,componentList)
				specResults<-c(specResults,specList)	
				
			}	# Check localEICPeakList 2
		}# Check localEICPeakList 1
	}#for loop 	
	
	# step 5: output
	cat(fileName,"finished deconvolution...\n")
	lib<-readMSP2Spec(filename=libFile,withRT=F)
	libmatch <- libMatching_Top10(lib,specResults,libMatch_Th) # target lib matching 
	libmatch$ET <- (libmatch$pkind-1)*params$ScanInterval+params$delaytime
	
	componentResults<-do.call(rbind,componentResults)	
	componentResults<-subset(componentResults,Intensity>0,select=c("mz","pkInd","lboundInd","rboundInd","Intensity","area","windowID","compoundID","isModel","isUnique","factor","compID"))
	write.csv(componentResults,row.names=FALSE,file=paste(params$WorkDir,"/output/decon/",fileName,"Decon.csv",sep=""))
	write.csv(subset(libmatch,select=c("ET","modelMass","compoundName","score")),file=paste(params$WorkDir,"/output/decon/",fileName,"_libmatch.csv",sep=""))
	#write(strOutput,file=outputPath)
	
	cat(fileName,"All deconvolution results are exported successfully!\n")
}

decomposition_ncdf4_test2<-function(inFilePath,params,cl,isDistParallel,clustingType)
{
	library(gdata)
	library(gtools)
	library(cluster)
	library("ncdf4")	
	#library("ncdf")
	#library(mclust)
	#options(expressions=1e5)
	
	DataFilelist<-params$DataFiles
	delaytime<-params$delaytime
	ScanInterval<-params$ScanInterval
	WorkDir<-params$WorkDir
	libFile <- params$libFile_decon
	StN_Th <- as.integer(params$StN_Th2)
	Shrp_Th <- as.numeric(params$Sharpness_Th)
	Gss_Th <- as.numeric(params$Gss_Th)	
	TotalScore_factor <- as.numeric(params$TotalScore_factor)
	ClusterDist_Th <- as.numeric(params$ClusteringDist_Th)
	SpecSimilarity_Th <- as.numeric(params$IntroSpecSimilarity_Th)
	libMatch_Th <- as.numeric(params$libMatching_Th1)
	
	
	fileName<-parseFileName(inFilePath)
	
	# read TIC data
#	TICfile<-paste(WorkDir,"output/TIC/denoised_",fileName,"_TIC.cdf",sep="") 
#	TICfile<-ifelse(file.exists(TICfile),TICfile,inFilePath)
#	ncid <- open.ncdf(TICfile)
#	TIC <- get.var.ncdf(ncid, varid="total_intensity")
#	close.ncdf(ncid)
#	remove(ncid)
#	
#	# read EIC data
#	rawEICFile<-paste(WorkDir,"output/EIC/",fileName,"EIC.cdf",sep="")
#	denoiseEICFile<-paste(WorkDir,"output/EIC/denoised_",fileName,"EIC.cdf",sep="")
#	EICFile<-ifelse(file.exists(denoiseEICFile),denoiseEICFile,rawEICFile)
#	if(file.exists(EICFile)==FALSE)next
#	ncid <- open.ncdf(EICFile)
#	vecInt <- get.var.ncdf(ncid, varid="intVec")
#	vectorMz<-get.var.ncdf(ncid, varid="mzVec")
#	close.ncdf(ncid)
#	remove(ncid)	
	
#	# read TIC data
	TICfile<-paste(WorkDir,"output/TIC/denoised_",fileName,"_TIC.cdf",sep="") 
	TICfile<-ifelse(file.exists(TICfile),TICfile,inFilePath)
	ncid <- nc_open(TICfile)
	TIC <- ncvar_get(ncid, varid="total_intensity" )
	nc_close(ncid)
	remove(ncid)
	
	# read EIC data
	rawEICFile<-paste(WorkDir,"output/EIC/",fileName,"EIC.cdf",sep="")
	denoiseEICFile<-paste(WorkDir,"output/EIC/denoised_",fileName,"EIC.cdf",sep="")
	EICFile<-ifelse(file.exists(denoiseEICFile),denoiseEICFile,rawEICFile)
	if(file.exists(EICFile)==FALSE)next
	ncid <- nc_open(EICFile)
	vecInt <- ncvar_get(ncid, varid="intVec")
	vectorMz<-ncvar_get(ncid, varid="mzVec")
	nc_close(ncid)
	remove(ncid)	
	
	# get decon window from TIC peak picking results
	cat(fileName,"reading TIC peak data...\n")
	TICPeakFile<-paste(params$WorkDir,"output/peakpicking/",fileName,"_TIC_PeakList.csv",sep="")
	TICpeaklist <- read.csv(TICPeakFile)
	TICpeaklist$RT<-(TICpeaklist$pkInd-1)*ScanInterval+delaytime
	TICpeaklist$window <- 0
	TICpeakIDs <- paste(TICpeaklist$lboundInd,TICpeaklist$rboundInd)
	uniqueID <- unique(TICpeakIDs)
	for (i in 1:length(uniqueID)) {
		curPos <- which(TICpeakIDs%in%uniqueID[i])
		TICpeaklist[curPos,]$window <- i
	}
	#winIDs<-1:length(uniqueID)
	
	# for target analysis of standards
	RTofInterest <- c(5.17,5.34,7.47,8.39,8.73,8.78,9.34,10.31,12.8,13.56,
			14.84,15.93,16.14,16.16,17.59,18.51,18.95,19.81,19.87,21.92,21.96,
			22.61,22.87,25.97,27.94,31.38,32.31)
	
	TICofInterest <- NULL
	LBET=(TICpeaklist$lboundInd-1)*params$ScanInterval+params$delaytime
	RBET=(TICpeaklist$rboundInd-1)*params$ScanInterval+params$delaytime
	for (i in 1:length(RTofInterest)) {
		curComponent <- RTofInterest[i]
		#TICofInterest <- rbind(TICofInterest,TICpeaklist[which.min(abs(TICpeaklist$RT-curComponent)),])
		TICofInterest <- rbind(TICofInterest,TICpeaklist[which(LBET<=curComponent&RBET>=curComponent),])
	}
	winIDs <- unique(TICofInterest$window)
	
	
	# read EIC peak picking result
	cat(fileName,"reading EIC peak data...\n")
	EICPeakFile<-paste(WorkDir,"/output/peakpicking/",fileName,"_EIC_PeakList.csv",sep="")
	if(file.exists(EICPeakFile)==FALSE)next
	EICpeaklist <- read.csv(EICPeakFile)
	#colnames(EICpeaklist)[which(colnames(EICpeaklist)=="curMz")]<-"mz"
	
	
	cat(fileName,"start deconvolution...\n")
	
	componentResults<-NULL
	specResults<-NULL
	mdlPkIDVec<-NULL
	#EICpeaklist$flag=0
	totalscan<-length(vecInt)/length(vectorMz)
	strOutput <- NULL
	compNum <- 0
	#outputPath<-paste(paste(params$WorkDir,"output/decon/",sep=""),paste(fileName,"_ModelPkProfiles.txt",sep=""),sep="")
	
	#for(windowID in c(1:29))
	for(windowID in winIDs)
	{

		cat("window:",windowID,"\n")
		curTICPk<-subset(TICpeaklist,window==windowID)[1,]
		minET=(curTICPk$lboundInd-1)*params$ScanInterval+params$delaytime
		maxET=(curTICPk$rboundInd-1)*params$ScanInterval+params$delaytime
		
		if(minET>=8)
		{
			CurNonUmassVec<-c(1:100,147,221)
		}else
		{
			CurNonUmassVec<-c(1:73,147,221)
		}
		
		# get all EIC peaks within current decon window
		localEICPeakList<-subset(EICpeaklist,pkInd>=curTICPk$lboundInd&pkInd<=curTICPk$rboundInd&lboundInd>0)
		#localEICPeakList <- subset(localEICPeakList,Intensity>=200&StN>=20) # done in peak picking
		#localEICPeakList<-subset(EICpeaklist,pkInd>=curTICPk$lboundInd&pkInd<=curTICPk$rboundInd)
		
		if (nrow(localEICPeakList)>=5) {
			noisyPeaks <- NULL
			# step 1: to remove noisy peaks within decon window
			clustResult1 <- hclust(dist(localEICPeakList$pkInd))
			#clustResult1 <- hclust(dist(localEICPeakList$pkInd),method="average")
			clustResult2 <- cutree(clustResult1,h=60)
			localEICPeakList$clust <- clustResult2
			
			# remove noisy peaks clustering less than two peaks
			#noisyPeaks <- which(clustResult2%in%which(table(clustResult2)<=2)==TRUE)
			# remove noisy peaks: small cluster with no significant mass
			for (i in unique(clustResult2)) {
#				curCluster <- subset(localEICPeakList,clust==i&Intensity>=500)
#				if (nrow(curCluster)==0) {
#					noisyPeaks <- c(noisyPeaks,which(localEICPeakList$clust==i))
#				}
				
				curCluster <- subset(localEICPeakList,clust==i)
				if (nrow(curCluster)<=2 & length(which(curCluster$Intensity>=500))==0) {
					
					noisyPeaks <- c(noisyPeaks,which(localEICPeakList$clust==i))
				}
				
			}
			
			noisyPeaks <- unique(noisyPeaks)
			# clean local EIC peaks
			if (length(noisyPeaks)>0) {
				localEICPeakList <- localEICPeakList[-noisyPeaks,]
			}
			
			
			if (nrow(localEICPeakList)>=5) {
				
				uniqueClust <- unique(localEICPeakList$clust)
				
				# step 2: get all EIC profiles
				allprofiles<-vector("list",nrow(localEICPeakList))
				for(i in 1:nrow(localEICPeakList))
				{
					curEICpk<-localEICPeakList[i,]			
					startInd<-curEICpk$offset+curEICpk$Lbound
					endInd<-curEICpk$offset+curEICpk$Rbound
					ET<-curEICpk$Lbound:curEICpk$Rbound
					curProfile<-data.frame(ET=ET,int=vecInt[startInd:endInd]) 
					allprofiles[[i]]<-curProfile
					names(allprofiles)[i]<-rownames(curEICpk)
				}	
				
				# step 3: get qualified unique mass feature for 2nd clustering and model peak selection
				#uniquePeakIDs <- which(localEICPeakList$isShared==0&(localEICPeakList$StN>=100|localEICPeakList$shrp>=10))
#				if (length(uniqueClust)>1) {
#					NonUniqueMassCheck <- 1
#					uniquePeakIDs <- which(localEICPeakList$isShared==0&localEICPeakList$StN>=100&localEICPeakList$shrp>=10&!localEICPeakList$mz%in%CurNonUmassVec)
#				}else {
#					NonUniqueMassCheck <- 0
				#shrpnessCutoff <- ifelse(diff(range(localEICPeakList$shrp,na.rm=TRUE))*0.01>50,diff(range(localEICPeakList$shrp,na.rm=TRUE))*0.01,50)
				#StNcutoff <- diff(range(localEICPeakList$StN))*0.1	
				#uniquePeakIDs <- which(localEICPeakList$isShared==0&localEICPeakList$StN>=StNcutoff&localEICPeakList$shrp>=shrpnessCutoff)
				uniquePeakIDs <- which(localEICPeakList$isShared==0&localEICPeakList$StN>=100&localEICPeakList$shrp>=10)
				
#				}
				
				if (length(uniquePeakIDs)>0) {
					modelPkList <- NULL				
					clust=0					
					uniquePeakList <- localEICPeakList[uniquePeakIDs,]
					
					for (i in uniqueClust) {
						# get qualified unique mass for clustering and model peak selection
						curClusterIDs <- which(uniquePeakList$clust==i)	
						#print(length(curClusterIDs))
						if (length(curClusterIDs)>0) {
							uniqueEICPeakList <- uniquePeakList[curClusterIDs,]
							uniqueEICPeakList$CompoundID <- 0
							
							if (length(curClusterIDs)==1) {
								clust <- clust + 1
								uniqueEICPeakList$CompoundID <- clust
								
								if (length(uniqueClust)==1) {
									# no limitations to the model peak if there's only one compound 
									modelPkList <- rbind(modelPkList,uniqueEICPeakList)
								}else {
									if(!uniqueEICPeakList$mz%in%CurNonUmassVec) {
										modelPkList <- rbind(modelPkList,uniqueEICPeakList)
									}
								}
								
							}else{
								curProfiles <- allprofiles[which(rownames(localEICPeakList)%in%rownames(uniqueEICPeakList))]			
								assign("Global.curProfs", value=curProfiles, envir = .GlobalEnv)			
								r<-DistCal4(isUnion=F)
								distance <- as.dist(r)
								#clustResut<-hclust(distance,method="average")
								clustResut<-hclust(distance)
								FinalClustResult<-cutree(clustResut,h=17)
								Clusters<-unique(FinalClustResult)	
								uniqueEICPeakList$CompoundID <- clust+FinalClustResult		
								clust <- clust + length(Clusters)
								
								for(curCluster in Clusters)
								{
									curUniqueEICPeakList <- uniqueEICPeakList[which(FinalClustResult==curCluster),]
									if (length(uniqueClust)==1&length(Clusters)==1) {
										scaledShrp <- curUniqueEICPeakList$shrp
										modelPkList <- rbind(modelPkList,curUniqueEICPeakList[which.max(scaledShrp),])
									}else {
										FilteredCurUniqueEICPeakList <- curUniqueEICPeakList[!curUniqueEICPeakList$mz%in%CurNonUmassVec,]
										if (nrow(FilteredCurUniqueEICPeakList)>0) {
											scaledShrp <- FilteredCurUniqueEICPeakList$shrp
											modelPkList <- rbind(modelPkList,FilteredCurUniqueEICPeakList[which.max(scaledShrp),])
										}
									}
									
								}
							}
						}	
					} # for loop
					
					
					uniqueClust <- uniqueClust[-which(uniqueClust%in%unique(uniquePeakList$clust)==TRUE)]
					
					# Step 4: decomposition on local EIC mass 
					if (!is.null(modelPkList)) { 
						decom = 1
						while (decom==1) {
							specList<-NULL
							componentList<-NULL
							mdlPkIDs<-rownames(modelPkList)							
							
							localEICPeakList.noModels <- localEICPeakList[-which(rownames(localEICPeakList)%in%rownames(modelPkList)),]
							
							# step 4.1: get all local EIC profiles within the window, X matrix format:column-EIC,row-scan
							lbound<-min(localEICPeakList.noModels$Lbound)
							rbound<-max(localEICPeakList.noModels$Rbound)
							mzVec<-unique(localEICPeakList.noModels$mz)
							
#							lbound<-min(localEICPeakList)
#							rbound<-max(localEICPeakList)
#							mzVec<-unique(localEICPeakList$mz)
							nMz<-length(mzVec)
							X<-NULL
							
							for(i in 1:nMz)
							{
								mz<-mzVec[i]
								mzInd<-which(vectorMz==mz)
								startInd<-(mzInd-1)*totalscan+1
								endInd<-startInd+totalscan-1
								curVecInt <- vecInt[startInd:endInd]
								X<-cbind(X,curVecInt[lbound:rbound])
							}
							colnames(X)<-mzVec	
							nScans<-nrow(X)
							
#							# testing *********************
#							mzVec1<-unique(localEICPeakList$mz)
#							nMz1<-length(mzVec1)
#							XX<-NULL
#							extracted_specList <- NULL
#							
#							for(i in 1:nMz1)
#							{
#								mz<-mzVec1[i]
#								mzInd<-which(vectorMz==mz)
#								startInd<-(mzInd-1)*totalscan+1
#								endInd<-startInd+totalscan-1
#								curVecInt <- vecInt[startInd:endInd]
#								XX<-cbind(XX,curVecInt[lbound:rbound])
#							}
#							
#							colnames(XX)<-mzVec1	
#							nScans<-nrow(XX)
#							
#							if (nrow(modelPkList)>1) {
#								for (i in 1:nrow(modelPkList)) {
#									M1 <- XX[which(c(lbound:rbound)==modelPkList[i,]$pkInd),]
#									extracted_specList[[i]] <- data.frame(mz=as.numeric(names(M1)),int=as.vector(M1))
#									#specList[[3]] <- rbind(specList[[3]],data.frame(mz=modelPkList[1,]$mz,int=modelPkList[1,]$Intensity))	
#								}
#							}
#							
#
#						
#							# testing *********************
							
							# step 4.2: Get model peak profiles as S matrix, each column is the EIC for one model peak	
							# add model mass information to spec and component list
							S<-NULL
							for(i in 1:nrow(modelPkList))
							{
								# model peak features are unique ones, directly get profiles from allprofiles
								curProfile<-merge(data.frame(ET=as.integer(lbound:rbound)),allprofiles[[rownames(modelPkList[i,])]],by.y="ET",all.x=T)
								curProfile$int[is.na(curProfile$int)]<-0 # expand the profile by filling zero values to those uncovered area
								S<-cbind(S,curProfile$int)	
								
								# init the component & spec list	
								specList[[i]]<-data.frame(mz=as.integer(),int=as.integer())
								names(specList)[i]<-paste(modelPkList$pkInd[i],windowID,i,modelPkList$mz[i],compNum+i)
								componentList[[i]]<-data.frame(mz=as.integer(),pkInd=as.integer(),Intensity=as.integer(),lboundInd=as.integer(),rboundInd=as.integer(),
										area=as.integer(),windowID=as.integer(),compoundID=as.integer(),isModel=as.integer(),isUnique=as.integer(),factor=as.numeric(),compID=as.integer())
								
								# add model peak mass information to spec list and component list
								curModelPeak <- modelPkList[i,]
								specList[[i]] <- rbind(specList[[i]],data.frame(mz=curModelPeak$mz,int=curModelPeak$Intensity))
								
								curModelPeak$windowID <- windowID
								curModelPeak$compoundID <- i
								curModelPeak$isModel <-	1	
								curModelPeak$isUnique <- 1 # keep for normalization format only, will be removed
								curModelPeak$factor <- 0
								curModelPeak$compID <-compNum+i				
								componentList[[i]]<-rbind(componentList[[i]],subset(curModelPeak,select=c("mz","pkInd","Intensity","lboundInd","rboundInd","area",
														"windowID","compoundID","isModel","isUnique","factor","compID")))
							}		
							colnames(S)<-mdlPkIDs
							
							# step 4.3: Residual minimization
							for(i in 1:ncol(X))
							{
								curMz<-colnames(X)[i]
								# to indicate unique or not, will be removed later
								isShared <- localEICPeakList[which(localEICPeakList$mz==as.integer(curMz))[1],]$isShared
								isUnique <- ifelse(isShared==0,1,0)
								
								M<-X[,i]
								A<-optim(par=rep(0,ncol(S)),fn=fr, gr = NULL,M,S,lower = 0, method="L-BFGS-B")$par
								A<-as.matrix(A)
								
								# restore each ion signals and save them into component list ans speclist
								
								for(m in 1:length(mdlPkIDs))
								{					
									curPkInd <- modelPkList[m,]$pkInd
									curSharedLbound <- modelPkList[m,]$lboundInd
									curSharedRbound <- modelPkList[m,]$rboundInd
									
									if(A[m]>0)
									{
										curIntensity <- as.integer(modelPkList[m,]$Intensity*A[m,])
										curArea <- as.integer(modelPkList[m,]$area*A[m,])
										
										specList[[m]]<-rbind(specList[[m]],data.frame(mz=as.integer(curMz),int=curIntensity))
										componentList[[m]]<-rbind(componentList[[m]],data.frame(mz=as.integer(curMz),pkInd=curPkInd,lboundInd=curSharedLbound,rboundInd=curSharedRbound,
														Intensity=curIntensity,area=curArea,windowID=windowID,compoundID=m,isModel=0,isUnique=isUnique,factor=A[m],compID=compNum+m))						
									}							
								}
							}
							
							# step 4.4: Component splitting correction: checking pairwise resovlved spectra similarity
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
										score <- specDistCal_check(specList[[index1]],specList[[index2]],isWeight=T,isNist=T) 
										#score2 <-specDistCal(extracted_specList[[index1]],extracted_specList[[index2]],isWeight=T,isNist=T)
										print(paste(index1,index2,score,0.9*score+(1-(abs(modelPKET[1]-modelPKET[2])/10))*100))
										if (0.9*score+(1-(abs(modelPKET[1]-modelPKET[2])/20))*100 >=750) {#SpecSimilarity_Th
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
								print(paste(compNum,Num,sep="  "))
								
#						## output model pk profiles
#						for (model in c(1:nrow(modelPkList))) {
#							curModel <- modelPkList[model,]
#							curModelInt <- S[,model]
#							strOutput <- paste(strOutput,getModelPkProfile(curModel,curModelInt,windowID,lbound,rbound,model,compNum),sep="\n\n")
#						}
							}
						}# decom while statement
					} # check model peak 
				} # check unique peaks for decon
				
				if (length(uniqueClust)>0) {
					
					if (is.null(modelPkList)) { # initialization if no component detected before
						specList<-NULL
						componentList<-NULL
					}
					
					CompoundNum <- ifelse(is.null(modelPkList),0,nrow(modelPkList))
					
					for (i in uniqueClust) {
						curCompPkInd <- round(median(subset(localEICPeakList,clust==i)$pkInd))
						curLeftbound <- round(median(subset(localEICPeakList,clust==i)$Lbound))
						curRightbound <- round(median(subset(localEICPeakList,clust==i)$Rbound))
						
						if (curLeftbound < curCompPkInd & curCompPkInd < curRightbound) {
							
							# initialization
							curSpecList<-data.frame(mz=as.integer(),int=as.integer())
							curComponentList<-data.frame(mz=as.integer(),pkInd=as.integer(),Intensity=as.integer(),lboundInd=as.integer(),rboundInd=as.integer(),
									area=as.integer(),windowID=as.integer(),compoundID=as.integer(),isModel=as.integer(),isUnique=as.integer(),factor=as.numeric(),compID=as.integer())
							
							# extract mass spectrum
							for (m in 1:length(allprofiles)) {
								curProfile <- allprofiles[[m]]
								
								curArea <- sum(subset(curProfile,ET>=curLeftbound&ET<=curRightbound)$int)	
								curPk <- subset(curProfile,ET==curCompPkInd)
								if (nrow(curPk)==1) {
									curIntensity <- as.integer(curPk$int)
									curMz <- as.integer(localEICPeakList$mz[m])
									
									curSpecList<-rbind(curSpecList,data.frame(mz=as.integer(curMz),int=curIntensity))
									curComponentList<-rbind(curComponentList,data.frame(mz=as.integer(curMz),pkInd=curCompPkInd,lboundInd=curLeftbound,rboundInd=curRightbound,
													Intensity=curIntensity,area=curArea,windowID=windowID,compoundID=CompoundNum+1,isModel=0,isUnique=0,factor=0,compID=compNum+1))	
								}
								
							}
							
							# check quality and import spectrum and component into list
							if (nrow(curSpecList) >5) {
								mostIntenseMassPos <- which.max(curComponentList$Intensity)
								curComponentList[mostIntenseMassPos,]$isModel <- 1
								curComponentList[mostIntenseMassPos,]$isUnique <- 1
								curModelmass <- curComponentList[mostIntenseMassPos,]$mz
								specList[[CompoundNum+1]] <- curSpecList
								names(specList)[CompoundNum+1]<-paste(curCompPkInd,windowID,CompoundNum+1,curModelmass,compNum+1)
								componentList[[CompoundNum+1]]<-curComponentList
								CompoundNum <- CompoundNum+1
								compNum <- compNum +1
								print(paste(compNum,"X",sep="  "))
							}
						} # check qualify
					}	# for loop
				} # end check missing compounds
				
				componentResults<-c(componentResults,componentList)
				specResults<-c(specResults,specList)	
				
			}	# Check localEICPeakList 2
		}# Check localEICPeakList 1
	}#for loop 	
	
	# step 5: output
	cat(fileName,"finished deconvolution...\n")
	lib<-readMSP2Spec(filename=libFile,withRT=F)
	libmatch <- libMatching_Top10(lib,specResults,libMatch_Th) # target lib matching 
	libmatch$ET <- (libmatch$pkind-1)*params$ScanInterval+params$delaytime
	
	componentResults<-do.call(rbind,componentResults)	
	componentResults<-subset(componentResults,Intensity>0,select=c("mz","pkInd","lboundInd","rboundInd","Intensity","area","windowID","compoundID","isModel","isUnique","factor","compID"))
	write.csv(componentResults,row.names=FALSE,file=paste(params$WorkDir,"/output/decon/",fileName,"Decon.csv",sep=""))
	write.csv(subset(libmatch,select=c("ET","modelMass","compoundName","score")),file=paste(params$WorkDir,"/output/decon/",fileName,"_libmatch.csv",sep=""))
	#write(strOutput,file=outputPath)
	
	cat(fileName,"All deconvolution results are exported successfully!\n")
}



