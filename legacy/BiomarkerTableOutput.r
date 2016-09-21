######################################################
#Input: aligned file
#Output: 4 biomarker tables & 1 biomarker list 
#Choice: PK intensity or area could be extracted
#Yan Ni, Aug 10, 2011
#Modified: Mar 6, 2012
######################################################

bioMarkerTableSQLite<-function(params,abundance="intensity")
{
	library(RSQLite)
	library(snow)
	JobName<-params$JobName
	workingDir <- params$WorkDir
	paraDB<- params$dbFile
	
	#========================================
	# read aligned file from database
	#========================================
	m <- dbDriver("SQLite")
	con <- dbConnect(m, dbname = paraDB)
	tname<-paste("t",JobName,"aligned",sep="_")
	##to write the produced csv file into database - modified by yan
	if(!dbExistsTable(con,tname))
	{
		aligndata_path <-paste(params$WorkDir,"output/",paste(paste(JobName,"aligned",sep="_"),"csv",sep="."),sep="")
		aligndata <- read.csv(aligndata_path)
		dbWriteTable(con,tname,aligndata)
	}
	##continue previous script
	if(!dbExistsTable(con,tname))
		stop(paste("Table '",tname,"' does not exist!",sep=""))
	featureSets<-dbReadTable(con,tname)
	featureSets<-subset(featureSets,int>0)
	
	#========================================
	# create biomarker list
	#========================================
	FileIDs<-unique(featureSets$FileID)
	nFile<-length(FileIDs)
	fileVsClust<-subset(featureSets,select=c("FileID","clustID"))#get common mother Ions which occur in 80% of datasets
	fileVsClust<-unique(fileVsClust)
	groupCount<-aggregate(fileVsClust$FileID,by=list(fileVsClust$clustID),length)#group count fileid by clustid
	nMajority<-round(as.numeric(params$SampleRatio)*nFile)
	CommonClustID<-subset(groupCount,x>=nMajority)$Group.1# get common clustid which occurs in over 4 datasets
	featureSets<-subset(featureSets,clustID%in%CommonClustID)#filter featureset by common clustid
	nCommonClustID <- length(CommonClustID)	
	#To make sure that common cluster ID has enough number or not
	params$nNode <- min(nCommonClustID,as.integer(params$nNode))
	c1<-makeCluster(as.numeric(params$nNode),type=params$clustType)	
	commonClustGroups<-clusterSplit(c1,CommonClustID)
	clusterExport(c1,"transFeatureSet")#assign the current function to all nodes
	clusterExport(c1,"weightSpec")#assign the current function to all nodes
	clusterExport(c1,"parseFileName")
	bioMarkerList<-clusterApplyLB(c1,commonClustGroups,transFeatureSet,FileIDs,featureSets,params,abundance)
	bioMarkerList<-do.call(rbind,bioMarkerList)
	bioMarkerList<-bioMarkerList[order(bioMarkerList$alignedET),]
	stopCluster(c1)
	##write into database and output csv file
	##Java code take care of inserting BioMarkerList table
	##database table's column name restriction, let R add X to 1_1, 1_2 
	dbWriteTable(con,"bioMarkerList",bioMarkerList,row.names=F,overwrite=T)	
	dbDisconnect(con)
	write.csv(bioMarkerList,file=paste(workingDir,"/output/",JobName,"_bioMarkerList.csv",sep=""),row.names=FALSE)
	
	#========================================
	# output 4 quali/quan tables
	#========================================
	QualQuanTableOutput(params,bioMarkerList) 
	cat("QualQuan tables have been output!\n")
}

transFeatureSet<-function(commonClustGroup,FileIDs,featureSets,params,abundance)
{
	nClust<-length(commonClustGroup)
	localTrainingSet<-NULL
	
	##To create raw file names list
	fileNames <- NULL
	fileNamesList <- params$DataFiles
	for (i in 1:length(fileNamesList)) {
		fileNames<-c(fileNames,as.character(parseFileName(fileNamesList[[i]])))
	}
	
	#========================================
	# Extract peak heights
	#========================================
	if (abundance=="intensity") {
		for(ind in 1:nClust)
		{
			print(ind)
			##for each clustID,transpose the fragments of different fileID by union merge 
			curClustID<-commonClustGroup[ind]
			curCompList<-subset(featureSets,clustID==curClustID)
			
			#========================================
			# get the most occurant model peak
			#========================================
			modelPk_list <- subset(curCompList,isModel==1)
			ModelPkMass <- as.integer(names(which.max(table(modelPk_list$mz))))		
			##avoid certain model peak with low intensity resulting worse quantitation performance
			if (max(table(subset(curCompList,isModel==1)$mz)) < round(nrow(modelPk_list)*0.4) & ModelPkMass>=320) {		
				modelPk_meanIntensity <- NULL
				for (mMz in unique(modelPk_list$mz)) {
					modelPk_meanIntensity <- c(modelPk_meanIntensity,mean(subset(curCompList,mz==mMz)$int))
				}
				ModelPkMass <- unique(modelPk_list$mz)[which.max(modelPk_meanIntensity)]
			}	
			
			#========================================
			# get the most occurant unique peak
			#========================================
			uniquePk_list <- subset(curCompList,isUnique==1)
			UniquePkMass <- as.integer(names(which.max(table(uniquePk_list$mz))))
			## select most intense unique mass if no one is most occurant
			if (max(table(uniquePk_list$mz))==1&nrow(uniquePk_list)>1) {
				UniquePkMass <- uniquePk_list[which.max(uniquePk_list$int),]$mz
			}
			
			#========================================
			# To extract PK intensity information
			#========================================
			##initialize the first file column
			curCmp<-subset(curCompList,FileID==FileIDs[1],select=c("mz","int"))
			curCmp <- curCmp[!duplicated(curCmp[,c("mz","int")]),]
			
			##merge the same fragments from the same file
			if(nrow(curCmp)>0)
				curCmp<-aggregate(curCmp$int,by=list(curCmp$mz),sum)
			names(curCmp)<-c("mz",fileNames[FileIDs[1]])
			
			##for loop to merge the rest of samples
			for(curFileID in FileIDs[-1])
			{
				curClust<-subset(curCompList,FileID==curFileID,select=c("mz","int"))
				curClust <- curClust[!duplicated(curClust[,c("mz","int")]),]
				if(nrow(curClust)>0)
					curClust<-aggregate(curClust$int,by=list(curClust$mz),sum)
				names(curClust)<-c("mz",fileNames[curFileID])
				curCmp<-merge(curCmp,curClust,by.x="mz",by.y="mz",all=T)
			}
			
			##define NAs as 0
			curCmp[is.na(curCmp)]<-0
			
			#========================================
			# organize peak list 
			#========================================
			#Extract ET,clustID, modelpk, uniqPk information
			curCmp$clustID<-curClustID
			curCmp$alignedET<-curCompList$alignedET[1]
			curCmp$isModel<-0
			curCmp[which(curCmp$mz==ModelPkMass),]$isModel=1
			curCmp$isUnique<-0
			curCmp[which(curCmp$mz==UniquePkMass),]$isUnique=1
			
			#========================================
			# get most intense peak
			#========================================
			if(nrow(curCmp)>0)
			{
				curCmp$int<-0
				curCmp$isMaxUmass<-0
				
				##get the median int of each ion intensity accross non-zero-value samples
				dataMatrix <- curCmp[,colnames(curCmp)%in%fileNames]			
				for(i in 1:nrow(curCmp))
				{
					nonZeroInd<-dataMatrix[i,]>0
					nonZeroInt<-dataMatrix[i,][nonZeroInd]
					curCmp[i,]$int<-median(nonZeroInt)
				}
				
				##add weight to median intensity
				curCmp<-weightSpec(curCmp)		
				curCmp_extract <- curCmp[!curCmp$mz%in%params$NonUmassVec,]
				
				##trying to get maximal mz which is not 73
				curCmp_extract <- curCmp_extract[order(curCmp_extract[,"int"],decreasing=T),]
				mzlist <- curCmp_extract$mz
				if (mzlist[1]==73&mzlist[2]<=250)
				{
					curCmp[which(curCmp$mz==mzlist[2]),]$isMaxUmass<-1
				}else
				{
					curCmp[which(curCmp$mz==mzlist[1]),]$isMaxUmass<-1
				}				
			}
			
			InfoMatrix <- subset(curCmp,select=c("clustID","alignedET","mz","int","isModel","isUnique","isMaxUmass"))
			curCmp <- cbind(InfoMatrix,dataMatrix)	##reformat the output and merge to biomarker list
			localTrainingSet<-rbind(localTrainingSet,curCmp)
		}#for loop
	}#if statement
	
	#========================================
	# Extract peak area
	#========================================
	if (abundance=="area") {	
		for(ind in 1:nClust)
		{
			##for each clustID,transpose the fragments of different fileID by union merge 
			curClustID<-commonClustGroup[ind]
			curCompList<-subset(featureSets,clustID==curClustID)
			
			##get the most occurant model peak
			modelPk_list <- subset(curCompList,isModel==1)
			ModelPkMass <- as.integer(names(which.max(table(modelPk_list$mz))))		
			##avoid certain model peak with low intensity resulting worse quantitation performance
			if (max(table(subset(curCompList,isModel==1)$mz)) < round(nrow(modelPk_list)*0.5) & ModelPkMass>=300) {
				ModelPkMass <- modelPk_list[which.max(modelPk_list$pkArea),]$mz
			}	
			
			#========================================
			# get the most occurant model peak
			#========================================
			uniquePk_list <- subset(curCompList,isUnique==1)
			UniquePkMass <- as.integer(names(which.max(table(uniquePk_list$mz))))
			## select most intense unique mass if no one is most occurant
			if (max(table(uniquePk_list$mz))==1&nrow(uniquePk_list)>1) {
				UniquePkMass <- uniquePk_list[which.max(uniquePk_list$pkArea),]$mz
			}
			
			#========================================
			# To extract PK intensity information
			#========================================
			
			##initialize the first file column
			curCmp <-subset(curCompList,FileID==FileIDs[1],select=c("mz","pkArea"))
			curCmp <- curCmp[!duplicated(curCmp[,c("mz","pkArea")]),]
			
			##merge the same fragments from the same file
			if(nrow(curCmp_area)>0)
				curCmp<-aggregate(curCmp$pkArea,by=list(curCmp$mz),sum)
			names(curCmp)<-c("mz",fileNames[FileIDs[1]])
			
			##for loop to merge the rest of samples
			for(curFileID in FileIDs[-1])
			{
				curClust<-subset(curCompList,FileID==curFileID,select=c("mz","pkArea"))
				curClust <- curClust[!duplicated(curClust[,c("mz","pkArea")]),]
				if(nrow(curClust)>0)
					curClust<-aggregate(curClust$pkArea,by=list(curClust$mz),sum)
				names(curClust)<-c("mz",fileNames[curFileID])
				curCmp<-merge(curCmp,curClust,by.x="mz",by.y="mz",all=T)
			}
			curCmp[is.na(curCmp)]<-0
			
			#========================================
			# organize peak list
			#========================================
			#Extract ET,clustID, modelpk, uniqPk information
			curCmp$clustID<-curClustID
			curCmp$alignedET<-curCompList$alignedET[1]
			curCmp$isModel<-0
			curCmp[which(curCmp$mz==ModelPkMass),]$isModel=1
			curCmp$isUnique<-0
			curCmp[which(curCmp$mz==UniquePkMass),]$isUnique=1
			
			#========================================
			# get most intense peak
			#========================================
			if(nrow(curCmp)>0)
			{
				curCmp$int<-0
				curCmp$isMaxUmass<-0
				
				##get the median area of each ion intensity accross non-zero-value samples
				dataMatrix <- curCmp[,colnames(curCmp)%in%fileNames]			
				for(i in 1:nrow(curCmp))
				{
					nonZeroInd<-dataMatrix[i,]>0
					nonZeroInt<-dataMatrix[i,][nonZeroInd]
					curCmp[i,]$int<-median(nonZeroInt)
				}
				
				##add weight to median intensity
				curCmp<-weightSpec(curCmp)		
				curCmp_extract <- curCmp[!curCmp$mz%in%params$NonUmassVec,]
				
				##trying to get maximal mz which is not 73
				curCmp_extract <- curCmp_extract[order(curCmp_extract[,"int"],decreasing=T),]
				mzlist <- curCmp_extract$mz
				if (mzlist[1]==73&mzlist[2]<=250)
				{
					curCmp[which(curCmp$mz==mzlist[2]),]$isMaxUmass<-1
				}else
				{
					curCmp[which(curCmp$mz==mzlist[1]),]$isMaxUmass<-1
				}				
			}
			
			InfoMatrix <- subset(curCmp,select=c("clustID","alignedET","mz","int","isModel","isUnique","isMaxUmass"))
			##reformat the output and merge to biomarker list
			curCmp <- cbind(InfoMatrix,dataMatrix)	
			localTrainingSet<-rbind(localTrainingSet,curCmp)
		}#for loop
	}#if statement
	localTrainingSet
}

QualQuanTableOutput <- function(params,bioMarkerList) 
{
	choice<-c("S","Q","U", "M")
	JobName<- params$JobName
	workDir <- params$WorkDir
	trainingSet<-bioMarkerList
	Ncol_Raw <- ncol(trainingSet)
	trainingSets <- NULL
	
	#========================================
	# Select 4 different Quant Mass
	#========================================
	
	# 1 #
	#to summarize all ions of each component, keep two columns:clustID and alignedET 
	trainingSet_Sum<-trainingSet[,c(1,2,8:Ncol_Raw)]
	Ncol_Sum <- ncol(trainingSet_Sum)
	trainingSets[[1]] <- aggregate(trainingSet_Sum[,3:Ncol_Sum],by=list(clustID=trainingSet_Sum$clustID,alignedET=trainingSet_Sum$alignedET),FUN=sum)
	
	# 2 #
	# Most intense mass, keep three columns-clustID,alignedET and MZ	
	trainingSet_Max<-subset(trainingSet,isMaxUmass==1)
	trainingSets[[2]]<-trainingSet_Max[,c(1:3,8:Ncol_Raw)]
	
	# 3 #
	# unique ions
	trainingSet_Unique<-subset(trainingSet,isUnique==1)
	trainingSets[[3]]<-trainingSet_Unique[,c(1:3,8:Ncol_Raw)]
	
	# 4 #
	# model peak
	trainingSet_Model<-subset(trainingSet,isModel==1)
	trainingSets[[4]]<-trainingSet_Model[,c(1:3,8:Ncol_Raw)]
	
	#========================================
	# Write 4 biomarker tables out as .csv
	#========================================
	
	# to select Qmass occuring more than nMarjority of samples
	nFile<-length(params$DataFiles)
	nMajority<-round(as.numeric(params$SampleRatio)*nFile)
	
	for (j in c(1:4)) {
		Cur_trainingSet <- trainingSets[[j]]
		Ncol <- ncol(Cur_trainingSet)
		CountNum <- NULL
		
		for (i in 1:nrow(Cur_trainingSet)) {
			count <- sum(Cur_trainingSet[i,][,c(4:Ncol)]>0)
			CountNum <- c(CountNum,count)
		}
		# add count number to filter and then remove it
		Cur_trainingSet$"count"<- CountNum
		Cur_trainingSet <- subset(Cur_trainingSet,count>=nMajority)
		Cur_trainingSet <- Cur_trainingSet[,-which(colnames(Cur_trainingSet)=="count")]
		write.csv(Cur_trainingSet,paste(workDir,"output/",JobName,"_",choice[j],"_","BiomarkerTable.csv",sep=""),row.names=F)
	}
	
	# Create a factor file for user to fill in sample information
	Sample <- colnames(trainingSet)[8:ncol(trainingSet)]
	factor <- as.data.frame(cbind(Sample,Y=rep(0,length(Sample))))
	write.csv(factor,paste(workDir,"output/factor.csv",sep=""),row.names=F)
}	

## to select either Qmass or Summation of all fragments
## to represent a component/potential compound

data_select <- function(params) 
{
	choice<-params$QuantitationType
	paraDB<- params$dbFile
	JobName<- params$JobName
	workDir <- params$WorkDir
	m <- dbDriver("SQLite")
	con <- dbConnect(m, dbname = paraDB)
	tname<-"bioMarkerList"

	if(!dbExistsTable(con,tname))
		stop(paste("Table '",tname,"' does not exist!",sep=""))
	trainingSet<-dbReadTable(con,tname)
	Ncol <- ncol(trainingSet)
	
	if (choice=="S") {
		
		##to summarize all ions of each component
		##Keep two columns:clustID and alignedET 
		trainingSet<-trainingSet[,c(1,2,8:Ncol)]
		Ncol <- ncol(trainingSet)
		trainingSet <- aggregate(trainingSet[,3:Ncol],by=list(clustID=trainingSet$clustID,alignedET=trainingSet$alignedET),FUN=sum)
		
		## to select Qmass occuring more than nMarjority of samples
		nFile<-length(params$DataFiles)
		nMajority<-round(as.numeric(params$sampleRatio)*nFile)
		
		
		CountNum <- NULL		
		for (i in 1:nrow(trainingSet)) {
			count <- sum(trainingSet[i,][,c(3:Ncol)]>0)
			CountNum <- c(CountNum,count)
		}
		
		##add count number to filter and then remove it
		trainingSet$"count"<- CountNum
		trainingSet <- subset(trainingSet,count>=nMajority)
		trainingSet <- trainingSet[,-which(colnames(trainingSet)=="count")]
	}else {
		
		##Other three conditions: keep three columns-clustID,alignedET and MZ	
		if (choice=="Q") {
			trainingSet<-subset(trainingSet,isMaxUmass==1)
			trainingSet<-trainingSet[,c(1:3,8:Ncol)]
		}
		
		##To select unique ions
		if (choice=="U") {
			trainingSet<-subset(trainingSet,isUnique==1)
			trainingSet<-trainingSet[,c(1:3,8:Ncol)]
		}
		
		##to select model peak for each compoenent
		if (choice=="M") {
			trainingSet<-subset(trainingSet,isModel==1)
			trainingSet<-trainingSet[,c(1:3,8:Ncol)]
		}
		
		## to select Qmass occuring more than nMarjority of samples
		nFile<-length(params$DataFiles)
		nMajority<-round(as.numeric(params$sampleRatio)*nFile)
		
		Ncol <- ncol(trainingSet)
		CountNum <- NULL
		
		for (i in 1:nrow(trainingSet)) {
			count <- sum(trainingSet[i,][,c(4:Ncol)]>0)
			CountNum <- c(CountNum,count)
		}
		# add count number to filter and then remove it
		trainingSet$"count"<- CountNum
		trainingSet <- subset(trainingSet,count>=nMajority)
		trainingSet <- trainingSet[,-which(colnames(trainingSet)=="count")]
	}
	print("selection is done!\n")
	# clean up
	dbDisconnect(con)
	# save the csv file under output dir.
	write.csv(trainingSet,paste(workDir,"output/",JobName,"_",choice,"_","BiomarkerTable.csv",sep=""),row.names=F)
}	
