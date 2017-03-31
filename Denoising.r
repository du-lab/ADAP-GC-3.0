#parTICDenoising<-function(params) 
#{
#	time1<-Sys.time()
#	library("snow")
#	###############
#	#get parameters
#	################
#	WorkDir <- params$WorkDir
#	DataFilelist <- params$DataFiles
#	codeDir<-params$codeDir
#	nNode<-as.integer(params$nNode)
#	clustType<-params$clustType
#	
#	c1<-makeCluster(nNode,type=clustType)
#	
#	clusterExport(c1,"TICDenoising")#assign the current function to all nodes
#	clusterApply(c1,DataFilelist,TICDenoising,WorkDir,codeDir,isBslnCrt=T)#parallel version
#	
#	stopCluster(c1)
#	time2<-Sys.time()
#	time2-time1
#}
#
#TICDenoising<-function(inFilePath,WorkDir,codeDir,isBslnCrt=F) 
#{
#	
#	library("ncdf")
#	source(paste(codeDir,"pipeline.r",sep="/"))
#	###############
#	#read EIC data
#	###############
#	#inFilePath <- DataFilelist[[fileindex]]
#	fileName<-parseFileName(inFilePath)
#	cat(fileName,"reading TIC data...\n")
#	
#	if(file.exists(inFilePath)==FALSE)next
#	ncid <- open.ncdf(inFilePath,write=TRUE)
#	vecInt <- get.var.ncdf(ncid, varid="total_intensity")
#	vecET<-get.var.ncdf(ncid, varid="scan_acquisition_time")
#	
#	totalscan<-length(vecInt)
#	
#	###########
#	#smoothing
#	###########
#	cat(fileName,"smoothing TIC data...\n")
#	
#	smoothingWindow<-10#50
#	ma = rep(1, smoothingWindow)/smoothingWindow
#	sn.ma=filter(vecInt,ma)
#	sn.ma[which(is.na(sn.ma))]<-0
#	vecInt<-sn.ma
#	
#	#####################
#	#baseline correction
#	#####################
#	if(isBslnCrt==T)
#	{
#		nbinSize<-50#240
#		isMin <- peaks(-vecInt, span=nbinSize)
#		minInd<-which(isMin==TRUE)	
#		
#		intervalOfScans<-c(minInd[1],minInd[-1]-minInd[-length(minInd)])
#		
#		bgs<-c(rep(vecInt[minInd],intervalOfScans),rep(0,totalscan-minInd[length(minInd)]))
#		f.lo <- loess(bgs ~ c(1:totalscan), span =0.05, degree = 2)
#		bsln <- f.lo$fitted
#		
#		## baseline subtraction
#		vecInt<-vecInt-bsln
#		vecInt[vecInt<0]<-0
#	}
#	
#	###################
#	#write to new  CDF
#	###################
#	cat(fileName,"write denoised TIC data...\n")
#	
#	ncid1<-create.ncdf(filename=paste(WorkDir,"/output/TIC/denoised_",fileName,"_TIC.cdf",sep=""),vars=ncid$var)
#	put.var.ncdf(ncid1,"total_intensity",vecInt)
#	close.ncdf(ncid1)
#	remove(ncid1)	
#	
#	close.ncdf(ncid)
#	remove(ncid)	
#}

parDenoising<-function(params) 
{
	library("snow")
	
	nNode<-params$nNode
	clustType<-params$clustType
	codeDir<-params$codeDir
	WorkDir<-params$WorkDir
	DataFilelist<-params$DataFiles
	
	cat("start denoising EIC data...\n")
	time1<-Sys.time()
	c1<-makeCluster(nNode,type=clustType)
	clusterExport(c1,"Denoising") ## assign the current function to all nodes
	clusterApply(c1,DataFilelist,Denoising,WorkDir,codeDir,isBslnCrt=T,isSm=T) ## parallel version
	stopCluster(c1)
	time2<-Sys.time()
	cat("end denoising EIC data...\n")
	time2-time1
}

Denoising<-function(inFilePath,WorkDir,codeDir,isBslnCrt=F,isSm=F) 
{
	library("ncdf")
	source(paste(codeDir,"pipeline.r",sep="/"))
	
	#====================================
	#read EIC data
	#====================================
	fileName<-parseFileName(inFilePath)
	
	cat(fileName,"reading EIC data...\n")
	EICFile<-paste(WorkDir,"output/EIC/",fileName,"EIC.cdf",sep="")
	ncid <- open.ncdf(EICFile)
	vecInt <- get.var.ncdf(ncid, varid="intVec")
	mzVec<-get.var.ncdf(ncid, varid="mzVec")
	totalscan<-length(vecInt)/length(mzVec)
	
	#====================================
	#smoothing and baseline correction
	#====================================
	cat(fileName,"smoothing EIC data...\n")
	
	DenoisedTIC <- 0
	time1<-Sys.time()
	
	for(mzInd in 1:length(mzVec))
	{
		mz<-mzVec[mzInd]	
		startInd<-(mzInd-1)*totalscan+1
		endInd<-startInd+totalscan-1
		curVecInt <- vecInt[startInd:endInd]
		
		#smoothing
		if(isSm==T)
		{
			smoothingWindow<-20
			ma = rep(1, smoothingWindow)/smoothingWindow
			sn.ma=filter(curVecInt,ma)
			sn.ma[which(is.na(sn.ma))]<-0
			curVecInt<-sn.ma
		}	
		
		#baseline correction
		if(isBslnCrt==T)
		{
			nbinSize<-240
			isMin <- peaks(-curVecInt, span=nbinSize)
			minInd<-which(isMin==TRUE)	
			
			if (length(minInd)==0) {bsln=0} # avoiding no minInd found - yan
			else {
				intervalOfScans<-c(minInd[1],minInd[-1]-minInd[-length(minInd)])
				
				bgs<-c(rep(curVecInt[minInd],intervalOfScans),rep(0,totalscan-minInd[length(minInd)]))
				f.lo <- loess(bgs ~ c(1:totalscan), span =0.05, degree = 2)
				bsln <- f.lo$fitted
			}
			
			#baseline subtraction
			curVecInt<-curVecInt-bsln
			curVecInt[curVecInt<0]<-0
		}
		
		#update vector of intensity
		vecInt[startInd:endInd]<-curVecInt
		DenoisedTIC <- DenoisedTIC + curVecInt
	}
	
	time2<-Sys.time()
	time2-time1
	
	#====================================
	#update EIC and TIC CDF
	#====================================
	cat(fileName,"write denoised EIC data...\n")
	ncid1<-create.ncdf(filename=paste(WorkDir,"/output/EIC/denoised_",fileName,"EIC.cdf",sep=""),vars=ncid$var)
	put.var.ncdf(ncid1,"intVec",vecInt)
	put.var.ncdf(ncid1,"mzVec",mzVec)
	close.ncdf(ncid1)
	remove(ncid1)	
	
	cat(fileName,"write denoised TIC data...\n")
	dimTIC<-dim.def.ncdf(name="totalInt", units="", vals=1:length(DenoisedTIC), unlim=FALSE, create_dimvar=TRUE )
	varTIC<-var.def.ncdf("total_intensity","" ,dimTIC,0 )
	ncid1<-create.ncdf(filename=paste(WorkDir,"/output/TIC/denoised_",fileName,"_TIC.cdf",sep=""),vars=varTIC)
	put.var.ncdf(ncid1,"total_intensity",DenoisedTIC)
	close.ncdf(ncid1)
	remove(ncid1)	
	
	close.ncdf(ncid)
	remove(ncid)	
}


Background_Denoising<-function(inFilePath,WorkDir,codeDir,isBslnCrt=F,isSm=F) 
{
	library("ncdf")
	source(paste(codeDir,"pipeline.r",sep="/"))
	
	#====================================
	#read EIC data
	#====================================
	fileName<-parseFileName(inFilePath)
	
	cat(fileName,"reading EIC data...\n")
	EICFile<-paste(WorkDir,"output/EIC/",fileName,"EIC.cdf",sep="")
	ncid <- open.ncdf(EICFile)
	vecInt <- get.var.ncdf(ncid, varid="intVec")
	mzVec<-get.var.ncdf(ncid, varid="mzVec")
	totalscan<-length(vecInt)/length(mzVec)
	

	maxInt_Vec <- NULL
	av_Vec <- NULL
	av_Vec_Nonzero <- NULL
	med_Vec <- NULL
	med_Vec_Nonzero <- NULL
	
	for(mzInd in 1:length(mzVec))
	{
		mz<-mzVec[mzInd]	
		startInd<-(mzInd-1)*totalscan+1
		endInd<-startInd+totalscan-1
		curVecInt <- vecInt[startInd:endInd]
		
		maxInt_Vec <- c(maxInt_Vec,max(curVecInt))
		av_Vec <- c(av_Vec,max(curVecInt)/mean(curVecInt))
		med_Vec <- c(med_Vec,max(curVecInt)/median(curVecInt))
		
		curVecInt[which(curVecInt==0)]=NA
		
		av_Vec_Nonzero <- c(av_Vec_Nonzero,mean(curVecInt,na.rm=TRUE))
		med_Vec_Nonzero <- c(med_Vec_Nonzero,median(curVecInt,na.rm=TRUE))
		#if (max(curVecInt) <= 300) print (mz)
	}
	
}


















