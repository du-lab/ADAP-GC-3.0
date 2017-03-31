codeDir <<-"code"
source(paste(codeDir,"pipeline.r",sep="/"))
#params<-readParaFromCmd("/home/yni/DataTest/U50Stds2013/U50Stds2013.db")
params<-readParaFromCmd("/home/yni/DataTest/LI102912/LI102912.db")
params$codeDir<-codeDir


library("snow")
library("ncdf")

WorkDir <- params$WorkDir
DataFilelist <- params$DataFiles

nNode<-as.integer(params$nNode)
clustType<-params$clustType
cl<-makeCluster(nNode,type=clustType)

args <- commandArgs(TRUE)
start <- as.integer(args[1])
end <- as.integer(args[2])

for(fileindex in start:end)
#for(fileindex in 1:length(DataFilelist))
{
	
	###############
	#read EIC data
	################
	inFilePath <- DataFilelist[[fileindex]]
	fileName<-parseFileName(inFilePath)
	
	cat(fileName,"reading EIC data...\n")
	#EICFile<-paste(WorkDir,"/output/EIC/",fileName,"EIC.cdf",sep="")
	
	#check if it exist denoised EIC file
	rawEICFile <-paste(WorkDir,"output/EIC/",fileName,"EIC.cdf",sep="")
	denoiseEICFile<-paste(WorkDir,"output/EIC/denoised_",fileName,"EIC.cdf",sep="")
	EICFile<-ifelse(file.exists(denoiseEICFile),denoiseEICFile,rawEICFile)
	
	if(file.exists(inFilePath)==FALSE)next
	ncid <- open.ncdf(EICFile)
	vecInt <- get.var.ncdf(ncid, varid="intVec")
	mzVec<-get.var.ncdf(ncid, varid="mzVec")
	
	#mzVec <- c(40:400) # for testing
	
	close.ncdf(ncid)
	remove(ncid)	
	
	###############
	#peak picking
	################
	cat(fileName,"peak picking...")
	#######peak picking
	mzGroups<-clusterSplit(cl,mzVec)
	
	clusterExport(cl,"peaks")#assign the current profiles to all nodes
	clusterExport(cl,"getPeaks")#assign the current profiles to all nodes
	clusterExport(cl,"GaussianDisimilarity")
	clusterExport(cl,"Sharpness")#assign the current profiles to all nodes
	clusterExport(cl,"Sharpness1")#assign the current profiles to all nodes
	clusterExport(cl,"dotProduct")#assign the current profiles to all nodes
	clusterExport(cl,"parseFileName")#assign the current profiles to all nodes
	clusterExport(cl,"MSWtoCWT")
	clusterExport(cl,"cwt")
	clusterExport(cl,"extendLength")
	clusterExport(cl,"extendNBase")
	clusterExport(cl,"wavCWTTree_ModifyByYan")
	clusterExport(cl,"wavCWTTree_ModifyByYan_Min")
	clusterExport(cl,"wavCWTPeaks_Modify")
	
	time1<-Sys.time()
	PeakList<-clusterApply(cl,mzGroups,getPeaksGroup,vecInt,mzVec,params,fileName)#parallel version
	time2<-Sys.time()
	time2-time1
	PeakList<-unlist(PeakList,recursive=F)
	PeakList<-do.call(rbind,PeakList)
	
	write.csv(PeakList,file=paste(WorkDir,"/output/peakpicking/",fileName,"_EIC_PeakList.csv",sep=""),row.names=F)
	
}
stopCluster(cl)
