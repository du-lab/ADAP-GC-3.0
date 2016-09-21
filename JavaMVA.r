#setwd("/Users/yanni/Documents/workspace/")
codeDir<-"code"
#codeDir<-"MSAP-Java"
source(paste(codeDir,"DataPretreatment_MSAP2.r",sep="/"))
source(paste(codeDir,"MVA_MSAP2.r",sep="/"))
source(paste(codeDir,"readPara.r",sep="/"))
#source(paste(codeDir,"MVAFunctions_MSAP2.r",sep="/"))
params<-readParaFromCmd(commandArgs(TRUE))
#params<-readParaFromCmd("~yanni/DataTest/stds/stds.db")
params$codeDir<- codeDir

if (params$isPretreatment==1) {
	pretreatment_data <- pretreatment_process(params)	
	#output preprocessed data
	processed_data <- do.call(cbind,pretreatment_data)
	#processed_data <- processed_data[,-which(colnames(processed_data)=="Y.Sample")]
	write.csv(processed_data,paste(params$WorkDir,"output/preprocessed_data.csv",sep=""))
	
	if (params$isPCA==1) {
		pca(params,pretreatment_data)
	}
	
	if (params$isPLSDA==1) {
		pls(params,pretreatment_data)
	}
}




