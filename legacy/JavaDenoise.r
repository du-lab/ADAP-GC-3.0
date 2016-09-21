codeDir<-"code"
source(paste(codeDir,"pipeline.r",sep="/"))
params<-readParaFromCmd(commandArgs(TRUE))
params$codeDir<-codeDir	
parDenoising(params)
#parTICDenoising(params) 
