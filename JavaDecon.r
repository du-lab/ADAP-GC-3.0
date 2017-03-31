codeDir<-"code"
source(paste(codeDir,"pipeline.r",sep="/"))
params<-readParaFromCmd(commandArgs(TRUE))
params$codeDir<-codeDir
fileInfo<-read.csv(paste(params$WorkDir,"output/",params$JobName,"_RT.csv",sep=""))
params$delaytime<-fileInfo$firstRT[1]
params$ScanInterval<-fileInfo$ScanInterval[1]
deconvolution(params,isDistParallel=F,clustingType="h",isDataParalell=T)



