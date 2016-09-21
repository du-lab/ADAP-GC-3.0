codeDir<-"code"
source(paste(codeDir,"pipeline.r",sep="/"))
params<-readParaFromCmd(commandArgs(TRUE))
params$codeDir<-codeDir
params$NonUmassVec<-c(1:50,147,221)
#params$NonUmassVec<-c(1:50,73,147,221)
bioMarkerTableSQLite(params)
