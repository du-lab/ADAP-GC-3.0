codeDir <<-"code"
source(paste(codeDir,"pipeline.r",sep="/"))
params<-readParaFromCmd(commandArgs(TRUE))
params$codeDir<-codeDir
parTICPeakpicking(params,denoised=F)
parEICpeakpicking(params)
#EICpeakpicking(params)


