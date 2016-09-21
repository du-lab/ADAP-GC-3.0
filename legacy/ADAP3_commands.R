# For ADAP3

# * please check "svn up" to get the most updated version of source code *

# (1) First step is to use ADAP2.0.4 software, configure params and run peak picking and deconvolution if the software works

# Otherwise, once the reconstructed EIC files and the working directory are created, we can terminate peak picking and then do command testing 

# (2) command testing 
# (2.1) Load project parameters
codeDir<-"code"
source(paste(codeDir,"pipeline.r",sep="/"))
params<-readParaFromCmd("/home/yni/DataTest/U50Stds2013/U50Stds2013.db")
params$codeDir<-codeDir
fileInfo<-read.csv(paste(params$WorkDir,"output/",params$JobName,"_RT.csv",sep=""))
params$delaytime<-fileInfo$firstRT[1]
params$ScanInterval<-fileInfo$ScanInterval[1]

# (2.2) parallel peak picking & deconvolution
parTICPeakpicking(params,denoised=F)
parEICpeakpicking(params)
deconvolution(params,isDistParallel=FALSE,clustingType="h",isDataParalell=TRUE)

# (2.3) sequential peak picking & deconvolution
DataFilelist <- params$DataFiles
for(fileindex in 1:length(DataFilelist)) {
  inFilePath <- DataFilelist[[fileindex]]
  TICPeakpicking(inFilePath,params,denoised=FALSE)
}

EICpeakpicking(params)

deconvolution(params,isDistParallel=TRUE,clustingType="h",isDataParalell=FALSE)
