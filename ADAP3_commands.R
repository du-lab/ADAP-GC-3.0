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



#########################################
# Xiuxia code 
#########################################

rm(list=ls())
graphics.off()




path_to_R_libs <- "/Users/xdu4/Documents/Duxiuxia/R_libs"





packageName <- "ncdf"
if (!is.element(el=packageName, set=installed.packages(lib.loc=path_to_R_libs)[,1])) {
    install.packages(pkgs=packageName, 
                     lib=path_to_R_libs, 
                     repos="http://cran.us.r-project.org",
                     dependencies=T)
}
library(package=packageName,
        lib.loc=path_to_R_libs,
        character.only=T)




packageName <- "splus2R"
if (!is.element(el=packageName, set=installed.packages(lib.loc=path_to_R_libs)[,1])) {
    install.packages(pkgs=packageName, 
                     lib=path_to_R_libs, 
                     repos="http://cran.us.r-project.org",
                     dependencies=T)
}
library(package=packageName,
        lib.loc=path_to_R_libs,
        character.only=T)





packageName <- "ifultools"
if (!is.element(el=packageName, set=installed.packages(lib.loc=path_to_R_libs)[,1])) {
    install.packages(pkgs=packageName, 
                     lib=path_to_R_libs, 
                     repos="http://cran.us.r-project.org",
                     dependencies=T)
}
library(package=packageName,
        lib.loc=path_to_R_libs,
        character.only=T)




packageName <- "wmtsa"
if (!is.element(el=packageName, set=installed.packages(lib.loc=path_to_R_libs)[,1])) {
    install.packages(pkgs=packageName, 
                     lib=path_to_R_libs, 
                     repos="http://cran.us.r-project.org",
                     dependencies=T)
}
library(package=packageName,
        lib.loc=path_to_R_libs,
        character.only=T)





codeDir <- "/Users/xdu4/Documents/Duxiuxia/Bitbucket/adap-gc_3"
workDir <- "/Users/xdu4/Documents/Duxiuxia/Temp/ADAP3_test"





source(paste(codeDir, "d_peakpicking.r", sep=.Platform$file.sep))




params <- list()
params$workDir <- workDir
params$peak_span <- 21
params$valley_span <- 11





EICpeakpicking(params)




path_to_R_libs <- "/Users/xdu4/Documents/Duxiuxia/R_Libs"





package_name <- "readMzXmlData"
if (!is.element(package_name, installed.packages()[,1])) {
    install.packages(package_nmae, 
                     lib=path_to_R_libs, 
                     repos="http://cran.us.r-project.org",
                     dependencies=T)
}

library(package_name, lib.loc=path_to_R_libs, character.only=T)





in_file_name <- "Pos_26.mzXML"
in_file_path <- "/Users/xdu4/Documents/Duxiuxia/Temp"
in_file_full_name <- paste(in_file_path, in_file_name, sep=.Platform$file.sep)

xs <- readMzXmlFile(in_file_full_name)

scan <- 4212

mz_vec <- xs[[scan]]$spectrum$mass
int_vec <- xs[[scan]]$spectrum$intensity

dataOut <- cbind(mz_vec, int_vecz)
write.table(x=dataOut, file="aEIC.txt", row.names=F, col.names=c("mass", "intensity"), append=F, sep="\t")

intensity_threshold <- 500
II <- which(y>intensity_threshold)
mz_vec <- mz_vec[II]
int_vec <- int_vec[II]

x_1 <- mz_vec[-length(mz_vec)]
x_2 <- mz_vec[-1]
mz_diff <- x_2 - x_1

#     y_interpolated <- spline(x, y, method="fmm")
int_interpolated <- approx(mz_vec, int_vec, method="linear", n=length(mz_vec))
layout(matrix(data=c(1,2), nrow=2, byrow=T))
plot(mz_vec, int_vec, type="l")
plot(int_interpolated$x, int_interpolated$y, type="l")






# ------------------fake signal-------------------
sampling_interval=0.01
f1 <- 1
f2 <- 4
t <- seq(from=0, to=24*pi, by=sampling_interval)
x <- rep(0, times=length(t))
II <- which(t >= pi & t <= 10*pi)
x[II] <- 2*sin(f1*t[II])
plot(t, x, type="l")

time_shift <-12*pi
y <- rep(0, times=length(t))
II <- which(t >= time_shift & t < time_shift+9*pi)
y[II] <- sin(f2*t[II])
plot(t, y, type="l", col="red")

z <- x+y
plot(t, z, pch=16, cex=0.2)
# ------------------------------------------------





params <- list()
params$path_to_R_libs <- path_to_R_libs
params$snr <- 10

source("d_peakpicking.R")
s <- getPeaks(x=t, y=z, params)


