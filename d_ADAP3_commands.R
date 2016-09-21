# Written by Xiuxia Du in September 2015





rm(list=ls())
graphics.off()





# set paths
path_to_R_libs <- paste("Users", "xdu4", "Documents", "Duxiuxia", "R_libs", sep=.Platform$file.sep)
path_to_R_libs <- paste(.Platform$file.sep, path_to_R_libs, .Platform$file.sep, sep="")

codeDir <- paste("Users", "xdu4", "Documents", "Duxiuxia", "Bitbucket", "adap-gc_3", sep=.Platform$file.sep)
codeDir <- paste(.Platform$file.sep, codeDir, .Platform$file.sep, sep="")

workDir <- paste("Users", "xdu4", "Documents", "Duxiuxia", "Temp", "ADAP3_test", sep=.Platform$file.sep)
workDir <- paste(.Platform$file.sep, workDir, .Platform$file.sep, sep="")
# 
# path_to_R_libs <- "/Users/xdu4/Documents/Duxiuxia/R_libs/"
# codeDir <- "/Users/xdu4/Documents/Duxiuxia/Bitbucket/adap-gc_3/"
# workDir <- "/Users/xdu4/Documents/Duxiuxia/Temp/ADAP3_test/"


setwd(workDir)





# install packages
packageName <- "ncdf4"
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





if (!dir.exists("output")) {
    dir.create(path="output")
}




# set parameters
params <- list()
params$workDir <- workDir

params$peak_span <- 11
params$valley_span <- 9

params$BHR_EIC <- 0.3
params$edgeHightDiffRatio_EIC <- 0.2
params$maxWindowlength_EIC <- 350

params$BHR_RIC <- 0.4
params$edgeHightDiffRatio_TIC <- 0.2
params$maxWindowlength_TIC <- 450

params$MaxPeakWidth <- 200
params$StN_Th <- 20





source(paste(codeDir, "d_peakpicking.r", sep=.Platform$file.sep))
EICpeakpicking(params)
