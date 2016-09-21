# test_wmtsa.R

# Written by Xiuxia Du.

# Started on August 25, 2015.
# Modified on September 14, 2015.



rm(list=ls())
graphics.off()






setwd("/Users/xdu4/Documents/Duxiuxia/Bitbucket/adap-gc_3")

path_to_R_libs <- "/Users/xdu4/Documents/Duxiuxia/R_libs"

if (!is.element("wmtsa", installed.packages(lib.loc=path_to_R_libs)[,1])) {
    install.packages("wmtsa", 
                     repos="http://cran.us.r-project.org",
                     lib=path_to_R_libs, 
                     dependencies=T)
}
require("splus2R", lib.loc=path_to_R_libs)
require("ifultools", lib.loc=path_to_R_libs)
library("wmtsa", lib.loc=path_to_R_libs)


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




source(paste(getwd(), "d_functions/wmtsa/d_wav_sig.R", sep=.Platform$file.sep))
source(paste(getwd(), "d_functions/ifultools/d_it_util.R", sep=.Platform$file.sep))
source(paste(getwd(), "d_functions/splus2R/Swrappers.R", sep=.Platform$file.sep))

s <- create.signalSeries(z, position=list(from=0, by=sampling_interval))

cwt_coefs <- d_wavCWT(s, 
                    scale.range=deltat(s)*c(1, length(s)),
                    n.scale=100,
                    wavelet="gaussian2",
                    shift=5,
                    variance=1)

list_of_attributes <- attributes(cwt_coefs)

def.par <- par(no.readonly = T)
nf <- layout(mat=matrix(c(1,2), nrow=2, ncol=1, byrow=T), heights=c(1,4), widths=c(1,1))
layout.show(nf)
par(mar=c(1,1,1,1))
plot(t, z, type="l", xlim=c(min(t), max(t)))
par(mar=c(1,1,1,1))
image(x=list_of_attributes$time, 
      y=log2(list_of_attributes$scale), 
      z=as.matrix(cwt_coefs), 
      col=heat.colors(24),
      xlab="Time", ylab="log2(scale)")
par(def.par)



plot(cwt_coefs, series=T)


## form CWT tree 
cwt_tree <- wavCWTTree(cwt_coefs)

tree_time_all <- vector(mode="numeric", length=0)
tree_scale_all <- vector(mode="numeric", length=0)


for (i in 1:length(cwt_tree)) {
    tree_time_all <- c(tree_time_all, cwt_tree[[i]]$time)
    tree_scale_all <- c(tree_scale_all, cwt_tree[[i]]$scale)
}

def.par <- par(no.readonly = T)

nf <- layout(mat=matrix(c(1,2), nrow=2, ncol=1, byrow=T))
layout.show(nf)

plot(t, z, cex=0.2)
plot(tree_time_all, log2(tree_scale_all), cex=0.2)
par(def.par)



plot(cwt_tree)







## estimate the peak locations using default 
## scale.range 
cwt_peaks <- wavCWTPeaks(cwt_tree)


plot(t, z, cex=0.2)
points(cwt_peaks$x, cwt_peaks$y, pch=16, col="red", cex=1.2)
