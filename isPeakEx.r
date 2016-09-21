peaks<-function (x, span = 3) 
{
    z <- embed(x, span)
    s <- span%/%2
    result <- max.col(z,ties.method="first") == 1 + s
    c(rep(FALSE, s), result, rep(FALSE, s))
}
embed <-function (x, dimension = 1) 
{
    
    
        n <- length(x)
        if ((dimension < 1) | (dimension > n)) 
            stop("wrong embedding dimension")
        m <- n - dimension + 1
        data <- x[1L:m + rep.int(dimension:1, rep.int(m, dimension)) - 
            1]
        dim(data) <- c(m, dimension)
        return(data)
   
}

parea <- function(f, pt, eps) 
				{
					 x <- f[, 1]
					sum(f[x <= pt + eps & x >= pt - eps, 2])
 				  }
sigma<-function (x, span = 5) 
{

    n <- length(x)
    z <- embed(x, span)
    s <- span%/%2
    msd <- apply(z, 1, sd)
    c(rep(msd[1], s), msd, rep(msd[n - 2 * s], s))
	
}

isPeakEx<-function (f, SoN = 2, span = 81, sm.span = 11, plot = FALSE, 
    add = FALSE, zerothrsh = 2, area.w = 0.003, ratio = 0.2, 
    ...) 
{
    #parea <- function(f, pt, eps = area.w) {
#        x <- f[, 1]
#        sum(f[x <= pt * (1 + eps) & x >= pt * (1 - eps), 2])
#    }
    n <- dim(f)[1]
    mz <- f[, 1]
    lo <- lnn(f[, 2], span = span, sm.span = sm.span)
    sm <- lo$fitted
    ispeak <- peaks(sm, span)
    sigmas <- lo$s
    peak <- ispeak & sm > zerothrsh & sm > SoN * sigmas
	
	if(length(peak[peak])>0)
	{
		#area <- sapply(mz[peak], parea, f = cbind(mz, sm)[!peak,])
	    #As <- rep(NA, length(peak))
	    #As[peak] <- area
		#peak[peak] <- peak[peak] & area/max(area) > ratio
		#pks <- data.frame(peak = peak, smooth = sm, mz = mz, sigmas = sigmas, area = As)
		pks <- data.frame(peak = peak, smooth = sm, mz = mz, sigmas = sigmas)
	}
	else
	{
		pks<-NULL	
	}
    if (plot) 
	{
        if (add) 
			{
				lines(mz, sm, col = "cyan")
				if(length(peak[peak])>0)points(mz[peak], sm[peak], col = "red")
			}
        else 
			if(length(peak[peak])>0)specZoom(pks, ...)
    }
    return(pks)
}

isPeakEx2<-function (f, SoN = 2, span = 81, sm.span = 11, plot = FALSE, 
    add = FALSE, zerothrsh = 2, area.w = 0.003, ratio = 0.2, 
    ...) 
{
   # parea <- function(f, pt, eps = area.w) {
#        x <- f[, 1]
#        sum(f[x <= pt * (1 + eps) & x >= pt * (1 - eps), 2])
#    }
    n <- dim(f)[1]
    mz <- f[, 1]
    #lo <- lnn(f[, 2], span = span, sm.span = sm.span)
    #sm <- lo$fitted
	sm<-f[, 2]
	isMax <- peaks(sm, span)
	noise<-noise(sm,span)
	isSoN <-  sm > SoN * noise
	isZeroThrsh<- sm > zerothrsh
    #sigmas <- lo$s
    isPeak <- isMax & isSoN & isZeroThrsh
	
	if(length(peak[peak])>0)
	{
		area <- sapply(mz[peak], parea, f = cbind(mz, sm)[!peak,])
	    As <- rep(NA, length(peak))
	    As[isPeak] <- area
		isPeak[isPeak] <- isPeak[isPeak] & area/min(area) > ratio
		pks <- data.frame(isPeak=isPeak,isMax = isMax, int = sm, mz = mz, noise = noise, isSoN=isSoN,area = As,isZeroThrsh=isZeroThrsh)
		#pks <- data.frame(peak = peak, smooth = sm, mz = mz, sigmas = sigmas)
	}
	else
	{
		pks<-NULL	
	}
    if (plot) 
	{
        if (add) 
			{
				lines(mz, sm, col = "cyan")
				if(length(peak[peak])>0)points(mz[peak], sm[peak], col = "red")
			}
        else 
			if(length(peak[peak])>0)specZoom(pks, ...)
    }
    return(pks)
}
