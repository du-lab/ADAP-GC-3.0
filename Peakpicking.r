cwt <- function(ms, scales=1, wavelet='mexh') {
	## Check for the wavelet format
	if (wavelet == 'mexh') {
		psi_xval <- seq(-8, 8, length=1024)
		psi <- (2/sqrt(3) * pi^(-0.25)) * (1 - psi_xval^2) *exp(-psi_xval^2/2)
		#plot(psi_xval, psi)
	} else if (is.matrix(wavelet)) {
		if (nrow(wavelet) == 2) {
			psi_xval <- wavelet[1,]
			psi <- wavelet[2,]
		} else if (ncol(wavelet) == 2) {
			psi_xval <- wavelet[,1]
			psi <- wavelet[,2]
		} else {
			stop('Unsupported wavelet format!')
		}
	} else {
		stop('Unsupported wavelet!')
	}
	
	oldLen <- length(ms)
	## To increase the computation effeciency of FFT, extend it as the power of 2
	## because of a strange signal length 21577 makes the FFT very slow!
	ms <- extendNBase(ms, nLevel=NULL, base=2)
	len <- length(ms)
	nbscales <- length(scales)
	wCoefs <- NULL
	
	psi_xval <- psi_xval - psi_xval[1]
	dxval <- psi_xval[2]
	xmax  <- psi_xval[length(psi_xval)]
	for (i in 1:length(scales)) {
		scale.i <- scales[i]
		f <- rep(0, len)
		j <- 1 + floor((0:(scale.i * xmax))/(scale.i * dxval))
		if (length(j) == 1)		j <- c(1, 1)
		lenWave <- length(j)
		f[1:lenWave] <- rev(psi[j]) - mean(psi[j])
		if (length(f) > len) stop(paste('scale', scale.i, 'is too large!'))
		wCoefs.i <- 1/sqrt(scale.i) * convolve(ms, f)
		## Shift the position with half wavelet width
		wCoefs.i <- c(wCoefs.i[(len-floor(lenWave/2) + 1) : len], wCoefs.i[1:(len-floor(lenWave/2))])
		wCoefs <- cbind(wCoefs, wCoefs.i)
	}
	if (length(scales) == 1) wCoefs <- matrix(wCoefs, ncol=1)
	colnames(wCoefs) <- scales
	wCoefs <- wCoefs[1:oldLen,,drop=FALSE]
	return(wCoefs)
}

extendNBase <- function(x, nLevel=1, base=2, ...) {
	
	if (!is.matrix(x)) {
		x <- matrix(x, ncol=1)
	} else if (min(dim(x)) == 1) {
		x <- matrix(x, ncol=1)
	}
	
	nR <- nrow(x)
	if (is.null(nLevel)) {
		nR1 <- nextn(nR, base)		
	} else {
		nR1 <- ceiling(nR / base^nLevel) * base^nLevel		
	}
	if (nR != nR1) {
		x <- extendLength(x, addLength=nR1-nR, ...)
	}
	
	return(x)
}

extendLength <- function(x, addLength=NULL, method=c('reflection', 'open', 'circular'), direction=c('right', 'left', 'both')) {
	if (is.null(addLength)) stop('Please provide the length to be added!')
	if (!is.matrix(x)) x <- matrix(x, ncol=1)	
	method <- match.arg(method)
	direction <- match.arg(direction)
	
	nR <- nrow(x)
	nR1 <- nR + addLength
	if (direction == 'both') {
		left <- right <- addLength
	} else if (direction == 'right') {
		left <- 0
		right <- addLength
	} else if (direction == 'left') {
		left <- addLength
		right <- 0
	}
	
	if (right > 0) {
		x <- switch(method,
				reflection =rbind(x, x[nR:(2 * nR - nR1 + 1), , drop=FALSE]),
				open = rbind(x, matrix(rep(x[nR,], addLength), ncol=ncol(x), byrow=TRUE)),
				circular = rbind(x, x[1:(nR1 - nR),, drop=FALSE]))
	}
	
	if (left > 0) {
		x <- switch(method,
				reflection =rbind(x[addLength:1, , drop=FALSE], x),
				open = rbind(matrix(rep(x[1,], addLength), ncol=ncol(x), byrow=TRUE), x),
				circular = rbind(x[(2 * nR - nR1 + 1):nR,, drop=FALSE], x))
	}
	if (ncol(x) == 1)  x <- as.vector(x)
	
	return(x)
}


MSWtoCWT <- function (curVecInt,scales) {

	z <- cwt(curVecInt,scales=scales,wavelet="mexh")
	z <- Re(z)
	attr(z, "scale") <- scales
	attr(z, "time") <- c(1:length(curVecInt))
	attr(z, "wavelet") <- "gaussian2"
	attr(z, "series") <- curVecInt
	attr(z, "sampling.interval") <- 1
	attr(z, "series.name") <- deparse(substitute(curVecInt))
	attr(z, "n.sample") <- length(curVecInt)
	attr(z, "n.scale") <- length(scales)
	attr(z, "filter.arg") <- 1
	oldClass(z) <- "wavCWT"
	z
}

wavCWT_ModifyByYan<-function (x, scale, wavelet = "gaussian2", shift = 5, variance = 1) 
{
	checkScalarType(wavelet, "character")
	checkScalarType(shift, "numeric")
	checkScalarType(variance, "numeric")
	
	series.name <- deparse(substitute(x))
	sampling.interval <- deltat(x)
	n.scale <- length(scale)
	
#	if (abs(min(scale) - sampling.interval) > .Machine$double.eps) 
#		stop("Minimum scale must be greater than or equal to sampling interval ", 
#				"of the time series")

	times <- time(x)
	
	x <- as.vector(x)
	storage.mode(x) <- "double"

	gauss1 <- c("gaussian1", "gauss1")
	gauss2 <- c("gaussian2", "gauss2", "mexican hat", "sombrero")
	supported.wavelets <- c("haar", gauss1, gauss2, "morlet")
	wavelet <- match.arg(lowerCase(wavelet), supported.wavelets)
	filter <- mutilsFilterTypeContinuous(wavelet)
	
	if (filter == 4) {
		filter.arg <- sqrt(variance)
		wavelet <- "gaussian1"
	}
	else if (filter == 5) {
		filter.arg <- sqrt(variance)
		wavelet <- "gaussian2"
	}
	else if (filter == 6) {
		filter.arg <- shift
		wavelet <- "morlet"
	}
	else if (filter == 7) {
		filter.arg <- 0
		wavelet <- "haar"
		scale <- sampling.interval * unique(round(scale/sampling.interval))
	}
	else stop("Unsupported filter type")
	
	z <- .Call("RS_wavelets_transform_continuous_wavelet", as.numeric(x), 
			as.numeric(sampling.interval), as.integer(filter), as.numeric(filter.arg), 
			as.numeric(scale), COPY = rep(FALSE, 5), CLASSES = c("numeric", 
					"numeric", "integer", "numeric", "numeric"), PACKAGE = "ifultools")
	
	if (wavelet != "morlet") 
		z <- Re(z)
	
	attr(z, "scale") <- scale
	attr(z, "time") <- as.vector(times)
	attr(z, "wavelet") <- wavelet
	attr(z, "series") <- x
	attr(z, "sampling.interval") <- sampling.interval
	attr(z, "series.name") <- series.name
	attr(z, "n.sample") <- length(x)
	attr(z, "n.scale") <- n.scale
	attr(z, "filter.arg") <- filter.arg
	oldClass(z) <- "wavCWT"
	z
}




wavCWT_Modify_V1<-function (x, scale.range = deltat(x) * c(1, length(x)), n.scale = 100, 
		wavelet = "gaussian2", shift = 5, variance = 1) 
{
	checkVectorType(scale.range, "numeric")
	checkScalarType(n.scale, "integer")
	checkScalarType(wavelet, "character")
	checkScalarType(shift, "numeric")
	checkScalarType(variance, "numeric")
	checkRange(n.scale, c(1, Inf))
	
	series.name <- deparse(substitute(x))
	if (length(scale.range) != 2) 
		stop("scale.range must be a two-element numeric vector")
	if (variance <= 0) 
		stop("variance input must be positive")
	
	sampling.interval <- deltat(x)
	
	# determine scales
	octave <- logb(scale.range, 2) # log2
	scale <- ifelse1(n.scale > 1, 2^c(octave[1] + seq(0, n.scale - 2) * diff(octave)/(floor(n.scale) - 1), octave[2]), scale.range[1])
	scale <- unique(round(scale/sampling.interval) * sampling.interval)
	n.scale <- length(scale)
	
	if (abs(min(scale) - sampling.interval) > .Machine$double.eps) 
		stop("Minimum scale must be greater than or equal to sampling interval ", 
				"of the time series")
	if (inherits(x, "signalSeries")) 
		times <- as(x@positions, "numeric")
	else times <- time(x)
	
	x <- as.vector(x)
	storage.mode(x) <- "double"
	gauss1 <- c("gaussian1", "gauss1")
	gauss2 <- c("gaussian2", "gauss2", "mexican hat", "sombrero")
	supported.wavelets <- c("haar", gauss1, gauss2, "morlet")
	wavelet <- match.arg(lowerCase(wavelet), supported.wavelets)
	filter <- mutilsFilterTypeContinuous(wavelet)
	
	if (filter == 4) {
		filter.arg <- sqrt(variance)
		wavelet <- "gaussian1"
	}
	else if (filter == 5) {
		filter.arg <- sqrt(variance)
		wavelet <- "gaussian2"
	}
	else if (filter == 6) {
		filter.arg <- shift
		wavelet <- "morlet"
	}
	else if (filter == 7) {
		filter.arg <- 0
		wavelet <- "haar"
		scale <- sampling.interval * unique(round(scale/sampling.interval))
	}
	else stop("Unsupported filter type")
	
	z <- .Call("RS_wavelets_transform_continuous_wavelet", as.numeric(x), 
			as.numeric(sampling.interval), as.integer(filter), as.numeric(filter.arg), 
			as.numeric(scale), COPY = rep(FALSE, 5), CLASSES = c("numeric", 
					"numeric", "integer", "numeric", "numeric"), PACKAGE = "ifultools")
	
	if (wavelet != "morlet") 
		z <- Re(z)
	
	attr(z, "scale") <- scale
	attr(z, "time") <- as.vector(times)
	attr(z, "wavelet") <- wavelet
	attr(z, "series") <- x
	attr(z, "sampling.interval") <- sampling.interval
	attr(z, "series.name") <- series.name
	attr(z, "n.sample") <- length(x)
	attr(z, "n.scale") <- n.scale
	attr(z, "filter.arg") <- filter.arg
	oldClass(z) <- "wavCWT"
	z
}



wavCWTTree_ModifyByYan<-function (x, n.octave.min = 1, tolerance = 0, type = "maxima",noise.fun = "quantile") 
{
	"WTMM" <- function(x, tolerance = NULL, type = "maxima") {
		
		if (!is(x, "wavCWT")) 
			stop("Input object must be of class wavCWT")
		
		x.attr <- attributes(x)
		times <- x.attr$time
		scales <- x.attr$scale
		n.sample <- x.attr$n.sample
		series <- x.attr$series
		
		if (is.null(tolerance)) {
			tolerance <- mad(Mod(x[, 1]))/scales
		}
		if (length(tolerance) < length(scales)) 
			tolerance <- tolerance[1]/sqrt(scales)
		
		wtmmz <- .Call("RS_wavelets_transform_continuous_wavelet_modulus_maxima", 
				as.matrix(x) + (0+0i), tolerance, mutilsTransformPeakType(type), 
				CLASSES = c("matrix", "numeric", "integer"), COPY = rep(FALSE, 
						3), PACKAGE = "ifultools")
		
		# wtmmz returns the extrema positions starting from 0
		z <- matrix(0, nrow = nrow(x), ncol = ncol(x))
		z[matrix(unlist(wtmmz), ncol = 2) + 1] <- 1
		z
	}# WTMM func
	
	"wtmmBranches" <- function(wtmm, extrema.mask, times, scales, 
			span.min = 5, gap.max = 3, skip = NULL, sampling.interval = 1) {
		
		scales <- as.integer(scales/sampling.interval)
		n.scale <- ncol(extrema.mask)
		n.sample <- nrow(extrema.mask)
		
		if (is.null(scales)) 
			scales <- 1:n.scale
		
		iwtmm <- which(extrema.mask[, n.scale] > 0) # the last scale ?
		
		if (length(iwtmm)!=0) {  # ?
			
			iscale <- seq(n.scale - 1, 1, -1)
			tree <- as.list(iwtmm)
			names(tree) <- iwtmm
			peakStatus <- as.list(rep(0, length(iwtmm)))
			names(peakStatus) <- iwtmm
			orphanRidgeList <- NULL
			orphanRidgeName <- NULL
			n.level <- length(iscale)
			
			for (j in seq(n.level)) {
				iscale.j <- iscale[j]
				scale.j <- scales[iscale.j]
				
				if (length(iwtmm) == 0) {
					iwtmm <- which(extrema.mask[, iscale.j] > 0)
					next
				}
				
				span <- scale.j * 2 + 1
				if (span < span.min) 
					span <- span.min
				
				remove.j <- selPeak.j <- NULL
				
				for (k in seq(along = iwtmm)) { # to count the length
					itime <- iwtmm[k]
					itime.start <- itime - span
					if (itime.start < 1) 
						itime.start <- 1
					itime.end <- itime + span
					if (itime.end > n.sample) 
						itime.end <- n.sample
					
					itime.candidates <- which(extrema.mask[itime.start:itime.end, 
									iscale.j] > 0) + itime.start - 1
					
					if (length(itime.candidates) == 0) {
						status.k <- peakStatus[[as.character(itime)]]
						
						if (length(status.k)>0&length(scale.j)>0) {
							if (status.k > gap.max & scale.j >= 2) {
								temp <- tree[[as.character(itime)]]
								orphanRidgeList <- c(orphanRidgeList, list(temp[1:(length(temp) - 
																	status.k)]))
								orphanRidgeName <- c(orphanRidgeName, paste(iscale.j + 
														status.k + 1, itime, sep = "_"))
								remove.j <- c(remove.j, as.character(itime))
								next
							}
							else {
								itime.candidates <- itime
								peakStatus[[as.character(itime)]] <- status.k + 
										1
							}
						}
					}
					else {
						peakStatus[[as.character(itime)]] <- 0
						
						if (length(itime.candidates) >= 2) 
							itime.candidates <- itime.candidates[which.min(abs(itime.candidates - 
															itime))]
					}
					
					tree[[as.character(itime)]] <- c(tree[[as.character(itime)]], itime.candidates)
					selPeak.j <- c(selPeak.j, itime.candidates)
				} # for loop
				
				if (length(remove.j) > 0) {
					bad.tree <- which(is.element(names(tree), remove.j))
					tree <- tree[-bad.tree]
					peakStatus <- peakStatus[-bad.tree]
				}
				
				dupPeak.j <- unique(selPeak.j[duplicated(selPeak.j)])
				
				if (length(dupPeak.j) > 0) {
					bad.tree <- NULL
					for (dupPeak.jk in dupPeak.j) {
						selInd <- which(selPeak.j == dupPeak.jk)
						selLen <- sapply(tree[selInd], length)
						bad.tree.jk <- which.max(selLen)
						bad.tree <- c(bad.tree, selInd[-bad.tree.jk])
						orphanRidgeList <- c(orphanRidgeList, tree[bad.tree.jk])
						orphanRidgeName <- c(orphanRidgeName, paste(iscale.j, 
										selPeak.j[bad.tree.jk], sep = "_"))
					}
					selPeak.j <- selPeak.j[-bad.tree]
					tree <- tree[-bad.tree]
					peakStatus <- peakStatus[-bad.tree]
				}
				
				names(tree) <- selPeak.j
				names(peakStatus) <- selPeak.j
				
				if (scale.j >= 2) {
					maxInd.next <- which(extrema.mask[, iscale.j] >0)
					unSelPeak.j <- maxInd.next[!is.element(maxInd.next, selPeak.j)]
					newPeak.j <- as.list(unSelPeak.j)
					names(newPeak.j) <- unSelPeak.j
					tree <- c(tree, newPeak.j)
					iwtmm <- c(selPeak.j, unSelPeak.j)
					newPeakStatus <- as.list(rep(0, length(newPeak.j)))
					names(newPeakStatus) <- newPeak.j
					peakStatus <- c(peakStatus, newPeakStatus)
				}
				else {
					iwtmm <- selPeak.j
				}
			}
			
			if (length(tree)!=0) {
				names(tree) <- paste(1, names(tree), sep = "_")
				names(orphanRidgeList) <- orphanRidgeName
				tree <- c(tree, orphanRidgeList)
				tree <- lapply(tree, rev)
				tree <- tree[unique(names(tree))]
				tree <- lapply(seq(along = tree), function(i, tree, iscale.min, times, scales, wtmm) {
							itime <- tree[[i]]
							iscale <- seq(iscale.min[i], length = length(itime))
							list(itime = itime, iscale = iscale, time = times[itime], 
									scale = scales[iscale], extrema = wtmm[cbind(itime, 
													iscale)])
						}, tree = tree, iscale.min = as.integer(gsub("_.*", "", 
										names(tree))), times = times, scales = scales * sampling.interval, 
						wtmm = wtmm)
				iflat <- lapply(tree, function(x, nr) (x$iscale - 1) *nr + x$itime, nr = nrow(wtmm))
				flatset <- iflat[[1]]
				bad <- NULL
				if(length(iflat)>1)
				{
					for (i in seq(2, length(iflat))) {
						if (any(is.element(iflat[[i]], flatset))) 
							bad <- c(bad, i)
						else flatset <- c(flatset, iflat[[i]])
					}
				}
				
				if (length(bad) > 0) 
					tree <- tree[-bad]
				tree
			}# tree length checking		
		}
	}# wtmmBranches function
	
	x.attr <- attributes(x)
	times <- x.attr$time  # scan
	scales <- x.attr$scale
	n.sample <- x.attr$n.sample # scan number
	sampling.interval <- x.attr$sampling.interval
	series <- x.attr$series # Intensity
	
	#border.times <- range(times) + sampling.interval * c(1, -1)
	extrema.mask <- WTMM(x, tolerance = tolerance, type = type)
	
	if (!identical(dim(x), dim(extrema.mask))) 
		stop("Input WTMM dimenions do not match those of the input CWT matrix")
	
	z <- wtmmBranches(ifelse1(is.complex(x), Mod(as.matrix(x)), 
					as.matrix(x)), extrema.mask, times, scales, sampling.interval = sampling.interval)
	
	if (!is.null(z)) {
		n.scale <- length(scales)
		n.octave <- log2(max(scales)/min(scales))
		n.voice <- (n.scale - 1)/n.octave
		n.scale.min <- as.integer(n.voice * n.octave.min)
		
		# filtering by the minimal scale number
		good <- which(unlist(lapply(z, function(x, n.scale.min) length(x[[1]]) > n.scale.min, n.scale.min = n.scale.min)))
		
		if (length(good)!=0) {
			z <- z[good]
			# modified by yan output the index @ highest intensity
			#endtime <- unlist(lapply(z, function(x, iscale) x$itime[iscale], iscale = which.min(scales)))
			endtime <- unlist(lapply(z, function(x) {ipos = which.max(series[x$itime]) 
								x$itime[ipos]}))
			isort <- order(endtime)
			endtime <- sort(endtime)
			z <- z[isort]
			#names(z) <- seq(z)
		
			noise <- x[, 1]
			noise.min <- quantile(abs(noise), prob = 0.05)
			noise.span <- max(0.01 * diff(range(times)), 5 * sampling.interval)
		
			noise.levels <- unlist(lapply(endtime, function(x, noise.fun, 
								times, times.range, noise, noise.min, noise.span) {
							time.start <- x - noise.span
							if (time.start < times.range[1]) 
								time.start <- times.range[1]
							time.end <- x + noise.span
							#if (time.end < times.range[2]) #?
							if (time.end > times.range[2]) 
								time.end <- times.range[2]
							ix <- which(times >= time.start & times <= time.end)
							noise.local <- noise.fun(abs(noise[ix]))
							if (noise.local < noise.min) 
								noise.local <- noise.min
							noise.local
						}, noise.fun = switch(noise.fun, quantile = function(x) {
									quantile(x, probs = 0.95)
								}, sd = sd, mad = function(x) {
									mad(x, center = 0)
								}), times = times, times.range = range(times), noise = noise, # noise is the coeff at smallest scale
						noise.min = noise.min, noise.span = noise.span))
		
		tmpargs <- lapply(z, function(x) unlist(lapply(x, function(x,imax) x[imax], imax = which.max(x$extrema))))
		peaks <- data.frame(do.call("rbind", tmpargs))
		peak.snr <- round(peaks[["extrema"]]/noise.levels)
		peaks <- data.frame(pkInd=endtime,StN=peak.snr)
		peaks
		}
	}
}

wavCWTTree_ModifyByYan_Min<-function (x, n.octave.min = 1, tolerance = 0, type = "minimal") 
{
	
	"WTMM" <- function(x, tolerance = NULL, type = "minimal") {
		
		if (!is(x, "wavCWT")) 
			stop("Input object must be of class wavCWT")
		
		x.attr <- attributes(x)
		times <- x.attr$time
		scales <- x.attr$scale
		n.sample <- x.attr$n.sample
		series <- x.attr$series
		
		if (is.null(tolerance)) {
			tolerance <- mad(Mod(x[, 1]))/scales
		}
		if (length(tolerance) < length(scales)) 
			tolerance <- tolerance[1]/sqrt(scales)
		
		wtmmz <- .Call("RS_wavelets_transform_continuous_wavelet_modulus_maxima", 
				as.matrix(x) + (0+0i), tolerance, mutilsTransformPeakType(type), 
				CLASSES = c("matrix", "numeric", "integer"), COPY = rep(FALSE, 
						3), PACKAGE = "ifultools")
		
		# wtmmz returns the extrema positions starting from 0
		z <- matrix(0, nrow = nrow(x), ncol = ncol(x))
		z[matrix(unlist(wtmmz), ncol = 2) + 1] <- 1
		z
	}# WTMM func
	
	"wtmmBranches" <- function(wtmm, extrema.mask, times, scales, 
			span.min = 5, gap.max = 3, skip = NULL, sampling.interval = 1) {
		
		scales <- as.integer(scales/sampling.interval)
		n.scale <- ncol(extrema.mask)
		n.sample <- nrow(extrema.mask)
		
		if (is.null(scales)) 
			scales <- 1:n.scale
		
		iwtmm <- which(extrema.mask[, n.scale] > 0) # the last scale ?
		
		if (length(iwtmm)!=0) {  # ?
			
			iscale <- seq(n.scale - 1, 1, -1)
			tree <- as.list(iwtmm)
			names(tree) <- iwtmm
			peakStatus <- as.list(rep(0, length(iwtmm)))
			names(peakStatus) <- iwtmm
			orphanRidgeList <- NULL
			orphanRidgeName <- NULL
			n.level <- length(iscale)
			
			for (j in seq(n.level)) {
				iscale.j <- iscale[j]
				scale.j <- scales[iscale.j]
				
				if (length(iwtmm) == 0) {
					iwtmm <- which(extrema.mask[, iscale.j] > 0)
					next
				}
				
				span <- scale.j * 2 + 1
				if (span < span.min) 
					span <- span.min
				
				remove.j <- selPeak.j <- NULL
				
				for (k in seq(along = iwtmm)) { # to count the length
					itime <- iwtmm[k]
					itime.start <- itime - span
					if (itime.start < 1) 
						itime.start <- 1
					itime.end <- itime + span
					if (itime.end > n.sample) 
						itime.end <- n.sample
					
					itime.candidates <- which(extrema.mask[itime.start:itime.end, 
									iscale.j] > 0) + itime.start - 1
					
					if (length(itime.candidates) == 0) {
						status.k <- peakStatus[[as.character(itime)]]
						
						if (length(status.k)>0&length(scale.j)>0) {
							if (status.k > gap.max & scale.j >= 2) {
								temp <- tree[[as.character(itime)]]
								orphanRidgeList <- c(orphanRidgeList, list(temp[1:(length(temp) - 
																	status.k)]))
								orphanRidgeName <- c(orphanRidgeName, paste(iscale.j + 
														status.k + 1, itime, sep = "_"))
								remove.j <- c(remove.j, as.character(itime))
								next
							}
							else {
								itime.candidates <- itime
								peakStatus[[as.character(itime)]] <- status.k + 
										1
							}
						}
					}
					else {
						peakStatus[[as.character(itime)]] <- 0
						
						if (length(itime.candidates) >= 2) 
							itime.candidates <- itime.candidates[which.min(abs(itime.candidates - 
															itime))]
					}
					
					tree[[as.character(itime)]] <- c(tree[[as.character(itime)]], itime.candidates)
					selPeak.j <- c(selPeak.j, itime.candidates)
				} # for loop
				
				if (length(remove.j) > 0) {
					bad.tree <- which(is.element(names(tree), remove.j))
					tree <- tree[-bad.tree]
					peakStatus <- peakStatus[-bad.tree]
				}
				
				dupPeak.j <- unique(selPeak.j[duplicated(selPeak.j)])
				
				if (length(dupPeak.j) > 0) {
					bad.tree <- NULL
					for (dupPeak.jk in dupPeak.j) {
						selInd <- which(selPeak.j == dupPeak.jk)
						selLen <- sapply(tree[selInd], length)
						bad.tree.jk <- which.max(selLen)
						bad.tree <- c(bad.tree, selInd[-bad.tree.jk])
						orphanRidgeList <- c(orphanRidgeList, tree[bad.tree.jk])
						orphanRidgeName <- c(orphanRidgeName, paste(iscale.j, 
										selPeak.j[bad.tree.jk], sep = "_"))
					}
					selPeak.j <- selPeak.j[-bad.tree]
					tree <- tree[-bad.tree]
					peakStatus <- peakStatus[-bad.tree]
				}
				
				names(tree) <- selPeak.j
				names(peakStatus) <- selPeak.j
				
				if (scale.j >= 2) {
					maxInd.next <- which(extrema.mask[, iscale.j] >0)
					unSelPeak.j <- maxInd.next[!is.element(maxInd.next, selPeak.j)]
					newPeak.j <- as.list(unSelPeak.j)
					names(newPeak.j) <- unSelPeak.j
					tree <- c(tree, newPeak.j)
					iwtmm <- c(selPeak.j, unSelPeak.j)
					newPeakStatus <- as.list(rep(0, length(newPeak.j)))
					names(newPeakStatus) <- newPeak.j
					peakStatus <- c(peakStatus, newPeakStatus)
				}
				else {
					iwtmm <- selPeak.j
				}
			}
			
			if (length(tree)!=0) {
				names(tree) <- paste(1, names(tree), sep = "_")
				names(orphanRidgeList) <- orphanRidgeName
				tree <- c(tree, orphanRidgeList)
				tree <- lapply(tree, rev)
				tree <- tree[unique(names(tree))]
				tree <- lapply(seq(along = tree), function(i, tree, iscale.min, times, scales, wtmm) {
							itime <- tree[[i]]
							iscale <- seq(iscale.min[i], length = length(itime))
							list(itime = itime, iscale = iscale, time = times[itime], 
									scale = scales[iscale], extrema = wtmm[cbind(itime, 
													iscale)])
						}, tree = tree, iscale.min = as.integer(gsub("_.*", "", 
										names(tree))), times = times, scales = scales * sampling.interval, 
						wtmm = wtmm)
				iflat <- lapply(tree, function(x, nr) (x$iscale - 1) *nr + x$itime, nr = nrow(wtmm))
				flatset <- iflat[[1]]
				bad <- NULL
				if(length(iflat)>1)
				{
					for (i in seq(2, length(iflat))) {
						if (any(is.element(iflat[[i]], flatset))) 
							bad <- c(bad, i)
						else flatset <- c(flatset, iflat[[i]])
					}
				}
				
				if (length(bad) > 0) 
					tree <- tree[-bad]
				tree
			}# tree length checking		
		}
	}# wtmmBranches function
	
	
	x.attr <- attributes(x)
	times <- x.attr$time
	scales <- x.attr$scale
	n.sample <- x.attr$n.sample
	sampling.interval <- x.attr$sampling.interval
	series <- x.attr$series
	
	border.times <- range(times) + sampling.interval * c(1, -1)
	extrema.mask <- WTMM(x, tolerance = tolerance, type = type)
	
	if (!identical(dim(x), dim(extrema.mask))) 
		stop("Input WTMM dimenions do not match those of the input CWT matrix")
	
	z <- wtmmBranches(ifelse1(is.complex(x), Mod(as.matrix(x)), 
					as.matrix(x)), extrema.mask, times, scales, sampling.interval = sampling.interval)
	
	if (!is.null(z)) {
		n.scale <- length(scales)
		n.octave <- log2(max(scales)/min(scales))
		n.voice <- (n.scale - 1)/n.octave
		n.scale.min <- as.integer(n.voice * n.octave.min)
		
		# filtering
		good <- which(unlist(lapply(z, function(x, n.scale.min) length(x[[1]]) > n.scale.min, n.scale.min = n.scale.min)))
		
		if (length(good)!=0) {
			z <- z[good]
			# modified by yan output the index @ highest intensity
			#endtime <- unlist(lapply(z, function(x, iscale) x$itime[iscale], iscale = which.min(scales)))
			endtime <- unlist(lapply(z, function(x) {ipos = which.min(series[x$itime]) 
								x$itime[ipos]}))
			isort <- order(endtime)
			z <- z[isort]
			names(z) <- seq(z)
			
			attr(z, "iendtime") <- endtime[isort]
			attr(z, "endtime") <- times[endtime[isort]] # earliest time
			attr(z, "time") <- times
			attr(z, "scale") <- scales
			attr(z, "extrema.mask") <- extrema.mask
			attr(z, "noise") <- x[, 1]  # the smallest scale
			attr(z, "branch.hist") <- colSums(extrema.mask * abs(x))
			attr(z, "wavelet") <- attr(x, "wavelet")
			attr(z, "filter.arg") <- attr(x, "filter.arg")
			attr(z, "series.name") <- attr(x, "series.name")
			attr(z, "series") <- attr(x, "series")
			attr(z, "sampling.interval") <- attr(x, "sampling.interval")
			oldClass(z) <- "wavCWTTree"
			z
		}
	}
}

wavCWTTree_Modify_V1<-function (x, n.octave.min = 1, tolerance = 0, type = "maxima") 
{
	
	"WTMM" <- function(x, tolerance = NULL, type = "maxima") {
		
		if (!is(x, "wavCWT")) 
			stop("Input object must be of class wavCWT")
		
		x.attr <- attributes(x)
		times <- x.attr$time
		scales <- x.attr$scale
		n.sample <- x.attr$n.sample
		series <- x.attr$series
		
		if (is.null(tolerance)) {
			tolerance <- mad(Mod(x[, 1]))/scales
		}
		if (length(tolerance) < length(scales)) 
			tolerance <- tolerance[1]/sqrt(scales)
		
		wtmmz <- .Call("RS_wavelets_transform_continuous_wavelet_modulus_maxima", 
				as.matrix(x) + (0+0i), tolerance, mutilsTransformPeakType(type), 
				CLASSES = c("matrix", "numeric", "integer"), COPY = rep(FALSE, 
						3), PACKAGE = "ifultools")
		
		# wtmmz returns the extrema positions starting from 0
		z <- matrix(0, nrow = nrow(x), ncol = ncol(x))
		z[matrix(unlist(wtmmz), ncol = 2) + 1] <- 1
		z
	}# WTMM func
	
	"wtmmBranches" <- function(wtmm, extrema.mask, times, scales, 
			span.min = 5, gap.max = 3, skip = NULL, sampling.interval = 1) {
		
		scales <- as.integer(scales/sampling.interval)
		n.scale <- ncol(extrema.mask)
		n.sample <- nrow(extrema.mask)
		
		if (is.null(scales)) 
			scales <- 1:n.scale
		
		iwtmm <- which(extrema.mask[, n.scale] > 0) # the last scale ?
		
		if (length(iwtmm)!=0) {  # ?
			
			iscale <- seq(n.scale - 1, 1, -1)
			tree <- as.list(iwtmm)
			names(tree) <- iwtmm
			peakStatus <- as.list(rep(0, length(iwtmm)))
			names(peakStatus) <- iwtmm
			orphanRidgeList <- NULL
			orphanRidgeName <- NULL
			n.level <- length(iscale)
			
			for (j in seq(n.level)) {
				iscale.j <- iscale[j]
				scale.j <- scales[iscale.j]
				
				if (length(iwtmm) == 0) {
					iwtmm <- which(extrema.mask[, iscale.j] > 0)
					next
				}
				
				span <- scale.j * 2 + 1
				if (span < span.min) 
					span <- span.min
				
				remove.j <- selPeak.j <- NULL
				
				for (k in seq(along = iwtmm)) { # to count the length
					itime <- iwtmm[k]
					itime.start <- itime - span
					if (itime.start < 1) 
						itime.start <- 1
					itime.end <- itime + span
					if (itime.end > n.sample) 
						itime.end <- n.sample
					
					itime.candidates <- which(extrema.mask[itime.start:itime.end, 
									iscale.j] > 0) + itime.start - 1
					
					if (length(itime.candidates) == 0) {
						status.k <- peakStatus[[as.character(itime)]]
						
						if (length(status.k)>0&length(scale.j)>0) {
							if (status.k > gap.max & scale.j >= 2) {
								temp <- tree[[as.character(itime)]]
								orphanRidgeList <- c(orphanRidgeList, list(temp[1:(length(temp) - 
																	status.k)]))
								orphanRidgeName <- c(orphanRidgeName, paste(iscale.j + 
														status.k + 1, itime, sep = "_"))
								remove.j <- c(remove.j, as.character(itime))
								next
							}
							else {
								itime.candidates <- itime
								peakStatus[[as.character(itime)]] <- status.k + 
										1
							}
						}
					}
					else {
						peakStatus[[as.character(itime)]] <- 0
						
						if (length(itime.candidates) >= 2) 
							itime.candidates <- itime.candidates[which.min(abs(itime.candidates - 
															itime))]
					}
					
					tree[[as.character(itime)]] <- c(tree[[as.character(itime)]], itime.candidates)
					selPeak.j <- c(selPeak.j, itime.candidates)
				} # for loop
				
				if (length(remove.j) > 0) {
					bad.tree <- which(is.element(names(tree), remove.j))
					tree <- tree[-bad.tree]
					peakStatus <- peakStatus[-bad.tree]
				}
				
				dupPeak.j <- unique(selPeak.j[duplicated(selPeak.j)])
				
				if (length(dupPeak.j) > 0) {
					bad.tree <- NULL
					for (dupPeak.jk in dupPeak.j) {
						selInd <- which(selPeak.j == dupPeak.jk)
						selLen <- sapply(tree[selInd], length)
						bad.tree.jk <- which.max(selLen)
						bad.tree <- c(bad.tree, selInd[-bad.tree.jk])
						orphanRidgeList <- c(orphanRidgeList, tree[bad.tree.jk])
						orphanRidgeName <- c(orphanRidgeName, paste(iscale.j, 
										selPeak.j[bad.tree.jk], sep = "_"))
					}
					selPeak.j <- selPeak.j[-bad.tree]
					tree <- tree[-bad.tree]
					peakStatus <- peakStatus[-bad.tree]
				}
				
				names(tree) <- selPeak.j
				names(peakStatus) <- selPeak.j
				
				if (scale.j >= 2) {
					maxInd.next <- which(extrema.mask[, iscale.j] >0)
					unSelPeak.j <- maxInd.next[!is.element(maxInd.next, selPeak.j)]
					newPeak.j <- as.list(unSelPeak.j)
					names(newPeak.j) <- unSelPeak.j
					tree <- c(tree, newPeak.j)
					iwtmm <- c(selPeak.j, unSelPeak.j)
					newPeakStatus <- as.list(rep(0, length(newPeak.j)))
					names(newPeakStatus) <- newPeak.j
					peakStatus <- c(peakStatus, newPeakStatus)
				}
				else {
					iwtmm <- selPeak.j
				}
			}
			
			if (length(tree)!=0) {
				names(tree) <- paste(1, names(tree), sep = "_")
				names(orphanRidgeList) <- orphanRidgeName
				tree <- c(tree, orphanRidgeList)
				tree <- lapply(tree, rev)
				tree <- tree[unique(names(tree))]
				tree <- lapply(seq(along = tree), function(i, tree, iscale.min, times, scales, wtmm) {
							itime <- tree[[i]]
							iscale <- seq(iscale.min[i], length = length(itime))
							list(itime = itime, iscale = iscale, time = times[itime], 
									scale = scales[iscale], extrema = wtmm[cbind(itime, 
													iscale)])
						}, tree = tree, iscale.min = as.integer(gsub("_.*", "", 
										names(tree))), times = times, scales = scales * sampling.interval, 
						wtmm = wtmm)
				iflat <- lapply(tree, function(x, nr) (x$iscale - 1) *nr + x$itime, nr = nrow(wtmm))
				flatset <- iflat[[1]]
				bad <- NULL
				if(length(iflat)>1)
				{
					for (i in seq(2, length(iflat))) {
						if (any(is.element(iflat[[i]], flatset))) 
							bad <- c(bad, i)
						else flatset <- c(flatset, iflat[[i]])
					}
				}
				
				if (length(bad) > 0) 
					tree <- tree[-bad]
				tree
			}# tree length checking		
		}
	}# wtmmBranches function
	
	
	x.attr <- attributes(x)
	times <- x.attr$time # Intensity
	scales <- x.attr$scale
	n.sample <- x.attr$n.sample # scan number
	sampling.interval <- x.attr$sampling.interval
	
	#border.times <- range(times) + sampling.interval * c(1, -1)
	extrema.mask <- WTMM(x, tolerance = tolerance, type = type)
	
	if (!identical(dim(x), dim(extrema.mask))) 
		stop("Input WTMM dimenions do not match those of the input CWT matrix")
	
	z <- wtmmBranches(ifelse1(is.complex(x), Mod(as.matrix(x)), 
					as.matrix(x)), extrema.mask, times, scales, sampling.interval = sampling.interval)
	
	if (!is.null(z)) {
		n.scale <- length(scales)
		n.octave <- log2(max(scales)/min(scales))
		n.voice <- (n.scale - 1)/n.octave
		n.scale.min <- as.integer(n.voice * n.octave.min)
		
		# filtering using n.scale.min
		good <- which(unlist(lapply(z, function(x, n.scale.min) length(x[[1]]) > n.scale.min, n.scale.min = n.scale.min)))
		
		if (length(good)!=0) {
			z <- z[good]
			endtime <- unlist(lapply(z, function(x, iscale) x$itime[iscale], iscale = which.min(scales)))
			isort <- order(endtime)
			z <- z[isort]
			names(z) <- seq(z)
			
			attr(z, "iendtime") <- endtime[isort]
			attr(z, "endtime") <- times[endtime[isort]] # earliest time
			attr(z, "time") <- times
			attr(z, "scale") <- scales
			attr(z, "extrema.mask") <- extrema.mask
			attr(z, "noise") <- x[, 1]  # the smallest scale
			attr(z, "branch.hist") <- colSums(extrema.mask * abs(x))
			attr(z, "wavelet") <- attr(x, "wavelet")
			attr(z, "filter.arg") <- attr(x, "filter.arg")
			attr(z, "series.name") <- attr(x, "series.name")
			attr(z, "series") <- attr(x, "series")
			attr(z, "sampling.interval") <- attr(x, "sampling.interval")
			oldClass(z) <- "wavCWTTree"
			z
		}
	}
}

wavCWTPeaks_Modify <- function( x, scale.range = NULL, noise.min = NULL, noise.span = NULL, noise.fun = "quantile")
{
	
	if (!is(x, "wavCWTTree")) 
		stop("Input must be an object of class wavCWTTree")
	
	xatt <- attributes(x)
	endtimes <- attr(x, "endtime")
	times <- attr(x, "time")
	scale <- attr(x, "scale")
	noise <- attr(x, "noise")
	wavelet <- attr(x, "wavelet")
	series <- attr(x, "series")
	branch.hist <- attr(x, "branch.hist")
	sampling.interval <- abs(diff(times[1:2]))
	
#	if (!is.element(wavelet, "gaussian2")) 
#		stop("Only CWT developed using the Mexican hat (gaussian2) filter are supported")
	if (is.null(noise.min)) 
		noise.min <- quantile(abs(attr(x, "noise")), prob = 0.05)
	#if (is.null(scale.range)) # no use currently
		#scale.range <- scale[range(which(branch.hist > quantile(branch.hist,prob = 0.8)))]
	if (is.null(noise.span)) 
		noise.span <- max(0.01 * diff(range(times)), 5 * sampling.interval)
	
	noise.levels <- unlist(lapply(endtimes, function(x, noise.fun, 
							times, times.range, noise, noise.min, noise.span) {
						time.start <- x - noise.span
						if (time.start < times.range[1]) 
							time.start <- times.range[1]
						time.end <- x + noise.span
						#if (time.end < times.range[2]) #?
						if (time.end > times.range[2]) 
							time.end <- times.range[2]
						ix <- which(times >= time.start & times <= time.end)
						noise.local <- noise.fun(abs(noise[ix]))
						if (noise.local < noise.min) 
							noise.local <- noise.min
						noise.local
					}, noise.fun = switch(noise.fun, quantile = function(x) {
								quantile(x, probs = 0.95)
							}, sd = sd, mad = function(x) {
								mad(x, center = 0)
							}), times = times, times.range = range(times), noise = noise, # noise is the coeff at smallest scale
					noise.min = noise.min, noise.span = noise.span))
	
	tmpargs <- lapply(x, function(x) unlist(lapply(x, function(x,imax) x[imax], imax = which.max(x$extrema))))
	peaks <- data.frame(do.call("rbind", tmpargs))
	peaks <- cbind(data.frame(branch = row.names(peaks)), peaks,data.frame(iendtime = attr(x, "iendtime")))
	peak.snr <- round(peaks[["extrema"]]/noise.levels)
	if(nrow(peaks)>0) row.names(peaks) <- as.character(seq(nrow(peaks)))
	
	z <- list(x = times[peaks$iendtime], y = series[peaks$iendtime])
	attr(z, "peaks") <- peaks
	attr(z, "scale.range") <- scale.range
	#attr(z, "length.min") <- length.min
	attr(z, "noise.span") <- noise.span
	attr(z, "noise.fun") <- noise.fun
	attr(z, "noise.min") <- noise.min
	attr(z, "snr") <- peak.snr
	z
}

wavCWTPeaks_Modify_V1 <- function( x, snr.min = 3, scale.range = NULL, length.min = 10, 
		noise.span = NULL, noise.fun = "quantile", noise.min = NULL)
{
	
	if (!is(x, "wavCWTTree")) 
		stop("Input must be an object of class wavCWTTree")
	
	xatt <- attributes(x)
	endtimes <- attr(x, "endtime")
	times <- attr(x, "time")
	scale <- attr(x, "scale")
	noise <- attr(x, "noise")
	wavelet <- attr(x, "wavelet")
	series <- attr(x, "series")
	branch.hist <- attr(x, "branch.hist")
	sampling.interval <- abs(diff(times[1:2]))
	
	if (!is.element(wavelet, "gaussian2")) 
		stop("Only CWT developed using the Mexican hat (gaussian2) filter are supported")
	if (is.null(noise.min)) 
		noise.min <- quantile(abs(attr(x, "noise")), prob = 0.05)
	if (is.null(scale.range)) # no use currently
		scale.range <- scale[range(which(branch.hist > quantile(branch.hist,prob = 0.8)))]
	if (is.null(noise.span)) 
		noise.span <- max(0.01 * diff(range(times)), 5 * sampling.interval)
	
	noise.levels <- unlist(lapply(endtimes, function(x, noise.fun, 
							times, times.range, noise, noise.min, noise.span) {
						time.start <- x - noise.span
						if (time.start < times.range[1]) 
							time.start <- times.range[1]
						time.end <- x + noise.span
						if (time.end < times.range[2]) #?
							time.end <- times.range[2]
						ix <- which(times >= time.start & times <= time.end)
						noise.local <- noise.fun(abs(noise[ix]))
						if (noise.local < noise.min) 
							noise.local <- noise.min
						noise.local
					}, noise.fun = switch(noise.fun, quantile = function(x) {
								quantile(x, probs = 0.95)
							}, sd = sd, mad = function(x) {
								mad(x, center = 0)
							}), times = times, times.range = range(times), noise = noise, 
					noise.min = noise.min, noise.span = noise.span))
	
	tmpargs <- lapply(x, function(x) unlist(lapply(x, function(x,imax) x[imax], imax = which.max(x$extrema))))
	peaks <- data.frame(do.call("rbind", tmpargs))
	peaks <- cbind(data.frame(branch = row.names(peaks)), peaks,data.frame(iendtime = attr(x, "iendtime")))
	
	peak.snr <- round(peaks[["extrema"]]/noise.levels)
	peak.scale <- peaks[["scale"]]
	
	branch.lengths <- unlist(lapply(x, function(x, scale.range) length(which(x$scale >= 
												scale.range[1] & x$scale <= scale.range[2])), scale.range = scale.range))
	
	#branch.lengths1 <- unlist(lapply(x, function(x, scale.range) length(which(x$scale >= scale[1] & x$scale <= scale[length(scale)])), scale.range = scale.range))
	
	#good.snr <- peak.snr >= snr.min
#	good.scale <- peak.scale >= scale.range[1]      
##	good.scale <- ((peak.scale >= max(scale.range[1],2))& (peak.scale <= min(scale.range[2], 25)))   #Modifed by Wenchao Zhang, Use lower scale boundary and Upper scale boundary to filter  
#	#good.length <- branch.lengths >= length.min
#	iendtime.min <- max(as.integer(noise.span/sampling.interval/4), 3)
#	iendtime.max <- length(times) - iendtime.min + 1
#	good.end <- peaks[["iendtime"]] > iendtime.min & peaks[["iendtime"]] < iendtime.max
	
	# Added some codes to filter out some very very small peaks. Added by wenchao zhang
	#peak_extrema <- peaks[["extrema"]]
	
	#good.peak_extrema <- peak_extrema> max(peak_extrema)/10.0
	
	#peaks <- peaks[which(good.snr & good.scale & good.length & good.end & good.peak_extrema), ]	
#	peaks <- peaks[which(good.scale  & good.end), ]	
#	snr <- peak.snr[which(good.scale & good.end)]
	snr <- peak.snr
	
	#Debug for none good peak cases, 
	if(nrow(peaks)>0) row.names(peaks) <- as.character(seq(nrow(peaks)))
	z <- list(x = times[peaks$iendtime], y = series[peaks$iendtime])
	attr(z, "peaks") <- peaks
	#attr(z, "snr.min") <- snr.min
	attr(z, "scale.range") <- scale.range
	attr(z, "length.min") <- length.min
	attr(z, "noise.span") <- noise.span
	attr(z, "noise.fun") <- noise.fun
	attr(z, "noise.min") <- noise.min
	attr(z, "snr") <- snr
	attr(z, "branch.lengths") <- branch.lengths
	z
}

parTICPeakpicking<-function(params,denoised=T)
{
	time1<-Sys.time()
	
	DataFilelist <- params$DataFiles
	nNode<-as.integer(params$nNode)
	clustType<-params$clustType
	
	library("snow")
	params$nNode <- min(as.integer(params$nNode),length(DataFilelist))
	cl<-makeCluster(params$nNode,type=clustType)
	
	clusterExport(cl,"peaks")#assign the current profiles to all nodes
	clusterExport(cl,"getPeaks")#assign the current profiles to all nodes
	clusterExport(cl,"parseFileName")#assign the current profiles to all nodes
	clusterExport(cl,"MSWtoCWT")
	clusterExport(cl,"cwt")
	clusterExport(cl,"extendLength")
	clusterExport(cl,"extendNBase")
	clusterExport(cl,"wavCWTTree_ModifyByYan")
	clusterExport(cl,"wavCWTTree_ModifyByYan_Min")
	clusterExport(cl,"wavCWTPeaks_Modify")
#	clusterExport(cl,"WTMM")
#	clusterExport(cl,"wtmmBranches")
	
	##calculate the distance between each pair of profiles
	clusterApply(cl,DataFilelist,TICPeakpicking,params,denoised)#parallel version
	stopCluster(cl)
	time2<-Sys.time()
	time2-time1
}

TICPeakpicking_ncdf4<-function(inFilePath,params,denoised=T)
{
	
	library("ncdf4")
	WorkDir <- params$WorkDir
	
	###############
	#read TIC data
	################
	fileName<-parseFileName(inFilePath)
	denoisedTICFile<-paste(WorkDir,"/output/TIC/denoised_",fileName,"_TIC.cdf",sep="")
	
	if(file.exists(inFilePath)==FALSE)next
	if(denoised)
	{
		cat(fileName,"reading denoised TIC data...\n")
		ncid <- nc_open(denoisedTICFile)
	}else
	{
		cat(fileName,"reading raw TIC data...\n")
		ncid <- nc_open(inFilePath)
	}
	
	vecInt <- ncvar_get(ncid, varid="total_intensity")
	nc_close(ncid)
	remove(ncid)	
	
	###############
	#peak picking
	###############
	cat(fileName,"TIC peak picking...")
	PeakList<-getPeaks(vecInt,params,mz=NULL,mzVec=NULL)
	## save peak picking results to cdf
	write.csv(PeakList,file=paste(WorkDir,"/output/peakpicking/",fileName,"_TIC_PeakList.csv",sep=""),row.names = F)
}

TICPeakpicking<-function(inFilePath,params,denoised=T)
{
	
	library("ncdf")
	WorkDir <- params$WorkDir
	
	###############
	#read TIC data
	################
	fileName<-parseFileName(inFilePath)
	denoisedTICFile<-paste(WorkDir,"/output/TIC/denoised_",fileName,"_TIC.cdf",sep="")
	
	if(file.exists(inFilePath)==FALSE)next
	if(denoised)
	{
		cat(fileName,"reading denoised TIC data...\n")
		ncid <- open.ncdf(denoisedTICFile)
	}else
	{
		cat(fileName,"reading raw TIC data...\n")
		ncid <- open.ncdf(inFilePath)
	}
	
	vecInt <- get.var.ncdf(ncid, varid="total_intensity")
	close.ncdf(ncid)
	remove(ncid)	
	
	###############
	#peak picking
	###############
	cat(fileName,"TIC peak picking...")
	PeakList<-getPeaks(vecInt,params,mz=NULL,mzVec=NULL)
	## save peak picking results to cdf
	write.csv(PeakList,file=paste(WorkDir,"/output/peakpicking/",fileName,"_TIC_PeakList.csv",sep=""),row.names = F)
}

parEICpeakpicking<-function(params)
{
	time1<-Sys.time()
	###############
	#get parameters
	###############
	#library(tcltk)
	library("snow")
	library("ncdf")
	#library("wmtsa")
	#library("MassSpecWavelet")
	
	WorkDir <- params$WorkDir
	DataFilelist <- params$DataFiles
	
	nNode<-as.integer(params$nNode)
	clustType<-params$clustType
	cl<-makeCluster(nNode,type=clustType)
	
	###############
	#peak picking
	###############
	
	for(fileindex in 1:length(DataFilelist))
	{
		
		###############
		#read EIC data
		################
		inFilePath <- DataFilelist[[fileindex]]
		fileName<-parseFileName(inFilePath)
		
		cat(fileName,"reading EIC data...\n")
		#EICFile<-paste(WorkDir,"/output/EIC/",fileName,"EIC.cdf",sep="")
		
		#check if it exist denoised EIC file
		rawEICFile <-paste(WorkDir,"output/EIC/",fileName,"EIC.cdf",sep="")
		denoiseEICFile<-paste(WorkDir,"output/EIC/denoised_",fileName,"EIC.cdf",sep="")
		EICFile<-ifelse(file.exists(denoiseEICFile),denoiseEICFile,rawEICFile)
		
		if(file.exists(inFilePath)==FALSE)next
		ncid <- open.ncdf(EICFile)
		vecInt <- get.var.ncdf(ncid, varid="intVec")
		mzVec<-get.var.ncdf(ncid, varid="mzVec")
		
		#mzVec <- c(40:400) # for testing
		
		close.ncdf(ncid)
		remove(ncid)	
		
		###############
		#peak picking
		################
		cat(fileName,"peak picking...")
		#######peak picking
		mzGroups<-clusterSplit(cl,mzVec)
		
		clusterExport(cl,"peaks")#assign the current profiles to all nodes
		clusterExport(cl,"getPeaks")#assign the current profiles to all nodes
		clusterExport(cl,"GaussianDisimilarity")
		clusterExport(cl,"Sharpness")#assign the current profiles to all nodes
		clusterExport(cl,"Sharpness1")#assign the current profiles to all nodes
		clusterExport(cl,"dotProduct")#assign the current profiles to all nodes
		clusterExport(cl,"parseFileName")#assign the current profiles to all nodes
		clusterExport(cl,"MSWtoCWT")
		clusterExport(cl,"cwt")
		clusterExport(cl,"extendLength")
		clusterExport(cl,"extendNBase")
		clusterExport(cl,"wavCWTTree_ModifyByYan")
		clusterExport(cl,"wavCWTTree_ModifyByYan_Min")
		clusterExport(cl,"wavCWTPeaks_Modify")
		
		time1<-Sys.time()
		PeakList<-clusterApply(cl,mzGroups,getPeaksGroup,vecInt,mzVec,params,fileName)#parallel version
		time2<-Sys.time()
		time2-time1
		PeakList<-unlist(PeakList,recursive=F)
		PeakList<-do.call(rbind,PeakList)
		
		write.csv(PeakList,file=paste(WorkDir,"/output/peakpicking/",fileName,"_EIC_PeakList.csv",sep=""),row.names=F)
		
	}
	stopCluster(cl)
	
	time2<-Sys.time()
	time2-time1
}


EICpeakpicking<-function(params)
{
	time1<-Sys.time()
	###############
	#get parameters
	###############

	library("ncdf")
	WorkDir <- params$WorkDir
	DataFilelist <- params$DataFiles
	
	for(fileindex in 1:length(DataFilelist))
	{
		
		###############
		#read EIC data
		################
		inFilePath <- DataFilelist[[fileindex]]
		fileName<-parseFileName(inFilePath)

		#check if it exist denoised EIC file
		rawEICFile <-paste(WorkDir,"output/EIC/",fileName,"EIC.cdf",sep="")
		denoiseEICFile<-paste(WorkDir,"output/EIC/denoised_",fileName,"EIC.cdf",sep="")
		EICFile<-ifelse(file.exists(denoiseEICFile),denoiseEICFile,rawEICFile)
		
		if(file.exists(inFilePath)==FALSE)next
		ncid <- open.ncdf(EICFile)
		vecInt <- get.var.ncdf(ncid, varid="intVec")
		mzVec<-get.var.ncdf(ncid, varid="mzVec")
		
		close.ncdf(ncid)
		remove(ncid)	
		
		cat(fileName,"start peak picking...\n")
		mzGroup <- mzVec
		time1<-Sys.time()
		PeakList <- getPeaksGroup(mzGroup,vecInt,mzVec,params,fileName)
		PeakList<-do.call(rbind,PeakList)
		write.csv(PeakList,file=paste(WorkDir,"/output/peakpicking/",fileName,"_EIC_PeakList.csv",sep=""),row.names=F)	
		time2<-Sys.time()
		time2-time1
		cat(fileName,"finish peak picking!\n")
	}
}

EICpeakpicking_ncdf4<-function(params)
{
	time1<-Sys.time()
	###############
	#get parameters
	###############
	
	library("ncdf4")
	WorkDir <- params$WorkDir
	DataFilelist <- params$DataFiles
	
	for(fileindex in 1:length(DataFilelist))
	{
		
		###############
		#read EIC data
		################
		inFilePath <- DataFilelist[[fileindex]]
		fileName<-parseFileName(inFilePath)
		
		#check if it exist denoised EIC file
		rawEICFile <-paste(WorkDir,"output/EIC/",fileName,"EIC.cdf",sep="")
		denoiseEICFile<-paste(WorkDir,"output/EIC/denoised_",fileName,"EIC.cdf",sep="")
		EICFile<-ifelse(file.exists(denoiseEICFile),denoiseEICFile,rawEICFile)
		
		if(file.exists(inFilePath)==FALSE)next
		ncid <- nc_open(EICFile)
		vecInt <- ncvar_get(ncid, varid="intVec")
		mzVec<-ncvar_get(ncid, varid="mzVec")
		
		nc_close(ncid)
		remove(ncid)
		
		cat(fileName,"start peak picking...\n")
		mzGroup <- mzVec
		time1<-Sys.time()
		PeakList <- getPeaksGroup(mzGroup,vecInt,mzVec,params,fileName)
		PeakList<-do.call(rbind,PeakList)
		write.csv(PeakList,file=paste(WorkDir,"/output/peakpicking/",fileName,"_EIC_PeakList.csv",sep=""),row.names=F)	
		time2<-Sys.time()
		time2-time1
		cat(fileName,"finish peak picking!\n")
	}
}

getPeaksGroup<-function(mzGroup,vecInt,mzVec,params,fileName)
{
	nMz<-length(mzGroup)
	groupResult<-vector("list",nMz)
	sink(paste(params$WorkDir,"/output/peakpicking/",fileName,".txt",sep=""),append=T)	
	
	#time1<-Sys.time()
	for(i in 1:nMz)
	{
		mz<-mzGroup[i]	
		print(mz)
		if (mz==0) next
		groupResult[[i]]<-getPeaks(vecInt,params,mz,mzVec)
	}
	#time2<-Sys.time()
	#time2-time1
	sink()
	groupResult
}

getPeaks<-function(vecInt,params,mz=NULL,mzVec=NULL)
{
	library("wmtsa")
	
	if(!is.null(mz))
	{
		#EIC, int vector from the mz position
		totalscan<-length(vecInt)/length(mzVec)
		mzInd<-which(mzVec==mz)
		startInd<-(mzInd-1)*totalscan+1
		endInd<-startInd+totalscan-1
		curVecInt <- vecInt[startInd:endInd]
#		BHR<-as.numeric(params$BHR_EIC)#0.3##boundary/height raio
#		edgeHightDiffRatio<-as.numeric(params$EHR_EIC)#0.2
		BHR<-0.3
		edgeHightDiffRatio<-0.2
		offset1<-startInd-1
		maxWindowlength<-350
		
	}else
	{
		#TIC, int vector directly from parameter
		curVecInt<-vecInt
		totalscan<-length(curVecInt)
		BHR<-0.4
		edgeHightDiffRatio<-0.2
		maxWindowlength<-450
	}
	
	
	WorkDir<-params$WorkDir
	nPoints<-length(curVecInt)
	delaytime<-params$delaytime
	ScanInterval<-params$ScanInterval

	# apex detection, using modified wavelet algorithms
	MaxPeakWidth <- 200 # scans, will be provided by ADAP user
	scalerange <- round(c(2,MaxPeakWidth)/2)
	scales <-  seq(from=scalerange[1], to=scalerange[2], by=2)

	CWTCoeff <- MSWtoCWT(curVecInt,scales)
	# identify peak apex using modified funcs from wtmsa package
	curPeakList <- wavCWTTree_ModifyByYan(CWTCoeff,type = "maxima") 
	#plot(curVecInt[curPeakList$pkInd])
	
	# local maximun calculation
	peakSpan<-as.integer(params$Peak_span)
	isPeak <- peaks(x=curVecInt, span=peakSpan)
	peakInd<-which(isPeak==TRUE)
	peakInd <- peakInd[which(curVecInt[peakInd]>=100)]

	valleySpan<-as.integer(params$Valley_span)
	isMin <- peaks(-curVecInt, valleySpan)
	NonZeros<-curVecInt>0
	isValley <- isMin&NonZeros
	valleyInd_LM<-which(isValley==TRUE)

	# filtering by StN, ignoring noisy peaks to avoid them merged to neighbouring peaks
	StN_Th <- 20
	curPeakList <- curPeakList[which(curVecInt[curPeakList$pkInd]>=200&curPeakList$StN>=StN_Th),]
	
	# identify boundaries, more accurate than local min but could miss severals
	tree_min <- wavCWTTree_ModifyByYan_Min(CWTCoeff,type = "minima") 
	valleyInd <- sort(wavCWTPeaks_Modify(tree_min)$x)


	#curPeakList <- subset(curPeakList,StN>=StN_Th)
	
	if (nrow(curPeakList) > 0) {
		# locate neighbouring peak boundaries initially
		curPeakList$Lbound <- 0
		curPeakList$Rbound <- 0
		curPeakList$Intensity <- curVecInt[curPeakList$pkInd]
		
		for (i in 1:nrow(curPeakList)) {
			curPeak <- curPeakList[i,]$pkInd 
			curBoundary_range <- findInterval(curPeak,valleyInd)
			if (curBoundary_range>0) {
				curPeakList[i,]$Lbound  <- valleyInd[curBoundary_range]
				
				if (curBoundary_range==length(valleyInd)) {
					curPeakList[i,]$Rbound <- nPoints
				}else {
					curPeakList[i,]$Rbound <- valleyInd[curBoundary_range+1]
				}
			}
		}
		
		curPeakList <- subset(curPeakList,Lbound>0&Rbound>0&Intensity>0)
		# exact boudnary detection and deviding windows for deconvolution 
		if (length(valleyInd) > nrow(curPeakList) & nrow(curPeakList) > 0 ) {
			#maxWindowlength<-as.integer(params$MaxWindow_length)
			
			curPeakList$lboundInd<-0
			curPeakList$rboundInd<-0
			curPeakList$isApex<-0 # flag for merging
			curPeakList$isShared<-0
			WaitingToMerge <- FALSE

			for(i in 1:nrow(curPeakList))
			#for(i in 1:8)
			{
				if(i==nrow(curPeakList)) # the last peak
				{
					#print(i)
					curPeakLbound <- curPeakList$Lbound[i]
					curPeakRbound <- curPeakList$Rbound[i]	
					curPeakHight<-curPeakList$Intensity[i]
					curLboundInt <- curVecInt[curPeakLbound]
					curRboundInt <- curVecInt[curPeakRbound]
					
					LeftToApex_Ratio <- curLboundInt/curPeakHight
					BoundToApex_Ratio <- abs(curLboundInt-curRboundInt)/curPeakHight
					
#					if (!is.null(mz)) {
#						decisition1 <- (LeftToApex_Ratio < BHR &BoundToApex_Ratio < edgeHightDiffRatio)|nrow(curPeakList)==1
#					}else{
#						decisition1 <- (BoundToApex_Ratio < edgeHightDiffRatio)|nrow(curPeakList)==1
#					}
#					
#					if (decisition1) {
					if ((LeftToApex_Ratio < BHR &BoundToApex_Ratio < edgeHightDiffRatio)|nrow(curPeakList)==1) { # include single peak exists
						
						if (WaitingToMerge == FALSE) {
							curPeakList[i,]$isApex=1
							curPeakList[i,]$lboundInd=curPeakLbound 
							curPeakList[i,]$rboundInd=curPeakRbound
							curPeakList[i,]$isShared<-0
						}else {
							localPkList <- which(curPeakList$isApex==-1)
							previousLbound<-curPeakList[localPkList,]$Lbound[1]
							
							if((curPeakRbound-previousLbound)<=maxWindowlength) {
								curPeakList[c(i,localPkList),]$isApex=1
								curPeakList[c(i,localPkList),]$lboundInd=previousLbound 
								curPeakList[c(i,localPkList),]$rboundInd=curPeakRbound
								curPeakList[c(i,localPkList),]$isShared<-1
								WaitingToMerge <- FALSE
							}else {
								curPeakList[c(i,localPkList),]$isApex=1
								curPeakList[i,]$lboundInd=curPeakLbound 
								curPeakList[i,]$rboundInd=curPeakRbound
								curPeakList[i,]$isShared<-0
								
								# not merging and force to split
								curPeakList[localPkList,]$lboundInd=previousLbound
								curPeakList[localPkList,]$rboundInd=max(curPeakList[localPkList,]$Rbound)
								
								if (length(localPkList)>1) {
									curPeakList[localPkList,]$isShared<-1	
								}else {
									curPeakList[localPkList,]$isShared<-0	
								}
							}
						}
							
						
					}else { # merge to the previous one
						previousLbound<-curPeakList[i-1,]$lboundInd			
						if((curPeakRbound-previousLbound)<=maxWindowlength)
						{								
							curPeakList[c(i-1,i),]$isApex=1
							curPeakList[i,]$lboundInd=previousLbound
							curPeakList[c(i-1,i),]$rboundInd=curPeakRbound
							curPeakList[c(i-1,i),]$isShared<-1	
						}else {
							curPeakList[i,]$isApex=1
							curPeakList[i,]$lboundInd=curPeakLbound 
							curPeakList[i,]$rboundInd=curPeakRbound
							curPeakList[i,]$isShared<-0
						}				
					}			
				}else {
					#print(i)
					curPeakInd<-curPeakList$pkInd[i]
					curPeakLbound <- curPeakList$Lbound[i]
					curPeakRbound <- curPeakList$Rbound[i]	
					NextPeakInd <- curPeakList$pkInd[i+1]
					
					if (curPeakRbound > NextPeakInd) {
						# find missing boundary between two neighbouring apex locations
						EIC_part <- curVecInt[c(curPeakInd:NextPeakInd)]
						add_BoundInd <- curPeakInd+which.min(EIC_part)-1
						curPeakRbound <- curPeakList[i,]$Rbound <- add_BoundInd
						curPeakList[i+1,]$Lbound <- add_BoundInd	
					}
					
					
					curLboundInt <- curVecInt[curPeakLbound]
					curRboundInt <- curVecInt[curPeakRbound]
					curPeakHight<-curPeakList$Intensity[i]
					LeftToApex_Ratio <- curLboundInt/curPeakHight
					RightToApex_Ratio <- curRboundInt/curPeakHight
					BoundToApex_Ratio <- abs(curLboundInt-curRboundInt)/curPeakHight
					
					
#					if (!is.null(mz)) {
#						decisition2 <- LeftToApex_Ratio < BHR & RightToApex_Ratio < BHR & BoundToApex_Ratio < edgeHightDiffRatio
#					}else{
#						decisition2 <- BoundToApex_Ratio < edgeHightDiffRatio
#					}
						
					
#					if (decisition2) {
					if (LeftToApex_Ratio < BHR & RightToApex_Ratio < BHR & BoundToApex_Ratio < edgeHightDiffRatio) {
						
						Rbound_Candidates <- valleyInd[which(valleyInd<=NextPeakInd&valleyInd>=curPeakInd)]
						if (length(Rbound_Candidates)>1) {
							Rbound_Candidates_Int <- curVecInt[Rbound_Candidates]
							RbtoApex <- Rbound_Candidates_Int/curPeakHight
							BdifftoApex <- abs(curLboundInt-Rbound_Candidates_Int)/curPeakHight
							qualify_candidates <- which(RbtoApex < BHR & BdifftoApex < edgeHightDiffRatio)
							if (length(qualify_candidates)>0) {
								curPeakRbound <- Rbound_Candidates[qualify_candidates[length(qualify_candidates)]]							
							
							}
						}
							
						if (WaitingToMerge == FALSE) {
							curPeakList[i,]$isApex=1
							curPeakList[i,]$lboundInd=curPeakLbound 
							curPeakList[i,]$rboundInd=curPeakRbound
							curPeakList[i,]$isShared<-0
						}else {
							localPkList <- which(curPeakList$isApex==-1)
							previousLbound<-curPeakList[localPkList,]$Lbound[1]
							curPeakList[c(i,localPkList),]$isApex=1
							curPeakList[c(i,localPkList),]$lboundInd=previousLbound 
							curPeakList[c(i,localPkList),]$rboundInd=curPeakRbound
							curPeakList[c(i,localPkList),]$isShared<-1
							WaitingToMerge <- FALSE
						}
						
						#preRbound <- curPeakRbound
					} else { # merging
						if (curLboundInt >= curRboundInt) { # merge to the left
							
							if (i==1) { # the first peak, can't merge to left 
								curPeakList[i,]$isApex=1
#								curPeakList[i,]$lboundInd=curPeakLbound 
#								curPeakList[i,]$rboundInd=curPeakRbound
								curPeakList[i,]$isShared<-0
							} else {
								if (WaitingToMerge == FALSE)  {
									previousLbound<-curPeakList[i-1,]$lboundInd
									if((curPeakRbound-previousLbound)<=maxWindowlength)
									{					
										
										curPeakList[i,]$isApex=1
										curPeakList[i,]$lboundInd=previousLbound
										PeaksToMerge <- which(curPeakList$lboundInd==previousLbound)
										curPeakList[PeaksToMerge,]$rboundInd=curPeakRbound
										curPeakList[c(i-1,i),]$isShared<-1	
									}else {
										curPeakList[i,]$isApex=1
										curPeakList[i,]$lboundInd=curPeakLbound 
										curPeakList[i,]$rboundInd=curPeakRbound
										curPeakList[i,]$isShared<-0
									}
								} else {
									localPkList <- which(curPeakList$isApex==-1)
									previousLbound<-curPeakList[localPkList,]$Lbound[1]
									
									curPeakList[c(i,localPkList),]$isApex=1
									curPeakList[c(i,localPkList),]$lboundInd=previousLbound
									curPeakList[c(i,localPkList),]$rboundInd=curPeakRbound
									curPeakList[c(i,localPkList),]$isShared<-1
									
									
#									if((curPeakRbound-previousLbound)<=maxWindowlength)
#									{		
#										curPeakList[c(i,localPkList),]$isApex=1
#										curPeakList[c(i,localPkList),]$lboundInd=previousLbound
##										PeaksToMerge <- which(curPeakList$lboundInd==previousLbound)
##										curPeakList[PeaksToMerge,]$rboundInd=curPeakRbound
#										curPeakList[c(i,localPkList),]$rboundInd=curPeakRbound
#										curPeakList[c(i,localPkList),]$isShared<-1
#										
#									}else {
##										curPeakList[c(i,localPkList),]$isApex=1
##										curPeakList[i,]$lboundInd=curPeakLbound 
##										curPeakList[i,]$rboundInd=curPeakRbound
##										curPeakList[i,]$isShared<-0
##										
##										# not merging and force to split
##										curPeakList[localPkList,]$lboundInd=previousLbound
##										curPeakList[localPkList,]$rboundInd=max(curPeakList[localPkList,]$Rbound)
##										
##										if (length(localPkList)>1) {
##											curPeakList[localPkList,]$isShared<-1	
##										}else {
##											curPeakList[localPkList,]$isShared<-0	
##										}
#										
#										curPeakList[c(i,localPkList),]$isApex=1
#										curPeakList[c(i,localPkList),]$lboundInd=previousLbound
#										curPeakList[c(i,localPkList),]$rboundInd=curPeakRbound
#										curPeakList[c(i,localPkList),]$isShared<-1
#										
#									}
								}
								
							}
							WaitingToMerge <- FALSE
						}else { # merging to the right
							Rbound_Candidates <- valleyInd[which(valleyInd<=NextPeakInd&valleyInd>=curPeakInd)]
							
							if (length(Rbound_Candidates)>1) {
								Rbound_Candidates_Int <- curVecInt[Rbound_Candidates]
								RbtoApex <- Rbound_Candidates_Int/curPeakHight
								BdifftoApex <- abs(curLboundInt-Rbound_Candidates_Int)/curPeakHight
								qualify_candidates <- which(RbtoApex < BHR & BdifftoApex < edgeHightDiffRatio)
								if (length(qualify_candidates)>0) {
									curPeakRbound <- Rbound_Candidates[qualify_candidates[length(qualify_candidates)]]
									curPeakList[i,]$isApex=1
									curPeakList[i,]$lboundInd=curPeakLbound 
									curPeakList[i,]$rboundInd=curPeakRbound
									curPeakList[i,]$isShared<-0
									WaitingToMerge <- FALSE
									
								}else{
									if ((curPeakRbound-curPeakLbound)<=maxWindowlength) {
										curPeakList[i,]$isApex=-1
										WaitingToMerge <- TRUE
										
									}else {
										curPeakList[i,]$isApex=1
										curPeakList[i,]$lboundInd=curPeakLbound 
										curPeakList[i,]$rboundInd=curPeakRbound
										curPeakList[i,]$isShared<-0
										WaitingToMerge <- FALSE
										
									}
								}
								
							}else {
								if ((curPeakRbound-curPeakLbound)<=maxWindowlength) {
									curPeakList[i,]$isApex=-1
									WaitingToMerge <- TRUE
								
								}else {
									curPeakList[i,]$isApex=1
									curPeakList[i,]$lboundInd=curPeakLbound 
									curPeakList[i,]$rboundInd=curPeakRbound
									curPeakList[i,]$isShared<-0
									WaitingToMerge <- FALSE
								}
							}							
						}
					}	
				}
			} # for loop
			#cat(paste("finish peak apex and boundary picking\n",sep=""))	
			

			# peak evaluation
			if(!is.null(mz))
			{
				
				if (nrow(curPeakList)>0) { 
					#curPeakList$gss <- 90
					curPeakList$shrp<- 0
					curPeakList$area <- 0	
					
					for(i in 1:nrow(curPeakList))
					{
						curPk<-curPeakList[i,]
						startInd<-curPk$Lbound
						endInd<-curPk$Rbound
						
						# combine local maximum to check shared peak feature
						pkCandidates <- peakInd[peakInd>=startInd&peakInd<=endInd]
						
						if (length(pkCandidates)>0) {
							candidate_ratio <- curVecInt[pkCandidates]/curPk$Intensity
							pk_num <- length(which(candidate_ratio>=0.1))
							
							if (pk_num>1) {
								curPeakList[i,]$isShared=1
							}
						}
						
#						if (length(pkCandidates)>1) {
#							coeluting_peaks <- pkCandidates[which(curVecInt[pkCandidates]!=curPk$Intensity)]
#							if (length(coeluting_peaks)>0) {
#								coeluting_peaks_ind <- NULL
#								for (j in 1:length(coeluting_peaks)) {
#									if (min(abs(coeluting_peaks[j]-valleyInd_LM))>5) {
#										coeluting_peaks_ind <- c(coeluting_peaks_ind,j)
#									}
#									
#								}
#								if (length(coeluting_peaks_ind)>0) {
#									candidate_ratio <- curVecInt[coeluting_peaks[coeluting_peaks_ind]]/curPk$Intensity
#									pk_num <- length(which(candidate_ratio>=0.1))
#									
#									if (pk_num>0) {
#										curPeakList[i,]$isShared=1
#									}
#								}
#							}
#						}
						
						curProfile <- curVecInt[startInd:endInd]
						curPeakList[i,]$area<-sum(curProfile)	
						
						pkPos <- which.max(curProfile)
						base <- c(curProfile[1],curProfile[length(curProfile)])
						time <- c(1,length(curProfile))
						linear.model <- lm(base~time)
						coeff <- as.numeric(linear.model$coefficients)
						pkPredict <- coeff[2]*pkPos+coeff[1]
						P25Height <- round((curProfile[pkPos]-pkPredict)*0.25)
						scan_diff <- abs(c(startInd:endInd)-curPk$pkInd)
						
						LSide_pair <- subset(data.frame(Int=curProfile[1:(pkPos-1)],scanDist=sort(c(1:(pkPos-1)),decreasing=TRUE)),Int>=P25Height)
						RSide_pair <- subset(data.frame(Int=curProfile[(pkPos+1):length(curProfile)],scanDist=c(1:(length(curProfile)-pkPos))),Int>=P25Height)
						
						curPeakList[i,]$shrp <- round(mean(c(median((curProfile[pkPos]-LSide_pair$Int)/LSide_pair$scanDist),
								median((curProfile[pkPos]-RSide_pair$Int)/RSide_pair$scanDist))))
		
		
						#curPeakList[i,]$gss<-round(GaussianDisimilarity(curProfile),digits=2)
		
					}
					
					curPeakList$mz=mz
					curPeakList$offset<-offset1		
					curPeakList<-subset(curPeakList,select=c("mz","pkInd","Lbound","Rbound","lboundInd","rboundInd","StN","Intensity","area","shrp","isShared","offset"))
					curPeakList
				}
			}else
			{
				# TIC peak picking output
				subset(curPeakList,select=c("pkInd","Lbound","Rbound","lboundInd","rboundInd","isShared"))
			}
		} # if statement
	}
}


copyShape<-function(p,apex,from) {
	xx<-p$ET
	yy<-p$int
	maxInd<-which(xx==apex)
	endInd<-length(yy)	
	if(from==1) 
	{
		yy<-yy[c(1:(maxInd-1),maxInd:1)]
		xx<-xx[1]:(xx[1]+length(yy)-1)
		
	}else
	{
		yy<-yy[c(endInd:(maxInd+1),maxInd:endInd)]	
		
		# in case the first several peaks starting almost at the 1st scan
		# if the mirror image comes from the right part, the resulting profile
		# would have negative left boundary
		if ((xx[endInd]-length(yy)+1)<0) {
			yy<-yy[c(1:(maxInd-1),maxInd:1)]
			xx<-xx[1]:(xx[1]+length(yy)-1)
		}else {
			xx<-(xx[endInd]-length(yy)+1):xx[endInd]
		}	
	}	
	data.frame(ET=xx,int=yy)
}


## gaussian curve fitting
GaussianDisimilarity<-function(yy) {
	res<-90
	xx<-1:length(yy)
	tryCatch({fit1<-nls(yy ~ (w/sqrt(2*pi*sigma1^2)) * exp(-(xx - mu)^2/(2 * sigma1^2)),control = list(maxiter = 50),  
						start = list(w = max(yy)*sqrt(2*pi*1^2), sigma1=1,mu = xx[which.max(yy)]),trace = F)
				res<-dotProduct(fitted(fit1),yy)},
			error=function(er)
			{
				
				#cat(mz,er);
			}
	)
	#ifelse(is.na(res),90,res
	res
}

## Produdce two mirror image: left/right half
## The cutting point is apex
GaussianDisimilarity1<-function(yy,apex) {
	
	maxInd<-apex
	endInd<-length(yy)
	yy1<-yy[c(1:(maxInd-1),maxInd:1)]
	yy2<-yy[c(endInd:(maxInd+1),maxInd:endInd)]
	
	dist1<-GaussianDisimilarity(yy1)
	dist2<-GaussianDisimilarity(yy2)
	c(dist1,dist2)
	
}



peaks<-function (x, span = 3) 
{
	z <- embed(as.vector(x), span)
	s <- span%/%2
	
	if (span%%2==1) {
		result <- max.col(z,ties.method="first") == 1 + s
		c(rep(FALSE, s), result, rep(FALSE, s))
	}else if (span%%2==0) {
		result <- max.col(z,ties.method="first") == s
		c(rep(FALSE, s), result, rep(FALSE, s-1))
	}
}

Sharpness<-function(yy) {
	res<-0
	localPeakInd<-which.max(yy)
	
	nLength<-length(yy)
	shrpVec<-c((yy[2:localPeakInd]-yy[1:(localPeakInd-1)])/yy[1:(localPeakInd-1)],(yy[localPeakInd:(nLength-1)]-yy[(localPeakInd+1):nLength])/yy[(localPeakInd+1):nLength])
	res<-sum(shrpVec[is.finite(shrpVec)])
	res
}

#Sharpness<-function(curProfile,rt,curPk) {
#	res<-0
#	localPeakInd<-which(rt)
#	
#	nLength<-length(yy)
#	shrpVec<-c((yy[2:localPeakInd]-yy[1:(localPeakInd-1)])/yy[1:(localPeakInd-1)],(yy[localPeakInd:(nLength-1)]-yy[(localPeakInd+1):nLength])/yy[(localPeakInd+1):nLength])
#	res<-sum(shrpVec[is.finite(shrpVec)])
#	res
#}
#

Sharpness1<-function(yy,apex) {
	maxInd<-apex
	endInd<-length(yy)
	yy1<-yy[c(1:(maxInd-1),maxInd:1)]
	yy2<-yy[c(endInd:(maxInd+1),maxInd:endInd)]
	
	shrp1<-Sharpness(yy1)
	shrp2<-Sharpness(yy2)
	c(shrp1,shrp2)
}

