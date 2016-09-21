
d_wavCWT <- function(x, scale.range=deltat(x) * c(1, length(x)), n.scale=100,
                     wavelet="gaussian2", shift=5, variance=1){
    
    # check input argument types and lengths
    checkVectorType(scale.range,"numeric")
    checkScalarType(n.scale,"integer")
    checkScalarType(wavelet,"character")
    checkScalarType(shift,"numeric")
    checkScalarType(variance,"numeric")
    checkRange(n.scale, c(1,Inf))
    
    # obtain series name
    series.name <- deparse(substitute(x))
    
    if (length(scale.range) != 2)
        stop("scale.range must be a two-element numeric vector")
    if (variance <= 0)
        stop("variance input must be positive")
    
    # obtain sampling interval
    sampling.interval <- deltat(x)
    
    # create a vector of log2-spaced scales over the specified range of scale
    octave <- logb(scale.range, 2)
    scale  <- ifelse1(n.scale > 1, 2^c(octave[1] + seq(0, n.scale - 2) * diff(octave) /
                                           (floor(n.scale) - 1), octave[2]), scale.range[1])
    
    # project the scale vector onto a uniform grid of sampling interval width
    scale <- unique(round(scale / sampling.interval) * sampling.interval)
    n.scale <- length(scale)
    
    # check scale range
    if (abs(min(scale) - sampling.interval) > .Machine$double.eps)
        stop("Minimum scale must be greater than or equal to sampling interval ",
             "of the time series")
    
    # map the time vector
    if (inherits(x, "signalSeries"))
    {
        times <- as(x@positions,"numeric")
        x <- x@data
    }
    else
    {
        
        times <- time(x)
        x <- as.vector(x)
    }
    
    # ensure that x is a double vector
    storage.mode(x) <- "double"
    
    # map the wavelet filter and corresponding argument
    gauss1 <- c("gaussian1", "gauss1")
    gauss2 <- c("gaussian2", "gauss2", "mexican hat", "sombrero")
    supported.wavelets <- c("haar", gauss1, gauss2, "morlet")
    wavelet <- match.arg(lowerCase(wavelet), supported.wavelets)
    
    # map the filter type to MUTILS type
    # 4: gaussian1
    # 5: gaussian2, sombrero, mexican hat
    # 6: morlet
    # 7: haar
    filter <- mutilsFilterTypeContinuous(wavelet)
    
    if (filter == 4)
    {
        filter.arg <- sqrt(variance)
        wavelet    <- "gaussian1"
    }
    else if (filter == 5)
    {
        filter.arg <- sqrt(variance)
        wavelet    <- "gaussian2"
    }
    else if (filter == 6)
    {
        filter.arg <- shift
        wavelet    <- "morlet"
    }
    else if (filter == 7)
    {
        filter.arg <- 0.0
        wavelet    <- "haar"
        
        # in the case of the Haar, the Euler-Macluarin approximation
        # to the DFT of the wavelet filter is not defined for non-integer
        # multiples of the sampling interval. therefore, we coerce the
        # scales accordingly in this case
        scale <- sampling.interval * unique(round(scale / sampling.interval))
    }
    else
    {
        stop("Unsupported filter type")
    }
    
    # calculate the CWT
    z <- itCall("RS_wavelets_transform_continuous_wavelet",
                as.numeric(x),
                as.numeric(sampling.interval),
                as.integer(filter),
                as.numeric(filter.arg),
                as.numeric(scale))
    #
    #COPY=rep(FALSE, 5),
    #CLASSES=c("numeric", "numeric", "integer", "numeric", "numeric"),
    #PACKAGE="ifultools")
    
    # if the impulse response of the waveleyt filter is purely real
    # then transform CWT coefficients to purely real as well
    if (wavelet != "morlet")
        z <- Re(z)
    
    # assign attributes
    attr(z, "scale")       <- scale
    attr(z, "time")        <- as.vector(times)
    attr(z, "wavelet")     <- wavelet
    attr(z, "series")      <- x
    attr(z, "sampling.interval") <- sampling.interval
    attr(z, "series.name") <- series.name
    attr(z, "n.sample")    <- length(x)
    attr(z, "n.scale")     <- n.scale
    attr(z, "filter.arg")  <- filter.arg
    
    oldClass(z) <- "wavCWT"
    
    z
}