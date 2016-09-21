# Modified by Xiuxia Du in September 2015.





EICpeakpicking<-function(params)
{
    
    workDir <- params$workDir
    eicFiles <- list.files(path=paste(workDir, "EIC", sep=.Platform$file.sep))
    
    
    
    
    peak_picking_dir <- paste(getwd(), .Platform$file.sep, "output", .Platform$file.sep, "peakpicking", sep="")
    if (!dir.exists(peak_picking_dir)) {
        dir.create(path=peak_picking_dir)
    }
    
    
    
    

    for (fileIndex in 1:length(eicFiles))
    {
        
      # 1. Import one EIC file
        
#         inFilePath <- DataFilelist[[fileindex]]
#         fileName<-parseFileName(inFilePath)
#         
#         #check if it exist denoised EIC file
#         rawEICFile <-paste(WorkDir,"output/EIC/",fileName,"EIC.cdf",sep="")
#         denoiseEICFile<-paste(WorkDir,"output/EIC/denoised_",fileName,"EIC.cdf",sep="")
#         EICFile<-ifelse(file.exists(denoiseEICFile),denoiseEICFile,rawEICFile)
#         
#         if(file.exists(inFilePath)==FALSE)next
        
        eicFileName <- paste(workDir, "EIC", eicFiles[fileIndex], sep=.Platform$file.sep)
        # ncid <- open.ncdf(eicFileName)
        # vecInt <- get.var.ncdf(ncid, varid="intVec")
        # mzVec <- get.var.ncdf(ncid, varid="mzVec")
        # close.ncdf(ncid)
        
        ncid <- nc_open(eicFileName)
        mzVec <- ncvar_get(ncid, varid="mzVec")
        intVec <- ncvar_get(ncid, varid="intVec")
        scan_acquisition_time <- ncvar_get(ncid, varid="scan_acquisition_time")
        nc_close(ncid)
        
        
        remove(ncid)	
        

        
        
        
        cat(paste(eicFiles[fileIndex], " peak picking starts ...\n", sep=""))
        time_start <- Sys.time()
        
        
        
        
        
        mzGroup <- mzVec
        # peakList <- getPeaksGroup(mzGroup, vecInt, mzVec, params, eicFiles[fileIndex])
        peakList <- getPeaksGroup(mzGroup, intVec, mzVec, params)
        peakList <- do.call(rbind, peakList)
        
        
        
        
        
        ind <- regexpr(pattern=".cdf", text=eicFiles[fileIndex])
        file_identifier <- substr(x=eicFiles[fileIndex], start=1, stop=ind-1)
        out_file_name <- paste(file_identifier, "_peaklist.csv", sep="") 
        out_file_full_name <- paste(workDir, "output", .Platform$file.sep, "peakpicking", .Platform$file.sep, 
                                    out_file_name, sep="")
        
        if (file.exists(out_file_full_name)) {
            file.remove(out_file_full_name)
        }
        write.csv(peakList,
                  file=out_file_full_name,
                  row.names=F)	
        
        
        
        
        
        time_end <- Sys.time()
        cat(paste(eicFiles[fileIndex], "peak picking finished! Time spent peak picking = ", 
                  round((time_end-time_start)/60, digits=2), " minutes, \n", sep=""))
    }
}





getPeaksGroup<-function(mzGroup, vecInt, mzVec, params)
{
    
#     mz_ind_start <- 1
#     mz_ind_end <- length(mzGroup)
    
    mz_ind_start <- 34
    mz_ind_end <- 35
    groupResult <- vector(mode="list", length=mz_ind_end-mz_ind_start+1)
    
    
    
    
    
    for (i in mz_ind_start:mz_ind_end)
    # for(i in 1:length(mzGroup))
    {
        cur_mz <- mzGroup[i]	
        print(cur_mz)
        
        
        
        
        
        if (cur_mz==0) {
            next
        } else if(!is.null(cur_mz)) {
            
            #EIC, int vector from the mz position
            totalscan <- length(vecInt)/length(mzVec)
            mzInd <- which(mzVec==cur_mz)
            startInd <- (mzInd-1)*totalscan + 1
            endInd <- startInd + totalscan - 1
            curVecInt <- vecInt[startInd:endInd]
            #		BHR<-as.numeric(params$BHR_EIC)#0.3##boundary/height raio
            #		edgeHightDiffRatio<-as.numeric(params$EHR_EIC)#0.2

            offset1 <- startInd-1
            
            groupResult[[i]] <- getPeaks(vecInt=curVecInt, mz=cur_mz, params, offset=offset1)
            
            
        } else {
            #TIC, int vector directly from parameter
            curVecInt <- vecInt
            totalscan <- length(curVecInt)
            BHR <- params$BHR_TIC
            edgeHightDiffRatio <- params$edgeHightDiffRatio_TIC
            maxWindowlength <- params$maxWindowlength_TIC
            
            groupResult[[i]] <- getPeaks(vecInt=curVecInt, params)
        }
    }
    
    return(groupResult)
}






getPeaks <- function(vecInt, mz=NULL, params, offset=NULL)
{

    MaxPeakWidth <- params$MaxPeakWidth # scans, will be provided by ADAP user
    StN_Th <- params$StN_Th
    
    BHR <- params$BHR_EIC
    edgeHightDiffRatio <- params$edgeHightDiffRatio_EIC
    maxWindowlength <- params$maxWindowlength_EIC
    
    
    
    
    
    nPoints <- length(vecInt) 
    # delaytime<-params$delaytime # not found used
    # ScanInterval<-params$ScanInterval # not found used
    
    
    
    
    
    # apex detection, using modified wavelet algorithms
    scalerange <- round(c(2,MaxPeakWidth)/2)
    scales <-  seq(from=scalerange[1], to=scalerange[2], by=2)
    
    CWTCoeff <- MSWtoCWT(vecInt, scales)
    # identify peak apex using modified funcs from wtmsa package
    curPeakList <- wavCWTTree_ModifyByYan(CWTCoeff, type="maxima") 
    #plot(curVecInt[curPeakList$pkInd])
    
    
    
    
    
    # plotting for debugging
#     xx <- vecInt[1:10000]
#     intensity_threshold <- 2.5e+4
#     II <- which(xx > intensity_threshold)
#     JJ <- intersect(curPeakList$pkInd, II)
#     plot(xx, type="l")
#     points(JJ, rep(intensity_threshold, times=length(JJ)), pch=16, cex=1, col="red")
    
    
    
    
    
    
    
    # local maximun calculation
    peakSpan <- as.integer(params$peak_span)
    isPeak <- peaks(x=vecInt, span=peakSpan)
    
    
    
    
    
    peakInd <- which(isPeak==TRUE)
    peakInd <- peakInd[which(vecInt[peakInd]>=100)]
    
    valleySpan <- as.integer(params$valley_span)
    isMin <- peaks(-vecInt, valleySpan)
    NonZeros <- vecInt>0
    isValley <- isMin&NonZeros
    valleyInd_LM <- which(isValley==TRUE)
    
    # filtering by StN, ignoring noisy peaks to avoid them merged to neighbouring peaks
    
    curPeakList <- curPeakList[which(vecInt[curPeakList$pkInd]>=200&curPeakList$StN>=StN_Th),]
    
    # identify boundaries, more accurate than local min but could miss severals
    tree_min <- wavCWTTree_ModifyByYan_Min(CWTCoeff,type = "minima") 
    valleyInd <- sort(wavCWTPeaks_Modify(tree_min)$x)
    
    
    #curPeakList <- subset(curPeakList,StN>=StN_Th)
    
    if (nrow(curPeakList) > 0) {
        # locate neighbouring peak boundaries initially
        curPeakList$Lbound <- 0
        curPeakList$Rbound <- 0
        curPeakList$Intensity <- vecInt[curPeakList$pkInd]
        
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
                    curLboundInt <- vecInt[curPeakLbound]
                    curRboundInt <- vecInt[curPeakRbound]
                    
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
                        EIC_part <- vecInt[c(curPeakInd:NextPeakInd)]
                        add_BoundInd <- curPeakInd+which.min(EIC_part)-1
                        curPeakRbound <- curPeakList[i,]$Rbound <- add_BoundInd
                        curPeakList[i+1,]$Lbound <- add_BoundInd	
                    }
                    
                    
                    curLboundInt <- vecInt[curPeakLbound]
                    curRboundInt <- vecInt[curPeakRbound]
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
                            Rbound_Candidates_Int <- vecInt[Rbound_Candidates]
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
                                Rbound_Candidates_Int <- vecInt[Rbound_Candidates]
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
                            candidate_ratio <- vecInt[pkCandidates]/curPk$Intensity
                            pk_num <- length(which(candidate_ratio>=0.1))
                            
                            if (pk_num>1) {
                                curPeakList[i,]$isShared=1
                            }
                        }
                        
                        #						if (length(pkCandidates)>1) {
                        #							coeluting_peaks <- pkCandidates[which(vecInt[pkCandidates]!=curPk$Intensity)]
                        #							if (length(coeluting_peaks)>0) {
                        #								coeluting_peaks_ind <- NULL
                        #								for (j in 1:length(coeluting_peaks)) {
                        #									if (min(abs(coeluting_peaks[j]-valleyInd_LM))>5) {
                        #										coeluting_peaks_ind <- c(coeluting_peaks_ind,j)
                        #									}
                        #									
                        #								}
                        #								if (length(coeluting_peaks_ind)>0) {
                        #									candidate_ratio <- vecInt[coeluting_peaks[coeluting_peaks_ind]]/curPk$Intensity
                        #									pk_num <- length(which(candidate_ratio>=0.1))
                        #									
                        #									if (pk_num>0) {
                        #										curPeakList[i,]$isShared=1
                        #									}
                        #								}
                        #							}
                        #						}
                        
                        curProfile <- vecInt[startInd:endInd]
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
                    curPeakList$offset <- offset		
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







MSWtoCWT <- function (vecInt,scales) {
    
    z <- cwt(vecInt,scales=scales,wavelet="mexh")
    z <- Re(z)
    attr(z, "scale") <- scales
    attr(z, "time") <- c(1:length(vecInt))
    attr(z, "wavelet") <- "gaussian2"
    attr(z, "series") <- vecInt
    attr(z, "sampling.interval") <- 1
    attr(z, "series.name") <- deparse(substitute(vecInt))
    attr(z, "n.sample") <- length(vecInt)
    attr(z, "n.scale") <- length(scales)
    attr(z, "filter.arg") <- 1
    oldClass(z) <- "wavCWT"
    z
}







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

