#=========================================================================
# To read unit GC-MS data
# split EIC and create re-constructed EIC files ready for peak picking
# Date: Feb 27 2013, written by Yan Ni
#=========================================================================
splitEIC_UnitMass <- function(params) {
	library("ncdf")
	
	WorkDir <- params$WorkDir
	DataFilelist <- params$DataFiles
	
	for(fileindex in 1:length(DataFilelist)) {
		inFilePath <- DataFilelist[[fileindex]]
		fileName<-parseFileName(inFilePath)
		
		#read raw data
		ncid <- open.ncdf(inFilePath)
		intensity <- get.var.ncdf(ncid, varid="intensity_values")
		mass <- get.var.ncdf(ncid, varid="mass_values")
		scan_index <- get.var.ncdf(ncid, varid="scan_index")
		scan_time <- get.var.ncdf(ncid, varid="scan_acquisition_time")
		#total_intensity <- get.var.ncdf(ncid, varid="total_intensity")
		
		# get unique integer mass (mz vector)
		mzVec <- sort(unique(mass))
		
		# get intensity vector
		intVec <- NULL
		total_scan <- length(scan_index)
		full_list <- data.frame(scan=c(1:total_scan)) # only used for filling zeros
		
		time1<-Sys.time()
		for (mz in mzVec) {
			print (mz)
			scan_vec <- NULL
			ind <- which(mass==mz)
			
			if (length(ind) >=1) {	
				for (j in ind) {
					#scan_vec <- c(scan_vec,which(scan_index>=j)[1]-1)
					#scan_vec <- c(scan_vec,binarysearch(j,scan_index,1,total_scan))
					scan_vec <- c(scan_vec,findInterval(j,scan_index)) # find the exact scan number
				}
			}
			
			# to save intensity for current mass
			int_vec <- intensity[ind]
			curmass <- data.frame(scan=scan_vec,int=int_vec)
			curmass <- curmass[order(curmass$scan),] 
			# fill zero
			curmass_block <- merge(full_list, curmass,by=1,all.x=T)
			curmass_block$int[is.na(curmass_block$int)]=0
			intVec <- c(intVec,curmass_block$int)
			#z <- with(curmass, merge(zoo(int, scan), zoo(, c(1:total_scan)), fill = 0))
			#intVec <- c(intVec,coredata(z))
		}
		
		time2<-Sys.time()
		time2 - time1
		
		
		# create EIC cdf file
		cat(fileName,"write out EIC data...\n")
		# create dimension vector and define netCDF variables
		dim_intVec <- dim.def.ncdf("int", units="count", vals=1:length(intVec), unlim=FALSE)
		dim_mzVec <- dim.def.ncdf("mz", units="count", vals=1:length(mzVec), unlim=FALSE)
		var_intVec <- var.def.ncdf("intVec",units="count",dim=dim_intVec,missval=NA,longname="intVec",prec="integer")
		var_mzVec <- var.def.ncdf("mzVec",units="count",dim=dim_mzVec,missval=NA,longname="mzVec",prec="integer")
		# define netcdf file and write data into it
		ncid1<-create.ncdf(filename=paste(WorkDir,"/output/EIC/",fileName,"EIC.cdf",sep=""),vars=list(var_intVec,var_mzVec))
		put.var.ncdf(ncid1,var_intVec,intVec)
		put.var.ncdf(ncid1,var_mzVec,mzVec)
		close.ncdf(ncid1)
		remove(ncid1)	
	}
}

binarysearch <- function(val,tab,L,H) {
		while (H>=L) { 
			M=L+(H-L) %/% 2 
			if (tab[M]>val) 
				H<-M-1 
			else if (tab[M]<val) 
				L<-M+1 
			else return(M)}
		return(L-1)
	} 
	


