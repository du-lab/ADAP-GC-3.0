
readParaSQLite<-function(paraDB)
{
	
	# create a SQLite instance and create one connection.
   m <- dbDriver("SQLite")
   con <- dbConnect(m, dbname = paraDB)
  # query
   rs <- dbGetQuery(con, "select * from options where pname='File'")
   # clean up
   dbDisconnect(con)

   params<-NULL

	for(i in 1:nrow(rs))
	{
		paraName<-rs$pname[i]
		paraValue<-rs$pvalue[i]
		expr<-paste("params$DataFiles$",paraName,i,"=paraValue",sep="")
		eval(parse(text=expr))
	}
	
	con <- dbConnect(m, dbname = paraDB)
	# query
   rs <- dbGetQuery(con, "select * from options where pname<>'File'")
   # clean up
   dbDisconnect(con)
	for(i in 1:nrow(rs))
	{
		paraName<-rs$pname[i]
		paraValue<-rs$pvalue[i]
		expr<-paste("params$",paraName,"=paraValue",sep="")
		eval(parse(text=expr))
	}
	

	return(params)

}
getDataFileFullPath<-function(params)
{
	DataFilelist=params$DataFiles
	for(i in 1:length(DataFilelist))
		DataFilelist[[i]]<-paste(params$DataDir,DataFilelist[[i]],sep="")
	DataFilelist
}

readParaFromCmd<-function(args)
{
	
	
	library(RSQLite)
#	dbFile<-"~/testing/KQC2/KQC.db"
	dbFile<-args[1] # the first arg is dbFile name (includes path)
	params<-readParaSQLite(dbFile)
	params$DataFiles<-getDataFileFullPath(params)
	params$clustType="MPI"
	params$dbFile<-dbFile
	params
	
}

split.matrix<-function(x,f) {
	#print('processing matrix')
	v=lapply(
			1:dim(x)[[2]]
			, function(i) {
				base:::split.default(x[,i],f)#the difference is here
			}
	)
	
	w=lapply(
			seq(along=v[[1]])
			, function(i) {
				result=do.call(
						cbind
						, lapply(v,
								function(vj) {
									vj[[i]]
								}
						)
				)
				colnames(result)=colnames(x)
				return(result)
			}
	)
	names(w)=names(v[[1]])
	return(w)
}

parseFileName<-function(filePath,isPath=T)
{
	if(isPath)
	{
		fileName<-strsplit(filePath,"/")[[1]]#parse the path
		fileName<-fileName[length(fileName)]#get last element as filename
	}else
	{
		fileName<-filePath
	}
	
	sepFilename<-strsplit(fileName,".",fixed=TRUE)[[1]]
	extName<-sepFilename[length(sepFilename)]
	fileName<-substr(fileName,0,nchar(fileName)-nchar(extName)-1)#remove extention name
	return(fileName)
}

dotProduct<- function(x, y) 
{
	if(length(x)<=10||length(y)<=10)
	{
		dotprod<-0
	}else
	{
		dotprod<-(x%*%y)/((x%*%x)^0.5*(y%*%y)^0.5)
		dotprod<-as.numeric(dotprod[1])
	}
	
	angle<-acos(dotprod)*180/pi
	if(is.na(angle))
		angle<-90
	angle
	#dotprod
}

f<- function(x, y) 
{
	(1-cor(x,y))*1000
}

distCal<-function (x) 
{
	
	
	ncy <- ncx <- ncol(x)
	if (ncx == 0) 
		stop("'x' is empty")
	r <- matrix(0, nrow = ncx, ncol = ncy)
	for (i in seq_len(ncx)) {
		for (j in seq_len(i)) {
			x2 <- x[, i]
			y2 <- x[, j]
			r[i, j] <- ifelse(i==j,0,dotProduct(x2, y2))
		}
	}
	r <- r + t(r) - diag(diag(r))
	rownames(r) <- colnames(x)
	colnames(r) <- colnames(x)
	r
	
	
	
}
