# Date: April 4, 2012
# To implement MSAP in ADAP-GC 2.0 GUI

pretreatment_process <- function(params) {
	# get all necessary parameters
	library(RSQLite)

	JobName<-params$JobName
	WorkDir <- params$WorkDir
	paraDB<- params$dbFile
	QualQuan <- params$QualQuanTable_Name
	
	# read it from database, it could be updated after checking manually
#	m <- dbDriver("SQLite")
#	con <- dbConnect(m, dbname = paraDB)
#	tname<-paste(JobName,QualQuan,"BioMarkerTable",sep="_")
#	
#	if(dbExistsTable(con,tname)) {
#		datafile <- dbReadTable(con,tname)
#	}else {
#		stop(paste("Table '",tname,"' does not exist!",sep=""))
#	}
	
	#read it directly from csv file
	tname<-paste(JobName,QualQuan,"BiomarkerTable.csv",sep="_")
	datafile <- read.csv(paste(WorkDir,"output/",tname,sep=""),check.names=F)
	 
	factor_file_name <- paste(WorkDir,"output/factor.csv",sep="")
	
	#temp
	#data_file_name <- paste(WorkDir,"output/stds_M_BiomarkerTable.csv",sep="")
	NormType <- params$Normalization
	ScaleType <- params$Scalling
	
	# combine datafile and factor file together as a list
	Combined_table <- data_org2(datafile,factor_file_name)
	# select normalization and scalling methods
	pretreatment_file <- pretreatment2(Combined_table,NormType,ScaleType)
}

data_org2 <- function(datafile,factor_file_name)
{
	# Input: datafile-raw data from data preprocessing
	# Format: row - samples, column - variables, containing "Y"
	# 1st row: primary name, 1 column: sample name
	# Return: organzied list 

	rawdata <- datafile
	temp.varName <- rawdata$clustID
	Xmatrix <- as.data.frame(t(rawdata[,-c(1:3)]))
	colnames(Xmatrix) <- temp.varName
	
	factor <- read.csv(factor_file_name)$Y
	#factor <- factor[order(rownames(factor)),]
	
	OrgData <- list(X=Xmatrix,Y=factor)
	print("Good start...your data has been produced!\n")
	return (OrgData)
}

pretreatment2 <- function(org.data,norm.type,ScaleType) { 
	datafile <- org.data$X
	if (norm.type == "Internal Standard") {
		datafile <- IsNorm(params,org.data)
	}
	if (norm.type == "Percentage") datafile <- PencentageNorm(org.data)
	if (norm.type == "none") datafile <- datafile	
	org.data$X <- datafile
	print("Normalization is done!\n")
	pretreatment_data <- scaling(org.data,ScaleType) 
	print("scaling is done!\n")
	pretreatment_data
}

##-----------Normalization -------------------------------

IsNorm <- function(params,org.data)
# datafile: data matrix
# is_name: name of most stable standard
{ 
	dataset <- org.data$X
	Nrow <- nrow(dataset)  
	Ncol <- ncol(dataset) 
	
	StdID <- as.character(params$StandardID)
	InStand<- dataset[,StdID]
#	expr<- paste("dataset",StdID,sep="$")
#	InStand<- eval(parse(text=expr)) 
	#mean_inStand <- mean(InStand)
	InStand[which(InStand==0)] = mean(InStand)
	#InStand<- datafile[,i]		
	ISMatrix <- matrix(rep(InStand,ncol(dataset)),ncol=ncol(dataset))
	dataset <- dataset/ISMatrix
	# remove specific internal standard finally
	dataset <- dataset[,-which(colnames(dataset)==StdID)]
	dataset
}

##  Percentage Normalization

PencentageNorm <- function(org.data) 
{
	datafile <- org.data$X
	Nrow <- nrow(datafile)  # the sample number
	Ncol <- ncol(datafile)  # the variable number
	sumlist<-rowSums(datafile)	## this calculation should not put into for loop, changing all the time
	for (i in 1:Nrow) {
		for (j in 1:Ncol) {		
			datafile[i,j] <- (datafile[i,j]/sumlist[i])*100
		}
	}
	datafile
}

##-------------scaling -------------------------------
# nor.data: data matrix after normalizaton
# scal.type: scalling type

scaling<-function(norm.data,scal.type) 
{
	datafile <- norm.data$X
	Nrow <- nrow(datafile)  # the sample number
	Ncol <- ncol(datafile)  # the variable number
	meanlist<-colMeans(datafile)	## this calculation should not put into for loop, changing all the time	
	for (i in 1:Ncol) 
	{
		# to avoid NAN when all column number equals zero
		if (sum(datafile[,i])==0) {
			datafile[,i] = datafile[,i] 
		}
		else 
		{
			std <- sd(datafile[,i])	
			if (scal.type=="centering")for (j in 1:Nrow) 
				{datafile[j,i] <- datafile[j,i] - meanlist[i]}
			if (scal.type=="auto scaling")for (j in 1:Nrow) 
				{datafile[j,i] <- (datafile[j,i] - meanlist[i])/std}	
			if (scal.type=="pareto scaling") for (j in 1:Nrow) 		
				{datafile[j,i] <- (datafile[j,i] - meanlist[i])/(std^0.5)}		
			if (scal.type=="vast scaling")  for (j in 1:Nrow) 			
				{datafile[j,i] <- ((datafile[j,i] - meanlist[i])*meanlist[i])/(std^2)}		
			if (scal.type=="range scaling") for (j in 1:Nrow) 		
				{datafile[j,i] <- (datafile[j,i] - meanlist[i])/(imax-imin)}
			if (scal.type=="level scaling") 
			{datafile[j,i] <- (datafile[j,i] - meanlist[i])/meanlist[i]}	
			if (scal.type=="none") for (j in 1:Nrow)  
				{datafile[j,i] <- datafile[j,i]}
		}
	}
	norm.data$X <- datafile
	norm.data
}
