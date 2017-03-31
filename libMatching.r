##############################
# Libray matching & evaluation
##############################

#add weight to intensity
weightSpec<-function(spec)
{

	n<-1
	m<-1
	sumInt<-sum(spec$int)
	# The weight formular refers to publised paper
	# "An integrated Method for Spectrum... S.E.Stein, 1999"
	w=1/(0.5+sumInt-1)
	weightFactorVec=1/(1+w*spec$int);
    # weighted spectrum considers original intensity, mz and weight factor	
	spec$int<-spec$int^n*spec$mz^m*weightFactorVec
	spec

}

# used for score calculation
RT_score<-function(RT1,RT2,maxRTShift)
{
	999*(1-abs(RT1-RT2)/maxRTShift)
}

# Two dot product calculation with no big difference
# usually keep the same with NIST method
regularDotProduct<-function(x,y)
{
	dotProd<-(x%*%y)/((x%*%x)^0.5*(y%*%y)^0.5)
	999*dotProd
	#angle<-acos(as.numeric(dotProd[1]))*180/pi
	#angle	
}

NISTDotProduct<-function(x,y)
{
	dotProdNew<-sum((x*y)^0.5)^2/(sum(x)*sum(y))
	999*dotProdNew
	#angle<-acos(as.numeric(dotProdNew[1]))*180/pi
	#angle
}

#use union set of two spectra
impureScore<-function(spec1,spec2,isNist)
{
	mzMergeSet<-merge(spec1,spec2,by.x="mz",by.y="mz",all.x=TRUE,all.y=TRUE)
	colnames(mzMergeSet)[2]="x"
	colnames(mzMergeSet)[3]="y"
	x<-mzMergeSet$x
	y<-mzMergeSet$y
	y[is.na(y)]<-0
	x[is.na(x)]<-0
	
	ifelse(isNist,NISTDotProduct(x,y),regularDotProduct(x,y))
}

#use intersecion set of two spectra
pureScore<-function(spec1,spec2,isNist)
{
	mzMergeSet<-merge(spec1,spec2,by.x="mz",by.y="mz")
	if(nrow(mzMergeSet)>=2)
	{
		colnames(mzMergeSet)[2]="x"
		colnames(mzMergeSet)[3]="y"
		x<-mzMergeSet$x
		y<-mzMergeSet$y
		y[is.na(y)]<-0
		x[is.na(x)]<-0
		ifelse(isNist,NISTDotProduct(x,y),regularDotProduct(x,y))
	}else{
		0#return maximum angle if no common fragments or the number of common fragments is less than 2
	}
	
}

# calculate pairwise spectra similarity
# only take two columns structure such as:(mz,int)
# The total score is calculated using both impure and pure score
# the ratio is 3:7

specDistCal<-function(spec1,spec2,isWeight,isNist)
{
	if(!(ncol(spec1)==2&&ncol(spec2)==2))stop("not valid spectrum!")
	#normalize
	spec1$int<-(spec1$int/max(spec1$int))*999
	spec2$int<-(spec2$int/max(spec2$int))*999
	
	if(isWeight)
	{	
		#scale and weight the intentisty
		spec1<-weightSpec(spec1)
		spec2<-weightSpec(spec2)
	}
	
	#remove zero values
	spec1<-subset(spec1,int>0)
	spec2<-subset(spec2,int>0)
	0.3*impureScore(spec1,spec2,isNist)+0.7*pureScore(spec1,spec2,isNist)

}

specDistCal_check2<-function(spec1,spec2,isWeight,isNist)
{
	if(!(ncol(spec1)==2&&ncol(spec2)==2))stop("not valid spectrum!")
	#normalize
	spec1$int<-(spec1$int/max(spec1$int))*999
	spec2$int<-(spec2$int/max(spec2$int))*999
	
	if(isWeight)
	{	
		#scale and weight the intentisty
		spec1<-weightSpec(spec1)
		spec2<-weightSpec(spec2)
	}
	
	#remove zero values
	spec1<-subset(spec1,int>0)
	spec2<-subset(spec2,int>0)
	1*impureScore(spec1,spec2,isNist)+0*pureScore(spec1,spec2,isNist)
	
}



specDistCal_check<-function(spec1,spec2,isWeight,isNist)
{
	if(!(ncol(spec1)==2&&ncol(spec2)==2))stop("not valid spectrum!")
	#normalize
	spec1$int<-(spec1$int/max(spec1$int))*999
	spec2$int<-(spec2$int/max(spec2$int))*999
	
	if(isWeight)
	{	
		#scale and weight the intentisty
		spec1<-weightSpec(spec1)
		spec2<-weightSpec(spec2)
	}
	
	#remove zero values
	spec1<-subset(spec1,int>0)
	spec2<-subset(spec2,int>0)
	pureScore(spec1,spec2,isNist)
	
}


#read msp file and save spectra in memory(mz,int,RT)
readMSP2Spec<-function(filename,withRT)	
{
	cat(paste("read ",filename,"...\n"))

	
	#read lib
	mspFile <- file(filename, "r")
	nCounter=0
	compList<-vector("list",length=0)
	compName<-NULL
	curComp<-NULL
	nFragments<-0
	RT<-0
	
	curLine<-readLines(mspFile,n=1)
	
	while(length(curLine<-readLines(mspFile,n=1))>0)
	{
		#read one line from file
		#extract lines starting with "Name","Num Peaks" or "Synon: RT/Retention time"
		#And extract spectrum data
		
		if(grepl(":",curLine)==TRUE)#if it is the meta data
		{
			#normalize the string by removing '#' and replacing '=' to ':'
			curLine<-sub('#+',"",curLine)
			curLine<-sub("=",":",curLine)
			splitArrays<-strsplit(curLine,":")[[1]]#parse the path
			field<-splitArrays[[1]]
			if(field=="Name")#new compound
			{
				#if not the first compound then save the previous comp into last list element
				if(nCounter>0)
				{
					compList[[nCounter]]<-curComp
					#compList[[nCounter]]<-curComp[-which(is.na(curComp$mz)==TRUE),]
					names(compList)[nCounter]<-compName
					curComp<-NULL#then clear current comp
				}
				#save compound name to current list element name
				nCounter=nCounter+1	
				compName<-trim(splitArrays[[2]])

			}else
			{
				if(field=="Num Peaks")
				{
					nPeaks<-trim(splitArrays[[2]])
					if(withRT)
					{
						curComp<-data.frame(mz=rep(0,nPeaks),int=rep(0,nPeaks),RT=rep(RT,nPeaks))#create new data frame
					}else
					{
						curComp<-data.frame(mz=rep(0,nPeaks),int=rep(0,nPeaks))#create new data frame
					}
					nFragments<-0#init the counter of fragments
				}
				if(withRT)
				{
					if(field=="Synon")
					{
						field2<-trim(splitArrays[[2]])
						if(field2=="RT(s)"&&!is.na(as.numeric(splitArrays[[3]]))) #if it is RT, add the RT column to component data frame
							RT<-as.numeric(splitArrays[[3]])
						if(field2=="RT(m)"&&!is.na(as.numeric(splitArrays[[3]])))
							RT<-as.numeric(splitArrays[[3]])*60
						if(field2=="Retention time"&&!is.na(as.numeric(splitArrays[[3]])))
							RT<-as.numeric(splitArrays[[3]])*60
						if(field2=="Retention second"&&!is.na(as.numeric(splitArrays[[3]])))
							RT<-as.numeric(splitArrays[[3]])
						
					}
				}
			}
		}else #spectrum data
		{
			SpecPoints<-strsplit(curLine,";")[[1]]#parse the pairs
			if(length(SpecPoints)>0)#skip empty line
			{
				#assign the mz and int from current line
				for(i in 1:length(SpecPoints))
				{
					nFragments=nFragments+1#increase counter of fragments
					curFragment<-trim(SpecPoints[i])#get current pair
					curFragment<-strsplit(curFragment," ")[[1]]#split the current pair
					if(length(curFragment)>0)
					{
						curComp[nFragments,]$mz<-as.numeric(curFragment[1])
						curComp[nFragments,]$int<-as.numeric(curFragment[length(curFragment)])
					}
				}
			}
			RT<-0###reset RT
		}

		#}
	}
	#add the last compound to the list
	compList[[nCounter]]<-curComp
	#compList[[nCounter]]<-curComp[-which(is.na(curComp$mz)==TRUE),]
	names(compList)[nCounter]<-compName

	close(mspFile)
	compList

	#close(pb)#close progress bar	
}

# That should be used for alignment to pairwisely compare the spectrum similarity
# which has been converted into C
# To extract each component from each file
# eg # of file = 15, then it will output 15*15 score matrix

specDistMatrix<-function(componentFeatures,addRTScore)
{
	library(gtools)
	FileIds<-sort(unique(componentFeatures$FileID))
	ncy <- ncx <- length(FileIds)
   	r <- matrix(0, nrow = ncx, ncol = ncy)
	indexPairs<-combinations(ncx,2)

	for(ind in 1:nrow(indexPairs))
	{
		indexPair<-indexPairs[ind,]
				
		i<-indexPair[1]
		j<-indexPair[2]
		spec1<-subset(componentFeatures,FileID==FileIds[i],select=c("mz","int","compoundET"))
		spec2<-subset(componentFeatures,FileID==FileIds[j],select=c("mz","int","compoundET"))
		
		RT1=spec1$compoundET[1]
		RT2=spec2$compoundET[1]
	
		spec1<-subset(spec1,select=c("mz","int"))
		spec2<-subset(spec2,select=c("mz","int"))	
		
		#calculate the distance between each pair of profiles
		v1<-specDistCal(spec1,spec2,isWeight=TRUE,isNist=TRUE)
		
		maxRTShift<-4/60
		v2<-RT_score(RT1,RT2,maxRTShift)
		
		#the final score has two options,considering RTscore or not
		v<-ifelse(addRTScore,0.9*v1+0.1*v2,v1)
		if(is.na(v))
		{
			v<-0
			warning(paste(i,j,v,":dot product is NA!replaced by 0"))
			
			}
		r[i,j]<-as.numeric(format(v, digits=2))
	}	


    r <- r + t(r) - diag(diag(r)) #fill another half symetric values of matrix

	columnNames<-FileIds
	rownames(r) <- columnNames# assign feature ID as col and row names
    colnames(r) <- columnNames
	r
}
##now java version is implemented and it is currently not used by JavaGUI
MSP2SQLite<-function(MSPFile,paraDB,withRT)
{
	library(RSQLite)
	components<-readMSP2Spec(MSPFile,withRT)
	tableName<-paste("spec",parseFileName(MSPFile),sep="_")
	nCounter<-0;
	m <- dbDriver("SQLite")
	con <- dbConnect(m, dbname = paraDB)

	res <- dbSendQuery(con, paste("drop table if exists ",tableName,sep=""))
	res <- dbSendQuery(con, paste("create table ",tableName,"(compoundName text",ifelse(withRT,",RT double",""),",spec text)",sep=""))
	total<-length(components)
	dbBeginTransaction(con)
	for(i in 1:total)
	{
		#get current cluster

		curComponent<-components[[i]]
		
		
		#output the best spectrum 

		componentName<-names(components[i])
		if(withRT)
			curRT<-curComponent$RT[1]
		
		
		nCounter<-nCounter+1;
		strOutput<-getSpecBlock(curComponent,isFormat=FALSE)
	
		res <- dbSendQuery(con, paste("insert into ",tableName," values('",componentName,"'",ifelse(withRT,paste(",'",curRT,"'",sep=""),""),",'",strOutput,"')",sep=""))
		print(paste(round(i/total*100, 0),  "% done"))
		flush.console()
	}
	dbCommit(con)
	# clean up
   	dbDisconnect(con)
	cat("Spectra saved! ","...\n")
}

##now java version is implemented and it is currently not used by JavaGUI

refSpec2SQLite<-function(params,paraDB)
{
	library(RSQLite)
	sampleRatio<-as.numeric(params$sampleRatio)
	JobName<-params$JobName
	filterFile<-params$filterFile
	Tolerance<-as.numeric(params$RT_Tolerance)
	filterType<-params$filterType
	WorkDir <- params$WorkDir
	alignedFile<-params$AlignedFile
	filterType<-ifelse(is.na(filterType),"none",filterType)

	components<-read.csv(file=alignedFile)
	######change int names to avoid reserved keyword in db system
	names(components)[names(components)=="int"]<-"intensity"
	#####save aligned results to db
	m <- dbDriver("SQLite")
	con <- dbConnect(m, dbname = paraDB)
	#res <- dbSendQuery(con, "drop table if exists AlignedPeaks")
	dbWriteTable(con,"AlignedPeaks",components,overwrite=T)
	
	FileIDs<-unique(components$FileID)
	nFile<-length(FileIDs)
	fileVsClust<-subset(components,select=c("FileID","clustID"))#get common mother Ions which occur in 80% of datasets
	fileVsClust<-unique(fileVsClust)
	groupCount<-aggregate(fileVsClust$FileID,by=list(fileVsClust$clustID),length)#group count fileid by clustid
	nMajority<-round(sampleRatio*nFile)
	CommonClustID<-subset(groupCount,x>=nMajority)$Group.1# get common clustid which occurs in over 4 datasets
	components<-subset(components,clustID%in%CommonClustID)#filter featureset by common clustid
	
	if(filterType=="RT")
	{
		cat(paste("filtering results by RT file:",filterFile,"...\n"))
		RTFilter<-read.csv(filterFile)
		RTFilter<-RTFilter$RT
		RTFilter<-RTFilter[!is.na(RTFilter)]
		RTFilter<-RTFilter*60
	
	}else{
		if(filterType=="clustID")
		{
			cat(paste("filtering results by clustID:",strClustID,"...\n"))
			clustIDvec<-unlist(strsplit(strClustID,split=","))
			components<-subset(components,clustID%in%clustIDvec)

		}
		
	}
	
	ClustIDVec=unique(components$clustID)
	total <- length(ClustIDVec)# create progress bar

	nCounter<-0;
		 # query
	res <- dbSendQuery(con, "drop table if exists ExpSpectra")
	res <- dbSendQuery(con, "create table ExpSpectra(clustID int,spec text)")
	####for NIST text output
	outputPath<-paste(WorkDir,"/output/",JobName,"_Nist.txt",sep="")
	strNistOutput<-NULL	

	dbBeginTransaction(con)
	for(i in 1:length(ClustIDVec))
	{
		#get current cluster
		curClustID<-ClustIDVec[i]
		curComponent<-subset(components,clustID==curClustID&isRef==1)
		
		
		#output the best spectrum 
		curCompoundET<-unique(curComponent$alignedET)
		componentName<-paste(format(curCompoundET,digit=4)," clust",curClustID,sep="")
		

		if(filterType=="RT")
			{
				if(length(RTFilter[abs(RTFilter-curCompoundET*60)<RT_Tolerance])>0)
				{
				nCounter<-nCounter+1;
				######for db output
				strOutput<-getMSPBlock(curComponent,isFormat=FALSE)
				####for NIST text output
				strNistOutput<-paste(strNistOutput,getMSPBlock(curComponent,curCompoundET*60,componentName,nCounter),sep="\n\n")
				}	
			}else{
		
			nCounter<-nCounter+1;
			######for db output
			strOutput<-getSpecBlock(curComponent,isFormat=FALSE)
			####for NIST text output
			strNistOutput<-paste(strNistOutput,getMSPBlock(curComponent,curCompoundET*60,componentName,nCounter),sep="\n\n")
			}
			res <- dbSendQuery(con, paste("insert into ExpSpectra values(",curClustID ,",'",strOutput,"')",sep=""))
			print(paste(round(i/total*100, 0),  "% done"))
			flush.console()
	}
	dbCommit(con)
	# clean up
   	dbDisconnect(con)
	####for NIST text output
	write(strNistOutput,file=outputPath)
	cat("Spectra saved! ","...\n")
}

#select reference components in the aligned data and save them into MSP file
refSpec2MSP<-function(sampleRatio,JobName,params,filterFile="None",RT_Tolerance,strClustID,filterType)
{
	library(tcltk)
	#JobName<-params$JobName
	WorkDir <- params$WorkDir
	alignedFile<-params$AlignedFile
	outputPath<-paste(WorkDir,"/output/",JobName,"_Nist.txt",sep="")

	components<-read.csv(file=alignedFile)
	
	FileIDs<-unique(components$FileID)
	nFile<-length(FileIDs)
	fileVsClust<-subset(components,select=c("FileID","clustID"))#get common mother Ions which occur in 80% of datasets
	fileVsClust<-unique(fileVsClust)
	groupCount<-aggregate(fileVsClust$FileID,by=list(fileVsClust$clustID),length)#group count fileid by clustid
	nMajority<-round(sampleRatio*nFile)
	CommonClustID<-subset(groupCount,x>=nMajority)$Group.1# get common clustid which occurs in over 4 datasets
	components<-subset(components,clustID%in%CommonClustID)#filter featureset by common clustid

	if(filterType=="RT")
	{
		cat(paste("filtering results by RT file:",filterFile,"...\n"))
		RTFilter<-read.csv(filterFile)
		RTFilter<-RTFilter$RT
		RTFilter<-RTFilter[!is.na(RTFilter)]
		RTFilter<-RTFilter*60
	
	}else{
		if(filterType=="clustID")
		{
			cat(paste("filtering results by clustID:",strClustID,"...\n"))
			clustIDvec<-unlist(strsplit(strClustID,split=","))
			components<-subset(components,clustID%in%clustIDvec)

		}
		
	}

	ClustIDVec=unique(components$clustID)
	total <- length(ClustIDVec)# create progress bar
	pb <- tkProgressBar(title = paste(JobName,":output to Nist file"), min = 0,
			                    max = total, width = 300)

	strOutput<-NULL
	nCounter<-0;
	for(i in 1:length(ClustIDVec))
	{
		#get current cluster
		curClustID<-ClustIDVec[i]
		curComponent<-subset(components,clustID==curClustID&isRef==1)
		
		#calculate and get the best spectrum
#		r<-specDistMatrix(curComponents,addRTScore=FALSE)
#		bestfileID<-colnames(r)[which.max(colMeans(r))]
#		curComponent<-subset(curComponents,FileID==bestfileID)
		#output the best spectrum 
		curCompoundET<-unique(curComponent$alignedET)
		componentName<-paste(format(curCompoundET,digit=4)," clust",curClustID,sep="")
		

		if(filterType=="RT")
			{
				if(length(RTFilter[abs(RTFilter-curCompoundET*60)<RT_Tolerance])>0)
				{
				nCounter<-nCounter+1;
				strOutput<-paste(strOutput,getMSPBlock(curComponent,curCompoundET*60,componentName,nCounter),sep="\n\n")
				}	
			}else{
		
			nCounter<-nCounter+1;

			strOutput<-paste(strOutput,getMSPBlock(curComponent,curCompoundET*60,componentName,nCounter),sep="\n\n")
			}
		setTkProgressBar(pb, i, label=paste(round(i/total*100, 0),  "% done"))
	}
	write(strOutput,file=outputPath)
	close(pb)#close progress bar	
	cat(paste(parseFileName(outputPath),"saved! ","...\n"))
}

refSpec2MSP_SQLite<-function(params)
{

	sampleRatio<-as.numeric(params$sampleRatio)
	JobName<-params$JobName
	filterFile<-params$filterFile
	Tolerance<-as.numeric(params$RT_Tolerance)
	filterType<-params$filterType
	WorkDir <- params$WorkDir
	alignedFile<-params$AlignedFile
	filterType<-ifelse(is.na(filterType),"none",filterType)

	outputPath<-paste(WorkDir,"/output/",JobName,"_Nist.txt",sep="")

	components<-read.csv(file=alignedFile)
	
	FileIDs<-unique(components$FileID)
	nFile<-length(FileIDs)
	fileVsClust<-subset(components,select=c("FileID","clustID"))#get common mother Ions which occur in 80% of datasets
	fileVsClust<-unique(fileVsClust)
	groupCount<-aggregate(fileVsClust$FileID,by=list(fileVsClust$clustID),length)#group count fileid by clustid
	nMajority<-round(sampleRatio*nFile)
	CommonClustID<-subset(groupCount,x>=nMajority)$Group.1# get common clustid which occurs in over 4 datasets
	components<-subset(components,clustID%in%CommonClustID)#filter featureset by common clustid
	
	if(filterType=="RT")
	{
		cat(paste("filtering results by RT file:",filterFile,"...\n"))
		RTFilter<-read.csv(filterFile)
		RTFilter<-RTFilter$RT
		RTFilter<-RTFilter[!is.na(RTFilter)]
		RTFilter<-RTFilter*60
	
	}else{
		if(filterType=="clustID")
		{
			cat(paste("filtering results by clustID:",strClustID,"...\n"))
			clustIDvec<-unlist(strsplit(strClustID,split=","))
			components<-subset(components,clustID%in%clustIDvec)

		}
		
	}
	
	ClustIDVec=unique(components$clustID)
	total <- length(ClustIDVec)# create progress bar
	strOutput<-NULL
	nCounter<-0;
	for(i in 1:length(ClustIDVec))
	{
		#get current cluster
		curClustID<-ClustIDVec[i]
		curComponent<-subset(components,clustID==curClustID&isRef==1)
		
		#calculate and get the best spectrum
#		r<-specDistMatrix(curComponents,addRTScore=FALSE)
#		bestfileID<-colnames(r)[which.max(colMeans(r))]
#		curComponent<-subset(curComponents,FileID==bestfileID)
		#output the best spectrum 
		curCompoundET<-unique(curComponent$alignedET)
		componentName<-paste(format(curCompoundET,digit=4)," clust",curClustID,sep="")
		

		if(filterType=="RT")
			{
				if(length(RTFilter[abs(RTFilter-curCompoundET*60)<RT_Tolerance])>0)
				{
				nCounter<-nCounter+1;
				strOutput<-paste(strOutput,getMSPBlock(curComponent,curCompoundET*60,componentName,nCounter),sep="\n\n")
				}	
			}else{
		
			nCounter<-nCounter+1;
			strOutput<-paste(strOutput,getMSPBlock(curComponent,curCompoundET*60,componentName,nCounter),sep="\n\n")

			}
			print(paste(round(i/total*100, 0),  "% done"))
			flush.console()
	}
	write(strOutput,file=outputPath)
	cat(paste(parseFileName(outputPath),"saved! ","...\n"))
}


#convert spec of dataframe format (mz,int) into Nist MSP format
getSpecBlock<-function(curComponent,isFormat=TRUE)
{

		strMzBlock<-NULL
		curComponent<-curComponent[order(curComponent$mz),]#order by mz
		curComponent$int<-as.numeric(curComponent$int)
		maxInt<-max(curComponent$int)
		curComponent$int<-round((curComponent$int/maxInt)*999)#normalize by maximun int
		curComponent<-subset(curComponent,int>0)#filter the mass with low intensity which is zero after normalization
		nPeak<-nrow(curComponent)
		
		for(j in 1:nrow(curComponent))
		{
			strPair<-paste(curComponent[j,]$mz," ",curComponent[j,]$int,";",sep="")
			strMzBlock<-paste(strMzBlock,strPair)
			
			if(isFormat&&j%%4==0)strMzBlock<-paste(strMzBlock,"\n",sep="")
		}
		list(Spec=strMzBlock,nPeak=nPeak)
}
#convert spec of dataframe format (mz,int) into Nist MSP format
getMSPBlock<-function(curComponent,RT,componentName,dbind)
{

		strBlock<-sprintf("Name:%s",componentName)
		#RT<-RT*60 #convert to seconds
		strBlock<-paste(strBlock,sprintf("Synon:RT(s):%.1f",RT),sep="\n")
		#strBlock<-paste(strBlock,sprintf("Synon:MCR-Reference:%.2f Clust%d",RT,curClustID),sep="\n")
		strBlock<-paste(strBlock,sprintf("DB#:%d",dbind),sep="\n")

		res<-getSpecBlock(curComponent,isFormat=TRUE)
		strMzBlock<-res[[1]]
		nPeak<-res[[2]]
		strBlock<-paste(strBlock,sprintf("Num Peaks:%d",nPeak),sep="\n")
				
		strBlock<-paste(strBlock,strMzBlock,sep="\n")
		strBlock

}
##read spec table from sqlite db and convert it into data.frame (mz,int) 
loadSpecSQL<-function(tableName,dbName)
{
	library(gdata)
	library(RSQLite)
	m <- dbDriver("SQLite")
	con <- dbConnect(m, dbname = dbName)
	rs <- dbGetQuery(con, paste("select * from",tableName))
	dbDisconnect(con)
	compList<-NULL

	for(i in 1:nrow(rs))
	{
		curComp<-NULL
		SpecPoints<-strsplit(rs$spec[[i]],";")[[1]]#parse the pairs
		if(length(SpecPoints)>0)#skip empty line
		{
			#assign the mz and int from current record
			for(j in 1:length(SpecPoints))
			{
				#nFragments=nFragments+1#increase counter of fragments
				curFragment<-trim(SpecPoints[j])#get current pair
				curFragment<-strsplit(curFragment," ")[[1]]#split the current pair
				if(length(curFragment)>0)
				{
					curComp<-rbind(curComp,data.frame(mz=as.numeric(curFragment[1]),int=as.numeric(curFragment[length(curFragment)])))
					#curComp[nFragments,]$mz<-as.numeric(curFragment[1])
					#curComp[nFragments,]$int<-as.numeric(curFragment[length(curFragment)])
				}
			}
		}
	compList[[i]]<-curComp
	names(compList)[i]<-rs$compoundName[[i]]

	}

	compList
}

libMatching_old<-function(refSpec,inputSpec,minSpecSimilarity)
{
	#library(tcltk)
	library(gdata)

	#IntPower<-1
#	massPower<-3
	
	if(is.character(refSpec))
	{ ######read from NIST file
		lib<-readMSP2Spec(refSpec,withRT=FALSE)
	}else
	{
		lib<-refSpec
		}
	
	####read user spectra
	if(is.character(inputSpec))
	{
		components<-readMSP2Spec(inputSpec,withRT=FALSE)
	}else
	{
		components<-inputSpec
	}		
#

	
	
	total <- length(components)# create progress bar
	#pb <- tkProgressBar(title = paste("matching lib"), min = 0,max = total, width = 400)
	
	candidates<-NULL
	for(i in 1:total)
	{
		curComponent<-components[i]#get current cluster which is one detected compound
		curComponentName<-names(curComponent)
		curComponentName<-strsplit(curComponentName," ")[[1]]
		pkind<-as.numeric(curComponentName[1])
		curWindowID <- curComponentName[2]
		curCompoundID <- curComponentName[3]
		mdpkMass<-curComponentName[4]
		curSpec<-curComponent[[1]]
		if(nrow(curSpec)>0)
		{
			curSpec<-subset(curSpec,select=c("mz","int"))
		
			#matching lib
			curMatch<-data.frame(pkind=pkind,windowID=curWindowID,compoundID=curCompoundID,modelMass=mdpkMass,compoundName="N/A",score=0)#init the match item
			for(j in 1:length(lib))
			{
				curCompound<-lib[j]
				refName<-names(curCompound)
				refSpec<-curCompound[[1]]
				refSpec<-subset(refSpec,select=c("mz","int"))
				
				score<-specDistCal(spec1=curSpec,spec2=refSpec,isWeight=TRUE,isNist=TRUE)
				#update the matched item with the top match
				if(score>curMatch$score)
				{
					curMatch$score<-score
					curMatch$compoundName<-refName
				}
					
			}	
			if(curMatch$score>minSpecSimilarity)
				candidates<-rbind(candidates,curMatch)
		
		}
		#setTkProgressBar(pb, i, label=paste(round(i/total*100, 0),  "% done"))		
	}
		

	#close(pb)#close progress bar	
	candidates

}


# modified libmathing based on libmatching_TopOne
# output all canidiates for each component once it 
# meet the score critera (minSpecSimilarity,eg. 700)

libMatching_Top10<-function(refSpec,inputSpec,minSpecSimilarity)
{
	library(gdata)
	
	if(is.character(refSpec))
	{ ######read from NIST file
		lib<-readMSP2Spec(refSpec,withRT=FALSE)
	}else
	{
		lib<-refSpec
	}
	
	####read user spectra
	if(is.character(inputSpec))
	{
		components<-readMSP2Spec(inputSpec,withRT=FALSE)
	}else
	{
		components<-inputSpec
	}		
	
	total <- length(components)# create progress bar
	candidates<-NULL
	
	# for loop for each component spectrum
	for(i in 1:total)
	{
		curComponent<-components[i]#get current cluster which is one detected compound
		curComponentName<-names(curComponent)
		curComponentName<-strsplit(curComponentName," ")[[1]]
		pkind<-as.numeric(curComponentName[1])
		curWindowID <- curComponentName[2]
		curCompoundID <- curComponentName[3]
		mdpkMass<-curComponentName[4]
		curCompID <- curComponentName[5]
		curSpec<-curComponent[[1]]
		
		if(nrow(curSpec)>0)
		{
			curSpec<-subset(curSpec,select=c("mz","int"))
			
			#initiate matching results (data frame format)
			curMatch<-data.frame(pkind=as.integer(),windowID=as.integer(),compoundID=as.integer(),compID=as.integer(),modelMass=as.integer(),compoundName=as.character(),score=as.integer())
			
			#matching lib from beginning to end
			for(j in 1:length(lib))
			{
				curCompound<-lib[j]
				refName<-names(curCompound)
				refSpec<-curCompound[[1]]
				#refSpec<-subset(refSpec,select=c("mz","int"))
				refSpec<-subset(refSpec,mz!=0)
				
				score<-round(specDistCal(spec1=curSpec,spec2=refSpec,isWeight=TRUE,isNist=TRUE))			
				curMatch <- rbind(curMatch,data.frame(pkind=pkind,windowID=curWindowID,compoundID=curCompoundID,compID=curCompID,modelMass=mdpkMass,compoundName=refName,score=score))
				
			}	
			
			#filtering: output top ten with score >= minSpecSimilarity (eg. 700)
			curMatch <- subset(curMatch,score>=minSpecSimilarity)
			if (nrow(curMatch)>10) {
				candidates <- rbind(candidates,curMatch[order(curMatch$score,decreasing=T)[1:10],])	
			}else{
				candidates <- rbind(candidates,curMatch)
			}	
		}
	}
	candidates
}



libMatching_TopOne<-function(refSpec,inputSpec,minSpecSimilarity)
{
	#library(tcltk)
	library(gdata)
	
	#IntPower<-1
#	massPower<-3
	
	if(is.character(refSpec))
	{ ######read from NIST file
		lib<-readMSP2Spec(refSpec,withRT=FALSE)
	}else
	{
		lib<-refSpec
	}
	
	####read user spectra
	if(is.character(inputSpec))
	{
		components<-readMSP2Spec(inputSpec,withRT=FALSE)
	}else
	{
		components<-inputSpec
	}		
#
	total <- length(components)# create progress bar
	#pb <- tkProgressBar(title = paste("matching lib"), min = 0,max = total, width = 400)
	
	candidates<-NULL
	for(i in 1:total)
	{
		curComponent<-components[i]#get current cluster which is one detected compound
		curComponentName<-names(curComponent)
		curComponentName<-strsplit(curComponentName," ")[[1]]
		pkind<-as.numeric(curComponentName[1])
		curWindowID <- curComponentName[2]
		curCompoundID <- curComponentName[3]
		mdpkMass<-curComponentName[4]
		curCompID <- curComponentName[5]
		curSpec<-curComponent[[1]]
		
		if(nrow(curSpec)>0)
		{
			curSpec<-subset(curSpec,select=c("mz","int"))
			
			#matching lib
			curMatch<-data.frame(pkind=pkind,windowID=curWindowID,compoundID=curCompoundID,compID=curCompID,modelMass=mdpkMass,compoundName="N/A",score=0)#init the match item
			for(j in 1:length(lib))
			{
				curCompound<-lib[j]
				refName<-names(curCompound)
				refSpec<-curCompound[[1]]
				refSpec<-subset(refSpec,select=c("mz","int"))
				
				score<-specDistCal(spec1=curSpec,spec2=refSpec,isWeight=TRUE,isNist=TRUE)
				#update the matched item with the top match
				if(score>curMatch$score)
				{
					curMatch$score<-score
					curMatch$compoundName<-refName
				}
				
			}	
			if(curMatch$score>minSpecSimilarity)
				candidates<-rbind(candidates,curMatch)
			
		}
		#setTkProgressBar(pb, i, label=paste(round(i/total*100, 0),  "% done"))		
	}
	
	
	#close(pb)#close progress bar	
	candidates
	
}

evaluateIdentification<-function(refSpec,inputSpec,RT_Tolerance,minSpecSimilarity,fileName,withRT,params)
{
#	library(tcltk)
	library(gdata)

	#IntPower<-1
#	massPower<-3
	
	if(is.character(refSpec))
	{ ######read from NIST file
		lib<-readMSP2Spec(refSpec,withRT)
	}else
	{
		lib<-refSpec
	}
	
	####read user spectra
	if(is.character(inputSpec))
	{
		components<-readMSP2Spec(inputSpec,withRT)
	}else
	{
		components<-inputSpec
	}		
#
	
	
	total <- length(lib)# create progress bar
#	pb <- tkProgressBar(title = paste("matching lib"), min = 0,
#			                    max = total, width = 400)
	
	candidates<-NULL
	for(i in 1:total)
	{
		curRef<-lib[i]#get current cluster which is one detected compound
		curRefName<-names(curRef)
		curRefName<-gsub(",",".",curRefName)
		#curRefName<-strsplit(curRefName," ")[[1]]
		curRef<-curRef[[1]]
		refRT<-ifelse(withRT,curRef$RT[1],0)
		
		curRefSpec<-curRef[,1:2]

		#matching lib
		curMatch<-data.frame(id=i,RT=refRT/60,compoundName=curRefName,componentName="N/A",score=0)#init the match item

		for(j in 1:length(components))
		{
			curComponent<-components[j]
			componentName<-names(curComponent)
			curComponent<-curComponent[[1]]
			if(nrow(curComponent)>0)
			{
				componentRT<-curComponent$RT[1]
				####if RT info is not in dataframe then get it from list name
				if(is.null(componentRT))
				{
					componentRT<-as.numeric(unlist(strsplit(names(components[j])," "))[1])*60
				}
				
				componentRT<-ifelse(withRT,componentRT,0)
				curSpec<-curComponent[,1:2]
				
				#only match the component within the RT range
				if(abs(componentRT-refRT)<=RT_Tolerance)
				{
					score<-specDistCal(spec1=curRefSpec,spec2=curSpec,isWeight=TRUE,isNist=TRUE)
					#update the matched item with the top match
					if(score>curMatch$score&&score>minSpecSimilarity)		
					{
						curMatch$score<-score
						curMatch$componentName<-componentName
					}
				}
			}
		}	
		#if(curMatch$angle<maxSpecAngle)
		curMatch$score<-round(curMatch$score)
		candidates<-rbind(candidates,curMatch)

		
		
#		setTkProgressBar(pb, i, label=paste(round(i/total*100, 0),  "% done"))		
	}
		
#	Sys.sleep(5)
#	close(pb)#close progress bar	

	#candidates<-as.data.frame(unlist(candidates))

	MatchResultName<-paste(params$WorkDir,"/output/evaluationResult_",fileName,".csv",sep="")
	write.csv(row.names=FALSE,candidates[order(candidates$RT),-1],file=MatchResultName)

}
