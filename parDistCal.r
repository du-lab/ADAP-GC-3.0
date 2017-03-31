parDistCal5<-function (cl,goodShapePeaklist,type="deconvolution",isUnion=TRUE) 
{
	ncy <- ncx <- nrow(goodShapePeaklist)
	
	if (ncx == 0) 
		stop("'x' is empty")
	r <- matrix(0, nrow = ncx, ncol = ncy)
	indexPairs<-combinations(ncx,2)
	indexPairLists<-vector("list",nrow(indexPairs))
	for(i in 1:nrow(indexPairs))
	{
		indexPairLists[[i]]<-indexPairs[i,]
	}
	indexPairsGroups<-clusterSplit(cl,indexPairLists)
	clusterExport(cl,"Global.vecInt")#assign the current profiles to all nodes
	clusterExport(cl,"dotProduct")#assign the current profiles to all nodes
	distResults<-clusterApply(cl,indexPairsGroups,pairDist5,goodShapePeaklist,type,isUnion)#calculate the distance between each pair of profiles
	distResults<-unlist(distResults,recursive=FALSE)
	distResults<-do.call(rbind,distResults)
	fillDistMatrix<-function(distResult)
	{
		
		i<-distResult[1]
		j<-distResult[2]
		v<-distResult[3]
		r[i,j]<<-v
		
	}
	
	res<-apply(distResults,1,fillDistMatrix)#fill the distance matrix with the results
	r <- r + t(r) - diag(diag(r)) #fill another half symetric values of matrix
	columnNames<-names(Global.curProfs)
	rownames(r) <- columnNames# assign feature ID as col and row names
	colnames(r) <- columnNames
	r  
}

DistCal4<-function(type="deconvolution",isUnion=TRUE) 
{
	ncy <- ncx <- length(Global.curProfs)
	
	if (ncx == 0) 
		stop("'x' is empty")
	r <- matrix(0, nrow = ncx, ncol = ncy)
	indexPairs<-combinations(ncx,2)
	indexPairLists<-vector("list",nrow(indexPairs))
	for(i in 1:nrow(indexPairs))
	{
		indexPairLists[[i]]<-indexPairs[i,]
	}
	
	#indexPairsGroups<-clusterSplit(c1,indexPairLists)
	#indexPairsGroups<-lapply(splitIndices(length(indexPairLists), 2), function(i) indexPairLists[i])
	
#	if (nrow(indexPairs)==1) {
#		indexPairsGroups<-lapply(splitIndices(length(indexPairLists), 1), function(i) indexPairLists[i])
#	}else {
#		indexPairsGroups<-lapply(splitIndices(length(indexPairLists), 2), function(i) indexPairLists[i])
#	}
	
	#indexPairsGroups <- indexPairLists
	distResults<- pairDist2(indexPairLists,type,isUnion)#
	
	
	#distResults<-lapply(indexPairsGroups,pairDist2,type,isUnion)#calculate the distance between each pair of profiles
	#distResults<-unlist(distResults,recursive=FALSE)
	distResults<-do.call(rbind,distResults)
	fillDistMatrix<-function(distResult)
	{
		
		i<-distResult[1]
		j<-distResult[2]
		v<-distResult[3]
		r[i,j]<<-v
		
	}
	
	apply(distResults,1,fillDistMatrix)#fill the distance matrix with the results
	r <- r + t(r) - diag(diag(r)) #fill another half symetric values of matrix
	columnNames<-names(Global.curProfs)
	rownames(r) <- columnNames# assign feature ID as col and row names
	colnames(r) <- columnNames
	r  
}
DistCal4_old<-function(type="deconvolution",isUnion=TRUE) 
{
	ncy <- ncx <- length(Global.curProfs)
	
	if (ncx == 0) 
		stop("'x' is empty")
	r <- matrix(0, nrow = ncx, ncol = ncy)
	indexPairs<-combinations(ncx,2)
	indexPairLists<-vector("list",nrow(indexPairs))
	for(i in 1:nrow(indexPairs))
	{
		indexPairLists[[i]]<-indexPairs[i,]
	}
	
	#indexPairsGroups<-clusterSplit(c1,indexPairLists)
	#indexPairsGroups<-lapply(splitIndices(length(indexPairLists), 2), function(i) indexPairLists[i])
	
	if (nrow(indexPairs)==1) {
		indexPairsGroups<-lapply(splitIndices(length(indexPairLists), 1), function(i) indexPairLists[i])
	}else {
		indexPairsGroups<-lapply(splitIndices(length(indexPairLists), 2), function(i) indexPairLists[i])
	}
	
	
	distResults<-lapply(indexPairsGroups,pairDist2,type,isUnion)#calculate the distance between each pair of profiles
	distResults<-unlist(distResults,recursive=FALSE)
	distResults<-do.call(rbind,distResults)
	fillDistMatrix<-function(distResult)
	{
		
		i<-distResult[1]
		j<-distResult[2]
		v<-distResult[3]
		r[i,j]<<-v
		
	}
	
	apply(distResults,1,fillDistMatrix)#fill the distance matrix with the results
	r <- r + t(r) - diag(diag(r)) #fill another half symetric values of matrix
	columnNames<-names(Global.curProfs)
	rownames(r) <- columnNames# assign feature ID as col and row names
	colnames(r) <- columnNames
	r  
}
parDistCal4<-function (cl,type="deconvolution",isUnion=TRUE) 
{
	ncy <- ncx <- length(Global.curProfs)
	
	if (ncx == 0) 
		stop("'x' is empty")
	r <- matrix(0, nrow = ncx, ncol = ncy)
	indexPairs<-combinations(ncx,2)
	indexPairLists<-vector("list",nrow(indexPairs))
	for(i in 1:nrow(indexPairs))
	{
		indexPairLists[[i]]<-indexPairs[i,]
	}
	indexPairsGroups<-clusterSplit(cl,indexPairLists)
	clusterExport(cl,"Global.curProfs")#assign the current profiles to all nodes
	clusterExport(cl,"dotProduct")#assign the current profiles to all nodes
	distResults<-clusterApply(cl,indexPairsGroups,pairDist2,type,isUnion)#calculate the distance between each pair of profiles
	distResults<-unlist(distResults,recursive=FALSE)
	distResults<-do.call(rbind,distResults)
	fillDistMatrix<-function(distResult)
	{
		
		i<-distResult[1]
		j<-distResult[2]
		v<-distResult[3]
		r[i,j]<<-v
		
	}
	
	res<-apply(distResults,1,fillDistMatrix)#fill the distance matrix with the results
	r <- r + t(r) - diag(diag(r)) #fill another half symetric values of matrix
	columnNames<-names(Global.curProfs)
	rownames(r) <- columnNames# assign feature ID as col and row names
	colnames(r) <- columnNames
	r  
}
parDistCal3<-function (c1,type="deconvolution",isUnion=TRUE) 
{
   	ncy <- ncx <- length(Global.curProfs)
	
   	if (ncx == 0) 
        stop("'x' is empty")
   	r <- matrix(0, nrow = ncx, ncol = ncy)
	indexPairs<-combinations(ncx,2)
	indexPairLists<-vector("list",nrow(indexPairs))
	for(i in 1:nrow(indexPairs))
	{
		indexPairLists[[i]]<-indexPairs[i,]
		}
	indexPairsGroups<-clusterSplit(c1,indexPairLists)
	clusterExport(c1,"Global.curProfs")#assign the current profiles to all nodes
	distResults<-clusterApply(c1,indexPairsGroups,pairDist2,type,isUnion)#calculate the distance between each pair of profiles
	distResults<-unlist(distResults,recursive=FALSE)
	distResults<-do.call(rbind,distResults)
	fillDistMatrix<-function(distResult)
	{
		
		i<-distResult[1]
		j<-distResult[2]
		v<-distResult[3]
		r[i,j]<<-v
	
	}

	apply(distResults,1,fillDistMatrix)#fill the distance matrix with the results
    r <- r + t(r) - diag(diag(r)) #fill another half symetric values of matrix
	getColName<-function(Global.curProfs)colnames(Global.curProfs)[2]
	columnNames<-do.call(rbind,lapply(Global.curProfs,getColName)) #get feature ID 

    rownames(r) <- columnNames# assign feature ID as col and row names
    colnames(r) <- columnNames
    r  
}
parDistCal2<-function (c1,type="deconvolution",isUnion=TRUE) 
{
   	ncy <- ncx <- length(Global.curProfs)
	
   	if (ncx == 0) 
        stop("'x' is empty")
   	r <- matrix(0, nrow = ncx, ncol = ncy)
	indexPairs<-combinations(ncx,2)
	indexPairLists<-vector("list",nrow(indexPairs))
	for(i in 1:nrow(indexPairs))
	{
		indexPairLists[[i]]<-indexPairs[i,]
		}
	indexPairsGroups<-clusterSplit(c1,indexPairLists)
	clusterExport(c1,"Global.curProfs")#assign the current profiles to all nodes
	distResults<-clusterApply(c1,indexPairsGroups,pairDist2,type,isUnion)#calculate the distance between each pair of profiles
	distResults<-unlist(distResults,recursive=FALSE)
	distResults<-do.call(rbind,distResults)
	fillDistMatrix<-function(distResult)
	{
		
		i<-distResult[1]
		j<-distResult[2]
		v<-distResult[3]
		
		if(is.na(v))
		{
			v<-0
			warning(paste(i,j,v,":dot product is NA!replaced by 0"))
			
			}
		r[i,j]<<-v
	}

	apply(distResults,1,fillDistMatrix)#fill the distance matrix with the results
    r <- r + t(r) - diag(diag(r)) #fill another half symetric values of matrix
	getColName<-function(Global.curProfs)colnames(Global.curProfs)[2]
	columnNames<-do.call(rbind,lapply(Global.curProfs,getColName)) #get feature ID 

    rownames(r) <- columnNames# assign feature ID as col and row names
    colnames(r) <- columnNames
    r  
}

pairDist5<-function(indexPairsGroup,goodShapePeaklist,type,isUnion=TRUE)
{
	
	
	if(length(indexPairsGroup)>0)
	{
		result<-vector("list",length(indexPairsGroup))
		for(ind in 1:length(indexPairsGroup))
		{
			i<-indexPairsGroup[[ind]][1]
			j<-indexPairsGroup[[ind]][2]
			
			ET.x<-(goodShapePeaklist[i,]$offset+goodShapePeaklist[i,]$lboundInd):(goodShapePeaklist[i,]$offset+goodShapePeaklist[i,]$rboundInd)
			x2 <- Global.vecInt[ET.x]
			ET.y<-(goodShapePeaklist[j,]$offset+goodShapePeaklist[j,]$lboundInd):(goodShapePeaklist[j,]$offset+goodShapePeaklist[j,]$rboundInd)					
			y2 <- Global.vecInt[ET.y]
			
			
			if(type=="Alignment")
			{
				x2.fileid<-strsplit(colnames(x2)[2],split="\\.")[[1]][1]
				y2.fileid<-strsplit(colnames(y2)[2],split="\\.")[[1]][1]
				if(x2.fileid==y2.fileid)
				{
					result[[ind]]<-c(i=i,j=j,dotProd=90)
					next#without merge the compound from the same file
				}
			}
			if(type=="deconvolution")
			{
				data1<-merge(x2,y2,by.x="ET",by.y="ET",all=isUnion)
			}else{
				data1<-merge(x2,y2,by.x="mz",by.y="mz",all=isUnion)
			}
			data1[is.na(data1)]<-0
			result[[ind]]<-c(i=i,j=j,dotProd=dotProduct(x=data1[,2], y=data1[,3]))
		}
		result
	}
}


pairDist2<-function(indexPairsGroup,type,isUnion=TRUE)
{


	if(length(indexPairsGroup)>0)
	{
		result<-vector("list",length(indexPairsGroup))
		for(ind in 1:length(indexPairsGroup))
		{
			i<-indexPairsGroup[[ind]][1]
			j<-indexPairsGroup[[ind]][2]
			x2 <- Global.curProfs[[i]]
		    y2 <- Global.curProfs[[j]]
			
			
			if(type=="Alignment")
			{
				x2.fileid<-strsplit(colnames(x2)[2],split="\\.")[[1]][1]
				y2.fileid<-strsplit(colnames(y2)[2],split="\\.")[[1]][1]
				if(x2.fileid==y2.fileid)
				{
					result[[ind]]<-c(i=i,j=j,dotProd=90)
					next#without merge the compound from the same file
				}
			}
			if(type=="deconvolution")
			{
				#x2$ET <-x2$ET- (which.max(x2$int)-which.max(y2$int))
				
				data1<-merge(x2,y2,by.x="ET",by.y="ET",all=isUnion)
				
				
#				x11()
#				plot(x2$int,type="l")
#				points(y2$int,type="l")
				
				
			}else{
				data1<-merge(x2,y2,by.x="mz",by.y="mz",all=isUnion)
			}
				data1[is.na(data1)]<-0
				result[[ind]]<-c(i=i,j=j,dotProd=dotProduct(x=data1[,2], y=data1[,3]))
		}
		result
	}
}

parDistCal<-function (c1) 
{
   	ncy <- ncx <- length(Global.curProfs)
	
   	if (ncx == 0) 
        stop("'x' is empty")
   	r <- matrix(0, nrow = ncx, ncol = ncy)
	indexPairs<-combinations(ncx,2)
	
	clusterExport(c1,"Global.curProfs")#assign the current profiles to all nodes
	distResults<-parApply(c1,indexPairs,1,pairDist)#calculate the distance between each pair of profiles
	
	fillDistMatrix<-function(distResult)
	{
		
		i<-distResult[1]
		j<-distResult[2]
		v<-distResult[3]
		r[i,j]<<-v
	
	}

	apply(distResults,2,fillDistMatrix)#fill the distance matrix with the results
    r <- r + t(r) - diag(diag(r)) #fill another half symetric values of matrix
	getColName<-function(Global.curProfs)colnames(Global.curProfs)[2]
	columnNames<-do.call(rbind,lapply(Global.curProfs,getColName)) #get feature ID 

    rownames(r) <- columnNames# assign feature ID as col and row names
    colnames(r) <- columnNames
    r  
}



pairDist<-function(indexPair)
{
	
	i<-indexPair[1]
	j<-indexPair[2]
	x2 <- Global.curProfs[i]
    y2 <- Global.curProfs[j]

	data1<-merge(x2,y2,by.x="ET",by.y="ET",all=TRUE)
	data1[is.na(data1)]<-0


	c(i=i,j=j,dotProd=dotProduct(data1[,2], data1[,3]))

}
matr2list<-function(mat)
{
	nlist<-nrow(mat)
	outlist<-vector("list",nlist)
	for(i in 1:nlist)outlist[[i]]<-mat[i,]
	outlist
	}
	
