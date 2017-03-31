library(ellipse)

#----- PCA-scores ----
# scal.data: data frame after scalling pretreatment
# a1,a2: pc number of score plot

# Make the PCA model
# scal.data: data from scalling, list format
# scal.data$X: save X data
# scal.data$Y: Y vector only
# output: a list of pca results

PCAModel <- function(scal.data)
{
	datafile <- scal.data$X
	#group <- scal.data$Y$"Y"+1
	group <- scal.data$Y+1
	pca<-prcomp(datafile)
	
	corr_matrix <- cor(datafile)
	myEig <- eigen(corr_matrix)   # eigenvalues % vectors
	variance <- 100*myEig$values/sum(myEig$values)
	pca$variance <- variance
	pca$Data <- datafile
	pca$group <- group
	pca
}

PCAScorePlot <- function(scal.data,a1,a2)
{
	pca <-PCAModel(scal.data)
	datafile <- pca$Data
	group <- pca$group
	scores_value <- pca$x  # scores
	variance <- pca$variance
	
	if (!missing(a2)) {
		pca.elps <- ellipse(cov(pca$x[,c(a1,a2)]),level = 0.95)
		
		plot(scores_value[,c(a1,a2)],main="PCA score plot",
				xlab=paste(paste("PC",a1,sep="")," (",format(variance[a1],digit=4),"%)",sep=""),
				ylab=paste(paste("PC",a2,sep="")," (",format(variance[a2],digit=4),"%)",sep=""),
				pch=16,cex=2,col=group,xlim=range(pca.elps[,1])*1.15,ylim=range(pca.elps[,2])*1.5)
		points(pca.elps,type="l",col="grey")
		abline(h=0,col="grey")
		abline(v=0,col="grey")
		text(pca$x[,a1]+2,pca$x[,a2],labels=rownames(datafile),cex=1)
	}
	if (missing(a2)) 
	{plot(scores_value[,a1],xlab="sample number",
				ylab = paste(paste("PC",a1,sep="")," (",format(var_pro[a1]*100,digit=4),"%)",sep=""),
				pch=16,cex=2,col=group)}
}

PCA3DPlot <- function(scal.data,a1,a2,a3) {
	pca <-PCAModel(scal.data)
	scores_value <- pca$x  # scores
	group <- pca$group
	
	library(rgl)  # checking 3D-PCA
	open3d()
	plot3d(scores_value[,c(a1,a2,a3)],type="s",size=3,main="3D-PCA",col=group)
	aspect3d(1,1,1)
}

PCALoadingPlot <- function(scal.data,a1,a2)
{
	pca <-PCAModel(scal.data)
	datafile <- pca$Data
	group <- pca$group
	
	pca<-prcomp(datafile)
	loadings_value <- pca$rotation # loadings
	
	if (!missing(a2)) {	
		plot(loadings_value[,c(a1,a2)],main="PCA loading plot",xlab=paste("P",a1,sep=""),
				ylab=paste("P",a2,sep=""), pch=3, col="blue")
		abline(h=0,col="grey")
		abline(v=0,col="grey")
		text(loadings_value[,c(a1,a2)],labels=colnames(datafile),cex=0.6,pos=1)
	}
	if (missing(a2)) {
		plot(loadings_value[,a1],xlab="variable number",ylab = paste("P",a1,sep=""), pch=3, col="blue")}
		abline(h=0,col="grey")
		text(loadings_value[,a1],labels=colnames(datafile),cex=0.6,pos=1)
}

pca_expvarcum <- function(scal.data)  
{
	pca <-PCAModel(scal.data)
	datafile <- pca$Data
	# outputing cummulative variance of n PCs
	Nrow <- nrow(datafile)
	Ncol <- ncol(datafile)
	n <- min(Nrow,Ncol)
	
	var_pro <- pca$variance/100	
	var_pro_cum <- var_pro
	for (i in 2:length(var_pro_cum)) {
		var_pro_cum[i] <- sum(var_pro_cum[i-1],var_pro_cum[i])
	}
	
	plot(var_pro_cum[1:n],main="PCA: Scree plot",xlab="Principal components",ylab="Proportion of variance",
			ylim=c(0,1),type="o",pch=19,col="red")
	lines(var_pro[1:n],col="blue",type="o",pch=19)
	legend(1,c("R2X(cum)","R2X"),cex=0.8,col=c("red","blue"),pch=19)
	abline(h=0.5,col="red") 
}

# checking single/cumulative explaining variance distribution
pca_expvarboth <- function(scal.data,p)
{
	pca <-PCAModel(scal.data)
	datafile <- pca$Data	
	var_pro <- pca$variance/100	
	var_pro_cum <- var_pro
	for (i in 2:length(var_pro_cum)) {
		var_pro_cum[i] <- sum(var_pro_cum[i-1],var_pro_cum[i])
	}
	#min(which(var_pro_cum>=1))	
	data_var <- rbind(var_pro[1:p],var_pro_cum[1:p])
	barplot(data_var,names.arg=c(1:p),xlab="Principal component number",beside=TRUE,col=c("green","blue"),ylim=c(0,1),main="PCA: explaining variance")
	legend("topleft",c("R2X","R2X(cum)"),fill=c("green","blue"))
}


PCAModelSummary<- function(scal.data,p) 
{ 
	library(pcaMethods) #pca,Q2 need pcaMethods library
	# default method is "ppca": probabilistic PCA
	# p: The amount of principal components to estimate Q2 for
	# fold: default no is 7
	datafile <- scal.data$X
	# name conflicts, to use :: avoiding the problem
	resPPCA <- pcaMethods::pca(datafile, method = "ppca", center = TRUE, nPcs = p)
	q2PPCA <- Q2(resPPCA,datafile,nruncv = 1,fold = 7)
	
	# TO combine R2x, Q2 variable together
	data_var <- rbind(resPPCA@R2,resPPCA@R2cum,t(q2PPCA))
	max_var <- max(data_var)
	min_var <- min(data_var)
	barplot(as.matrix(data_var),xlab="Principal component number",ylab="Variance",ylim=c(min_var-0.1,max_var+0.1),beside=TRUE,col=c("yellow","green","blue"))
	legend("topleft",c("R2X","R2X(cum)","Q2X(cum)"),fill=c("yellow","green","blue"))
	abline(h=0)
}

R2VX_sin <- function(scal.data,f) { # f: PC component number
	pca <- PCAModel(scal.data)	
	datafile <- pca$Data		
	scores_value <- round(pca$x,2)
	loadings_value <- round(pca$rotation,2)	
	loadings_value <- t(loadings_value)  # Be careful about matrix size
	if (f==1) {
		error_value <- datafile - scores_value[,1:f] %*% t(as.matrix(loadings_value[1:f,]))}
	if (f!=1) {
		error_value <- datafile - scores_value[,1:f] %*% loadings_value[1:f,]}	
	R2V_list <- c()	
	for (j in 1:ncol(error_value)) {
		S2VXj <- sum(error_value[,j]^2)/nrow(error_value)
		S2VXj0 <- sum(datafile[,j]^2)/nrow(datafile)
		R2V_list <- c(R2V_list,1-S2VXj/S2VXj0)
	}
	return(R2V_list)
}

R2VX_showmore <- function(scal.data,f,workdir) { 	
	datafile <- scal.data$X		
	f <- as.integer(f)
	data_matrix <- NULL
	
	if (f==1) {
		data_matrix <- rbind(data_matrix,R2VX_sin(scal.data,f))
		colnames(data_matrix) <- colnames(datafile)		
		data_matrix  <- data_matrix[,which(data_matrix>= 0.7)]	
		name_list <- colnames(data_matrix)
		colnames(data_matrix) <- NULL	
		barplot(data_matrix,beside=F,xlab="Variable Num",ylab="Variance", ylim=c(0,1),col=rainbow(f)[f],main="Well Explained variables")
		legend("top",ncol=f,"PC1",fill=rainbow(f))
	}
	
	else {			
		for (i in 1:f) {
			data_matrix <- rbind(data_matrix,R2VX_sin(scal.data,i))
		}		
		colnames(data_matrix) <- colnames(datafile)
		# the cummulative explained variance > 0.7
		data_matrix  <- data_matrix[,which(data_matrix[f,]>= 0.7)]	
		name_list <- colnames(data_matrix)
		
		# output csv file
		#work.dir <- getwd()
		tempDir <- paste(workdir,"output/PCA/WellExpVar/",sep="")
		dir.create(tempDir)
		write.csv(t(data_matrix),paste(tempDir,"well.explained.variables_PCA.csv",sep=""))
		
		# bar plots of individual component variance
		colnames(data_matrix) <- NULL	
		position <-barplot(data_matrix[f,],beside=F,xlab="Variable Num",ylab="Variance", ylim=c(0,1.2),col=rainbow(f)[f],main="Well Explained variables")
		axis(1,at=position,labels=name_list)	
		
		# text(test,par("usr")[3] - 0.2,labels=name_list,srt=90,pos=3,offset=0)	
		PCs <- rev(c(1:(f-1)))
		for (i in PCs)
		{
			barplot(data_matrix[i,],add=T,col=rainbow(f)[i],main="Well Explained variables")}	
		#barplot(data_matrix,beside=F,col=rainbow(f),main="Well Explained variables")	
		legend.name <- NULL
		for (i in 1:f) {
			legend.name <- c(legend.name,paste("PC",i,sep=""))
		}
		legend("top",ncol=f,legend.name,fill=rainbow(f))
	}
}

R2VXadj <- function(scal.data,p) { 
	pca <- PCAModel(scal.data)
	datafile <- pca$Data		
	scores_value <- round(pca$x,2)
	loadings_value <- round(pca$rotation,2)	
	loadings_value <- t(loadings_value)  # Be careful about matrix size
	#error_value <- datafile - scores_value[,1:f] %*% loadings_value[1:f,]	
	if (p==1) {
		error_value <- datafile - scores_value[,1:p] %*% t(as.matrix(loadings_value[1:p,]))}
	if (p!=1) {
		error_value <- datafile - scores_value[,1:p] %*% loadings_value[1:p,]}	
	R2V_list <- c()	
	for (j in 1:ncol(error_value)) {
		S2VXj <- sum(error_value[,j]^2)/nrow(error_value)
		S2VXj0 <- sum(datafile[,j]^2)/nrow(datafile)
		R2V_list <- c(R2V_list,1-S2VXj/S2VXj0)
	}
	
	names(R2V_list) <- colnames(datafile)
	
	plot(R2V_list,type="h",xlab="variable Num",ylab="R2VX(adj,cum)",main=paste("PCA overview (R2VX(adj,cum), ",p,"PCs)"))	
	if (length(which(R2V_list>0.7))>=1) {
		highvalue_list <- R2V_list[which(R2V_list>0.7)]	
		text(which(R2V_list>0.7),highvalue_list,pos=4,col="red",offset =0.5,labels=names(highvalue_list))}
}


R2VXadj_plot <- function(scal.data,p){	
	# the P value should be an integer
	par(mfrow=c(as.integer(p),1))	
	for (i in 1:p) {	
		R2VXadj(scal.data,i)}}


## ---------------------------------------------------------------------------	

Dmodx <- function(scal.data,f) {
	pca <- PCAModel(scal.data)
	datafile <- pca$Data	
	group <- as.numeric(pca$group)	
	scores_value <- round(pca$x,2)
	loadings_value <- round(pca$rotation,2)
	loadings_value <- t(loadings_value)  # Be careful about matrix size
	if (f==1) {                                                                                                                
		error_value <- datafile - scores_value[,1:f] %*% t(as.matrix(loadings_value[1:f,]))}
	if (f!=1) {
		error_value <- datafile - scores_value[,1:f] %*% loadings_value[1:f,]}
	error_list <- c()
	for (i in 1:nrow(error_value)) {
		variance <- var(as.numeric(error_value[i,]))
		error_list <- c(error_list,sqrt(sum(error_value[i,]^2)/(ncol(error_value)*variance)))
	}		
	plot(error_list,col=group,pch=group,xlab="Num",ylab=paste("DModx[",f,"](norm)",sep=""),main="Distance to the model")
	text(error_list,pos=3,offset = 1,labels=rownames(error_value))
	lines(error_list)	
}

# Providing Pearson relationship between group and selected score values 
# t-test shows if score is significant different between groups

pca_corr <- function(scal.data,a1)
{
	pca <- PCAModel(scal.data)
	datafile <- pca$Data	
	group <- as.numeric(pca$group)
	scores_value <- pca$x
	corr <- cor(group,scores_value[,a1])
	ttest <- pairwise.t.test(scores_value[,a1],group)
	#pvalue <- as.numeric(ttest$p.value)
	pvalue <- round(ttest$p.value,2)
	boxplot(scores_value[,a1]~group,main=paste("Group VS. Score[",a1,"]",sep=""),xlab="group",
			ylab=paste("score [",a1,"]",sep=""),cex.axis=1.5,cex.names=3,col="light blue")
	legend("bottomright",c(paste("Pearson corr: ",round(corr,2)),paste("P-ttest:",round(pvalue,6))))
}


#-----------Pearson regression relationship between selected scores and variables-----------

# calculating Pearson correlation coefficients between selected scores and variables
# A: seleted PC number

pca_regression <- function(scal.data,A)
{
	pca <- PCAModel(scal.data)
	datafile <- pca$Data	
	group <- as.numeric(pca$group)
	corr_matrix <- c(cor(pca$x[,A],datafile))
	
	names(corr_matrix) <- colnames(datafile)
	barplot(corr_matrix,col="green",ylab=paste("Coefficients [",A,"]",sep=""),
			ylim=c(-1,1),main="Pearson: Scores vs. Variables")
	abline(h=0.5,col="red")
	abline(h=-0.5,col="red")
	#highvalue_list <- corr_matrix[which(abs(corr_matrix)>0.5)]
	#text(highvalue_list,pos=4,col="red",offset =0.5,labels=names(highvalue_list))
}	

pca <- function(params,scal.data)
{
	workdir <- params$WorkDir
	f <- as.integer(params$PCA_PCs)
	a1 <- as.integer(params$PCA_CID)
			
			
	#setwd(workdir)
	#work.dir <- getwd()
	create.dir <- paste(workdir,"output/PCA/",sep="")
	dir.create(create.dir)
	
	pca <- PCAModel(scal.data)
	write.csv(pca$x,paste(create.dir,"pca-scores.csv",sep=""))
	write.csv(pca$rotation,paste(create.dir,"pca-loadings.csv",sep="/"))
	
	# f: # of PCs to show
#	pdf(paste(create.dir,"PCA.pdf",sep=""))
#	rown <- as.integer(f-1)
#	par(mfrow=c(rown,2)) 
#	for (i in c(2:f)) {
#		PCAScorePlot(scal.data,1,i)
#		PCALoadingPlot(scal.data,1,i)
#		
#	}
#	dev.off()
	
	pdf(paste(create.dir,"PCA-scores.pdf",sep=""))
	for (i in c(2:f)) {
		PCAScorePlot(scal.data,1,i)
	}
	dev.off()
	
	pdf(paste(create.dir,"PCA-loadings.pdf",sep=""))
	for (i in c(2:f)) {
		PCALoadingPlot(scal.data,1,i)
	}
	dev.off()
	
	pdf(paste(create.dir,"PCA-ModelSummary.pdf",sep=""))
	par(mfrow=c(2,1)) 
	pca_expvarcum(scal.data)
	pca_expvarboth(scal.data,f)	
	dev.off()

	pdf(paste(create.dir,"PCA-VariableVariance.pdf",sep=""))
	R2VX_showmore(scal.data,f,workdir)
	dev.off()
	
	pdf(paste(create.dir,"PCA-DistanceToModel.pdf",sep=""))
	par(mfrow=c(f,1)) 
	for (i in c(1:f)) {
	Dmodx(scal.data,i)
	}
	dev.off()
	
	pdf(paste(create.dir,"PCA-scoreVsgroup.pdf",sep=""))
	pca_corr (scal.data,a1)
	dev.off()	
	
	print("PCA Done!")
}

############################################################################################
############################################################################################

################### 
#----- PLS -------- 
###################

library(pls)

## parameters ----------------------
## A: principal component number
## a1,a2: ploting axis
## choice: cross validation choice
## "v" variable name, "y" group information
## parameters ----------------------

pls.model <- function(scal.data,choice,A) {
	# creating a data frame	
	#data_com <- data.frame(y=scal.data$Y$"Y")
	data_com <- data.frame(y=scal.data$Y+1)
	data_com$v <- as.matrix(scal.data$X)
	# Leave-one-out
	if (choice=="LOO") 
		pls_fit<- plsr(y~v, ncomp = A, data = data_com, validation = choice)
	# the default is 7 fold cross validation
	else if (choice=="CV") 
		pls_fit<- plsr(y~v, ncomp = A, data = data_com, validation = choice,segments=7)
	return(pls_fit) 
}

#----- PLS-scores ----------------------------------
library(ellipse)

# a2 could be missing
pls_score.x <- function(scal.data,choice,A,a1,a2) {  
	pls_data <- pls.model(scal.data,choice,A)
	group <- pls_data$model$y
	
	scores <- pls_data$scores	
	### Ploting ellipse in PLS plot
	if (!missing(a2)) {
		pls.elps<- ellipse(cov(scores[,c(a1,a2)]),level = 0.95)
		
		plot(scores[,c(a1,a2)],xlab=paste("t [",a1,"]",sep=""),
				ylab=paste("t [",a2,"]",sep=""), main="PLS score plot (X)",pch=16,cex=2,
				col=group,xlim=range(pls.elps[,1])*1.15,ylim=range(pls.elps[,2])*1.15)
		points(pls.elps,type="l",col="grey")
		abline(h=0,col="grey")
		abline(v=0,col="grey")
		### labeling choice
		text(scores[,a1]+2,scores[,a2],labels=rownames(scal.data$X),cex=0.8)
	}
	if (missing(a2)) {
		plot(scores[,a1],main="PLS score plot (X)",xlab="sample number",ylab = paste("t [",a1,"]",sep=""),pch=16,cex=2,col=group)
		text(scores[,a1]+0.5,labels=rownames(scal.data$X),cex=0.8)
		abline(h=0,col="green")}
}

#pls_score.x(s_data,15,10,1,2,"LOO",10)

pls_score.y <- function(scal.data,choice,A,a1,a2) {  
	pls_data <- pls.model(scal.data,choice,A)
	group <- pls_data$model$y
	Yscores <- pls_data$Yscores
	
	if (!missing(a2)) {
		pls.elps<- ellipse(cov(Yscores[,c(a1,a2)]),level = 0.95)
		
		plot(Yscores[,c(a1,a2)], xlab=paste("u [",a1,"]",sep=""),
				ylab=paste("u [",a2,"]",sep=""),main="PLS score plot (Y)",pch=16,cex=2,col=group,xlim=range(pls.elps[,1])*1.15,ylim=range(pls.elps[,2])*1.15)
		points(pls.elps,type="l",col="grey")
		abline(h=0,col="grey")
		abline(v=0,col="grey")
		text(Yscores[,c(a1,a2)],labels=rownames(scal.data$X),cex=0.8)
	}
	if (missing(a2)) {
		plot(Yscores[,a1],xlab="sample number",ylab = paste("u [",a1,"]",sep=""),pch=16,cex=2,col=group)
		text(Yscores[,a1],labels=rownames(scal.data$X),cex=0.8,pos=3)
	}
		
}

pls_score.xy <- function(scal.data,choice,A,a1,a2) {  
	pls_data <- pls.model(scal.data,choice,A)
	group <- pls_data$model$y
	
	scores <- pls_data$scores	
	Yscores <- pls_data$Yscores
	
	if (!missing(a2)) {
		com_matrix <- cbind(scores[,a1],Yscores[,a2])
		pls.elps<- ellipse(cov(com_matrix[,c(1,2)]),level = 0.95)
		
		plot(scores[,a1],Yscores[,a2],xlab=paste("t [",a1,"]",sep=""),
				ylab=paste("u [",a2,"]",sep=""),pch=16,cex=2,main="PLS score plot (X,Y)",col=group,xlim=range(pls.elps[,1])*1.15,ylim=range(pls.elps[,2])*1.15)
		points(pls.elps,type="l",col="grey")
		abline(h=0,col="grey")
		abline(v=0,col="grey")
		text(scores[,a1],Yscores[,a2],labels=rownames(scal.data$X),cex=0.8)
	}
	if (missing(a2)) 
	{
		plot(scores[,a1],xlab="sample number",ylab = paste("t [",a1,"]",sep=""),pch=16,col=group)
		text(scores[,a1],labels=rownames(scal.data$X),cex=0.8,pos=3)
		abline(h=0,col="grey")
	}
}



pls_3D <- function(scal.data,choice,A,a1,a2,a3) { 
	pls_data <- pls.model(scal.data,choice,A)
	group <- pls_data$model$y+1
	library(rgl)  # checking 3D-PCA
	open3d()
	plot3d(pls_data$scores[,c(a1,a2,a3)],type="s",size=3,main="PLS-3D scores",col=group)
	aspect3d(1,1,1)
}

pls_3D_loading <- function(scal.data,choice,A,a1,a2,a3) { 
	pls_data <- pls.model(scal.data,choice,A)
	group <- pls_data$model$y+1
	library(rgl)  # checking 3D-PCA
	open3d()
	plot3d(pls_data$loadings[,c(a1,a2,a3)],main="PLS-3D loadings",col="blue")
	aspect3d(1,1,1)
}	



#-----  PLS-loadings ----------------------------------

pls_loading <- function(scal.data,choice,A,a1,a2) {  
	pls_data <- pls.model(scal.data,choice,A)
	if (!missing(a2)) {
		
		plot(pls_data$loadings[,c(a1,a2)], main="PLS loading plot",
				xlab=paste("PC",a1,sep=""),ylab=paste("PC",a2,sep=""), pch=3, col="blue")
		abline(h=0,col="grey")
		abline(v=0,col="grey")
		text(pls_data$loadings[,c(a1,a2)],labels=colnames(scal.data$X),cex=0.6,pos=1)
	}
	if (missing(a2)) 
		plot(pls_data$loadings[,a1],main="PLS loading plot",xlab="variable number",ylab = paste("p [",a1,"]",sep=""),pch=3, col="blue")
		abline(h=0,col="grey")
		text(pls_data$loadings[,a1],labels=colnames(scal.data$X),cex=0.6,pos=1)
}	

#----- Top n PCs-----------------------------------------
#--------------------------------------------------------

#----- PLS-scores -----
pls_score.top <- function(scal.data,choice,A,f) { 
	#f :  # of PCs to show scores
	pls_data <- pls.model(scal.data,choice,A)
	group <- pls_data$model$y
	
	plot(pls_data, plottype = "scores",pch=16,cex=2,comps = 1:f,col=group)
}

#----- PLS-loadings----	
pls_loading.top <- function(scal.data,choice,A,f) {  
	pls_data <- pls.model(scal.data,choice,A)
	group <- pls_data$model$y
	plot(pls_data, plottype="loadings", comps = 1:f,legendpos = "topleft", xlab = "variable")
	abline(h = 0)
}


#----- PLS-evaluation 1 -------------------------------------	
## checking optimal number of PCs
pls_validation <- function(scal.data,choice,A) { 
	pls_data <- pls.model(scal.data,choice,A)
	
	plot(pls_data,legendpos= "topright",plottype="validation")
	#or: plot(RMSEP(pls_data))
}

pls_varcum <- function(scal.data,choice,A) { 	
	pls_data <- pls.model(scal.data,choice,A)
	y <- as.numeric(pls_data$model$y)
	datafile <- pls_data$model$v
	
	# explained variance for x	
	var_pro <- explvar(pls_data)
	var_pro_cum <- var_pro
	for (i in 2:length(var_pro_cum)) {
		var_pro_cum[i] <- sum(var_pro_cum[i-1],var_pro_cum[i])
	}
	
	# explained variance for Y, cv
	yve <- 100 * drop(R2(pls_data, estimate = "CV", 
					intercept = FALSE)$val)
	
	# predicted variance for Y, cv, self definition	
	Q2.value <- 100*pls.Q2(datafile,y,nc=A,cv=TRUE)
	Q2 <- Q2.value[1,]
	Q2cum <- Q2.value[2,]	
	data_var <- rbind(var_pro,var_pro_cum,yve,Q2,Q2cum)	
	barplot(data_var,names.arg=c(1:A),xlab="Principal component number",ylab="explained variance (%)",
			beside=TRUE,col=c("yellow","green","blue","grey","red"),main="PLS: explaining variance",
			ylim=c(-50,100))
	
	#barplot(explvar(pls_data)) # explaining variance
	legend("topleft",c("R2X","R2X(cum)","R2Y(cum)","Q2","Q2Y(cum)"),
			fill=c("yellow","green","blue","grey","red"))
}	

#----- PLS-DModX,Y---------------------------	
## f: PC number	
## m,n: group number
## choice: cross validation choice
#
#PLS_Dmodx <- function(scal.data,f,choice,A) {
#	pls_data <- pls.model(scal.data,choice,A)
#	group <- pls_data$model$y+1
#		
#	scores_value <- round(pls_data$scores,2)
#	loadings_value <- round(pls_data$loadings,2)
# Be careful about matrix size
#	loadings_value <- t(loadings_value)  
#	if (f==1) {
#		error_value <- scal.data$v - scores_value[,1:f] %*% t(as.matrix(loadings_value[1:f,]))}
#	if (f!=1) {
#		error_value <- scal.data$v - scores_value[,1:f] %*% loadings_value[1:f,]}
#	error_list <- c()
#	for (i in 1:nrow(error_value)) {
#		variance <- var(as.numeric(error_value[i,]))
#		error_list <- c(error_list,sqrt(sum(error_value[i,]^2)/(ncol(error_value)*variance)))
#		}			
#	plot(error_list,col=group,pch=group,xlab="Num",ylab=paste("DModx[",f,"](norm)",sep=""),main="Distance to the model (X)")
#	text(error_list,pos=3,offset = 1,labels=rownames(error_value))
#	lines(error_list)	
#	}
#

PLS_Dmodx <- function(scal.data,f,choice,A) {
	# A: # of components to make pls model
	# f: the component ID to show distance
	
	pls_data <- pls.model(scal.data,choice,A)
	group <- pls_data$model$y
	datafile <- pls_data$model$v
	scores_value <- round(pls_data$scores,2)
	loadings_value <- round(pls_data$loadings,2)
	# Be careful about matrix size
	loadings_value <- t(loadings_value)  
	if (f==1) {
		error_value <- datafile - scores_value[,1:f] %*% t(as.matrix(loadings_value[1:f,]))}
	if (f!=1) {
		error_value <- datafile - scores_value[,1:f] %*% loadings_value[1:f,]}
	error_list <- c()
	# two parameters
	a <- ncol(datafile)-A
	b <- nrow(datafile)-A-1
	# So calculation
	C <- sqrt(sum(rowSums(error_value^2))/(a*b))
	
	# Distance calculation : absolute distance to model (normalization)
	for (i in 1:nrow(error_value)) {
		#variance <- var(as.numeric(error_value[i,]))
		error_list <- c(error_list,sqrt(sum(error_value[i,]^2)/a)/C)
	}			
	plot(error_list,col=group,pch=16,xlab="Num",ylab=paste("DModx[",f,"](norm)",sep=""),main="Distance to the model (X)",ylim=c(0,1.2))
	text(error_list,pos=3,offset = 1,labels=rownames(error_value))
	lines(error_list)	
}

#----------------------------------------------

#PLS_DmodY <- function(scal.data,f,choice,A) {
#	# "v" variable name, "y" group information
#	# group info has 1 vector does not work actually
#	# pls_data$loadings = zero
#	
#	pls_data <- pls.model(scal.data,choice,A)
#	group <- pls_data$model$y	
#	datafile <- pls_data$model$v
#	scores_value <- round(pls_data$Yscores,2)
#	loadings_value <- round(pls_data$Yloadings,2)
#	loadings_value <- t(loadings_value)  # Be careful about matrix size
#	if (f==1) {
#		error_value <- datafile - scores_value[,1:f] %*% t(as.matrix(loadings_value[1:f,]))
#	}
#	if (f!=1) {
#		error_value <- datafile - scores_value[,1:f] %*% loadings_value[1:f,]
#	}
#	error_list <- c()
#	for (i in 1:nrow(error_value)) {
#		variance <- var(as.numeric(error_value[i,]))
#		error_list <- c(error_list,sqrt(sum(error_value[i,]^2)/(ncol(error_value)*variance)))
#	}	
#	
#	plot(error_list,col=group,pch=16,xlab="Num",ylab=paste("DModY[",f,"](norm)",sep=""),main="Distance to the model (Y)")
#	text(error_list,pos=3,offset = 1,labels=rownames(error_value))
#	lines(error_list)	
#}

## ~~~~~~~PLS_R2VX variable explaining variance~~~~~~~~~~~~------------------

PLS_R2VXadj <- function(scal.data,f,choice,A) { 
	pls_data <- pls.model(scal.data,choice,A)
	group <- pls_data$model$y
	datafile <- pls_data$model$v
	scores_value <- round(pls_data$scores,2)
	loadings_value <- round(pls_data$loadings,2)
	loadings_value <- t(loadings_value)  # Be careful about matrix size
	#error_value <- datafile - scores_value[,1:f] %*% loadings_value[1:f,]	
	if (f==1) {
		error_value <- datafile - scores_value[,1:f] %*% t(as.matrix(loadings_value[1:f,]))}
	if (f!=1) {
		error_value <- datafile - scores_value[,1:f] %*% loadings_value[1:f,]}
	error_list <- c()	
	for (j in 1:ncol(error_value)) {
		S2VXj <- sum(error_value[,j]^2)/nrow(error_value)
		S2VXj0 <- sum(datafile[,j]^2)/nrow(datafile)
		error_list <- c(error_list,1-S2VXj/S2VXj0)
	}
	names(error_list) <- colnames(datafile)	
	plot(error_list,type="h",xlab="variable Num",ylab="R2VX(adj,cum)",main=paste("PLS(R2VX(adj,cum) [PC",f,"]",sep=""))	
	if (length(which(error_list>0.7))>=1) {
		highvalue_list <- error_list[which(error_list>0.7)]	
		text(which(error_list>0.7),highvalue_list,pos=4,col="red",offset =0.5,labels=names(highvalue_list))}
}

PLS_R2VXadj_plot <- function(scal.data,f,choice,A){
	
	par(mfrow=c(as.integer(f),1))
	for (i in 1:f) {
		PLS_R2VXadj(scal.data,i,choice,A)}}

PLS.R2VX_sin <- function(scal.data,f,choice,A){ # f: PC component number
	pls_data <- pls.model(scal.data,choice,A)
	datafile <- pls_data$model$v
	group <- pls_data$model$y	
	
	scores_value <- round(pls_data$scores,2)
	loadings_value <- round(pls_data$loadings,2)
	loadings_value <- t(loadings_value) 
	if (f==1) {
		error_value <- datafile - scores_value[,1:f] %*% t(as.matrix(loadings_value[1:f,]))}
	if (f!=1) {
		error_value <- datafile - scores_value[,1:f] %*% loadings_value[1:f,]}	
	R2V_list <- c()	
	for (j in 1:ncol(error_value)) {
		S2VXj <- sum(error_value[,j]^2)/nrow(error_value)
		S2VXj0 <- sum(datafile[,j]^2)/nrow(datafile)
		R2V_list <- c(R2V_list,1-S2VXj/S2VXj0)
	}
	return(R2V_list)
}


PLS.R2VX.showmore <- function(scal.data,f,choice,A) { 	
	# A: PC number to make PLS model
	# f: PC ID
	datafile <- scal.data$X
	f <- as.integer(f)
	data_matrix <- NULL
	
	if (f==1) {
		data_matrix <- rbind(data_matrix,PLS.R2VX_sin(scal.data,f,choice,A))
		colnames(data_matrix) <- colnames(datafile)
		data_matrix  <- data_matrix[,which(data_matrix[f,]>= 0.4)]	
		name_list <- colnames(data_matrix)
		colnames(data_matrix) <- NULL	
		barplot(data_matrix,beside=F,xlab="Variable Num",ylab="Variance", ylim=c(0,1),col=rainbow(f)[f],main="Well Explained variables")
		legend("top",ncol=f,"PC1",fill=rainbow(f))
	}
	
	# f >= 2
	else {
		for (i in 1:f) {
			data_matrix <- rbind(data_matrix,PLS.R2VX_sin(scal.data,i,choice,A))
		}
		
		colnames(data_matrix) <- colnames(datafile)
		data_matrix  <- data_matrix[,which(data_matrix[f,]>= 0.4)]	
		name_list <- colnames(data_matrix)
		
		#  create saving dir
		work.dir <- getwd()
		create.dir <- paste(work.dir,"PLS/WellExpVar/",sep="/")
		dir.create(create.dir)
		write.csv(t(data_matrix),paste(create.dir,"well.explained.variables_PLS.csv",sep=""))
		
		colnames(data_matrix) <- NULL	
		position <-barplot(data_matrix[f,],beside=F,xlab="Variable Num",ylab="Variance", ylim=c(0,1),col=rainbow(f)[f],main="Well Explained variables")
		axis(1,at=position,labels=name_list)	
		
		PCs <- rev(c(1:(f-1)))
		for (i in PCs)
		{
			barplot(data_matrix[i,],add=T,col=rainbow(f)[i],main="Well Explained variables")}
		
		#barplot(data_matrix,col=rainbow(f),main="Well Explained variables")
		
		legend.name <- NULL
		for (i in 1:f) {
			legend.name <- c(legend.name,paste("PC",i,sep=""))
		}
		legend("topright",ncol=f,legend.name,fill=rainbow(f))
	}	
}

#----- PLS-evaluation 2 -------------------------------------	
pls_coef <- function(scal.data,choice,A,f) { 
	pls_data <- pls.model(scal.data,choice,A)
	
	coefplot(pls_data, legendpos="bottomright",ncomp = 1:f)
}	

pls_corr <- function(scal.data,choice,A,f) {
	pls_data <- pls.model(scal.data,choice,A)
	
	corrplot(pls_data,legendpos="topright",comps=1:f,col="blue",pch=4)
}	



pls_prediction <- function(scal.data,choice,f) { 
	# A: PCs to model and show prediction data	
	pls_data <- pls.model(scal.data,choice,f)
	group <- pls_data$model$y
	#plot(pls_data,plottype="prediction",col=group)
	
	par(mfrow=c(1,2))
	plot(pls_data, ncomp = f, asp = 1, line = TRUE,col=group) # predicted vs measured, asp=y/x axis ratio
	ttest <- pairwise.t.test(pls_data$validation$pred[,,f],group)
	#pvalue <- as.numeric(ttest$p.value)
	pvalue <- round(ttest$p.value,2)
	boxplot(pls_data$validation$pred[,,f]~group,col="lightblue",ylab="predicted scores")
	legend("bottomright",paste("P-ttest:",round(pvalue,2)))
	# Or: plot(y,pls_data$validation$pred[,,2],col=group)
}	


#pls.feature.select <- function(scal.data,choice,A,comp,n1)
#{
# A: component num to establish pls model
# comp: the pc number to show variables
# n1: top n1 number of variables
#	
#	pls_data <- pls.model(scal.data,choice,A)
#	
# order equals sort funtion, finding top n1 variables
#	fstore<- order(-abs(pls_data$coeff[,,comp]))[1:n1]
#	value<-  abs(pls_data$coeff[,,comp])[fstore]
#names(value)<- fstore
#	
#    bplot<- barplot(value,main = paste("PLS coeffcients (Top ",n1," variables) @ ","PC",comp,sep=""),
#		col="green",ylim=c(0,max(value)*1.03))
#text(bplot,value*1.02,labels=fstore)
#}
library(mixOmics)

pls.feature.select <- function(workdir,scal.data,choice,A,comp)
{
	# A: component num to establish pls model
	# comp: the pc number to show variables
	# n1: top n1 number of variables
	
	pls_data <- pls.model(scal.data,choice,A)
	datafile <- pls_data$model$v
	y <- pls_data$model$y
	# regression coefficients: B
	pls_reg_coef <- pls_data$coeff[,,comp]
	# pearson relation: variable vs Y
	cor_y <- cor(datafile,y)
	# pearson relation: variable vs score
	cor_score <- cor(datafile,pls_data$scores[,comp])
	# pearson covarance: variable vs score
	cov_score <- cov(datafile,pls_data$scores[,comp])
	# loadings
	loading <- pls_data$loadings[,comp]
	# weights
	loading.weights <- pls_data$loading.weights[,comp]
	# VIP
	x <- as.matrix(datafile,nrow=nrow(datafile))
	y <- as.factor(y)	
	pls <- plsda(x,y,ncomp=A)
	vip_value <- vip(pls)[,comp]
	
	#  create saving dir
	#work.dir <- getwd()
	create.dir <- paste(workdir,"output/PLS/MultiVarEval",sep="/")
	dir.create(create.dir)
	
	pdf(paste(create.dir,"/VarPlot1.pdf",sep=""))
	#x11()
	par(mfrow=c(3,1))
	plot(pls_reg_coef,type="h",col="red",xlab="Sample Num",ylab="regression coef")
	plot(loading,type="h",xlab="Sample Num",col="purple")
	plot(loading.weights,type="h",xlab="Sample Num")
	dev.off()
	
	pdf(paste(create.dir,"/VarPlot2.pdf",sep=""))
	#pdf("VarPlot2.pdf")
	#x11()
	par(mfrow=c(4,1))
	plot(cor_y,type="h",col="blue",xlab="Sample Num",ylab="coef_Y")
	plot(cor_score,type="h",col="green",xlab="Sample Num",ylab=paste("corr_Score[",comp,"]"))
	plot(cov_score,type="h",col="green",xlab="Sample Num",ylab=paste("cov_Score[",comp,"]"))
	plot(vip_value,type="h",xlab="Sample Num",ylab=paste("VIP[",comp,"]",sep=""))	
	dev.off()
	
	data_list <- rbind(pls_reg_coef,t(cor_y),t(cor_score),t(cov_score),t(loading),loading.weights,t(vip_value)) 
	data_list <- t(data_list)
	colnames(data_list) <- c("coef","cor(y)","cor(score)","cov(score)","loading","weights",paste("VIP",comp,sep=""))	
	write.csv(data_list,paste(create.dir,"/var_evaluation.csv",sep=""))
	print ("Output is done!")
}

pls.VarPlot <- function(scal.data,choice,A,comp,type)
{
	# A: component num to establish pls model
	# comp: the pc number to show variables
	# n1: top n1 number of variables
	
	pls_data <- pls.model(scal.data,choice,A)
	# pearson relation: variable vs score
	datafile <- pls_data$model$v
	cor_score <- cor(datafile,pls_data$scores[,comp])
	# pearson covarance: variable vs score
	cov_score <- cov(datafile,pls_data$scores[,comp])
	# VIP
	x <- as.matrix(datafile,nrow=nrow(datafile))
	y <- as.factor(pls_data$model$y)	
	pls <- plsda(x,y,ncomp=A)
	vip_value <-vip(pls)[,comp]
	
	imp_vip1 <- length(vip_value[vip_value >=1])
	imp_vip2 <- length(vip_value[vip_value >=1.5])	
	imp_vip3 <- length(vip_value[vip_value >=2])
	
	filterVIP <- datafile[,names(vip_value[vip_value >=1])]
	
	#Pcorr
	y <- as.numeric(pls_data$scores[,comp])	
	pvalue_list <- c()	
	for (i in 1:ncol(datafile)) 
	{
		#calculating the correlation coefficients and sig. level
		result <- cor.test(datafile[,i],y)
		#result <- pairwise.t.test(scal.data$v[,i], scal.data$y)
		pvalue_list <- c(pvalue_list,round(as.numeric(result$p.value),2))
	}
	
	datamatrix <- data.frame(cbind(cor_score,cov_score,vip_value,pvalue_list))
	# P value calculates significant level of corr
	colnames(datamatrix) <- c("cor","cov","vip","p")
	
	fiter_data1 <- subset(datamatrix,datamatrix$vip>1)
	fiter_data2 <- subset(datamatrix,datamatrix$p<0.05)
	
	#x11()	# the plot is large, so x11() is applied
	if (type=="VIP") 
	{
		vip_value <- rev(sort(vip(pls)[,comp]))
		plot(vip_value,type="h",col="green",xlab="Sample Num",ylab=paste("VIP[",comp,"]",sep=""),main="PLS-VIP")
		legend("topright",c(paste(imp_vip1," variables >= 1",sep=""),paste(imp_vip2," variables >= 1.5",sep=""),
						paste(imp_vip3," variables >= 2",sep="")))}
	
	if (type=="SPlot") 
	{
		plot(cov_score,cor_score,pch=3,col="blue",xlab="Cov(t,x)",ylab="Corr(t,x)",main="PLS-SPlot")
		points(fiter_data1$cov,fiter_data1$cor,pch=18,col="red")
		points(fiter_data2$cov,fiter_data2$cor,pch=2,col="green")
		legend("topleft",c("VIP >1","p<0.05"),pch=c(18,2),col=c("red","green")) 
	}
	
	if (type=="VPlot") {
		plot(cor_score,vip_value,pch=3,col="blue",xlab="Corr(t,x)",ylab=paste("VIP[",comp,"]",sep=""),main="PLS-VPlot")
		points(fiter_data1$cor,fiter_data1$vip,pch=18,col="red")
		points(fiter_data2$cor,fiter_data2$vip,pch=2,col="green")
		legend("topleft",c("VIP >1","p<0.05"),pch=c(18,2),col=c("red","green")) 
	}
	
	if (type=="HeatMap") {
		heatmap(filterVIP,col=rainbow(256),xlab="Heat Map")
	}
}




# calculating Q2 value (plspm)
# x: n*m data matrix
# y: group information
# nc: pcs
# cv: TRUE or FALSE

pls.Q2 <- function(x, y, nc=2, cv=FALSE)
{
	# ============ checking arguments ============
	X <- as.matrix(x)
	n <- nrow(X)
	p <- ncol(X)
	Y <- as.matrix(y)
	if (any(is.na(X))) na.miss<-TRUE else na.miss<-FALSE       
	if (na.miss) cv <- FALSE
	# ============ setting inputs ==============
	Xx <- scale(X) 
	Yy <- scale(Y)
	X.old <- Xx
	Y.old <- Yy
	# RSS first value is n-1
	RSS <- c(n-1, rep(NA, nc))
	PRESS <- rep(NA, nc)
	Q2 <- NULL
	# ============ pls regression algorithm ==============
	w.old <- rep(1,p)
	t.new <- rep(1,n)
	p.new <- rep(NA, p)
	h <- 1
	repeat
	{
		if (!na.miss) # no missing data
		{
			w.old <- t(X.old) %*% Y.old / sum(Y.old^2)
			w.new <- w.old / sqrt(sum(w.old^2)) # normalization
			t.new <- X.old %*% w.new
			p.new <- t(X.old) %*% t.new / sum(t.new^2)
		}
		
		c.new <- t(Y.old) %*% t.new / sum(t.new^2)
		u.new <- Y.old / as.vector(c.new)
		
		if (!na.miss) 
		{
			# cross validation "leave-one-out"
			RSS[h+1] <-  sum((Y.old - t.new%*%c.new)^2)
			press <- rep(0,n)
			for (i in 1:n)
			{
				Xy.aux <- t(X.old[-i,]) %*% Y.old[-i]
				wh.si <- Xy.aux %*% sqrt(solve(t(Xy.aux)%*%Xy.aux))
				th.si <- X.old[-i,] %*% wh.si
				ch.si <- t(Y.old[-i]) %*% th.si %*% solve(t(th.si)%*%th.si)
				ch.si <- as.vector(ch.si)
				Yhat.si <- ch.si * X.old[i,] %*% wh.si
				press[i] <- (Y.old[i] - Yhat.si)^2
			}
			PRESS[h] = sum(press)
			Q2[h] = 1 - PRESS[h]/RSS[h]
		}
		
		Y.old <- Y.old - t.new%*%c.new# deflate y.old
		X.old <- X.old - (t.new %*% t(p.new))# deflate X.old
		if (!na.miss)
			if (cv) 
				if (Q2[h]<0 || h==nc) break
		if (!cv)
			if (h==nc) break
		h <- h + 1
	}
	
	q2cum <- rep(NA,h)
	for (k in 1:h)
		q2cum[k] <- prod(PRESS[1:k]) / prod(RSS[1:k])
	Q2cum <- round((1 - q2cum), 4)
	result <- rbind(Q2,Q2cum)
	return(result)
}

PLS.pred <- function(m.data,t.data,choice,A,a1=1,a2=2) {
	# m.data: training set after scalling
	# t.data: test set after same processing procedure
	
	model <- pls.model(m.data,choice,A) # model
	p.model <- predict(model,ncomp=c(a1,a2),newdata=t.data,type="scores")	
	p.model2 <- predict(model,ncomp=c(a1,a2),newdata=t.data)
	par(mfrow=c(1,2))
	pls_score.x(m.data,choice,A,a1,a2) # score plot	
	points(p.model,col="green",pch=19)
	
	plot(model$validation$pred[,,a1],col=m.data$y+1,main="Predicted Response",
			xlab="Sample Num",ylab=paste("response@PC",a1,sep="")) # response plot
	points(p.model2[,,a1],col="green",pch=19)
}


pls <- function(params,scal.data)
{
	workdir <- params$WorkDir
	f <- as.integer(params$PLS_PCs)
	A <- as.integer(params$PLS_CID)
	#work.dir <- getwd()
#	setwd(workdir)
#	work.dir <- getwd()
	create.dir <- paste(workdir,"output/PLS/",sep="")
	dir.create(create.dir)
	
	pls_data <- pls.model(scal.data,"LOO",f)
	pls_scores <- pls_data$scores
	pls_loadings <- pls_data$loadings
	write.csv(pls_scores,paste(create.dir,"pls-scores.csv",sep="/"))
	write.csv(pls_loadings,paste(create.dir,"pls-loadings.csv",sep="/"))

#	pdf(paste(create.dir,"PLSDA.pdf",sep=""))
#	rownum <- as.integer(f-1)
#	par(mfrow=c(rownum,2)) 
#	for (i in c(2:f)) {
#		pls_score.x(scal.data,"LOO",f,1,i)
#		pls_loading(scal.data,"LOO",f,1,i)
#	}
#	dev.off()
		
	pdf(paste(create.dir,"PLS-scores.pdf",sep=""))
	for (i in c(2:f)) {
		pls_score.x(scal.data,"LOO",f,1,i)
	}
	dev.off()
	
	
	pdf(paste(create.dir,"PLS-loadings.pdf",sep=""))
	for (i in c(2:f)) {
		pls_loading(scal.data,"LOO",f,1,i)
	}
	dev.off()
	
	pdf(paste(create.dir,"PLS-modelSummary.pdf",sep=""))	
	pls_varcum (scal.data,"LOO",f)	
	dev.off()	
	
	pdf(paste(create.dir,"PLS-DistanceToModel.pdf",sep=""))
	par(mfrow=c(f,1)) 	
	for (i in c(1:f)) {
		PLS_Dmodx (scal.data,i,"LOO",f)
		}
	dev.off()
		
	# 1st component
	pls.feature.select (workdir,scal.data,"LOO",f,1)
	
	pdf(paste(create.dir,"PLS-VarPlots.pdf",sep=""))	
	par(mfrow=c(2,2))
	pls.VarPlot (scal.data,"LOO",f,A,"SPlot")
	pls.VarPlot (scal.data,"LOO",f,A,"VPlot")
	pls.VarPlot (scal.data,"LOO",f,A,"VIP")
	dev.off()
	
	pdf(paste(create.dir,"PLS-Heatmap.pdf",sep=""))	
	pls.VarPlot (scal.data,"LOO",f,A,"HeatMap")
	dev.off()
	
	pdf(paste(create.dir,"PLS-prediction.pdf",sep=""))	
	pls_prediction (scal.data,"LOO",f)
	dev.off()
	
	# Default is PC 2
	pdf(paste(create.dir,"PLS-VariableVariance.pdf",sep=""))	
	PLS.R2VX.showmore (scal.data,A,"LOO",f)
	dev.off()
	
	print("PLS Done!")
}
