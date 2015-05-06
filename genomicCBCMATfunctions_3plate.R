packageLoad = function(x) {
	# returns a matrix of the installed packages
	inst = installed.packages()
	matchVals = match(x, inst[,"Package"])
	missing = x[is.na(matchVals)]
	if (length(missing)>0) 
	{
		cat("Installing packages...", "\n")
		lapply(missing, install.packages, repos="http://www.revolution-computing.com/cran/")
	}
	if ("limma" %in% missing)
	{
		source("http://bioconductor.org/biocLite.R")
		biocLite("limma")
	}
	# then load libraries
	for (i in 1:length(x))
	{
		require(x[i], character.only=TRUE)
	}
}

readData<-function(path, dataFile)
{
	dataPath<-paste(path, dataFile, sep="")
	checkFile<-unlist(strsplit(dataFile, "\\."))
	fileType<-tolower(checkFile[length(checkFile)])
	if (fileType=="csv")
	  	curData<-read.csv(dataPath)
	if (fileType=="txt" | fileType=="tsv")
		curData<-read.delim(dataPath)
	return(curData)
}

# subset samples and genes to only those in design and moduleAssign
subsetExpData<-function(expData, designData, moduleAssignData)
{
	subExpData<-expData
	# first remove any columns (samples that are not in designData)
	matchIndex<-match(colnames(subExpData), designData$sampleID)
	if (any(is.na(matchIndex)))
	{
		remIndex<-which(is.na(matchIndex))
		subExpData<-subExpData[,-remIndex]
	}
	if (any(!(designData$groupID %in% c(0,1))))
	{
		# remove these samples that are not controls or cases (may be NTC or reference)
		keepIds<-designData$sampleID[c(which(designData$groupID==0), which(designData$groupID==1))]
		# check - yes this is correct
		#table(designData$groupID[match(keepIds, designData$sampleID)])		
		matchIndex<-match(keepIds, colnames(subExpData))
		keepIndex<-matchIndex[which(!is.na(matchIndex))]
		subExpData<-subExpData[,keepIndex]
		# check that this is correct - yes
		#table(designData$groupID[match(colnames(subExpData), designData$sampleID)])
	}
	
	# next remove rows (such as housekeeping genes)
	matchIndex<-match(rownames(subExpData), moduleAssignData$geneID)
	if (any(is.na(matchIndex)))
	{
		remIndex<-which(is.na(matchIndex))
		subExpData<-subExpData[-remIndex,]
	}	
	
	return(subExpData)
}

# calculate FC values compared to controls rather than reference samples
FCtocontrols<-function(subExpData, designData)
{
	controlIndex<-which(!is.na(match(colnames(subExpData), designData$sampleID[which(designData$groupID==0)])))
	# convert to delta delta Ct values before standardizing to control means
#	controlDelta<- -log2(subExpData[,controlIndex])
	controlDelta<-apply(subExpData[,controlIndex], 2, function(x)
	{
		-log2(as.numeric(as.character(x)))
	})
	controlMeans<-apply(controlDelta, 1, mean, na.rm=T)	
#	allDelta<- -log2(subExpData)
	allDelta<-apply(subExpData, 2, function(x)
	{
		-log2(as.numeric(as.character(x)))
	})
	allDeltaToControls<-apply(allDelta, 2, function(x)
	{
		x-controlMeans
	})
#	allFCtoControls<-2^(-allDeltaToControls)
	allFCtoControls<-apply(allDeltaToControls, 2, function(x)
	{
		2^(-as.numeric(as.character(x)))
	})
	rownames(allFCtoControls)<-rownames(subExpData)
	return(allFCtoControls)
}

# convert the FC values from gene level to module level
convertFCtoModule<-function(FCvals, moduleAssignData)
{
	# loop through all modules
	uniModules<-sort(as.character(unique(moduleAssignData$module)))
	FCtoModule<-c()
	for (i in 1:length(uniModules))
	{
		curMod<-uniModules[i]
		curGeneIDs<-moduleAssignData$geneID[which(curMod==moduleAssignData$module)]
		matchIndex<-match(curGeneIDs, rownames(FCvals))
		if (any(is.na(matchIndex)))
			matchIndex<-matchIndex[-which(is.na(matchIndex))]
		curMeans<-2^apply(log2(FCvals[matchIndex,]), 2, mean, na.rm=T)
		FCtoModule<-rbind(FCtoModule, curMeans)
	}	
	rownames(FCtoModule)<-uniModules
	#apply(FCtoModule, 1, function(x) {sum(is.na(x))})
	return(FCtoModule)
}

# run model fits at the gene level
fitLimmaToGenes<-function(FCvals, designData, moduleAssignData)
{
	# make a design matrix - first order the columns in FCvals to match those in designData
	matchIndex<-match(colnames(FCvals), designData$sampleID)
	curGroups<-designData$groupID[matchIndex]
	if (all(curGroups %in% c(0,1)))
	{
		curMM<-model.matrix(~1+factor(curGroups))
		FCvals2<-apply(FCvals, 2, function(x)
		{
			log2(as.numeric(as.character(x)))
		})
		rownames(FCvals2)<-rownames(FCvals)
#		curFit<-lmFit(log2(FCvals), curMM)
		curFit<-lmFit(FCvals2, curMM)
		curFit<-eBayes(curFit)
		curTT<-topTable(curFit, coef=2, adjust.method="BH", number=nrow(FCvals))
		# output: gene name, module, test statistic, p-value, and adjusted p-value
		# need to add module to curTT
		curTT<-curTT[,1:5]
		curTT$module<-moduleAssignData$module[match(rownames(curTT), moduleAssignData$geneID)]
		curTT$geneName<-moduleAssignData$geneName[match(rownames(curTT), moduleAssignData$geneID)]
		return(curTT)
	}
	else
		print("ERROR: Columns in FC matrix belong to a group that is not case or control.")
}

# assumes the FC values are already at the module level
# starts with module level FC and compares case to control
CaseVsControlModuleLevelFC<-function(FCvals, designData)
{
	controlMatch<-match(designData$sampleID[which(designData$groupID==0)], colnames(FCvals))
	if (any(is.na(controlMatch)))
		controlMatch<-controlMatch[-which(is.na(controlMatch))]
	controlVals<-FCvals[,controlMatch]
	caseMatch<-match(designData$sampleID[which(designData$groupID==1)], colnames(FCvals))
	if (any(is.na(caseMatch)))
		caseMatch<-caseMatch[-which(is.na(caseMatch))]
	caseVals<-FCvals[,caseMatch]

	# take mean on log2 scale
	controlMeans<-2^apply(log2(controlVals), 1, mean, na.rm=T)
	caseMeans<-2^apply(log2(caseVals), 1, mean, na.rm=T)
	
	modMeans<-caseMeans/controlMeans
	return(modMeans)
}

# also try second method of calculating average
# here start with FC values at the gene level
# calculate case vs. control for each gene and then roll up to module level 
CaseVsControlModuleLevelFC2<-function(FCvals, designData, moduleAssignData)
{
	controlMatch<-match(designData$sampleID[which(designData$groupID==0)], colnames(FCvals))
	if (any(is.na(controlMatch)))
		controlMatch<-controlMatch[-which(is.na(controlMatch))]
	controlVals<-FCvals[,controlMatch]
	caseMatch<-match(designData$sampleID[which(designData$groupID==1)], colnames(FCvals))
	if (any(is.na(caseMatch)))
		caseMatch<-caseMatch[-which(is.na(caseMatch))]
	caseVals<-FCvals[,caseMatch]
	
	# take mean on log2 scale
	controlMeans<-2^apply(log2(controlVals), 1, mean, na.rm=T)
	caseMeans<-2^apply(log2(caseVals), 1, mean, na.rm=T)
	
	geneMeans<-caseMeans/controlMeans
	
	# now roll up to module level
	uniModules<-sort(as.character(unique(moduleAssignData$module)))
	modMeans<-c()
	for (i in 1:length(uniModules))
	{
		curMod<-uniModules[i]
		curIDs<-moduleAssignData$geneID[which(curMod==moduleAssignData$module)]
		matchIndex<-match(curIDs, names(geneMeans))
		if (any(is.na(matchIndex)))
			matchIndex<-matchIndex[-which(is.na(matchIndex))]
		modMeans<-c(modMeans, 2^mean(log2(geneMeans[matchIndex]), na.rm=T))
	}
	names(modMeans)<-uniModules
	return(modMeans)
}

# look at correlation of genes within modules
modGeneCorrelation<-function(FCvals, moduleAssignData, designData, fileName)
{
	uniModules<-unique(as.character(moduleAssignData$module))
	modCor<-list()
	for (i in 1:length(uniModules))
	{
		curMod<-uniModules[i]
		curGenes<-moduleAssignData$geneID[which(curMod==moduleAssignData$module)]
		matchIndex<-match(curGenes, rownames(FCvals))
		if (any(is.na(matchIndex)))
			matchIndex<-matchIndex[-which(is.na(matchIndex))]
		# subset FC to only cases
		cases<-as.character(designData$sampleID[which(designData$groupID==1)])
		caseIndex<-match(cases, colnames(FCvals))
		corMatrix<-cor(t(FCvals[matchIndex,caseIndex]), use="complete.obs")
		modCor[[i]]<-corMatrix
	}
	yvals<-c()
	xvals<-c()
	for (i in 1:length(modCor))
	{	
		newYvals<-modCor[[i]][upper.tri(modCor[[i]])]
		yvals<-c(yvals, newYvals)
		xvals<-c(xvals, rep(i, length(newYvals)))
	}
	# reorder the modules based on median correlation
	orderIndex<-order(tapply(yvals, xvals, median))

	pdf(file=fileName, width=11)
	plot(yvals~factor(xvals, levels=orderIndex), ylab="Correlation", xlab="", axes=FALSE, ylim=c(min(yvals, na.rm=T), 1), main="Correlation of Cases Only")
	box()
	# color by plate
	for (i in 1:length(uniModules))
	{
		curModule<-uniModules[orderIndex][i]
		if (moduleAssignData$plate3[which(curModule==moduleAssignData$module)][1]=="No")
			curCol<-"black"
		else
			curCol<-"red"
		axis(1, at=i, labels=uniModules[orderIndex][i], las=2, col=curCol, col.axis=curCol)
	}
	#axis(1, at=1:length(uniModules), labels=uniModules[orderIndex], las=2)
	axis(2, at=seq(-1,1,0.2), labels=seq(-1,1,0.2), las=2)
	for (i in seq(-1,1,0.2))
		abline(h=i, col="gray", lty=2)		
	dev.off()
}

# fit module level model - use gene and individual as random effect
moduleLevelModelWithRandomEffect<-function(FCcomparedToControls, designData, moduleAssignData)
{	
	designSub<-designData[which(designData$groupID %in% c(0,1)),]
	uniModules<-as.character(unique(moduleAssignData$module))	
	
	# need to rearrange designSub to match column names in designSub
	designSub<-designSub[match(colnames(FCcomparedToControls), designSub$sampleID),]
	all(designSub$sampleID==colnames(FCcomparedToControls))
	
	# fit a model for each module
	groupPvals<-c()
	groupEst<-c()
	for (i in 1:length(uniModules))
	{
		# get the genes
		curMod<-uniModules[i]
		geneID<-as.character(moduleAssignData$geneID[which(moduleAssignData$module==curMod)])
		modelData<-c()
		for (j in 1:length(geneID))
		{
			curGene<-geneID[j]
			expIndex<-which(curGene==rownames(FCcomparedToControls))
			modelData<-rbind(modelData, cbind(as.vector(as.matrix(FCcomparedToControls[expIndex,])), designSub$groupID, rep(curGene, nrow(designSub)), colnames(FCcomparedToControls)))
		}
		modelData<-as.data.frame(modelData)
		colnames(modelData)<-c("exp","group","gene","individual")
		modelData$exp<-as.numeric(as.character(modelData$exp))
		modelData$group<-as.factor(as.character(modelData$group))
		modelData$gene<-as.factor(as.character(modelData$gene))
		modelData$individual<-as.factor(as.character(modelData$individual))
		
		# now make model
		largeMod<-lmer(log2(exp) ~ group + (1|gene) + (1|individual), data=modelData, REML=FALSE)
		smallMod<-lmer(log2(exp) ~ 1 + (1|gene) + (1|individual), data=modelData, REML=FALSE)
		groupPvals<-c(groupPvals, anova(largeMod, smallMod)[2,8])
		groupEst<-c(groupEst, summary(largeMod)$coef[2,1])
	}
#	hist(groupPvals, breaks=20)
	groupPvalsAdj<-p.adjust(groupPvals, method="fdr")
#	hist(groupPvalsAdj, breaks=20)
	modLevelResults<-cbind(groupPvals, groupPvalsAdj, groupEst)
	modLevelResults<-as.data.frame(modLevelResults)
	rownames(modLevelResults)<-uniModules
	
	colnames(modLevelResults)<-c("pval","fdr","groupParameter")
	modLevelResults<-modLevelResults[order(modLevelResults$fdr),]
	return(modLevelResults)
}







################
# PLOTTING FUNCTIONS
################

#####
# PCA plot
#####
makePCAplots<-function(FCvals, designData, resultsDirectory)
{
	# remove samples and genes that have > 10 missing values before running PCA
	FCorig<-FCvals
	rowCount<-apply(FCvals, 1, function(x) {sum(is.na(x))})
	if (any(rowCount > 10))
	{
		remRows<-which(apply(FCvals, 1, function(x) {sum(is.na(x))}) > 10)
		FCvals<-FCvals[-remRows,]
	}
	colCount<-apply(FCvals, 2, function(x) {sum(is.na(x))})
	if (any(colCount > 10))
	{
		remCols<-which(apply(FCvals, 2, function(x) {sum(is.na(x))}) > 10)
		FCvals<-FCvals[,-remCols]
	}	
	if (nrow(designData) > ncol(FCvals))
	{
		remIndex<-which(!(designData$sampleID %in% colnames(FCvals)))
		tempDesign<-designData[-remIndex,]
		# now reorder
		tempDesign<-tempDesign[match(colnames(FCvals), tempDesign$sampleID),]
		# all(tempDesign$sampleID==colnames(FCvals))
	}
	else
	{
		tempDesign<-designData
		tempDesign<-tempDesign[match(colnames(FCvals), tempDesign$sampleID),]
	}
	# now make 2 plots - one with all remaining genes and one with all remaining samples
	######## plot 1: remove samples with NAs (leave in all remaining genes)
	samPCAdat<-as.data.frame(t(log2(FCvals)))
	remSamsPCA<-prcomp(~ ., data=samPCAdat, na.action=na.omit, scale=TRUE, center=TRUE)
	numPCAsams<-sum(apply(FCvals, 2, function(x) {sum(is.na(x))})==0)
	numOrigSams<-ncol(FCorig)
	PC1var<-(remSamsPCA$sdev[1]^2)/(sum(remSamsPCA$sdev^2))*100
	PC2var<-(remSamsPCA$sdev[2]^2)/(sum(remSamsPCA$sdev^2))*100
	possColors<-rainbow(length(unique(tempDesign$groupName)))
	uniGroups<-unique(tempDesign$groupName)
	# color by groupName
	curCols<-rep("", nrow(remSamsPCA$x))
	for (i in 1:nrow(remSamsPCA$x))
	{
		tempRN<-rownames(remSamsPCA$x)[i]
		tempGroup<-tempDesign$groupName[which(tempRN==tempDesign$sampleID)]
		curCols[i]<-possColors[which(tempGroup==uniGroups)]
	}
		
	pdf(file=paste(resultsDirectory,"PCA_showingGroups_AllGenes.pdf",sep=""))
	plot(x=remSamsPCA$x[,1], y=remSamsPCA$x[,2], xlab=paste("PC1 (",round(PC1var,1),"%)",sep=""), ylab=paste("PC2 (",round(PC2var,1),"%)",sep=""), main=paste("Original # of Samples:", numOrigSams, "\n# Samples in PCA:", numPCAsams), pch=19, col=curCols)
	# add a legend
	legend(x=max(remSamsPCA$x[,1]), y=max(remSamsPCA$x[,2]), xjust=1, yjust=1, pch=19, col=possColors, legend=uniGroups, cex=0.8)
	dev.off()
	
	######## plot 2: remove genes with NAs (leave in all remaining samples)
	remIndex<-which(apply(FCvals, 1, function(x) {sum(is.na(x))}) > 0)
	numOrigGenes<-nrow(FCorig)
	if (length(remIndex) > 0)
		FCvalsSub<-FCvals[-remIndex,]
	else
		FCvalsSub<-FCvals
	numPCAgenes<-nrow(FCvalsSub)
		
	genePCAdat<-as.data.frame(t(log2(FCvalsSub)))
	remGenesPCA<-prcomp(~ ., data=genePCAdat, na.action=na.omit, scale=TRUE, center=TRUE)
	gPC1var<-(remGenesPCA$sdev[1]^2)/(sum(remGenesPCA$sdev^2))*100
	gPC2var<-(remGenesPCA$sdev[2]^2)/(sum(remGenesPCA$sdev^2))*100
	# color by groupName
	gcurCols<-rep("", nrow(remGenesPCA$x))
	for (i in 1:nrow(remGenesPCA$x))
	{
		tempRN<-rownames(remGenesPCA$x)[i]
		tempGroup<-tempDesign$groupName[which(tempRN==tempDesign$sampleID)]
		gcurCols[i]<-possColors[which(tempGroup==uniGroups)]
	}

	pdf(file=paste(resultsDirectory,"PCA_showingPlates_AllSamples.pdf",sep=""))
	plot(x=remGenesPCA$x[,1], y=remGenesPCA$x[,2], xlab=paste("PC1 (",round(gPC1var,1),"%)",sep=""), ylab=paste("PC2 (",round(gPC2var,1),"%)",sep=""), main=paste("Original # of Genes:", numOrigGenes, "\n# Genes in PCA:", numPCAgenes), pch=19, col=gcurCols)
	# add a legend
	legend(x=max(remGenesPCA$x[,1]), y=max(remGenesPCA$x[,2]), xjust=1, yjust=1, pch=19, col=possColors, legend=uniGroups, cex=0.8)
	dev.off()
}

#####
# set up the coordinates to plot the polygon inside a circle
# assume we always start from the top vertical line in a circle
#####
calcFillPolygon<-function(centerx, centery, radius, percentPos, percentNeg, units="npc", toplot="y")
{
	# calculate angle
	curAnglePos<-percentPos*360
	curLength<-100
	allx<-c()
	ally<-c()
	allxn<-c()
	allyn<-c()
	# calculate this first part for all angles
	if (percentPos > 0)
	{
		# then in the top right quadrant of the circle
		# find the x-value to stop at
		opAnglePos<-90-curAnglePos
		if (opAnglePos < 0)
			opAnglePos<-0
		# cosine needs the parameter in radians
		stopAtx<-radius*cos(opAnglePos*pi/180)
		
		# now calculate the x and y values for the polygon
		# start at x value of zero for this quadrant of the circle
		startX<-0
		xvals<-seq(startX, stopAtx, length=curLength)
		yvals<-sqrt(radius^2-xvals^2)
		allx<-c(centerx, centerx, centerx+xvals)
		ally<-c(centery, centery+radius, centery+yvals)
	}
	if (percentPos > 0.25)
	{
		# in the bottom right quadrant of the circle
		# find the x-value to stop at
		opAnglePos<-curAnglePos-90
		if (opAnglePos > 90)
			opAnglePos<-90
		# cosine needs the parameter in radians
		stopAtx<-radius*cos(opAnglePos*pi/180)
		
		# now calculate the x and y values for the polygon
		# start at x value equal to the radius for this quadrant of the circle
		startX<-radius
		xvals<-seq(startX, stopAtx, length=curLength)
		yvals<-sqrt(radius^2-xvals^2)
		allx<-c(allx, centerx+xvals)
		ally<-c(ally, centery-yvals)
	}
	if (percentPos > 0.5)
	{
		# in the bottom left quadrant of the circle
		# find the x-value to stop at
		opAnglePos<-270-curAnglePos
		if (opAnglePos < 0)
			opAnglePos<-0
		stopAtx<--radius*cos(opAnglePos*pi/180)
		
		# now calculate the x and y values for the polygon
		# start at x value of zero for this quadrant of the circle
		startX<-0
		xvals<-seq(startX, stopAtx, length=curLength)
		yvals<--sqrt(radius^2-xvals^2)
		allx<-c(allx, centerx+xvals)
		ally<-c(ally, centery+yvals)
	}
	if (percentPos > 0.75)
	{
		# in the top left quadrant of the circle
		# find the x-value to stop at
		opAnglePos<-curAnglePos-270
		stopAtx<--radius*cos(opAnglePos*pi/180)
		
		# now calculate the x and y values for the polygon
		# start at x value equal to negative radius for this quadrant of the circle
		startX<- -radius
		xvals<-seq(startX, stopAtx, length=curLength)
		yvals<-sqrt(radius^2-xvals^2)
		allx<-c(allx, centerx+xvals)
		ally<-c(ally, centery+yvals)
	}
	if (toplot=="y" & percentPos > 0)
		grid.polygon(x=allx, y=ally, gp=gpar(fill="red"), default.units=units)
	
	# now work on the negative percent
	curAngleNeg<-percentNeg*360
	# calculate this first part for all angles
	if (percentNeg > 0)
	{
		# then in the top left quadrant of the circle (going backwards for negative percent)
		# find the x-value to stop at
		opAngleNeg<-90-curAngleNeg
		if (opAngleNeg < 0)
			opAngleNeg<-0
		# cosine needs the parameter in radians
		stopAtx<--radius*cos(opAngleNeg*pi/180)
		
		# now calculate the x and y values for the polygon
		# start at x value of zero for this quadrant of the circle
		startX<-0
		xvals<-seq(startX, stopAtx, length=curLength)
		yvals<-sqrt(radius^2-xvals^2)
		allxn<-c(centerx, centerx, centerx+xvals)
		allyn<-c(centery, centery+radius, centery+yvals)
	}
	if (percentNeg > 0.25)
	{
		# in the bottom left quadrant of the circle (going backwards for negative percent)
		# find the x-value to stop at
		opAngleNeg<-curAngleNeg-90
		if (opAngleNeg > 90)
			opAngleNeg<-90
		# cosine needs the parameter in radians
		stopAtx<--radius*cos(opAngleNeg*pi/180)
		
		# now calculate the x and y values for the polygon
		# start at x value equal to the radius for this quadrant of the circle
		startX<--radius
		xvals<-seq(startX, stopAtx, length=curLength)
		yvals<--sqrt(radius^2-xvals^2)
		allxn<-c(allxn, centerx+xvals)
		allyn<-c(allyn, centery+yvals)
	}
	if (percentNeg > 0.5)
	{
		# in the bottom right quadrant of the circle (going backwards for negative percent)
		# find the x-value to stop at
		opAngleNeg<-270-curAngleNeg
		if (opAngleNeg < 0)
			opAngleNeg<-0
		stopAtx<-radius*cos(opAngleNeg*pi/180)
		
		# now calculate the x and y values for the polygon
		# start at x value of zero for this quadrant of the circle
		startX<-0
		xvals<-seq(startX, stopAtx, length=curLength)
		yvals<--sqrt(radius^2-xvals^2)
		allxn<-c(allxn, centerx+xvals)
		allyn<-c(allyn, centery+yvals)

	}
	if (percentNeg > 0.75)
	{
		# in the top right quadrant of the circle
		# find the x-value to stop at
		opAngleNeg<-curAngleNeg-270
		stopAtx<-radius*cos(opAngleNeg*pi/180)
		
		# now calculate the x and y values for the polygon
		# start at x value equal to negative radius for this quadrant of the circle
		startX<-radius
		xvals<-seq(startX, stopAtx, length=curLength)
		yvals<-sqrt(radius^2-xvals^2)
		allxn<-c(allxn, centerx+xvals)
		allyn<-c(allyn, centery+yvals)
	}
	if (toplot=="y" & percentNeg > 0)
		grid.polygon(x=allxn, y=allyn, gp=gpar(fill="blue"), default.units=units)
	
	return(list(allx=allx, ally=ally, allxn=allxn, allyn=allyn))
}


make6x11gridPercent<-function(modPercent, fileName, scale=0.33)
{
	if (nrow(modPercent)==66)
	{
		# draw the 3 TFA plates
		png(fileName, width=14*scale, height=8*scale, units="in", res=200, pointsize=20)
		pushViewport(viewport(width=unit(11*scale, "inches"), height=unit(6*scale, "inches"), 
			x=unit(2*scale, "inches"), y=unit(1*scale, "inches"), default.units="inches", 
			just=c("left", "bottom")))
		vpWidth<-11*scale
		vpHeight<-6*scale
		numCol<-11
		numRow<-6
		
		lineType<-1
		lineWidth<-0.5
		lineColor<-"gray"
		for (i in seq(0,1,1/numRow))
			grid.lines(x=seq(0,1,1/numCol)*vpWidth, y=i*vpHeight, gp=gpar(lwd=lineWidth, lty=lineType, col=lineColor), default.units="inches")
		for (i in seq(0,1,1/numCol))
			grid.lines(x=i*vpWidth, y=seq(0,1,1/numRow)*vpHeight, gp=gpar(lwd=lineWidth, lty=lineType, col=lineColor), default.units="inches")
			
		# 1/23/14 add a dark line to separate plates 1 + 2 and plate 3
		grid.lines(x=seq(0,1,1/numCol)*vpWidth, y=2/6*vpHeight, gp=gpar(lwd=1, lty=lineType, col="black"), default.units="inches")

		# now need to add circles and color
		x<-seq(1/(numCol*2), 1, 1/numCol)
		y<-seq(1/(numRow*2), 1, 1/numRow)
		allX<-c(x,x,x,x,x,x)
		# start at top and go down from left to right
		allY<-c(rep(y[6],11), rep(y[5],11), rep(y[4],11), rep(y[3],11), rep(y[2],11), rep(y[1],11))
		
		# convert values to inches
		nallX<-allX*numCol*scale
		nallY<-allY*numRow*scale
		
		curPosData<-modPercent[,1]
		curNegData<-modPercent[,2]

		# now create polygons		
		k<-1
		posx<-c()
		posid<-c()
		posy<-c()
		negx<-c()
		negid<-c()
		negy<-c()

		for (i in 1:nrow(modPercent))
		{
			grid.circle(nallX[i], nallY[i], r=(1/2.4)*scale, default.units="inches")
			if (curPosData[i] != 0 | curNegData[i] != 0)
			{
				retList<-calcFillPolygon(centerx=nallX[i], centery=nallY[i], radius=(1/2.4)*scale,
					percentPos=curPosData[i], percentNeg=curNegData[i],
					units="inches", toplot="n")
							
				posx<-c(posx, c(retList$allx))
				posid<-c(posid, rep(i, length(retList$allx)))
				posy<-c(posy, c(retList$ally))
				negx<-c(negx, c(retList$allxn))
				negid<-c(negid, rep(i, length(retList$allxn)))
				negy<-c(negy, c(retList$allyn))
			}
		}
		# now plot the polygons if there is data to plot
		if (length(posx)>0)
			grid.polygon(x=posx, y=posy, id=posid, gp=gpar(fill="red"), default.units="inches")
		if (length(negx)>0)
			grid.polygon(x=negx, y=negy, id=negid, gp=gpar(fill="blue"), default.units="inches")
		
		dev.off()
	}
}

# return a matrix of percent up and down by module
statSigCountByModule<-function(geneLevelStats)
{
	uniModules<-sort(as.character(unique(geneLevelStats$module)))
	modLevelPercents<-c()
	for (i in 1:length(uniModules))
	{
		curMod<-uniModules[i]
		matchIndex<-which(curMod==geneLevelStats$module)
		countUp<-sum(geneLevelStats$adj.P.Val[matchIndex] < 0.05 & geneLevelStats$t[matchIndex] > 0, na.rm=T)
		countDown<-sum(geneLevelStats$adj.P.Val[matchIndex] < 0.05 & geneLevelStats$t[matchIndex] < 0, na.rm=T)
		modLevelPercents<-rbind(modLevelPercents, c(countUp/4, countDown/4))
	}
	rownames(modLevelPercents)<-uniModules
	colnames(modLevelPercents)<-c("PercentUp","PercentDown")
	return(modLevelPercents)
}

# order the modules by round and number (may change to reorder based on annotation)
reorderModules<-function(modFC, moduleAssignData)
{
	curRN<-rownames(modFC)
	# now split into plates 1 and 2 and plate 3 (put plate 3 modules at the end)
	RNplate12<-curRN[curRN %in% moduleAssignData$module[which(moduleAssignData$plate3=="No")]]
	RNplate3<-curRN[curRN %in% moduleAssignData$module[which(moduleAssignData$plate3=="Yes")]]
	
	modRound12<-unlist(lapply(strsplit(RNplate12, "\\."), function(x) {x[1]}))
	modNumber12<-as.numeric(unlist(lapply(strsplit(RNplate12, "\\."), function(x) {x[2]})))
	modRound12<-as.numeric(substr(modRound12, 2, nchar(modRound12)))
	orderIndex12<-order(modRound12, modNumber12)

	modRound3<-unlist(lapply(strsplit(RNplate3, "\\."), function(x) {x[1]}))
	modNumber3<-as.numeric(unlist(lapply(strsplit(RNplate3, "\\."), function(x) {x[2]})))
	modRound3<-as.numeric(substr(modRound3, 2, nchar(modRound3)))
	orderIndex3<-order(modRound3, modNumber3)

	modFC12<-modFC[match(RNplate12, rownames(modFC)),]
	modFC3<-modFC[match(RNplate3, rownames(modFC)),]
	
	newModFC12<-modFC12[orderIndex12,]
	newModFC3<-modFC3[orderIndex3,]
	newModFC<-as.data.frame(rbind(newModFC12, newModFC3))

	return(newModFC)
}

# make the case vs. control grid of 66 modules
# maxVal is the maximum log2(FC) to show in the color key
make6x11grid<-function(modFC, fileName, scale=0.33, maxVal=2, sigCutoff=1, modLevRes)
{
	if (nrow(modFC)==66)
	{
		png(fileName, width=14*scale, height=8*scale, units="in", res=200, pointsize=20)
		pushViewport(viewport(width=unit(11*scale, "inches"), height=unit(6*scale, "inches"), 
			x=unit(2*scale, "inches"), y=unit(1*scale, "inches"), default.units="inches", 
			just=c("left", "bottom")))
		vpWidth<-11*scale
		vpHeight<-6*scale
		numCol<-11
		numRow<-6
		
		lineType<-1
		lineWidth<-0.5
		lineColor<-"gray"
		for (i in seq(0,1,1/numRow))
			grid.lines(x=seq(0,1,1/numCol)*vpWidth, y=i*vpHeight, gp=gpar(lwd=lineWidth, lty=lineType, col=lineColor), default.units="inches")
		for (i in seq(0,1,1/numCol))
			grid.lines(x=i*vpWidth, y=seq(0,1,1/numRow)*vpHeight, gp=gpar(lwd=lineWidth, lty=lineType, col=lineColor), default.units="inches")
		
		# 1/23/14 add a dark line to separate plates 1 + 2 and plate 3
		grid.lines(x=seq(0,1,1/numCol)*vpWidth, y=2/6*vpHeight, gp=gpar(lwd=1, lty=lineType, col="black"), default.units="inches")

		# now need to add circles and color
		x<-seq(1/(numCol*2), 1, 1/numCol)
		y<-seq(1/(numRow*2), 1, 1/numRow)
		allX<-c(x,x,x,x,x,x)
		# start at top and go down from left to right
		allY<-c(rep(y[6],11), rep(y[5],11), rep(y[4],11), rep(y[3],11), rep(y[2],11), rep(y[1],11))
		
		# before setting the color set FC to 1 for those modules that do not meet with sig cutoff
		if (sigCutoff < 1)
		{
			for (i in 1:nrow(modFC))
			{
				curMod<-rownames(modFC)[i]
				curFDR<-modLevRes$fdr[which(curMod==rownames(modLevRes))]
				if (curFDR > sigCutoff)
					modFC$FC[i]<-1
			}
		}
		
		setFillAlpha<-function(curVals)
		{
			allFill<-c()
			allAlpha<-c()
			allRadius<-c()
			for (i in 1:length(curVals))
			{
				allRadius<-c(allRadius, 1)
				if (curVals[i] < 0)
				{
					allFill<-c(allFill, "blue")
					allAlpha<-c(allAlpha, abs(curVals[i])/maxVal)
				}
				if (curVals[i] > 0)
				{
					allFill<-c(allFill, "red")
					allAlpha<-c(allAlpha, curVals[i]/maxVal)
				}
				if (curVals[i]==0)
				{
					allFill<-c(allFill, "transparent")
					allAlpha<-c(allAlpha, 1)
				}
			}
			if (any(allAlpha > 1))
				allAlpha[which(allAlpha > 1)]<-1
			return(list(allFill=allFill, allAlpha=allAlpha, allRadius=allRadius))
		}
		retList<-setFillAlpha(curVals=log2(modFC$FC))
		allFill<-retList$allFill
		allAlpha<-retList$allAlpha	
		allRadius<-retList$allRadius
		
		grid.circle(allX*vpWidth, allY*vpHeight, r=1/(numRow*2.4)*allRadius*vpHeight, gp=gpar(col="white", fill=allFill, alpha=allAlpha, lwd=0), default.units="inches")

		popViewport()
		pushViewport(viewport(width=unit(0.5*scale, "inches"), height=unit(6*scale, "inches"), 
			x=unit(1*scale, "inches"), y=unit(4*scale, "inches"), default.units="inches", 
			just=c("left", "bottom")))
			
		# need to add a key for the color

		# now loop through the color key
		curVal<- -0.5
		# blue values
		for (i in 100:1)
		{
			grid.rect(x=0, y=0+curVal, width=1, height=1/201, gp=gpar(col="transparent", fill="blue", alpha=i/100))
			curVal<-curVal+1/201
		}
		# now set the zero value
		grid.rect(x=0, y=0+curVal, width=1, height=1/201, gp=gpar(col="transparent", fill="white"))
		curVal<-curVal+1/201
		# red values
		for (i in 1:100)
		{
			grid.rect(x=0, y=0+curVal, width=1, height=1/201, gp=gpar(col="transparent", fill="red", alpha=i/100))
			curVal<-curVal+1/201
		}		
		
		grid.rect(x=0,y=0,width=1, height=1, gp=gpar(fill="transparent"))

		# now add text
		curCex<-0.35
		grid.text(label=1/(2^maxVal), x=1.1, y=-0.5, gp=gpar(cex=curCex))
		grid.text(label=1/(2^(maxVal/2)), x=1, y=-0.5+0.25, gp=gpar(cex=curCex))
		grid.text(label=1, x=1, y=-0.5+0.5, gp=gpar(cex=curCex))
		grid.text(label=2^(maxVal/2), x=1, y=-0.5+0.75, gp=gpar(cex=curCex))
		grid.text(label=2^(maxVal), x=1, y=-0.5+1, gp=gpar(cex=curCex))
				
		grid.text(label="FC", x=0, y=0.55, gp=gpar(cex=0.4))
		
		dev.off()
	}
}




# make the case vs. control grid of 66 modules
# here there are 3 options: sig up, sig down, not sig for each module
make6x11gridModuleLevel<-function(modLevelResults, fileName, scale=0.33, maxVal=1)
{
	if (nrow(modLevelResults)==66)
	{	
		png(fileName, width=14*scale, height=8*scale, units="in", res=200, pointsize=20)
		pushViewport(viewport(width=unit(11*scale, "inches"), height=unit(6*scale, "inches"), 
			x=unit(2*scale, "inches"), y=unit(1*scale, "inches"), default.units="inches", 
			just=c("left", "bottom")))
		vpWidth<-11*scale
		vpHeight<-6*scale
		numCol<-11
		numRow<-6
		
		lineType<-1
		lineWidth<-0.5
		lineColor<-"gray"
		for (i in seq(0,1,1/numRow))
			grid.lines(x=seq(0,1,1/numCol)*vpWidth, y=i*vpHeight, gp=gpar(lwd=lineWidth, lty=lineType, col=lineColor), default.units="inches")
		for (i in seq(0,1,1/numCol))
			grid.lines(x=i*vpWidth, y=seq(0,1,1/numRow)*vpHeight, gp=gpar(lwd=lineWidth, lty=lineType, col=lineColor), default.units="inches")
		
		# 1/23/14 add a dark line to separate plates 1 + 2 and plate 3
		grid.lines(x=seq(0,1,1/numCol)*vpWidth, y=2/6*vpHeight, gp=gpar(lwd=1, lty=lineType, col="black"), default.units="inches")

		# now need to add circles and color
		x<-seq(1/(numCol*2), 1, 1/numCol)
		y<-seq(1/(numRow*2), 1, 1/numRow)
		allX<-c(x,x,x,x,x,x)
		# start at top and go down from left to right
		allY<-c(rep(y[6],11), rep(y[5],11), rep(y[4],11), rep(y[3],11), rep(y[2],11), rep(y[1],11))
		
		setFillAlpha<-function(curVals)
		{
			allFill<-c()
			allAlpha<-c()
			allRadius<-c()
			for (i in 1:length(curVals))
			{
				allRadius<-c(allRadius, 1)
				if (curVals[i] < 0)
				{
					allFill<-c(allFill, "blue")
					allAlpha<-c(allAlpha, abs(curVals[i])/maxVal)
				}
				if (curVals[i] > 0)
				{
					allFill<-c(allFill, "red")
					allAlpha<-c(allAlpha, curVals[i]/maxVal)
				}
				if (curVals[i]==0)
				{
					allFill<-c(allFill, "transparent")
					allAlpha<-c(allAlpha, 1)
				}
			}
			if (any(allAlpha > 1))
				allAlpha[which(allAlpha > 1)]<-1
			return(list(allFill=allFill, allAlpha=allAlpha, allRadius=allRadius))
		}
		retList<-setFillAlpha(curVals=modLevelResults$sig)
		allFill<-retList$allFill
		allAlpha<-retList$allAlpha	
		allRadius<-retList$allRadius
		
		grid.circle(allX*vpWidth, allY*vpHeight, r=1/(numRow*2.4)*allRadius*vpHeight, gp=gpar(col="white", fill=allFill, alpha=allAlpha, lwd=0), default.units="inches")

		popViewport()
		pushViewport(viewport(width=unit(0.5*scale, "inches"), height=unit(6*scale, "inches"), 
			x=unit(1*scale, "inches"), y=unit(4*scale, "inches"), default.units="inches", 
			just=c("left", "bottom")))
					
		dev.off()
	}
}






# show where the modules are on the grid
make6x11gridkey<-function(modFC, fileName, scale=0.33)
{
	if (nrow(modFC)==66)
	{
		png(fileName, width=14*scale, height=8*scale, units="in", res=200, pointsize=20)
		pushViewport(viewport(width=unit(11*scale, "inches"), height=unit(6*scale, "inches"), 
			x=unit(2*scale, "inches"), y=unit(1*scale, "inches"), default.units="inches", 
			just=c("left", "bottom")))
		vpWidth<-11*scale
		vpHeight<-6*scale
		numCol<-11
		numRow<-6
		
		lineType<-1
		lineWidth<-0.5
		lineColor<-"gray"
		for (i in seq(0,1,1/numRow))
			grid.lines(x=seq(0,1,1/numCol)*vpWidth, y=i*vpHeight, gp=gpar(lwd=lineWidth, lty=lineType, col=lineColor), default.units="inches")
		for (i in seq(0,1,1/numCol))
			grid.lines(x=i*vpWidth, y=seq(0,1,1/numRow)*vpHeight, gp=gpar(lwd=lineWidth, lty=lineType, col=lineColor), default.units="inches")

		# 1/23/14 add a dark line to separate plates 1 + 2 and plate 3
		grid.lines(x=seq(0,1,1/numCol)*vpWidth, y=2/6*vpHeight, gp=gpar(lwd=1, lty=lineType, col="black"), default.units="inches")

		# now need to add circles and color
		x<-seq(1/(numCol*2), 1, 1/numCol)
		y<-seq(1/(numRow*2), 1, 1/numRow)
		allX<-c(x,x,x,x,x,x)
		# start at top and go down from left to right
		allY<-c(rep(y[6],11), rep(y[5],11), rep(y[4],11), rep(y[3],11), rep(y[2],11), rep(y[1],11))
		
		grid.text(label=rownames(modFC), x=allX, y=allY, gp=gpar(cex=0.25))

		dev.off()
	}
}

makeClusteredHeatmapWithSampleInfo<-function(curData, fileName, curBreaks, clusterSams=FALSE, clusterMods=TRUE, sampleInfo, allDesign, includeControl=TRUE, modLevelResults, fdr=1)
{
	# first see if we should remove any samples (based on whether to include control)
	if (includeControl==FALSE)
	{
		remIndex<-which(allDesign$groupID==0)
		# need to remove these from the expression data (curData) and from sample information (sampleInfo)
		curData<-curData[-remIndex,]
		if (ncol(sampleInfo) > 1)
			sampleInfo<-sampleInfo[-remIndex,]
		else
			sampleInfo<-as.data.frame(sampleInfo[-remIndex,])
	}
	# then see if we need to remove any modules (based on fdr in modLevelResults)
	continue<-TRUE
	if (fdr < 1)
	{
		if (fdr > 1 | fdr < 0)
			fdr<-0.05
		keepMods<-rownames(modLevelResults)[which(modLevelResults$fdr < fdr)]
		if (length(keepMods)==0)
		{
			continue<-FALSE
			print("There were no modules that met that level of fdr so no plot was created.")
		}
		else
		{
			if (length(keepMods)==1)
			{
				continue<-FALSE
				print("Only one module was significant at that level of fdr so no plot was created.")
			}
			else
			{
				continue<-TRUE
				# now subset to those modules
				matchIndex<-match(keepMods, colnames(curData))
				curData<-curData[, matchIndex]
			}
		}
	}
		
	if (continue==TRUE)
	{
		# cluster based on samples
		samdd<-as.dendrogram(hclust(dist(curData)))	
		samdd_ord<-order.dendrogram(samdd)
	
		# cluster based on modules
		moddd<-as.dendrogram(hclust(dist(t(curData))))
		moddd_ord<-order.dendrogram(moddd)

		canvas_size = setDimensions(nrow(curData), ncol(curData))
		#print(canvas_size)

		color_palette = colorRampPalette(c("blue", "white", "red"))(1000)
		color_palette[c(495:505)] = "#FFFFFF"
	
		# assume we want to color with only one column from sampleInfo
		palette_fun<-setPalettes(ncol(sampleInfo))
	
		factors = unlist(apply(sampleInfo, 2, function(x) as.data.frame(factor(x))), recursive=F)
  		factor_levels = lapply(factors, function(x){ as.character(levels(x)) }) 
  		factor_lengths = lapply(factor_levels, function(x) length(x))

		# to find out more about add parameter in dendrogramGrob look at gpar
		# fill sets the fill color, col sets the color of borders, and alpha (0 to 1) sets transparency
		color_gen<-list()
		add_list<-list()
		for (i in 1:length(factor_levels))
		{
			tempColor<-palette_fun[[i]](as.numeric(factor_lengths[i]))
			color_gen[[i]]<-tempColor
			add_list[[i]]<-list(rect=list(col="transparent", fill=tempColor[factors[[i]]]))
		}
		add_list<-unlist(add_list, recursive=F)

		curKey<-list(space="bottom", border=TRUE, rep=FALSE, cex=2.0)	
		for (i in 1:length(factor_levels))
		{
			curKeyLen<-length(curKey)
			curKey[[curKeyLen+1]]<-list(factor_levels[[i]])
			curKey[[curKeyLen+2]]<-list(pch=15, col=color_gen[[i]])
			names(curKey)[(curKeyLen+1):(curKeyLen+2)]<-c("text","points")
		}
	
		# set any values greater than curBreaks to max value
		curData<-t(apply(curData, 1, function(x)
		{
			if (any(abs(x) > curBreaks, na.rm=T))
			{
				if (any(x > curBreaks, na.rm=T))
					x[which(x > curBreaks)]<-curBreaks
				if (any(x < -curBreaks, na.rm=T))
					x[which(x < -curBreaks)]<- -curBreaks
			}
			return(x)
		}))
	
		# now create the plots - only show dendrogram for the samples (too messy for modules)
		# first samples only
		png(fileName, width=canvas_size[[1]], height=canvas_size[[2]], units="in", res=200, pointsize=20)
	    
		# color key is the blue/red continuous color rectangle
		# legend is at top that shows the dendrogram and the independent variable coloring
		# key is the left side grob that shows the variable coloring
		cexVal=1.5
		if (clusterSams==TRUE & clusterMods==FALSE)
			lp<-levelplot(curData[samdd_ord,], xlab = NULL, ylab = NULL, 
				at = do.breaks(c(-curBreaks, curBreaks), 100), panel = panel.corrgram.2, 
				scales = list(x = list(rot = 90, cex=cexVal), y=list(cex=cexVal)), 
				colorkey = list(space = "right", labels=list(cex=cexVal)), col.regions = color_palette) 
		if (clusterSams==TRUE & clusterMods==TRUE)
		{
#			lp<-levelplot(curData[samdd_ord,moddd_ord], xlab = NULL, ylab = NULL, 
#				at = do.breaks(c(-curBreaks, curBreaks), 100), panel = panel.corrgram.2, 
#				scales = list(x = list(rot = 90, cex=cexVal), y=list(cex=cexVal)), 
#				colorkey = list(space = "right", labels=list(cex=cexVal)), col.regions = color_palette) 

			lp<-levelplot(curData[samdd_ord,moddd_ord], xlab = NULL, ylab = NULL, 
				at = do.breaks(c(-curBreaks, curBreaks), 100), panel = panel.corrgram.2,  
				scales = list(x = list(rot = 90, cex=cexVal), y=list(cex=cexVal)), key=curKey, 
				legend = list(top = list(fun = dendrogramGrob, args = list(x = samdd, ord = samdd_ord, 
				side = "top", size = 5, size.add = 0.5, add=add_list, type = "rectangle"))), 
				colorkey = list(space = "right", labels=list(cex=cexVal)), col.regions = color_palette, 
				par.settings=list(axis.components=list(right=list(tck=0), left=list(tck=0), bottom=list(tck=0), 
				top=list(tck=0)))) 
		}
		if (clusterSams==FALSE & clusterMods==TRUE)
			lp<-levelplot(curData[,moddd_ord], xlab = NULL, ylab = NULL, 
				at = do.breaks(c(-curBreaks, curBreaks), 100), panel = panel.corrgram.2, 
				scales = list(x = list(rot = 90, cex=cexVal), y=list(cex=cexVal)), 
				colorkey = list(space = "right", labels=list(cex=cexVal)), col.regions = color_palette) 
		if (clusterSams==FALSE & clusterMods==FALSE)
			lp<-levelplot(curData, xlab = NULL, ylab = NULL, 
				at = do.breaks(c(-curBreaks, curBreaks), 100), panel = panel.corrgram.2, 
				scales = list(x = list(rot = 90, cex=cexVal), y=list(cex=cexVal)), 
				colorkey = list(space = "right", labels=list(cex=cexVal)), col.regions = color_palette) 
		print(lp)
		dev.off()			
		return(samdd)
	}
}

setDimensions = function(x, y) 
{
	width = x * 0.22 + 5.57
	height = y * 0.22 + 5.57

	return(list(width, height))
} 

setPalettes = function(n) {
  # these can be seen in ?brewer.pal
  pal = list("PuOr", "Spectral", "PiYG", "BrBG", "RdGy", "PRGn", "RdYlGn")
  sel = pal[1:n]
  lapply(sel, getColor)
}

panel.corrgram.2 <- function(x, y, z, subscripts, at = pretty(z), scale = 0.8, ...) 
{ 
	# x and y are the integer locations for the circle centers
	x <- as.numeric(x)[subscripts]
	y <- as.numeric(y)[subscripts]
	
	# z is the numeric value to indicate what color the circle should be (percent positive, percent negative, zscore, etc.)
	z <- as.numeric(z)[subscripts]
	zcol <- level.colors(z, at = at, ...)
	
	# do not need to loop through the values
	grid.circle(x = x, y = y, r = .5 * scale, default.units = "native", gp = gpar(fill = zcol, col="white"))
} 

getColor = function(x) {
 colorRampPalette(brewer.pal(7, x), interpolate = "spline")
}
