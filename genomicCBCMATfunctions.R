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
		curTT$module<-moduleAssignData$module[match(curTT$ID, moduleAssignData$geneID)]
		curTT$geneName<-moduleAssignData$geneName[match(curTT$ID, moduleAssignData$geneID)]
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
