######
# function from MAT to load packages
######
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
	# then load libraries
	for (i in 1:length(x))
	{
		require(x[i], character.only=TRUE)
	}
}


######
# same function used in MAT
######
returnParameter<-function(parameters, key, type, split=T)
{
	matchIndex<-which(toupper(trimString(parameters$Key))==toupper(trimString(key)))
	curKey<-parameters$Value[matchIndex]
	# now change the data type to return
	if (type=="numeric")
		output<-as.numeric(as.character(curKey))
	if (type=="character")
		output<-as.character(curKey)
	if (type=="statement")
	{
		output = eval(parse(text=as.character(curKey)))
    	if (split) {
      		output = trimString(splitString(output))
    	}
	}
	return(output)
}

#####
# trim any leading or trailing spaces in a string when doing a comparison
#####
trimString = function (string) {
  gsub("^\\s+|\\s+$", "", string)
}

######
# 08/15/12 mainly used to take the place of combineDataFromMultiplePlates function (so rest of code works)
# before was actually combine the data from multiple files into one matrix, now just getting the max assay number
######
prepareDataForPalx<-function(ctVal, ctCall)
{
	assayVals<-rownames(ctVal)
	maxAssay<-max(as.numeric(substr(unique(unlist(lapply(strsplit(assayVals, "_"), function(x) {x[1]}))), 2, 2)))
	return(list(ctValDataForPalx=ctVal, ctCallDataForPalx=ctCall, maxAssay=maxAssay))
}

######
# make matrix of Ct.Call, Ct.Quality, and Ct.Value - change format into a matrix of samples (Sample.Name) by assays (EvaGreen.Name)
# fluidData = the original fluidigm data (where each gene and sample has one row)
# colName = the column name to convert to a matrix of samples (columns) by assays (rows)
# note: the output will have chamber IDs rather than sample names because the sample name is often not unique (also well 
#  IDs for assays since they may also not be unique)
######
changeFluidFormat<-function(fluidData, colName="Ct.Value")
{
	# find number of rows and columns
	numRow<-length(unique(fluidData$AssayID))
	numCol<-length(unique(fluidData$SampleID))
	curRowNames<-as.character(unique(fluidData$AssayID))
	curColNames<-as.character(unique(fluidData$SampleID))
	
	colIndex<-match(colName, colnames(fluidData))
	
	retMatrix<-c()
	for (i in 1:length(curColNames))
	{
		tempCol<-curColNames[i]
		matchIndex<-which(tempCol==fluidData$SampleID)
		orderIndex<-match(fluidData$AssayID[matchIndex], curRowNames)
		if (colName=="Ct.Call")
			retMatrix<-cbind(retMatrix, as.character(fluidData[matchIndex[orderIndex],colIndex]))		
		else
			retMatrix<-cbind(retMatrix, fluidData[matchIndex[orderIndex],colIndex])
	}
	colnames(retMatrix)<-curColNames
	rownames(retMatrix)<-curRowNames
	
	# return a character or numeric matrix
	return(retMatrix)
}

######
# change pass criteria based on different quality number (everything below or equal to cutoff will be fail)
# qualMatrix = the Ct quality matrix (samples by genes)
# cutoff = the threshold to use (default is 0.65); the value must be between 0 and 1
######
changePassCriteria<-function(qualMatrix, cutoff=0.65)
{
	newCallMatrix<-apply(qualMatrix, 1, function(x)
	{
		y<-rep("Pass", length(x))
		y[which(as.numeric(x) <= cutoff)]<-"Fail"
		return(y)
	})	
	return(t(newCallMatrix))
}

######
# PALX on Ct Call matrix (look at genes that are not valid)
# ctCallMatrix = the matrix of "Pass" and "Fail" or "No Call" values (samples by genes)
# ctValueMatrix = the matrix of Ct values (samples by genes)
# percentage = what level of PALX should be performed (default is 10%)
# note: "No Call" is treated as a failure
######
performPALXGenes<-function(ctCallMatrix, ctValueMatrix, percentage=0.1) 
{
	passCounts<-apply(ctCallMatrix, 1, function(x)
	{
		sum(x=="Pass")
	})
	keepIndex<-which(passCounts >= ceiling(percentage*ncol(ctCallMatrix)))
	return(list(ctValueMatrix[keepIndex,], ctCallMatrix[keepIndex,]))
}
performPALXSamples<-function(ctCallMatrix, ctValueMatrix, percentage=0.1)
{
	passCounts<-apply(ctCallMatrix, 2, function(x)
	{
		sum(x=="Pass")
	})
	keepIndex<-which(passCounts >= ceiling(percentage*nrow(ctCallMatrix)))
	return(list(ctValueMatrix[, keepIndex], ctCallMatrix[, keepIndex]))
}
performPALXBoth<-function(ctCallMatrix, ctValueMatrix, percentage=0.1)
{
	passCounts<-apply(ctCallMatrix, 1, function(x)
	{
		sum(x=="Pass")
	})
	keepIndex1<-which(passCounts >= ceiling(percentage*ncol(ctCallMatrix)))
	
	passCounts<-apply(ctCallMatrix, 2, function(x)
	{
		sum(x=="Pass")
	})
	keepIndex2<-which(passCounts >= ceiling(percentage*nrow(ctCallMatrix)))
	return(list(ctValueMatrix[keepIndex1, keepIndex2], ctCallMatrix[keepIndex1, keepIndex2]))
}

######
# convert 999 values in ct values to NA since they are not representative values
# ctValueMatrix = the matrix of Ct values (rows are genes and columns are samples)
# floorValue = the value to floor 999 values to (can be NA or numeric)
######
convert999<-function(ctValueMatrix, floorValue)
{
	newCtVal<-apply(ctValueMatrix, 1, function(x)
	{
		x[which(x==999)]<-floorValue
		return(x)
	})
	return(t(newCtVal))
}

######
# set Ct values to NA when Ct Call = "Fail" and Ct value is not 999
# ctValueMatrix = the matrix of Ct values (rows are genes and columns are samples)
# ctCallMatrix = the matrix of Ct call values, which are "Pass" or "Fail" or "No Call" (rows are genes and columns are samples)
# note: No Call is treated as a failure
######
failToNA<-function(ctValueMatrix, ctCallMatrix)
{
	for (i in 1:ncol(ctCallMatrix))
	{
		callCheck<-(ctCallMatrix[,i]!="Pass")
		valCheck<-(ctValueMatrix[,i]!="999")
		if (any(is.na(valCheck)))
			valCheck[which(is.na(valCheck))]<-TRUE	
		if (any(callCheck & valCheck))
		{
			setIndex<-which(callCheck & valCheck)
			ctValueMatrix[setIndex,i]<-NA
		}
	}
	return(ctValueMatrix)
}

######
# calculate delta Ct using HK genes
# ctValueMatrix = the matrix of Ct values (rows are genes and columns are samples)
# HKgenes = character vector of HK gene well IDs
######
calculateDeltaCt<-function(ctValueMatrix, HKgenes)
{
	matchIndex<-c()
	for (i in 1:length(HKgenes))
	{
		matchIndex<-c(matchIndex, grep(HKgenes[i], rownames(ctValueMatrix)))
	}
	# now split based on plate
	curAssays<-sort(unique(unlist(lapply(strsplit(rownames(ctValueMatrix)[matchIndex], "_"), function(x){x[1]}))))
	deltaCtList<-list()
	for (i in 1:length(curAssays))
	{
		tempAssay<-curAssays[i]
		rowIndex<-grep(paste(tempAssay,"_",sep=""), rownames(ctValueMatrix))
		tempMatrix<-ctValueMatrix[rowIndex,]
		#curHKgenes<-paste(tempAssay, "_", HKgenes, sep="")
		curHKgenes<-HKgenes[grep(paste(tempAssay, "_", sep=""), HKgenes)]
		matchIndex<-match(curHKgenes, rownames(tempMatrix))
		HKgeneMatrix<-tempMatrix[matchIndex,]
		HKprod<-apply(HKgeneMatrix, 2, function(x)
		{
			prod(as.numeric(as.character(x)), na.rm=T)
		})
		allLen<-apply(HKgeneMatrix, 2, function(x) {sum(!is.na(x))})
		HKgeoMean<-HKprod^(1/allLen)
		tempDeltaCt<-t(apply(tempMatrix, 1, function(x)
		{
			as.numeric(as.character(x)) - HKgeoMean
		}))
		# do not remove HK genes for now
		deltaCtList[[i]]<-tempDeltaCt
	}
	# now combine all list elements together (HK genes have already been removed)
	deltaCt<-c()
	for (i in 1:length(deltaCtList))
	{
		deltaCt<-rbind(deltaCt, deltaCtList[[i]])
	}
	return(deltaCt)
}

######
# calculate delta delta Ct when there is a reference group (e.g. genomic CBC)
# deltaCt = the delta Ct matrix created in calculateDeltaCt (rows are genes/assays and columns are samples) - values can be positive or negative
# sampleDesignData = the sample design data (contains sample IDs, sample names, and type)
######
calculateDeltaDeltaCtWithRefGroup<-function(deltaCt, sampleDesignData)
{
	# need to convert the reference samples to chamber IDs
	# also need to split into multiple plates
	uniPlates<-sort(unique(unlist(lapply(strsplit(colnames(deltaCt), "_"), function(x) {x[1]}))))
	allDeltaDeltaCt<-c()
	for (i in 1:length(uniPlates))
	{
		# trim leading and trailing spaces from ref sample names
		refSampleNames<-unique(as.character(sampleDesignData$SampleName[which(sampleDesignData$Type=="Ref")]))
		
		allChIDmatrix<-unique(sampleDesignData[which(sampleDesignData$Type=="Ref"),1:2])
		curChIDmatrix<-allChIDmatrix[grep(paste(uniPlates[i], "_", sep=""),allChIDmatrix$SampleID),]
		
		# need to take average of multiple wells of same sample (will not know ahead of time how many wells)
		refMatrix<-c()
		for (j in 1:length(refSampleNames))
		{
			curChIDs<-as.character(curChIDmatrix[which(refSampleNames[j]==curChIDmatrix[,2]),1])
			tempDeltaCt<-deltaCt[,match(curChIDs, colnames(deltaCt))]
			# now take average
			if (length(curChIDs)>1)
				tempDeltaCt<-apply(tempDeltaCt, 1, mean, na.rm=T)
			refMatrix<-cbind(refMatrix, tempDeltaCt)
		}
		colnames(refMatrix)<-refSampleNames
		# now get the average of these values
		refMeans<-apply(refMatrix, 1, mean, na.rm=T)
		# now use on non-reference samples for this plate
		curPlateDeltaCt<-deltaCt[,grep(paste(uniPlates[i],"_",sep=""), colnames(deltaCt))]
		tempDeltaDeltaCt<-apply(curPlateDeltaCt, 2, function(x)
		{
			x - refMeans
		})
		allDeltaDeltaCt<-cbind(allDeltaDeltaCt, tempDeltaDeltaCt)
	}
	# do not remove reference samples
	return(allDeltaDeltaCt)
}


