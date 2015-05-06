######
# function to read in data and set up the column names for the Fluidigm data
# fileName = the directory plus file name of the Fluidigm csv file
# skipRows = the number of rows to skip in the csv file (standard file has 10 rows that should be skipped)
######
readFluidData<-function(fileName, skipRows=10)
{
	# read in data
	fluidData<-read.csv(file=fileName, skip=skipRows)
	# now alter the column names
	oldCN<-colnames(fluidData)
	rem1<-grep(".1", oldCN)
	for (i in 1:length(rem1))
	{
		oldCN[rem1[i]]<-substr(oldCN[rem1[i]], 1, nchar(oldCN[rem1[i]])-2)
	}
	rem2<-grep(".2", oldCN)
	for (i in 1:length(rem1))
	{
		oldCN[rem2[i]]<-substr(oldCN[rem2[i]], 1, nchar(oldCN[rem2[i]])-2)
	}
	rem3<-grep(".3", oldCN)
	for (i in 1:length(rem1))
	{
		oldCN[rem3[i]]<-substr(oldCN[rem3[i]], 1, nchar(oldCN[rem3[i]])-2)
	}
	newCN<-paste(oldCN, as.matrix(fluidData[1,]), sep=".")
	newCN<-gsub(" ", ".", newCN)
	colnames(fluidData)<-newCN
	fluidData<-fluidData[2:nrow(fluidData),]
	return(fluidData)
}

######
# 5/20/13
# reorder file names so that they are ordered based on the second element in the file name (S1A1, S1A2, S2A1, ...)
######
reorderFileNames<-function(fileNames)
{
	# sample assay names should be in the second position (after the first underscore)
	sampleAssayNames<-unlist(lapply(strsplit(fileNames, "_"), function(x) {x[[2]]}))
	# can not just order using these names because could have S1 and S10
	samNum<-unlist(lapply(strsplit(sampleAssayNames, "A"), function(x) {x[[1]]}))
	samNum<-as.numeric(substr(samNum, 2, nchar(samNum)))
	assayNum<-as.numeric(unlist(lapply(strsplit(sampleAssayNames, "A"), function(x) {x[[2]]})))
	newFileNames<-fileNames[order(samNum, assayNum)]
	return(newFileNames)
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
# 1/3/13
# check that samples and genes are in the same chamber IDs across plates
# this function assumes we are using 96 by 96 well plates!!!
######
checkSampleGeneChamberIDs<-function(dataList)
{
	anyErrors<-FALSE
	if (length(dataList) > 1)
	{
		# need the names to know which comparisons to make
		dataNames<-names(dataList)
		plateInfo<-unlist(lapply(strsplit(dataNames, "_"), function(x) {x[2]}))
		sampleNumber<-unlist(lapply(strsplit(plateInfo, "A"), function(x) {x[1]}))
		assayNumber<-unlist(lapply(strsplit(plateInfo, "A"), function(x) {x[2]}))
		
		# first check samples across plates 
		sampleResults<-c()
		uniSamples<-unique(sampleNumber)
		for (i in 1:length(uniSamples))
		{
			curIndex<-which(sampleNumber==uniSamples[i])
			if (length(curIndex) > 1)
			{
				for (j in 2:length(curIndex))
				{
					# this assumes we are using a 96 by 96 well plate!!
					sampleResults<-c(sampleResults, all(as.character(dataList[[curIndex[j-1]]]$Sample.Name[seq(1,9216,96)])==as.character(dataList[[curIndex[j]]]$Sample.Name[seq(1,9216,96)])))
				}
			}
		}
		if (any(sampleResults==FALSE))
		{
			print("ERROR: sample names do not match across some plates. Will stop processing.")
			anyErrors<-TRUE
		}
				
		# next check genes within one plate
		geneResultsOnePlate<-c()
		for (i in 1:length(dataList))
		{
			geneResultsOnePlate<-c(geneResultsOnePlate, all(rep(dataList[[i]]$EvaGreen.Name[1:96], 96)==dataList[[i]]$EvaGreen.Name))
		}
		if (any(geneResultsOnePlate==FALSE))
		{
			print("ERROR: gene names do not match on same plate. Will stop processing.")
			anyErrors<-TRUE
		}	
			
		# next check genes across plates
		geneResults<-c()
		uniAssays<-unique(assayNumber)
		for (i in 1:length(uniAssays))
		{
			curIndex<-which(assayNumber==uniAssays[i])
			if (length(curIndex) > 1)
			{
				for (j in 2:length(curIndex))
				{
					geneResults<-c(geneResults, all(as.character(dataList[[curIndex[j-1]]]$EvaGreen.Name)==as.character(dataList[[curIndex[j]]]$EvaGreen.Name)))
				}
			}
		}
		if (any(geneResults==FALSE))
		{
			print("ERROR: gene names do not match across some plates. Will stop processing.")
			anyErrors<-TRUE
		}
	}
	if (anyErrors==FALSE)
		print("Gene and sample names matched across multiple plates.")
	return(anyErrors)
}

######
# combine Ct call and value data from multiple plates (may have multiple sample and assay plates)
# dataList = a list of Fluidigm data from multiple plates
# ctValList = a list of Ct value matrices (from multiple plates)
# ctCallList = a list of Ct call matrices (from multiple plates)
######
combineDataFromMultiplePlates<-function(dataList, ctValList, ctCallList)
{
	# need to determine how to combine using the file names
	sampleAssay<-unlist(lapply(strsplit(names(dataList), "_"), function(x) {x[2]}))	
	numSam<-unlist(lapply(strsplit(unlist(lapply(strsplit(substr(sampleAssay, 2, nchar(sampleAssay)),"S"), function(x) {x[1]})), "A"), function(x) {x[[1]]}))
	numAssay<-unlist(lapply(strsplit(substr(sampleAssay, 2, nchar(sampleAssay)),"A"), function(x) {x[2]}))	
	maxSam<-max(as.numeric(numSam))
	maxAssay<-max(as.numeric(numAssay))
	minSam<-min(as.numeric(numSam))
	minAssay<-min(as.numeric(numAssay))

	ctValDataForPalx<-c()
	ctCallDataForPalx<-c()
	for (i in minSam:maxSam)
	{
		for (j in minAssay:maxAssay)
		{
			curSamAssay<-paste("S",i,"A",j,sep="")
			curIndex<-grep(curSamAssay, names(dataList))
			if (j==minAssay)
			{
				tempVal<-ctValList[[curIndex]]
				tempCall<-ctCallList[[curIndex]]
				rownames(tempVal)<-paste("A",j,"_",rownames(tempVal),sep="")
				rownames(tempCall)<-paste("A",j,"_",rownames(tempCall),sep="")
			}	
			else
			{
				curVal<-ctValList[[curIndex]]
				curCall<-ctCallList[[curIndex]]
				rownames(curVal)<-paste("A",j,"_",rownames(curVal),sep="")
				rownames(curCall)<-paste("A",j,"_",rownames(curCall),sep="")
				tempVal<-rbind(tempVal, curVal)
				tempCall<-rbind(tempCall, curCall)
			}	
		}
		colnames(tempVal)<-paste("S",i,"_",colnames(tempVal),sep="")
		colnames(tempCall)<-paste("S",i,"_",colnames(tempCall),sep="")
		ctValDataForPalx<-cbind(ctValDataForPalx, tempVal)
		ctCallDataForPalx<-cbind(ctCallDataForPalx, tempCall)
	}
	return(list(ctValDataForPalx=ctValDataForPalx, ctCallDataForPalx=ctCallDataForPalx, maxAssay=maxAssay))
}

######
# make matrix of Ct.Call, Ct.Quality, and Ct.Value - change format into a matrix of samples (Sample.Name) by genes (EvaGreen.Name)
# fluidData = the original fluidigm data (where each gene and sample has one row)
# colName = the column name to convert to a matrix of samples (columns) by genes (rows)
# note: the output will have chamber IDs rather than sample names because the sample name is often not unique
######
changeFluidFormat<-function(fluidData, colName="Ct.Value")
{
	# problem is that sample name is often not unique - need to use chamber.id
	curChamber<-as.character(fluidData$Chamber.ID)
	sample<-unlist(lapply(strsplit(curChamber, "-"), function(x) {x[1]}))
	uniSample<-unique(sample)
#	uniGenes<-unique(as.character(fluidData$EvaGreen.Name))
	assays<-unlist(lapply(strsplit(curChamber, "-"), function(x) {x[2]}))
	uniGenes<-unique(assays)
	colIndex<-match(colName, colnames(fluidData))

	# make a matrix of genes (rows) by samples (columns)
	newMatrix<-matrix(as.character(fluidData[,colIndex]), nrow=length(uniGenes), ncol=length(uniSample))		
	rownames(newMatrix)<-uniGenes
	colnames(newMatrix)<-uniSample
	# return a character matrix
	return(newMatrix)
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
# ctCallMatrix = the matrix of "Pass" and "Fail" values (samples by genes)
# ctValueMatrix = the matrix of Ct values (samples by genes)
# percentage = what level of PALX should be performed (default is 10%)
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
# ctCallMatrix = the matrix of Ct call values, which are "Pass" or "Fail" (rows are genes and columns are samples)
######
failToNA<-function(ctValueMatrix, ctCallMatrix, curFloor, workdir)
{
	origVals<-c()
	floorVals<-c()
	for (i in 1:ncol(ctCallMatrix))
	{
		callCheck<-(ctCallMatrix[,i]=="Fail")
		valCheck<-(ctValueMatrix[,i]!="999")
		if (any(is.na(valCheck)))
			valCheck[which(is.na(valCheck))]<-TRUE	
		if (any(callCheck & valCheck))
		{
			setIndex<-which(callCheck & valCheck)
			if (is.na(curFloor))
				ctValueMatrix[setIndex,i]<-NA
			else
			{
				# then the value should be floored rather than set to NA
				#print(i)
				# show the values that will be replaced
				#print(ctValueMatrix[setIndex,i])
				for (j in 1:length(setIndex))
				{
					origVals<-c(origVals, as.numeric(as.matrix(ctValueMatrix[setIndex[j],i])))
					if (as.numeric(ctValueMatrix[setIndex[j],i]) < 5)
					{
						ctValueMatrix[setIndex[j],i]<-5
					}
					else
					{
						ctValueMatrix[setIndex[j],i]<-curFloor
					}
					floorVals<-c(floorVals, as.numeric(as.matrix(ctValueMatrix[setIndex[j],i])))
				}
			}
		}
	}
	if (length(origVals) > 0)
	{
		pdf(file=paste(workdir,"Results/Floored_vs_Original_Values_Scatterplot.pdf",sep=""))
		plot(x=origVals, y=floorVals, xlab="Original Values", ylab="Floored Values", main="Values not equal to 999 that had a Call of Fail")
		dev.off()
	}
	return(ctValueMatrix)
}

######
# calculate delta Ct using HK genes
# ctValueMatrix = the matrix of Ct values (rows are genes and columns are samples)
# HKgenes = character vector of the housekeeping gene names
# ctCallMatrix = the matrix of call values (Pass and Fail - rows are genes and columns are samples); needed to not use fails in this normalization step
######
calculateDeltaCt<-function(ctValueMatrix, HKgenes, ctCallMatrix)
{
	matchIndex<-c()
	matchCallIndex<-c()
	for (i in 1:length(HKgenes))
	{
		matchIndex<-c(matchIndex, grep(HKgenes[i], rownames(ctValueMatrix)))
		matchCallIndex<-c(matchCallIndex, grep(HKgenes[i], rownames(ctCallMatrix)))
	}
	# now split based on plate
	curAssays<-sort(unique(unlist(lapply(strsplit(rownames(ctValueMatrix)[matchIndex], "_"), function(x){x[1]}))))
	deltaCtList<-list()
	for (i in 1:length(curAssays))
	{
		tempAssay<-curAssays[i]
		rowIndex<-grep(paste(tempAssay,"_",sep=""), rownames(ctValueMatrix))
		tempMatrix<-ctValueMatrix[rowIndex,]
		
		rowIndex<-grep(paste(tempAssay,"_",sep=""), rownames(ctCallMatrix))
		tempCallMatrix<-ctCallMatrix[rowIndex,]
		
		curHKgenes<-HKgenes[grep(paste(tempAssay, "_", sep=""), HKgenes)]
		matchIndex<-match(curHKgenes, rownames(tempMatrix))
		matchCallIndex<-match(curHKgenes, rownames(tempCallMatrix))
		
		HKgeneMatrix<-tempMatrix[matchIndex,]
		HKgeneCallMatrix<-tempCallMatrix[matchIndex,]
		# set any failed values to NA before normalizing to HK genes (will not normalize to NA unless all chosen HK genes are failures)
		# 5/21/13 this step is necessary in case the values have been floored (rather than set to NA)
		if (any(toupper(HKgeneCallMatrix)==toupper("Fail")))
		{
			# then set these values to NA before normalizing
			for (j in 1:ncol(HKgeneMatrix))
			{
				if (any(toupper(HKgeneCallMatrix[,j])==toupper("Fail")))
				{
					HKgeneMatrix[which(toupper(HKgeneCallMatrix[,j])==toupper("Fail")), j]<-NA
				}
			}
		}
		
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
# calculate delta delta Ct when there is one reference sample for each non-reference sample (e.g. TollGen, Epigen)
# deltaCt = the delta Ct matrix created in calculateDeltaCt (rows are genes and columns are samples) - values can be positive or negative
# nonRefSample = a vector of the stimulated (non-reference) samples
# refSample = a vector of the non-stimulated (reference) samples; refSample must have the same length as nonRefSample
######
calculateDeltaDeltaCt<-function(deltaCt, nonRefSample, refSample)
{
	if (length(nonRefSample) != length(refSample))
		print("The sample list lengths must be the same")
	else
	{
		matchNonRef<-match(nonRefSample, colnames(deltaCt))
		matchRef<-match(refSample, colnames(deltaCt))
		deltadeltaCt<-deltaCt[,matchNonRef] - deltaCt[,matchRef]
		return(deltadeltaCt)
	}
}

#######
# make plots of reference sample delta Ct values compared to non-reference sample delta Ct values
# ctDeltaMatrix = the matrix of delta Ct values (rows are genes and columns are samples)
# RefSampleChamberIDs = the chamber IDs of the reference samples
#######
RefvsNonRefgenesBoxplot<-function(ctDeltaMatrix, RefSampleChamberIDs)
{
	Refmat<-ctDeltaMatrix[, match(RefSampleChamberIDs, colnames(ctDeltaMatrix))]
	nonRefmat<-ctDeltaMatrix[, which(!((1:ncol(ctDeltaMatrix)) %in% match(RefSampleChamberIDs, colnames(ctDeltaMatrix))))]
	boxplot(c(as.numeric(Refmat), as.numeric(nonRefmat))~c(rep(1,length(as.numeric(Refmat))),rep(2,length(as.numeric(nonRefmat)))), axes=FALSE,
		ylab="Delta Ct")
	box()
	axis(1, at=1:2, labels=c("Ref samples", "non Ref samples"))
	axis(2)
}


######
# determine the chamberIDs for the stim and nonstim samples
# fluidData = the original fluidigm data (where each gene and sample has one row)
# sampleName = the sample names to convert to chamber IDs (can have multiple sample names)
# note: a sample name may not be unique so need to return all chamberIDs for that sample name
######
determineChamberIDs<-function(fluidData, sampleName)
{
	curChamber<-as.character(fluidData$Chamber.ID)
	sample<-unlist(lapply(strsplit(curChamber, "-"), function(x) {x[1]}))	
	sampleNames<-trimString(as.character(fluidData$Sample.Name))
	chamberNames<-unique(cbind(sample, sampleNames))
	if (nrow(chamberNames)!=96)
		print("The number of sample names is not 96.")
	# 5/10/12 - problem: the sample names may have replicates - return a matrix of sample names and chamber IDs
	sampleName<-unique(sampleName)
	
	matchIndex<-c()
	retSamNames<-c()
	for (i in 1:length(sampleName))
	{
		curIndex<-which(sampleName[i]==chamberNames[,2])
		if (length(curIndex)==1)
		{
			matchIndex<-c(matchIndex, curIndex)
			retSamNames<-c(retSamNames, sampleName[i])
		}
		else
		{
			# take all chamberIDs
			tempIndex<-which(!(curIndex %in% matchIndex))
			matchIndex<-c(matchIndex, curIndex[tempIndex])
			retSamNames<-c(retSamNames, rep(sampleName[i],length(tempIndex)))
		}
	}
	if (any(is.na(matchIndex)))
	{
		print("Some of the sample names do not have matches.")
		matchIndex<-matchIndex[which(!is.na(matchIndex))]
	}
	if (length(matchIndex)>0)
	{
		chamberIDs<-as.character(chamberNames[matchIndex,1])
		retMat<-cbind(chamberIDs, retSamNames)
		return(retMat)	
	}
	else
	{
		return(c())
	}
}

######
# use a subset of the sample name to determine actual sample name
# fluidData = the original fluidigm data (where each gene and sample has one row)
# partialName = a part of the sample name, but not the full name (can be a vector of length > 1)
######
determineSampleNamesFromPartName<-function(fluidData, partialName)
{
	uniSampleNames<-as.character(unique(fluidData$Sample.Name))
	curNames<-c()
	for (i in 1:length(partialName))
	{
		curIndex<-grep(partialName[i], uniSampleNames)
		if (length(curIndex)==0)
			print(paste("No match", i))
		else
			curNames<-c(curNames, uniSampleNames[curIndex])
	}
	return(curNames)
}

######
# use chamberID to get sample names
# fluidData = the original fluidigm data (where each gene and sample has one row)
# chamberID = the chamber ID(s)
######
determineSampleNamesFromChamberID<-function(fluidData, chamberID)
{
	chIDs<-unlist(lapply(strsplit(as.character(fluidData$Chamber.ID[seq(1,nrow(fluidData),96)]),"-"), function(x){x[1]}))
	sNames<-as.character(fluidData$Sample.Name[seq(1,nrow(fluidData),96)])
	matchIndex<-match(chamberID, chIDs)
	if (any(is.na(matchIndex)))
		matchIndex<-matchIndex[-which(is.na(matchIndex))]
	if (length(matchIndex)>0)
		return(sNames[matchIndex])
}

######
# 12/24/12
# get the sample names from the chamber IDs, but use those that include plate number
######
determineSampleNamesFromChamberIDInDataList<-function(dataList, chamberID, maxAssay)
{
	samplePlateNum<-unlist(lapply(strsplit(chamberID, "_"), function(x) {x[1]}))
	plateNum<-as.numeric(substr(samplePlateNum, 2, nchar(samplePlateNum)))
	# find the correct dataList element (need to take into account that may have more than one assay plate in a panel)
	dataListElement<-seq(1,99,maxAssay)[plateNum]
	samNames<-c()
	for (i in 1:length(dataListElement))
	{
		curIndex<-grep(paste(samplePlateNum[i],"A",sep=""), names(dataList))[1]
		samNames<-c(samNames, determineSampleNamesFromChamberID(dataList[[curIndex]], unlist(lapply(strsplit(chamberID, "_"), function(x) {x[2]}))[i]))
	}
	return(samNames)
}

######
# make a heatmap of FC - assume they want to cluster both genes and samples
# fcVals = a matrix of FC values on original scale (columns are samples and rows are genes)
# maxVal = the max log2(FC) to be shown on the plot (often a matrix of FC vals may have a few large values
# 		that would wash out all the color on the heatmap so set a reasonable max value; default is 4, which 
#		gives a FC range of 1/16 to 16)
# fileName = the directory and file name to place the png file
# curWidth = the width of the png file in inches
# curHeight = the height of the png file in inches
# aspect = this relates to the ratio of height to width of each cell square on the matrix
######
FCheatmap<-function(fcVals, maxVal=4, fileName, curWidth=8, curHeight=8, aspect=0.8, clusterRows=TRUE, clusterCols=TRUE)
{
	library(lattice)
	library(latticeExtra)
	library(grid)
	
	color_palette<-colorRampPalette(c("blue","yellow","red"))(1000)
	
	newData<-t(log2(fcVals))
	
	# use correlation for distance metric
	samCor<-as.dist(1-cor(t(newData), use="complete.obs"))
	geneCor<-as.dist(1-cor(newData, use="complete.obs"))
	samdd<-as.dendrogram(hclust(samCor))
	genedd<-as.dendrogram(hclust(geneCor))
	samdd_ord<-order.dendrogram(samdd)
	genedd_ord<-order.dendrogram(genedd)
	
	if (any(newData > maxVal))
	{
		newData[which(newData > maxVal)]<-maxVal
	}
	if (any(newData < -maxVal))
	{
		newData[which(newData < -maxVal)]<- -maxVal
	}
	
	png(fileName, width=curWidth, height=curHeight, units="in", res=200, pointsize=15)

	if (clusterRows==TRUE & clusterCols==TRUE)
	{
		lp<-levelplot(newData[samdd_ord,genedd_ord], xlab = NULL, ylab = NULL, 
				at = do.breaks(c(-maxVal, maxVal), 100),  
				scales = list(x = list(rot = 90, cex=0.8), y=list(cex=0.7)), 
				legend = list(top = list(fun = dendrogramGrob, args = list(x = samdd, ord = samdd_ord, 
				side = "top", size = 5, size.add = 0.5, type = "rectangle"))), 
				colorkey = list(space = "right", labels=list(cex=0.8)), col.regions = color_palette, aspect=aspect,
				par.settings=list(axis.components=list(right=list(tck=0), left=list(tck=0), bottom=list(tck=0), top=list(tck=0)))) 
	}
	else if (clusterRows==TRUE & clusterCols==FALSE)
	{
			lp<-levelplot(newData[,genedd_ord], xlab = NULL, ylab = NULL, 
				at = do.breaks(c(-maxVal, maxVal), 100),  
				scales = list(x = list(rot = 90, cex=0.8), y=list(cex=0.7)), 
				colorkey = list(space = "right", labels=list(cex=0.8)), col.regions = color_palette, aspect=aspect,
				par.settings=list(axis.components=list(right=list(tck=0), left=list(tck=0), bottom=list(tck=0), top=list(tck=0)))) 
	}
	else if (clusterRows==FALSE & clusterCols==TRUE)
	{
			lp<-levelplot(newData[samdd_ord,], xlab = NULL, ylab = NULL, 
				at = do.breaks(c(-maxVal, maxVal), 100),  
				scales = list(x = list(rot = 90, cex=0.8), y=list(cex=0.7)), 
				legend = list(top = list(fun = dendrogramGrob, args = list(x = samdd, ord = samdd_ord, 
				side = "top", size = 5, size.add = 0.5, type = "rectangle"))), 
				colorkey = list(space = "right", labels=list(cex=0.8)), col.regions = color_palette, aspect=aspect,
				par.settings=list(axis.components=list(right=list(tck=0), left=list(tck=0), bottom=list(tck=0), top=list(tck=0)))) 
	}
	else
	{
			lp<-levelplot(newData, xlab = NULL, ylab = NULL, 
				at = do.breaks(c(-maxVal, maxVal), 100),  
				scales = list(x = list(rot = 90, cex=0.8), y=list(cex=0.7)), 
				colorkey = list(space = "right", labels=list(cex=0.8)), col.regions = color_palette, aspect=aspect,
				par.settings=list(axis.components=list(right=list(tck=0), left=list(tck=0), bottom=list(tck=0), top=list(tck=0)))) 
	}
	print(lp)
	dev.off()
}



######
# make a heatmap of FC with sample info key - assume they want to cluster both genes and samples
# fcVals = a matrix of FC values on original scale (columns are samples and rows are genes)
# maxVal = the max log2(FC) to be shown on the plot (often a matrix of FC vals may have a few large values
# 		that would wash out all the color on the heatmap so set a reasonable max value; default is 4, which 
#		gives a FC range of 1/16 to 16)
# fileName = the directory and file name to place the png file
# curWidth = the width of the png file in inches
# curHeight = the height of the png file in inches
# aspect = this relates to the ratio of height to width of each cell square on the matrix
# sampleInfo = a matrix of variables about the samples (assumes the samples (rows) are ordered same as samples (columns) in FC matrix)
######
FCheatmapWithSamInfo<-function(fcVals, maxVal=4, fileName, curWidth=8, curHeight=8, aspect=0.8, clusterRows=TRUE, clusterCols=TRUE, 
	sampleInfo)
{
	library(lattice)
	library(latticeExtra)
	library(grid)
	
	color_palette<-colorRampPalette(c("blue","yellow","red"))(1000)
	
	newData<-t(log2(fcVals))
	
	# use correlation for distance metric
	samCor<-as.dist(1-cor(t(newData), use="complete.obs"))
	geneCor<-as.dist(1-cor(newData, use="complete.obs"))
	samdd<-as.dendrogram(hclust(samCor))
	genedd<-as.dendrogram(hclust(geneCor))
	samdd_ord<-order.dendrogram(samdd)
	genedd_ord<-order.dendrogram(genedd)
	
	if (any(newData > maxVal))
	{
		newData[which(newData > maxVal)]<-maxVal
	}
	if (any(newData < -maxVal))
	{
		newData[which(newData < -maxVal)]<- -maxVal
	}
	
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

	curKey<-list(space="bottom", border=TRUE, rep=FALSE, cex=0.7)	
	for (i in 1:length(factor_levels))
	{
		curKeyLen<-length(curKey)
		curKey[[curKeyLen+1]]<-list(factor_levels[[i]])
		curKey[[curKeyLen+2]]<-list(pch=15, col=color_gen[[i]])
		names(curKey)[(curKeyLen+1):(curKeyLen+2)]<-c("text","points")
	}
	png(fileName, width=curWidth, height=curHeight, units="in", res=200, pointsize=15)

	if (clusterRows==TRUE & clusterCols==TRUE)
	{
		lp<-levelplot(newData[samdd_ord,genedd_ord], xlab = NULL, ylab = NULL, 
				at = do.breaks(c(-maxVal, maxVal), 100),  
				scales = list(x = list(rot = 90, cex=0.8), y=list(cex=0.7)), key=curKey, 
				legend = list(top = list(fun = dendrogramGrob, args = list(x = samdd, ord = samdd_ord, 
				side = "top", size = 5, size.add = 0.5, add=add_list, type = "rectangle"))), 
				colorkey = list(space = "right", labels=list(cex=0.8)), col.regions = color_palette, aspect=aspect,
				par.settings=list(axis.components=list(right=list(tck=0), left=list(tck=0), bottom=list(tck=0), top=list(tck=0)))) 
	}
	print(lp)
	dev.off()
}


# functions from MAT
getColor = function(x) {
 colorRampPalette(brewer.pal(7, x), interpolate = "spline")
}

setPalettes = function(n) {
  # these can be seen in ?brewer.pal
  pal = list("PuOr", "Spectral", "PiYG", "BrBG", "RdGy", "PRGn", "RdYlGn")
  sel = pal[1:n]
  lapply(sel, getColor)
}



######
# remove HK genes
# curMatrix = any matrix that has genes for rows (could be deltaCt or deltadeltaCt)
# HKgenes = vector of HK gene names
######
removeHKgenes<-function(curMatrix, HKgenes)
{
	matchIndex<-match(HKgenes, rownames(curMatrix))
	curMatrix<-curMatrix[-matchIndex,]
	return(curMatrix)
}

######
# calculate delta delta Ct when there is a reference group (e.g. genomic CBC)
# deltaCt = the delta Ct matrix created in calculateDeltaCt (rows are genes and columns are samples) - values can be positive or negative
# refSampleNames = the name of the reference samples (e.g. healthy controls)
# dataList - the data list directly from the Fluidigm files
# ctCallMatrix - the matrix of call values (Pass and Fail)
######
calculateDeltaDeltaCtWithRefGroup<-function(deltaCt, refSampleNames, dataList, ctCallMatrix)
{
	# need to convert the reference samples to chamber IDs
	# also need to split into multiple plates
	uniPlates<-sort(unique(unlist(lapply(strsplit(colnames(deltaCt), "_"), function(x) {x[1]}))))
	allDeltaDeltaCt<-c()
	remChIDs<-c()
	for (i in 1:length(uniPlates))
	{
		# take the first plate (should not matter which one is picked as long as the sample order is consistent)
		curIndex<-grep(paste("_",uniPlates[i],"A",sep=""), names(dataList))[1]
		# trim leading and trailing spaces from ref sample names
		refSampleNames<-trimString(refSampleNames)
		curChIDmatrix<-determineChamberIDs(dataList[[curIndex]], refSampleNames)
		# need to take average of multiple wells of same sample (will not know ahead of time how many wells)
		refMatrix<-c()
		for (j in 1:length(refSampleNames))
		{
			curChIDs<-curChIDmatrix[which(refSampleNames[j]==curChIDmatrix[,2]),1]
			remChIDs<-c(remChIDs, paste(uniPlates[i],"_",curChIDs,sep=""))
			tempDeltaCt<-deltaCt[,match(paste(uniPlates[i],"_",curChIDs,sep=""), colnames(deltaCt))]
			# set any reference samples that are failed to NA
			tempCallMatrix<-ctCallMatrix[,match(paste(uniPlates[i],"_",curChIDs,sep=""), colnames(ctCallMatrix))]
			if (any(toupper(tempCallMatrix)==toupper("Fail")))
			{
				if (length(curChIDs) > 1)
				{
					for (k in 1:length(curChIDs))
					{
						if (any(toupper(tempCallMatrix[,k])==toupper("Fail")))
							tempDeltaCt[which(toupper(tempCallMatrix[,k])==toupper("Fail")), k]<-NA
					}
				}
				else
				{
					tempDeltaCt[which(toupper(tempCallMatrix)==toupper("Fail"))]<-NA
				}
			}
			
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
#	# then remove reference samples before returning value
#	remIndex<-match(remChIDs, colnames(allDeltaDeltaCt))
#	if (any(is.na(remIndex)))
#		remIndex<-remIndex[-which(is.na(remIndex))]
#	if (length(remIndex)>0)
#		allDeltaDeltaCt<-allDeltaDeltaCt[,-remIndex]
	return(allDeltaDeltaCt)
}

######
# change the column names from chamber ID to sample name
# curMatrix = matrix where genes are rows and samples are columns (delta Ct, delta delta Ct, FC)
# fluidData = the original fluidigm data (where each gene and sample has one row)
######
changeColumnNames<-function(curMatrix, fluidData)
{
	curCNames<-colnames(curMatrix)
	curChamber<-as.character(fluidData$Chamber.ID)
	sample<-unlist(lapply(strsplit(curChamber, "-"), function(x) {x[1]}))
	chamberSample<-unique(cbind(sample, as.character(fluidData$Sample.Name)))
	newCNames<-chamberSample[match(curCNames, chamberSample[,1]),2]
	colnames(curMatrix)<-newCNames
	return(curMatrix)
}

######
# 12/26/12
# make a conversion matrix between gene names and wells on the plate
# wellID = the well IDs for gene names
# dataList = the fluidigm data (list of matrices with 9216 rows)
######
convertWellToGeneName<-function(wellID, dataList)
{
	returnMapWellGeneName<-c()
	for (i in 1:length(dataList))
	{
		curChamber<-as.character(dataList[[i]]$Chamber.ID)
		assays<-unlist(lapply(strsplit(curChamber, "-"), function(x) {x[2]}))
		curGenes<-as.character(dataList[[i]]$EvaGreen.Name)
		
		sampleAssay<-substr(unlist(lapply(strsplit(names(dataList)[i], "_"), function(x) {x[2]})),3,4)	
		completeAssay<-paste(sampleAssay,"_", assays, sep="")
		mapGeneWell<-unique(cbind(completeAssay, curGenes))
		keepIndex<-match(wellID, mapGeneWell[,1])
		# remove any NAs
		if (any(is.na(keepIndex)))
			keepIndex<-keepIndex[-which(is.na(keepIndex))]
		if (length(keepIndex) > 0)
			returnMapWellGeneName<-rbind(returnMapWellGeneName, mapGeneWell[keepIndex,])
	}
	# get unique rows - may have repeats due to multiple sample plates
	returnMapWellGeneName<-unique(returnMapWellGeneName)
	# sort the HK genes to ensure that they are in the same order
	if (nrow(returnMapWellGeneName) > 1)
		returnMapWellGeneName<-returnMapWellGeneName[order(returnMapWellGeneName[,1]),]
	colnames(returnMapWellGeneName)<-c("WellID", "GeneName")
	return(returnMapWellGeneName)
}

#######
# 4/17/14
# make PCA plot to look for plate batch effect
#######
makePCAplot<-function(FCvals, workdir, retMap, refSamples, SampleChamberID)
{
	curCN<-colnames(FCvals)
	curPlates<-unlist(lapply(strsplit(curCN, "_"), function(x) {x[1]}))
	uniPlates<-unique(curPlates)
	numPlates<-length(uniPlates)
	
	# remove HK genes 
	if (any(retMap[,1] %in% rownames(FCvals)))
	{
		matchIndex<-match(retMap[,1], rownames(FCvals))
		if (any(is.na(matchIndex)))
			matchIndex<-matchIndex[-which(is.na(matchIndex))]
		FCvals<-FCvals[-matchIndex,]
	}
	# remove ref samples
	refSampleNames<-trimString(refSamples)
	remRefs<-c()
	for (i in 1:length(refSampleNames))
	{
		remRefs<-c(remRefs, as.character(SampleChamberID[which(refSampleNames[i]==SampleChamberID[,2]),1]))
	}
	if (any(remRefs %in% colnames(FCvals)))
	{
		matchIndex<-match(remRefs, colnames(FCvals))
		if (any(is.na(matchIndex)))
			matchIndex<-matchIndex[-which(is.na(matchIndex))]		
		FCvals<-FCvals[,-matchIndex]
	}

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
	
	if (numPlates > 1)
	{
		# then make PCA plot to look for plate batch effect
		###### 1st leave in all genes so samples with NAs will be removed
		samPCAdat<-as.data.frame(t(log2(FCvals)))
		remSamsPCA<-prcomp(~ ., data=samPCAdat, na.action=na.omit, scale=TRUE, center=TRUE)
		numPCAsams<-sum(apply(FCvals, 2, function(x) {sum(is.na(x))})==0)
		numOrigSams<-ncol(FCorig)
		PC1var<-(remSamsPCA$sdev[1]^2)/(sum(remSamsPCA$sdev^2))*100
		PC2var<-(remSamsPCA$sdev[2]^2)/(sum(remSamsPCA$sdev^2))*100
		possColors<-rainbow(numPlates)
		curCols<-rep("", nrow(remSamsPCA$x))
		for (i in 1:nrow(remSamsPCA$x))
		{
			tempRN<-strsplit(rownames(remSamsPCA$x)[i], "_")[[1]][1]
			curCols[i]<-possColors[which(tempRN==uniPlates)]
		}
		
		pdf(file=paste(workdir,"Results/PCA_showingPlates_AllGenes.pdf",sep=""))
		plot(x=remSamsPCA$x[,1], y=remSamsPCA$x[,2], xlab=paste("PC1 (",round(PC1var,1),"%)",sep=""), ylab=paste("PC2 (",round(PC2var,1),"%)",sep=""), main=paste("Original # of Samples:", numOrigSams, "\n# Samples in PCA:", numPCAsams), pch=19, col=curCols)
		# add a legend
		legend(x=max(remSamsPCA$x[,1]), y=max(remSamsPCA$x[,2]), xjust=1, yjust=1, pch=19, col=possColors, legend=uniPlates, cex=0.8)
		dev.off()
		
		###### 2nd remove all genes that have NAs so we get all samples 
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
		gCurCols<-rep("", nrow(remGenesPCA$x))
		for (i in 1:nrow(remGenesPCA$x))
		{
			tempRN<-strsplit(rownames(remGenesPCA$x)[i], "_")[[1]][1]
			gCurCols[i]<-possColors[which(tempRN==uniPlates)]
		}
		pdf(file=paste(workdir,"Results/PCA_showingPlates_AllSamples.pdf",sep=""))
		plot(x=remGenesPCA$x[,1], y=remGenesPCA$x[,2], xlab=paste("PC1 (",round(gPC1var,1),"%)",sep=""), ylab=paste("PC2 (",round(gPC2var,1),"%)",sep=""), main=paste("Original # of Genes:", numOrigGenes, "\n# Genes in PCA:", numPCAgenes), pch=19, col=gCurCols)
		# add a legend
		legend(x=max(remGenesPCA$x[,1]), y=max(remGenesPCA$x[,2]), xjust=1, yjust=1, pch=19, col=possColors, legend=uniPlates, cex=0.8)
		dev.off()
	}
}

makePCAplotForAssays<-function(FCvals, workdir, retMap, refSamples, SampleChamberID)
{	
	# remove HK genes 
	if (any(retMap[,1] %in% rownames(FCvals)))
	{
		matchIndex<-match(retMap[,1], rownames(FCvals))
		if (any(is.na(matchIndex)))
			matchIndex<-matchIndex[-which(is.na(matchIndex))]
		FCvals<-FCvals[-matchIndex,]
	}
	# remove ref samples
	refSampleNames<-trimString(refSamples)
	remRefs<-c()
	for (i in 1:length(refSampleNames))
	{
		remRefs<-c(remRefs, as.character(SampleChamberID[which(refSampleNames[i]==SampleChamberID[,2]),1]))
	}
	if (any(remRefs %in% colnames(FCvals)))
	{
		matchIndex<-match(remRefs, colnames(FCvals))
		if (any(is.na(matchIndex)))
			matchIndex<-matchIndex[-which(is.na(matchIndex))]		
		FCvals<-FCvals[,-matchIndex]
	}

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
	
	# remove samples with failures so all assays are plotted
	remCols<-which(apply(FCvals, 2, function(x) {sum(is.na(x))}) > 0)
	FCvalsSub<-FCvals[,-remCols]
	assayPCAdat<-as.data.frame(log2(FCvalsSub))
	assayPCA<-prcomp(~ ., data=assayPCAdat, na.action=na.omit, scale=TRUE, center=TRUE)
	
	aPC1var<-(assayPCA$sdev[1]^2)/(sum(assayPCA$sdev^2))*100
	aPC2var<-(assayPCA$sdev[2]^2)/(sum(assayPCA$sdev^2))*100

	pdf(file=paste(workdir,"Results/PCA_studyAssays.pdf",sep=""))
	plot(x=assayPCA$x[,1], y=assayPCA$x[,2], pch=19, xlab=paste("PC1 (", round(aPC1var, 1), "%)", sep=""), ylab=paste("PC2 (", round(aPC2var, 1), "%)", sep=""), main=paste("PCA of ", nrow(assayPCA$x), " assays", sep=""))
	dev.off()
}
