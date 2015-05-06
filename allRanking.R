######
# Function to perform all 2-way comparisons and then combine to overall ranking
# Steps:
# 1. perform preprocessing (2 ways: 1 with all data and 1 with PALO)
# 2. call ranking function for each 2-way comparison (needs to be called twice for each comparison: once
#    for all data and once for PALO) - include plotting option
# 3. combine into overall ranking for full data set
# Input Parameters:
# 1. curData - the data frame that contains the data (should include both signal intensity and detection p-values)
# 2. designData - the design data frame that matches samples to groups
# 2. caseNames - a vector of all the case names (must be the same length as controlNames)
# 3. controlNames - a vector of all the control names (must be the same length as caseNames)
#    will loop through the vector of caseNames and controlNames for all the 2-way comparisons so they need to
#    match up in terms of what 2-way comparisons to perform (for example, a controlName could appear twice)
# 4. toPlot - a boolean indicating whether plots of Diff vs. FC should be created
# 5. fileDirectory - the directory where the plot files should be placed
######
performAllRanking<-function(curData, designData, caseNames, controlNames, toPlot=FALSE, fileDirectory)
{
	if (length(caseNames)!=length(controlNames))
		print("ERROR: Must have the same number of cases and controls.")

	colnames(curData)[1]<-"PROBE_ID"

	curData[,2:ncol(curData)]<-apply(curData[,2:ncol(curData)], 2, function(x)
	{
		as.numeric(as.character(x))
	})
	
	# check that data are on the original intensity scale
	expSignal<-curData[,grep("Signal", colnames(curData))]
	if (max(expSignal, na.rm=TRUE) < 20)
	{
		# then assume the data has already been log2 transformed
		# want to put back on original intensity scale
		expSignal<-2^expSignal
		curData[,grep("Signal", colnames(curData))]<-expSignal
	}
	
	# perform preprocessing without and with PALO
	curDataPre<-preprocess(expDataFile=curData, outputFile=NA, sigColName="Signal", 		pvalCol="Detection", allCols=FALSE, performNormalize=TRUE, performFloor=TRUE, 		performLog2=FALSE, performPALO=FALSE)
	curDataPrePalo<-preprocess(expDataFile=curData, outputFile=NA, sigColName="Signal", 		pvalCol="Detection", allCols=FALSE, performNormalize=TRUE, performFloor=TRUE, 		performLog2=FALSE, performPALO=TRUE)

	# need to call the perform one ranking for each case vs. control comparison
	allList<-list()
	PALOList<-list()
	comparisons<-c()
	for (i in 1:length(caseNames))
	{
		retListAll<-performOneRanking(curDataPre, designData, caseNames[i], controlNames[i], toPlot, fileDirectory, type="All")
		retListPALO<-performOneRanking(curDataPrePalo, designData, caseNames[i], controlNames[i], toPlot, fileDirectory, type="PALO")
		allList[[i]]<-retListAll
		PALOList[[i]]<-retListPALO
		comparisons<-c(comparisons, paste(caseNames[i], "vs", controlNames[i], sep="_"))
	}
	comparisons<-gsub(" ", "_", comparisons)
	
	# if there is more than one comparison, then need to combine into one ranking for All and one ranking for PALO
	# but perform the same steps regardless of the number of comparisons
	# get rankings based on all comparisons
	aBothOrder<-c()
	for (i in 1:length(allList))
	{
		curBothOrder<-as.numeric(as.character(allList[[i]][[1]]$BothOrder))
		aBothOrder<-cbind(aBothOrder, curBothOrder)
	}
	# now take average of the order
	rankBothAll<-order(apply(aBothOrder, 1, mean, na.rm=TRUE), decreasing=FALSE)
	allBothOrder<-match(1:nrow(curDataPre), rankBothAll)
		
	pBothOrder<-c()
	for (i in 1:length(PALOList))
	{
		curBothOrder<-as.numeric(as.character(PALOList[[i]][[1]]$BothOrder))
		pBothOrder<-cbind(pBothOrder, curBothOrder)
	}
	# now take average of the order
	rankBothPALO<-order(apply(pBothOrder, 1, mean, na.rm=TRUE), decreasing=FALSE)
	paloBothOrder<-match(1:nrow(curDataPrePalo), rankBothPALO)

	# now create the returnList
	keepTempOutAll<-c()
	for (i in 1:length(allList))
	{
		curTempOut<-allList[[i]][[1]]
		colnames(curTempOut)[2:ncol(curTempOut)]<-paste(comparisons[i], colnames(curTempOut)[2:ncol(curTempOut)], sep="_")
		if (i == 1)
			keepTempOutAll<-curTempOut
		else
			keepTempOutAll<-cbind(keepTempOutAll, curTempOut[,2:ncol(curTempOut)])
		# write out curTempOut to file - now write out 3 files for each
		curTempOutFC<-curTempOut[,c(1,grep("FC", colnames(curTempOut)), grep("Mean",colnames(curTempOut)), grep("Pval",colnames(curTempOut)))]
		curTempOutDiff<-curTempOut[,c(1,grep("Diff", colnames(curTempOut)), grep("Mean",colnames(curTempOut)), grep("Pval",colnames(curTempOut)))]
		curTempOutBoth<-curTempOut[,c(1,grep("FC$", colnames(curTempOut)), grep("Diff$", colnames(curTempOut)),grep("Both", colnames(curTempOut)),
				grep("Mean",colnames(curTempOut)), ncol(curTempOut))]

		csvFileNameFC<-paste(fileDirectory, "/", "fc_All_", caseNames[i], "_vs_", controlNames[i], ".csv",sep="")
		csvFileNameFC<-gsub(" ", "_", csvFileNameFC)
		csvFileNameDiff<-paste(fileDirectory, "/", "diff_All_", caseNames[i], "_vs_", controlNames[i], ".csv",sep="")
		csvFileNameDiff<-gsub(" ", "_", csvFileNameDiff)
		csvFileNameBoth<-paste(fileDirectory, "/", "diffFc_All_", caseNames[i], "_vs_", controlNames[i], ".csv",sep="")
		csvFileNameBoth<-gsub(" ", "_", csvFileNameBoth)

		write.csv(curTempOutFC, file=csvFileNameFC, row.names=FALSE)
		write.csv(curTempOutDiff, file=csvFileNameDiff, row.names=FALSE)
		write.csv(curTempOutBoth, file=csvFileNameBoth, row.names=FALSE)
	}
	keepTempOutAll<-cbind(keepTempOutAll, allBothOrder)
	colnames(keepTempOutAll)[ncol(keepTempOutAll)]<-"AllBothOrder"
		
	keepTempOutPALO<-c()
	for (i in 1:length(PALOList))
	{
		curTempOut<-PALOList[[i]][[1]]
		colnames(curTempOut)[2:ncol(curTempOut)]<-paste(comparisons[i], colnames(curTempOut)[2:ncol(curTempOut)], sep="_")
		if (i == 1)
			keepTempOutPALO<-curTempOut
		else
			keepTempOutPALO<-cbind(keepTempOutPALO, curTempOut[,2:ncol(curTempOut)])
		# write out curTempOut to file - now write out 3 files for each
		curTempOutFC<-curTempOut[,c(1,grep("FC", colnames(curTempOut)), grep("Mean",colnames(curTempOut)), grep("Pval",colnames(curTempOut)))]
		curTempOutDiff<-curTempOut[,c(1,grep("Diff", colnames(curTempOut)), grep("Mean",colnames(curTempOut)), grep("Pval",colnames(curTempOut)))]
		curTempOutBoth<-curTempOut[,c(1,grep("FC$", colnames(curTempOut)), grep("Diff$", colnames(curTempOut)),grep("Both", colnames(curTempOut)),
				grep("Mean",colnames(curTempOut)), ncol(curTempOut))]

		csvFileNameFC<-paste(fileDirectory, "/", "fc_PALO_", caseNames[i], "_vs_", controlNames[i], ".csv",sep="")
		csvFileNameFC<-gsub(" ", "_", csvFileNameFC)
		csvFileNameDiff<-paste(fileDirectory, "/", "diff_PALO_", caseNames[i], "_vs_", controlNames[i], ".csv",sep="")
		csvFileNameDiff<-gsub(" ", "_", csvFileNameDiff)
		csvFileNameBoth<-paste(fileDirectory, "/", "diffFc_PALO_", caseNames[i], "_vs_", controlNames[i], ".csv",sep="")
		csvFileNameBoth<-gsub(" ", "_", csvFileNameBoth)

		write.csv(curTempOutFC, file=csvFileNameFC, row.names=FALSE)
		write.csv(curTempOutDiff, file=csvFileNameDiff, row.names=FALSE)
		write.csv(curTempOutBoth, file=csvFileNameBoth, row.names=FALSE)
	}
	keepTempOutPALO<-cbind(keepTempOutPALO, paloBothOrder)
	colnames(keepTempOutPALO)[ncol(keepTempOutPALO)]<-"AllBothOrder"
		
	# now make the outBothAll and outBothPALO data frames
	colKeepAll<-sort(c(1,grep("FC$", colnames(keepTempOutAll)), grep("Diff$", colnames(keepTempOutAll))))
	colKeepPALO<-sort(c(1,grep("FC$", colnames(keepTempOutPALO)), grep("Diff$", colnames(keepTempOutPALO))))
		
	# note that the ordering can be slightly different between these two at the beginning if some of the probes 
	# that don't pass PALO end up slightly changing the ordering
	outBothAll<-keepTempOutAll[rankBothAll, colKeepAll]
	outBothPALO<-keepTempOutPALO[rankBothPALO, colKeepPALO]

	# write out a file for all ordering - now just include PROBE_ID and Order (rather than all FC and Diff)
	colKeepAll<-sort(c(1, grep("AllBothOrder", colnames(keepTempOutAll))))
	colKeepPALO<-sort(c(1,grep("AllBothOrder", colnames(keepTempOutPALO))))
	
	writeAll<-keepTempOutAll[,colKeepAll]
	writePALO<-keepTempOutPALO[,colKeepPALO]

	csvFileName<-paste(fileDirectory, "/", "All", ".csv",sep="")
	csvFileName<-gsub(" ", "_", csvFileName)
	write.csv(writeAll, file=csvFileName, row.names=FALSE)	
	csvFileName<-paste(fileDirectory, "/", "PALO", ".csv",sep="")
	csvFileName<-gsub(" ", "_", csvFileName)
	write.csv(writePALO, file=csvFileName, row.names=FALSE)	
	retList<-list(dataOutAll=keepTempOutAll, outBothAll=outBothAll, 
				dataOutPALO=keepTempOutPALO, outBothPALO=outBothPALO)

	return(retList)
}

######
# perform ranking for one case vs. control comparison
# Input Parameters:
# 1. tempData - the data set that has been preprocessed (all or PALO data)
# 2. designData - the design data set
# 3. caseName - the name of the case group
# 4. controlName - the name of the control group
# 5. toPlot - a boolean indicating whether plot of Diff vs. FC should be created
# 6. fileDirectory - the directory to place the plot
# 7. type - "All" or "PALO" to indicate whether this is all data or PALO passing data (to include in plot file name)
######
performOneRanking<-function(tempData, designData, caseName, controlName, toPlot, fileDirectory, type="All")
{
	# figure out case and control columns
	caseBarcodes<-as.character(designData$barcode[which(as.character(designData$group_label)==as.character(caseName))])
	controlBarcodes<-as.character(designData$barcode[which(as.character(designData$group_label)==as.character(controlName))])
	caseIndex<-c()
	for (i in 1:length(caseBarcodes))
		caseIndex<-c(caseIndex, grep(caseBarcodes[i], colnames(tempData)))
	controlIndex<-c()
	for (i in 1:length(controlBarcodes))
		controlIndex<-c(controlIndex, grep(controlBarcodes[i], colnames(tempData)))

	controls<-tempData[,controlIndex]
	cases<-tempData[,caseIndex]

	# make sure that the results of subsetting are matrices
	if (length(controlIndex)==1)
		controls<-matrix(controls, ncol=1)
	if (length(caseIndex)==1)
		cases<-matrix(cases, ncol=1)

	caseMean<-apply(cases, 1, function(x)
	{
		tempX<-as.numeric(as.character(x))
		mean(tempX, na.rm=TRUE)
	})
	controlMean<-apply(controls, 1, function(x)
	{
		tempX<-as.numeric(as.character(x))
		mean(tempX, na.rm=TRUE)
	})
	
	# if both groups have sample size of 5 or greater then perform t.test
	if (ncol(cases)>=5 & ncol(controls)>5)
	{
		# use limma
		library(limma)
		all<-cbind(cases, controls)
		row.names(all)<-tempData$PROBE_ID
		designMat<-cbind(rep(1,ncol(all)), c(rep(1,ncol(cases)),rep(0,ncol(controls))))
		colnames(designMat)<-c("Control","CasevsControl")
		curFit<-lmFit(all, designMat)
		curFit<-eBayes(curFit)
		curTT<-topTable(curFit, coef="CasevsControl", adjust="BH", number=nrow(all))
		allPvalsAdjust<-curTT$adj.P.Val[match(row.names(all), curTT$ID)]

#		allPvals<-apply(all, 1, function(x)
#		{
#			if (all(x==min(x)))
#				return(1)
#			else
#				return(t.test(x[1:ncol(controls)], x[(ncol(controls)+1):ncol(all)])$p.value)
#		})
#		# now adjust for multiple testing
#		allPvalsAdjust<-p.adjust(allPvals, method="fdr")
	}
	else
		allPvalsAdjust<-rep(NA, nrow(controls))

	# now calculate FC and diff
	tempDataFC<-caseMean/controlMean
	tempDataDiff<-caseMean-controlMean

	# rank based on FC, then diff (use diff to distinguish identical FC)
	# rank based on diff, then FC (use FC to distinguish identical diff)
	rankFC<-order(abs(log2(as.numeric(tempDataFC))), abs(as.numeric(tempDataDiff)), decreasing=TRUE)
	rankDiff<-order(abs(as.numeric(tempDataDiff)), abs(log2(as.numeric(tempDataFC))), decreasing=TRUE)
	
	# correlation between FC and diff is low
	allCor<-cor(abs(log2(as.numeric(tempDataFC))), abs(as.numeric(tempDataDiff)))

	# rank based on both FC and diff
	FCorder<-match(1:nrow(tempData), rankFC)
	DiffOrder<-match(1:nrow(tempData), rankDiff)
	# take average of FC and diff rankings
	rankBoth<-order((FCorder+DiffOrder)/2, decreasing=FALSE)
	BothOrder<-match(1:nrow(tempData), rankBoth)

	# create 4 outputs for all data
	# 1. all info
	# 2. ranking by FC
	# 3. ranking by diff
	# 4. ranking by both
	# output 1
	tempDataOut<-cbind(as.character(tempData$PROBE_ID), tempDataFC, tempDataDiff, FCorder, DiffOrder, 		BothOrder, caseMean, controlMean, allPvalsAdjust)
	tempDataOut<-as.data.frame(tempDataOut)
	colnames(tempDataOut)<-c("PROBE_ID","FC", "Diff", "FCOrder", "DiffOrder", "BothOrder", 		"caseMean", "controlMean", "AdjustedPvals")
	# output 2
	outFC<-cbind(as.character(tempData$PROBE_ID)[rankFC], tempDataFC[rankFC])
	outFC<-as.data.frame(outFC)
	colnames(outFC)<-c("PROBE_ID", "FC")
	# output 3
	outDiff<-cbind(as.character(tempData$PROBE_ID)[rankDiff], tempDataDiff[rankDiff])
	outDiff<-as.data.frame(outDiff)
	colnames(outDiff)<-c("PROBE_ID", "Diff")
	# output 4
	outBoth<-cbind(as.character(tempData$PROBE_ID)[rankBoth], tempDataFC[rankBoth], 
		tempDataDiff[rankBoth])
	outBoth<-as.data.frame(outBoth)
	colnames(outBoth)<-c("PROBE_ID", "FC", "Diff")

	if (toPlot==TRUE)
	{
		plotFileName<-paste(fileDirectory, "/", caseName, "_vs_", controlName, "_", type, ".png",sep="")
		plotFileName<-gsub(" ", "_", plotFileName)
		png(plotFileName, width=600, height=400)
		# make plot of 3 different ranking results
		par(mfrow=c(1,3))
		curCol<-rep("black", length(rankFC))
		curCol[rankFC[1:10]]<-"red"
		curCol[rankFC[11:20]]<-"orange"
		curCol[rankFC[21:30]]<-"yellow"
		curCol[rankFC[31:40]]<-"green"
		curCol[rankFC[41:50]]<-"blue"
		curPch<-rep(1, length(rankFC))
		curPch[rankFC[1:50]]<-19
		plot(x=log2(as.numeric(as.character(tempDataOut[,2]))), y=as.numeric(as.character(tempDataOut[,3])), 
			xlab="log2(FC)", ylab="Diff", main=paste(type, "ranked by FC\nCor =", round(allCor, 3)), 
			col=curCol, pch=curPch)
		abline(h=0)

		# add a legend
		legend(x=min(log2(as.numeric(as.character(tempDataOut[,2]))), na.rm=TRUE), 
			y=max(as.numeric(as.character(tempDataOut[,3])), na.rm=TRUE), 
			legend=c("top 10", "top 11-20", "top 21-30", "top 31-40", "top 41-50"),
			pch=19, col=c("red","orange","yellow", "green","blue"))

		curCol<-rep("black", length(rankDiff))
		curCol[rankDiff[1:10]]<-"red"
		curCol[rankDiff[11:20]]<-"orange"
		curCol[rankDiff[21:30]]<-"yellow"
		curCol[rankDiff[31:40]]<-"green"
		curCol[rankDiff[41:50]]<-"blue"
		curPch<-rep(1, length(rankDiff))
		curPch[rankDiff[1:50]]<-19
		plot(x=log2(as.numeric(as.character(tempDataOut[,2]))), y=as.numeric(as.character(tempDataOut[,3])), 
			xlab="log2(FC)", ylab="Diff", main=paste(type, "ranked by Diff"), col=curCol, pch=curPch)
		abline(h=0)

		curCol<-rep("black", length(rankBoth))
		curCol[rankBoth[1:10]]<-"red"
		curCol[rankBoth[11:20]]<-"orange"
		curCol[rankBoth[21:30]]<-"yellow"
		curCol[rankBoth[31:40]]<-"green"
		curCol[rankBoth[41:50]]<-"blue"
		curPch<-rep(1, length(rankBoth))
		curPch[rankBoth[1:50]]<-19
		plot(x=log2(as.numeric(as.character(tempDataOut[,2]))), y=as.numeric(as.character(tempDataOut[,3])), 
			xlab="log2(FC)", ylab="Diff", main=paste(type, "ranked by Both:\nFC & Diff"), col=curCol, pch=curPch)
		abline(h=0)

		dev.off()		
	}
	return(list(tempDataOut=tempDataOut, outFC=outFC, outDiff=outDiff, outBoth=outBoth))
}
