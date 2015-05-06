######
# calculate M for HK genes (see 2002 Genome Biology paper: "Accurate normalization of real-time quantitate RT-PCR data 
# by geometric averaging of multiple internal control genes")
# ctValueMatrix = the matrix of Ct values
# mapHK = the map between sample ID and sample name for HK genes
######
calculateMforHK<-function(ctValueMatrix, mapHK)
{
	HKgenes<-mapHK[,1]
	matchIndex<-c()
	remHK<-c()
	for (i in 1:length(HKgenes))
	{
		if (length(grep(HKgenes[i], rownames(ctValueMatrix)))>0)				
			matchIndex<-c(matchIndex, grep(HKgenes[i], rownames(ctValueMatrix)))
		else
		{
			print(paste("This housekeeping gene:", HKgenes[i], "is not in the Fluidigm gene list."))
			remHK<-c(remHK, i)
		}
	}
	if (length(remHK)>0)
	{
		HKgenes<-HKgenes[-i]
		matchIndex<-c()
		for (i in 1:length(HKgenes))
		{
			matchIndex<-c(matchIndex, grep(HKgenes[i], rownames(ctValueMatrix)))
		}
	}

	curAssays<-sort(unique(unlist(lapply(strsplit(rownames(ctValueMatrix)[matchIndex], "_"), function(x){x[1]}))))
	# calculate M values for each plate
	Mvalues<-list()
	for (i in 1:length(curAssays))
	{
		tempAssay<-curAssays[i]
		curHKgenes<-HKgenes[grep(paste(tempAssay, "_", sep=""), HKgenes)]
		
		# make a list of all possible combinations between 3 and number of HK genes
		comboList<-list()
		for (z in 3:(length(curHKgenes)-1))
		{
			curCombMatrix<-combn(length(curHKgenes), z)
			for (y in 1:ncol(curCombMatrix))
			{
				curLen<-length(comboList)
				comboList[[curLen+1]]<-curCombMatrix[,y]
			}
		}
		curLen<-length(comboList)
		comboList[[curLen+1]]<-1:length(curHKgenes)

		assayMvalues<-c()
		for (j in 1:length(comboList))
		{
			retValue<-calculateMForOneCombo(comboList[[j]], ctValueMatrix, curHKgenes)
			assayMvalues<-c(assayMvalues, retValue)
		}
		Mvalues[[i]]<-assayMvalues			
	}

	# now sum the values from multiple plates and see which combination of HK gives the lowest average SD
	sumMvalues<-c()
	for (i in 1:length(Mvalues))
	{
		if (i==1)
			sumMvalues<-Mvalues[[i]]
		else
			sumMvalues<-sumMvalues + Mvalues[[i]]
	}
	# now find genes with smallest M value
	curIndex<-comboList[[which(sumMvalues==min(sumMvalues))[1]]]
	# need to use a subset of the HK genes to get the right values (because only want one plate's worth of HK genes)
	keepHKgenes<-curHKgenes[curIndex]
	# now find all matching genes on other plates
	actualHK<-mapHK[match(keepHKgenes, mapHK[,1]),2]
	# now find all matches in the mapping
	retHK<-c()
	for (i in 1:length(actualHK))
	{
		retHK<-c(retHK, as.character(mapHK[which(actualHK[i]==mapHK[,2]),1]))
	}
	return(retHK)
}

######
# calculate the M value for one possible combination of HK genes (see Genome Biology 2002)
# oneCombo = one possible combination of HK genes (must have at least 3 HK genes)
# ctValueMatrix = the matrix of Ct values 
# curHKgenes = the current HK genes on one plate
######
calculateMForOneCombo<-function(oneCombo, ctValueMatrix, curHKgenes)
{
	matchIndex<-match(curHKgenes, rownames(ctValueMatrix))
	curLen<-length(oneCombo)
	# find all 2-way ratios for this combo
	curComb<-combn(curLen, 2)
	# get ratio for each 2-way option
	allSD<-c()
	for (i in 1:ncol(curComb))
	{
		index1<-matchIndex[oneCombo[curComb[,i]]][1]
		index2<-matchIndex[oneCombo[curComb[,i]]][2]
		curSD<-sd(log2(as.numeric(as.character(ctValueMatrix[index1,]))/as.numeric(as.character(ctValueMatrix[index2,]))), na.rm=T)
		allSD<-c(allSD, curSD)
	}
	curMvalue<-mean(allSD)
	return(curMvalue)
}

########
# remove questionable samples (based on being outliers for comparing HK genes across plates)
# ctValueMatrix = the matrix of Ct values
# mapHK = the matrix of wells and housekeeping genes
# maxAssay = the number of plates for the assay (e.g. genomic CBC has 2 plates)
########
remQuestSam<-function(ctValueMatrix, mapHK, maxAssay)
{
	setList<-list()
	for (i in 1:maxAssay)
	{
		setList[[i]]<-ctValueMatrix[grep(paste("A",i,"_",sep=""), rownames(ctValueMatrix)),]
	}
	# now loop through all 2-way combinations
	questSam<-c()
	for (i in 1:(length(setList)-1))
	{
		for (j in (i+1):length(setList))
		{
			corVals<-corHKacrossPlates(setList[[i]], setList[[j]], mapHK)
			if (any(corVals < 0.9))
			{
				lowCor<-names(corVals)[which(corVals < 0.9)]
				for (k in 1:length(lowCor))
				{
					curHK<-lowCor[k]
					curAssays<-mapHK[which(curHK==mapHK[,2]),1]
					index1<-match(curAssays, rownames(setList[[i]]))
					index1<-index1[-which(is.na(index1))]
					index2<-match(curAssays, rownames(setList[[j]]))
					index2<-index2[-which(is.na(index2))]
					questSam<-c(questSam, plotHKGeneAcrossPlate(setList[[i]], setList[[j]], lowCor[k], index1, index2))
				}
				questSam<-unique(questSam)
			}
		}
	}
	if (length(questSam)>0)
		ctValueMatrix<-ctValueMatrix[,-match(questSam, colnames(ctValueMatrix))]
	return(ctValueMatrix)
}

########
# look at the correlation of HK genes across plates (requires more than one plate for the assay)
# ctValueMatrix1 = the matrix of Ct values for plate 1 (rows are genes and columns are samples)
# ctValueMatrix2 = the matrix of Ct values for plate 2 (rows are genes and columns are samples)
# mapHK = the matrix of wells and housekeeping genes
########
corHKacrossPlates<-function(ctValueMatrix1, ctValueMatrix2, mapHK)
{
	corVals<-c()
	if (ncol(ctValueMatrix1) != ncol(ctValueMatrix2))
	{
		keepCols<-colnames(ctValueMatrix1)[which(colnames(ctValueMatrix1) %in% colnames(ctValueMatrix2))]
		ctValueMatrix1<-ctValueMatrix1[,match(keepCols, colnames(ctValueMatrix1))]
		ctValueMatrix2<-ctValueMatrix2[,match(keepCols, colnames(ctValueMatrix2))]
	}
	if (!all(colnames(ctValueMatrix1)==colnames(ctValueMatrix2)))
	{
		keepCols<-colnames(ctValueMatrix1)[which(colnames(ctValueMatrix1) %in% colnames(ctValueMatrix2))]
		ctValueMatrix1<-ctValueMatrix1[,match(keepCols, ctValueMatrix1)]
		ctValueMatrix2<-ctValueMatrix2[,match(keepCols, ctValueMatrix2)]
	}
	HKgenes<-unique(mapHK[,2])	
	for (i in 1:length(HKgenes))
	{
		curHK<-HKgenes[i]
		curAssays<-mapHK[which(curHK==mapHK[,2]),1]
		curIndex1<-match(curAssays, rownames(ctValueMatrix1))
		curIndex1<-curIndex1[-which(is.na(curIndex1))]
		curIndex2<-match(curAssays, rownames(ctValueMatrix2))
		curIndex2<-curIndex2[-which(is.na(curIndex2))]
		curCor<-cor(as.numeric(as.character(ctValueMatrix1[curIndex1,])), 
			as.numeric(as.character(ctValueMatrix2[curIndex2,])), use="complete.obs")
		corVals<-c(corVals, curCor)
	}
	allCols<-c("black","red","orange","yellow","green","blue","purple","gray")
	if (min(corVals) > 0.5)
	{
		plot(y=corVals, x=1:length(HKgenes), pch=19, col=allCols[1:length(HKgenes)], ylim=c(0.5,1), axes=FALSE, ylab="Correlation", xlab="",
			main="Correlation Across Plates")
		box()
		axis(1, at=1:length(HKgenes), labels=HKgenes, las=2, cex.axis=0.7)
		axis(2)
		legend(x=1,y=0.7, legend=HKgenes, pch=19, col=allCols[1:length(HKgenes)])
	}	
	else
	{
		plot(y=corVals, x=1:length(HKgenes), pch=19, col=allCols[1:length(HKgenes)], ylim=c(min(corVals),1), axes=FALSE, ylab="Correlation", 
			xlab="", main="Correlation Across Plates")
		box()
		axis(1, at=1:length(HKgenes), labels=HKgenes, las=2, cex.axis=0.7)
		axis(2)
		legend(x=1, y=0.7, legend=HKgenes, pch=19, col=allCols[1:length(HKgenes)])
	}
	if (min(corVals) < 0.9)
		print("WARNING: The minimum correlation is lower than expected, indicating that there may be a plate batch effect.")
	names(corVals)<-HKgenes
	return(corVals)
}

#####
# look at correlation of HK gene across plates in a plot to see if any samples are outliers
# ctValueMat1 = the Ct value matrix for plate 1
# ctValueMat2 = the Ct value matrix for plate 2
# curHKgene = the current HK gene to plot
# index1 = the row index for that particular HK gene on plate 1
# index2 = the row index for that particular HK gene on plate 2
#####
plotHKGeneAcrossPlate<-function(ctValueMat1, ctValueMat2, curHKgene, index1, index2)
{
	plot(x=ctValueMat1[index1,], y=ctValueMat2[index2,], xlab="Plate 1", ylab="Plate 2", main=curHKgene)
	curLm<-lm(as.numeric(ctValueMat2[index2,]) ~ as.numeric(ctValueMat1[index1,]))
	abline(a=curLm$coef[1], b=curLm$coef[2], col="red")
	# potential outliers
	return(colnames(ctValueMat1)[as.numeric(names(which(abs(curLm$residuals)>4)))])
}
