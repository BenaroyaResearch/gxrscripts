######
# study correlation of housekeeping genes
# ctValueMatrix = the matrix of Ct values (rows are genes and columns are samples)
# HKgenes = character vector of the housekeeping gene names
######
HKgeneCorrelation<-function(ctValueMatrix, HKgenes)
{
	if (!(all(HKgenes %in% rownames(ctValueMatrix))))
	{
		matchIndex<-match(HKgenes, rownames(ctValueMatrix))
		print(paste("This housekeeping gene:", HKgenes[which(is.na(matchIndex))], "is not in the Fluidigm gene list."))
	}
	else
	{
		matchIndex<-match(HKgenes, rownames(ctValueMatrix))
		corVals<-c()
		for (i in 1:length(matchIndex))
		{
			corVals<-cbind(corVals, as.numeric(as.character(ctValueMatrix[matchIndex[i],])))
		}
		colnames(corVals)<-HKgenes
		rownames(corVals)<-colnames(ctValueMatrix)
		HKcor<-cor(corVals, use="complete.obs")
		return(HKcor)
	}
}

######
#
######
StudyHKCorrelationIncludePlates<-function(ctValueMatrix, HKgenes, maxAssay)
{
	# split into multiple plates
	setList<-list()
	for (i in 1:maxAssay)
	{
		setList[[i]]<-ctValueMatrix[grep(paste("A",i,"_",sep=""), rownames(ctValueMatrix)),]
		rownames(setList[[i]])<-unlist(lapply(strsplit(rownames(setList[[i]]), "_"), function(x) {x[2]}))
	}
	# look into M values from Genome Biology 2002 paper
	

	# then get the correlation among HK genes for each plates
	corList<-list()
	for (i in 1:maxAssay)
	{
		corList[[i]]<-HKgeneCorrelation(setList[[i]], HKgenes)
	}
	# now how to decide on the HK genes?
	corResults<-lapply(corList, corThreshold, threshold=0.8)

}

######
# count number of times a HK gene has a correlation that doesn't pass a threshold, 
# indicating it doesn't correlate well with other HK genes
# corVals = the correlation matrix obtained from the function, HKgeneCorrelation
# threshold = the correlation threshold that we would like the HK genes to have (default is 0.9)
######
corThreshold<-function(corVals, threshold = 0.9)
{
	apply(corVals, 1, function(x)
	{
		sum(x >= threshold)
	})
}

######
# make a conversion matrix between HK genes and wells on the plate
# HKgenes = the housekeeping genes
# dataList = the fluidigm data (list of matrices with 9216 rows)
######
convertHKtoWell<-function(HKgenes, dataList)
{
	returnMapHKWell<-c()
	for (i in 1:length(dataList))
	{
		curChamber<-as.character(dataList[[i]]$Chamber.ID)
		assays<-unlist(lapply(strsplit(curChamber, "-"), function(x) {x[2]}))
		curHKgenes<-as.character(dataList[[i]]$EvaGreen.Name)
		
		tempAssay<-unlist(lapply(strsplit(names(dataList)[i], "_"), function(x) {x[2]}))
		sampleAssay<-substr(tempAssay,regexpr("A",tempAssay)[1],nchar(tempAssay))	
		completeAssay<-paste(sampleAssay,"_", assays, sep="")
		mapHKWell<-unique(cbind(completeAssay, curHKgenes))
		keepIndex<-match(HKgenes, mapHKWell[,2])
		returnMapHKWell<-rbind(returnMapHKWell, mapHKWell[keepIndex,])
	}
	# get unique rows - may have repeats due to multiple sample plates
	returnMapHKWell<-unique(returnMapHKWell)
	# sort the HK genes to ensure that they are in the same order
	returnMapHKWell<-returnMapHKWell[order(returnMapHKWell[,2]),]
	return(returnMapHKWell)
}

######
# calculate M for HK genes (see 2002 Genome Biology paper: "Accurate normalization of real-time quantitate RT-PCR data 
# by geometric averaging of multiple internal control genes")
# 5/14/13 if there are 3 or less HK genes, then use all of them
######
calculateMforHK<-function(ctValueMatrix, mapHK)
{
	HKgenes<-mapHK[,1]
	numHKgenes<-length(unique(mapHK[,2]))
	# if there are greater than 3 HK genes, then pick the best ones
	if (numHKgenes > 3)
	{
		matchIndex<-c()
		remHK<-c()
		for (i in 1:numHKgenes)
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
#			curHKgenes<-paste(tempAssay, "_", HKgenes, sep="")
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
			retHK<-c(retHK, mapHK[which(actualHK[i]==mapHK[,2]),1])
		}
	}
	else
	{
		# use all input HK genes because there are 3 or less
		retHK<-HKgenes
	}
	return(retHK)
}


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

######
# calculate the geometric mean for a set of HK genes
# ctValueMatrix = the matrix of Ct values
# HKgenes = the housekeeping genes
######
calculateGeoMean<-function(ctValueMatrix, HKgenes)
{
	matchIndex<-c()
	for (i in 1:length(HKgenes))
	{
		matchIndex<-c(matchIndex, grep(HKgenes[i], rownames(ctValueMatrix)))
	}
	curAssays<-sort(unique(unlist(lapply(strsplit(rownames(ctValueMatrix)[matchIndex], "_"), function(x){x[1]}))))
	# calculate geometric mean for each plate
	geoMeans<-list()
	for (i in 1:length(curAssays))
	{
		tempAssay<-curAssays[i]
		curHKgenes<-paste(tempAssay, "_", HKgenes, sep="")
		matchIndex<-match(curHKgenes, rownames(ctValueMatrix))
		# now calculate geometric mean
		HKprod<-apply(ctValueMatrix[matchIndex,], 2, function(x)
		{
			prod(as.numeric(as.character(x)), na.rm=T)
		})
		allLen<-apply(ctValueMatrix[matchIndex,], 2, function(x){sum(!is.na(x))})
		HKgeoMean<-HKprod^(1/allLen)
		geoMeans[[i]]<-HKgeoMean
	}
	return(geoMeans)
}

######
# plots of HKgenes
# ctValueMatrix = the matrix of Ct values (rows are genes and columns are samples)
# HKgenes = character vector of the housekeeping gene names
######
HKgenePlots<-function(ctValueMatrix, HKgenes)
{
	if (!(all(HKgenes %in% rownames(ctValueMatrix))))
	{
		matchIndex<-match(HKgenes, rownames(ctValueMatrix))
		print(paste("This housekeeping gene:", HKgenes[which(is.na(matchIndex))], "is not in the Fluidigm gene list."))
	}
	else
	{
		matchIndex<-match(HKgenes, rownames(ctValueMatrix))
		corVals<-c()
		for (i in 1:length(matchIndex))
		{
			corVals<-cbind(corVals, as.numeric(as.character(ctValueMatrix[matchIndex[i],])))
		}
		colnames(corVals)<-HKgenes
		rownames(corVals)<-colnames(ctValueMatrix)
		if (length(matchIndex)==3)
			par(mfrow=c(1,3))
		if (length(matchIndex)==4)
			par(mfrow=c(2,3))
		if (length(matchIndex)==5)
			par(mfrow=c(2,5))
		if (length(matchIndex)==6)
			par(mfrow=c(3,5))
		if (length(matchIndex)==7)
			par(mfrow=c(3,7))
		if (length(matchIndex)==8)
			par(mfrow=c(4,7))
		for (i in 1:(length(matchIndex)-1))
		{
			for (j in (i+1):length(matchIndex))
			{
				plot(x=corVals[,i], y=corVals[,j], xlab=colnames(corVals)[i], ylab=colnames(corVals)[j], 
					main=round(cor(corVals[,i], corVals[,j], use="complete.obs"),3))
				curLm<-lm(corVals[,j]~corVals[,i])
				abline(a=curLm$coef[1], b=curLm$coef[2], col="red")
			}
		}
	}
}

########
# make line plots to look at the correlation among HK genes and to see the range of HK genes
########
HKgeneLinePlots<-function(ctValueMatrix, HKgenes)
{
	HKmat<-ctValueMatrix[match(HKgenes, rownames(ctValueMatrix)),]
	allCols<-c("black","red","orange","yellow","green","blue","purple","gray")
#	orderIndex<-order(as.numeric(HKmat[1,]))
	orderIndex<-1:ncol(HKmat)
	plot(x=1:ncol(HKmat), y=as.numeric(HKmat[1,])[orderIndex], ylim=c(min(as.numeric(HKmat), na.rm=T), max(as.numeric(HKmat), na.rm=T)), 
		xlab="Sample", ylab="Ct", col=allCols[1], type="l")
	for (i in 2:nrow(HKmat))
	{
		lines(x=1:ncol(HKmat), y=as.numeric(HKmat[i,])[orderIndex], col=allCols[i])
	}
	# add legend
	legend(x=(ncol(HKmat)/2), y=max(as.numeric(HKmat), na.rm=T), lty=1, col=allCols[1:nrow(HKmat)], legend=rownames(HKmat), cex=0.8)
}

########
# look at the range of HK genes compared to non-HK genes
########
HKvsNonHKgenesBoxplot<-function(ctValueMatrix, HKgenes)
{
	HKmat<-ctValueMatrix[match(HKgenes, rownames(ctValueMatrix)),]
	nonHKmat<-ctValueMatrix[which(!((1:nrow(ctValueMatrix)) %in% match(HKgenes, rownames(ctValueMatrix)))),]
	boxplot(c(as.numeric(HKmat), as.numeric(nonHKmat))~c(rep(1,length(as.numeric(HKmat))),rep(2,length(as.numeric(nonHKmat)))), axes=FALSE,
		ylab="Ct")
	box()
	axis(1, at=1:2, labels=c("HK genes", "non HK genes"))
	axis(2)
}

########
# plot of CV (coefficient of variation)
########
plotCV<-function(ctValueMatrix, HKgenes)
{
	HKmat<-ctValueMatrix[match(HKgenes, rownames(ctValueMatrix)),]
	CVvals<-c()
	meanVals<-c()
	for (i in 1:nrow(HKmat))
	{
		curVals<-as.numeric(HKmat[i,])
		curMean<-mean(curVals, na.rm=T)
		curSD<-sd(curVals, na.rm=T)
		meanVals<-c(meanVals, curMean)
		CVvals<-c(CVvals, curSD/curMean)
	}
	plot(x=1:length(HKgenes), y=CVvals, ylab="CV", axes=FALSE, xlab="", pch=19)
	box()
	axis(1, at=1:length(HKgenes), labels=HKgenes)
	axis(2)
}

########
# if there is more than one plate in the assay, determine the HK genes by looking for overlap genes between the 2 plates
# ctValueMatrix1 = the matrix of Ct values for plate 1 (rows are genes and columns are samples)
# ctValueMatrix2 = the matrix of Ct values for plate 2 (rows are genes and columns are samples)
########
determineHKgenes<-function(ctValueMatrix1, ctValueMatrix2)
{
	HKgenes<-rownames(ctValueMatrix1)[which(rownames(ctValueMatrix1) %in% rownames(ctValueMatrix2))]
	return(HKgenes)
}

########
# remove questionable samples (based on being outliers for comparing HK genes across plates)
# ctValueMatrix = the matrix of Ct values
# mapHK = the matrix of wells and housekeeping genes
# maxAssay = the number of plates for the assay (e.g. genomic CBC has 2 plates)
########
remQuestSam<-function(ctValueMatrix, mapHK, maxAssay, workdir)
{
	setList<-list()
	listNames<-c()
	for (i in 1:maxAssay)
	{
		setList[[i]]<-ctValueMatrix[grep(paste("A",i,"_",sep=""), rownames(ctValueMatrix)),]
		listNames<-c(listNames, paste("A",i,sep=""))
		#rownames(setList[[i]])<-unlist(lapply(strsplit(rownames(setList[[i]]), "_"), function(x) {x[2]}))
	}
	# now loop through all 2-way combinations
	questSam<-c()
	for (i in 1:(length(setList)-1))
	{
		for (j in (i+1):length(setList))
		{	
			corVals<-corHKacrossPlates(setList[[i]], setList[[j]], mapHK, listNames[i], listNames[j], workdir)
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
					questSam<-c(questSam, plotHKGeneAcrossPlate(setList[[i]], setList[[j]], lowCor[k], index1, index2, workdir, listNames[i], listNames[j]))
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
# aname1 = name of assay 1 (to print plot file)
# aname2 = name of assay 2
# workdir = directory to print out plot
########
corHKacrossPlates<-function(ctValueMatrix1, ctValueMatrix2, mapHK, aname1, aname2, workdir)
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
	pdf(file=paste(workdir,"Results/HKgenesAcross_",aname1,"_",aname2,".pdf",sep=""))
	if (min(corVals) > 0.5)
	{
		plot(y=corVals, x=1:length(HKgenes), pch=19, col=allCols[1:length(HKgenes)], ylim=c(0.5,1), axes=FALSE, ylab="Correlation", xlab="",
			main="Correlation Across Plates")
		abline(h=0.9, lty=2)
		box()
		axis(1, at=1:length(HKgenes), labels=HKgenes, las=2, cex.axis=0.7)
		axis(2)
		legend(x=1,y=0.7, legend=HKgenes, pch=19, col=allCols[1:length(HKgenes)], cex=0.8)
	}	
	else
	{
		plot(y=corVals, x=1:length(HKgenes), pch=19, col=allCols[1:length(HKgenes)], ylim=c(min(corVals),1), axes=FALSE, ylab="Correlation", 
			xlab="", main="Correlation Across Plates")
		abline(h=0.9, lty=2)	
		box()
		axis(1, at=1:length(HKgenes), labels=HKgenes, las=2, cex.axis=0.7)
		axis(2)
		legend(x=1, y=0.7, legend=HKgenes, pch=19, col=allCols[1:length(HKgenes)], cex=0.8)
	}
	dev.off()
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
# aname1 = name of assay 1 (to print plot file)
# aname2 = name of assay 2
# workdir = directory to print out plot
#####
plotHKGeneAcrossPlate<-function(ctValueMat1, ctValueMat2, curHKgene, index1, index2, workdir, aname1, aname2)
{
	# remove any backslashes in HK gene name because it messed up the file name
	curHKgene<-gsub("/","_",curHKgene)
	pdf(file=paste(workdir,"Results/LowCorrelationAcross_",aname1,"_",aname2, "_ForHKGene_", curHKgene,".pdf",sep=""))
	plot(x=ctValueMat1[index1,], y=ctValueMat2[index2,], xlab="Plate 1", ylab="Plate 2", main=curHKgene)
	curLm<-lm(as.numeric(ctValueMat2[index2,]) ~ as.numeric(ctValueMat1[index1,]))
	abline(a=curLm$coef[1], b=curLm$coef[2], col="red")
	dev.off()
	# potential outliers
	return(colnames(ctValueMat1)[as.numeric(names(which(abs(curLm$residuals)>4)))])
}
