# R code to start with Fluidigm Ct data and return delta delta Ct values
path<-"/Users/ewhalen/Projects/Fluidigm/AutomateTests/"
directory<-"P14_U01Mexico_3Plates_Phase2"

workdir<-paste(path, directory, "/", sep="")

library(combinat)

# read in scripts
#source(file=paste(workdir,"scripts/","Fluidigm_R_Functions.R",sep=""))
#source(file=paste(workdir,"scripts/","Fluidigm_R_HKgenes_Functions.R",sep=""))
# use these temporarily while I'm working on the files
source(file="/Users/ewhalen/Rscripts/Fluidigm_R_Functions.R")
source(file="/Users/ewhalen/Rscripts/Fluidigm_R_HKgenes_Functions.R")

# read in parameters
parameters<-read.csv(file=paste(workdir,"Parameters/parameters.csv",sep=""))

# read in all data files stored in Data directory - do not know how many there will be 
fileNames<-list.files(path=paste(workdir,"Data/",sep=""))
fileNames<-reorderFileNames(fileNames)

# read in data
dataList<-list()
for (i in 1:length(fileNames))
{
	dataList[[i]]<-readFluidData(paste(workdir,"Data/",fileNames[i],sep=""))
}
names(dataList)<-fileNames

# now check chamber IDs to make sure samples and genes match across plates
# if they do not match, then stop processing - need to fix the input data
anyError<-checkSampleGeneChamberIDs(dataList)

if (anyError==FALSE)
{
	# convert into matrix of rows (assays/genes) and columns (samples represented by chamber IDs)
	ctValList<-list()
	ctCallList<-list()
	for (i in 1:length(dataList))
	{
		ctValList[[i]]<-changeFluidFormat(dataList[[i]], colName="Ct.Value")
		ctCallList[[i]]<-changeFluidFormat(dataList[[i]], colName="Ct.Call")
	}

	# do PALX across all genes and samples
	# note the row names will not be unique (will have overlap from HK genes)
	retList<-combineDataFromMultiplePlates(dataList, ctValList, ctCallList)

	curPercentage<-returnParameter(parameters, "PALX", "numeric")

	ctValPalx<-performPALXBoth(retList$ctCallDataForPalx, retList$ctValDataForPalx, percentage=curPercentage)[[1]]
	ctCallPalx<-performPALXBoth(retList$ctCallDataForPalx, retList$ctValDataForPalx, percentage=curPercentage)[[2]]

	# check which columns and rows were removed by PALX
	colRemove<-colnames(retList$ctValDataForPalx)[which(!(colnames(retList$ctValDataForPalx) %in% colnames(ctValPalx)))]
	samRemove<-c()
	if (length(colRemove)>0)
	{
		print("Columns (samples) removed during PALX:")
		colsToPrint<-cbind(colRemove, determineSampleNamesFromChamberIDInDataList(dataList, colRemove, retList$maxAssay))
		colnames(colsToPrint)<-c("ChamberID","SampleName")
		print(colsToPrint)

		samRemove<-cbind(colRemove, determineSampleNamesFromChamberIDInDataList(dataList, colRemove, retList$maxAssay), rep("Failed PALX", length(colRemove)))
		colnames(samRemove)<-c("ChamberID", "SampleName","ReasonRemoved")
	}
	geneRemove<-c()
	rowRemove<-rownames(retList$ctValDataForPalx)[which(!(rownames(retList$ctValDataForPalx) %in% rownames(ctValPalx)))]
	if (length(rowRemove)>0)
	{
		print("Rows (genes/assays) removed during PALX:")
		rowsToPrint<-convertWellToGeneName(rowRemove, dataList)
		print(rowsToPrint)
		
		geneRemove<-cbind(rowsToPrint, rep("Failed PALX", length(rowRemove)))
		colnames(geneRemove)[3]<-"ReasonRemoved"
	}

	# set 999 and Fail to NA (may want to eventually change 999 to a floored value)
	# need to perform failToNA first (only set fails when Ct!=999) - these are typically very small Ct values (high expression)
	curFloor<-returnParameter(parameters, "floor", "numeric")
	ctValCorrected<-failToNA(ctValPalx, ctCallPalx, curFloor, workdir)
	ctValCorrect<-convert999(ctValCorrected, curFloor)

	# pick HK genes (if there are > 1 plates per assay, then also do QC on HK genes to see if samples show poor correlation)
	# get HK genes from parameter file?
	HKgenes<-unlist(strsplit(returnParameter(parameters, "HKgenes","character"), " "))
	# convert HK genes to well IDs
	retMap<-convertHKtoWell(HKgenes, dataList)
	if (length(rowRemove) > 0)
	{
		# check if any HK genes were removed during PALX
		if (any(rowRemove %in% retMap[,1]))
		{
			# need to remove that HK gene from consideration
			remIndex<-c()
			for (i in 1:nrow(geneRemove))
			{
				if (geneRemove[i,2] %in% retMap[,2])
					remIndex<-c(remIndex, which(geneRemove[i,2]==retMap[,2]))
			}
			retMap<-retMap[-remIndex,]
		}
	}

	# do QC on samples using the HK genes if there is more than one plate for the assay
	{
		if (retList$maxAssay > 1)
		{
			ctValUpdate<-remQuestSam(ctValCorrect, retMap, retList$maxAssay, workdir)
		}
		else
		{
			ctValUpdate<-ctValCorrect
		}
	}
	colRemove<-colnames(ctValCorrect)[which(!(colnames(ctValCorrect) %in% colnames(ctValUpdate)))]
	if (length(colRemove)>0)
	{
		print("Columns (samples) removed because of low correlation across plates for housekeeping genes:")
		colsToPrint<-cbind(colRemove, determineSampleNamesFromChamberIDInDataList(dataList, colRemove, retList$maxAssay))
		colnames(colsToPrint)<-c("ChamberID","SampleName")
		print(colsToPrint)
		
		samRemove2<-cbind(colRemove, determineSampleNamesFromChamberIDInDataList(dataList, colRemove, retList$maxAssay), rep("Low Correlation Across Plates for HK genes", length(colRemove)))
		colnames(samRemove2)<-c("ChamberID", "SampleName","ReasonRemoved")
		samRemove<-rbind(samRemove, samRemove2)
	}

	# now pick the HK genes
	keepHKgenes<-calculateMforHK(ctValUpdate, retMap)
	# print HK genes to the console
	print("The housekeeping genes used in Delta Ct calculations:")
	print(unique(retMap[match(keepHKgenes, retMap[,1]),2]))
	HKnames<-as.matrix(unique(retMap[match(keepHKgenes, retMap[,1]),2]))
	colnames(HKnames)<-"ChosenHK"
	
	# print out samples and genes removed as well as HK genes chosen
	if (length(samRemove) > 0)
		write.csv(samRemove, file=paste(workdir,"Results/removedSamples.csv",sep=""), row.names=FALSE)
	write.csv(HKnames, file=paste(workdir, "Results/chosenHKgenes.csv",sep=""), row.names=FALSE)

	# calculate delta Ct using HK genes that were picked using M value (from Genome Biology paper)
	deltaCt<-calculateDeltaCt(ctValUpdate, keepHKgenes, ctCallPalx)
	
	# next get the reference samples to calculate delta delta Ct - these will be stored in the parameters file
	refSamples<-unlist(strsplit(returnParameter(parameters, "RefSamples","character"), ";"))
	deltadeltaCt<-calculateDeltaDeltaCtWithRefGroup(deltaCt, refSamples, dataList, ctCallPalx)

	# 5/5/14 do extra check to see if using reference samples creates a lot more NAs (if ref samples failed)
	if (any(apply(deltadeltaCt, 1, function(x) {sum(is.na(x))}) > (1-curPercentage)*ncol(deltadeltaCt)))
	{
		# remove this assay
		remIndex<-which(apply(deltadeltaCt, 1, function(x) {sum(is.na(x))}) > (1-curPercentage)*ncol(deltadeltaCt))
		deltadeltaCt<-deltadeltaCt[-remIndex,]
		rowRemove<-rownames(deltaCt)[which(!(rownames(deltaCt) %in% rownames(deltadeltaCt)))]
		
		print("Rows (genes/assays) removed because reference sample failed:")
		rowsToPrint<-convertWellToGeneName(rowRemove, dataList)
		print(rowsToPrint)
		
		geneRemove2<-cbind(rowsToPrint, rep("Failed reference sample", length(rowRemove)))
		colnames(geneRemove2)[3]<-"ReasonRemoved"
		geneRemove<-rbind(geneRemove, geneRemove2)
	}
	if (length(geneRemove) > 0)
		write.csv(geneRemove, file=paste(workdir, "Results/removedGenes.csv",sep=""), row.names=FALSE)

	# now calculate FC
	FC<-2^(-deltadeltaCt)

	# do I need to change the data before writing to csv (like column names)?
	# could probably strip the assay number before the row name
	# for now leave as is and decide after discussion
	write.csv(FC, file=paste(workdir,"Results/FC.csv",sep=""))	

	# write out design files for matching chamber IDs to gene/assay names and sample names
	GeneChamberID<-convertWellToGeneName(rownames(FC), dataList)
	# need to rearrange the names to match the order in FC
	GeneChamberID<-GeneChamberID[match(rownames(FC), GeneChamberID[,1]),]
	#all(GeneChamberID[,1]==rownames(FC))
	# now add a column for number of NAs per gene
	GeneChamberID<-as.data.frame(GeneChamberID)
	GeneChamberID$numNAs<-apply(FC, 1, function(x) {sum(is.na(x))})
	GeneChamberID$numSamples<-rep(ncol(FC), nrow(GeneChamberID))
	GeneChamberID$percentNAs<-as.numeric(as.character(GeneChamberID$numNAs))/as.numeric(as.character(GeneChamberID$numSamples))

	SampleChamberID<-cbind(colnames(FC), determineSampleNamesFromChamberIDInDataList(dataList, colnames(FC), retList$maxAssay))
	colnames(SampleChamberID)<-c("ChamberID","SampleName")
	SampleChamberID<-as.data.frame(SampleChamberID)
	SampleChamberID$numNAs<-apply(FC, 2, function(x) {sum(is.na(x))})
	SampleChamberID$numGenes<-rep(nrow(FC), nrow(SampleChamberID))
	SampleChamberID$percentNAs<-as.numeric(as.character(SampleChamberID$numNAs))/as.numeric(as.character(SampleChamberID$numGenes))

	write.csv(GeneChamberID, file=paste(workdir,"Results/designForGeneChamberID.csv",sep=""), row.names=FALSE)
	write.csv(SampleChamberID, file=paste(workdir,"Results/designForSampleChamberID.csv",sep=""), row.names=FALSE)

	# finally make PCA plot to check for plate batch effect
	makePCAplot(FCvals=FC, workdir, retMap, refSamples, SampleChamberID)
	
	# make PCA plot to look at assays as well
	makePCAplotForAssays(FCvals=FC, workdir, retMap, refSamples, SampleChamberID)
}