# R code to start with Fluidigm Ct data and return delta delta Ct values
#path<-"/Users/ewhalen/Projects/Fluidigm/AutomateTests/"
#directory<-"NewNoemieLupusGenomicCBC"

#workdir<-paste(path, directory, "/", sep="")
params<-commandArgs(trailingOnly = TRUE)
workdir<-params[1]

# read in scripts
source(file=paste(workdir,"scripts/","NewDataDesign_Fluidigm_R_Functions.R",sep=""))
source(file=paste(workdir,"scripts/","NewDataDesign_Fluidigm_R_HKgenes_Functions.R",sep=""))
# use these temporarily while I'm working on the files
#source(file="/Users/ewhalen/Rscripts/NewDataDesign_Fluidigm_R_Functions.R")
#source(file="/Users/ewhalen/Rscripts/NewDataDesign_Fluidigm_R_HKgenes_Functions.R")

packageLoad("combinat")

# read in parameters
parameters<-read.csv(file=paste(workdir,"Parameters/parameters.csv",sep=""))

# now expect there to be 3 files in the Data directory
# read in data
fluidData<-read.csv(paste(workdir,"Data/dataFile.csv",sep=""))
sampleDesignData<-read.csv(paste(workdir,"Data/sampleDesignFile.csv",sep=""))
assayDesignData<-read.csv(paste(workdir,"Data/assayDesignFile.csv",sep=""))

ctVal<-changeFluidFormat(fluidData, colName="Ct.Value")
ctCall<-changeFluidFormat(fluidData, colName="Ct.Call")

# do PALX across all genes and samples
retList<-prepareDataForPalx(ctVal, ctCall)

curPercentage<-returnParameter(parameters, "PALX", "numeric")

ctValPalx<-performPALXBoth(retList$ctCallDataForPalx, retList$ctValDataForPalx, percentage=curPercentage)[[1]]
ctCallPalx<-performPALXBoth(retList$ctCallDataForPalx, retList$ctValDataForPalx, percentage=curPercentage)[[2]]

# check which columns and rows were removed by PALX
colRemove<-colnames(retList$ctValDataForPalx)[which(!(colnames(retList$ctValDataForPalx) %in% colnames(ctValPalx)))]
if (length(colRemove)>0)
{
	print("Columns (samples) removed during PALX:")
	print(colRemove)
}
rowRemove<-rownames(retList$ctValDataForPalx)[which(!(rownames(retList$ctValDataForPalx) %in% rownames(ctValPalx)))]
if (length(rowRemove)>0)
{
	print("Rows (genes/assays) removed during PALX:")
	print(rowRemove)
}

# set 999 and Fail to NA (may want to eventually change 999 to a floored value)
# need to perform failToNA first (only set fails when Ct!=999) - these are typically very small Ct values (high expression)
ctValCorrected<-failToNA(ctValPalx, ctCallPalx)
curFloor<-returnParameter(parameters, "floor", "numeric")
ctValCorrect<-convert999(ctValCorrected, curFloor)

# pick HK genes (if there are > 1 plates per assay, then also do QC on HK genes to see if samples show poor correlation)
# get HK genes from parameter file?
HKgenes<-unique(as.character(assayDesignData$Target[which(assayDesignData$Type=="HK")]))
# convert HK genes to well IDs
retMap<-unique(assayDesignData[which(assayDesignData$Type=="HK"), 1:2])

# do QC on samples using the HK genes if there is more than one plate for the assay
{
	if (retList$maxAssay > 1)
	{
		ctValUpdate<-remQuestSam(ctValCorrect, retMap, retList$maxAssay)
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
	print(colRemove)
}

# now pick the HK genes
keepHKgenes<-calculateMforHK(ctValUpdate, retMap)
# print HK genes to the console
print("The housekeeping genes used in Delta Ct calculations:")
print(as.character(unique(retMap[match(keepHKgenes, retMap[,1]),2])))

# calculate delta Ct using HK genes that were picked using M value (from Genome Biology paper)
deltaCt<-calculateDeltaCt(ctValUpdate, keepHKgenes)
	
# calculate delta delta Ct using the reference samples in the sample design data
deltadeltaCt<-calculateDeltaDeltaCtWithRefGroup(deltaCt, sampleDesignData)

# now calculate FC
FC<-2^(-deltadeltaCt)

# do I need to change the data before writing to csv (like column names)?
# could probably strip the assay number before the row name
# for now leave as is and decide after discussion
write.csv(FC, file=paste(workdir,"Results/FC.csv",sep=""))	

	


