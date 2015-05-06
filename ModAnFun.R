#####
# test that the list of functions has been installed
# if not, then install from CRAN
# finally load libraries
#####
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
	if ("preprocessCore" %in% missing)
	{
		source("http://bioconductor.org/biocLite.R")
		biocLite("preprocessCore")
	}
	if ("impute" %in% missing)
	{
		source("http://bioconductor.org/biocLite.R")
		biocLite("impute")
	}
	if ("edgeR" %in% missing)
	{
		source("http://bioconductor.org/biocLite.R")
		biocLite("edgeR")
	}
	# then load libraries
	for (i in 1:length(x))
	{
		require(x[i], character.only=TRUE)
	}
}

#####
# read the parameters file
#####
readParameters<-function(path)
{
	parameters_name = list.files(paste(path, "/Parameters/", sep=""), pattern=".csv")[1]
	parameters = read.csv(paste(path, "Parameters/", parameters_name, sep=""))
	return(parameters)
}

#####
# trim any leading or trailing spaces in a string when doing a comparison
#####
trimString = function (string) {
  gsub("^\\s+|\\s+$", "", string)
}

#####
# split a string based on a comma
#####
splitString = function(string) {
  unlist(strsplit(as.character(string), "\\,"))
}

#####
# return a specific row/value from the parameters data frame 
#####
returnParameter<-function(parameters, key, type, split=T)
{
	matchIndex<-which(toupper(trimString(parameters$Parameter))==toupper(trimString(key)))
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
# read in csv files
#####
readData<-function(path, project_name, dataFile)
{
	dataPath<-paste(path, "Data/", project_name, "/", dataFile, sep="")
	checkFile<-unlist(strsplit(dataFile, "\\."))
	fileType<-tolower(checkFile[length(checkFile)])
	if (fileType=="csv")
	  	curData<-read.csv(dataPath)
	if (fileType=="txt" | fileType=="tsv")
		curData<-read.delim(dataPath)
	return(curData)
}

#####
# create results directory
#####
createResults<-function(curPath, project_name)
{
	fullPath<-paste(curPath, "Results/", sep="")
	lFiles<-list.files(fullPath)
	fullPath<-paste(fullPath, trimString(project_name), "/", sep="")
	# replace any spaces with underscore in the directory name
	fullPath<-gsub(" ", "_", fullPath)
	# only create the directory if it's not currently there
	if (!(trimString(project_name) %in% trimString(lFiles)))
		dir.create(fullPath)
	return(fullPath)
}

#####
# check that the design data has the right columns
#####
designDataCheck<-function(designData, parameters)
{
	passValue<-TRUE
	missingCols<-c()
	# there are 2 required columns
	requireCols<-c("barcode","group")
	recommendedCols<-c("sample_id", "group_label")
	if (!all(requireCols %in% colnames(designData)))
	{
		missingCols<-c(missingCols, requireCols[!(requireCols %in% colnames(designData))])
		passValue<-FALSE
	}

	# make sure the barcode column will match the column headers from the expression file
	if ("barcode" %in% colnames(designData))
	{
		designData$origBarcode<-designData$barcode
		designData$barcode<-as.character(designData$barcode)
		designData$barcode<-make.names(designData$barcode, unique=TRUE)
	}
	
	# now see if any optional columns are missing
	orderVals<-returnParameter(parameters, "order", "statement", split=T)
	color_groups<-returnParameter(parameters, "color_groups", "statement", split=T)
	checkVals<-unique(c(orderVals, color_groups))
	if (any(checkVals==""))
		checkVals<-checkVals[-which(checkVals=="")]
	if (length(checkVals) > 0)
	{
		if (!all(checkVals %in% colnames(designData)))
		{
			missingCols<-c(missingCols, checkVals[!(checkVals %in% colnames(designData))])
			passValue<-FALSE
		}
	}	
	if (length(missingCols) > 0)
		print(paste("ERROR: The following columns are missing in the design data set (column names are case sensitive): ",
			paste(missingCols, collapse=", "), ".", sep=""))	
	return(list(passValue=passValue, designData=designData))
}

#####
# check that the expression data has the correct columns
#####
expDataCheck<-function(expData, parameters)
{
	platform<-returnParameter(parameters, "platform_type", type="character")
	if (length(platform)==0)
		platform<-"Microarray"
	passValue<-TRUE
	if (length(grep(toupper("Signal"), toupper(colnames(expData))))==0)
	{
		print("ERROR: The expression data should have column names that include 'Signal' to indicate signal intensity.")
		passValue<-FALSE
	}
	if (length(grep(toupper("Pval"), toupper(colnames(expData))))==0 & toupper(platform)==toupper("Microarray"))
	{
		print("ERROR: The expression data should have column names that include 'Pval' to indicate detection p-values.")
		passValue<-FALSE		
	}
	return(passValue)	
}

#####
# decide which method of synchronizing to use based on platform
#####
syncData<-function(expData, designData, moduleAssignData, parameter, rearrange=FALSE)
{
	platform<-returnParameter(parameters, "platform_type", type="character")
	if (length(platform)==0)
		platform<-"Microarray"
	if (toupper(platform)=="MICROARRAY")
		retData<-syncMicroarrayData(expData, designData, moduleAssignData, parameters, rearrange)
	if (toupper(platform)=="RNASEQ")
		retData<-syncRNAseqData(expData, designData, moduleAssignData, parameters)
	return(retData)
}

#####
# synchronize the experimental, design and module data for microarrays
# if rearrange = TRUE, then rearrange the columns so that cases occur first, then controls (and within rearrange by sample_id)
#####
syncMicroarrayData<-function(expData, designData, moduleAssignData, parameters, rearrange=FALSE)
{
	# get the preprocess parameters
	performNormalize<-returnParameter(parameters, "performNormalize", "character")
	performFloor<-returnParameter(parameters, "performFloor", "character")
	performLog2<-returnParameter(parameters, "performLog2", "character")
	performPALO<-returnParameter(parameters, "performPALO", "character")
	PALX<-returnParameter(parameters, "PALX", "numeric")
	dataPass<-TRUE

	# check the data to see if log2 should be performed (override the value stored in parameters)
	sigCols<-grep(toupper("Signal"), toupper(colnames(expData)))
	# 4/17/12 do not override the performLog2 value if it is set to FALSE (the default is TRUE)
	floorValue<-returnParameter(parameters, "min_expr", "numeric")
	expData[,sigCols]<-apply(expData[,sigCols],2, function(x) {as.numeric(as.character(x))})
	if (performLog2==TRUE)
	{
		if (max(expData[,sigCols], na.rm=TRUE) < 50)
		{
			print("WARNING: The maximum intensity is less than 50 so the log2 calculation will not be performed.")
			performLog2<-FALSE
			floorValue<-log2(returnParameter(parameters, "min_expr", "numeric"))
		}
	}
	else
	{
		# performLog2 has been set to FALSE (actively changed from the default), but should warn them if there are any large values
		if (any(expData[,sigCols] > 100, na.rm=TRUE))
			print("WARNING: Some of the intensity values are greater than 100.  Should consider performing log2 calculation.")
	}
		
	# do not reset the floor boolean variable

	# if at least one of the 4 preprocessing steps is TRUE, then call preprocess
	if (any(c(performNormalize, performFloor, performLog2, performPALO)=="TRUE"))
	{
		if (PALX==0)
		{
			PALO<-TRUE
			PALX<-NA
		}
		else
			PALO<-FALSE
		# 10/10/11 for now do not have a range requirement (may want to incorporate into parameter file)
		# but setting the range is available in the preprocess function
		# added range to parameters file 12/5/11
		rangeCut<-returnParameter(parameters, "range_cut", "numeric")
		subexpData<-preprocess(expData, outputFile=NA, floorValue=floorValue, PALO=PALO, PALX=PALX, 
			performNormalize=performNormalize, performFloor=performFloor, performLog2=performLog2, 
			performPALO=performPALO, rangeCutoff=rangeCut)
		if (nrow(subexpData)==0)
		{
			print(paste("ERROR: No data was returned from preprocessing.  Either no probes passed PALX or no probes passed the range cutoff of ", rangeCut, ".", sep=""))
			dataPass<-FALSE
			
		}
	}
	else
		subexpData<-expData
			
	# subset the columns of the experimental data based on the design
	matchIndex<-match(toupper(designData$uniqueID), toupper(colnames(subexpData)))
	if (any(is.na(matchIndex)))
		matchIndex<-matchIndex[-which(is.na(matchIndex))]
	# add column one to include the probe id
	subexpData<-subexpData[,c(1,matchIndex)]
	# make sure the first column is labeled PROBE_ID
	colnames(subexpData)[1]<-"PROBE_ID"

	# 4/16/12 after subsetting to only samples in the design file, check the probes variability again
	# remove any probe with no variability
	checkProbeVar<-apply(subexpData[,2:ncol(subexpData)], 1, function(x) {length(unique(x))})
	# length of 1 means the entire row has only one value
	remIndex<-which(checkProbeVar==1)
	if (length(remIndex)>0)
		subexpData<-subexpData[-remIndex,]

	if (rearrange==TRUE)
	{
		ssubexpData<-subexpData[,-1]
		matchIndex<-match(toupper(colnames(ssubexpData)), toupper(designData$uniqueID))
		if (any(is.na(matchIndex)))
			matchIndex<-matchIndex[-which(is.na(matchIndex))]
		curGroups<-designData$group[matchIndex]
		curSampleIDs<-designData$sample_id[matchIndex]
		orderIndex<-order(curGroups, curSampleIDs)
		ssubexpData<-ssubexpData[,orderIndex]
		subexpData<-cbind(subexpData$PROBE_ID, ssubexpData)
		colnames(subexpData)[1]<-"PROBE_ID"
	}
	
	# keep all the probes for running limma
	origSubExpData<-subexpData
		
	# subset experimental data based on the transcripts in the modules
	matchIndex<-match(toupper(moduleAssignData$Systematic.Name), toupper(subexpData$PROBE_ID))
	if (any(is.na(matchIndex)))
		matchIndex<-matchIndex[-which(is.na(matchIndex))]
	if (dataPass==TRUE & length(matchIndex)==0)
		print("ERROR: The module version file does not match the probe IDs in the expression file.")
	subexpData<-subexpData[matchIndex,]

	# add the module information to the experimental data
	# may not have all the probes
	matchIndex<-match(toupper(subexpData$PROBE_ID), toupper(moduleAssignData$Systematic.Name))
	if (any(is.na(matchIndex)))
		matchIndex<-matchIndex[-which(is.na(matchIndex))]
	subexpData$Module<-as.character(moduleAssignData$module[matchIndex])
	
	return(list(subexpData=subexpData, origSubExpData=origSubExpData))	
}

#####
# synchronize the experimental, design and module data for RNAseq data
#####
syncRNAseqData<-function(expData, designData, moduleAssignData, parameters)
{
	countData<-expData[,grep("Signal", colnames(expData))]
	rownames(countData)<-expData$ProbeID
	roundCountData<-t(apply(countData, 1, round))
	
	if (any(roundCountData != countData))
	{
		print("WARNING: the RNA seq count data is not integers. The values will be rounded.")
		countData<-roundCountData
	}
		
	# create DGE object
	d<-DGEList(counts=countData)

	# perform TMM normalization (perform normalization before filtering)
	d<-calcNormFactors(d)

	# subset to samples of interest
	dSub<-d[,match(designData$uniqueID, colnames(d))]
	
	# now filter rows (genes) based on those that have a count of 1 in at least PALX% sample in normalized data
	PALX<-returnParameter(parameters, "PALX", "numeric")
	keepRows<-rowSums(round(cpm(dSub$counts)) > 1) >= ceiling(PALX*ncol(dSub))
	dSubFilter<-dSub[keepRows,]
	
	# now subset to only genes in the modules
	matchIndex<-match(moduleAssignData$Systematic.Name, rownames(dSubFilter))
	if (any(is.na(matchIndex)))
		matchIndex<-matchIndex[-which(is.na(matchIndex))]
	
	dSubFilterMod<-dSubFilter[matchIndex,]		
	
	# test
	# all(rownames(dSubFilterMod) %in% moduleAssignData$Systematic.Name)
	
	return(list(subexpData=dSubFilterMod, origSubExpData=dSubFilter))		
}

#####
# remove modules that don't have any probes in the expression data
#####
updateModuleAssignData<-function(subexpData, moduleAssignData, parameters)
{
	platform<-returnParameter(parameters, "platform_type", type="character")
	if (length(platform)==0)
		platform<-"Microarray"
		
	if (toupper(platform)=="MICROARRAY")
	{
		subexpData$Module<-as.character(subexpData$Module)
		modsInData<-unique(subexpData$Module)
	}
	if (toupper(platform)=="RNASEQ")
	{
		modsInData<-unique(as.character(moduleAssignData$module[match(rownames(subexpData), moduleAssignData$Systematic.Name)]))
	}
	keepIndex<-match(moduleAssignData$module, modsInData)
	keepIndex<-which(!is.na(keepIndex))
	moduleAssignData<-moduleAssignData[keepIndex,]
	moduleAssignData$module<-as.character(moduleAssignData$module)
	return(moduleAssignData)
}

#####
# remove modules that don't have any probes in the expression data
#####
updateModuleAnnotateData<-function(subexpData, moduleAnnotateData, parameters)
{
	platform<-returnParameter(parameters, "platform_type", type="character")
	if (length(platform)==0)
		platform<-"Microarray"

	if (toupper(platform)=="MICROARRAY")
	{
		subexpData$Module<-as.character(subexpData$Module)
		modsInData<-unique(subexpData$Module)
	}
	if (toupper(platform)=="RNASEQ")
	{
		modsInData<-unique(as.character(moduleAssignData$module[match(rownames(subexpData), moduleAssignData$Systematic.Name)]))
	}
	keepIndex<-match(moduleAnnotateData$Module, modsInData)
	keepIndex<-which(!is.na(keepIndex))
	moduleAnnotateData<-moduleAnnotateData[keepIndex,]
	moduleAnnotateData$Module<-as.character(moduleAnnotateData$Module)
	return(moduleAnnotateData)
}

#####
# calculate summary statistics (mean, sd, median, and mad) per sample and module
#####
calcSummaryStats<-function(subexpData, parameters, moduleAssignData)
{
	platform<-returnParameter(parameters, "platform_type", type="character")
	if (length(platform)==0)
		platform<-"Microarray"

	if (toupper(platform)=="MICROARRAY")
	{
		# remove the probe and module columns
		remIndex<-match(toupper(c("PROBE_ID", "Module")), toupper(colnames(subexpData)))
		# put the values back on the raw intensity scale before calculating summary stats
		moduleMean<-apply(subexpData[,-remIndex], 2, function(x) {tapply(2^x, subexpData$Module, mean, na.rm=TRUE)})
		moduleSd<-apply(subexpData[,-remIndex], 2, function(x) {tapply(2^x, subexpData$Module, sd, na.rm=TRUE)})
		moduleMedian<-apply(subexpData[,-remIndex], 2, function(x) {tapply(2^x, subexpData$Module, median, na.rm=TRUE)})
		moduleMad<-apply(subexpData[,-remIndex], 2, function(x) {tapply(2^x, subexpData$Module, mad, na.rm=TRUE)})
	}
	if (toupper(platform)=="RNASEQ")
	{
		# get normalized counts	
		seData<-cpm(subexpData$counts)
		curMods<-as.character(moduleAssignData$module[match(rownames(subexpData), moduleAssignData$Systematic.Name)])
		moduleMean<-apply(seData, 2, function(x) {tapply(x, curMods, mean, na.rm=TRUE)})
		moduleSd<-apply(seData, 2, function(x) {tapply(x, curMods, sd, na.rm=TRUE)})
		moduleMedian<-apply(seData, 2, function(x) {tapply(x, curMods, median, na.rm=TRUE)})
		moduleMad<-apply(seData, 2, function(x) {tapply(x, curMods, mad, na.rm=TRUE)})
	}
	
	return(list(moduleMean=moduleMean, moduleSd=moduleSd, moduleMedian=moduleMedian, moduleMad=moduleMad))
}

#####
# compare case to control (for now just assume 2 groups and 1 time point)
# 4/13/11: now can have multiple reference groups - NEED TO ADD ERROR CHECKING!!
#####
compareCaseControl<-function(parameters, subexpData, designData, moduleAssignData, origSubExpData, geneSymbolData)
{
	platform<-returnParameter(parameters, "platform_type", type="character")
	if (length(platform)==0)
		platform<-"Microarray"

	if (toupper(platform)=="MICROARRAY")
	{
		retList<-compareCaseControlMicroarray(parameters, subexpData, designData, moduleAssignData, origSubExpData, geneSymbolData)
	}
	if (toupper(platform)=="RNASEQ")
	{
		retList<-compareCaseControlRNAseq(parameters, subexpData, designData, moduleAssignData, origSubExpData, geneSymbolData)
	}
	return(retList)
}

#####
# get results for RNAseq - should mainly match microarray
#####
compareCaseControlRNAseq<-function(parameters, subexpData, designData, moduleAssignData, origSubExpData, geneSymbolData)
{
	# assume we just have one (2 sample) group comparison to make
	indivRes<-calcIndivMetricsRNAseq(parameters, subexpData, designData)
	
	# run edgeR model at gene level
	modelRes<-calcEdgeRmodel(parameters, origSubExpData, designData, moduleAssignData)

	# match order of microarray results
	DEG<-modelRes[,c(1,3:7)]
	DEG$probeID<-rownames(DEG)
	DEG$DEGs<-rep(0, nrow(DEG))
	alpha<-returnParameter(parameters, "FDR", "numeric")	
	MTC<-returnParameter(parameters, "MTC", "character")
	if (MTC==TRUE)
		DEG$DEGs[which(DEG$FDR < alpha)]<-1
	if (MTC==FALSE)
	{
		DEG$DEGs[which(DEG$PValue < alpha)]<-1
		DEG$FDR<-rep(NA, nrow(DEG))
	}
	# need to match the gene symbol to the data
	matIndex<-match(DEG$probeID, geneSymbolData$probe_id)
	if (all(is.na(matIndex)))
		DEG$geneSymbol<-rep(NA, nrow(DEG))
	else
		DEG$geneSymbol<-geneSymbolData$gene_symbol[matIndex]
	# now rearrange, rename, rescale columns
	DEG$logFC<-2^DEG$logFC
	DEG<-DEG[,c(7,6,9,2,1,5,3,4,8)]
	colnames(DEG)<-c("probeID","module","geneSymbol","LR","FC","Diff","originalPval","adjustedPval","DEGs")
	
	# now need to summarize this data
	moduleResults<-calcModuleMetrics(indivRes, parameters, moduleAssignData, modelRes, subexpData)
	
	percentMetrics<-calculatePercentMetric(moduleResults$posPercentZscore, moduleResults$negPercentZscore, moduleResults$posPercentFC, moduleResults$negPercentFC, parameters)
	
	# only missing the GSA results
	return(list(delta=indivRes$delta, zscore=indivRes$zscore, fc=indivRes$fc, posPercentZscore=moduleResults$posPercentZscore,
		negPercentZscore=moduleResults$negPercentZscore, posPercentFC=moduleResults$posPercentFC, 
		negPercentFC=moduleResults$negPercentFC, medianZscore=moduleResults$medianZscore, medianFC=moduleResults$medianFC,
		posCountZscore=moduleResults$posCountZscore, negCountZscore=moduleResults$negCountZscore, 
		posCountFC=moduleResults$posCountFC, negCountFC=moduleResults$negCountFC, 
		maxPercentZscore=percentMetrics$maxPercentZscore, diffPercentZscore=percentMetrics$diffPercentZscore,
		valenceZscore=percentMetrics$valenceZscore, maxPercentFC=percentMetrics$maxPercentFC, 
		diffPercentFC=percentMetrics$diffPercentFC, valenceFC=percentMetrics$valenceFC, moduleOrder=list(moduleResults$topMod),
		tPercentList=list(moduleResults$tPercent), tBothPercentList=list(moduleResults$tBothPercent), 
		tMeanList=list(moduleResults$origTmean), DEGs=list(DEG)))
}

######
# roll gene level stats into module level
######
calcModuleMetrics<-function(indivRes, parameters, moduleAssignData, modelRes, subexpData)
{
	zscore<-indivRes$zscore
	delta<-indivRes$delta
	fc<-indivRes$fc
	zscoreThresh<-returnParameter(parameters, "zscore_cut", "numeric")
	deltaThresh<-returnParameter(parameters, "delta_cut", "numeric")
	fcThresh<-returnParameter(parameters, "fold_cut", "numeric")
		
	# apply thresholds to score and delta to count the positive and negative transcripts
	posMatZscore<-((zscore > zscoreThresh) & (abs(delta) > deltaThresh))
	negMatZscore<-((zscore < -zscoreThresh) & (abs(delta) > deltaThresh))
	posMatFC<-((fc > fcThresh) & (abs(delta) > deltaThresh))
	negMatFC<-((fc < 1/fcThresh) & (abs(delta) > deltaThresh))

	groupFC<-2^modelRes$logFC[which(!is.na(modelRes$module))]
	names(groupFC)<-rownames(modelRes)[which(!is.na(modelRes$module))]
	# need to reorder groupFC to match those of subexpData
	matchIndex<-match(rownames(subexpData), names(groupFC))
	groupFC<-groupFC[matchIndex]
	all(names(groupFC)==rownames(subexpData))	# test
	allPosVecFC<-(groupFC > fcThresh)
	allNegVecFC<-(groupFC < 1/fcThresh)
	
	# now repeat with FDR
	MTC<-returnParameter(parameters, "MTC", "character")
	if (MTC==TRUE)	
		modFDR<-modelRes$FDR[which(!is.na(modelRes$module))]
	if (MTC==FALSE)
		modFDR<-modelRes$PValue[which(!is.na(modelRes$module))]
	names(modFDR)<-rownames(modelRes)[which(!is.na(modelRes$module))]
	modFDR<-modFDR[matchIndex]
	
	posMatZscore<-replace(posMatZscore, which(posMatZscore==TRUE), 1)
	posMatZscore<-replace(posMatZscore, which(posMatZscore==FALSE), 0)
	posMatFC<-replace(posMatFC, which(posMatFC==TRUE), 1)
	posMatFC<-replace(posMatFC, which(posMatFC==FALSE), 0)
	
	negMatZscore<-replace(negMatZscore, which(negMatZscore==TRUE), 1)
	negMatZscore<-replace(negMatZscore, which(negMatZscore==FALSE), 0)	
	negMatFC<-replace(negMatFC, which(negMatFC==TRUE), 1)
	negMatFC<-replace(negMatFC, which(negMatFC==FALSE), 0)
	
	# combine counts into modules
	all(rownames(subexpData)==rownames(negMatFC))	# make sure the order is the same for modules

	curMods<-moduleAssignData$module[match(rownames(subexpData$counts), moduleAssignData$Systematic.Name)]
	posCountZscore<-apply(posMatZscore, 2, function(x) {tapply(x, curMods, sum, na.rm=TRUE)})
	negCountZscore<-apply(negMatZscore, 2, function(x) {tapply(x, curMods, sum, na.rm=TRUE)})
	posCountFC<-apply(posMatFC, 2, function(x) {tapply(x, curMods, sum, na.rm=TRUE)})
	negCountFC<-apply(negMatFC, 2, function(x) {tapply(x, curMods, sum, na.rm=TRUE)})

	# combine all counts into modules
	allPosCountFC<-apply(as.matrix(allPosVecFC), 2, function(x) {tapply(x, curMods, sum, na.rm=TRUE)})
	colnames(allPosCountFC)<-"allFCPos"
	allNegCountFC<-apply(as.matrix(allNegVecFC), 2, function(x) {tapply(x, curMods, sum, na.rm=TRUE)})
	colnames(allNegCountFC)<-"allFCNeg"	
	
	# FIX EVERYTHING BELOW - fixed
	tstats<-log2(groupFC)	# use logFC to see if t-stat would have been up or down
	paraFDR<-returnParameter(parameters, "FDR", "numeric")
	if (paraFDR < 0 | paraFDR > 1)
		paraFDR<-0.05
	tstats[which(modFDR > paraFDR)]<-0
		
	tPosCount<-tapply(tstats, curMods, function(x){sum(x>0, na.rm=TRUE)})
	tNegCount<-tapply(tstats, curMods, function(x){sum(x<0, na.rm=TRUE)})
	# now combine into one tCount
	tBoth<-cbind(tPosCount, -tNegCount)
	colnames(tBoth)<-c("posPercent","negPercent")
	# take the difference rather than the maximum
	tCount<-tPosCount - tNegCount

	# use LR to replace t-stats, though LR is only positive
	origTstats<-modelRes$LR[which(!is.na(modelRes$module))]
	names(origTstats)<-rownames(modelRes)[which(!is.na(modelRes$module))]
	# need to reorder groupFC to match those of subexpData
	matchIndex<-match(rownames(subexpData), names(origTstats))
	origTstats<-origTstats[matchIndex]
	origTsum<-tapply(origTstats, curMods, sum, na.rm=TRUE)
	
	# now get the module counts
	# get the module counts from the probe assignment dataset
	# if used the expression data, could possibly miss some probes that were removed
	modCount<-table(moduleAssignData$module)
	
	if (all(modCount>0))
	{
	 	# calculte the positive and negative percentages
		posPercentZscore<-apply(posCountZscore, 2, function(x) {x/modCount})*100
		negPercentZscore<-apply(negCountZscore, 2, function(x) {x/modCount})*100
		posPercentFC<-apply(posCountFC, 2, function(x) {x/modCount})*100
		negPercentFC<-apply(negCountFC, 2, function(x) {x/modCount})*100
		allPosPercentFC<-apply(allPosCountFC, 2, function(x) {x/modCount})*100
		allNegPercentFC<-apply(allNegCountFC, 2, function(x) {x/modCount})*100
		tPercent<-tCount/modCount*100
		tBothPercent<-apply(tBoth, 2, function(x) {x/modCount})*100
		origTmean<-origTsum/modCount
		
		# now bind t sum with modules to find top modules
		topMod<-cbind(names(modCount), modCount, origTsum, origTmean, match(1:length(origTmean), order(abs(origTmean),decreasing=T)), allPosPercentFC, allNegPercentFC, allPosPercentFC-allNegPercentFC)	
		colnames(topMod)<-c("Module", "moduleCount", "LRStatisticSum", "LRStatisticMean", "Rank", "allFCPercentPos", "allFCPercentNeg", "allFCPercentDiff")
		topMod<-as.data.frame(topMod)
		topMod$Module<-as.character(topMod$Module)
		for (j in 2:ncol(topMod))
		{
			topMod[,j]<-as.numeric(as.character(topMod[,j]))
		}
	}
	
	# also calculate the median zscore & fold change for each module
	medianZscore<-apply(zscore, 2, function(x) {tapply(x, curMods, median, na.rm=TRUE)})
	medianFC<-apply(fc, 2, function(x) {tapply(x, curMods, median, na.rm=TRUE)})
	
	return(list(posPercentZscore=posPercentZscore, negPercentZscore=negPercentZscore, medianZscore=medianZscore,
			posPercentFC=posPercentFC, negPercentFC=negPercentFC, medianFC=medianFC,
	        posCountZscore=posCountZscore, negCountZscore=negCountZscore, posCountFC=posCountFC,
	        negCountFC=negCountFC, topMod=topMod, tPercent=tPercent, 
	        tBothPercent=tBothPercent, origTmean=origTmean))
}

######
# run the edgeR model to compare case vs. control
######
calcEdgeRmodel<-function(parameters, origSubExpData, designData, moduleAssignData)
{
	# need to add design data to the DGEList object
	if (all(designData$uniqueID == rownames(origSubExpData$samples)))
	{
		# should be true since we subset the original data based on design
		origSubExpData$samples$group<-designData$group
		
		# calculate difference using normalized counts
		normVals<-cpm(origSubExpData)
		caseSams<-designData$uniqueID[which(designData$group==1)]
		controlSams<-designData$uniqueID[which(designData$group==0)]
		caseNorms<-normVals[,match(caseSams, colnames(normVals))]
		controlNorms<-normVals[,match(controlSams, colnames(normVals))]
		normDiff<-apply(caseNorms, 1, mean, na.rm=T)-apply(controlNorms, 1, mean, na.rm=T)
		
		# make model matrix
		curMM<-model.matrix(~group, data=designData)
		
		# estimate common dispersion
		origSubExpData<-estimateGLMCommonDisp(origSubExpData, curMM)
		origSubExpData<-estimateGLMTrendedDisp(origSubExpData, curMM)
		origSubExpData<-estimateGLMTagwiseDisp(origSubExpData, curMM)

		curFit<-glmFit(origSubExpData, curMM)
		# case vs. control is the last column in the design matrix
		curLRT<-glmLRT(curFit, coef=ncol(curMM))
		
		tt<-topTags(curLRT, adjust.method="BH", n=nrow(origSubExpData$counts))$table
		
		# now add module information and difference values
		normDiff<-normDiff[match(rownames(tt), names(normDiff))]
		tt$normDiff<-normDiff
		
		# now add module information
		tt$module<-moduleAssignData$module[match(rownames(tt), moduleAssignData$Systematic.Name)]
		
		return(tt)
	}
}

######
# these calculations will need further thought to determine reasonably ways to calculate stable metrics
######
calcIndivMetricsRNAseq<-function(parameters, subexpData, designData)
{
	# get normalized data with genes assigned to modules only
	normSEData<-cpm(subexpData)
	
	# now get the controls
	contSams<-designData$uniqueID[which(designData$group==0)]
	contData<-normSEData[,match(contSams, colnames(normSEData))]
	
	# zmeasure is either mean or median
	zmeasure<-returnParameter(parameters, "zmeasure", "character")
	
	# take measure of central tendency and variation
	if (tolower(zmeasure)=="mean")
	{
		contMeans<-apply(contData, 1, mean, na.rm=T)
		contSD<-apply(contData, 1, sd, na.rm=T)
	}
	if (tolower(zmeasure)=="median")
	{
		contMeans<-apply(contData, 1, median, na.rm=T)
		contSD<-apply(contData, 1, mad, na.rm=T)
	}
	if (any(contSD==0))
	{
		# replace with minimum non-zero SD value
		replaceVal<-min(contSD[which(contSD!=0)], na.rm=T)
		contSD[which(contSD==0)]<-replaceVal
	}

	# now can calculate delta, zscore and FC
	# delta and FC are on the original scale 
	# previously (with microarray) calculated zscore on log2 scale, but it's not stable on RNAseq data because it's not floored
	delta<-apply(normSEData, 2, function(x)
	{
		x-contMeans
	})
	# calculated on original scale - not a good metric for RNAseq data
	zscore<-apply(normSEData, 2, function(x)
	{
		(x-contMeans)/contSD
	})
	# don't allow zero means for FC calculation
	if (any(contMeans==0))
	{
		# replace with minimum non-zero SD value
		replaceVal<-min(contMeans[which(contMeans!=0)], na.rm=T)
		contMeans[which(contMeans==0)]<-replaceVal
	}
	fc<-apply(normSEData, 2, function(x)
	{
		x/contMeans
	})
	return(list(delta=delta, zscore=zscore, fc=fc))
}

#####
# code to compare case and control in microarray data
#####
compareCaseControlMicroarray<-function(parameters, subexpData, designData, moduleAssignData, origSubExpData, geneSymbolData)
{
	# remove the trasncript and module information before doing calculations on the expression values
	remIndex<-match(toupper(c("PROBE_ID", "Module")), toupper(colnames(subexpData)))
	seData<-subexpData[,-remIndex]
	# but set the row names to be PROBE IDs
	rownames(seData)<-as.character(subexpData$PROBE_ID)
	
	# repeat the same steps on the original probe data
	remIndex<-match("PROBE_ID", toupper(colnames(origSubExpData)))
	oseData<-origSubExpData[,-remIndex]
	rownames(oseData)<-as.character(origSubExpData$PROBE_ID)
	
	# get statistics to calculate
	# zmeasure is either mean or median
	zmeasure<-returnParameter(parameters, "zmeasure", "character")
		
	# get thresholds
	curDeltaThresh<-returnParameter(parameters, "delta_cut", "numeric")
	curZscoreThresh<-returnParameter(parameters, "zscore_cut", "numeric")
	curFCThresh<-returnParameter(parameters, "fold_cut", "numeric")
		
	# check for a Subset column in designData (add if not there)
	if (!(toupper("Subset") %in% toupper(colnames(designData))))
	{
		designData$Subset<-rep(1, nrow(designData))
	}
	
	# for now check that there are only 2 groups to compare and that they represent case or non-reference (1) and reference (0)
	if (all(sort(unique(designData$group))==c(0,1)))
	{
		Iter<-length(unique(designData$Subset))
		
		allDelta<-c()
		allZscore<-c()
		allFC<-c()
		allPosPercentZscore<-c()
		allNegPercentZscore<-c()
		allPosPercentFC<-c()
		allNegPercentFC<-c()
		allMedianZscore<-c()
		allMedianFC<-c()
		allPosCountZscore<-c()
		allNegCountZscore<-c()
		allPosCountFC<-c()
		allNegCountFC<-c()		
		allTopMod<-list()
		allMaxPercentZscore<-c()
		allDiffPercentZscore<-c()
		allValenceZscore<-c()
		allMaxPercentFC<-c()
		allDiffPercentFC<-c()
		allValenceFC<-c()

		allTPercents<-list()
		allTBothPercents<-list()
		allTMeans<-list()
		allDEGs<-list()
		allGSA<-list()
		for (i in 1:Iter)
		{		
			curVal<-unique(designData$Subset)[i]
			# now subset to just that group
			subDdata<-designData[which(designData$Subset==curVal),]
			subSEdata<-seData[,match(toupper(subDdata$uniqueID), toupper(colnames(seData)))]
			subOSEdata<-oseData[,match(toupper(subDdata$uniqueID), toupper(colnames(oseData)))]
						
			if (all(sort(unique(subDdata$group))==c(0,1)))
			{
				# calculate the score (zscore or fold change) and delta measures
				deltascore<-calculateDeltaScore(subSEdata, subDdata, zmeasure, subOSEdata)
				allDelta<-cbind(allDelta, deltascore$delta)
				allZscore<-cbind(allZscore, deltascore$zscore)
				allFC<-cbind(allFC, deltascore$fc)
				DEGmat<-cbind(deltascore$origTstats, deltascore$groupFC, deltascore$groupDiff, deltascore$origPvals, deltascore$adjustPvals)
				DEGmat<-as.data.frame(DEGmat)
				colnames(DEGmat)<-c("tStatistic","FC","Diff","originalPval","adjustedPval")
				# add column to indicate statistical significance
				DEGmat$DEGs<-deltascore$tstats
				DEGmat$DEGs[which(DEGmat$DEGs!=0)]<-1
				# add module assignment to the matrix
				DEGmat$module<-moduleAssignData$module[match(rownames(DEGmat), moduleAssignData$Systematic.Name)]
				DEGmat$probeID<-rownames(DEGmat)
				# need to match the gene symbol to the data
				matIndex<-match(DEGmat$probeID, geneSymbolData$probe_id)
				if (all(is.na(matIndex)))
					DEGmat$geneSymbol<-rep(NA, nrow(DEGmat))
				else
					DEGmat$geneSymbol<-geneSymbolData$gene_symbol[matIndex]

				DEGmat<-DEGmat[,c(8,7,9,1:6)]
				allDEGs[[i]]<-DEGmat
				
				# subset results to only probes assigned to modules
				matchIndex<-match(rownames(deltascore$delta), names(deltascore$tstats))
						
				# now apply thresholds to score and delta
				posneg<-calculatePosNeg(deltascore$delta, deltascore$zscore, deltascore$fc, deltascore$tstats[matchIndex], curDeltaThresh, curZscoreThresh, curFCThresh, subexpData, moduleAssignData, deltascore$origTstats[matchIndex], deltascore$groupFC[matchIndex], parameters)
				allPosPercentZscore<-cbind(allPosPercentZscore, posneg$posPercentZscore)
				allNegPercentZscore<-cbind(allNegPercentZscore, posneg$negPercentZscore)
				allPosPercentFC<-cbind(allPosPercentFC, posneg$posPercentFC)
				allNegPercentFC<-cbind(allNegPercentFC, posneg$negPercentFC)
				allMedianZscore<-cbind(allMedianZscore, posneg$medianZscore)
				allMedianFC<-cbind(allMedianFC, posneg$medianFC)
				allPosCountZscore<-cbind(allPosCountZscore, posneg$posCountZscore)
				allNegCountZscore<-cbind(allNegCountZscore, posneg$negCountZscore)
				allPosCountFC<-cbind(allPosCountFC, posneg$posCountFC)
				allNegCountFC<-cbind(allNegCountFC, posneg$negCountFC)				
				# these are different in how to combine multiple subsets
				allTopMod[[i]]<-posneg$topMod
				allTPercents[[i]]<-posneg$tPercent
				allTBothPercents[[i]]<-posneg$tBothPercent
				allTMeans[[i]]<-posneg$origTmean
			
				# calculate difference, valence, and maximum percentages
				percentMetrics<-calculatePercentMetric(posneg$posPercentZscore, posneg$negPercentZscore, posneg$posPercentFC,
					posneg$negPercentFC, parameters)

				allMaxPercentZscore<-cbind(allMaxPercentZscore, percentMetrics$maxPercentZscore)
				allDiffPercentZscore<-cbind(allDiffPercentZscore, percentMetrics$diffPercentZscore)
				allValenceZscore<-cbind(allValenceZscore, percentMetrics$valenceZscore)
				allMaxPercentFC<-cbind(allMaxPercentFC, percentMetrics$maxPercentFC)
				allDiffPercentFC<-cbind(allDiffPercentFC, percentMetrics$diffPercentFC)
				allValenceFC<-cbind(allValenceFC, percentMetrics$valenceFC)
				
				# run GSA
				uniModules<-unique(moduleAssignData$module)
				moduleList<-list()
				for (j in 1:length(uniModules))
				{
					curModule<-uniModules[j]
					curProbes<-moduleAssignData$Systematic.Name[which(moduleAssignData$module==curModule)]
					moduleList[[j]]<-curProbes
				}
				curGSA<-GSA(x=as.matrix(subSEdata), y=subDdata$group+1, genenames=subexpData$PROBE_ID, genesets=moduleList, resp.type="Two class unpaired", 
					nperms=1000, random.seed=81809, minsize=9)
				curGSAmat<-cbind(curGSA$GSA.scores, curGSA$pvalues.lo, curGSA$pvalues.hi)
				curGSAmat<-as.data.frame(curGSAmat)
				rownames(curGSAmat)<-as.character(uniModules)
				colnames(curGSAmat)<-c("GSAscore", "pvalLo", "pvalHi")
				curGSAmat$smallPval<-curGSAmat$pvalHi
				curGSAmat$smallPval[which(curGSAmat$pvalLo < curGSAmat$pvalHi)]<-curGSAmat$pvalLo[which(curGSAmat$pvalLo < curGSAmat$pvalHi)]
				allGSA[[i]]<-curGSAmat
			}
			else
				print(paste("ERROR: Subset group", curVal, "did not have case vs. control groups."))
		}
		# return delta, score, posPercent, negPercent, medianscore, maxPercent, diffPercent, valence
		return(list(delta=allDelta, zscore=allZscore, fc=allFC, posPercentZscore=allPosPercentZscore, 
				negPercentZscore=allNegPercentZscore, posPercentFC=allPosPercentFC, negPercentFC=allNegPercentFC, 
				medianZscore=allMedianZscore, medianFC=allMedianFC, posCountZscore=allPosCountZscore, 
				negCountZscore=allNegCountZscore, posCountFC=allPosCountFC, negCountFC=allNegCountFC,
				maxPercentZscore=allMaxPercentZscore, diffPercentZscore=allDiffPercentZscore, valenceZscore=allValenceZscore,
				maxPercentFC=allMaxPercentFC, diffPercentFC=allDiffPercentFC, valenceFC=allValenceFC, 
				moduleOrder=allTopMod, tPercentList=allTPercents,
				tBothPercentList=allTBothPercents, tMeanList=allTMeans, DEGs=allDEGs, GSAList=allGSA))
	}
}

#####
# calculate delta and score (zscore/fold change)
#####
calculateDeltaScore<-function(seData, designData, zmeasure, oseData)
{
	# split into case (1) and control (0)
	# calculate the fold change on the original expression values
	matchIndex<-match(toupper(designData$uniqueID[which(designData$group==0)]), toupper(colnames(seData)))
	if (any(is.na(matchIndex)))
		matchIndex<-matchIndex[-which(is.na(matchIndex))]
	control<-seData[,matchIndex]
	matchIndex<-match(toupper(designData$uniqueID[which(designData$group==1)]), toupper(colnames(seData)))
	if (any(is.na(matchIndex)))
		matchIndex<-matchIndex[-which(is.na(matchIndex))]
	case<-seData[,matchIndex]
	
	# repeat for all probes
	matchIndex<-match(toupper(designData$uniqueID[which(designData$group==0)]), toupper(colnames(oseData)))
	if (any(is.na(matchIndex)))
		matchIndex<-matchIndex[-which(is.na(matchIndex))]
	controlO<-oseData[,matchIndex]
	matchIndex<-match(toupper(designData$uniqueID[which(designData$group==1)]), toupper(colnames(oseData)))
	if (any(is.na(matchIndex)))
		matchIndex<-matchIndex[-which(is.na(matchIndex))]
	caseO<-oseData[,matchIndex]

	# perform t-tests to compare case and control
	nco<-ncol(control)
	nca<-ncol(case)
	all<-cbind(case, control)
	ncoO<-ncol(controlO)
	ncaO<-ncol(caseO)
	allO<-cbind(caseO, controlO)
	
	# use all probes that pass PALX, not just probes assigned to modules
	designMat<-cbind(rep(1,ncol(allO)), c(rep(1,ncaO),rep(0,ncoO)))
	colnames(designMat)<-c("Control","CasevsControl")
	curFit<-lmFit(allO, designMat)
	curFit<-eBayes(curFit)
	
	origTstats<-curFit$t[,2]
	origPvals<-curFit$p.value[,2]
	curTT<-topTable(curFit, coef="CasevsControl", adjust="BH", number=nrow(allO))
	# need to reorder the top table to match the original ordering
	adjustP<-curTT$adj.P.Val[match(row.names(allO), rownames(curTT))]
	names(adjustP)<-row.names(allO)
	
	# now do multiple testing correction if MTC is TRUE
	MTC<-returnParameter(parameters, "MTC", "character")
	if (MTC==TRUE)
	{
		FDR<-returnParameter(parameters, "FDR", "numeric")
		# make sure value is valid
		if (FDR < 0 | FDR > 1)
			FDR<-0.05
		tstats<-origTstats
		tstats[which(adjustP > FDR)]<-0
		adjustPvals<-adjustP
	}
	else
	{
		tstats<-origTstats
		tstats[which(origPvals > 0.05)]<-0
		adjustPvals<-NA
	}
	
	# first calculate zscore on log2 scale
	# zmeasure will be either "mean" or "median"
	if (tolower(zmeasure)=="mean")
	{
		comean<-apply(control, 1, zmeasure, na.rm=TRUE)
		cosd<-apply(control, 1, sd, na.rm=TRUE)
	}
	if (tolower(zmeasure)=="median")
	{
		comean<-apply(control, 1, zmeasure, na.rm=TRUE)
		cosd<-apply(control, 1, mad, na.rm=TRUE)
	}
	# some control sd are zero - need to fix
	if (any(cosd==0))
	{
		replaceVal<-min(cosd[which(cosd!=0)], na.rm=T)
		cosd[which(cosd==0)]<-replaceVal
	}
	curZscore<-apply(seData, 2, function(x){ (x-comean)/cosd})

	# then go back to original scale for fold change
	seData<-2^seData
	control<-2^control
	case<-2^case	
	comean<-apply(control, 1, zmeasure, na.rm=TRUE)	
	curFC<-apply(seData, 2, function(x) {x/comean})

	# now calculate delta on the original scale also		
	curDelta<-apply(seData, 2, function(x) {x-comean})	

	# calculate groupFC and groupDiff from all probes that pass PALX
	caseO<-2^caseO
	controlO<-2^controlO
	controlMean<-apply(controlO, 1, zmeasure, na.rm=TRUE)
	caseMean<-apply(caseO, 1, zmeasure, na.rm=TRUE)
	
	groupFC<-caseMean / controlMean
	groupDiff<-caseMean - controlMean 
	
	return(list(zscore=curZscore, fc=curFC, delta=curDelta, tstats=tstats, origTstats=origTstats, groupFC=groupFC, groupDiff=groupDiff, 
			adjustPvals=adjustPvals, origPvals=origPvals))
}

#####
# calculate the positive and negative percentages from delta and score (zscore/fold change)
# also calculate the median score (zscore or fold change) for each module
#####
calculatePosNeg<-function(delta, zscore, fc, tstats, deltaThresh, zscoreThresh, fcThresh, subexpData, moduleAssignData, origTstats, groupFC, parameters)
{
	# apply thresholds to score and delta to count the positive and negative transcripts
	posMatZscore<-((zscore > zscoreThresh) & (abs(delta) > deltaThresh))
	negMatZscore<-((zscore < -zscoreThresh) & (abs(delta) > deltaThresh))
	posMatFC<-((fc > fcThresh) & (abs(delta) > deltaThresh))
	negMatFC<-((fc < 1/fcThresh) & (abs(delta) > deltaThresh))

	curFCThresh<-returnParameter(parameters, "fold_cut", "numeric")
	allPosVecFC<-(groupFC > fcThresh)
	allNegVecFC<-(groupFC < 1/fcThresh)
	
	posMatZscore<-replace(posMatZscore, which(posMatZscore==TRUE), 1)
	posMatZscore<-replace(posMatZscore, which(posMatZscore==FALSE), 0)
	posMatFC<-replace(posMatFC, which(posMatFC==TRUE), 1)
	posMatFC<-replace(posMatFC, which(posMatFC==FALSE), 0)
	
	negMatZscore<-replace(negMatZscore, which(negMatZscore==TRUE), 1)
	negMatZscore<-replace(negMatZscore, which(negMatZscore==FALSE), 0)	
	negMatFC<-replace(negMatFC, which(negMatFC==TRUE), 1)
	negMatFC<-replace(negMatFC, which(negMatFC==FALSE), 0)
	
	# combine counts into modules
	posCountZscore<-apply(posMatZscore, 2, function(x) {tapply(x, subexpData$Module, sum, na.rm=TRUE)})
	negCountZscore<-apply(negMatZscore, 2, function(x) {tapply(x, subexpData$Module, sum, na.rm=TRUE)})
	posCountFC<-apply(posMatFC, 2, function(x) {tapply(x, subexpData$Module, sum, na.rm=TRUE)})
	negCountFC<-apply(negMatFC, 2, function(x) {tapply(x, subexpData$Module, sum, na.rm=TRUE)})
	# combine all counts into modules
	allPosCountFC<-apply(as.matrix(allPosVecFC), 2, function(x) {tapply(x, subexpData$Module, sum, na.rm=TRUE)})
	colnames(allPosCountFC)<-"allFCPos"
	allNegCountFC<-apply(as.matrix(allNegVecFC), 2, function(x) {tapply(x, subexpData$Module, sum, na.rm=TRUE)})
	colnames(allNegCountFC)<-"allFCNeg"	
	
	tPosCount<-tapply(tstats, subexpData$Module, function(x){sum(x>0, na.rm=TRUE)})
	tNegCount<-tapply(tstats, subexpData$Module, function(x){sum(x<0, na.rm=TRUE)})
	# now combine into one tCount
	tBoth<-cbind(tPosCount, -tNegCount)
	colnames(tBoth)<-c("posPercent","negPercent")
	# take the difference rather than the maximum
	tCount<-tPosCount - tNegCount
	origTsum<-tapply(origTstats, subexpData$Module, sum, na.rm=TRUE)
	
	# now get the module counts
	# get the module counts from the probe assignment dataset
	# if used the expression data, could possibly miss some probes that were removed
	modCount<-table(moduleAssignData$module)
	
	if (all(modCount>0))
	{
	 	# calculte the positive and negative percentages
		posPercentZscore<-apply(posCountZscore, 2, function(x) {x/modCount})*100
		negPercentZscore<-apply(negCountZscore, 2, function(x) {x/modCount})*100
		posPercentFC<-apply(posCountFC, 2, function(x) {x/modCount})*100
		negPercentFC<-apply(negCountFC, 2, function(x) {x/modCount})*100
		allPosPercentFC<-apply(allPosCountFC, 2, function(x) {x/modCount})*100
		allNegPercentFC<-apply(allNegCountFC, 2, function(x) {x/modCount})*100
		tPercent<-tCount/modCount*100
		tBothPercent<-apply(tBoth, 2, function(x) {x/modCount})*100
		origTmean<-origTsum/modCount
		
		# now bind t sum with modules to find top modules
		topMod<-cbind(names(modCount), modCount, origTsum, origTmean, match(1:length(origTmean), order(abs(origTmean),decreasing=T)), allPosPercentFC, allNegPercentFC, allPosPercentFC-allNegPercentFC)	
		colnames(topMod)<-c("Module", "moduleCount", "tStatisticSum", "tStatisticMean", "Rank", "allFCPercentPos", "allFCPercentNeg", "allFCPercentDiff")
		topMod<-as.data.frame(topMod)
		topMod$Module<-as.character(topMod$Module)
		for (j in 2:ncol(topMod))
		{
			topMod[,j]<-as.numeric(as.character(topMod[,j]))
		}
	}
	
	# also calculate the median zscore & fold change for each module
	medianZscore<-apply(zscore, 2, function(x) {tapply(x, subexpData$Module, median, na.rm=TRUE)})
	medianFC<-apply(fc, 2, function(x) {tapply(x, subexpData$Module, median, na.rm=TRUE)})
	
	return(list(posPercentZscore=posPercentZscore, negPercentZscore=negPercentZscore, medianZscore=medianZscore,
				posPercentFC=posPercentFC, negPercentFC=negPercentFC, medianFC=medianFC,
	            posCountZscore=posCountZscore, negCountZscore=negCountZscore, posCountFC=posCountFC,
	            negCountFC=negCountFC, topMod=topMod, tPercent=tPercent, 
	            tBothPercent=tBothPercent, origTmean=origTmean))
}

#####
# calculate metrics to combine the positive and negative percentages into one value (and apply thresholds to some)
# e.g. difference, valence, and maximum
#####
calculatePercentMetric<-function(posPercentZscore, negPercentZscore, posPercentFC, negPercentFC, parameters)
{
	# maximum
	nnegPercentZscore<- -negPercentZscore
	nnegPercentFC<- -negPercentFC
	
	maxPercentZscore<-posPercentZscore
	maxPercentZscore[which(negPercentZscore > posPercentZscore)]<-nnegPercentZscore[which(negPercentZscore > posPercentZscore)]
	maxPercentFC<-posPercentFC
	maxPercentFC[which(negPercentFC > posPercentFC)]<-nnegPercentFC[which(negPercentFC > posPercentFC)]
		
	# difference 
	diffPercentZscore<- posPercentZscore - negPercentZscore
	diffPercentFC<- posPercentFC - negPercentFC
		
	# valence
	valenceZscore<-sqrt(abs(posPercentZscore^2-negPercentZscore^2))
	valenceZscore[which(negPercentZscore > posPercentZscore)]<- -valenceZscore[which(negPercentZscore > posPercentZscore)]
	valenceFC<-sqrt(abs(posPercentFC^2-negPercentFC^2))
	valenceFC[which(negPercentFC > posPercentFC)]<- -valenceFC[which(negPercentFC > posPercentFC)]

	return(list(maxPercentZscore=maxPercentZscore, diffPercentZscore=diffPercentZscore, valenceZscore=valenceZscore,
				maxPercentFC=maxPercentFC, diffPercentFC=diffPercentFC, valenceFC=valenceFC))
}

######
# calculate summaries over the modules within a sample
######
sampleSummary<-function(score, posCount, negCount, parameters)
{
	if (toupper(returnParameter(parameters, "delta_type", "character"))=="FOLD")
	{
#		samScore<-apply(log2(score), 2, sum, na.rm=TRUE)
		if (any(score==0))
		{
			replaceVal<-min(score[which(score!=0)])
		}
		samScore<-apply(score, 2, function(x)
		{
			if (any(x==0))
			{
				x[which(x==0)]<-replaceVal
			}
			sum(log2(x), na.rm=T)
		})
	}
	else
		samScore<-apply(score, 2, sum, na.rm=TRUE)
		
	samPosCount<-apply(posCount, 2, sum, na.rm=TRUE)
	samNegCount<-apply(negCount, 2, sum, na.rm=TRUE)
	samTotalCount<-samPosCount+samNegCount
	
	return(list(samScore=samScore, samPosCount=samPosCount, samNegCount=samNegCount, 
				samTotalCount=samTotalCount))
}


#####
# write multiple output tables
#####
writeOutputFiles<-function(ccc, samsum, summaryStats, path_results, project_name, designData, parameters)
{
	platform<-returnParameter(parameters, "platform_type", type="character")
	if (length(platform)==0)
		platform<-"Microarray"

	# output the probe by individual (sample) delta (difference), zscore, and fold change
	outputOneFile(curOut=ccc$delta, path_results, project_name, curTitle="Individual Delta ", designData)
	outputOneFile(curOut=ccc$zscore, path_results, project_name, curTitle="Individual Zscore ", designData)
	outputOneFile(curOut=ccc$fc, path_results, project_name, curTitle="Individual FC ", designData)
		
	# difference in positive and negative percentages by sample and module
	outputOneFile(curOut=ccc$diffPercentZscore, path_results, project_name, 
		curTitle="Percent Difference Zscore ", designData)
	outputOneFile(curOut=ccc$diffPercentFC, path_results, project_name, 
		curTitle="Percent Difference FC ", designData)
		
	# positive percentage by sample and module
	outputOneFile(curOut=ccc$posPercentZscore, path_results, project_name,
		curTitle="Positive Percent Zscore ", designData)		
	outputOneFile(curOut=ccc$posPercentFC, path_results, project_name,
		curTitle="Positive Percent FC ", designData)
	# negative percentage by sample and module
	outputOneFile(curOut=ccc$negPercentZscore, path_results, project_name,
		curTitle="Negative Percent Zscore ", designData)
	outputOneFile(curOut=ccc$negPercentFC, path_results, project_name,
		curTitle="Negative Percent FC ", designData)

	# the median zscore or fold change by sample and module
	outputOneFile(curOut=ccc$medianZscore, path_results, project_name,
		curTitle="Median Zscore ", designData)
	outputOneFile(curOut=ccc$medianFC, path_results, project_name,
		curTitle="Median FC ", designData)
				
	# next output summary statistics
	outputOneFile(curOut=summaryStats$moduleMean, path_results, project_name,
		curTitle="Module Mean ", designData)

	outputOneFile(curOut=summaryStats$moduleSd, path_results, project_name,
		curTitle="Module Standard Deviation ", designData)
	
	outputOneFile(curOut=summaryStats$moduleMedian, path_results, project_name,
		curTitle="Module Median ", designData)

	outputOneFile(curOut=summaryStats$moduleMad, path_results, project_name,
		curTitle="Module Median Absolute Deviation ", designData)
		
	# write DEGs out slightly differently because it is a list
	for (i in 1:length(ccc$DEGs))
	{
		if (length(ccc$DEGs)==1)
		{
			fileName<-paste(path_results, project_name, "_", "Probe Level Statistics ",
			    format(Sys.time(), "%a %b %d %H%M%S %Y"), ".csv", sep="")
		}
		else
		{
			curSubset<-sort(unique(designData$Subset))[i]
			fileName<-paste(path_results, project_name, "_", "Probe Level Statistics Subset ",
				curSubset, "_", format(Sys.time(), "%a %b %d %H%M%S %Y"), ".csv", sep="")
		}
		# replace spaces with underscores in the file name
		fileName<-gsub(" ", "_", fileName)
		dataToWrite<-as.matrix(ccc$DEGs[[i]])
		write.csv(dataToWrite, fileName, row.names=FALSE)
	}
	
	# also write out the top modules
	for (i in 1:length(ccc$moduleOrder))
	{
		if (length(ccc$moduleOrder)==1)
		{
			fileName<-paste(path_results, project_name, "_", "Module Level Statistics ",
			    format(Sys.time(), "%a %b %d %H%M%S %Y"), ".csv", sep="")
		}
		else
		{
			curSubset<-sort(unique(designData$Subset))[i]
			fileName<-paste(path_results, project_name, "_", "Module Level Statistics Subset ",
				curSubset, "_", format(Sys.time(), "%a %b %d %H%M%S %Y"), ".csv", sep="")
		}
		# replace spaces with underscores in the file name
		fileName<-gsub(" ", "_", fileName)
		dataToWrite<-as.matrix(ccc$moduleOrder[[i]])
		write.csv(dataToWrite, fileName, row.names=FALSE)
	}
	
	# write out module group comparison for differential and GSA
	for (i in 1:length(ccc$tPercentList))
	{
		# write the results to a csv file
		if (length(ccc$tPercentList)==1)
			fileName<-paste(path_results, project_name, "_", "260_Module_Group_Comparison_",
			    format(Sys.time(), "%a %b %d %H%M%S %Y"), ".csv", sep="")
		else
		{
			curSubset<-sort(unique(designData$Subset))[i]
			fileName<-paste(path_results, project_name, "_", "260_Module_Group_Comparison_Subset_",
				curSubset, "_",
			    format(Sys.time(), "%a %b %d %H%M%S %Y"), ".csv", sep="")		}	
		# replace spaces with underscores in the file name
		fileName<-gsub(" ", "_", fileName)
		dataToWrite<-ccc$tBothPercentList[[i]]
		colnames(dataToWrite)<-c("positivePercent", "negativePercent")
		write.csv(dataToWrite, fileName)
		
		if (toupper(platform)==toupper("Microarray"))
		{
			# also write out the results for GSAList
			if (length(ccc$GSAList)==1)
				fileName<-paste(path_results, project_name, "_", "Gene_Set_Analysis_260_Module_Group_Comparison_",
				    format(Sys.time(), "%a %b %d %H%M%S %Y"), ".csv", sep="")
			else
			{
				curSubset<-sort(unique(designData$Subset))[i]
				fileName<-paste(path_results, project_name, "_", "Gene_Set_Analysis_260_Module_Group_Comparison_Subset_",
					curSubset, "_",
			    	format(Sys.time(), "%a %b %d %H%M%S %Y"), ".csv", sep="")		}	
			# replace spaces with underscores in the file name
			fileName<-gsub(" ", "_", fileName)
			dataToWrite<-as.matrix(ccc$GSAList[[i]][,c(1,4)])
			colnames(dataToWrite)[2]<-"pvalue"
			write.csv(dataToWrite, fileName)		
		}
	}
		
	# next need to output sample summaries - bind together the 4 list elements
	curOut<-rbind(samsum$samScore, samsum$samPosCount, samsum$samNegCount, samsum$samTotalCount)
	rownames(curOut)<-c("SampleScore","SamplePositiveCount","SampleNegativeCount","SampleTotalCount")
	outputOneFile(curOut=curOut, path_results, project_name,
		curTitle="Sample Summaries ", designData)
}

#####
# write one output file
#####
outputOneFile<-function(curOut, path_results, project_name, curTitle, designData)
{
	fName<-paste(path_results, project_name, "_", curTitle,
	    format(Sys.time(), "%a %b %d %H%M%S %Y"), ".csv", sep="")
	# replace any spaces with underscores in the file name
	fName<-gsub(" ", "_", fName)
	# change column names from unique id to sample id
	if (colnames(curOut)[1]=="All")
	{
		newCNames<-c("All", designData$sample_id[match(colnames(curOut)[2:ncol(curOut)], designData$uniqueID)])
		if (length(newCNames)==length(unique(newCNames)))
			colnames(curOut)<-newCNames	
	}
	else
	{
		newCNames<-designData$sample_id[match(colnames(curOut), designData$uniqueID)]
		if (length(newCNames)==length(unique(newCNames)))
			colnames(curOut)<-newCNames
	}
	# this will write out both the row names and column names
	write.csv(x=curOut, file=fName)
}


#######
# same as p.adjust (only uses the BH or fdr method)
# only difference is how it treats missing values (this function is the same as how SAS calculates FDR)
# include the full length of the vector (including missing values) when adjusting the p-value
# this change makes the algorithm more lenient (more probes will pass after adjusting)
#######
my.p.adjust<-function(p)
{
	p0 <- p
    if (all(nna <- !is.na(p))) 
        nna <- TRUE
    p <- as.vector(p[nna])
    n<-length(p0)
    stopVal<-sum(is.na(p0))+1
    p0[nna] <-
	{
		i <- n:stopVal
        o <- order(p, decreasing = TRUE)
        ro <- order(o)
        pmin(1, cummin(n/i * p[o]))[ro]
	}
	return(p0)
}


