## source functions for analysis
#source("/Users/ewhalen/Rscripts/ModAnFun.R")
#source("/Users/ewhalen/Rscripts/Preprocessing.R")
## path for file locations
#path<-"/Users/ewhalen/Software/ModuleMap/ModuleAnalysis/"

params <- commandArgs(trailingOnly = T)
workdir <- params[1]
scriptdir <- paste(workdir, "/scripts", sep="/")

# source functions for analysis
source(paste(scriptdir, "ModAnFun.R", sep="/"))
source(paste(scriptdir, "Preprocessing.R", sep="/"))
path<-workdir


# for now read parameters from the parameters file
parameters<-readParameters(path)

if (toupper(returnParameter(parameters, "platform_type", "character"))==toupper("RNA-seq"))
{
	parameters$Value<-as.character(parameters$Value)
	parameters$Value[which(parameters$Parameter=="platform_type")]<-"RNAseq"
}

# load libraries - check that all of these libraries are needed
packageLoad(c("lattice", "ggplot2", "grid", "preprocessCore", "latticeExtra", "limma", "GSA", "impute","edgeR"))
 
# create results directory, if necessary
project_name <- returnParameter(parameters, key="project", type="character")
path_results<-createResults(path, project_name)

# read in data
expDataFile<-returnParameter(parameters, "expression", "character")
designDataFile<-returnParameter(parameters, "design", "character")
moduleAssignFile<-returnParameter(parameters, "module_version", "character")
moduleAnnotateFile<-returnParameter(parameters, "module_function", "character")

expData<-readData(path, project_name, expDataFile)
designData<-readData(path, project_name, designDataFile)
moduleAssignData<-readData(path, "ModuleData", moduleAssignFile)
# there are duplicates in the module file so remove them
moduleAssignData<-unique(moduleAssignData)
moduleAnnotateData<-readData(path, "ModuleData", moduleAnnotateFile)
# read in gene symbol information
geneSymbolData<-readData(path, "ModuleData", "ProbeToSymbol.csv")
geneSymbolData<-unique(geneSymbolData)
# reorder module annotations based on function (hopefully makes the ordering more meaningful)
moduleAnnotateData<-moduleAnnotateData[order(moduleAnnotateData$Function),]

continue<-TRUE

designPassCheck<-designDataCheck(designData, parameters)
expPassCheck<-expDataCheck(expData, parameters)
designData<-designPassCheck$designData

if (designPassCheck$passValue==FALSE | expPassCheck==FALSE)
	continue<-FALSE

# write out information to log file to help with error checking
print(paste("The design file has ", nrow(designData), " samples.", sep=""))
print(paste("The dimensions of the expression data are ", nrow(expData), " rows and ", ncol(expData), " columns.", sep=""))

if (continue)
{
	# match the sample names to the column names in the signal intensity file
	if (toupper(trimString(returnParameter(parameters, "input", "character")))==toupper("BeadStudio"))
		designData$uniqueID<-paste("X", designData$barcode, "_", designData$position, ".", 
			returnParameter(parameters, key="signal_pattern", type="character"), sep="")
	else
	{
		# check if the first character in barcode is numeric (R will add an X if this is true)
		firstChar<-unique(substr(designData$barcode, 1, 1))
		if (all(firstChar %in% as.character(0:9)))
			designData$uniqueID<-paste("X", designData$barcode, ".", 
				returnParameter(parameters, key="signal_pattern", type="character"), sep="")
		else
			designData$uniqueID<-paste(designData$barcode, ".", 
				returnParameter(parameters, key="signal_pattern", type="character"), sep="")
	}
	
	# sample_id is used in all of the plots to label the samples
	if (!("sample_id" %in% colnames(designData)))
		designData$sample_id<-designData$barcode

	# synchronize the experimental, design, and module data
	# preprocessing will most likely result in less probes that are greater than the threshold
	retSubData<-syncData(expData, designData, moduleAssignData, parameters, rearrange=FALSE)
	subexpData<-retSubData$subexpData
	origSubExpData<-retSubData$origSubExpData
	if (nrow(subexpData)==0 | nrow(origSubExpData)==0)
		continue<-FALSE
	# no longer need the original expression data
	rm(expData)
	
	# write out information to log file to help with error checking
	print(paste("After preprocessing, the dimensions of the expression data are ", nrow(origSubExpData), " rows and ", ncol(origSubExpData), " columns.", sep=""))
	print(paste("After preprocessing and matching to probes in modules, the dimensions of the expression data are ", nrow(subexpData), " rows and ", ncol(subexpData), " columns.", sep=""))	
	
	if (continue)
	{
		# test if any modules are missing from the expression data 
		moduleAssignData<-updateModuleAssignData(subexpData, moduleAssignData, parameters)
		moduleAnnotateData<-updateModuleAnnotateData(subexpData, moduleAnnotateData, parameters)

		# calculate summary statistics per sample and module
		summaryStats<-calcSummaryStats(subexpData, parameters, moduleAssignData)

		# next compare case to control - contains the 5 data sets used to make the chiaroscuro sketches (plus other data sets)
		ccc<-compareCaseControl(parameters, subexpData, designData, moduleAssignData, origSubExpData, geneSymbolData)
	
		print(paste("The probe level data has ", nrow(ccc$DEGs[[1]]), " rows. This should match the rows in the preprocessed expression data: ", nrow(origSubExpData), ".", sep=""))
	
		# make sample summaries (for barplots)
		if (returnParameter(parameters, "delta_type", "character")=="fold")
			samsum<-sampleSummary(ccc$fc, ccc$posCountFC, ccc$negCountFC, parameters)
		else if (returnParameter(parameters, "delta_type", "character")=="zscore")
			samsum<-sampleSummary(ccc$zscore, ccc$posCountZscore, ccc$negCountZscore, parameters)
		
		# write the csv files
		writeOutputFiles(ccc, samsum, summaryStats, path_results, project_name, designData, parameters)

		print("Run complete.")
	}
}
