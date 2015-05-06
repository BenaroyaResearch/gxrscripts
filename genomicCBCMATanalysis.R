# read in input files
currentDirectory<-"/Users/ewhalen/Projects/genomicCBC/MATanalysis/P14_U01Mexico/"


dataDirectory<-paste(currentDirectory, "Data/", sep="")
resultsDirectory<-paste(currentDirectory, "Results/", sep="")
scriptsDirectory<-paste(currentDirectory, "scripts/", sep="")

source(paste(scriptsDirectory, "genomicCBCMATfunctions.R", sep=""))

# load required packages
packageLoad("limma")

# read in data
expData<-readData(dataDirectory, "FC.csv")
designData<-readData(dataDirectory, "designFile.csv")
moduleAssignData<-readData(dataDirectory, "moduleAssignFile.csv")

# FC values have row names in the first column
rownames(expData)<-as.character(expData[,1])
expData<-expData[,-1]

subExpData<-subsetExpData(expData, designData, moduleAssignData)

# next calculate FC values compared to controls (rather than reference samples)
FCcomparedToControls<-FCtocontrols(subExpData, designData)

# these are the individual module level results
FCmoduleLevelControls<-convertFCtoModule(FCcomparedToControls, moduleAssignData)

# write to csv
write.csv(FCmoduleLevelControls, file=paste(resultsDirectory, "individualModuleLevelFC.csv", sep=""))

# next run gene level analysis using limma to find statistically significant genes
# use the original FC values (compared to reference sample) - get the exact same results if use data compared to controls
geneLevel<-fitLimmaToGenes(subExpData, designData, moduleAssignData)

write.csv(geneLevel, file=paste(resultsDirectory, "geneLevelStatistics.csv", sep=""), row.names=FALSE)

# now calculate an average FC for each module
modLevelFCavg<-CaseVsControlModuleLevelFC(FCmoduleLevelControls, designData)

modLevelFCavg<-as.data.frame(modLevelFCavg)
colnames(modLevelFCavg)<-"FC"

write.csv(modLevelFCavg, file=paste(resultsDirectory, "groupModuleLevelFC.csv", sep=""))


