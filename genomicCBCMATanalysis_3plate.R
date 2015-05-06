# read in input files
currentDirectory<-"/Users/ewhalen/Projects/genomicCBC/MATanalysis/P14_U01Mexico_3Plates_Phase2_ShigellaVsHC/"

dataDirectory<-paste(currentDirectory, "Data/", sep="")
resultsDirectory<-paste(currentDirectory, "Results/", sep="")
scriptsDirectory<-paste(currentDirectory, "scripts/", sep="")

source(paste(scriptsDirectory, "genomicCBCMATfunctions_3plate.R", sep=""))

# load required packages
packageLoad(c("limma","grid","lattice","latticeExtra","lme4"))

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

# add plate 3 information
geneLevel$plate3<-"No"
for (i in 1:nrow(geneLevel))
{
	curID<-as.character(rownames(geneLevel)[i])
	geneLevel$plate3[i]<-as.character(moduleAssignData$plate3[which(curID==as.character(moduleAssignData$geneID))])
}
write.csv(geneLevel, file=paste(resultsDirectory, "geneLevelStatistics.csv", sep=""), row.names=FALSE)

# now calculate an average FC for each module
modLevelFCavg<-CaseVsControlModuleLevelFC(FCmoduleLevelControls, designData)

orderIndex<-order(abs(log2(modLevelFCavg)), decreasing=TRUE)
modLevelFCavgOrd<-as.data.frame(modLevelFCavg[orderIndex])
colnames(modLevelFCavgOrd)<-"FC"
modLevelFCavgOrd$plate3<-"No"
modLevelFCavgOrd$numSigGenes<-rep(0, nrow(modLevelFCavgOrd))
for (i in 1:nrow(modLevelFCavgOrd))
{
	curMod<-rownames(modLevelFCavgOrd)[i]
	modLevelFCavgOrd$plate3[i]<-as.character(moduleAssignData$plate3[which(moduleAssignData$module==curMod)[1]])
	modLevelFCavgOrd$numSigGenes[i]<-sum(geneLevel$adj.P.Val[which(geneLevel$module==curMod)] < 0.05, na.rm=T)
}

write.csv(modLevelFCavgOrd, file=paste(resultsDirectory, "groupModuleLevelFC.csv", sep=""))

# fit module level model with random effects for individual and gene
modLevelResults<-moduleLevelModelWithRandomEffect(FCcomparedToControls, designData, moduleAssignData)

write.csv(modLevelResults, file=paste(resultsDirectory, "moduleLevelModelFit.csv", sep=""))

##### Plots

# look at correlation of genes within a module
modGeneCorrelation(FCvals=subExpData, moduleAssignData, designData, fileName=paste(resultsDirectory, "moduleGeneCorrelation.pdf",sep=""))

# stat sig by module pie chart
modLevelPercents<-reorderModules(statSigCountByModule(geneLevelStats=geneLevel), moduleAssignData)
colnames(modLevelPercents)<-c("PercentUp", "PercentDown")

# make plot
make6x11gridPercent(modLevelPercents, file=paste(resultsDirectory, "moduleGridPercents.png", sep=""))

# now make grid key
newModFC<-reorderModules(modFC=modLevelFCavgOrd, moduleAssignData)
make6x11grid(modFC=newModFC, fileName=paste(resultsDirectory,"moduleGrid.png",sep=""), maxVal=2, sigCutoff=0.001, modLevRes=modLevelResults)
make6x11gridkey(modFC=newModFC, fileName=paste(resultsDirectory,"moduleGridKey.png",sep=""))

# 6/19/14 add extra column for significance for module level model fit
modLevelResults$sig<-0
modLevelResults$sig[which(modLevelResults$fdr < 0.05 & modLevelResults$groupParameter > 0)]<-1
modLevelResults$sig[which(modLevelResults$fdr < 0.05 & modLevelResults$groupParameter < 0)]<- -1

newModLevelResults<-reorderModules(modFC=modLevelResults, moduleAssignData)

make6x11gridModuleLevel(newModLevelResults, fileName=paste(resultsDirectory,"moduleLevelModelGrid.png",sep=""), maxVal=1)


# now make heatmap for individual values
# need to make sample info matrix
designSub<-designData[match(colnames(FCmoduleLevelControls), designData$sampleID),]
# check that they match
all(colnames(FCmoduleLevelControls)==designSub$sampleID)

# change the column names to sample names
FCforHeatmap<-FCmoduleLevelControls
colnames(FCforHeatmap)<-designSub$sampleName

# for generic code, just use groupName in samInfo
keepIndex<-which(colnames(designSub)=="groupName")
samInfo<-as.data.frame(designSub[,keepIndex])

makeClusteredHeatmapWithSampleInfo(t(log2(FCforHeatmap)), fileName=paste(resultsDirectory, "individualModuleLevel_WithSampleInfo.png", sep=""), curBreaks=3, clusterSams=TRUE, clusterMods=TRUE, sampleInfo=samInfo, allDesign=designSub, includeControl=TRUE, modLevelResults=modLevelResults, fdr=1)

makeClusteredHeatmapWithSampleInfo(t(log2(FCforHeatmap)), fileName=paste(resultsDirectory, "individualModuleLevel_WithSampleInfo_noControls.png", sep=""), curBreaks=3, clusterSams=TRUE, clusterMods=TRUE, sampleInfo=samInfo, allDesign=designSub, includeControl=FALSE, modLevelResults=modLevelResults, fdr=1)

makeClusteredHeatmapWithSampleInfo(t(log2(FCforHeatmap)), fileName=paste(resultsDirectory, "individualModuleLevel_WithSampleInfo_sigMods.png", sep=""), curBreaks=3, clusterSams=TRUE, clusterMods=TRUE, sampleInfo=samInfo, allDesign=designSub, includeControl=TRUE, modLevelResults=modLevelResults, fdr=0.01)

# add a PCA plot 
makePCAplots(FCcomparedToControls, designData, resultsDirectory)

