params <- commandArgs(trailingOnly = T)
path<- params[1]
packageLoad = function(x) {

inst = installed.packages()
        matchVals = match(x, inst[,"Package"])
        missing = x[is.na(matchVals)]
        if (length(missing)>0)
        {
                cat("Installing packages...", "\n")
                lapply(missing, install.packages, repos="http://www.revolution-computing.com/cran/")
        }
        # then load libraries
        for (i in 1:length(x))
        {
                require(x[i], character.only=TRUE)
        }


 if ("limma" %in% missing)
        {
                source("http://bioconductor.org/biocLite.R")
                biocLite("limma")
                require(limma)
        }
}

packageLoad("limma")


# load the 3 input files
expData<-read.delim(paste(path,"Data/expression.tsv", sep=""))

designData<-read.csv(paste(path,"Data/design.csv", sep=""))
comparisonData<-read.csv(paste(path,"Data/comparison.csv",sep=""))

# load the preprocessing file
source(paste(path, "scripts/Preprocessing.R", sep=""))
source(paste(path, "scripts/allRanking.R", sep=""))

caseNames<-comparisonData$cases
controlNames<-comparisonData$controls

retValue<-performAllRanking(curData=expData, designData=designData, caseNames=caseNames, controlNames=controlNames, toPlot=FALSE, fileDirectory=paste(path,"Results/",sep=""))
