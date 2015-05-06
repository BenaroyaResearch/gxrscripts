### Preprocessing steps:
# 1. assume input is background subtracted data
# 2. quantile normalization
# 3. floor to 10
# 4. log2
# 5. PALO

## ASSUMPTION: array outliers have been removed before normalization!!

##inputs:
# expDataFile - the background subtracted data file name (either txt or csv file) or the actual data set
#		if this is a character string, then will try to load file; otherwise expect a dataframe
# outputFile - where to place the normalized data
#		if this is NA, then it will just return the data set rather than writing it to file
# floorValue - the value to floor to (default is 10)
# PALO - a boolean (T or F) - if F, then assume performing PALX
# PALX - if PALO=F, then this is a percentage of the samples that must have p-value < 0.01 (a value between 0 and 1)
#		as an example: if PALX is 0.5 and there are 37 arrays, then 19 arrays must have a p-value < 0.01 (so more stringent
#		than expecting 18 arrays to pass) - basically, I take PALX*num_columns and then round up to the integer
# allCols - a boolean (T or F) to indicate whether all columns should be included in the output
#		if T, then include p-value columns; if F, p-values are not included in the output
# ComBat - a boolean (T or F) to indicate whether ComBat code should be called to correct batch effects
# ComBatFile - the file location of the ComBat.R code (will only be used if ComBat=T)
# sampleInfoFile - the sample information file that will be used by ComBat (will only be used if ComBat=T)
#		it is assumed that the user knows how to set up the sample information file correctly (see ComBat), must be .txt file
# writeNormalizeFile - the file location to write the normalized file (will only be used if ComBat=T), must be .txt file (this is temporary file - will be deleted)
# par.prior - a boolean (T or F) for the ComBat function to indicate whether to have a parametric prior
# added parameters: performNormalize, performFloor, performLog2, performPALO: all are booleans (T or F) to indicate whether each 
#  step in the preprocessing algorithm should be followed (default is TRUE to follow all steps)
# sigColName - the string found in the column names that represents the column as signal intensity (and not detection p-value)
# pvalCol - the string found in the column names that represents the columns as detection p-values (not signal intensity)
# rangeCutoff - the minimum range a probe must show to pass the filter (to filter probes based on variability)
preprocess<-function(expDataFile, outputFile, floorValue=10, PALO=FALSE, PALX=0.1, allCols=TRUE, ComBat=FALSE, ComBatFile, sampleInfoFile,
	writeNormalizeFile, par.prior=TRUE, performNormalize=TRUE, performFloor=TRUE, performLog2=TRUE, performPALO=TRUE, sigColName="Signal", 	pvalCol="Pval", rangeCutoff=0)
{
	# require preprocessCore package to run quantile normalization
	packageBioLoad("preprocessCore")
	
	if (is.character(expDataFile))
	{
		fileType<-tolower(strsplit(expDataFile, "\\.")[[1]][2])
		if (fileType=="txt")
			expData<-read.delim(expDataFile)
		else if (fileType=="csv")
			expData<-read.csv(expDataFile)
	} else {expData<-expDataFile}

	# run quantile normalization - this returns intensity columns only 
	if (performNormalize==TRUE)
		subexpData<-normalization(expData, sigColName)
	else
		subexpData<-expData[,grep(toupper(sigColName), toupper(colnames(expData)))]
	
	if (ComBat==TRUE)
	{
		# need to add the probe id column before calling ComBat
		temp<-cbind(expData[,1], subexpData)
		colnames(temp)[1]<-"PROBE_ID"
		write.table(temp, writeNormalizeFile, sep="\t", row.names=F)
		source(file=ComBatFile)

		# next perform ComBat
		tmp<-ComBat(expression_xls=writeNormalizeFile, sample_info_file=sampleInfoFile, write=F, par.prior=par.prior, skip=1)
		
		# if tmp is a character string, then an error occurred
		if (!is.character(tmp))
		{
			# remove the probe id column
			subexpData<-tmp[,2:ncol(tmp)]
		}
		else
			print(tmp)
			
		# want to remove the file that was created with the normalized data
		unlink(writeNormalizeFile)
	}
	
	# then perform flooring
	if (performFloor==TRUE)
		subexpData<-replace(subexpData, subexpData < floorValue, floorValue)
	
	# then take log2
	if (performLog2==TRUE)
		subexpData<-log2(subexpData)
	
	# finally perform PALO
	if (PALO==TRUE)
		percentage<-NA
	else
	{
		if (PALX < 0 | PALX > 1)
			print("ERROR: PALX must fall between 0 and 1.")
		else
			percentage<-PALX
	}
	if (performPALO==TRUE)
		keepRows<-performPALOfun(expData, pvalCut=0.01, percentage=percentage, pvalCol)
	else
		keepRows<-rep(TRUE, nrow(expData))
		
	if (any(is.na(keepRows)))
	{
		keepRows[which(is.na(keepRows))]<-FALSE
	}
	
	# remove any probes that do not pass range filter
	if (rangeCutoff >= 0)	
	{
		if (performLog2==TRUE)
		{
			dataRanges<-apply(2^subexpData, 1, function(x)
			{
				if (all(is.na(x)))
					return(NA)
				else
					return(diff(range(x, na.rm=TRUE)))
			})
		}
		else
		{
			# if there are no large values, then assume it's been log2 transformed already
			if (max(subexpData, na.rm=TRUE) < 50)
			{
				dataRanges<-apply(2^subexpData, 1, function(x)
				{
					if (all(is.na(x)))
						return(NA)
					else
						return(diff(range(x, na.rm=TRUE)))
				})
			}
			else
			{
				dataRanges<-apply(subexpData, 1, function(x)
				{
					if (all(is.na(x)))
						return(NA)
					else
						return(diff(range(x, na.rm=TRUE)))
				})
			}
		}
		keepRowsByRange<-(dataRanges > rangeCutoff)
		# NA for range indicates that all the values for that row are missing
		if (any(is.na(keepRowsByRange)))
			keepRowsByRange[which(is.na(keepRowsByRange))]<-FALSE
		keepRows<-(keepRows & keepRowsByRange)
	}
		
	# now get data ready to write to file
	if (allCols==TRUE)
	{
		newExpData<-expData
		changeCols<-grep(toupper(sigColName), toupper(colnames(expData)))
		newExpData[,changeCols]<-as.matrix(subexpData)
		# now subset to those that pass PALO
		finalData<-newExpData[keepRows,]
	}
	else
	{
		subexpData<-subexpData[keepRows,]
		probeCol<-grep("PROBE_ID", toupper(colnames(expData)))
		subexpData<-cbind(expData[keepRows,probeCol], subexpData)
		finalData<-subexpData
		colnames(finalData)[1]<-"PROBE_ID"
	}
	
	if (!is.na(outputFile))
	{
		fileType<-tolower(strsplit(outputFile, "\\.")[[1]][2])
		if (fileType=="txt")
			write.table(finalData, file=outputFile, sep="\t", row.names=F)
		else if (fileType=="csv")
			write.csv(finalData, file=outputFile, row.names=F)
	} else (return(finalData))
}


normalization<-function(expData, sigColName)
{
	# assume that the column names with intensity values contain "AVG" or whatever the signal column name string is
	keepIndex<-grep(toupper(sigColName), toupper(colnames(expData)))
	x.matrix<-expData[,keepIndex]
	
	if (class(x.matrix)!="matrix")
		x.matrix<-as.matrix(x.matrix)
	
	# do quantile normalization - will probably have some negative values because background corrected
	xPre<-as.data.frame(normalize.quantiles(x.matrix))
	colnames(xPre)<-colnames(x.matrix)

	return(xPre)
}


performPALOfun<-function(expData, pvalCut=0.01, percentage=NA, pvalCol)
{
	# subset to just p-value data
	keepIndex<-grep(toupper(pvalCol), toupper(colnames(expData)))
	useData<-expData[,keepIndex]
	
	if(is.na(percentage))
	{
		# then assume only need 1 p-value to pass the filter
		percentage<-1/ncol(useData)
	}
	keepRowIndex<-apply(useData, 1, function(x)
	{	
		# if all values are missing (NA), then do not filter out that row
		if (all(is.na(x)))
			return(TRUE)
		else
		{
			# if some of the pvalues are NA, then set the pvalues temporarily to zero as though they should pass the filter (be on the safe side)
			if (any(is.na(x)))
				x[-which(is.na(x))]<-0
			return((sum(x < pvalCut) >= ceiling(percentage*ncol(useData))))
		}	
	})
	return(keepRowIndex)
}

packageBioLoad = function(x) {
	# returns a matrix of the installed packages
	inst = installed.packages()
	matchVals = match(x, inst[,"Package"])
	missing = x[is.na(matchVals)]
	if (length(missing)>0) 
	{
		cat("Installing packages...", "\n")
		#lapply(missing, install.packages, repos="http://www.revolution-computing.com/cran/")
		source("http://bioconductor.org/biocLite.R")
		biocLite("preprocessCore")
	}
	# then load libraries
	for (i in 1:length(x))
	{
		require(x[i], character.only=TRUE)
	}
}

