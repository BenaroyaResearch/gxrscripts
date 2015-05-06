#####
# function to create the bar plots for the sample summaries
#####
makeBarPlot<-function(curData, designData, parameters, curTitle, path_results)
{
	# change the names on curData to be the sample IDs (rather than chip code)
	matchIndex<-match(toupper(names(curData)), toupper(designData$uniqueID))
	if (any(is.na(matchIndex)))
		matchIndex<-matchIndex[-which(is.na(matchIndex))]
	samID<-as.character(designData$sample_id[matchIndex])

#	Groups<-factor(designData$group[matchIndex])
#	# 0 refers to control and 1 refers to case
#	levels(Groups)<-c("Control", "Case")
	if ("group_label" %in% colnames(designData))
		Groups<-factor(designData$group_label[matchIndex])
	else
	{
		Groups<-factor(designData$group[matchIndex])
		# 0 refers to control and 1 refers to case
		levels(Groups)<-c("Control", "Case")
	}

	names(curData)<-samID

	# now make plot
	# use factor on the names to ensure that the order is not rearranged (otherwise, will do alphabetical order)
	curPlot<-qplot(factor(names(curData),levels=names(curData)), curData, geom="bar", stat="identity", fill=Groups) + opts(title=curTitle, axis.text.x=theme_text(angle=-90, hjust=0)) + xlab("") + ylab("") + scale_fill_brewer(palette="Set2")
	
	# write the plot to a file
	# turn the warnings off
	options(warn=-1)
	project_name<-returnParameter(parameters, "project", "character")
	curFileName<-paste(path_results,  project_name, "_", curTitle, "_", format(Sys.time(), "%a %b %d %H%M%S %Y"), ".png", sep="")
	# replace any spaces with underscore in the file name
	curFileName<-gsub(" ", "_", curFileName)
	png(filename=curFileName, width=10, height=8, units="in", res=200, pointsize=20)
	print(curPlot)
	dev.off()
	options(warn=0)
}

#####
# make heatmap plots of samples vs. modules
#####
makeHeatmap<-function(curData, parameters, designData, moduleAnnotateData, dataName, path_results, moduleOrder)
{
	# if the first column name is All, then remove for plotting
	if (colnames(curData)[1]=="All")
	{
		curData<-curData[,-1]
	}
	# make plots of all modules
	newData<-alterModule(curData, parameters, designData, moduleAnnotateData, type="All", NA)
	# need to get the full range of the zscores for coloring the heatmap
	if (toupper(dataName) == toupper("Zscore") | toupper(dataName) == toupper("Fold"))
	{
		if (toupper(dataName)==toupper("Zscore"))
			curBreaks<-max(abs(newData))
		else
		{
			# need to set the fold change to negative when it's less than one to show it's downregulated
			#less1<-which(newData<1)
			#newData[less1]<--(1/newData[less1])
			newData<-log2(newData)
			curBreaks<-max(abs(newData))
			# do not let the largest value get too large for fold change or it will remove all color from the plot
			if (curBreaks > log2(20))
				curBreaks<-log2(20)
			# then need to reset any data that are larger than this value
			newData<-apply(newData, 2, function(x)
			{
				if (any(x > curBreaks))
					x[which(x>curBreaks)]<-curBreaks
				if (any(x < -curBreaks))
					x[which(x< -curBreaks)]<- -curBreaks
				return(x)
			})
			dataName<-paste("Log2", dataName, sep="")
		}
	}
	else
		curBreaks<-101.01
	plotModuleSketch(newData, dataName, path_results, selection="All", parameters, curBreaks)
	# then make clustered heatmap from newData
	makeClusteredHeatmap(newData, parameters, designData, moduleAnnotateData, dataName, path_results, type="All", curBreaks)
	
	# make plots of annotated modules, if the annotated data exist
	if (nrow(moduleAnnotateData)>0)
	{
		newData<-alterModule(curData, parameters, designData, moduleAnnotateData, type="Annotated", NA)
		if (length(grep(toupper("Fold"), toupper(dataName)))>0)
		{
			newData<-log2(newData)
			newData<-apply(newData, 2, function(x)
			{
				if (any(x > curBreaks))
					x[which(x>curBreaks)]<-curBreaks
				if (any(x < -curBreaks))
					x[which(x< -curBreaks)]<- -curBreaks
				return(x)
			})	
		}
		plotModuleSketch(newData, dataName, path_results, selection="Annotated", parameters, curBreaks)
		makeClusteredHeatmap(newData, parameters, designData, moduleAnnotateData, dataName, path_results, type="Annotated", curBreaks)
	}
	
	# make plots of top changing modules (show most difference between case and control)
	newData<-alterModule(curData, parameters, designData, moduleAnnotateData, type="Changing", moduleOrder)
	if (length(grep(toupper("Fold"), toupper(dataName)))>0)
	{
		newData<-log2(newData)
		newData<-apply(newData, 2, function(x)
		{
			if (any(x > curBreaks))
				x[which(x>curBreaks)]<-curBreaks
			if (any(x < -curBreaks))
				x[which(x< -curBreaks)]<- -curBreaks
			return(x)
		})	
	}
	plotModuleSketch(newData, dataName, path_results, selection="Top Changing", parameters, curBreaks)
	# then make clustered heatmap from newData
	makeClusteredHeatmap(newData, parameters, designData, moduleAnnotateData, dataName, path_results, type="Top Changing", curBreaks)
}

#####
# plot heatmap with clustering
#####
makeClusteredHeatmap<-function(curData, parameters, designData, moduleAnnotateData, dataName, path_results, type="All", curBreaks)
{
	# 2 sets of clusters: samples only, and samples and modules
	if (dataName=="Percent_Negative")
		curData<-curData*(-1)	
	
	# get the color groups
	color_groups<-returnParameter(parameters, "color_groups", "statement", split=T)
		
	# cluster based on samples
	samdd<-as.dendrogram(hclust(dist(curData)))	
	samdd_ord<-order.dendrogram(samdd)
	# reorder the dendrogram based on the color groups so it's easiest to visualize grouping
	samddn<-reorder(samdd, designData[samdd_ord, color_groups[1]])
	samddn_ord<-order.dendrogram(samddn)
	
	# cluster based on modules
	moddd<-as.dendrogram(hclust(dist(t(curData))))
	moddd_ord<-order.dendrogram(moddd)

	# could have multiple color groups (showing independent variables for samples)
	groups<-as.data.frame(designData[,color_groups])
	n_cgroups<-ncol(groups)
	palette_fun<-setPalettes(n_cgroups)

	# make sure each group is of class type, factor
	factors = unlist(apply(groups, 2, function(x) as.data.frame(factor(x))), recursive=F)
  	factor_levels = lapply(factors, function(x){ as.character(levels(x)) }) 
  	factor_lengths = lapply(factor_levels, function(x) length(x))

	color_gen<-list()
	add_list<-list()
	for (i in 1:n_cgroups)
	{
		tempColor<-palette_fun[[i]](as.numeric(factor_lengths[i]))
		color_gen[[i]]<-tempColor
		add_list[[i]]<-list(rect=list(col="transparent", fill=tempColor[factors[[i]]]))
	}
	add_list<-unlist(add_list, recursive=F)
	
	# now set up the key - points = 15 means the points are boxes
	# can have multiple text and points elements (each represents one column)
	# unfortunately, cannot have more than one title 
	curKey<-list(space="left", border=TRUE, rep=FALSE)
	for (i in 1:length(factor_levels))
	{
		curKeyLen<-length(curKey)
		curKey[[curKeyLen+1]]<-list(factor_levels[[i]])
		curKey[[curKeyLen+2]]<-list(pch=15, col=color_gen[[i]])
		names(curKey)[(curKeyLen+1):(curKeyLen+2)]<-c("text","points")
	}

	canvas_size = setDimensions(nrow(curData), ncol(curData), parameters)
	#print(canvas_size)

	color_palette = colorRampPalette(c("blue", "white", "red"))(1000)
	color_palette[c(495:505)] = "#FFFFFF"
	
	# now create the plots - only show dendrogram for the samples (too messy for modules)
	# first samples only
	curFileName<-paste(path_results, project_name, "_", dataName, "_",  type, 
		"_Module_Sketch_Clustered_Samples_", format(Sys.time(), "%a %b %d %H%M%S %Y"), ".png", sep="")
	# replace any spaces with underscore in the file name
	curFileName<-gsub(" ", "_", curFileName)
	png(curFileName, width=canvas_size[[1]], height=canvas_size[[2]], units="in", res=200, pointsize=20)
	    
	# color key is the blue/red continuous color rectangle
	# legend is at top that shows the dendrogram and the independent variable coloring
	# key is the left side grob that shows the variable coloring
	lp<-levelplot(curData[samddn_ord,], xlab = NULL, ylab = NULL, 
		at = do.breaks(c(-curBreaks, curBreaks), 100), panel = panel.corrgram.2, 
		scales = list(x = list(rot = 90)),  key=curKey,
		legend = list(top = list(fun = dendrogramGrob, args = list(x = samddn, ord = samddn_ord, 
		side = "top", size = 5, size.add = 0.5, add = add_list, type = "rectangle"))), 
		colorkey = list(space = "right"), col.regions = color_palette) 
	print(lp)
	dev.off()	
	
	
	# then both samples and modules
	curFileName<-paste(path_results, project_name, "_", dataName, "_",  type, 
		"_Module_Sketch_Clustered_",
	    format(Sys.time(), "%a %b %d %H%M%S %Y"), ".png", sep="")
	# replace any spaces with underscore in the file name
	curFileName<-gsub(" ", "_", curFileName)
	png(curFileName, width=canvas_size[[1]], height=canvas_size[[2]], units="in", res=200, pointsize=20)
	    
	print(levelplot(curData[samddn_ord,moddd_ord], xlab = NULL, ylab = NULL, 
		at = do.breaks(c(-curBreaks, curBreaks), 100), panel = panel.corrgram.2, 
		scales = list(x = list(rot = 90)), key = curKey, 
		legend = list(top = list(fun = dendrogramGrob, args = list(x = samddn, ord = samddn_ord, 
		side = "top", size = 5, size.add = 0.5, add = add_list, type = "rectangle"))), 
		colorkey = list(space = "right"), col.regions = color_palette) )

	dev.off()			
}

getColor = function(x) {
 colorRampPalette(brewer.pal(7, x), interpolate = "spline")
}

setPalettes = function(n) {
  pal = list("Spectral", "PuOr", "PiYG")
  sel = pal[1:n]
  lapply(sel, getColor)
}

#####
# get the module format set up for plotting heatmap
#####
alterModule<-function(curData, parameters, designData, moduleAnnotateData, type="All", moduleOrder)
{
	# remove columns that are not expression values
	# assume the column names that we want to keep include Signal
	if (length(grep(toupper("Signal"), toupper(colnames(curData)), invert=TRUE))>0)
	{
		rmIndex<-grep(toupper("Signal"), toupper(colnames(curData)), invert=TRUE)
		curData<-curData[,-rmIndex]
	}
	# take the transpose of the matrix (so samples are rows and modules are columns)
	newData<-t(curData)
	
	# change the sample names from unique ID to sample ID
	matchIndex<-match(toupper(rownames(newData)), toupper(designData$uniqueID))
	if (any(is.na(matchIndex)))
		matchIndex<-matchIndex[-which(is.na(matchIndex))]
	newRNames<-as.character(designData$sample_id[matchIndex])
	rownames(newData)<-newRNames
	
	# change the column names to module names, including annotation
	origCNames<-colnames(newData)
	newCNames<-colnames(newData)
	addIndex<-match(toupper(newCNames), toupper(moduleAnnotateData$Module))
	newCNames[which(!is.na(addIndex))]<-paste(newCNames[which(!is.na(addIndex))], moduleAnnotateData$Function[addIndex[which(!is.na(addIndex))]])
	colnames(newData)<-newCNames
	
	# now decide if the data need to be subset
	if (type=="Annotated")
	{
		newData<-newData[,match(toupper(moduleAnnotateData$Module), toupper(origCNames))]
	}
	if (type=="Changing")
	{
		# how many of the top modules to retain
		# take the top 10% of modules
		topModule<-ceiling(0.10*nrow(moduleOrder[[1]]))
		# need to think about how this will be accomplished (previous tool uses sd)
#		keepModules<-rownames(moduleOrder)[1:topModule]
		if (length(moduleOrder)==1)
		{
			keepModules<-rownames(moduleOrder[[1]])[match(1:topModule, moduleOrder[[1]]$Rank)]
		}
		else
		{
			sumMeans<-rep(0, nrow(moduleOrder[[1]]))
			for (i in 1:length(moduleOrder))
			{
				sumMeans<-sumMeans+moduleOrder[[i]]$tStatisticMean
			}
			keepModules<-rownames(moduleOrder[[1]])[order(abs(sumMeans), decreasing=T)][1:topModule]
		}
		newData<-newData[,match(toupper(keepModules), toupper(origCNames))]
	}
	return(newData)
}

#####
# make the module sketch (module vs. sample heatmap) in a png file
#####
plotModuleSketch<-function(curData, dataName, path_results, selection="All", parameters, curBreaks)
{
	# currently, do not plan to rescale zscore to be a range of -100 to 100
	if (dataName=="Percent_Negative")
		curData<-curData*(-1)
	
	# set up the png file
	canvas_size = setDimensions(nrow(curData), ncol(curData), parameters)
	#print(canvas_size)

	color_palette = colorRampPalette(c("blue", "white", "red"))(1000)
	color_palette[c(495:505)] = "#FFFFFF"

	curFileName<-paste(path_results, project_name, "_", dataName, "_",  selection, "_Module_Sketch_",
	    format(Sys.time(), "%a %b %d %H%M%S %Y"), ".png", sep="")
	# replace spaces with underscores in the file name
	curFileName<-gsub(" ", "_", curFileName)
	png(curFileName, width=canvas_size[[1]], height=canvas_size[[2]], units="in", res=200, pointsize=20)
	    	
	curPlot<-levelplot(curData, xlab = NULL, ylab = NULL, at = do.breaks(c(-curBreaks, curBreaks), 100), 
			panel = panel.corrgram.2, scales = list(x = list(rot = 90)), colorkey = list(space = "top"), 
			col.regions = color_palette)

	print(curPlot)
	dev.off()
}

#####
# set the dimensions of the plotting device
#####
setDimensions = function(x, y, parameters) 
{
	width = x * 0.22 + 5.57
	height = y * 0.22 + 5.57

	return(list(width, height))
} 

#####
# function to set up the panel; uses grid library to set up the circles in the heatmap
#####
panel.corrgram.2 <- function(x, y, z, subscripts, at = pretty(z), scale = 0.8, ...) 
{ 
	# x and y are the integer locations for the circle centers
	x <- as.numeric(x)[subscripts]
	y <- as.numeric(y)[subscripts]
	
	# z is the numeric value to indicate what color the circle should be (percent positive, percent negative, zscore, etc.)
	z <- as.numeric(z)[subscripts]
	zcol <- level.colors(z, at = at, ...)
	
	# do not need to loop through the values
	grid.circle(x = x, y = y, r = .5 * scale, default.units = "native", gp = gpar(fill = zcol, col="white"))
} 


#####
# make pie charts of the positive and negative percentages
# for the testing just pdf
#####
plotPieChart<-function(curPosData, curNegData, path_results, selection="All", designData, moduleAnnotateData, project_name, moduleOrder)
{
	# originally scale was 1 in
	scale<-0.5
	
	# remove the first column if it is All
	if (colnames(curPosData)[1]=="All")
	{
		curPosData<-curPosData[,-1]
	}
	if (colnames(curNegData)[1]=="All")
	{
		curNegData<-curNegData[,-1]
	}	

	if (any(curPosData>1))
		curPosData<-curPosData/100
	if (any(curNegData<0))
		curNegData<-abs(curNegData)
	if (any(curNegData>1))
		curNegData<-curNegData/100
	
	if (any(round((curPosData+curNegData),10)>1))
		print("ERROR: sum of positive and negative percentages is greater than 1")
	
	# set up the pdf file: make the width 3 inches wider and the height 5 inches taller
	curWidth<-ncol(curPosData)*scale+3
	curHeight<-nrow(curPosData)*scale+5
		
	if (selection=="Annotated")
	{
		# subset curPosData and curNegData to only those modules that are annotated
		matchIndex<-match(toupper(moduleAnnotateData$Module), toupper(rownames(curPosData)))
		if (any(is.na(matchIndex)))
			matchIndex<-matchIndex[-which(is.na(matchIndex))]
		curPosData<-curPosData[matchIndex,]
		curNegData<-curNegData[matchIndex,]		
		
		# now resort from top to bottom to get the correct ordering on the plot
		curPosData<-curPosData[nrow(curPosData):1,]
		curNegData<-curNegData[nrow(curNegData):1,]
		
		yLabels<-paste(as.character(rownames(curPosData)), 
			as.character(moduleAnnotateData$Function[nrow(moduleAnnotateData):1]))
		curWidth<-ncol(curPosData)*scale+7
		curHeight<-nrow(curPosData)*scale+5
	}
	if (selection=="Top Changing")
	{
		# how many of the top modules to retain
		# take the top 10% of modules
		topModule<-ceiling(0.10*nrow(moduleOrder[[1]]))
		# need to think about how this will be accomplished (previous tool uses sd)
#		keepModules<-rownames(moduleOrder)[1:topModule]
		
		if (length(moduleOrder)==1)
		{
			keepModules<-rownames(moduleOrder[[1]])[match(1:topModule, moduleOrder[[1]]$Rank)]
		}
		else
		{
			sumMeans<-rep(0, nrow(moduleOrder[[1]]))
			for (i in 1:length(moduleOrder))
			{
				sumMeans<-sumMeans+moduleOrder[[i]]$tStatisticMean
			}
			keepModules<-rownames(moduleOrder[[1]])[order(abs(sumMeans), decreasing=T)][1:topModule]
		}

		matchIndex<-match(toupper(keepModules), toupper(rownames(curPosData)))
		if (any(is.na(matchIndex)))
			matchIndex<-matchIndex[-which(is.na(matchIndex))]
		curPosData<-curPosData[matchIndex,]
		curNegData<-curNegData[matchIndex,]	
		
		# now resort from top to bottom to get the correct ordering on the plot
		curPosData<-curPosData[nrow(curPosData):1,]
		curNegData<-curNegData[nrow(curNegData):1,]
		
		yLabels<-as.character(rownames(curPosData))
		# set up the pdf file: make the width 3 inches wider and the height 5 inches taller
		curWidth<-ncol(curPosData)*scale+3
		curHeight<-nrow(curPosData)*scale+5
	}
	if (selection=="All")
	{
		# rearrange the order so that it is ordered based on number rather than character for module
		ppModules<-rownames(curPosData)
		npModules<-rownames(curNegData)
		ppRound<-substr(ppModules, 1, 2)
		npRound<-substr(npModules, 1, 2)
		ppSecondNum<-as.numeric(unlist(lapply(strsplit(ppModules, "\\."), function(x) {x[2]})))
		npSecondNum<-as.numeric(unlist(lapply(strsplit(npModules, "\\."), function(x) {x[2]})))
		curPosData<-curPosData[order(ppRound, ppSecondNum),]
		curNegData<-curNegData[order(npRound, npSecondNum),]
		# now resort from top to bottom
		curPosData<-curPosData[nrow(curPosData):1,]
		curNegData<-curNegData[nrow(curNegData):1,]
		yLabels<-as.character(rownames(curPosData))
	}
	# special case for PK's data
	if (is.null(designData))
	{
		curWidth<-ncol(curPosData)*scale+7
	}

	curFileName<-paste(path_results, project_name, "_", "PieChart", "_",  selection, "_Module_Sketch_",
	    format(Sys.time(), "%a %b %d %H%M%S %Y"), ".png", sep="")
	# replace spaces with underscores in the file name 
	curFileName<-gsub(" ", "_", curFileName)
	
	png(curFileName, width=curWidth, height=curHeight, units="in", res=200, pointsize=20)
	# set up the viewport plot region (need more space for text when annotated)
	if (selection!="Annotated"  & !is.null(designData))
		pushViewport(viewport(width=unit(ncol(curPosData)*scale,"inches"), 
				height=unit(nrow(curPosData)*scale,"inches"), 
				x=unit(2,"inches"), y=unit(3,"inches"),
				default.units="inches", just=c("left","bottom"), name="plotRegion"))
	else if (selection=="Annotated" | is.null(designData))
		pushViewport(viewport(width=unit(ncol(curPosData)*scale,"inches"), 
				height=unit(nrow(curPosData)*scale,"inches"), 
				x=unit(6,"inches"), y=unit(3,"inches"),
				default.units="inches", just=c("left","bottom"), name="plotRegion"))
	grid.rect()
	y<-rep(seq(0.5*scale, nrow(curPosData)*scale, 1*scale), ncol(curPosData))
	x<-sort(rep(seq(0.5*scale, ncol(curPosData)*scale, 1*scale), nrow(curPosData)))
	grid.circle(x=x, y=y, r=0.4*scale, default.units="inches")
	
	# need to figure out how to rotate the text for the x-axis labels
	# need to make the text vertical for x-axis
	grid.xaxis(at=seq(0.5,ncol(curPosData),1)/ncol(curPosData), label=rep("",ncol(curPosData)))
	
	if (!is.null(dim(designData)))
	{
		matchIndex<-match(toupper(colnames(curPosData)), toupper(designData$uniqueID))
		if (any(is.na(matchIndex)))
			matchIndex<-matchIndex[-which(is.na(matchIndex))]
		grid.text(label=as.character(designData$sample_id[matchIndex]), x=seq(0.5*scale,length(matchIndex)*scale,1*scale), 
			y=rep(-1.5*scale,length(matchIndex)), rot=90, default.units="inches", gp=gpar(cex=1))
	}		
	else
	{
		grid.text(label=as.character(colnames(curPosData)), x=seq(0.5*scale,ncol(curPosData)*scale,1*scale), 
			y=rep(-1.5*scale,ncol(curPosData)), rot=90, default.units="inches", gp=gpar(cex=1))
	}

	grid.yaxis(at=seq(0.5,nrow(curPosData),1)/nrow(curPosData), label=yLabels, gp=gpar(cex=1))
	
	k<-1
	posx<-c()
	posid<-c()
	posy<-c()
	negx<-c()
	negid<-c()
	negy<-c()
	for (i in 1:ncol(curPosData))
	{
		for (j in 1:nrow(curPosData))
		{
			if (curPosData[j,i] != 0 | curNegData[j,i] != 0)
			{
				retList<-calcFillPolygon(centerx=x[k], centery=y[k], radius=0.4*scale, 
					percentPos=curPosData[j,i], 
					percentNeg=curNegData[j,i],
					units="inches", toplot="n")

				if (i==1 & j==1)
				{
					posx<-c(retList$allx)
					posid<-c(rep(k,length(retList$allx)))
					posy<-c(retList$ally)
					negx<-c(retList$allxn)
					negid<-c(rep(k,length(retList$allxn)))
					negy<-c(retList$allyn)
				}
				else
				{
					posx<-c(posx, c(retList$allx))
					posid<-c(posid, rep(k, length(retList$allx)))
					posy<-c(posy, c(retList$ally))
					negx<-c(negx, c(retList$allxn))
					negid<-c(negid, rep(k, length(retList$allxn)))
					negy<-c(negy, c(retList$allyn))
				}
			}
			k<-k+1
		}
	}
	# now plot the polygons
	if (!all(curPosData==0))
		grid.polygon(x=posx, y=posy, id=posid, gp=gpar(fill="red"), default.units="inches")
	if (!all(curNegData==0))
		grid.polygon(x=negx, y=negy, id=negid, gp=gpar(fill="blue"), default.units="inches")

	dev.off()
}

#####
# set up the coordinates to plot the polygon inside a circle
# assume we always start from the top vertical line in a circle
#####
calcFillPolygon<-function(centerx, centery, radius, percentPos, percentNeg, units="npc", toplot="y")
{
	# calculate angle
	curAnglePos<-percentPos*360
	curLength<-100
	allx<-c()
	ally<-c()
	allxn<-c()
	allyn<-c()
	# calculate this first part for all angles
	if (percentPos > 0)
	{
		# then in the top right quadrant of the circle
		# find the x-value to stop at
		opAnglePos<-90-curAnglePos
		if (opAnglePos < 0)
			opAnglePos<-0
		# cosine needs the parameter in radians
		stopAtx<-radius*cos(opAnglePos*pi/180)
		
		# now calculate the x and y values for the polygon
		# start at x value of zero for this quadrant of the circle
		startX<-0
		xvals<-seq(startX, stopAtx, length=curLength)
		yvals<-sqrt(radius^2-xvals^2)
		allx<-c(centerx, centerx, centerx+xvals)
		ally<-c(centery, centery+radius, centery+yvals)
	}
	if (percentPos > 0.25)
	{
		# in the bottom right quadrant of the circle
		# find the x-value to stop at
		opAnglePos<-curAnglePos-90
		if (opAnglePos > 90)
			opAnglePos<-90
		# cosine needs the parameter in radians
		stopAtx<-radius*cos(opAnglePos*pi/180)
		
		# now calculate the x and y values for the polygon
		# start at x value equal to the radius for this quadrant of the circle
		startX<-radius
		xvals<-seq(startX, stopAtx, length=curLength)
		yvals<-sqrt(radius^2-xvals^2)
		allx<-c(allx, centerx+xvals)
		ally<-c(ally, centery-yvals)
	}
	if (percentPos > 0.5)
	{
		# in the bottom left quadrant of the circle
		# find the x-value to stop at
		opAnglePos<-270-curAnglePos
		if (opAnglePos < 0)
			opAnglePos<-0
		stopAtx<--radius*cos(opAnglePos*pi/180)
		
		# now calculate the x and y values for the polygon
		# start at x value of zero for this quadrant of the circle
		startX<-0
		xvals<-seq(startX, stopAtx, length=curLength)
		yvals<--sqrt(radius^2-xvals^2)
		allx<-c(allx, centerx+xvals)
		ally<-c(ally, centery+yvals)
	}
	if (percentPos > 0.75)
	{
		# in the top left quadrant of the circle
		# find the x-value to stop at
		opAnglePos<-curAnglePos-270
		stopAtx<--radius*cos(opAnglePos*pi/180)
		
		# now calculate the x and y values for the polygon
		# start at x value equal to negative radius for this quadrant of the circle
		startX<- -radius
		xvals<-seq(startX, stopAtx, length=curLength)
		yvals<-sqrt(radius^2-xvals^2)
		allx<-c(allx, centerx+xvals)
		ally<-c(ally, centery+yvals)
	}
	if (toplot=="y" & percentPos > 0)
		grid.polygon(x=allx, y=ally, gp=gpar(fill="red"), default.units=units)
	
	# now work on the negative percent
	curAngleNeg<-percentNeg*360
	# calculate this first part for all angles
	if (percentNeg > 0)
	{
		# then in the top left quadrant of the circle (going backwards for negative percent)
		# find the x-value to stop at
		opAngleNeg<-90-curAngleNeg
		if (opAngleNeg < 0)
			opAngleNeg<-0
		# cosine needs the parameter in radians
		stopAtx<--radius*cos(opAngleNeg*pi/180)
		
		# now calculate the x and y values for the polygon
		# start at x value of zero for this quadrant of the circle
		startX<-0
		xvals<-seq(startX, stopAtx, length=curLength)
		yvals<-sqrt(radius^2-xvals^2)
		allxn<-c(centerx, centerx, centerx+xvals)
		allyn<-c(centery, centery+radius, centery+yvals)
	}
	if (percentNeg > 0.25)
	{
		# in the bottom left quadrant of the circle (going backwards for negative percent)
		# find the x-value to stop at
		opAngleNeg<-curAngleNeg-90
		if (opAngleNeg > 90)
			opAngleNeg<-90
		# cosine needs the parameter in radians
		stopAtx<--radius*cos(opAngleNeg*pi/180)
		
		# now calculate the x and y values for the polygon
		# start at x value equal to the radius for this quadrant of the circle
		startX<--radius
		xvals<-seq(startX, stopAtx, length=curLength)
		yvals<--sqrt(radius^2-xvals^2)
		allxn<-c(allxn, centerx+xvals)
		allyn<-c(allyn, centery+yvals)
	}
	if (percentNeg > 0.5)
	{
		# in the bottom right quadrant of the circle (going backwards for negative percent)
		# find the x-value to stop at
		opAngleNeg<-270-curAngleNeg
		if (opAngleNeg < 0)
			opAngleNeg<-0
		stopAtx<-radius*cos(opAngleNeg*pi/180)
		
		# now calculate the x and y values for the polygon
		# start at x value of zero for this quadrant of the circle
		startX<-0
		xvals<-seq(startX, stopAtx, length=curLength)
		yvals<--sqrt(radius^2-xvals^2)
		allxn<-c(allxn, centerx+xvals)
		allyn<-c(allyn, centery+yvals)

	}
	if (percentNeg > 0.75)
	{
		# in the top right quadrant of the circle
		# find the x-value to stop at
		opAngleNeg<-curAngleNeg-270
		stopAtx<-radius*cos(opAngleNeg*pi/180)
		
		# now calculate the x and y values for the polygon
		# start at x value equal to negative radius for this quadrant of the circle
		startX<-radius
		xvals<-seq(startX, stopAtx, length=curLength)
		yvals<-sqrt(radius^2-xvals^2)
		allxn<-c(allxn, centerx+xvals)
		allyn<-c(allyn, centery+yvals)
	}
	if (toplot=="y" & percentNeg > 0)
		grid.polygon(x=allxn, y=allyn, gp=gpar(fill="blue"), default.units=units)
	
	return(list(allx=allx, ally=ally, allxn=allxn, allyn=allyn))
}

#####
# create either the full 260 module plot or the 62 module plot (first 6 rounds)
# this only works for gen 2 of modules
# radiusModSize indicates whether the circle radius should be based on the number of probes assigned to the module
#####
create260ModulePlot<-function(path_results, project_name, allPercent, designData, allModules=TRUE, curSubset=NA, polygon=FALSE, moduleAssignData, radiusModSize=FALSE, type="percent")
{
	scale<-0.33
	if (!("Subset" %in% colnames(designData)))
		designData$Subset<-1
	else
	{
		curSubset<-sort(unique(designData$Subset))[curSubset]
	}
	if (polygon==FALSE)
	{
		if (type=="percent")
		{
			if (allModules==TRUE)
				fName<-"260 Module Plot "
			else
				fName<-"62 Module Plot "
		}
		if (type=="meant")
		{
			if (allModules==TRUE)
				fName<-"260 Module Plot Mean Test Statistic "
			else
				fName<-"62 Module Plot Mean Test Statistic "
		}
		if (type=="gsa")
		{
			if (allModules==TRUE)
				fName<-"260 Module Plot GSA "
			else
				fName<-"62 Module Plot GSA "
		}
	}
	else
	{
		if (allModules==TRUE)
			fName<-"260 Module Plot Pie Chart "
		else
			fName<-"62 Module Plot Pie Chart "
	}
	if (length(unique(designData$Subset))==1)
	{
		fileName<-paste(path_results, project_name, "_", fName,
		    	format(Sys.time(), "%a %b %d %H%M%S %Y"), ".png", sep="")
	}
	else
	{
		fileName<-paste(path_results, project_name, "_", fName, "_Subset_",
				curSubset, "_", 
		    	format(Sys.time(), "%a %b %d %H%M%S %Y"), ".png", sep="")
	}
	# replace spaces with underscores in the file name
	fileName<-gsub(" ", "_", fileName)
	if (allModules==TRUE)
	{
		png(fileName, width=23*scale, height=19*scale, units="in", res=200, pointsize=20)
		pushViewport(viewport(width=unit(20*scale, "inches"), height=unit(17*scale, "inches"), 
			x=unit(2*scale, "inches"), y=unit(1*scale, "inches"), default.units="inches", 
			just=c("left", "bottom")))
		vpWidth<-20*scale
		vpHeight<-17*scale
		numCol<-20
		numRow<-17
	}
	else
	{
		# draw the first 6 rounds of modules
		png(fileName, width=23*scale, height=8*scale, units="in", res=200, pointsize=20)
		pushViewport(viewport(width=unit(20*scale, "inches"), height=unit(6*scale, "inches"), 
			x=unit(2*scale, "inches"), y=unit(1*scale, "inches"), default.units="inches", 
			just=c("left", "bottom")))
		vpWidth<-20*scale
		vpHeight<-6*scale
		numCol<-20
		numRow<-6
	}

	# first draw horizontal lines
	lineType<-1
	lineWidth<-0.5
	lineColor<-"gray"
	for (i in seq(0,1,1/numRow))
		grid.lines(x=seq(0,1,1/numCol)*vpWidth, y=i*vpHeight, gp=gpar(lwd=lineWidth, lty=lineType, col=lineColor), default.units="inches")
	for (i in seq(0,1,1/numCol))
		grid.lines(x=i*vpWidth, y=seq(0,1,1/numRow)*vpHeight, gp=gpar(lwd=lineWidth, lty=lineType, col=lineColor), default.units="inches")
	
	# add rectangles of background color
	fillBackColor<-"lightgray"
	lineType<-1
	grid.rect(x=2*vpWidth/numCol, y=(numRow-1)*vpHeight/numRow, width=18*vpWidth/numCol, height=1*vpHeight/numRow, just=c("left", "bottom"), 
		gp=gpar(fill=fillBackColor, lty=lineType, col=fillBackColor), default.units="inches")
	grid.rect(x=3*vpWidth/numCol, y=(numRow-2)*vpHeight/numRow, width=17*vpWidth/numCol, height=1*vpHeight/numRow, just=c("left", "bottom"), 
		gp=gpar(fill=fillBackColor, lty=lineType, col=fillBackColor), default.units="inches")
	grid.rect(x=6*vpWidth/numCol, y=(numRow-3)*vpHeight/numRow, width=14*vpWidth/numCol, height=1*vpHeight/numRow, just=c("left", "bottom"), 
		gp=gpar(fill=fillBackColor, lty=lineType, col=fillBackColor), default.units="inches")
	grid.rect(x=16*vpWidth/numCol, y=(numRow-4)*vpHeight/numRow, width=4*vpWidth/numCol, height=1*vpHeight/numRow, just=c("left", "bottom"), 
		gp=gpar(fill=fillBackColor, lty=lineType, col=fillBackColor), default.units="inches")
	grid.rect(x=15*vpWidth/numCol, y=(numRow-5)*vpHeight/numRow, width=5*vpWidth/numCol, height=1*vpHeight/numRow, just=c("left", "bottom"), 
		gp=gpar(fill=fillBackColor, lty=lineType, col=fillBackColor), default.units="inches")
	if (allModules==TRUE)
	{	
		grid.rect(x=15*vpWidth/numCol, y=(numRow-8)*vpHeight/numRow, width=5*vpWidth/numCol, height=1*vpHeight/numRow, just=c("left", "bottom"), 
			gp=gpar(fill=fillBackColor, lty=lineType, col=fillBackColor), default.units="inches")
		grid.rect(x=11*vpWidth/numCol, y=(numRow-14)*vpHeight/numRow, width=9*vpWidth/numCol, height=1*vpHeight/numRow, just=c("left", "bottom"), 
			gp=gpar(fill=fillBackColor, lty=lineType, col=fillBackColor), default.units="inches")
		grid.rect(x=12*vpWidth/numCol, y=0*vpHeight, width=8*vpWidth/numCol, height=1*vpHeight/numRow, just=c("left", "bottom"), 
			gp=gpar(fill=fillBackColor, lty=lineType, col=fillBackColor), default.units="inches")
	}

	lineWidth<-0.8
	if (allModules==TRUE)
		lineColor<-"black"
	else
		lineColor<-"gray"
	grid.rect(x=0*vpWidth,y=(numRow-1)*vpHeight/numRow,width=2*vpWidth/numCol,height=1*vpHeight/numRow, just=c("left","bottom"), gp=gpar(lwd=lineWidth, col=lineColor, fill="transparent"), default.units="inches")
	grid.rect(x=0*vpWidth,y=(numRow-2)*vpHeight/numRow,width=3*vpWidth/numCol,height=1*vpHeight/numRow, just=c("left","bottom"), gp=gpar(lwd=lineWidth, col=lineColor, fill="transparent"), default.units="inches")
	grid.rect(x=0*vpWidth,y=(numRow-3)*vpHeight/numRow,width=6*vpWidth/numCol,height=1*vpHeight/numRow, just=c("left","bottom"), gp=gpar(lwd=lineWidth, col=lineColor, fill="transparent"), default.units="inches")
	grid.rect(x=0*vpWidth,y=(numRow-4)*vpHeight/numRow,width=16*vpWidth/numCol,height=1*vpHeight/numRow, just=c("left","bottom"), gp=gpar(lwd=lineWidth, col=lineColor, fill="transparent"), default.units="inches")
	grid.rect(x=0*vpWidth,y=(numRow-5)*vpHeight/numRow,width=15*vpWidth/numCol,height=1*vpHeight/numRow, just=c("left","bottom"), gp=gpar(lwd=lineWidth, col=lineColor, fill="transparent"), default.units="inches")
	grid.rect(x=0*vpWidth,y=(numRow-6)*vpHeight/numRow,width=20*vpWidth/numCol,height=1*vpHeight/numRow, just=c("left","bottom"), gp=gpar(lwd=lineWidth, col=lineColor, fill="transparent"), default.units="inches")
	if (allModules==TRUE)
	{
		# draw lines for modules 7, 8, and 9
		grid.lines(x=0*vpWidth, y=c(0,11/numRow)*vpHeight, gp=gpar(lwd=lineWidth, col=lineColor), default.units="inches")
		# module 7
		grid.lines(y=9*vpHeight/numRow,x=c(0,15/numCol)*vpWidth, gp=gpar(lwd=lineWidth, col=lineColor), default.units="inches")
		grid.lines(x=15*vpWidth/numCol, y=c(9/numRow,10/numRow)*vpHeight, gp=gpar(lwd=lineWidth, col=lineColor), default.units="inches")
		grid.lines(y=10*vpHeight/numRow, x=c(15/numCol,1)*vpWidth, gp=gpar(lwd=lineWidth, col=lineColor), default.units="inches")
		grid.lines(x=1*vpWidth, y=c(10/numRow,11/numRow)*vpHeight, gp=gpar(lwd=lineWidth, col=lineColor), default.units="inches")
		# module 8
		grid.lines(y=3*vpHeight/numRow, x=c(0,11/numCol)*vpWidth, gp=gpar(lwd=lineWidth, col=lineColor), default.units="inches")
		grid.lines(x=11*vpWidth/numCol, y=c(3/numRow,4/numRow)*vpHeight, gp=gpar(lwd=lineWidth, col=lineColor), default.units="inches")
		grid.lines(y=4*vpHeight/numRow, x=c(11/numCol,1)*vpWidth, gp=gpar(lwd=lineWidth, col=lineColor), default.units="inches")
		grid.lines(x=1*vpWidth, y=c(4/numRow,9/numRow)*vpHeight, gp=gpar(lwd=lineWidth, col=lineColor), default.units="inches")
		grid.lines(y=9*vpHeight/numRow, x=c(15/numCol,1)*vpWidth, gp=gpar(lwd=lineWidth, col=lineColor), default.units="inches")
		# module 9
		grid.lines(y=3*vpHeight/numRow, x=c(11/numCol,1)*vpWidth, gp=gpar(lwd=lineWidth, col=lineColor), default.units="inches")
		grid.lines(y=0*vpHeight, x=c(0,12/numCol)*vpWidth, gp=gpar(lwd=lineWidth, col=lineColor), default.units="inches")
		grid.lines(y=1*vpHeight/numRow, x=c(12/numCol,1)*vpWidth, gp=gpar(lwd=lineWidth, col=lineColor), default.units="inches")
		grid.lines(x=12*vpWidth/numCol, y=c(0,1/numRow)*vpHeight, gp=gpar(lwd=lineWidth, col=lineColor), default.units="inches")
		grid.lines(x=1*vpWidth, y=c(1/numRow,2/numRow)*vpHeight, gp=gpar(lwd=lineWidth, col=lineColor), default.units="inches")
	}
	
	# now add text - don't really want the tickmarks
	# new x-axes
	grid.text(label=1:numCol, x=seq(1/(numCol*2),1,1/numCol)*vpWidth, y=-0.3*scale, default.units="inches", gp=gpar(cex=0.6))
	grid.text(label=1:numCol, x=seq(1/(numCol*2),1,1/numCol)*vpWidth, y=(numRow+0.3)*scale, default.units="inches", gp=gpar(cex=0.6))
	# y-axis
	if (allModules==TRUE)
	{
		grid.text(label=c("41-52","21-40","1-20","101-111","81-100","61-80","41-60","21-40","1-20",
			"21-35","1-20","1-20","1-15","1-16","1-6","1-3","1-2"), y=seq(1/(numRow*2),1,1/numRow)*vpHeight, 
			x=-0.1*scale, gp=gpar(cex=0.35), just="right", default.units="inches")
		grid.text(label=c("M9","M8","M7","M6","M5","M4","M3","M2","M1"), x=-1.8*scale,
			y=c(1.5/numRow, 6/numRow,10/numRow,11.5/numRow,12.5/numRow,13.5/numRow,14.5/numRow,15.5/numRow,
			16.5/numRow)*vpHeight, just="left", default.units="inches", gp=gpar(cex=0.6))
	}
	else
	{
		grid.text(label=c("1-20","1-15","1-16","1-6","1-3","1-2"), y=seq(1/(numRow*2),1,1/numRow)*vpHeight, 
			x=-0.1*scale, gp=gpar(cex=0.35), just="right", default.units="inches")
		grid.text(label=c("M6","M5","M4","M3","M2","M1"), x=-1.8*scale, gp=gpar(cex=0.6),
			y=c(0.5/numRow,1.5/numRow,2.5/numRow,3.5/numRow,4.5/numRow,5.5/numRow)*vpHeight, just="left", default.units="inches")
	}
	grid.rect(gp=gpar(lwd=1.5,fill="transparent"))

	# now need to add circles and color
	x<-seq(1/(numCol*2), 1, 1/numCol)
	y<-seq(1/(numRow*2), 1, 1/numRow)
	if (allModules==TRUE)
	{
		allX<-c(x[1:12],x,x,x[1:11],x,x,x,x,x,x[1:15],x,x,x[1:15],x[1:16],x[1:6],x[1:3],x[1:2])
		allY<-c(rep(y[1],12),rep(y[2],20),rep(y[3],20),rep(y[4],11),rep(y[5],20),rep(y[6],20),
			rep(y[7],20),rep(y[8],20),rep(y[9],20),rep(y[10],15),rep(y[11],20),rep(y[12],20),
			rep(y[13],15),rep(y[14],16),rep(y[15],6),rep(y[16],3),rep(y[17],2))
	}
	else
	{
		allX<-c(x,x[1:15],x[1:16],x[1:6],x[1:3],x[1:2])
		allY<-c(rep(y[1],20),rep(y[2],15),rep(y[3],16),rep(y[4],6),rep(y[5],3),rep(y[6],2))
	}

	# set the fill colors if not drawing pie charts
	setFillAlpha<-function(allFill, allAlpha, curVals, allRadius, type, designData)
	{	
 		for (i in 1:length(curVals))
		{
			if (type=="meant")
			{
				# calculate the cutoffs (.1%, 1%, 5%, 10%, 20%)
				co1<-qt(0.9995, df=nrow(designData)-1)
				co2<-qt(0.995, df=nrow(designData)-1)
				co3<-qt(0.975, df=nrow(designData)-1)
				co4<-qt(0.95, df=nrow(designData)-1)
				co5<-qt(0.9, df=nrow(designData)-1)
			}
			if (type=="gsa")
			{
				co1<-0.001
				co2<-0.01
				co3<-0.03
				co4<-0.05
			}
			if (is.na(curVals[i]))
			{
				allFill<-c(allFill, "transparent")
				allAlpha<-c(allAlpha, 1)
				allRadius<-c(allRadius, 1)
			}
			else if (curVals[i]<0)
			{
				allFill<-c(allFill, "blue")
				if (type=="percent")
					allAlpha<-c(allAlpha, abs(curVals[i])/100)
				if (type=="gsa")
				{
					if (abs(curVals[i]) < co1)
						allAlpha<-c(allAlpha, 0.8)
					else
					{
						if (abs(curVals[i]) < co2)
							allAlpha<-c(allAlpha, 0.6)
						else
						{
							if (abs(curVals[i]) < co3)
								allAlpha<-c(allAlpha, 0.4)
							else
							{
								if (abs(curVals[i]) < co4)
									allAlpha<-c(allAlpha, 0.2)
								else
									allAlpha<-c(allAlpha, 0)
							}
						}
					}
				}
				if (type=="meant")
				{
					if (curVals[i] < -co1)
						allAlpha<-c(allAlpha, 1)
					else
					{
						if (curVals[i] < -co2)
							allAlpha<-c(allAlpha, 0.8)
						else
						{
							if (curVals[i] < -co3)
								allAlpha<-c(allAlpha, 0.6)
							else
							{
								if (curVals[i] < -co4)
									allAlpha<-c(allAlpha, 0.4)
								else
								{
									if (curVals[i] < -co5)
										allAlpha<-c(allAlpha, 0.2)
									else
										allAlpha<-c(allAlpha, 0)
								}
							}
						}
					}
				}
				if (abs(curVals[i]) < 20)
					allRadius<-c(allRadius, 0.8)
				else
					allRadius<-c(allRadius, .75+abs(curVals[i])*.2/80)
			}
			else if (curVals[i]>0)
			{
				allFill<-c(allFill, "red")
				if (type=="percent")
					allAlpha<-c(allAlpha, curVals[i]/100)
				if (type=="gsa")
				{
					if (abs(curVals[i]) < co1)
						allAlpha<-c(allAlpha, 0.8)
					else
					{
						if (abs(curVals[i]) < co2)
							allAlpha<-c(allAlpha, 0.6)
						else
						{
							if (abs(curVals[i]) < co3)
								allAlpha<-c(allAlpha, 0.4)
							else
							{
								if (abs(curVals[i]) < co4)
									allAlpha<-c(allAlpha, 0.2)
								else
									allAlpha<-c(allAlpha, 0)
							}
						}
					}
				}
				if (type=="meant")
				{
					if (curVals[i] > co1)
						allAlpha<-c(allAlpha, 1)
					else
					{
						if (curVals[i] > co2)
							allAlpha<-c(allAlpha, 0.8)
						else
						{
							if (curVals[i] > co3)
								allAlpha<-c(allAlpha, 0.6)
							else
							{
								if (curVals[i] > co4)
									allAlpha<-c(allAlpha, 0.4)
								else
								{
									if (curVals[i] > co5)
										allAlpha<-c(allAlpha, 0.2)
									else
										allAlpha<-c(allAlpha, 0)
								}
							}
						}
					}
				}
				if (curVals[i] < 20)
					allRadius<-c(allRadius, 0.8)
				else
					allRadius<-c(allRadius, .75+curVals[i]*.2/80)
			}
			else
			{
				allFill<-c(allFill, "white")
				allAlpha<-c(allAlpha, 0)
				allRadius<-c(allRadius, 1)
			}
		}
		return(list(allFill=allFill, allAlpha=allAlpha, allRadius=allRadius))
	}

	moduleOrder<-c()
	if (allModules==TRUE)
	{
		moduleOrder<-c(paste("M9.",41:52, sep=""),
				paste("M9.",21:40, sep=""), paste("M9.",1:20, sep=""),
				paste("M8.",101:111, sep=""), paste("M8.",81:100, sep=""),
				paste("M8.",61:80, sep=""), paste("M8.",41:60, sep=""),
				paste("M8.",21:40, sep=""), paste("M8.",1:20, sep=""),
				paste("M7.",21:35, sep=""), paste("M7.",1:20, sep=""))
	}
	moduleOrder<-c(moduleOrder, paste("M6.",1:20, sep=""),
				paste("M5.",1:15, sep=""), paste("M4.",1:16, sep=""),
				paste("M3.",1:6, sep=""), paste("M2.",1:3, sep=""), paste("M1.",1:2, sep=""))


	XlineColor<-"gray"
	# check if there are any missing modules	
	if (allModules==TRUE)
	{
		if (any(!(moduleOrder %in% moduleAssignData$module)))
		{
			# want to cross out the modules that are missing
			missModules<-moduleOrder[which(!(moduleOrder %in% moduleAssignData$module))]
			for (i in 1:length(missModules))
			{
				curMiss<-missModules[i]
				# find out where this module is on the grid
				modRound<-strsplit(curMiss, "\\.")[[1]][1]
				modNum<-as.numeric(strsplit(curMiss, "\\.")[[1]][2])
				if (modRound=="M1")
				{
					curTop<-(17/numRow)*vpHeight
					curBottom<-(16/numRow)*vpHeight
					curLeft<-(modNum-1)*vpWidth/numCol
					curRight<-((modNum-1)+1)*vpWidth/numCol
				}
				if (modRound=="M2")
				{
					curTop<-(16/numRow)*vpHeight
					curBottom<-(15/numRow)*vpHeight
					curLeft<-(modNum-1)*vpWidth/numCol
					curRight<-((modNum-1)+1)*vpWidth/numCol
				}
				if (modRound=="M3")
				{
					curTop<-(15/numRow)*vpHeight
					curBottom<-(14/numRow)*vpHeight
					curLeft<-(modNum-1)*vpWidth/numCol
					curRight<-((modNum-1)+1)*vpWidth/numCol
				}
				if (modRound=="M4")
				{
					curTop<-(14/numRow)*vpHeight
					curBottom<-(13/numRow)*vpHeight
					curLeft<-(modNum-1)*vpWidth/numCol
					curRight<-((modNum-1)+1)*vpWidth/numCol
				}
				if (modRound=="M5")
				{
					curTop<-(13/numRow)*vpHeight
					curBottom<-(12/numRow)*vpHeight
					curLeft<-(modNum-1)*vpWidth/numCol
					curRight<-((modNum-1)+1)*vpWidth/numCol
				}
				if (modRound=="M6")
				{
					curTop<-(12/numRow)*vpHeight
					curBottom<-(11/numRow)*vpHeight
					curLeft<-(modNum-1)*vpWidth/numCol
					curRight<-((modNum-1)+1)*vpWidth/numCol
				}
				if (modRound=="M7")
				{
					if (modNum<=35 & modNum>=21)
					{
						curTop<-(10/numRow)*vpHeight
						curBottom<-(9/numRow)*vpHeight
						curLeft<-(modNum-21)*vpWidth/numCol
						curRight<-((modNum-21)+1)*vpWidth/numCol
					}
					if (modNum<=20 & modNum>=1)
					{
						curTop<-(11/numRow)*vpHeight
						curBottom<-(10/numRow)*vpHeight
						curLeft<-(modNum-1)*vpWidth/numCol
						curRight<-((modNum-1)+1)*vpWidth/numCol
					}
				}
				if (modRound=="M8")
				{
					if (modNum<=111 & modNum>=101)
					{
						curTop<-(4/numRow)*vpHeight
						curBottom<-(3/numRow)*vpHeight
						curLeft<-(modNum-101)*vpWidth/numCol
						curRight<-((modNum-101)+1)*vpWidth/numCol
					}
					if (modNum<=100 & modNum>=81)
					{
						curTop<-(5/numRow)*vpHeight
						curBottom<-(4/numRow)*vpHeight
						curLeft<-(modNum-81)*vpWidth/numCol
						curRight<-((modNum-81)+1)*vpWidth/numCol
					}
					if (modNum<=80 & modNum>=61)
					{
						curTop<-(6/numRow)*vpHeight
						curBottom<-(5/numRow)*vpHeight
						curLeft<-(modNum-61)*vpWidth/numCol
						curRight<-((modNum-61)+1)*vpWidth/numCol
					}
					if (modNum<=60 & modNum>=41)
					{
						curTop<-(7/numRow)*vpHeight
						curBottom<-(6/numRow)*vpHeight
						curLeft<-(modNum-41)*vpWidth/numCol
						curRight<-((modNum-41)+1)*vpWidth/numCol
					}
					if (modNum<=40 & modNum>=21)
					{
						curTop<-(8/numRow)*vpHeight
						curBottom<-(7/numRow)*vpHeight
						curLeft<-(modNum-21)*vpWidth/numCol
						curRight<-((modNum-21)+1)*vpWidth/numCol
					}
					if (modNum<=20 & modNum>=1)
					{
						curTop<-(9/numRow)*vpHeight
						curBottom<-(8/numRow)*vpHeight
						curLeft<-(modNum-1)*vpWidth/numCol
						curRight<-((modNum-1)+1)*vpWidth/numCol
					}
				}
				if (modRound=="M9")
				{
					if (modNum<=52 & modNum>=41)
					{
						curTop<-(1/numRow)*vpHeight
						curBottom<-0*vpHeight
						curLeft<-(modNum-41)*vpWidth/numCol
						curRight<-((modNum-41)+1)*vpWidth/numCol
					}
					if (modNum<=40 & modNum>=21)
					{
						curTop<-(2/numRow)*vpHeight
						curBottom<-(1/numRow)*vpHeight
						curLeft<-(modNum-21)*vpWidth/numCol
						curRight<-((modNum-21)+1)*vpWidth/numCol
					}
					if (modNum<=20 & modNum>=1)
					{
						curTop<-(3/numRow)*vpHeight
						curBottom<-(2/numRow)*vpHeight
						curLeft<-(modNum-1)*vpWidth/numCol
						curRight<-((modNum-1)+1)*vpWidth/numCol
					}
				}
				grid.lines(x=c(curLeft, curRight), y=c(curTop, curBottom), gp=gpar(lwd=0.8, col=XlineColor), default.units="inches")
				grid.lines(x=c(curRight, curLeft), y=c(curTop, curBottom), gp=gpar(lwd=0.8, col=XlineColor), default.units="inches")
			}
		}
	}
	if (allModules==FALSE)
	{
		if (any(!(moduleOrder %in% moduleAssignData$module)))
		{
			# want to cross out the modules that are missing
			missModules<-moduleOrder[which(!(moduleOrder %in% moduleAssignData$module))]
			for (i in 1:length(missModules))
			{
				curMiss<-missModules[i]
				# find out where this module is on the grid
				modRound<-strsplit(curMiss, "\\.")[[1]][1]
				modNum<-as.numeric(strsplit(curMiss, "\\.")[[1]][2])
				if (modRound=="M1")
				{
					curTop<-(6/numRow)*vpHeight
					curBottom<-(5/numRow)*vpHeight
					curLeft<-(modNum-1)*vpWidth/numCol
					curRight<-((modNum-1)+1)*vpWidth/numCol
				}
				if (modRound=="M2")
				{
					curTop<-(5/numRow)*vpHeight
					curBottom<-(4/numRow)*vpHeight
					curLeft<-(modNum-1)*vpWidth/numCol
					curRight<-((modNum-1)+1)*vpWidth/numCol
				}
				if (modRound=="M3")
				{
					curTop<-(4/numRow)*vpHeight
					curBottom<-(3/numRow)*vpHeight
					curLeft<-(modNum-1)*vpWidth/numCol
					curRight<-((modNum-1)+1)*vpWidth/numCol
				}
				if (modRound=="M4")
				{
					curTop<-(3/numRow)*vpHeight
					curBottom<-(2/numRow)*vpHeight
					curLeft<-(modNum-1)*vpWidth/numCol
					curRight<-((modNum-1)+1)*vpWidth/numCol
				}
				if (modRound=="M5")
				{
					curTop<-(2/numRow)*vpHeight
					curBottom<-(1/numRow)*vpHeight
					curLeft<-(modNum-1)*vpWidth/numCol
					curRight<-((modNum-1)+1)*vpWidth/numCol
				}
				if (modRound=="M6")
				{
					curTop<-(1/numRow)*vpHeight
					curBottom<-(0/numRow)*vpHeight
					curLeft<-(modNum-1)*vpWidth/numCol
					curRight<-((modNum-1)+1)*vpWidth/numCol
				}
				grid.lines(x=c(curLeft, curRight), y=c(curTop, curBottom), gp=gpar(lwd=0.8, col=XlineColor), default.units="inches")
				grid.lines(x=c(curRight, curLeft), y=c(curTop, curBottom), gp=gpar(lwd=0.8, col=XlineColor), default.units="inches")
			}
		}	
	}	

	if (polygon==FALSE)
	{
		allFill<-c()
		allAlpha<-c()
		allRadius<-c()
		curNames<-names(allPercent)
		curVals<-allPercent[match(moduleOrder, curNames)]
		retList<-setFillAlpha(allFill, allAlpha, curVals, allRadius, type, designData)
		allFill<-retList$allFill
		allAlpha<-retList$allAlpha	
		allRadius<-retList$allRadius
		
		if (radiusModSize==TRUE)
		{
			# calculate the radius based on module size
			modSize<-log2(table(moduleAssignData$module))
			relSize<-modSize/max(modSize)
			allRadius<-relSize[match(moduleOrder, names(relSize))]
		}
		# for now use the same circle radius
		allRadius<-rep(1, length(allFill))
		grid.circle(allX*vpWidth, allY*vpHeight, r=1/(numRow*2.4)*allRadius*vpHeight, gp=gpar(col="white", fill=allFill, alpha=allAlpha, lwd=0), default.units="inches")
	}
	else
	{
		# change all values into inches so that the circle and polygon match (with npc they don't match!!)
		nallX<-allX*numCol*scale
		nallY<-allY*numRow*scale
		
#		grid.circle(nallX, nallY, r=(1/2.4)*scale, default.units="inches")

		curPosData<-allPercent[,1]
		curNegData<- -allPercent[,2]
		# now create polygons
		k<-1
		posx<-c()
		posid<-c()
		posy<-c()
		negx<-c()
		negid<-c()
		negy<-c()
		
		# only include modules that are not missing
		if (any(!(moduleOrder %in% unique(moduleAssignData$module))))
		{	
			missIndex<-which(!(moduleOrder %in% unique(moduleAssignData$module)))
			missModules<-moduleOrder[missIndex]
			keepIndex<-moduleOrder %in% unique(moduleAssignData$module)
			moduleOrder<-moduleOrder[keepIndex]
			missX<-nallX[missIndex]
			missY<-nallY[missIndex]
			# also need to subset the x and y positions
			nallX<-nallX[keepIndex]
			nallY<-nallY[keepIndex]
			# for missing modules, fill them in as black
#			grid.circle(missX, missY, r=(1/2.4)*scale, gp=gpar(col="transparent", fill=rep("transparent", length(missIndex)), alpha=rep(1, length(missIndex))), default.units="inches")
		}	
		for (i in 1:length(moduleOrder))
		{
			grid.circle(nallX, nallY, r=(1/2.4)*scale, default.units="inches")
			# find the correct module
			curModule<-moduleOrder[i]
			curIndex<-which(names(curPosData)==curModule)
			if (curPosData[curIndex] != 0 | curNegData[curIndex] != 0)
			{
				retList<-calcFillPolygon(centerx=nallX[i], centery=nallY[i], radius=(1/2.4)*scale,
					percentPos=curPosData[curIndex]/100, percentNeg=curNegData[curIndex]/100,
					units="inches", toplot="n")
							
				posx<-c(posx, c(retList$allx))
				posid<-c(posid, rep(i, length(retList$allx)))
				posy<-c(posy, c(retList$ally))
				negx<-c(negx, c(retList$allxn))
				negid<-c(negid, rep(i, length(retList$allxn)))
				negy<-c(negy, c(retList$allyn))
			}
		}
		# now plot the polygons if there is data to plot
		if (length(posx)>0)
			grid.polygon(x=posx, y=posy, id=posid, gp=gpar(fill="red"), default.units="inches")
		if (length(negx)>0)
			grid.polygon(x=negx, y=negy, id=negid, gp=gpar(fill="blue"), default.units="inches")
	}
	dev.off()
}

#####
# create the key that goes with the module plots
#####
createModuleKey<-function(path_results, project_name, designData)
{
	scale<-0.25
	fileName<-paste(path_results, project_name, "_", "Module_Plot_Key_",
			    format(Sys.time(), "%a %b %d %H%M%S %Y"), ".png", sep="")
	# replace underscores with spaces in the file name
	fileName<-gsub(" ", "_", fileName)
	png(fileName, width=5*scale, height=13*scale, units="in", res=200, pointsize=20)

	numCol<-2
	numRow<-11
	curWidth<-(numCol+1.5)*scale
	curHeight<-numRow*scale
	pushViewport(viewport(width=curWidth, height=curHeight, x=1.2*scale, y=1*scale, default.units="inches", 
		just=c("left", "bottom")))

	allFill<-c(rep("red", 11), rep("blue", 11))
	allAlpha<-rep(seq(0,1,.1),2)

	x<-seq(1/(numCol*2), 1, 1/numCol)
	y<-seq(1/(numRow*2), 1, 1/numRow)

	allX<-c(rep(x[1],numRow), rep(x[2],numRow))*curWidth
	allY<-c(y,y)*curHeight

	allRadius<-rep(1, length(allX))
	
	grid.circle(allX, allY, r=(1/(5*numCol))*allRadius, gp=gpar(col="white", fill=allFill, alpha=allAlpha), default.units="inches")
	
	grid.text(label=seq(0,100,10), y=seq(1/(numRow*2),1,1/numRow)*curHeight, x=-0.5*scale, default.units="inches", gp=gpar(cex=0.4))
	if (returnParameter(parameters, "MTC", "character")==TRUE)
		grid.text(label=paste("% probe sets with p < ", returnParameter(parameters, "FDR", "numeric"), sep=""), x=0.35*curWidth, y=-0.5*scale, default.units="inches", gp=gpar(cex=0.35))
	else
		grid.text(label="% probe sets with nominal p < 0.05", x=0.35*curWidth, y=-0.5*scale, default.units="inches", gp=gpar(cex=0.35))
	grid.text(label="OVER-XP", y=1.03*curHeight, x=0*curWidth, just="left", default.units="inches", gp=gpar(cex=0.3))
	grid.text(label="UNDER-XP", y=1.03*curHeight, x=0.55*curWidth, just="left", default.units="inches", gp=gpar(cex=0.3))

	dev.off()
	
	
#	# now make a separate key for mean test statistic
#	fileName<-paste(path_results, project_name, "_", "Module_Plot_Mean_Test_Statistic_Key_",
#			    format(Sys.time(), "%a %b %d %H%M%S %Y"), ".png", sep="")
#	# replace underscores with spaces in the file name
#	fileName<-gsub(" ", "_", fileName)
#	png(fileName, width=6*scale, height=6*scale, units="in", res=200, pointsize=20)
#
#	numCol<-2
#	numRow<-5
#	curWidth<-(numCol+1.5)*scale
#	curHeight<-numRow*scale
#	pushViewport(viewport(width=curWidth, height=curHeight, x=2*scale, y=0.5*scale, default.units="inches", 
#		just=c("left", "bottom")))
#	allFill<-c(rep("red",5),rep("blue",5))
#	allAlpha<-rep(seq(0.2,1,0.2),2)
#
#	x<-seq(1/(numCol*2), 1, 1/numCol)
#	y<-seq(1/(numRow*2), 1, 1/numRow)
#
#	allX<-c(rep(x[1],numRow), rep(x[2],numRow))*curWidth
#	allY<-c(y,y)*curHeight
#
#	allRadius<-rep(1, length(allX))
#	
#	grid.circle(allX, allY, r=(1/(5*numCol))*allRadius, gp=gpar(col="white", fill=allFill, alpha=allAlpha), default.units="inches")
#
#	cutoffs<-c(qt(0.90, df=nrow(designData)), qt(0.95, df=nrow(designData)), qt(0.975, df=nrow(designData)),
#		qt(0.995, df=nrow(designData)), qt(0.9995, df=nrow(designData)))
#	cutoffs<-round(cutoffs, 2)
#	grid.text(label=paste("p-value < ", c(0.20, 0.10, 0.05,0.01,0.001), "\n(t-stat > ", cutoffs, ")", sep=""), y=seq(1/(numRow*2),1,1/numRow)*curHeight, x=-0.6*scale, default.units="inches", gp=gpar(cex=0.25))
##	grid.text(label=paste("p-value < ", c(0.20, 0.10, 0.05,0.01,0.001), sep=""), y=seq(1/(numRow*2),1,1/numRow)*curHeight, x=-0.5*scale, default.units="inches", gp=gpar(cex=0.3))
##	grid.text(label=paste("t-stat > ", cutoffs, sep=""), y=seq(1/(numRow*2),1,1/numRow)*curHeight, x=(4)*scale, default.units="inches", gp=gpar(cex=0.3))
#	grid.text(label="OVER-XP", y=1.03*curHeight, x=0*curWidth, just="left", default.units="inches", gp=gpar(cex=0.3))
#	grid.text(label="UNDER-XP", y=1.03*curHeight, x=0.55*curWidth, just="left", default.units="inches", gp=gpar(cex=0.3))
#
#	dev.off()


	# now make a separate key for mean test statistic
	fileName<-paste(path_results, project_name, "_", "Module_Plot_Gene_Set_Analysis_Key_",
			    format(Sys.time(), "%a %b %d %H%M%S %Y"), ".png", sep="")
	# replace underscores with spaces in the file name
	fileName<-gsub(" ", "_", fileName)
	png(fileName, width=6*scale, height=5*scale, units="in", res=200, pointsize=20)

	numCol<-2
	numRow<-4
	curWidth<-(numCol+1.5)*scale
	curHeight<-numRow*scale
	pushViewport(viewport(width=curWidth, height=curHeight, x=2*scale, y=0.5*scale, default.units="inches", 
		just=c("left", "bottom")))
	allFill<-c(rep("red",4),rep("blue",4))
	allAlpha<-rep(seq(0.2,0.8,0.2),2)

	x<-seq(1/(numCol*2), 1, 1/numCol)
	y<-seq(1/(numRow*2), 1, 1/numRow)

	allX<-c(rep(x[1],numRow), rep(x[2],numRow))*curWidth
	allY<-c(y,y)*curHeight

	allRadius<-rep(1, length(allX))
	
	grid.circle(allX, allY, r=(1/(5*numCol))*allRadius, gp=gpar(col="white", fill=allFill, alpha=allAlpha), default.units="inches")

	grid.text(label=paste("p-value < ", c(0.05, 0.03, 0.01, 0.001), sep=""), y=seq(1/(numRow*2),1,1/numRow)*curHeight, x=-0.6*scale, default.units="inches", gp=gpar(cex=0.25))
	grid.text(label="OVER-XP", y=1.03*curHeight, x=0*curWidth, just="left", default.units="inches", gp=gpar(cex=0.3))
	grid.text(label="UNDER-XP", y=1.03*curHeight, x=0.55*curWidth, just="left", default.units="inches", gp=gpar(cex=0.3))

	dev.off()
}

createAll260ModulePlots<-function(tPercentList, tBothPercentList, path_results, project_name, designData, moduleAssignData, tMeanList, GSAList)
{
	# make module plots - assumes using generation 2 of modules!!!
	for (i in 1:length(tPercentList))
	{
		create260ModulePlot(path_results, project_name, allPercent=tPercentList[[i]],
			designData, allModules=TRUE, curSubset=i, polygon=FALSE, moduleAssignData)
		create260ModulePlot(path_results, project_name, allPercent=tPercentList[[i]],
			designData, allModules=FALSE, curSubset=i, polygon=FALSE, moduleAssignData)
		create260ModulePlot(path_results, project_name, allPercent=tBothPercentList[[i]],
			designData, allModules=TRUE, curSubset=i, polygon=TRUE, moduleAssignData)
		create260ModulePlot(path_results, project_name, allPercent=tBothPercentList[[i]],
			designData, allModules=FALSE, curSubset=i, polygon=TRUE, moduleAssignData)
		
#		create260ModulePlot(path_results, project_name, allPercent=tMeanList[[i]],
#			designData, allModules=TRUE, curSubset=i, polygon=FALSE, moduleAssignData, type="meant")
#		create260ModulePlot(path_results, project_name, allPercent=tMeanList[[i]],
#			designData, allModules=FALSE, curSubset=i, polygon=FALSE, moduleAssignData, type="meant")

		gsaPercent<-GSAList[[i]]$smallPval
		if (any(gsaPercent==0))
			gsaPercent[which(gsaPercent==0)]<-0.0001
		gsaPercent[which(GSAList[[i]]$GSAscore < 0)]<- -gsaPercent[which(GSAList[[i]]$GSAscore < 0)]
		names(gsaPercent)<-rownames(GSAList[[i]])
		create260ModulePlot(path_results, project_name, allPercent=gsaPercent,
			designData, allModules=TRUE, curSubset=i, polygon=FALSE, moduleAssignData, type="gsa")
		create260ModulePlot(path_results, project_name, allPercent=gsaPercent,
			designData, allModules=FALSE, curSubset=i, polygon=FALSE, moduleAssignData, type="gsa")		
		# write output to csv elsewhere (don't want it in a plotting function)
#		# write the results to a csv file
#		if (length(tPercentList)==1)
#			fileName<-paste(path_results, project_name, "_", "260_Module_Group_Comparison_",
#			    format(Sys.time(), "%a %b %d %H%M%S %Y"), ".csv", sep="")
#		else
#		{
#			curSubset<-sort(unique(designData$Subset))[i]
#			fileName<-paste(path_results, project_name, "_", "260_Module_Group_Comparison_Subset_",
#				curSubset, "_",
#			    format(Sys.time(), "%a %b %d %H%M%S %Y"), ".csv", sep="")		}	
#		# replace spaces with underscores in the file name
#		fileName<-gsub(" ", "_", fileName)
#		dataToWrite<-tBothPercentList[[i]]
#		colnames(dataToWrite)<-c("positivePercent", "negativePercent")
#		write.csv(dataToWrite, fileName)
#		
##		# also write out the results for the tMeanList
##		if (length(tMeanList)==1)
##			fileName<-paste(path_results, project_name, "_", "Mean_tStat_260_Module_Group_Comparison_",
##			    format(Sys.time(), "%a %b %d %H%M%S %Y"), ".csv", sep="")
##		else
##		{
##			curSubset<-sort(unique(designData$Subset))[i]
##			fileName<-paste(path_results, project_name, "_", "Mean_tStat_260_Module_Group_Comparison_Subset_",
##				curSubset, "_",
##			    format(Sys.time(), "%a %b %d %H%M%S %Y"), ".csv", sep="")
##		}	
##		# replace spaces with underscores in the file name
##		fileName<-gsub(" ", "_", fileName)
##		dataToWrite<-as.matrix(tMeanList[[i]])
##		colnames(dataToWrite)<-c("meanTStat")
##		write.csv(dataToWrite, fileName)		
#
#		# also write out the results for GSAList
#		if (length(GSAList)==1)
#			fileName<-paste(path_results, project_name, "_", "Gene_Set_Analysis_260_Module_Group_Comparison_",
#			    format(Sys.time(), "%a %b %d %H%M%S %Y"), ".csv", sep="")
#		else
#		{
#			curSubset<-sort(unique(designData$Subset))[i]
#			fileName<-paste(path_results, project_name, "_", "Gene_Set_Analysis_260_Module_Group_Comparison_Subset_",
#				curSubset, "_",
#			    format(Sys.time(), "%a %b %d %H%M%S %Y"), ".csv", sep="")		}	
#		# replace spaces with underscores in the file name
#		fileName<-gsub(" ", "_", fileName)
#		dataToWrite<-as.matrix(GSAList[[i]][,c(1,4)])
#		colnames(dataToWrite)[2]<-"pvalue"
#		write.csv(dataToWrite, fileName)
	}
	createModuleKey(path_results, project_name, designData)
}

#######
# create a full grid for the gen3 whole blood modules
#######
create382ModulePlot<-function(path_results, project_name, allPercent, designData, allModules=TRUE, curSubset=NA, polygon=FALSE, moduleAssignData)
{
	scale<-0.33
	if (!("Subset" %in% colnames(designData)))
		designData$Subset<-1
	else
	{
		curSubset<-sort(unique(designData$Subset))[curSubset]
	}
	if (polygon==FALSE)
	{
		if (allModules==TRUE)
			fName<-"382 Module Plot "
		else
			fName<-"29 Module Plot "
	}
	else
	{
		if (allModules==TRUE)
			fName<-"382 Module Plot Pie Chart "
		else
			fName<-"29 Module Plot Pie Chart "
	}
	if (length(unique(designData$Subset))==1)
	{
		fileName<-paste(path_results, project_name, "_", fName,
		    	format(Sys.time(), "%a %b %d %H%M%S %Y"), ".png", sep="")
	}
	else
	{
		fileName<-paste(path_results, project_name, "_", fName, "_Subset_",
				curSubset, "_", 
		    	format(Sys.time(), "%a %b %d %H%M%S %Y"), ".png", sep="")
	}
	# replace spaces with underscores in the file name
	fileName<-gsub(" ", "_", fileName)
	if (allModules==TRUE)
	{
		png(fileName, width=23*scale, height=28*scale, units="in", res=200, pointsize=20)
		pushViewport(viewport(width=unit(20*scale, "inches"), height=unit(26*scale, "inches"), 
			x=unit(2*scale, "inches"), y=unit(1*scale, "inches"), default.units="inches", 
			just=c("left", "bottom")))
		numCol<-20
		numRow<-26
		vpWidth<-numCol*scale
		vpHeight<-numRow*scale
	}
	else
	{
		# draw the first 6 rounds of modules (rounds 3, 8, 9, 10, 11, and 12)
		png(fileName, width=23*scale, height=8*scale, units="in", res=200, pointsize=20)
		pushViewport(viewport(width=unit(20*scale, "inches"), height=unit(6*scale, "inches"), 
			x=unit(2*scale, "inches"), y=unit(1*scale, "inches"), default.units="inches", 
			just=c("left", "bottom")))
		numCol<-20
		numRow<-6
		vpWidth<-numCol*scale
		vpHeight<-numRow*scale
	}

	# first draw horizontal lines
	lineType<-1
	lineColor<-"gray"
	for (i in seq(0,1,1/numRow))
		grid.lines(x=seq(0,1,1/numCol)*vpWidth, y=i*vpHeight, gp=gpar(lty=lineType, col=lineColor), default.units="inches")
	for (i in seq(0,1,1/numCol))
		grid.lines(x=i*vpWidth, y=seq(0,1,1/numRow)*vpHeight, gp=gpar(lty=lineType, col=lineColor), default.units="inches")
	
	# add rectangles of background color
	fillBackColor<-"lightgray"
	lineType<-1
	grid.rect(x=1*vpWidth/numCol, y=(numRow-1)*vpHeight/numRow, width=19*vpWidth/numCol, height=1*vpHeight/numRow, just=c("left", "bottom"), 
		gp=gpar(fill=fillBackColor, lty=lineType, col=fillBackColor), default.units="inches")
	grid.rect(x=3*vpWidth/numCol, y=(numRow-2)*vpHeight/numRow, width=17*vpWidth/numCol, height=1*vpHeight/numRow, just=c("left", "bottom"), 
		gp=gpar(fill=fillBackColor, lty=lineType, col=fillBackColor), default.units="inches")
	grid.rect(x=2*vpWidth/numCol, y=(numRow-3)*vpHeight/numRow, width=18*vpWidth/numCol, height=1*vpHeight/numRow, just=c("left", "bottom"), 
		gp=gpar(fill=fillBackColor, lty=lineType, col=fillBackColor), default.units="inches")
	grid.rect(x=4*vpWidth/numCol, y=(numRow-4)*vpHeight/numRow, width=16*vpWidth/numCol, height=1*vpHeight/numRow, just=c("left", "bottom"), 
		gp=gpar(fill=fillBackColor, lty=lineType, col=fillBackColor), default.units="inches")
	grid.rect(x=4*vpWidth/numCol, y=(numRow-5)*vpHeight/numRow, width=16*vpWidth/numCol, height=1*vpHeight/numRow, just=c("left", "bottom"), 
		gp=gpar(fill=fillBackColor, lty=lineType, col=fillBackColor), default.units="inches")
	grid.rect(x=15*vpWidth/numCol, y=(numRow-6)*vpHeight/numRow, width=5*vpWidth/numCol, height=1*vpHeight/numRow, just=c("left", "bottom"),
		gp=gpar(fill=fillBackColor, lty=lineType, col=fillBackColor), default.units="inches")
	if (allModules==TRUE)
	{	
		grid.rect(x=12*vpWidth/numCol, y=(numRow-8)*vpHeight/numRow, width=8*vpWidth/numCol, height=1*vpHeight/numRow, just=c("left", "bottom"), 
			gp=gpar(fill=fillBackColor, lty=lineType, col=fillBackColor), default.units="inches")
		grid.rect(x=3*vpWidth/numCol, y=(numRow-13)*vpHeight/numRow, width=17*vpWidth/numCol, height=1*vpHeight/numRow, just=c("left", "bottom"), 
			gp=gpar(fill=fillBackColor, lty=lineType, col=fillBackColor), default.units="inches")
		grid.rect(x=7*vpWidth/numCol, y=(numRow-20)*vpHeight/numRow, width=13*vpWidth/numCol, height=1*vpHeight/numRow, just=c("left", "bottom"), 
			gp=gpar(fill=fillBackColor, lty=lineType, col=fillBackColor), default.units="inches")
		grid.rect(x=11*vpWidth/numCol, y=0*vpHeight, width=9*vpWidth/numCol, height=1*vpHeight/numRow, just=c("left", "bottom"), 
			gp=gpar(fill=fillBackColor, lty=lineType, col=fillBackColor), default.units="inches")
	}

	lineWidth<-1
	if (allModules==TRUE)
		lineColor<-"black"
	else
		lineColor<-"gray"
	grid.rect(x=0*vpWidth,y=(numRow-1)*vpHeight/numRow,width=1*vpWidth/numCol,height=1*vpHeight/numRow, just=c("left","bottom"), gp=gpar(lwd=lineWidth, col=lineColor, fill="transparent"), default.units="inches")
	grid.rect(x=0*vpWidth,y=(numRow-2)*vpHeight/numRow,width=3*vpWidth/numCol,height=1*vpHeight/numRow, just=c("left","bottom"), gp=gpar(lwd=lineWidth, col=lineColor, fill="transparent"), default.units="inches")
	grid.rect(x=0*vpWidth,y=(numRow-3)*vpHeight/numRow,width=2*vpWidth/numCol,height=1*vpHeight/numRow, just=c("left","bottom"), gp=gpar(lwd=lineWidth, col=lineColor, fill="transparent"), default.units="inches")
	grid.rect(x=0*vpWidth,y=(numRow-4)*vpHeight/numRow,width=4*vpWidth/numCol,height=1*vpHeight/numRow, just=c("left","bottom"), gp=gpar(lwd=lineWidth, col=lineColor, fill="transparent"), default.units="inches")
	grid.rect(x=0*vpWidth,y=(numRow-5)*vpHeight/numRow,width=4*vpWidth/numCol,height=1*vpHeight/numRow, just=c("left","bottom"), gp=gpar(lwd=lineWidth, col=lineColor, fill="transparent"), default.units="inches")
	grid.rect(x=0*vpWidth,y=(numRow-6)*vpHeight/numRow,width=15*vpWidth/numCol,height=1*vpHeight/numRow, just=c("left","bottom"), gp=gpar(lwd=lineWidth, col=lineColor, fill="transparent"), default.units="inches")
	if (allModules==TRUE)
	{
		# draw lines for modules 13, 14, 15, and 16
		grid.lines(x=0*vpWidth, y=c(0,20/numRow)*vpHeight, gp=gpar(lwd=lineWidth, col=lineColor), default.units="inches")
		# module 13
		grid.lines(y=18*vpHeight/numRow, x=c(0,12/numCol)*vpWidth, gp=gpar(lwd=lineWidth, col=lineColor), default.units="inches")
		grid.lines(x=12*vpWidth/numCol, y=c(18/numRow, 19/numRow)*vpHeight, gp=gpar(lwd=lineWidth, col=lineColor), default.units="inches")
		grid.lines(y=19*vpHeight/numRow, x=c(12/numCol,1)*vpWidth, gp=gpar(lwd=lineWidth, col=lineColor), default.units="inches")
		grid.lines(x=1*vpWidth, y=c(12/numCol, 13/numCol)*vpHeight, gp=gpar(lwd=lineWidth, col=lineColor), default.units="inches")
		grid.lines(y=20*vpHeight/numRow, x=c(15/numCol,1)*vpWidth, gp=gpar(lwd=lineWidth, col=lineColor), default.units="inches")
		# module 14
		grid.lines(y=13*vpHeight/numRow, x=c(0,3/numCol)*vpWidth, gp=gpar(lwd=lineWidth, col=lineColor), default.units="inches")
		grid.lines(x=3*vpWidth/numCol, y=c(13/numRow, 14/numRow)*vpHeight, gp=gpar(lwd=lineWidth, col=lineColor), default.units="inches")
		grid.lines(y=14*vpHeight/numRow, x=c(3/numCol, 1)*vpWidth, gp=gpar(lwd=lineWidth, col=lineColor), default.units="inches")
		grid.lines(x=1*vpWidth, y=c(14/numRow, 18/numRow)*vpHeight, gp=gpar(lwd=lineWidth, col=lineColor), default.units="inches")
		grid.lines(y=18*vpHeight/numRow, x=c(12/numCol, 1)*vpWidth, gp=gpar(lwd=lineWidth, col=lineColor), default.units="inches")
		# module 15
		grid.lines(y=6*vpHeight/numRow, x=c(0,7/numCol)*vpWidth, gp=gpar(lwd=lineWidth, col=lineColor), default.units="inches")
		grid.lines(x=7*vpWidth/numCol, y=c(6/numRow, 7/numRow)*vpHeight, gp=gpar(lwd=lineWidth, col=lineColor), default.units="inches")
		grid.lines(y=7*vpHeight/numRow, x=c(7/numCol, 1)*vpWidth, gp=gpar(lwd=lineWidth, col=lineColor), default.units="inches")
		grid.lines(x=1*vpWidth, y=c(7/numRow, 13/numRow)*vpHeight, gp=gpar(lwd=lineWidth, col=lineColor), default.units="inches")
		grid.lines(y=13*vpHeight/numRow, x=c(3/numCol, 1)*vpWidth, gp=gpar(lwd=lineWidth, col=lineColor), default.units="inches")
		# module 16
		grid.lines(y=0*vpHeight, x=c(0,11/numCol)*vpWidth, gp=gpar(lwd=lineWidth, col=lineColor), default.units="inches")
		grid.lines(x=11*vpWidth/numCol, y=c(0, 1/numRow)*vpHeight, gp=gpar(lwd=lineWidth, col=lineColor), default.units="inches")
		grid.lines(y=1*vpHeight/numRow, x=c(11/numCol, 1)*vpWidth, gp=gpar(lwd=lineWidth, col=lineColor), default.units="inches")
		grid.lines(x=1*vpWidth, y=c(1/numRow, 6/numRow)*vpHeight, gp=gpar(lwd=lineWidth, col=lineColor), default.units="inches")
		grid.lines(y=6*vpHeight/numRow, x=c(7/numCol, 1)*vpWidth, gp=gpar(lwd=lineWidth, col=lineColor), default.units="inches")
	}
	
	# now add text - don't really want the tickmarks
	#grid.xaxis(at=seq(1/(numCol*2),1,1/numCol), label=1:numCol)
	#grid.xaxis(at=seq(1/(numCol*2),1,1/numCol), label=1:numCol, main=FALSE)
	# new x-axes
	grid.text(label=1:numCol, x=seq(1/(numCol*2),1,1/numCol)*vpWidth, y=-0.3*scale, default.units="inches", gp=gpar(cex=0.6))
	grid.text(label=1:numCol, x=seq(1/(numCol*2),1,1/numCol)*vpWidth, y=(numRow+0.3)*scale, default.units="inches", gp=gpar(cex=0.6))
	# y-axis
	if (allModules==TRUE)
	{
		grid.text(label=c("101-111","81-100","61-80","41-60","21-40","1-20", "121-127","101-120",
			"81-100","61-80","41-60","21-40","1-20", "81-83","61-80","41-60","21-40","1-20", 
			"21-32","1-20", "1-15", "1-4", "1-4", "1-2", "1-3", "1"), y=seq(1/(numRow*2),1,1/numRow)*vpHeight, 
			x=-0.1*scale, default.units="inches", gp=gpar(cex=0.35), just="right")
		grid.text(label=c("M16","M15","M14","M13","M12","M11","M10","M9","M8","M3"), x=-1.8*scale, default.units="inches",
			y=c(3/numRow, 9.5/numRow, 15.5/numRow, 19/numRow, 20.5/numRow, 21.5/numRow, 22.5/numRow, 23.5/numRow, 
			24.5/numRow, 25.5/numRow)*vpHeight, just="left", gp=gpar(cex=0.6))
	}
	else
	{
		grid.text(label=c("1-15", "1-4", "1-4", "1-2", "1-3", "1"), y=seq(1/(numRow*2),1,1/numRow)*vpHeight, 
			x=-0.1*scale, default.units="inches", gp=gpar(cex=0.35), just="right")
		grid.text(label=c("M12","M11","M10","M9","M8","M3"), x=-1.8*scale, default.units="inches",
			y=c(0.5/numRow,1.5/numRow,2.5/numRow,3.5/numRow,4.5/numRow,5.5/numRow)*vpHeight, just="left", gp=gpar(cex=0.6))
	}
	grid.rect(gp=gpar(lwd=1.5), fill="transparent")

	# now need to add circles and color
	x<-seq(1/(numCol*2), 1, 1/numCol)
	y<-seq(1/(numRow*2), 1, 1/numRow)
	if (allModules==TRUE)
	{
		allX<-c(x[1:11],x,x,x,x,x, x[1:7],x,x,x,x,x,x, x[1:3],x,x,x,x, x[1:12],x, 
			x[1:15], x[1:4], x[1:4], x[1:2], x[1:3], x[1])
		allY<-c(rep(y[1],11), rep(y[2],20), rep(y[3],20), rep(y[4],20), rep(y[5],20), rep(y[6],20),
			rep(y[7],7), rep(y[8],20), rep(y[9],20), rep(y[10],20), rep(y[11],20), rep(y[12],20), rep(y[13],20),
			rep(y[14],3), rep(y[15],20), rep(y[16],20), rep(y[17],20), rep(y[18],20),
			rep(y[19],12), rep(y[20],20),
			rep(y[21],15), rep(y[22],4), rep(y[23],4), rep(y[24],2), rep(y[25],3), rep(y[26],1))
	}
	else
	{
		allX<-c(x[1:15], x[1:4], x[1:4], x[1:2], x[1:3], x[1])
		allY<-c(rep(y[1],15),rep(y[2],4),rep(y[3],4),rep(y[4],2),rep(y[5],3),rep(y[6],1))
	}


	# set the fill colors if not drawing pie charts
	setFillAlpha<-function(allFill, allAlpha, curVals, allRadius)
	{	
 		for (i in 1:length(curVals))
		{
			if (is.na(curVals[i]))
			{
				allFill<-c(allFill, "black")
				allAlpha<-c(allAlpha, 1)
			}
			else if (curVals[i]<0)
			{
				allFill<-c(allFill, "blue")
				allAlpha<-c(allAlpha, abs(curVals[i])/100)
			}
			else if (curVals[i]>0)
			{
				allFill<-c(allFill, "red")
				allAlpha<-c(allAlpha, curVals[i]/100)
			}
			else
			{
				allFill<-c(allFill, "white")
				allAlpha<-c(allAlpha, 0)
			}
		}
		return(list(allFill=allFill, allAlpha=allAlpha))
	}

	moduleOrder<-c()
	if (allModules==TRUE)
	{
		moduleOrder<-c(paste("M16.",101:111,sep=""),
				paste("M16.",81:100,sep=""), paste("M16.",61:80,sep=""),
				paste("M16.",41:60,sep=""), paste("M16.",21:40,sep=""),
				paste("M16.",1:20,sep=""), paste("M15.",121:127,sep=""),
				paste("M15.",101:120,sep=""), paste("M15.",81:100,sep=""), 
				paste("M15.",61:80,sep=""), paste("M15.",41:60,sep=""), 
				paste("M15.",21:40,sep=""), paste("M15.",1:20,sep=""),
				paste("M14.",81:83,sep=""), paste("M14.",61:80,sep=""),
				paste("M14.",41:60,sep=""), paste("M14.",21:40,sep=""),
				paste("M14.",1:20,sep=""), paste("M13.",21:32,sep=""),
				paste("M13.",1:20,sep=""))
	}
	moduleOrder<-c(moduleOrder, paste("M12.",1:15, sep=""),
				paste("M11.",1:4, sep=""), paste("M10.",1:4, sep=""),
				paste("M9.",1:2, sep=""), paste("M8.",1:3, sep=""), paste("M3.",1, sep=""))
	if (polygon==FALSE)
	{
		allFill<-c()
		allAlpha<-c()
		curNames<-names(allPercent)
		curVals<-allPercent[match(moduleOrder, curNames)]
		retList<-setFillAlpha(allFill, allAlpha, curVals)
		allFill<-retList$allFill
		allAlpha<-retList$allAlpha	
		
		# calculate the radius based on module size
		modSize<-log2(table(moduleAssignData$module))
		relSize<-modSize/max(modSize)
		allRadius<-relSize[match(moduleOrder, names(relSize))]
		
		# for now use the same circle radius
		allRadius<-rep(1, length(allFill))
		
		if (allModules==TRUE)
			grid.circle(allX*vpWidth, allY*vpHeight, r=1/(numRow*1.8)*allRadius*vpHeight, gp=gpar(col="white", fill=allFill, alpha=allAlpha, lwd=0), default.units="inches")
		else
			grid.circle(allX*vpWidth, allY*vpHeight, r=1/(numRow*2.4)*allRadius*vpHeight, gp=gpar(col="white", fill=allFill, alpha=allAlpha, lwd=0), default.units="inches")
	}
	else
	{
		# change all values into inches so that the circle and polygon match (with npc they don't match!!)
		nallX<-allX*numCol*scale
		nallY<-allY*numRow*scale
		
		grid.circle(nallX, nallY, r=(1/2.4)*scale, default.units="inches")

		curPosData<-allPercent[,1]
		curNegData<- -allPercent[,2]
		# now create polygons
		k<-1
		posx<-c()
		posid<-c()
		posy<-c()
		negx<-c()
		negid<-c()
		negy<-c()
		
		# only include modules that are not missing
		if (any(!(moduleOrder %in% unique(moduleAssignData$module))))
		{	
			missIndex<-which(!(moduleOrder %in% unique(moduleAssignData$module)))
			missModules<-moduleOrder[missIndex]
			keepIndex<-moduleOrder %in% unique(moduleAssignData$module)
			moduleOrder<-moduleOrder[keepIndex]
			missX<-nallX[missIndex]
			missY<-nallY[missIndex]
			# need to subset the x and y positions for the polygons
			nallX<-nallX[keepIndex]
			nallY<-nallY[keepIndex]
			# for missing modules, want to fill them in as black
			grid.circle(missX, missY, r=(1/2.4)*scale, gp=gpar(col="white", fill=rep("black", length(missIndex)), alpha=rep(1, length(missIndex))), default.units="inches")
		}	
		for (i in 1:length(moduleOrder))
		{
			# find the correct module
			curModule<-moduleOrder[i]
			curIndex<-which(names(curPosData)==curModule)
			if (curPosData[curIndex] != 0 | curNegData[curIndex] != 0)
			{
				retList<-calcFillPolygon(centerx=nallX[i], centery=nallY[i], radius=(1/2.4)*scale,
					percentPos=curPosData[curIndex]/100, percentNeg=curNegData[curIndex]/100,
					units="inches", toplot="n")
							
				posx<-c(posx, c(retList$allx))
				posid<-c(posid, rep(i, length(retList$allx)))
				posy<-c(posy, c(retList$ally))
				negx<-c(negx, c(retList$allxn))
				negid<-c(negid, rep(i, length(retList$allxn)))
				negy<-c(negy, c(retList$allyn))
			}
		}
		# now plot the polygons if there is data to plot
		if (length(posx)>0)
			grid.polygon(x=posx, y=posy, id=posid, gp=gpar(fill="red"), default.units="inches")
		if (length(negx)>0)
			grid.polygon(x=negx, y=negy, id=negid, gp=gpar(fill="blue"), default.units="inches")
	}
	dev.off()
}

