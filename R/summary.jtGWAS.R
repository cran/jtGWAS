summary.jtGWAS <- function(object, marker2Print = 1:10, SNP2Print =1:10, printP = TRUE, outTopN=NA, subObj=FALSE,...)
{
	#-----------------------------------------------------------------------#
	# Check if the user set SNP2Print and marker2Pring 
	SNP.is.default.for.top.hit <- .checkSNPDefault(SNP2Print)
	marker.is.default <- .checkMarkerDefault(marker2Print)
	# get the colnames 
	markerNames <- colnames(object$J)
	
	#-----------------------------------------------------------------------#
	# check the validaty of marker range
	marker2Print <- .markerRangeCheck(object, marker2Print)
	# check the validaty of SNP range
	SNP2Print <-.SNPRangeCheck(object,SNP2Print)
	
	#-----------------------------------------------------------------------#
	# compute the p-values
	if(attr(object,'standardized'))
		p.values <- 2*pnorm(-abs(object$J))
	else
		printP = FALSE
	#-----------------------------------------------------------------------#
	# print out note
	.note.title.print(object, outTopN,printP, marker.is.default, SNP.is.default.for.top.hit)
	
	#-----------------------------------------------------------------------#
	# check print out the summary
	#-----------------------------------------------------------------------#
	if(is.na(attr(object,'outTopN'))&&is.na(outTopN))
	{
		# Determine the SNPs that are required for out print 
		# for special case that the user provide the SNP names.
		if(is.character(SNP2Print))
		{
			# SNP name not in the range warning
			if(length(SNP2Print[!(SNP2Print %in% object$gSnipID)]))
			{	
				cat("\nWARNING: specified SNP IDs/range: ")	
				cat(SNP2Print[!(SNP2Print %in% object$gSnipID)]) 
				cat("\nare not found, results for them are omitted.\n")
			}
			# keep SNP that appears in the input data.
			SNP2Print <- intersect(SNP2Print, object$gSnipID)
		}
		if(!attr(object,'standardized'))
			 print(object$J[SNP2Print,marker2Print])
		else{
			if(printP)
           		print(p.values[SNP2Print,marker2Print])
       		else
           		print(object$J[SNP2Print,marker2Print])
		}
		object$J <-object$J[SNP2Print,marker2Print]
		object$gSnipID <- object$gSnipID[SNP2Print]
		colnames(object$J) <-colnames(object$J)[marker2Print]
		if(subObj)	
			return(object)
	}
	
	#-----------------------------------------------------------------------#	
	if(is.na(attr(object,'outTopN'))&& !is.na(outTopN))		
	{
		if(is.numeric(marker2Print))
			marker.name.to.print <- markerNames[marker2Print]
		else
			marker.name.to.print <- marker2Print
		
		summary.mat <- .sortTopN(object,marker2Print, outTopN, printP)
		if(!attr(object,'standardized'))
			.printJStatisticsFormat(summary.mat, marker.name.to.print)	
		else{
			if(printP)
       	    	.printPvaluesFormat(summary.mat, marker.name.to.print)
       		else
           		.printJstarFormat(summary.mat, marker.name.to.print)
		}
		subobject <- .sort2saveObj(object, marker2Print,outTopN)
		colnames(subobject$J) <- marker.name.to.print
		colnames(subobject$gSnipID) <- marker.name.to.print
		if(subObj)	
			return(subobject)
	}		
	#-----------------------------------------------------------------------#
    if(!is.na(attr(object,'outTopN')))	
	{
		if(SNP.is.default.for.top.hit)
			SNP2Print = seq(1:nrow(object$gSnipID))
		
		if(is.numeric(marker2Print))
			marker.name.to.print <- markerNames[marker2Print]
		else
			marker.name.to.print <- marker2Print
		
		summary.mat <- .summary.mat(object, marker2Print, SNP2Print, printP)	
		
		if(!attr(object,'standardized'))
			.printJStatisticsFormat(summary.mat, marker.name.to.print)	
		else{
			if(printP)
       	    	.printPvaluesFormat(summary.mat, marker.name.to.print)
       		else
           		.printJstarFormat(summary.mat, marker.name.to.print)
		}
		object$J <- object$J[SNP2Print,marker2Print]
		object$gSnipID <- object$gSnipID[SNP2Print,marker2Print]	
		colnames(object$J) <-colnames(object$J)[marker2Print]
		colnames(object$gSnipID) <-colnames(object$gSnipID)[marker2Print] 
		if(subObj)
			return(object)
	}
}
#------------------------------------------------------------------------#
.checkSNPDefault <- function(SNP2Print)
{
	SNP.is.default.for.top.hit =T
	if(is.numeric(SNP2Print) && !(SNP2Print[length(SNP2Print)]==10&&length(SNP2Print)==10))
		SNP.is.default.for.top.hit = F
	if(length(SNP2Print)==1)
	{
		if(is.na(SNP2Print))
			SNP.is.default.for.top.hit =F
	}
  if(is.character(SNP2Print))
		SNP.is.default.for.top.hit = F
	return(SNP.is.default.for.top.hit)
}
#------------------------------------------------------------------------#
.checkMarkerDefault <- function(marker2Print)
{
	marker.is.default = T
	if(is.numeric(marker2Print) && !(marker2Print[length(marker2Print)]==10&&length(marker2Print)==10))
        marker.is.default = F
	if(is.character(marker2Print))
		marker.is.default = F
	return(marker.is.default)
}

#------------------------------------------------------------------------#
.sortTopN <- function(object, marker2Print, outTopN, printP)
{
 	markerNames <- colnames(object$J)
	col.name <- NULL
	summary.mat <-NULL
	nSNP <- nrow(object$J)
    SNPIDs <-rownames(object$J)
	pvalues <- 2*pnorm(-abs(object$J))
	
	for(i in marker2Print )
    {
        # generate colname names for the dataframe
        col.name <- cbind(col.name, "SNPID")
        if(printP)
            col.name <- cbind(col.name, "P-value")
        else{
            if(attr(object,'standardized'))
				col.name <- cbind(col.name, "J*")
			else
				col.name <- cbind(col.name, "J")
		}
        # create summary data frames    
        temp <- object$J[,i]
        temp <- cbind(temp,pvalues[,i])

        rownames(temp) <- SNPIDs
        colnames(temp) <- colnames(temp) <- c("J","pvalues")
        temp <-as.data.frame(temp)
        temp <- temp[order(-abs(temp$J)),]

        #rownames(order(object$J[,i]))
        summary.mat <- cbind(summary.mat, rownames(temp)[1:outTopN])
        if(printP)
            summary.mat <- cbind(summary.mat,temp$pvalues[1:outTopN])
        else
           summary.mat <- cbind(summary.mat,temp$J[1:outTopN])

    }
	colnames(summary.mat) <- col.name
    summary.mat = as.data.frame(summary.mat, stringsAsFactors = FALSE)
	return(summary.mat)
}

#-------------------------------------------------------------------------#
.sort2saveObj <- function(object,marker2Print,outTopN)
{
    markerNames <- colnames(object$J)
    col.name <- markerNames[marker2Print]
    gSnipID <-NULL
    Js <-NULL

    for(i in marker2Print )
    {
        # create summary data frames    
        temp <- object$J[,i]
        temp <- temp[order(-abs(temp))]

        gSnipID <- cbind(gSnipID, names(temp)[1:outTopN])
        Js <- cbind(Js,temp[1:outTopN])

    }
    colnames(Js) <- col.name
    colnames(gSnipID) <- col.name

    res <- NULL
    res$J <- as.matrix(Js)
    rownames(res$J) <- NULL
    res$gSnipID <- as.matrix(gSnipID)
    # set class for the result
    class(res) <- "jtGWAS"

    # add attribute for the result
    attr(res, 'outTopN') <- outTopN
    attr(res, 'standardized') <- attr(object, 'standardized')
    return(res)
}

#-------------------------------------------------------------------------#
# build summary matrix
.summary.mat <- function(object, marker2Print, SNP2Print, printP)
{
	markerNames <- colnames(object$J)
	p.values <- 2*pnorm(-abs(object$J))
	summary.mat <-NULL
	col.name <- NULL
	for(i in marker2Print )
	{
		# generate colname names for the dataframe
		col.name <- cbind(col.name, "SNPID")
		if(printP)
			col.name <- cbind(col.name, "P-value")
		else{
			if(attr(object,'standardized'))
				col.name <- cbind(col.name, "J*")
			else
				col.name <- cbind(col.name, "J")
		}

		# create summary data frames	
		summary.mat <- cbind(summary.mat,object$gSnipID[SNP2Print,i])
		if(printP)
			summary.mat <- cbind(summary.mat,p.values[SNP2Print,i])
		else
			summary.mat <- cbind(summary.mat,object$J[SNP2Print,i])
	} 
   
	colnames(summary.mat) <- col.name
	summary.mat = as.data.frame(summary.mat, stringsAsFactors = FALSE)
	return(summary.mat)
}

## format printing function for p-values
.printPvaluesFormat <- function(object, markerNames)
{
	
	# process the input
    summary.mat <-object
	col.name <- colnames(summary.mat)

	# printing parameters
	# width: maximum number of charactor allow to print in each line default to be 100
	# n.print.repeat: number of table to be cutted due to the width 
	
	n = floor(getOption("width")/20)
	width = n*20
    n.marker = ncol(summary.mat)/2
    n.print.repeat  <- ceiling(n.marker/n)
		
	# some basic element printint functions 
    marker.name.print <- function(x) cat(sprintf("%19s|",x,quote=F))
    space.print <- function(x) cat(" ")
    dash.print <- function(x) cat("-")
    double.dash.print <- function(x) cat("=")
    col.name.print <- function(i, col.name)
                    {
                        cat(sprintf("%11s", col.name[2*i-1]))
                        cat(sprintf("%8s", col.name[2*i]))
                        cat("|")
                    }
    SNPid.jstat.print <- function(j, i, summary.mat)
                    {
                        cat(sprintf("%11s", summary.mat[i, 2*j-1]))
                        cat(sprintf("%8.1e", as.numeric(summary.mat[i, 2*j])))
                        cat("|")
                    }

	# start of printing
    cat("\n\n")
    # Print out " ***Johckheere-Terpstra Test for Large Matrices****"
    title.string <- "Johckheere-Terpstra Test for Large Matrices\n"
    dummy <- sapply(1:(floor((width - nchar(title.string))/2)),space.print)
    cat(title.string)

    # Print out " *** Top Standardized Statistics and The Corresponding Variables IDs ***"
    subtitle.string <- "P-values for Top Standardized Statistics\n"
    dummy <- sapply(1:(ceiling((width - nchar(subtitle.string))/2)),space.print)
    cat(subtitle.string)

    # Print out "=======...="    
    dummy <- sapply(1:width,double.dash.print)
    cat("\n\n")
    for(n.repeat in 1:n.print.repeat)
    {
        start.marker = (n.repeat-1)*n + 1
        end.marker = min((n.repeat)*n, n.marker)

        # Print out "|  marker: x | marker: x| ...|"
    #    cat("|")
        dummy <- sapply(markerNames[start.marker:end.marker],marker.name.print)
        cat("\n")

        # Print out "------...-"    
        dummy <- sapply(1:((end.marker-start.marker+1)*20),dash.print)
        cat("\n")

        # Print out "|  ID   J* | SNPID   J*| ....|"
     #   cat("|")
        dummy <- sapply(start.marker:end.marker, col.name.print, col.name = col.name)
        cat("\n")

        # Print out "------...-"    
        dummy <- sapply(1:((end.marker-start.marker+1)*20),dash.print)
        cat("\n")

        # Print out "|"rnID"     "numeric value of J*" | ... |" 
        for(i in 1: nrow(summary.mat)){
     #       cat("|")
            dummy <- sapply(start.marker:end.marker, SNPid.jstat.print, i=i,summary.mat = summary.mat)
            cat("\n")
        }
        cat("\n\n")
    }
}

# print out statistics.

.printJStatisticsFormat <- function(summary.mat, markerNames)
{
	
	col.name <- colnames(summary.mat)
	
	# printing parameters
	# width: maximum number of charactor allow to print in each line default to be 100
	# n.print.repeat: number of table to be cutted due to the width 
	n = floor(getOption("width")/20)
	width = n*20
    n.marker = ncol(summary.mat)/2
    n.print.repeat  <- ceiling(n.marker/n)
   
	# some basic element printint functions  
	marker.name.print <- function(x) cat(sprintf("%19s|",x,quote=F))
    space.print <- function(x) cat(" ")
    dash.print <- function(x) cat("-")
    double.dash.print <- function(x) cat("=")
    col.name.print <- function(i, col.name)
                    {
                        cat(sprintf("%11s", col.name[2*i-1]))
                        cat(sprintf("%8s", col.name[2*i]))
                        cat("|")
                    }
    SNPid.jstat.print <- function(j, i, summary.mat)
                    {
                        cat(sprintf("%11s", summary.mat[i, 2*j-1]))
			            cat(sprintf("%8.1e", as.numeric(summary.mat[i, 2*j])))
						cat("|")
                    }

    cat("\n\n")
    # Print out " ***Johckheere-Terpstra Test for Large Matrices****"
    title.string <- "Johckheere-Terpstra Test for Large Matrices\n"
    dummy <- sapply(1:(floor((width - nchar(title.string))/2)),space.print)
    cat(title.string)

    # Print out " *** Top Standardized Statistics and The Corresponding Variables IDs ***"
    subtitle.string <- "Top Statistics\n"
    dummy <- sapply(1:(ceiling((width - nchar(subtitle.string))/2)),space.print)
    cat(subtitle.string)

    # Print out "=======...="    
    dummy <- sapply(1:width,double.dash.print)
    cat("\n\n")
    for(n.repeat in 1:n.print.repeat)
    {
        start.marker = (n.repeat-1)*n + 1
        end.marker = min((n.repeat)*n, n.marker)

        # Print out "|  marker: x | marker: x| ...|"
 #       cat("|")
        dummy <- sapply(markerNames[start.marker:end.marker],marker.name.print)
        cat("\n")

        # Print out "------...-"    
        dummy <- sapply(1:((end.marker-start.marker+1)*20),dash.print)
        cat("\n")

        # Print out "|  ID   J* | SNPID   J*| ....|"
  #      cat("|")
        dummy <- sapply(start.marker:end.marker, col.name.print, col.name = col.name)
        cat("\n")

        # Print out "------...-"    
        dummy <- sapply(1:((end.marker-start.marker+1)*20),dash.print)
        cat("\n")

        # Print out "|"rnID"     "numeric value of J*" | ... |" 
        for(i in 1:nrow(summary.mat)){
 #           cat("|")
            dummy <- sapply(start.marker:end.marker, SNPid.jstat.print, i=i,summary.mat = summary.mat)
            cat("\n")
        }
        cat("\n\n")
    }
}
# format printing function for standardized statistics 

.printJstarFormat <- function(object, markerNames)
{
	
	# process input
    summary.mat <-object
	col.name <- colnames(summary.mat)
	
	# printing parameters
	# width: maximum number of charactor allow to print in each line default to be 100
	# n.print.repeat: number of table to be cutted due to the width 
	n = floor(getOption("width")/20)
	width = n*20
    n.marker = ncol(summary.mat)/2
    n.print.repeat  <- ceiling(n.marker/n)
   
	# some basic element printint functions  
	marker.name.print <- function(x) cat(sprintf("%19s|",x,quote=F))
    space.print <- function(x) cat(" ")
    dash.print <- function(x) cat("-")
    double.dash.print <- function(x) cat("=")
    col.name.print <- function(i, col.name)
                    {
                        cat(sprintf("%11s", col.name[2*i-1]))
                        cat(sprintf("%8s", col.name[2*i]))
                        cat("|")
                    }
    SNPid.jstat.print <- function(j, i, summary.mat)
                    {
                        cat(sprintf("%11s", summary.mat[i, 2*j-1]))
			            cat(sprintf("%8.3f", as.numeric(summary.mat[i, 2*j])))
						cat("|")
                    }

    cat("\n\n")
    # Print out " ***Johckheere-Terpstra Test for Large Matrices****"
    title.string <- "Johckheere-Terpstra Test for Large Matrices\n"
    dummy <- sapply(1:(floor((width - nchar(title.string))/2)),space.print)
    cat(title.string)

    # Print out " *** Top Standardized Statistics and The Corresponding Variables IDs ***"
    subtitle.string <- "Top Standardized Statistics\n"
    dummy <- sapply(1:(ceiling((width - nchar(subtitle.string))/2)),space.print)
    cat(subtitle.string)

    # Print out "=======...="    
    dummy <- sapply(1:width,double.dash.print)
    cat("\n\n")
    for(n.repeat in 1:n.print.repeat)
    {
        start.marker = (n.repeat-1)*n + 1
        end.marker = min((n.repeat)*n, n.marker)

        # Print out "|  marker: x | marker: x| ...|"
 #       cat("|")
        dummy <- sapply(markerNames[start.marker:end.marker],marker.name.print)
        cat("\n")

        # Print out "------...-"    
        dummy <- sapply(1:((end.marker-start.marker+1)*20),dash.print)
        cat("\n")

        # Print out "|  ID   J* | SNPID   J*| ....|"
  #      cat("|")
        dummy <- sapply(start.marker:end.marker, col.name.print, col.name = col.name)
        cat("\n")

        # Print out "------...-"    
        dummy <- sapply(1:((end.marker-start.marker+1)*20),dash.print)
        cat("\n")

        # Print out "|"rnID"     "numeric value of J*" | ... |" 
        for(i in 1:nrow(summary.mat)){
   #         cat("|")
            dummy <- sapply(start.marker:end.marker, SNPid.jstat.print, i=i,summary.mat = summary.mat)
            cat("\n")
        }
        cat("\n\n")
    }
}


# function to check valid the requested marker2Print
.markerRangeCheck <- function(object,marker2Print)	
{
	# get the colnames and rownames 
	markerNames <- colnames(object$J)	
	
	# for the NA cases, which means full print out the results.
	if(is.na(marker2Print[1]))
		marker2Print = seq(1:ncol(object$J))

	# out of range
	if(is.character(marker2Print))
	{
		# if the user specify an wired names for the marker, we will tell
		# the user we ignore it.
		if(length(marker2Print[!(marker2Print %in% markerNames)]))
		{	
			cat("\nWARNING: specified marker ID/range: ")
			cat(marker2Print[!(marker2Print %in% markerNames)]) 
			cat(" are not found, \nresults for them are omitted.\n")
		}
		marker2Print <- intersect(marker2Print, markerNames)
	}
	
	if(is.numeric(marker2Print))
	{
		# if the user specify an wired numbers for the marker, we will tell
        # the user we ignore it.
		
        if(length(marker2Print[!(marker2Print %in% seq( 1: ncol(object$J)) )]))
        {   
			if( length(marker2Print) == 10&& marker2Print[length(marker2Print)] == 10)
			{}
			else
			{	
				cat("\nWARNING: specified marker ID/range: ")
            	cat(marker2Print[!(marker2Print %in% seq( 1: ncol(object$J)))])
            	cat(" are not found, \nresults for them are omitted.\n")
			}
        }
		marker2Print <- intersect(marker2Print, seq( 1: ncol(object$J)) )
	}
	return(marker2Print)	
}

# check the validaty of the requested SNP range
.SNPRangeCheck <- function(object, SNP2Print)
{	
    if(is.na(SNP2Print[1]))
        SNP2Print <- seq(1:nrow(object$J))

	if(is.numeric(SNP2Print))
    {
        # if the user specify an wired range for the SNP, we will tell
        # the user we ignore it.
        if(length(SNP2Print[!(SNP2Print %in% seq( 1: nrow(object$J)) )]))
        {
            cat("WARNING: specified SNP ID/range: ")
            cat(SNP2Print[!(SNP2Print %in% seq( 1: nrow(object$J)))])
            cat(" are not found, results for them are omitted.\n")
        }
        SNP2Print <- intersect(SNP2Print, seq(1: nrow(object$J)) )
    }
	return(SNP2Print)
}


# warning and tile print function

.note.title.print <- function(object, outTopN, printP, marker.is.default, SNP.is.default.for.top.hit)
{
	# print out note
	if(marker.is.default && SNP.is.default.for.top.hit)
	{
		if(is.na(attr(object,'outTopN'))&&is.na(outTopN))
		{
			cat("\n\nNote: By default, only part of the results may be printed!\n")
			cat("Please specify names/ranges using 'marker2Print' and 'SNP2Print',\n")
			cat("or set to 'NA' to print the full results.\n\n")
		}
		else
		{
			cat("\n\nNote: By default, all requested top hits are printed for up to 10 markers.\n")
			cat("Please specify names/ranges using 'marker2Print' and 'SNP2Print',\n")
			cat("or set to 'NA' to print the full results.\n\n")
			
		}
	}
	
	# print out titles
	if(is.na(attr(object,'outTopN'))&&is.na(outTopN))
	{           
		if(printP){
         	title.String <- "\n      Johckheere-Terpstra Test \n    P-values Based on Standardized Statistics\n\n"
     	}else{
			if(attr(object,'standardized'))
				title.String <- "\n    Johckheere-Terpstra Test Standardized Statistics\n\n"
			else
				title.String <- "\n    Johckheere-Terpstra Test Statistics\n\n"
     	}
     	cat(title.String)
	}

}




















