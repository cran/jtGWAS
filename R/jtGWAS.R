jtGWAS <- function(X, G, outTopN = 15L, numThreads = 1L, standardized = TRUE) {
	markerNames <- colnames(X)
	SNPNames <- colnames(G)

	# provide default col and row names for the data frame if not provided.
    if (is.null(markerNames))
        markerNames <- paste("Mrk:", 1:ncol(X), sep="")
    if (is.null(SNPNames))
        SNPNames <- paste("SNP:", 1:ncol(G), sep="")
	
	if(numThreads == 1) # calling single treaded function
	{	
	# compute the statistics
		if(is.na(outTopN))
    	{
			res <- .Call('jtGWAS_jtGWAS', PACKAGE = 'jtGWAS', X, G, !is.na(outTopN), 15, standardized)
			rownames(res$J) <- SNPNames
			res$gSnipID <- SNPNames
    	}
		else
		{	
			res<-.Call('jtGWAS_jtGWAS', PACKAGE = 'jtGWAS', X, G, !is.na(outTopN), outTopN, standardized)
		
			# using the SNP names instead of there matrix index for the result
			temp <- (res$J)	
			rows <- min(outTopN,ncol(G))
			for(i in 1:rows)
			{		for(j in 1:ncol(X))
					{temp[i,j] <- SNPNames[res$gSnipID[i,j]+1]}
			}
			res$gSnipID <-temp
			colnames(res$gSnipID) <- markerNames
		}
	}
	
	if(numThreads != 1) # calling parallel version function
	{	
		# compute the statistics
		if(is.na(outTopN))
    	{
			res<-.Call('jtGWAS_jtGWASmp', PACKAGE = 'jtGWAS', X, G, !is.na(outTopN), numThreads, 15, standardized)
        	rownames(res$J) <- SNPNames
			res$gSnipID <- SNPNames
		}
		else
		{	
			res<-.Call('jtGWAS_jtGWASmp', PACKAGE = 'jtGWAS', X, G, !is.na(outTopN), numThreads, outTopN, standardized)
			# using the SNP names instead of there matrix index for the result
			temp_mp <- (res$J)
			rowsmp <- min(outTopN,ncol(G))	
			for(i in 1: rowsmp)
			{	for(j in 1:ncol(X))
					{temp_mp[i,j] <- SNPNames[res$gSnipID[i,j]+1]}
			}
			res$gSnipID <-temp_mp
			colnames(res$gSnipID) <- markerNames
		}
	}	
	# set the col and row names for the results matrix.
	colnames(res$J) <- markerNames
	
	# set class for the result
	class(res) <- "jtGWAS"	
	
	# add attribute for the result
	attr(res, 'outTopN') <- outTopN
	attr(res, 'standardized') <- standardized	
	return(res)
}




