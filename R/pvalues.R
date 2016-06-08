pvalues <- function(jtGWAS.object)
{
	if(attr(jtGWAS.object,'standardized'))
	 	return(2*pnorm(-abs(jtGWAS.object$J)))
	else
		 stop("This function is only to be used with standardized results.")
}

