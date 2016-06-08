\name{summary.jtGWAS}
\alias{summary.jtGWAS}
\docType{package}
\title{
  Summarize Jonckheere-Terpstra Test Statistics and P-Values
}
\description{
  Summary method for jtGWAS results.
}

\usage{
\method{summary}{jtGWAS}(object, marker2Print=1:10, SNP2Print=1:10, 
				 printP=TRUE, outTopN=NA, subObj=FALSE, \dots)	
}

\arguments{
  \item{object}{A jtGWAS object that is the return of method \code{jtGWAS}. Required.
}
  \item{marker2Print}{Either a numeric or character vector that indicates the desired dependent variables to print. The default is 1:10. Set to \code{NA} to print full results. Optional.
}
  \item{SNP2Print}{Either a numeric or character vector that indicates the desired independent variables to print. If \code{outTopN=NA} in function \code{jtGWAS}, i.e., the results are not sorted, both numeric and character vectors can be used to set the print range. The default range is 1:10. Set to \code{NA} to print full results. If \code{outTopN!=NA} in function \code{jtGWAS}, the range of \code{SNP2Print} refers to the range of top normalized statistics computed in \code{jtGWAS} and only numeric vectors are processed. Optional.
}
  \item{printP}{A boolean indicating whether to print the p-value (\code{TRUE}) or the standardized statistics (\code{FALSE}). The default value is \code{TRUE}.
}
  \item{outTopN}{An integer specifying the number of top hits to print in the summary, if the statistics were not sorted during the statistics evaluation (i.e., if \code{outTopN=NA} in function \code{jtGWAS}). Optional, and the default value is \code{NA}.
}
  \item{subObj}{A boolean indicating whether to return a \code{jtGWAS} object subset per the requested summary (\code{TRUE}). Optional. The default value is \code{FALSE} (nothing is returned).
}
  \item{...}{Additional arguments affecting the summary produced.
}
}

\value{
  If \code{subObj=TRUE}, this method returns a \code{jtGWAS} object matching the subset of statistics from \code{object} that are being printed. For example, if \code{object} is not sorted by top hits, \code{summary(object, outTopN=10, subObj=TRUE)} will print the summary and return a subset of \code{object} that contains only the top 10 SNPs for each marker. If \code{subObj=FALSE}, nothing is returned.
}
\note{
  This function prints a matrix or paired columns of SNP IDs and statistics/p-values to the log.
}
\seealso{
  \code{\link{jtGWAS}}, \code{\link{pvalues}}
}

\examples{
# Generate dummy data	
num_patient <- 100
num_marker 	<- 10
num_SNP     <- 500
set.seed(12345);
X_pat_mak <- matrix(rnorm(num_patient*num_marker), num_patient, num_marker)
G_pat_SNP <- matrix(rbinom(num_patient*num_SNP, 2, 0.5), num_patient, num_SNP)
colnames(X_pat_mak) <- colnames(X_pat_mak, do.NULL=FALSE, prefix="Mrk:")
colnames(G_pat_SNP) <- colnames(G_pat_SNP, do.NULL=FALSE, prefix="SNP:")

res <- jtGWAS(X_pat_mak, G_pat_SNP, outTopN=5)
summary(res, marker2Print=c("Mrk:1","Mrk:2"), SNP2Print=1:5, printP=FALSE)
summary(res, marker2Print=NA, SNP2Print=1:5, printP=FALSE)

res <- jtGWAS(X_pat_mak, G_pat_SNP, outTopN=NA)
summary(res, marker2Print=1:10, SNP2Print=1:10, printP=TRUE)
summary(res, marker2Print=c("Mrk:1","Mrk:2"), SNP2Print=c("SNP:1","SNP:2"), printP=TRUE) 

res <- jtGWAS(X_pat_mak, G_pat_SNP, outTopN=NA, standardized=FALSE)
summary(res, outTopN=10)
}

