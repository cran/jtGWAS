\name{jtGWAS-package}
\alias{jtGWAS-package}
\alias{jtGWAS-package}
\docType{package}
\title{
	Efficient Jonckheere-Terpstra Test Statistics
}
\description{The core of this his \code{Rcpp}-based package is a function to compute standardized Jonckheere-Terpstra test statistics for large numbers of dependent and independent variables, e.g. genome-wide analysis. It implements \code{OpenMP}, allowing the option of computing on multiple threads. Supporting functions are also provided to calculate p-values and summarize results.
}
\details{
\tabular{ll}{
Package: \tab jtGWAS\cr
Type: \tab Package\cr
Version: \tab 1.5.1\cr
Date: \tab 2017-08-14\cr
License: \tab GPL-3\cr
}
Please see the example function calls below, or refer to the individual function documentation or the included vignette for more information. The package vignette serves as a tutorial for using this package. The technical details are provided in the reference cited below. Specifically, the calculation of the standardized test statistic employs the null variance equation as defined by Hollander and Wolfe (1999, eq. 6.19) to account for ties in the data.
}

\author{
Jiaxing Lin, Alexander Sibley, Ivo Shterev, and Kouros Owzar

Maintainer: Jiaxing Lin <jiaxing.lin@duke.edu>
}

\references{
    Hollander, M. and Wolfe, D. A. (1999) \emph{Nonparametric Statistical Methods}. New York: Wiley, 2nd edition.
}

\keyword{ package }
\seealso{
    \code{\link{Rcpp}}
}

\examples{
# Generate dummy data	
num_patient <- 100
num_marker  <- 10
num_SNP     <- 500
set.seed(12345)
X_pat_mak <- matrix(rnorm(num_patient*num_marker), num_patient, num_marker)
G_pat_SNP <- matrix(rbinom(num_patient*num_SNP, 2, 0.5), num_patient, num_SNP)
colnames(X_pat_mak) <- colnames(X_pat_mak, do.NULL=FALSE, prefix="Mrk:")
colnames(G_pat_SNP) <- colnames(G_pat_SNP, do.NULL=FALSE, prefix="SNP:")
   
res <- jtGWAS(X_pat_mak, G_pat_SNP, outTopN=15)
res
}
