\name{jtGWAS}
\alias{jtGWAS}
\docType{package}
\title{
  Compute the Jonckheere-Terpstra Test Statistics
}
\description{
  A method to compute the Jonckheere-Terpstra test statistics for large numbers of dependent and independent variables, with optional multi-threaded execution. The calculation of the standardized test statistic employs the null variance equation as defined by Hollander and Wolfe (1999, eq. 6.19) to account for ties in the data.
}

\usage{
  jtGWAS(X, G, outTopN=15L, numThreads=1L, standardized=TRUE)  	
}

\arguments{
  \item{X}{A matrix of dependent variables, e.g. marker levels. Row names are sample IDs, and column names are marker IDs. Required.
}
  \item{G}{A matrix of independent variables, e.g. SNP counts. Row names are sample IDs, and column names are SNP IDs. Required.
}
  \item{outTopN}{An integer to indicate the number of top hits to be reported for each marker, based on the standardized Jonckheere-Terpstra test statistics. Optional, and the default value is \code{15L}. If set to \code{NA}, all results are returned.}
  \item{numThreads}{A integer to indicate the number of threads used in the computation. Optional, and the default value is \code{1L}(sequential computation).
}
  \item{standardized}{A boolean to specify whether to return standardized statistics (\code{TRUE}) or non-standardized statistics (\code{FALSE}). Optional, the default value is \code{TRUE}.
} 
}

\note{
  Rows (samples) are assumed to be in the same order in \code{X} and \code{G}.
}

\references{
    Hollander, M. and Wolfe, D. A. (1999) \emph{Nonparametric Statistical Methods}. New York: Wiley, 2nd edition.
}

\value{
  A \code{list} with two objects
  \item{J}{A matrix of the standardized/non-standardized Jonckheere-Terpstra test statistics,
depending on the value of the \code{standardized} argument.
}
  \item{gSnipID}{If \code{outTopN} was specified, this object is a matrix of the column IDs of \code{G} corresponding to the top standardized Jonckheere-Terpstra test statistics for each marker. Otherwise this is a vector of column IDs of \code{G}.
}
}


\examples{
# Generate dummy data	
num_patient <- 100
num_marker 	<- 10
num_SNP     <- 500
set.seed(12345)
X_pat_mak <- matrix(rnorm(num_patient*num_marker), num_patient, num_marker)
G_pat_SNP <- matrix(rbinom(num_patient*num_SNP, 2, 0.5), num_patient, num_SNP)
colnames(X_pat_mak) <- colnames(X_pat_mak, do.NULL=FALSE, prefix="Mrk:")
colnames(G_pat_SNP) <- colnames(G_pat_SNP, do.NULL=FALSE, prefix="SNP:")

res <- jtGWAS(X_pat_mak, G_pat_SNP, outTopN=5)
res
res <- jtGWAS(X_pat_mak, G_pat_SNP, outTopN=NA)
head(res) 
}

