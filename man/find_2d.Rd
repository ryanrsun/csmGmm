% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_incongruous.R
\name{find_2d}
\alias{find_2d}
\title{Tells if row x if allTestStats is an incongruous result (has a higher lfdr than a set of
test statistics with lower magnitudes). For K=2 case.}
\usage{
find_2d(x, allTestStats)
}
\arguments{
\item{x}{Scalar, which row of allTestStats to check.}

\item{allTestStats}{J*K vector of all test statistics.}
}
\value{
A scalar denoting the number of sets with lower lfdr and test statistics of lower magnitude. 0 means congruous result.
}
\description{
Tells if row x if allTestStats is an incongruous result (has a higher lfdr than a set of
test statistics with lower magnitudes). For K=2 case.
}
\examples{
zMatrix <- cbind(rnorm(10^5), rnorm(10^5))
find_2d(x = 5, allTestStats = allTestStats)

}
