\name{getcmat}
\alias{getcmat}
\title{A matrix of point process configurations in p-dimensional covariate space}
\description{
This function returns a matrix where the rows correspond to each possible
point process specification in p-dimensional covariate space; this may be
helpful in specifying the argument \code{settozero} in the functions \code{\link{monoreg}} and \code{\link{monosurv}}.
}
\usage{
getcmat(p)
}
\arguments{
\item{p}{Number of covariate axes.}
}
\value{
A zero-one matrix with 2^p - 1 rows and p columns.
}
\author{Olli Saarela <olli.saarela@utoronto.ca>}
