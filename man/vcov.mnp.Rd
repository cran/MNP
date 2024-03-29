% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcov.mnp.R
\name{vcov.mnp}
\alias{vcov.mnp}
\title{Extract Multinomial Probit Model Covariance Matrix}
\usage{
\method{vcov}{mnp}(object, subset = NULL, ...)
}
\arguments{
\item{object}{An output object from \code{mnp}.}

\item{subset}{A scalar or a numerical vector specifying the row number(s) of
\code{param} in the output object from \code{mnp}. If specified, the
posterior draws of covariance matrix for those rows are extracted. The
default is \code{NULL} where all the posterior draws are extracted.}

\item{...}{further arguments passed to or from other methods.}
}
\value{
When a numerical vector or \code{NULL} is specified for
\code{subset} argument, \code{vcov.mnp} returns a three dimensional array
where the third dimension indexes posterior draws. When a scalar is
specified for \code{subset} arugment, \code{vcov.mnp} returns a matrix.
}
\description{
\code{vcov.mnp} is a function which extracts the posterior draws of
covariance matrix from objects returned by \code{mnp}.
}
\seealso{
\code{mnp}, \code{coef.mnp};
}
\author{
Kosuke Imai, Department of Government and Department of Statistics, Harvard University
\email{imai@harvard.edu}
}
\keyword{methods}
