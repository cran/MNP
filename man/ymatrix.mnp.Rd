\name{ymatrix.mnp}

\alias{ymatrix.mnp}

\title{ymatrix.mnp: Generating the response matrix for the Multinomial
  Probit Models}

\description{
  \code{ymatrix.mnp.R} is a function used to generate the response matrix
  for the Multinomial Probit models.
}

\usage{
  ymatrix.mnp(data, base=NULL, extra=FALSE, verbose=verbose)
}

\arguments{
  \item{data}{A model frame.}
  \item{base}{The name of the base category. The default is the lowest
    level of the response variable. For the Multinomial ordered Probit
    model, the default base category will be the last column in the
    response variables. 
  }
  \item{extra}{
  logical. \code{TRUE} for generating the names for the coefficients as
  well as the model matrix. The default is \code{FALSE}.  
  }
  \item{verbose}{logical. If \code{TRUE}, helpful messages
    are printed on the screen. The default is \code{FALSE}.
  }
}

\value{
  If \code{extra=TRUE}, a list  containing the following elements:
  \item{Y}{The response matrix.}
  \item{MoP}{logical. \code{TRUE} for the Multinomial ordered Probit
    model, and \code{FALSE} for the standard Multinomial Probit model.}
  \item{lev}{The factor level of the response variable.}
  \item{p}{The total number of alternatives.}
  \item{base}{The base category name.}
  If \code{extra=FALSE}, the response matrix will be returned.
}

\seealso{\code{mprobit}, \code{xmatrix.mnp}}

\author{
  Kosuke Imai, Department of Politics, Princeton University
  \email{kimai@Princeton.Edu}
}

\keyword{misc}
