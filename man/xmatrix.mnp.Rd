\name{xmatrix.mnp}

\alias{xmatrix.mnp}

\title{xmatrix.mnp: Generating the model matrix for the Multinomial
  Probit Models}

\description{
  \code{xmatrix.mnp.R} is a function used to generate the model matrix
  for the Multinomial Probit models.
}

\usage{
  xmatrix.mnp(formula, data = sys.parent(), choiceX=NULL,
              cXnames=NULL, base=NULL, n.dim, lev,
              MoP=FALSE, verbose=FALSE, extra=FALSE)
}

\arguments{
  \item{formula}{A symbolic description of the model to be fit
    specifying the response variable and covariates. The formula should 
    not include the choice-specific covariates. Details and
    specific examples are given below. 
  }
  \item{data}{An optional data frame in which to interpret the variables
    in \code{formula} and \code{choiceX}. The default is the environment
    in which \code{mprobit} is called.
  }
  \item{choiceX}{An optional list containing a matrix of choice-specific
    covariates for each category. Details and examples are provided
    below.
  }
  \item{cXnames}{A vector of the names for the choice-specific covariates
    specified in \code{choiceX}. The details and examples are provided
    below.
  }
  \item{base}{
    Character. The base category name.
  }
  \item{n.dim}{
    A positive integer. The number of dimensions, i.e., the total number
    of alternatives minus one.
  }
  \item{lev}{
    The factor level for the response variable.
  }
  \item{MoP}{
    logical. \code{TRUE} for the Multinomial ordered Probit model,
    and \code{FALSE} for the standard Multinomial Probit model. The
    default is \code{FALSE}.
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
  \item{X}{The model matrix.}
  \item{coefnames}{A vector of the names for the coefficients.}
  If \code{extra=FALSE}, the model matrix will be returned.
}

\seealso{\code{mprobit}, \code{ymatrix.mnp}}

\author{
  Kosuke Imai, Department of Politics, Princeton University
  \email{kimai@Princeton.Edu}
}

\keyword{misc}
