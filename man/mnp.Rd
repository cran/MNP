\name{mnp}

\alias{mnp}
\alias{MNP}

\title{Fitting the Multinomial Probit Models via Markov chain Monte Carlo} 

\description{
  \code{mnp} is used to fit the (Bayesian) Multinomial Probit
  models via Markov chain Monte Carlo. Along with the standard
  Multinomial Probit model, it can also fit models with different choice
  sets for each observation, and complete or partial ordering of all the
  available alternatives. The computation uses the efficient marginal
  data augmentation algorithm that is developed by Imai and van Dyk (2004).
}

\usage{
mnp(formula, data = parent.frame(), choiceX = NULL, cXnames = NULL,
        base = NULL, n.draws = 5000, p.var = "Inf",  p.df = n.dim+1, p.scale = 1,
        coef.start = 0, cov.start = 1, burnin = 0, thin = 0, verbose = FALSE) 
}

\arguments{
  \item{formula}{A symbolic description of the model to be fit
    specifying the response variable and covariates. The formula should 
    not include the choice-specific covariates. Details and
    specific examples are given below. 
  }
  \item{data}{An optional data frame in which to interpret the variables
    in \code{formula} and \code{choiceX}. The default is the environment
    in which \code{mnp} is called.
  }
  \item{choiceX}{An optional list containing a matrix of choice-specific
    covariates for each category. Details and examples are provided
    below.
  }
  \item{cXnames}{A vector of the names for the choice-specific covariates
    specified in \code{choiceX}. The details and examples are provided
    below.
  }
  \item{base}{The name of the base category. For the standard
    Multinomial Probit model, the default is the lowest level of the
    response variable. For the Multinomial ordered Probit model, the
    default base category will be the last column in the matrix of
    response variables.
  }
  \item{n.draws}{A positive integer. The number of MCMC draws. The
    default is \code{5000}.
  }
  \item{p.var}{A positive definite matrix. The prior variance of the
    coefficients.  A scalar input can set the prior variance to the
    diagonal matrix whose diagonal element is equal to that value. The
    default is \code{"Inf"}, which represents an improper noninformative
    prior on the coefficients.   
  }
  \item{p.df}{A positive integer greater than \code{n.dim-1}. The prior
  degree of freedom parameter for the covariance matrix. The default is
    \code{n.dim+1}, which is equal to the total number of alternatives.
  }
  \item{p.scale}{A positive definite matrix whose first diagonal element
    is set to \code{1}. The prior scale matrix for the covariance
    matrix. The first diagonal element will be set to 1 if
    it is not equal to 1 already. A scalar input can be used to set the
    scale matrix to a diagonal matrix with diagonal elements equal to
    the scalar input value, except that the first diagonal element will
    be set to one. The default is \code{1}.
  }
  \item{coef.start}{A vector. The starting values for the coefficients.
    A scalar input will set the starting values for all the coefficients
    equal to that value.  The default is \code{0}.
  }
  \item{cov.start}{A positive definite matrix whose first diagonal
    element is set to \code{1}. The starting values for the covariance
    matrix. The first diagonal element will be set to 1 if it is not
    equal to 1 already. A scalar input can be used to set the starting
    value to a diagonal matrix with diagonal elements equal to the scalar
    input value, except that the first diagonal element will be set to
    one. The default is \code{1}. 
  }
  \item{burnin}{A positive integer. The burnin interval for the Markov
    chain; i.e., the number of initial Gibbs draws that should not be
    stored. The default is \code{0}. 
  }
  \item{thin}{A positive integer. The thinning interval for the Markov
    chain; i.e., the number of Gibbs draws between the recorded values
    that are skipped. The default is \code{0}. 
  }
  \item{verbose}{logical. If \code{TRUE}, helpful messages along with
    a progress report (every 10\%) of the Gibbs sampling are printed
    on the screen. The default is \code{FALSE}.
  }
}

\details{
  For \strong{the standard Multinomial Probit model} where only the most
  preferred choice is observed, use the syntax, \code{y ~ x1 +
    x2}, where \code{y} is a factor variable indicating the most
  preferred choice and \code{x1} and \code{x2} are individual-specific
  covariates. The interactions of individual-specific variables with each
  of the choice indicator variables will be fit.
  
  To specify \strong{choice specific covariates}, use the syntax,
  \code{choiceX=list(A=cbind(z1, z2), B=cbind(z3, z4), C=cbind(z5,
  z6))}, where \code{A}, \code{B}, and \code{C} represent the choice
  names of the response variable, and \code{z1} and \code{z2} are each
  vectors of length \eqn{n} that record the values of the two
  choice-specific covariates for each individual for choice A, likewise
  for \code{z3}, ..., \code{z6}. The corresponding variable names via
  \code{cXnames=c("price", "quantity")}
  need to be specified, where \code{price} refers to the coefficient
  name for \code{z1}, \code{z3}, and \code{z5}, and \code{quantity}
  refers to that for \code{z2}, \code{z4}, and \code{z6}.

  If \strong{the choice set} varies from one observation to another, use the
  syntax, \code{cbind(y1, y2, y3) ~ x1 + x2}, in the case of a
  three choice problem, and indicate unavailable alternatives by
  \code{NA}. If only the most preferred choice is observed, y1, y2, and y3 are
   indicator variables that take on the value one for individuals who
   prefer that choice and zero otherwise. The last column of the
   response matrix, \code{y3} in this particular example syntax, is
   used as the base category.
  
  For \strong{the Multinomial ordered Probit model} where the complete
  or partial ordering of the available alternatives is recorded, use the
  same syntax as when the choice set varies (i.e., \code{cbind(y1, y2,
  y3, y4) ~ x1 + x2}). For each observation, all the available
  alternatives in the response variables should be numerically ordered
  in terms of preferences such as \code{1 2 2 3}. Ties are allowed. The
  missing values in the response variable should be denoted by
  \code{NA}. The software will impute these missing values using the
  specified covariates. The resulting uncertainty estimates of the
  parameters will properly reflect the amount of missing data. For
  example, we expect the standard errors to be larger when there is more
  missing data.
}

\examples{
## load the detergent data
data(detergent)
## run the standard Multinomial Probit model with intercepts and the price
res1 <- mnp(choice ~ 1, choiceX=list(Surf=Surf, Tide=Tide, Wisk=Wisk,
                                     EraPlus=EraPlus, Solo=Solo, All=All),
                cXnames=c("price"), data=detergent, n.draws=5000,
                p.df=6, burnin=1000, thin=3, verbose=TRUE) 
summary(res1)

## load the Japanese election data
data(japan)
## run the Multinomial ordered Probit model
res2 <- mnp(cbind(LDP, NFP, SKG, JCP) ~ sex + education + age,
            data = japan, n.draws=5000, verbose = TRUE)
summary(res2)
}

\value{
  An object of class \code{mnp} containing the following elements:
  \item{param}{A matrix of the Gibbs draws for each parameter; i.e., the
    coefficients and covariance matrix. For the covariance matrix, the
    elements on or above the diagonal are returned.
  }
  \item{call}{The matched call.}
  \item{x}{The matrix of covariates.}
  \item{y}{The vector or matrix of the response variable.}
  \item{n.alt}{The total number of alternatives.}
  \item{p.var}{The prior variance for the coefficients.}
  \item{p.df}{The prior degrees of freedom parameter for the covariance matrix.}
  \item{p.scale}{The prior scale matrix for the covariance matrix.}
  \item{burnin}{The number of initial burnin draws.}
  \item{thin}{The thinning interval.}
  }
}

\references{
  Imai, Kosuke and David A. van Dyk. (2004) \dQuote{A Bayesian Analysis of the
  Multinomial Probit Model Using the Marginal Data Augmentation,}
  \emph{Journal of Econometrics},
  Forthcoming. \url{http://www.princeton.edu/~kimai/research/mnp.html} 
}

\author{
  Kosuke Imai, Department of Politics, Princeton University
  \email{kimai@Princeton.Edu}, \url{http://www.princeton.edu/~kimai};
  Jordan R. Vance, Princeton University; David A. van Dyk, Department of
  Statistics, University of California, Irvine \email{dvd@uci.edu},
  \url{http://www.ics.uci.edu/~dvd}.  
}

\seealso{\code{summary.mnp}; MNP home page at
  \url{http://www.princeton.edu/~kimai/research/MNP.html}}

\keyword{models}





