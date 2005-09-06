\name{predict.mnp}

\alias{predict.mnp}

\title{Posterior Prediction under the Bayesian Multinomial Probit Models}

\description{
  Obtains posterior predictions under a fitted (Bayesian) multinomial
  probit model. \code{predict} method for class \code{mnp}.
}

\usage{
  \method{predict}{mnp}(object, newdata = NULL, newdraw = NULL,
          type = c("prob", "choice", "order", "latent"), verbose = FALSE, ...)
}

\arguments{
  \item{object}{An output object from \code{mnp}.}
  \item{newdata}{An optional data frame containing the values of the
    predictor variables. Predictions for multiple values of the
    predictor variables can be made simultaneously if \code{newdata} has
    multiple rows. The default is the original data frame used
    for fitting the model.
  }
  \item{newdraw}{An optional matrix of MCMC draws to be used for
    posterior predictions. The default is the original MCMC draws stored
    in \code{object}.
  }
  \item{type}{The type of posterior predictions required. There are
    four options:
    \code{type = "prob"} returns the predictive probabilities of being
    the most preferred choice among the choice set.
    \code{type = "choice"} returns the Monte Carlo sample of the 
    most preferred choice,
    \code{type = "order"} returns the Monte Carlo sample of the
    ordered preferences,
    and \code{type = "latent"} returns the Monte Carlo sample of the
    predictive values of the latent variable. The default is to return
    all four types of posterior predictions.
  }
  \item{verbose}{logical. If \code{TRUE}, helpful messages along with a
    progress report on the Monte Carlo sampling from the posterior 
    predictive distributions are printed on the screen. The default is
    \code{FALSE}.
  }
  \item{...}{further arguments passed to or from other methods.}
}

\details{The posterior predictive values are computed using the
  Monte Carlo sample stored in the \code{mnp} output (or other sample if
  \code{newdraw} is specified). Given each Monte Carlo sample of the
  parameters and each vector of predictor variables, we sample the
  vector-valued latent variable from the appropriate multivariate Normal
  distribution. Then, using the sampled predictive values of the latent
  variable, we construct the most preferred choice as well as the
  ordered preferences. Averaging over the Monte Carlo sample
  of the preferred choice, we obtain the predictive probabilities of
  each choice being most preferred given the values of the predictor 
  variables.  Since the predictive values are computed via Monte Carlo
  simulations, each run may produce somewhat different values. The
  computation may be slow if predictions with many values of the
  predictor variables are required and/or if a large Monte Carlo
  sample of the model parameters is used. In either case, setting
  \code{verbose = TRUE} may be helpful in monitoring the progress of the
  code.
}

\value{
  \code{predict.mnp} yields a list containing at least one of the
  following elements:
  \item{o}{A three dimensional array of the Monte Carlo sample from the
    posterior predictive distribution of the ordered preferences. The
    first dimension corresponds to the alternatives in the choice set,
    the second dimension corresponds to the rows of \code{newdata} (or
    the original data set if \code{newdata} is left unspecified), and
    the third dimension indexes the Monte Carlo sample.}
  \item{p}{A matrix of the posterior predictive probabilities for each
  alternative in the choice set being most preferred. The
  rows correspond to the rows of \code{newdata} (or the original data
  set if \code{newdata} is left unspecified) and the columns correspond
  to the alternatives in the choice set.
  }
  \item{y}{A matrix of the Monte Carlo sample from the posterior predictive
   distribution of the most preferred choice. The rows
   correspond to the rows of \code{newdata} (or the original data set if
   \code{newdata} is left unspecified) and the columns index the Monte
   Carlo sample.
  }
  \item{w}{A three dimensional array of the Monte Carlo sample from the 
    posterior predictive distribution of the latent
    variable. The first dimension corresponds to the alternatives in the
     choice set, the second dimension corresponds to the rows of
     \code{newdata} (or the original data set if \code{newdata} is left
     unspecified), and the third dimension indexes the Monte Carlo sample.
  }
}

\seealso{\code{mnp}; MNP home page at
  \url{http://imai.princeton.edu/research/MNP.html}}

\author{
  Kosuke Imai, Department of Politics, Princeton University
  \email{kimai@Princeton.Edu}
}

\keyword{methods}
