\name{predict.mnp}

\alias{predict.mnp}

\title{Predict Method for the Multinomial Probit Models}

\description{
  Obtains predictions from a
  fitted  multinomial  probit model  object.  \code{predict} method  for
  class \code{mnp}.
}

\usage{
  \method{predict}{mnp}(object, newdata = NULL, type = c("prob", "choice", "latent"),
          verbose = FALSE, ...)
}

\arguments{
  \item{object}{An output object from \code{mnp}.}
  \item{newdata}{Optionally, a data frame in which to look for variables
    with which to predict. The default is the original data frame used
    for fitting the model.
  }
  \item{type}{The type of predictions required. There are three options:
    \code{type = "prob"} returns the predicted probabilities,
    \code{type = "choice"} returns the predicted choices, and
    \code{type = "latent"} returns the
    predicted values of the latent variable.
    The default is to return
    all three types of predictions.
  }
  \item{verbose}{logical. If \code{TRUE}, helpful messages along with a
    progress report of drawing from the predictive distributions are
    printed on the screen. The default is \code{FALSE}.
  }
  \item{...}{further arguments passed to or from other methods.}
}

\value{
  \code{predict.mnp} yields a list containing at least one of the
  following elements:
  \item{p}{A matrix of the predicted probabilities where the row
    indicate the observations and the columns indicate the alternatives,
  }
  \item{y}{A matrix of the predicted choices where the rows represent
    the observations and the columns indexes the predictive draws.
  }
  \item{w}{A 3 dimensional arrays of predicted values for the latent
    variable where the first dimension represents the alternatives,
     the second dimension indexes the observations, and the third
    dimension indexes the predictive draws.
  }
}

\seealso{\code{mnp}; MNP home page at
  \url{http://www.princeton.edu/~kimai/research/MNP.html}}

\author{
  Kosuke Imai, Department of Politics, Princeton University
  \email{kimai@Princeton.Edu}
}

\keyword{methods}
