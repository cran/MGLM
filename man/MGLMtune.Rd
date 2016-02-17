\name{MGLMtune}
\alias{MGLMtune}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Choose the tuning parameter value in sparse regression.
}
\description{
Finds the tuning parameter value that yields the smallest BIC.
}
\usage{
MGLMtune(formula, data, dist, penalty, lambdas, ngridpt, warm.start=TRUE, 
  keep.path=FALSE, display=FALSE, init, weight, penidx, ridgedelta, 
  maxiters=150, epsilon=1e-5, regBeta=FALSE, overdisp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{formula}{an object of class \code{formula} (or one that can be coerced
  to that class): a symbolic description of the model to be fitted.
  The response has to be on the left hand side of ~. 
}
  \item{data}{
  an optional data frame, list or environment (or object coercible by 
  \code{as.data.frame} to a data frame) containing the variables in the model.
  If not found in \code{data} when using function \code{MGLMtune}, the variables 
  are taken from \code{environment(formula)}, typically the environment from
  which \code{MGLMtune} is called.
}
  \item{dist}{
	a description of the distribution to fit. See \code{\link{dist}} for the details.
}
  \item{penalty}{
	penalty type for the regularization term. Can be chosen from \code{"sweep"}, 
	\code{"group"}, or \code{"nuclear"}. See \link{MGLMsparsereg} for the description of each penalty type. 
}
  \item{lambdas}{
	an optional vector of the penalty values to tune.  If missing, the vector of 
	penalty values will be set inside the function.  \code{ngridpt} must
	be provided if \code{lambdas} is missing.
}
  \item{ngridpt}{
	an optional numeric variable specifying the number of grid points to tune.  If \code{lambdas} is given, \code{ngridpt}
	will be ignored.  Otherwise, the maximum \eqn{\lambda} is
	determined from the data.  The smallest \eqn{\lambda}
	is set to \eqn{1/n}, where \eqn{n} is the sample size.
}
  \item{warm.start}{
	an optional logical variable to specify whether to give warm start at each
	tuning grid point.  If \code{warm.start=TRUE}, the fitted sparse regression 
	coefficients will be used as the initial value when fitting the sparse
	regression with the next tuning grid. 
}
  \item{keep.path}{
  an optional logical variable controling whether to output the whole solution path. 
  The default is \code{keep.path=FALSE}. If \code{keep.path=TRUE}, the sparse regression result at each grid point will be kept, and saved in the output object 
  \code{select.list}.
  }
  \item{display}{
	an optional logical variable to specify whether to show each tuning 
	step.
  }
  \item{weight}{
  an optional vector of weights assigned to each row of the data. 
	Should be \code{NULL} or a numeric vector. Could be a variable from the \code{data}, 
	or a variable from \code{environment(formula)} with the length equal to
	the number of rows of the data.
	If \code{weight=NULL}, equal weights of ones will be assigned.
  }
  \item{init}{
	an optional matrix of initial value of the parameter estimates.
	Should have the compatible dimension with the data. See \code{\link{dist}} for
	details of dimensions in each distribution.
  }
  \item{penidx}{
	a logical vector indicating the variables to be penalized. The default value is \code{rep(TRUE, p)}, which means all predictors are subject to regularization. If \code{X} contains intercept, use \code{penidx=c(FALSE,rep(TRUE,p-1))}.  
  }

  \item{maxiters}{
	an optional numeric controlling the maximum number of iterations. The default
	value is \code{maxiters=150}.
  }
  \item{ridgedelta}{
  an optional numeric controlling the behavior of the Nesterov's accelerated proximal gradient method. The default value is \eqn{\frac{1}{pd}}.
}
  \item{epsilon}{
  an optional numeric controlling the stopping criterion. The algorithm terminates when the relative change in the objective values of two successive iterates is less then \code{epsilon}. The default value is \code{epsilon=1e-5}.
  }
  \item{regBeta}{
	an optional logical variable used when running negative multinomial regression (\code{dist="NegMN"}). \code{regBeta} controls whether to run regression on the over-dispersion parameter.
	The default is \code{regBeta=FALSE}.
  }
  \item{overdisp}{
	an optional numerical variable used only when fitting sparse negative multinomial 
  model and \code{regBeta=FALSE}.  \code{overdisp} gives the over-dispersion value
	for all the observations.  The default value is estimated using negative-multinomial
	regression.  When \code{dist="MN", "DM", "GDM"} or \code{regBeta=TRUE}, the value of 
	\code{overdisp} is ignored.
  }
}
\value{
\item{select}{
	the final sparse regression result, using the optimal
	tuning parameter.
}
\item{path}{
	a data frame with degrees of freedom and BICs at each lambda.
}
}
\author{
Yiwen Zhang and Hua Zhou
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{MGLMsparsereg}}
}

\examples{
\dontrun{
set.seed(345)
n <- 200
d <- 4
p <- 30
beta <- matrix(0, p, d-1)
beta[c(1,3,5), ] <- 1
beta <- matrix(0, p, d-1)
beta[c(1,3,5), ] <- 1
X <- matrix(rnorm(n*p),n, p)
Alpha <- exp(X\%*\%alpha)
Beta <- exp(X\%*\%beta)
m <- rbinom(n, 80, 0.5)
Y <- rgdirm(m, Alpha, Beta)

select <- MGLMtune(Y~0+X, dist="GDM", penalty="group", ngridpt=30, display=TRUE)
}
}
\keyword{Sparse regression}

