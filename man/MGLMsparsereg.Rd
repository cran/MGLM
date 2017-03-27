\name{MGLMsparsereg}
\alias{MGLMsparsereg}
\alias{MGLMsparsereg.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit multivariate GLM sparse regression
}
\description{
Fit sparse regression in multivariate generalized linear models.
}
\usage{
MGLMsparsereg(formula, data, dist, lambda, penalty, weight, 
  init, penidx, maxiters = 150, ridgedelta, epsilon = 1e-5, regBeta=FALSE, overdisp)
MGLMsparsereg.fit(Y, X, dist, lambda, penalty, weight, 
  init, penidx, maxiters = 150, ridgedelta, epsilon = 1e-5, regBeta=FALSE, overdisp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{formula}{An object of class \code{formula} (or one that can be coerced
  to that class): a symbolic description of the model to be fitted.
  The response has to be on the left hand side of ~. 
}
  \item{data}{
  An optional data frame, list or environment (or object coercible by 
  \code{as.data.frame} to a data frame) containing the variables in the model.
  If not found in \code{data} when using function \code{MGLMsparsereg}, the variables 
  are taken from \code{environment(formula)}, typically the environment from
  which \code{MGLMsparsereg} is called.
}
  \item{Y}{
  A matrix containing the multivariate categorical response data. 
  Rows of the matrix represent observations, while columns are the different
  categories.  Rows and columns of all zeros are automatically removed when
  \code{dist="MN"}, \code{"DM"}, or \code{"GDM"}.
  }
  \item{X}{
   Design matrix (including intercept).  
   Number of rows of the matrix should match that of \code{Y}. 
  }
  \item{dist}{
	A description of the error distribution to fit. See \code{\link{dist}} for details.
  }
  \item{weight}{
    An optional vector of weights assigned to each row of the data. 
	Should be \code{NULL} or a numeric vector. Could be a variable from \code{data}, 
	or a variable from \code{environment(formula)} with the length equal to
	the number of rows of the data.
	If \code{weight=NULL}, equal weights of ones will be assigned.
  }
  \item{init}{
	An optional matrix of initial value of the parameter estimates.
	Should have the compatible dimension with the data. See \code{\link{dist}} for
	details of the dimensions in each distribution.
}
  \item{lambda}{
	Penalty parameter.  
  }
  \item{penalty}{
	Penalty type for the regularization term. Can be chosen from \code{"sweep"}, 
	\code{"group"}, or \code{"nuclear"}. See Details for the description of each penalty type. 
  }
  \item{penidx}{
	A logical vector indicating the variables to be penalized. The default value is \code{rep(TRUE, p)}, which means all predictors are subject to regularization. If \code{X} contains intercept, use \code{penidx=c(FALSE,rep(TRUE,p-1))}.
  }
  \item{maxiters}{
	An optional numeric controlling the maximum number of iterations. The default
	value is maxiters=150.
  }
  \item{ridgedelta}{
	An optional numeric controlling the behavior of the Nesterov's accelerated proximal gradient method. The default value is \eqn{\frac{1}{pd}}.
  }
  \item{epsilon}{
  An optional numeric controlling the stopping criterion. The algorithm terminates when the relative change in the objective values of two successive iterates is less then \code{epsilon}.
  	The default value is \code{epsilon=1e-5}.
  }
  \item{regBeta}{
	An optional logical variable used when running negative multinomial regression (\code{dist="NegMN"}).  
	\code{regBeta} controls whether to run regression on the over-dispersion parameter.
	The default is \code{regBeta=FALSE}.
  }
  \item{overdisp}{
	An optional numerical variable used only when fitting sparse negative multinomial 
	model \code{dist="NegMN"} and \code{regBeta=FALSE}.  \code{overdisp} gives the over dispersion value
	for all the observations.  The default value is estimated using negative-multinomial regression.  When \code{dist="MN", "DM", "GDM"} or \code{regBeta=TRUE}, the value of \code{overdisp} is ignored.
  }
}
\details{
In general, we consider regularization problem
\deqn{
\min_B h(B) = -l(B)+ J(B),
}
where \eqn{l(B)} is the loglikelihood function and \eqn{J(B)} is the 
regularization function.  

Sparsity in the individual elements of the parameter matrix \eqn{B} is achieved 
by the lasso penalty (\code{dist="sweep"})
\deqn{
J(B) = \lambda \sum_{k\in penidx} \sum_{j=1}^d \|B_{kj}\|
}

Sparsity in the rows of the regression parameter matrix \eqn{B} is achieved
by the group penalty (\code{dist="group"})
\deqn{
J(B) = \lambda \sum_{k \in penidx} \|B_{k \cdot}\|_2,
}
where \eqn{\|v\|_2} is the \eqn{l_2} norm of a vector \eqn{v}. In other words, 
\eqn{\|B_{k\cdot}\|_2} is the \eqn{l_2} norm of the \eqn{k}-th row of the 
parameter matrix \eqn{B}.

Sparsity in the rank of the parameter matrix \eqn{B} is achieved by the nuclear norm penalty (\code{dist="nuclear"})
\deqn{
J(B) = \lambda \|B\|_*= \lambda \sum_{i=1}^{min(p, d)} \sigma_i(B),
}
where \eqn{\sigma_i(B)} are the singular values of the parameter matrix \eqn{B}. 
The nuclear norm \eqn{\|B\|_*} is a convex relaxation of \eqn{rank(B)=\|\sigma(B)\|_0}.

See \code{\link{dist}} for details about distributions. 
}

\value{
Returns an object of class \code{"MGLMsparsereg"}. 
An object of class \code{"MGLMsparsereg"}
is a list containing at least the following components:  
\item{coefficients}{the estimated matrix of regression coefficients.}
\item{logL}{The final loglikelihood value.}
\item{AIC}{Akaike information criterion.}
\item{BIC}{Bayesian information criterion.}
\item{Dof}{Degrees of freedom of the estimated parameter.}
\item{iter}{Number of iterations used. }
\item{maxlambda}{The maxmum tuning parameter such that 
	the estimated coefficients are not all zero.  This value is returned only
	when the tuning parameter \code{lambda} given to the function is large enough 
	such that all the parameter estimates are zero; otherwise, \code{maxlambda}
	is not computed.}
\item{call}{A matched call.}
\item{data}{The data used to fit the model: a list of the predictor matrix
and the response matrix.}
\item{penalty}{The penalty chosen when running the penalized regression.}
}
\author{
Yiwen Zhang and Hua Zhou
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
\dontrun{
##----------------------------------------##
## Generate Dirichlet Multinomial data
dist <- "DM"
n <- 100
p <- 15
d <- 5
m <- runif(n, min=0, max=25) + 25
set.seed(134)
X <- matrix(rnorm(n*p),n, p)
alpha <- matrix(0, p, d)
alpha[c(1,3, 5), ] <- 1
Alpha <- exp(X\%*\%alpha)
Y <- rdirm(size=m, alpha=Alpha)
##----------------------------------------##

## Tuning
ngridpt <- 10
p <- ncol(X)
d <- ncol(Y)

pen <- 'nuclear'
	fit <- MGLMsparsereg(formula=Y~0+X, dist=dist, lambda=Inf, penalty=pen)
	maxlambda <- fit$maxlambda
	lambdas <- exp(seq(from=log(maxlambda), to=log(maxlambda/100), 
					length.out=ngridpt))
	BICs <- rep(0, ngridpt)
	AICs <- rep(0, ngridpt)
	LogLs <- rep(0, ngridpt)
	Dofs <- rep(0, ngridpt)
	
	
	ptm <- proc.time()
	for(j in 1:ngridpt){
		print(j)
		if(j ==1 ) B0 <- matrix(0, p, ncol(fit$coefficients)) else B0 <- B_hat
		select.fit <- MGLMsparsereg(formula=Y~0+X, dist=dist,lambda=lambdas[j], 
				penalty=pen, init=B0)
		B_hat <- select.fit$coefficients
		BICs[j] <- select.fit$BIC
		LogLs[j] <- select.fit$logL
		AICs[j] <- select.fit$AIC
		Dofs[j] <- select.fit$Dof
	}
	proc.time() - ptm
	
	chosen.lambda <- lambdas[which.min(BICs)]
	select <- MGLMsparsereg(formula=Y~0+X, dist=dist, lambda=chosen.lambda, 
			penalty=pen)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Sparse regression}
