\name{MGLMreg}
\alias{MGLMreg}
\alias{DMD.MN.reg}
\alias{DMD.DM.reg}
\alias{DMD.GDM.reg}
\alias{DMD.NegMN.reg}
\alias{DMD.NegMN.Alpha.reg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit multivariate response GLM regression
}
\description{
\code{MGLMreg} fits the multivariate response generalized linear models,
specified by a symbolic description of the linear predictor and
a description of the error distribution.
}
\usage{
MGLMreg(formula, data, dist,init, weight,  epsilon=1e-8, 
maxiters=150, display=FALSE, LRT=FALSE, parallel=FALSE, cores, regBeta=FALSE)

DMD.MN.reg(Y, init, X, weight, epsilon, maxiters, display, parallel, cores, cl, sys)
DMD.DM.reg(Y, init, X, weight,  epsilon, maxiters, display, parallel, cores, cl, sys)
DMD.GDM.reg(Y,  init, X, weight, epsilon, maxiters, display, parallel, cores, cl, sys)
DMD.NegMN.reg(Y, init, X, weight, epsilon, maxiters, display, parallel, cores, cl, sys)
DMD.NegMN.Alpha.reg(
	Y, init, X, weight, epsilon, maxiters, display, parallel, cores, cl, sys)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{an object of class \code{formula} (or one that can be coerced
  to that class): a symbolic description of the model to be fitted.
  The response has to be on the left hand side of ~.
}
  \item{data}{
  an optional data frame, list or environment (or object coercible by 
  as.data.frame to a data frame) containing the variables in the model.
  If not found in \code{data} when using function \code{MGLMreg}, the variables 
  are taken from \code{environment(formula)}, typically the environment from
  which \code{MGLMreg} is called.
}
  \item{Y, X}{
  for \code{DMD.MN.reg}, \code{DMD.DM.reg}, \code{DMD.GDM.reg}, and
  \code{DMD.NegMN.reg}, \code{X} is a design matrix of dimension 
  \code{n*(p+1)} and \code{Y} is the response matrix of dimension 
  \code{n*d}.
}
  \item{dist}{
	a description of the error distribution to fit. See \code{\link{dist}} for details. 
}
  \item{weight}{
	an optional vector of weights assigned to each row of the data. 
	Should be \code{NULL} or a numeric vector. Could be a variable from \code{data}, 
	or a variable from \code{environment(formula)} with the length equal to
	the number of rows of the data.
	If \code{weight=NULL}, equal weights of ones will be assigned.
}
  \item{init}{
	an optional matrix of initial value of the parameter estimates.
	Should have the compatible dimension with \code{data}. See \code{\link{dist}} for
	details of the dimensions in each distribution.
}
  \item{epsilon}{
	an optional numeric controlling the stopping criterion. The algorithm terminates when the relative change in the loglikelihoods of two successive iterates is less than \code{epsilon}.  The default value is \code{epsilon=1e-8}.
}
  \item{maxiters}{
	an optional numeric controlling the maximum number of iterations. The default
	value is \code{maxiters=150}.
}
  \item{display}{
	an optional logical variable controlling the display of iterations. The default
	value is \code{display=FALSE}.
}
\item{LRT}{
	an optional logical variable controlling whether to perform likelihood ratio test on
	each predictor. The default value is \code{LRT=FALSE}, in which case only the Wald test is performed.
}
\item{parallel}{
   an optional logical variable controlling whether to perform parallel computing. 
   On a multi-core Windows machine, a cluster is created based on socket; 
	on a multi-core Linux/Mac machine, a cluster is created based on forking. The default
  value is \code{parallel=FALSE}.
}
\item{cores}{
	an optional value specifying the number of cores to use.  
	Default value is half of the logical cores.
}
\item{cl}{
	a cluster object, created by the package \pkg{parallel} or by package \pkg{snow}. 
	If \code{parallel=TRUE}, use the registered default cluster; if \code{parallel=FALSE}, any given value to \code{cl} will be ignored.
}
\item{sys}{
	the operating system.  Will be used when choosing parallel type. 
}	
\item{regBeta}{
	an optional logical variable.  When \code{dist="NegMN"}, the user can decide 
	whether to run regression on the over dispersion parameter \eqn{\beta}.  The 
	default is \code{regBeta=FALSE}. 
}
}

\details{
The formula should be in the form responses ~ covariates where the 
responses are the multivariate count matrix or a few columns from
a data frame which is specified by \code{data}.  
The covariates are either matrices or from the data
frame.  The covariates can be numeric or character or factor. 

See \code{\link{dist}} for details about distributions. 
}

\value{
Returns an object of class \code{"MGLMreg"}. 
An object of class \code{"MGLMreg"} is a list 
containing at least the following components: 

\item{coefficients}{the estimated regression coefficients.}
\item{SE}{the standard errors of the estimates.}
\item{Hessian}{the Hessian at the estimated parameter values.}
\item{gradient}{the gradient at the estimated parameter values.}
\item{wald.value}{the Wald statistics.}
\item{wald.p}{the p values of Wald test.}
\item{test}{test statistic and the corresponding p-value. if 
\code{LRT=FALSE}, only returns test results
from Wald test; if \code{LRT=TRUE}, returns the test results from
both Wald test and likelihood ratio test.}
\item{logL}{the final loglikelihood.}
\item{BIC}{Bayesian information criterion. }
\item{AIC}{Akaike information criterion.}
\item{fitted}{the fitted values from the regression model}
\item{iter}{the number of iterations used.}
\item{call}{the matched call.}
\item{distribution}{the distribution fitted.}
\item{data}{the data used to fit the model.}
}
\author{
Yiwen Zhang and Hua Zhou
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link{MGLMfit}} for distribution fitting.
}
\examples{
##----------------------------------------##
## Generate data
n <- 2000
p <- 5
d <- 4
m <- rep(20, n)
set.seed(1234)
X <- 0.1* matrix(rnorm(n*p),n, p)
alpha <- matrix(1, p, d-1)
beta <- matrix(1, p, d-1)
Alpha <- exp(X\%*\%alpha)
Beta <- exp(X\%*\%beta)
gdm.Y <- rgdirm(m, Alpha, Beta)

##----------------------------------------##
## Regression
\donttest{
gdm.reg <- MGLMreg(gdm.Y~X, dist="GDM", LRT=FALSE)
plot(gdm.reg)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Models}
\keyword{Regression}
