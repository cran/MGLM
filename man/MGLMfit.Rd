\name{MGLMfit}
\alias{MGLMfit}
\alias{DMD.DM.fit}
\alias{DMD.GDM.fit}
\alias{DMD.NegMN.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit multivariate discrete distributions
}
\description{
Fit the specified multivariate discrete distribution.
}
\usage{
MGLMfit(data, dist, init, weight, epsilon=1e-8, maxiters=150, display=FALSE)

DMD.DM.fit(data, init, weight, epsilon=1e-8, maxiters=150, display=FALSE)
DMD.GDM.fit(data, init, weight, epsilon=1e-8, maxiters=150, display=FALSE)
DMD.NegMN.fit(data, init, weight, epsilon=1e-8, maxiters=150, display=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame or matrix containing the count data. 
	Rows of the matrix represent observations and columns are the categories. 
	Rows and columns of all zeros are automatically removed.
}
  \item{dist}{
A description of the distribution to fit. Choose from \code{"MN"}, \code{"DM"}, \code{"GDM"}, \code{"NegMN"}.  See \code{\link{dist}} for details.
}
  \item{weight}{
An optional vector of weights assigned to each row of the data. 
	Should be Null or a numeric vector with the length equal to
	the number of rows of \code{data}.
	If \code{weight=NULL}, equal weights of all ones will be assigned.
}
  \item{init}{
An optional vector of initial value of the parameter estimates.
	Should have the same dimension as the estimated parameters. See \code{\link{dist}} for
	details.
}
  \item{epsilon}{
An optional numeric controlling the stopping criterion. The algorithm terminates when the relative change in the log-likelihoods of two successive iterates is less than \code{epsilon}. The default value is \code{epsilon=1e-8}.
}
  \item{maxiters}{
An optional number controlling the maximum number of iterations. The default
value is \code{maxiters=150}.
}
  \item{display}{
An optional logical variable controlling the display of iterations. The default
value is FALSE.
}
}
\details{
See \code{\link{dist}} for details about model parameterization.
}
\value{
Returns an object of class \code{"MGLMfit"}.  An object of class \code{"MGLMfit"}
is a list containing at least the following components: 
\item{estimate}{The vector of the distribution prameter estimates.}
\item{SE}{The vector of standard errors of the estimates.}
\item{vcov}{The variance-covariance matrix of the estimates.}
\item{logL}{The loglikelihood value.}
\item{iter}{The number of iterations used.}
\item{BIC}{Bayesian information criterion.}
\item{AIC}{Akaike information criterion.}
\item{distribution}{The distribution fitted.}
\item{LRT}{When \code{dist="DM"} or \code{"GDM"}, it is the likelihood ratio test statistic for
comparing the current model to the multinomial model.  No
LRT provided when \code{dist="NegMN"}.}
\item{LRTpvalue}{The likelihood ratio test P value.}
\item{DoF}{The degrees of freedom of the model.}
}
\author{
Yiwen Zhang and Hua Zhou
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(rnaseq)
Y <- as.matrix(rnaseq[, 1:6])
fit <- MGLMfit(data=Y, dist="GDM") 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Models}
\keyword{Distribution fitting}