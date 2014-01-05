\name{MGLMfit-class}
\Rdversion{1.1}
\docType{class}
\alias{MGLMfit-class}

\title{Class \code{"MGLMfit"}}

\description{
A list containing the model fitting results from the \code{MGLMfit}.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MGLMfit", ...)}.
}

\section{Slots}{
  \describe{
	\item{\code{estimate}:}{Object of class \code{"vector"},
	containing the parameter estimates.}
	\item{\code{SE}:}{Object of class \code{"vector"},
	containing the standard errors of the estimates.}
	\item{\code{vcov}:}{Object of class \code{"matrix"},
	the variance covariance matrix of the parameter estimates.}
    \item{\code{logL}:}{Object of class \code{"numeric"}, 
    the fitted log likelihood. }
    \item{\code{BIC}:}{Object of class \code{"numeric"}, 
    Bayesian information criterion.}
	\item{\code{AIC}:}{Object of class \code{"numeric"},
	Akaike information criterion.}
	\item{\code{LRTpvalue}:}{Object of class \code{"numeric"},
    likelihood ratio test p value.}
	\item{\code{gradient}:}{Object of class \code{"numeric"},
	containing the gradient.}
    \item{\code{iter}:}{Object of class \code{"numeric"}, 
    number of iteration used.}
    \item{\code{distribution}:}{Object of class \code{"character"},
    the distribution fitted.}
	\item{\code{fitted}:}{Object of class \code{"vector"},
	the fitted mean of each category}
  }
}
\section{Methods}{
\describe{
    \item{print}{\code{signature(x = "MGLMfit")}: 
		print out the results from the class MGLMfit.}
}
}
\author{
Yiwen Zhang and Hua Zhou
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
showClass("MGLMfit")
}
\keyword{classes}
