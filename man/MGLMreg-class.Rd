\name{MGLMreg-class}
\Rdversion{1.1}
\docType{class}
\alias{MGLMreg-class}

\title{Class \code{"MGLMreg"}}
\description{
Objects can be created by calls of the form \code{new("MGLMreg", ...)}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MGLMreg", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{call}:}{Object of class \code{"function"}.}
    \item{\code{data}:}{Object of class \code{"list"} ,
    consists of both the predictor matrix and the response matrix. }
    \item{\code{coefficients}:}{Object of class \code{"matrix"},
    the estimated parameters. }
    \item{\code{SE}:}{Object of class \code{"matrix"},
    the standard errors of the parameters. }
    \item{\code{test}:}{Object of class \code{"matrix"},
    the test statistics and p-values.}
    \item{\code{Hessian}:}{Object of class \code{"matrix"},
    the Hessian matrix.}
    \item{\code{logL}:}{Object of class \code{"numeric"},
    the loglikelihood.}
    \item{\code{BIC}:}{Object of class \code{"numeric"},
    Bayesian information criterion.}
	\item{\code{AIC}:}{Object of class \code{"numeric"},
	Akaike information criterion.}
    \item{\code{iter}:}{Object of class \code{"numeric"},
    the number of iterations used.}
    \item{\code{LRT}:}{Object of class \code{"numeric"},
    likelihood ratio test statistics.}
    \item{\code{distribution}:}{Object of class \code{"character"},
    the error distribution fitted.}
	\item{\code{fitted}:}{Object of class \code{"matrix"},
	  the fitted value.}
  }
}
\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "MGLMreg")}: 
		print out regression results from the class \code{MGLMreg}.}
    \item{plot}{\code{signature(x = "MGLMreg")}: 
  	plot regression results from the class \code{MGLMreg}.}
	 }
}
\author{
Yiwen Zhang and Hua Zhou
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
showClass("MGLMreg")
}
\keyword{classes}
