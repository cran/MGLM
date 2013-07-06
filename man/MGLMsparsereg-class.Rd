\name{MGLMsparsereg-class}
\Rdversion{1.1}
\docType{class}
\alias{MGLMsparsereg-class}

\title{Class \code{"MGLMsparsereg"}}

\description{
A list containing the results from the \code{MGLMsparsereg}.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MGLMsparsereg", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{call}:}{Object of class \code{"function"}. }
    \item{\code{data}:}{Object of class \code{"list"},
    consists of both the predictor matrix and the response matrix}
    \item{\code{coefficients}:}{Object of class \code{"matrix"},
    the estimated regression coefficients.}
    \item{\code{logL}:}{Object of class \code{"numeric"},
    the loglikelihood. }
    \item{\code{BIC}:}{Object of class \code{"numeric"},
    Bayesian information criterion. }
    \item{\code{AIC}:}{Object of class \code{"numeric"}, 
    Akaike information criterion.}
	\item{\code{Beta}:}{Object of class \code{"numeric"},
	the over dispersion parameter of the negative multinomial regression.}
    \item{\code{Dof}:}{Object of class \code{"numeric"},
    the degrees of freedom.}
    \item{\code{iter}:}{Object of class \code{"numeric"},
    number of iterations used. }
    \item{\code{maxlambda}:}{Object of class \code{"numeric"},
    the maximum tuning parameter that ensures the estimated regression coefficients are not all zero.}
    \item{\code{lambda}:}{Object of class \code{"numeric"},
    the tuning parameter used. }
    \item{\code{distribution}:}{Object of class \code{"character"},
    the distribution fitted. }
    \item{\code{penalty}:}{Object of class \code{"character"},
    the chosen penalty when running penalized regression. }
  }
}
\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "MGLMsparsereg")}: 
	print out sparse regression results from class \code{"MGLMsparsereg"}.}
	 }
}
\author{
Yiwen Zhang and Hua Zhou
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
showClass("MGLMsparsereg")
}
\keyword{classes}
