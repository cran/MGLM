\name{MGLMtune-class}
\Rdversion{1.1}
\docType{class}
\alias{MGLMtune-class}
\alias{print, MGLMtune-method}

\title{Class \code{"MGLMtune"}}

\description{
A list containing the results from the \code{MGLMsparsereg}.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MGLMsparsereg", ...)}.
}
\section{Slots}{
  \describe{  
    \item{\code{call}:}{Object of class \code{"function"}. }
    \item{\code{select}:}{Object of class \code{"MGLMsparsereg"},
    regularized regression results given by the optimal tuning parameter.}
    \item{\code{path}:}{Object of class \code{"data.frame"},
    the BIC, AIC, log-likelihood and degrees of freedom given each tuning parameter.}
    \item{\code{select.list}:}{Object of class \code{"list"},
    the regularized regression results at each tuning grid point. }
  }
}
\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "MGLMtune")}: 
	Print out sparse regression solution path from class \code{"MGLMtune"}.}
	}
}
\author{
Yiwen Zhang and Hua Zhou
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
showClass("MGLMtune")
}
\keyword{classes}
