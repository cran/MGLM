\name{print-methods}
\docType{methods}
\alias{print-methods}
\alias{print,ANY-method}
\alias{print,MGLMfit-method}
\alias{print,MGLMreg-method}
\alias{print,MGLMsparsereg-method}
\alias{print,MGLMtune-method}
\title{Print out the defined objects }
\description{
Print out the objects of class \code{"MGLMfit"}, \code{"MGLMreg"},
\code{"MGLMsparsereg"}, \code{"MGLMtune"} in a nice way. 
}
\section{Methods}{
\describe{

\item{\code{signature(x = "ANY")}}{
Print main function output.
}

\item{\code{signature(x = "MGLMfit")}}{
Print main function output for the fitted distribution.
}

\item{\code{signature(x = "MGLMreg")}}{
Print main function output for regression.
}

\item{\code{signature(x = "MGLMsparsereg")}}{
Print main function output for sparse regression.
}

\item{\code{signature(x = "MGLMtune")}}{
Print main function output for sparse regression tuning process.
}

}}
\keyword{methods}

