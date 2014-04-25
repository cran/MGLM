\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot,ANY-method}
\alias{plot,MGLMreg-method}
\alias{plot,MGLMsparsereg-method}
\alias{plot,MGLMtune-method}
\title{plot the defined objects }
\description{
Plot the objects of class \code{"MGLMreg"},
\code{"MGLMsparsereg"}, \code{"MGLMtune"} in a nice way. 
}
\section{Methods}{
\describe{

\item{\code{signature(x = "ANY")}}{
Plot main function output.
}

\item{\code{signature(x = "MGLMreg")}}{
Plot main function output for regression.
}

\item{\code{signature(x = "MGLMsparsereg")}}{
Plot main function output for sparse regression.
}

\item{\code{signature(x = "MGLMtune")}}{
Plot main function output for sparse regression tuning process.
}

}}
\keyword{methods}
