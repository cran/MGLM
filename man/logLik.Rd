% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logLik.MGLM.R
\name{logLik}
\alias{logLik}
\alias{logLik,MGLMfit-method}
\alias{logLik,MGLMreg-method}
\alias{logLik,MGLMsparsereg-method}
\title{Extract log-likelihood}
\usage{
\S4method{logLik}{MGLMfit}(object)

\S4method{logLik}{MGLMreg}(object)

\S4method{logLik}{MGLMsparsereg}(object)
}
\arguments{
\item{object}{an object from which a log-likelihood value can be extracted.}
}
\value{
Returns a log-likelihood value of \code{object}.
}
\description{
\code{logLik} extracts log-likelihood for classes \code{"MGLMfit"}, 
\code{"MGLMreg"}, \code{"MGLMsparsereg"}.
}
\examples{
library("MGLM")
data("rnaseq")
data <- rnaseq[, 1:6]
dmFit <- MGLMfit(data, dist = "DM")
logLik(dmFit)
}
