# logLik method to extract log-likelihood 
# Author: Juhyun Kim 
##============================================================## 

#' @name logLik
#' @title Extract log-likelihood 
#' @description \code{logLik} extracts log-likelihood for classes \code{"MGLMfit"}, 
#' \code{"MGLMreg"}, \code{"MGLMsparsereg"}. 
#' @param object an object from which a log-likelihood value can be extracted.
#' @importFrom stats4 logLik
#' @return Returns a log-likelihood value of \code{object}.
#' @examples
#' library("MGLM")
#' data("rnaseq")
#' data <- rnaseq[, 1:6]
#' dmFit <- MGLMfit(data, dist = "DM")
#' logLik(dmFit)
NULL 


logLikMGLM <- function(object) {
	object@logL
}


#' @rdname logLik
#' @exportMethod logLik
setMethod("logLik", "MGLMfit", logLikMGLM)

#' @rdname logLik
#' @exportMethod logLik
setMethod("logLik", "MGLMreg", logLikMGLM)

#' @rdname logLik
#' @exportMethod logLik
setMethod("logLik", "MGLMsparsereg", logLikMGLM)