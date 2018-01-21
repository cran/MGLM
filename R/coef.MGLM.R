# Coef method to extract coefficients 
# Author: Juhyun Kim 
##============================================================## 

#' Method coef.
#' @name coef
#' @title Extract Model Coefficients 
#' @description \code{coef} extracts estimated model coefficients of class. \code{coefficients} is an \emph{alias} for it. 
#' @param object an object for which the extraction of model coefficients is meaningful. 
#' One of the following classes \code{"MGLMfit"}, \code{"MGLMreg"},
#' \code{"MGLMsparsereg"}, \code{"MGLMtune"}
#' @importFrom stats4 coef
#' @return Coefficients extracted from the model object \code{object}.
#' 
#' For the class \code{"MGLMtune"}, the function returns model coefficients 
#' based on the optimal tuning parameter.
#' @exportMethod coef
#' @examples
#' library("MGLM")
#' data("rnaseq")
#' data <- rnaseq[, 1:6]
#' mnreg <- MGLMreg(formula = cbind(X1, X2, X3, X4, X5, X6) ~ log(totalReads) + 
#' treatment + age + gender, data = rnaseq, dist = "MN")
#' coef(mnreg)
NULL 


coefMGLM <- function(object, tune = FALSE) {
  if (tune) object@select@coefficients 
  else object@coefficients 
}

 
setMethod("coefficients", "MGLMfit", function(object) coefMGLM(object))
#' @rdname coef
setMethod("coef", "MGLMfit", function(object) coefMGLM(object))
setMethod("coefficients", "MGLMreg", function(object) coefMGLM(object))
#' @rdname coef
setMethod("coef", "MGLMreg", function(object) coefMGLM(object))
setMethod("coefficients", "MGLMsparsereg", function(object) coefMGLM(object))
#' @rdname coef
setMethod("coef", "MGLMsparsereg", function(object) coefMGLM(object))
setMethod("coefficients", "MGLMtune", function(object) coefMGLM(object, TRUE))
#' @rdname coef
setMethod("coef", "MGLMtune", function(object) coefMGLM(object, TRUE))