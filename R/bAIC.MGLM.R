# extract AIC or BIC 
# Author: Juhyun Kim 
##============================================================## 

#' @name AIC
#' @title Akaike's Information Criterion (AIC)
#' @description Calculates the Akaike's information criterion (AIC) for a fitted model object. 
#' @param object MGLM object. \code{"MGLMfit"}, \code{"MGLMreg"},
#' \code{"MGLMsparsereg"}, or \code{"MGLMtune"}
#' @return Returns a numeric value with the corresponding AIC.
#' 
#' For the class \code{"MGLMtune"}, the function returns AIC 
#' based on the optimal tuning parameter.
#'
#' @examples 
#' set.seed(124)
#' n <- 200
#' d <- 4
#' alpha <- rep(1, d-1)
#' beta <- rep(1, d-1)
#' m <- 50
#' Y <- rgdirmn(n, m, alpha, beta)
#' gdmFit <- MGLMfit(Y, dist="GDM")
#' AIC(gdmFit)
NULL 

#' @name BIC
#' @title Bayesian information criterion (BIC)   
#' @description Calculates the Bayesian information criterion (BIC) for a fitted model object. 
#' @param object MGLM object. \code{"MGLMfit"}, \code{"MGLMreg"},
#' \code{"MGLMsparsereg"}, or \code{"MGLMtune"}
#' @return Returns a numeric value with the corresponding BIC.
#' 
#' For the class \code{"MGLMtune"}, the function returns BIC 
#' based on the optimal tuning parameter.
#' @examples 
#' set.seed(124)
#' n <- 200
#' d <- 4
#' alpha <- rep(1, d-1)
#' beta <- rep(1, d-1)
#' m <- 50
#' Y <- rgdirmn(n, m, alpha, beta)
#' gdmFit <- MGLMfit(Y, dist="GDM")
#' BIC(gdmFit)
NULL 



aicMGLM <- function(object, tune = FALSE) {
  if (tune) object@select@AIC 
  else object@AIC
}

bicMGLM <- function(object, tune = FALSE) {
  if (tune) object@select@BIC 
  else object@BIC
}

#' @rdname AIC 
#' @exportMethod AIC
setMethod("AIC", "MGLMfit", function(object) aicMGLM(object))
#' @rdname AIC
#' @exportMethod AIC
setMethod("AIC", "MGLMreg", function(object) aicMGLM(object))
#' @rdname AIC
#' @exportMethod AIC
setMethod("AIC", "MGLMsparsereg", function(object) aicMGLM(object))
#' @rdname AIC
#' @exportMethod AIC
setMethod("AIC", "MGLMtune", function(object) aicMGLM(object, TRUE))

#' @rdname BIC
#' @exportMethod BIC
setMethod("BIC", "MGLMfit", function(object) bicMGLM(object))
#' @rdname BIC
#' @exportMethod BIC
setMethod("BIC", "MGLMreg", function(object) bicMGLM(object))
#' @rdname BIC
#' @exportMethod BIC
setMethod("BIC", "MGLMsparsereg", function(object) bicMGLM(object))
#' @rdname BIC
#' @exportMethod BIC
setMethod("BIC", "MGLMtune", function(object) bicMGLM(object, TRUE))