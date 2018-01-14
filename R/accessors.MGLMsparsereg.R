# extract maxlambda from MGLMsparsereg class
##============================================================## 

#' @name maxlambda 
#' @title Extract maximum lambda 
#' @description \code{maxlambda} extracts the maximum tuning parameter that ensures 
#' the estimated regression coefficients are not all zero for the object of class \code{MGLMsparsereg}. 
#' @param object an object of class \code{MGLMsparsereg} from which
#' maximum lambda value can be extracted.
#' @return Returns a maximum lambda value of \code{object}.
#' @examples
#' library("MGLM")
#' dist <- "DM"
#' n <- 100
#' p <- 10
#' d <- 5
#' set.seed(118)
#' m <- rbinom(n, 200, 0.8)
#' X <- matrix(rnorm(n * p), n, p)
#' alpha <- matrix(0, p, d)
#' alpha[c(1, 3, 5), ] <- 1
#' Alpha <- exp(X %*% alpha)
#' Y <- rdirmn(size = m, alpha = Alpha)
#' pen <- "group"
#' ngridpt <- 30
#' spmodelfit <- MGLMsparsereg(formula = Y ~ 0 + X, dist = dist, 
#'                             lambda = Inf, penalty = pen)
#' maxlambda <- maxlambda(spmodelfit)
NULL 


maxlambdaMGLM <- function(object) {
  object@maxlambda
}


#' @rdname maxlambda
#' @exportMethod maxlambda
setMethod("maxlambda", "MGLMsparsereg", function(object) maxlambdaMGLM(object))





# extract degrees of freedom from MGLMsparsereg class
##============================================================## 

#' @name dof 
#' @title Extract degrees of freedom  
#' @description \code{dof} extracts the degrees of freedom of the estimated parameter 
#' from the object of class \code{MGLMsparsereg}. 
#' @param object an object of class \code{MGLMsparsereg} 
#' @return Returns degrees of freedom of \code{object}.
#' @examples
#' library("MGLM")
#' dist <- "DM"
#' n <- 100
#' p <- 10
#' d <- 5
#' set.seed(118)
#' m <- rbinom(n, 200, 0.8)
#' X <- matrix(rnorm(n * p), n, p)
#' alpha <- matrix(0, p, d)
#' alpha[c(1, 3, 5), ] <- 1
#' Alpha <- exp(X %*% alpha)
#' Y <- rdirmn(size = m, alpha = Alpha)
#' pen <- "group"
#' ngridpt <- 30
#' spmodelfit <- MGLMsparsereg(formula = Y ~ 0 + X, dist = dist, 
#'                             lambda = Inf, penalty = pen)
#' df <- dof(spmodelfit)
NULL 


dofMGLM <- function(object) {
  object@Dof
}


#' @rdname dof
#' @exportMethod dof
setMethod("dof", "MGLMsparsereg", function(object) dofMGLM(object))