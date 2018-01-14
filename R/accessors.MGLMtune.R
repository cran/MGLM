# extract path from MGLMtune class
# Author: Juhyun Kim 
##============================================================## 

#' @name path 
#' @title Extract path 
#' @description \code{path} extracts from object of the class \code{MGLMtune} the path of 
#' BIC, AIC, log-likelihood and degrees of freedom given each tuning parameter.
#' @param object an object of class \code{MGLMtune} from which
#' path can be extracted.
#' @return Returns a path of \code{object}.
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
#' select <- MGLMtune(Y ~ 0 + X, dist = "DM", penalty = "nuclear", 
#' ngridpt = 10, display = FALSE)
#' select_path <- path(select)
NULL 


pathMGLM <- function(object) {
  object@path
}


#' @rdname path
#' @exportMethod path
setMethod("path", "MGLMtune", function(object) pathMGLM(object))