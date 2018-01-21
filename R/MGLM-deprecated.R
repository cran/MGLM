#' Deprecated function(s) in the MGLM package
#' 
#' These functions are provided for compatibility with older version of
#' the yourPackageName package.  They may eventually be completely
#' removed.
#' @name MGLM-deprecated
#' @param ... parameters to be passed to the modern version of the function
#' @aliases ddirm dgdirm dneg rdirm rgdirm 
#' @section Details:
#' \tabular{rl}{
#'   \code{ddirm} \tab now a synonym for \code{\link{ddirmn}}\cr
#'   \code{dgdirm} \tab now a synonym for \code{\link{dgdirmn}}\cr
#'   \code{dneg} \tab now a synonym for \code{\link{dnegmn}}\cr
#'   \code{rdirm} \tab now a synonym for \code{\link{rdirmn}}\cr
#'   \code{rgdirm} \tab now a synonym for \code{\link{rgdirmn}}\cr
#' }
#' 
#' The function \code{dneg} has been deprecated. Use \code{dnegmn} instead.
#' 
#' Note the change in argument order: 
#' \code{dneg(Y, prob, beta)} and \code{dnegmn(Y, alpha, beta)} from MGLM_0.0.8 have been deprecated;
#' use \code{dnegmn(Y, beta, prob = alpha/(rowSums(alpha)+1), alpha=NULL)} instead.
NULL 

#' @rdname MGLM-deprecated
#' @export 
ddirm <- function(...) { 
  .Deprecated("ddirmn")
  ddirmn(...)
}
#' @rdname MGLM-deprecated
#' @export 
rdirm <- function(...) {
  .Deprecated("rdirmn")
  rdirmn(...)
}
#' @rdname MGLM-deprecated
#' @export 
dgdirm <- function(...) { 
  .Deprecated("dgdirmn")
  dgdirmn(...)
}
#' @rdname MGLM-deprecated
#' @export 
rgdirm <- function(...) {
  .Deprecated("rgdirmn")
  rgdirmn(...)
}
#' @rdname MGLM-deprecated
#' @param Y,alpha,beta for functions \code{dnegmn}, note the change in argument order. See Details.  
#' @export  
dneg <- function(Y, alpha, beta) { 
  .Deprecated("dnegmn")
  warning(" note the deprecated argument order;\n dneg(Y, prob, beta) and dnegmn(Y, alpha, beta) from MGLM_0.0.8 have been deprecated;\n use dnegmn(Y, beta, prob = alpha/(rowSums(alpha)+1), alpha=NULL) instead", 
          call. = FALSE)
  dnegmn(Y, beta, alpha=alpha)
}