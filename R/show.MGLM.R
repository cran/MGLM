## Set method to show the output 
## Author: Yiwen Zhang
##============================================================## 

#' @name show
#' @aliases show
#' @title Show an object  
#' @docType methods
#' @description Display the object by printing its class. 
#' @param object an object to be printed. Should be of class \code{"MGLMfit"}, \code{"MGLMreg"},
#' \code{"MGLMsparsereg"} or \code{"MGLMtune"}. 
#' @examples 
#' library("MGLM")
#' data("rnaseq")
#' data <- rnaseq[, 1:6]
#' gdmFit <- MGLMfit(data, dist = "GDM")
#' show(gdmFit)
#' @importFrom methods show
#' @exportMethod show 
NULL

#' @rdname show 
setMethod("show", signature = "MGLMfit", function(object) {
    digits = 6
    print(cbind(estimate = round(object@estimate, digits), SE = round(object@SE, digits)))
    cat("\n")
    cat("Distribution: ", object@distribution, "\n", sep = "")
    cat("Log-likelihood: ", object@logL, "\n", sep = "")
    cat("BIC: ", object@BIC, "\n", sep = "")
    cat("AIC: ", object@AIC, "\n", sep = "")
    cat("LRT test p value: ", ifelse(object@LRTpvalue > 1e-04, sprintf("%.3f", object@LRTpvalue), 
        "<0.0001"), "\n", sep = "")
    cat("Iterations: ", object@iter, "\n", sep = "")
})


#' @rdname show 
setMethod("show", signature = "MGLMreg", function(object) {
    digits = 6
    cat("Call: ")
    print(object@call)
    cat("\n")
    cat("Coefficients:\n")
    if (is.numeric(object@coefficients)) {
      print(round(object@coefficients, digits))
    } else if (is.list(object@coefficients)) {
      print(lapply(object@coefficients, round, digits))
    }
    cat("\n")
    
    
    cat("Hypothesis test: \n")
    if (all(is.numeric(object@test[, 2]) & object@test[, 2] >= 1e-6)) {
        test <- cbind(object@test[, 1], round(object@test[, 2], digits))
    }
    else {
      test <- cbind(object@test[, 1], object@test[, 2])
    }
    colnames(test) <- colnames(object@test)
    print(test)
    cat("\n")
    
    cat("Distribution: ", object@distribution, "\n", sep = "")
    cat("Log-likelihood: ", object@logL, "\n", sep = "")
    cat("BIC: ", object@BIC, "\n", sep = "")
    cat("AIC: ", object@AIC, "\n", sep = "")
    cat("Iterations: ", object@iter, "\n", sep = "")
})

#' @rdname show 
setMethod("show", signature = "MGLMsparsereg", function(object) {
    cat("Call: ")
    print(object@call)
    cat("\n")
    cat("Distribution: ", object@distribution, "\n", sep = "")
    cat("Log-likelihood: ", object@logL, "\n", sep = "")
    cat("BIC: ", object@BIC, "\n", sep = "")
    cat("AIC: ", object@AIC, "\n", sep = "")
    cat("Degrees of freedom: ", object@Dof, "\n", sep = "")
    cat("Lambda: ", object@lambda, "\n", sep = "")
    if (!is.null(object@maxlambda)) {
        cat("Max lambda: ", object@maxlambda, "\n", sep = "")
    }
    cat("Iterations: ", object@iter, "\n", sep = "")
})

#' @rdname show 
setMethod("show", signature = "MGLMtune", function(object) {
    cat("Call: ")
    print(object@select@call)
    cat("\n")
    cat("Distribution: ", object@select@distribution, "\n", sep = "")
    cat("Log-likelihood: ", object@select@logL, "\n", sep = "")
    cat("BIC: ", object@select@BIC, "\n", sep = "")
    cat("AIC: ", object@select@AIC, "\n", sep = "")
    cat("Degrees of freedom: ", object@select@Dof, "\n", sep = "")
    cat("Lambda: ", object@select@lambda, "\n", sep = "")
    cat("Number of grid points: ", nrow(object@path), "\n", sep = "")
}) 
