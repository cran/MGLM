# Fit multivariate response distribution 
# Author: Hua Zhou and Yiwen Zhang
#============================================================##


##============================================================## 
## Distribution fitting function 
##============================================================##
#' @name MGLMfit
#' @aliases MGLMfit
#' @title Fit multivariate discrete distributions
#'
#' @description Fit the specified multivariate discrete distribution.
#'
#' @param data a data frame or matrix containing the count data. 
#' Rows of the matrix represent observations and columns are the categories. 
#' Rows and columns of all zeros are automatically removed.
#' @param dist a description of the distribution to fit. Choose from \code{"MN"}, \code{"DM"}, \code{"GDM"}, \code{"NegMN"}. See \code{\link{dist}} for details.
#' @param weight an optional vector of weights assigned to each row of the data. Should be Null or a numeric vector with the length equal to the number of rows of \code{data}. 
#' If \code{weight=NULL}, equal weights of all ones will be assigned.
#' @param init an optional vector of initial value of the parameter estimates. Should have the same dimension as the estimated parameters. See \code{\link{dist}} for details.
#' @param epsilon an optional numeric controlling the stopping criterion. The algorithm terminates when the relative change in the log-likelihoods of two successive iterates is less than \code{epsilon}. The default value is \code{epsilon=1e-8}.
#' @param maxiters an optional number controlling the maximum number of iterations. The default value is \code{maxiters=150}.
#' @param display an optional logical variable controlling the display of iterations. The default value is FALSE.
#'
#' @details See \code{\link{dist}} for details about model parameterization.
#' @return Returns an object of S4 class \code{"MGLMfit"}. An object of class \code{"MGLMfit"} is a list containing at least the following components: \itemize{
#' \item{\code{estimate}}{ the vector of the distribution prameter estimates.}
#' \item{\code{SE}}{ the vector of standard errors of the estimates.}
#' \item{\code{vcov}}{ the variance-covariance matrix of the estimates.}
#' \item{\code{logL}}{ the loglikelihood value.}
#' \item{\code{iter}}{ the number of iterations used.}
#' \item{\code{BIC}}{ Bayesian information criterion.}
#' \item{\code{AIC}}{ Akaike information criterion.}
#' \item{\code{distribution}}{ the distribution fitted.}
#' \item{\code{LRT}}{ when \code{dist="DM"} or \code{"GDM"}, it is the likelihood ratio test statistic for comparing the current model to the multinomial model. No LRT provided when \code{dist="NegMN"}.} 
#' \item{\code{LRTpvalue}}{ the likelihood ratio test P value.}
#' \item{\code{gradient}}{ the gradient at the estimated parameter values.}
#' \item{\code{DoF}}{ the degrees of freedom of the model.}
#' }
#'	
#' @author Yiwen Zhang and Hua Zhou
#'
#' @examples
#' data(rnaseq)
#' Y <- as.matrix(rnaseq[, 1:6])
#' fit <- MGLMfit(data=Y, dist="GDM") 
#' 
#' 
#' @keywords Models Distribution fitting 
#' 
#' 
#' @export 
MGLMfit <- function(data, dist, init, weight, epsilon = 1e-08, maxiters = 150, display = FALSE) {
    
    N <- nrow(data)
    d <- ncol(data)
    
    ## ----------------------------------------## 
    ## Check weights
    ## ----------------------------------------##
    if (!missing(weight) && length(weight) != N) 
        stop("Length of weights doesn't match with the sample size.")
    if (!missing(weight) && any(weight < 0)) 
        stop("Negative weights are not allowed.")
    if (!is.element(dist, c("MN", "DM", "GDM", "NegMN"))) 
        stop(paste("Dist '", dist, "' is not valid. \n", sep = ""))
    
    ##----------------------------------------## 
    ## Give warnings about zero rows
    ##----------------------------------------##
    if (dist != "NegMN") {
        if (any(rowSums(data) == 0)) {
            rmv <- sum(rowSums(data) == 0)
            warning(paste(rmv, " rows are removed because the row sums are 0.", sep = ""))
        }
        if (any(colSums(data) == 0)) {
            rmv <- sum(colSums(data) == 0)
            warning(paste(rmv, " columns are removed because the column sums are 0.", 
                sep = ""))
        }
    }
    
    ##----------------------------------------## 
    ## Fit distribution
    ##----------------------------------------##
    if (dist == "MN") {
        if (missing(weight)) weight <- rep(1, N)
        est <- DMD.MN.fit(data = data, weight = weight)
    } else if (dist == "DM") {
        if (!missing(init) && length(init) != d) 
            stop("Dimension of the initial values is not compatible with the data.")
        est <- DMD.DM.fit(data = data, init = init, weight = weight, epsilon = epsilon, 
            maxiters = maxiters, display = display)
        
    } else if (dist == "GDM") {
        if (!missing(init) && length(init) != 2 * (d - 1)) 
            stop("Dimension of the initial values is not compatible with the data.")
        est <- DMD.GDM.fit(data = data, init = init, weight = weight, epsilon = epsilon, 
            maxiters = maxiters, display = display)
        
    } else if (dist == "NegMN") {
        if (!missing(init) && length(init) != (d + 2)) 
            stop("Dimension of the initial values is not compatible with the data.")
        est <- DMD.NegMN.fit(data = data, init = init, weight = weight, epsilon = epsilon, 
            maxiters = maxiters, display = display)
    }
    
    ##----------------------------------------##
    
    new("MGLMfit", estimate = est$estimate, SE = est$SE, vcov = est$vcov, 
    	logL = est$logL, BIC = est$BIC, AIC = est$AIC, LRT = est$LRT,
    	LRTpvalue = est$LRTpvalue, iter = est$iter, 
    	distribution = est$distribution, gradient = est$gradient, fitted = est$fitted)
}

