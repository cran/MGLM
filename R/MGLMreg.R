# Fitting multivariate generalized linear model.  
# Author: Hua Zhou and Yiwen Zhang 
##============================================================## 

##============================================================## 
## Regression function 
##============================================================##
#' @title Fit multivariate response GLM regression
#'
#' @description \code{MGLMreg} fits multivariate response generalized linear models, specified by a symbolic description of the linear predictor and a description of the error distribution.
#'
#' @param formula an object of class \code{formula} (or one that can be coerced to that class): a symbolic description of the model to be fitted. The response has to be on the left hand side of ~.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in \code{data} when using function \code{MGLMreg}, the variables are taken from \code{environment(formula)}, typically the environment from which \code{MGLMreg} is called.
#' @param Y,X for \code{MGLMreg.fit}, \code{X} is a design matrix of dimension \code{n*(p+1)} and \code{Y} is the response matrix of dimension \code{n*d}.
#' @param dist a description of the error distribution to fit. See \code{\link{dist}} for details. 
#' @param weight an optional vector of weights assigned to each row of the data. Should be \code{NULL} or a numeric vector. Could be a variable from \code{data}, or a variable from \code{environment(formula)} with the length equal to the number of rows of the data. If \code{weight=NULL}, equal weights of ones will be assigned. Default is \code{NULL}.
#' @param init an optional matrix of initial value of the parameter estimates. Should have the compatible dimension with \code{data}. See \code{\link{dist}} for details of the dimensions in each distribution.
#' @param epsilon an optional numeric controlling the stopping criterion. The algorithm terminates when the relative change in the loglikelihoods of two successive iterates is less than \code{epsilon}.  The default value is \code{epsilon=1e-8}.
#' @param maxiters an optional numeric controlling the maximum number of iterations. The default value is \code{maxiters=150}.
#' @param display an optional logical variable controlling the display of iterations. The default value is \code{display=FALSE}.
#' @param LRT an optional logical variable controlling whether to perform likelihood ratio test on each predictor. The default value is \code{LRT=FALSE}, in which case only the Wald test is performed.
#' @param parallel an optional logical variable controlling whether to perform parallel computing. On a multi-core Windows machine, a cluster is created based on socket; on a multi-core Linux/Mac machine, a cluster is created based on forking. The default value is \code{parallel=FALSE}.
#' @param cores an optional value specifying the number of cores to use. Default value is half of the logical cores.
#' @param cl a cluster object, created by the package \pkg{parallel} or by package \pkg{snow}. If \code{parallel=TRUE}, use the registered default cluster; if \code{parallel=FALSE}, any given value to \code{cl} will be ignored.
#' @param sys the operating system.  Will be used when choosing parallel type. 	
#' @param regBeta an optional logical variable.  When \code{dist="NegMN"}, the user can decide whether to run regression on the overdispersion parameter \eqn{\beta}.  The default is \code{regBeta=FALSE}.
#' 
#' 
#' @details The formula should be in the form responses ~ covariates where the responses are the multivariate count matrix or a few columns from a data frame which is specified by \code{data}. The covariates are either matrices or from the data frame.  The covariates can be numeric or character or factor. 
#'   See \code{\link{dist}} for details about distributions. 
#'   
#'   Instead of using the formula, the user can directly input the design matrix and the response vector using \code{MGLMreg.fit} function.
#' 
#' @return Returns an object of class \code{"MGLMreg"}. An object of class \code{"MGLMreg"} is a list containing the following components: \itemize{
#'  \item{\code{coefficients}}{ the estimated regression coefficients.}
#'  \item{\code{SE}}{ the standard errors of the estimates.}
#'  \item{\code{Hessian}}{ the Hessian at the estimated parameter values.}
#'  \item{\code{gradient}}{ the gradient at the estimated parameter values.}
#'  \item{\code{wald.value}}{ the Wald statistics.}
#'  \item{\code{wald.p}}{ the p values of Wald test.}
#'  \item{\code{test}}{ test statistic and the corresponding p-value. If \code{LRT=FALSE}, only returns test resultsfrom Wald test; if \code{LRT=TRUE}, returns the test results from both Wald test and likelihood ratio test.}
#'  \item{\code{logL}}{ the final loglikelihood.}
#'  \item{\code{BIC}}{ Bayesian information criterion. }
#'  \item{\code{AIC}}{ Akaike information criterion.}
#'  \item{\code{fitted}}{ the fitted values from the regression model}
#'  \item{\code{iter}}{ the number of iterations used.}
#'  \item{\code{call}}{ the matched call.}
#'  \item{\code{distribution}}{ the distribution fitted.}
#'  \item{\code{data}}{ the data used to fit the model.}
#'  \item{\code{Dof}}{ degrees of freedom.}
#'  }
#'
#' @author Yiwen Zhang and Hua Zhou
#'
#' @seealso See also \code{\link{MGLMfit}} for distribution fitting.
#' 
#' 
#' @examples
#' ##----------------------------------------##
#' ## Generate data
#' n <- 2000
#' p <- 5
#' d <- 4
#' m <- rep(20, n)
#' set.seed(1234)
#' X <- 0.1* matrix(rnorm(n*p),n, p)
#' alpha <- matrix(1, p, d-1)
#' beta <- matrix(1, p, d-1)
#' Alpha <- exp(X %*% alpha)
#' Beta <- exp(X %*% beta)
#' gdm.Y <- rgdirmn(n, m, Alpha, Beta)
#' 
#' ##----------------------------------------##
#' ## Regression
#' gdm.reg <- MGLMreg(gdm.Y~X, dist="GDM", LRT=FALSE)
#'
#' 
#' @keywords models regression  
#' 
#' @export
MGLMreg <- function(formula, data, dist, init = NULL, weight = NULL, epsilon = 1e-08, 
                    maxiters = 150, display = FALSE, LRT = FALSE, 
                    parallel = FALSE, cores = NULL, cl = NULL, sys = NULL, 
                    regBeta = FALSE) {
  
  ##----------------------------------------## 
  ## Creating the environment
  ##----------------------------------------## 
  if (missing(data)) 
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "weight"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame(n = 1))
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  weight <- mf$`(weight)`
  X <- model.matrix(mt, mf, contrasts)
  if (is.null(colnames(Y))) 
    colnames(Y) <- paste("Col", 1:ncol(Y), sep = "_")
  
  ##----------------------------------------## 
  ## Call MGLMreg.fit  
  ##----------------------------------------## 
  est <- eval(call("MGLMreg.fit", Y = Y, X = X, dist = dist, weight = weight, 
                   epsilon = epsilon, maxiters = maxiters, display = display, 
                   LRT = LRT, parallel = parallel, cores = cores, cl = cl, 
                   sys = sys, regBeta = regBeta))
  est@call <- match.call()
  return(est)
  
}

## ============================================================## 
## The MGLM reg function 
## ============================================================##
#' @rdname MGLMreg
#' @export 
MGLMreg.fit <- function(Y, init = NULL, X, dist, weight = NULL, epsilon = 1e-08, 
                        maxiters = 150, display = FALSE, LRT = FALSE, 
                        parallel = FALSE, cores = NULL, cl = NULL, sys = NULL,
                        regBeta = FALSE) {
  
  ow <- getOption("warn")
  
  ##----------------------------------------## 
  ## Give warnings about zero rows
  ##----------------------------------------## 
  if (dist != "NegMN") {
    if (any(rowSums(Y) == 0)) {
      rmv <- rowSums(Y) == 0
      warning(paste(sum(rmv), " rows are removed because the row sums are 0.", 
                    "\n", sep = ""))
      Y <- Y[!rmv, ]
      X <- X[!rmv, ]
    }
    if (any(colSums(Y) == 0)) {
      rmv <- colSums(Y) == 0
      warning(paste(sum(rmv), " columns are removed because the column sums are 0.", 
                    "\n", sep = ""))
      Y <- Y[, !rmv]
    }
  }
  
  ## ----------------------------------------## 
  ## Check dimensions
  N <- nrow(Y)
  d <- ncol(Y)
  p <- ncol(X)
  
  ## ----------------------------------------## 
  ## Check weight length and values
  if (nrow(X) != N) 
    stop("Unequal numbers of observations in X and Y.")
  if (is.null(weight)) 
    weight <- rep(1, N)
  if (!is.null(weight) && length(weight) != N) 
    stop("Length of weights doesn't match with the sample size.")
  if (!is.null(weight) && any(weight < 0)) 
    stop("Negative weights are not allowed.")
   #if (missing(weight)) 
   #  weight <- rep(1, N)
  
  ## ----------------------------------------## 
  ## Check distribution
  if (!is.element(dist, c("MN", "DM", "NegMN", "GDM"))) 
    stop(paste("Dist '", dist, "' is not valid.", sep = ""))
  
  ## ----------------------------------------## 
  ## If parallel with unspecified cores, use half of the logical cores
  if (is.null(cores)) {
    if (parallel) 
      cores <- floor(detectCores()/2) 
  }
  
  ## ----------------------------------------## 
  ## Fit the model
  if (dist != "GDM") {
    if (parallel) {
      cl <- makeCluster(getOption("cl.cores", cores))
      clusterExport(cl, "glm.private")
      sys <- Sys.info()
      sys <- sys[1]
    } else {
      cl <- NULL
      sys <- NULL
    }
    if (dist == "MN") {
      ## ----------------------------------------## 
      ## MN
      if (N < p * (d - 1)) 
        warning(paste("Sample size is smaller than the number of parameters.", 
                      "\n", sep = ""))
      
      if (is.null(init)) {
        init <- matrix(0, p, (d - 1))
        for (i in 1:(d - 1)) {
          init[, i] <- glm.fit(X, Y[, i], family = poisson(link = "log"), 
                               weights = weight)$coefficients
        }
      } else if (any(dim(init) != c(p, (d - 1)))) 
        stop("Dimension of the initial values is not compatible with the data.")
      
      est <- eval(call("DMD.MN.reg", Y = Y, X = X, weight = weight, init = init, 
                       epsilon = epsilon, maxiters = maxiters, display = display, parallel = parallel, 
                       cores = cores, cl = cl, sys = sys))
    } else if (dist == "DM") {
      ## ----------------------------------------## 
      ## DM
      if (N < p * d) 
        warning(paste("Sample size is smaller than the number of parameters.", 
                      "\n", sep = ""))
      
      if (is.null(init)) {
        # alphahat <- as.vector( DMD.DM.fit(data=Y, weight=weight,
        # epsilon=epsilon)$estimate) alphahat <- alphahat/sum(alphahat) init <-
        # rbind(alphahat,matrix(0,(p-1),d))
        init <- matrix(0.1, p, d)
        options(warn = -1)
        for (j in 1:d) {
          fit <- glm.fit(x = X, y = Y[, j]/rowSums(Y), family = binomial(link = "logit"))
          init[, j] <- fit$coefficients
        }
        options(warn = ow)
      } else if (any(dim(init) != c(p, d))) {
        stop("Dimension of the initial values is not compatible with the data.")
      }
      est <- eval(call("DMD.DM.reg", Y = Y, X = X, weight = weight, init = init, 
                       epsilon = epsilon, maxiters = maxiters, display = display, parallel = parallel, 
                       cores = cores, cl = cl, sys = sys))
    } else if (dist == "NegMN") {
      ## ----------------------------------------## 
      ## NegMN
      if (N < p * (d + 1)) 
        warning(paste("Sample size is smaller than the number of parameters.", 
                      "\n", sep = ""))
      
      if (regBeta) {
        if (is.null(init)) {
          alpha_init <- matrix(0, p, d)
          for (i in 1:d) {
            alpha_init[, i] <- glm.fit(X, Y[, i], family = poisson(link = "log"), 
                                       weights = weight)$coefficients
          }
          beta_init <- glm.fit(X, (rowSums(Y) + 10), family = quasipoisson(link = "log"), 
                               weights = weight)$coefficients
          init <- cbind(alpha_init, beta_init)
        } else {
          if (any(dim(init) != c(p, (d + 1)))) 
            stop("Dimension of the initial values is not compatible with the data.") else {
              alpha_init <- init[, 1:d]
              beta_init <- init[, (d + 1)]
            }
        }
        est <- eval(call("DMD.NegMN.reg", Y = Y, init = init, X = X, weight = weight, 
                         epsilon = epsilon, maxiters = maxiters, display = display, parallel = parallel, 
                         cores = cores, cl = cl, sys = sys))
      } else {
        if (is.null(init)) {
          init <- matrix(0, p, d)
          for (i in 1:d) {
            init[, i] <- glm.fit(X, Y[, i], family = poisson(link = "log"), 
                                 weights = weight)$coefficients
          }
        } else {
          if (any(dim(init) != c(p, d))) 
            stop("Dimension of the initial values is not compatible with the data.")
        }
        est <- eval(call("DMD.NegMN.Alpha.reg", Y = Y, init = init, X = X, 
                         weight = weight, epsilon = epsilon, maxiters = maxiters, display = display, 
                         parallel = parallel, cores = cores, cl = cl, sys = sys))
      }
    }
    if (parallel) 
      stopCluster(cl)
    ## ----------------------------------------## 
    ## GDM
  } else if (dist == "GDM") {
    if (d == 2) 
      stop("When d=2, GDM model is equivilant to DM model, please use dist='DM'.")
    if (parallel) {
      cl <- makeCluster(getOption("cl.cores", cores))
      clusterExport(cl, "DMD.DM.reg")
      clusterExport(cl, "ddirmn")
      sys <- Sys.info()
      sys <- sys[1]
    } else {
      cl <- NULL
      sys <- NULL
    }
    
    # if (is.null(init)) {
    #   # init <- as.vector( DMD.GDM.fit(data=Y, weight=weight,
    #   # epsilon=epsilon)$estimate) init <- init/sum(init) init <-
    #   # rbind(rep(init[1:(d-1)],2), matrix(0, (p-1),2*(d-1)))
    #   init <- NULL
    #   
    if (any(dim(init) != c(p, 2 * (d - 1)))) {
      stop("Dimension of the initial values is not compatible with the data")
    } else if (N < p * (d - 1) * 2) 
      warning(paste("Sample size is smaller than the number of parameters.", 
                    "\n", sep = ""))
    
    est <- eval(call("DMD.GDM.reg", Y = Y, X = X, weight = weight, init = init, 
                     epsilon = epsilon, maxiters = maxiters, display = display, parallel = parallel, 
                     cores = cores, cl = cl, sys = sys))
  }
  
  ## ----------------------------------------## 
  ## Clean up the results
  ## ----------------------------------------## 
  wald.value <- c(est$wald.value)
  wald.p <- c(est$wald.p)
  est$test <- cbind(wald.value, wald.p)
  rownames(est$test) <- colnames(X)
  colnames(est$test) <- c("wald value", "Pr(>wald)")
  
  
  ## ----------------------------------------## 
  ## LRT test the hypothesis beta_j=0
  ## ----------------------------------------## 
  logL <- est$logL
  Dof <- numeric()
  if (LRT) {
    options(warn = -1)
    LRT.value <- rep(NA, p)
    LRT.p <- rep(NA, p)
    
    for (t in 1:p) {
      subX <- X[, -t]
      if (dist == "MN") {
        subest <- eval(call("DMD.MN.reg", Y = Y, X = X[, -t], weight = weight, 
                            init = init[-t, ], epsilon = epsilon, maxiters = maxiters, display = FALSE, 
                            parallel = FALSE, cores = cores, cl = cl, sys = sys))
      } else if (dist == "DM") {
        subest <- eval(call("DMD.DM.reg", Y = Y, X = X[, -t], weight = weight, 
                            init = init[-t, ], epsilon = epsilon, maxiters = maxiters, display = FALSE, 
                            parallel = FALSE, cores = cores, cl = cl, sys = sys))
      } else if (dist == "NegMN" & regBeta) {
        subest <- eval(call("DMD.NegMN.reg", Y = Y, X = X[, -t], weight = weight, 
                            init = init[-t, ], epsilon = epsilon, maxiters = maxiters, display = FALSE, 
                            parallel = FALSE, cores = cores, cl = cl, sys = sys))
      } else if (dist == "NegMN" & (!regBeta)) {
        subest <- eval(call("DMD.NegMN.Alpha.reg", Y = Y, X = X[, -t], weight = weight, 
                            init = init[-t, ], epsilon = epsilon, maxiters = maxiters, display = FALSE, 
                            parallel = FALSE, cores = cores, cl = cl, sys = sys))
      } else if (dist == "GDM") {
        subest <- eval(call("DMD.GDM.reg", Y = Y, X = subX, weight = weight, 
                            init = init[-p, ], epsilon = epsilon, maxiters = maxiters, display = FALSE, 
                            parallel = FALSE, cores = cores, cl = cl, sys = sys))
      }
      sublogL <- subest$logL
      LRT.value[t] <- 2 * (logL - sublogL)
    }
    ## ----------------------------------------## 
    ## Calculate the degrees of freedom
    ## ----------------------------------------## 
    Dof <- ifelse(dist == "MN", d - 1, 
           ifelse(dist == "DM", d, 
           ifelse(dist == "GDM", 2 * (d - 1), 
           ifelse(regBeta, d + 1, d))))
    LRT.p <- pchisq(LRT.value, Dof, lower.tail = FALSE)
    test <- cbind(LRT.value, LRT.p)
    colnames(test) <- c("LRT value", "Pr(>LRT)")
    est$test <- cbind(est$test, test)
    options(warn = ow)
  }
  
  ## ----------------------------------------## 
  ## More things to output
  ## ----------------------------------------## 
  est$call <- match.call()
  est$data <- list(Y = Y, X = X)
  
  ## ----------------------------------------## 
  ## Set class
  ## ----------------------------------------## 
  new("MGLMreg", coefficients = est$coefficients, SE = est$SE, Hessian = est$Hessian, 
      gradient = est$gradient, wald.value = est$wald.value, wald.p = est$wald.p, 
      test = est$test, logL = est$logL, BIC = est$BIC, AIC = est$AIC, 
      fitted = est$fitted, call = est$call, distribution = est$distribution, 
      data = est$data, iter = est$iter, Dof = Dof)
  
}



## ============================================================## 
## re-organize the glm.fit function 
## ============================================================##
glm.private <- function(Y, start = NULL, weights, X, family) {
    fit <- glm.fit(x = X, y = Y, weights = weights, family = poisson(link = "log"), 
        start = start)
    return(fit$coefficients)
} 
