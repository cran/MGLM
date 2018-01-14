# A function to tune the optimal penalty threshold, and perform sparse
# regression.  
# Author: Zhang
## ============================================================## 


##============================================================## 
## Tuning function
##============================================================##
#' @title Choose the tuning parameter value in sparse regression
#'
#' @description Finds the tuning parameter value that yields the smallest BIC.
#'
#' @param formula an object of class \code{formula} (or one that can be coerced to that class): a symbolic description of the model to be fitted. The response has to be on the left hand side of ~. 
#' @param data an optional data frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) containing the variables in the model. If not found in \code{data} when using function \code{MGLMtune}, the variables are taken from \code{environment(formula)}, typically the environment from which \code{MGLMtune} is called.
#' @param dist a description of the distribution to fit. See \code{\link{dist}} for the details.
#' @param penalty penalty type for the regularization term. Can be chosen from \code{"sweep"}, \code{"group"}, or \code{"nuclear"}. See \link{MGLMsparsereg} for the description of each penalty type. 
#' @param lambdas an optional vector of the penalty values to tune.  If missing, the vector of penalty values will be set inside the function.  \code{ngridpt} must be provided if \code{lambdas} is missing.
#' @param ngridpt an optional numeric variable specifying the number of grid points to tune.  If \code{lambdas} is given, \code{ngridpt} will be ignored.  Otherwise, the maximum \eqn{\lambda} is determined from the data.  The smallest \eqn{\lambda}is set to \eqn{1/n}, where \eqn{n} is the sample size.
#' @param warm.start an optional logical variable to specify whether to give warm start at each tuning grid point.  If \code{warm.start=TRUE}, the fitted sparse regression coefficients will be used as the initial value when fitting the sparseregression with the next tuning grid. 
#' @param keep.path an optional logical variable controling whether to output the whole solution path. The default is \code{keep.path=FALSE}. If \code{keep.path=TRUE}, the sparse regression result at each grid point will be kept, and saved in the output object \code{select.list}.
#' @param display an optional logical variable to specify whether to show each tuning step.
#' @param weight an optional vector of weights assigned to each row of the data. Should be \code{NULL} or a numeric vector. Could be a variable from the \code{data}, or a variable from \code{environment(formula)} with the length equal to the number of rows of the data. If \code{weight=NULL}, equal weights of ones will be assigned.
#' @param init an optional matrix of initial value of the parameter estimates. Should have the compatible dimension with the data. See \code{\link{dist}} for details of dimensions in each distribution. 
#' @param penidx a logical vector indicating the variables to be penalized. The default value is \code{rep(TRUE, p)}, which means all predictors are subject to regularization. If \code{X} contains intercept, use \code{penidx=c(FALSE,rep(TRUE,p-1))}.
#' @param maxiters an optional numeric controlling the maximum number of iterations. The default value is \code{maxiters=150}.
#' @param ridgedelta an optional numeric controlling the behavior of the Nesterov's accelerated proximal gradient method. The default value is \eqn{\frac{1}{pd}}{1/(pd)}.
#' @param epsilon an optional numeric controlling the stopping criterion. The algorithm terminates when the relative change in the objective values of two successive iterates is less then \code{epsilon}. The default value is \code{epsilon=1e-5}.
#' @param regBeta an optional logical variable used when running negative multinomial regression (\code{dist="NegMN"}). \code{regBeta} controls whether to run regression on the over-dispersion parameter. The default is \code{regBeta=FALSE}.
#' @param overdisp an optional numerical variable used only when fitting sparse negative multinomial model and \code{regBeta=FALSE}.  \code{overdisp} gives the over-dispersion value for all the observations.  The default value is estimated using negative-multinomial regression.  When \code{dist="MN", "DM", "GDM"} or \code{regBeta=TRUE}, the value of \code{overdisp} is ignored.
#'
#' @return \itemize{
#' 	\item{\code{select}}{ the final sparse regression result, using the optimal tuning parameter.}
#'  \item{\code{path}}{ a data frame with degrees of freedom and BICs at each lambda.}
#' }
#'
#' @author Yiwen Zhang and Hua Zhou
#'
#' @seealso \code{\link{MGLMsparsereg}}
#'
#' @examples 
#' set.seed(118)
#' n <- 50
#' p <- 10
#' d <- 5
#' m <- rbinom(n, 100, 0.8)
#' X <- matrix(rnorm(n * p), n, p)
#' alpha <- matrix(0, p, d)
#' alpha[c(1, 3, 5), ] <- 1
#' Alpha <- exp(X %*% alpha)
#' Y <- rdirmn(size=m, alpha=Alpha)
#' sweep <- MGLMtune(Y ~ 0 + X, dist="DM", penalty="sweep", ngridpt=10)
#' show(sweep)
#' 
#' 
#' @export
MGLMtune <- function(formula, data, dist, penalty, lambdas, ngridpt, warm.start = TRUE, 
    keep.path = FALSE, display = FALSE, init, weight, penidx, ridgedelta, maxiters = 150, 
    epsilon = 1e-05, regBeta = FALSE, overdisp) {
    
    call <- match.call()
    if (missing(data)) 
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    mf <- mf[c(1L, match(c("formula", "data", "weight"), names(mf), 0L))]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame(n = 1))
    mt <- attr(mf, "terms")
    Y <- model.response(mf, "any")
    X <- model.matrix(mt, mf, contrasts)
    ow <- getOption("warn")
    
    if (!missing(weight)) 
        weight <- weight[rowSums(Y) != 0]
    X <- as.matrix(X[rowSums(Y) != 0, ])
    Y <- as.matrix(Y[rowSums(Y) != 0, colSums(Y) != 0])
    
    d <- ncol(Y)
    m <- rowSums(Y)
    p <- ncol(X)
    N <- nrow(Y)
    if (missing(weight)) 
        weight <- rep(1, N)
    ## ----------------------------------------## 
    ## Check distribution and d
    ## ----------------------------------------## 
    if (dist == "GDM" && d == 2) 
        stop("When d=2, GDM model is equivilant to DM model, please use dist='DM'.")
    ## ----------------------------------------## 
    ## Pennalty type
    ## ----------------------------------------## 
    if (penalty == "group") 
        penalty <- "group_row"
    if (!penalty %in% c("sweep", "group_row", "nuclear")) 
        stop("Penalty type can only be sweep, group, or nuclear.")
    ## ----------------------------------------##
    ## regularization set
    ## ----------------------------------------##
    if (missing(penidx)) 
        penidx <- rep(TRUE, p)
    ## ----------------------------------------## 
    ## Starting point
    ## ----------------------------------------## 
    if (missing(init)) {
        if (dist == "MN") 
            init <- matrix(0, p, (d - 1)) else if (dist == "DM") 
            init <- matrix(0, p, d) else if (dist == "GDM") 
            init <- matrix(0, p, 2 * (d - 1)) else if (dist == "NegMN") {
            if (regBeta) 
                init <- matrix(0, p, (d + 1)) else init <- matrix(0, p, d)
        }
    }
    if (dist == "NegMN" && regBeta == FALSE && missing(overdisp)) {
        options(warn = -1)
        est <- eval(call("DMD.NegMN.Alpha.reg", Y = Y, init = init, X = X, weight = weight, 
            epsilon = epsilon, maxiters = maxiters, parallel = FALSE, cores = 1, 
            cl = NULL, sys = NULL, display = FALSE))
        overdisp <- est$coefficients$phi
        options(warn = ow)
    } else {
        overdisp <- NULL
    }
    ## ----------------------------------------## 
    ## Ridgedelta
    ## ----------------------------------------## 
    if (missing(ridgedelta)) 
        ridgedelta <- 1/(p * d)
    ## ----------------------------------------## 
    ## Pennalty values
    ## ----------------------------------------## 
    fit.max <- NULL
    if (missing(lambdas)) {
        if (missing(ngridpt)) 
            ngridpt <- 10
        ## ----------------------------------------## 
        ## Find maximum lambda
        ## ----------------------------------------## 
        fit.max <- eval(call("MGLMsparsereg.fit", Y = Y, X = X, dist = dist, lambda = Inf, 
            penalty = penalty, weight = weight, penidx = penidx, init = init, ridgedelta = ridgedelta, 
            maxiters = maxiters, epsilon = epsilon, regBeta = regBeta, overdisp = overdisp))
        
        maxlambda <- fit.max@maxlambda
        # lambdas <- exp(seq(from=log(maxlambda/N), to=log(maxlambda), length.out=15))
        # ----------------------------------------## Find minimum lambda
        # if(penalty=='group_row'){ for(j in 1:15){ if(j==1) B0 <- fit.max$coefficients
        # else B0 <- B_hat temp <- eval(call('MGLMsparsereg.fit', Y=Y, X=X, dist=dist,
        # lambda=lambdas[j], penalty=penalty, weight=weight, penidx=penidx, init=init,
        # ridgedelta=ridgedelta, maxiters=maxiters, epsilon=epsilon, regBeta=regBeta,
        # overdisp=overdisp)) B_hat <- temp$coefficients nz <- sum(rowSums(B_hat^2)>0)
        # if(nz==p){ next }else{ if(j>1) minlambda <- lambdas[j-1] else minlambda <-
        # lambdas[1]/10 break } } }else{
        minlambda <- maxlambda/N
        # }
        
        lambdas <- exp(seq(from = log(maxlambda), to = log(minlambda), length.out = ngridpt))
    } else {
        ngridpt <- length(lambdas)
    }
    
    BICs <- rep(NA, ngridpt)
    AICs <- rep(NA, ngridpt)
    logL <- rep(NA, ngridpt)
    Dof <- rep(NA, ngridpt)
    select.list <- list()
    
    for (j in 1:ngridpt) {
        if (j == 1 & !is.null(fit.max)) {
            temp <- fit.max
            select.list[[j]] <- temp
            B_hat <- temp@coefficients
            BICs[j] <- temp@BIC
            AICs[j] <- temp@AIC
            logL[j] <- temp@logL
            Dof[j] <- temp@Dof
            if (display) 
                print(paste(j, " lamda=", sprintf("%.2f", lambdas[j]), "  BIC=", 
                  sprintf("%.2f", BICs[j]), " AIC=", sprintf("%.2f", AICs[j]), " logL=", 
                  sprintf("%.2f", logL[j]), " Dof=", sprintf("%.2f", Dof[j]), sep = ""))
            next
        }
        
        if (warm.start) {
            if (j == 1) 
                B0 <- init else B0 <- B_hat
        } else {
            B0 <- init
        }
        temp <- eval(call("MGLMsparsereg.fit", Y = Y, X = X, dist = dist, lambda = lambdas[j], 
            penalty = penalty, weight = weight, penidx = penidx, init = B0, ridgedelta = ridgedelta, 
            maxiters = maxiters, epsilon = epsilon, regBeta = regBeta, overdisp = overdisp))
        select.list[[j]] <- temp
        B_hat <- temp@coefficients
        BICs[j] <- temp@BIC
        AICs[j] <- temp@AIC
        logL[j] <- temp@logL
        Dof[j] <- temp@Dof
        
        if (display) 
            print(paste(j, " lamda=", sprintf("%.2f", lambdas[j]), "  BIC=", sprintf("%.2f", 
                BICs[j]), " AIC=", sprintf("%.2f", AICs[j]), " logL=", sprintf("%.2f", 
                logL[j]), " Dof=", sprintf("%.2f", Dof[j]), sep = ""))
        
    }
    chosen.lambda <- lambdas[which.min(BICs)]
    select <- select.list[[which.min(BICs)]]
    select@call <- match.call()
    select@data <- list(Y = Y, X = X)
    select@distribution <- ifelse(dist == "MN", "Multinomial", ifelse(dist == "DM", 
        "Dirichlet Multinomial", ifelse(dist == "GDM", "Generalized Dirichlet Multinomial", 
            "Negative Multinomial")))
    select@penalty <- ifelse(penalty == "group_row", "group", ifelse(penalty == "nuclear", 
              "nuclear", "sweep"))
    select@lambda <- chosen.lambda
  
    outDf <- data.frame(Dof = Dof, Lambda = lambdas, BIC = BICs, AIC = AICs, logL = logL)
    
    
    ## ----------------------------------------## 
    ## Set class
    ## ----------------------------------------## 
    if (keep.path) new("MGLMtune", call = select@call, select = select, path = outDf, select.list = select.list)
    else new("MGLMtune", call = select@call, select = select, path = outDf, select.list = list())
  
} 
