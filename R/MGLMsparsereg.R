# Variable selection 
# Author: Zhang
#============================================================## 

##============================================================## 
## Sparse reg function 
##============================================================##
#' @title Fit multivariate GLM sparse regression
#'
#' @description Fit sparse regression in multivariate generalized linear models.
#'
#' @param formula an object of class \code{formula} (or one that can be coerced
#'   to that class): a symbolic description of the model to be fitted.
#'   The response has to be on the left hand side of ~. 
#' @param data an optional data frame, list or environment (or object coercible by 
#'   \code{as.data.frame} to a data frame) containing the variables in the model.
#'   If not found in \code{data} when using function \code{MGLMsparsereg}, the variables 
#'   are taken from \code{environment(formula)}, typically the environment from
#'   which \code{MGLMsparsereg} is called.
#' @param Y a matrix containing the multivariate categorical response data. 
#'   Rows of the matrix represent observations, while columns are the different
#'   categories.  Rows and columns of all zeros are automatically removed when
#'   \code{dist="MN"}, \code{"DM"}, or \code{"GDM"}.  
#' @param X design matrix (including intercept).  
#'    Number of rows of the matrix should match that of \code{Y}.  
#' @param dist a description of the error distribution to fit. See \code{\link{dist}} for details. 
#' @param weight an optional vector of weights assigned to each row of the data. 
#' 	Should be \code{NULL} or a numeric vector. Could be a variable from \code{data}, 
#' 	or a variable from \code{environment(formula)} with the length equal to
#' 	the number of rows of the data.
#' 	If \code{weight=NULL}, equal weights of ones will be assigned.
#' @param init an optional matrix of initial value of the parameter estimates.
#' 	Should have the compatible dimension with the data. See \code{\link{dist}} for
#' 	details of the dimensions in each distribution. 
#' @param lambda penalty parameter. 
#' @param penalty penalty type for the regularization term. Can be chosen from \code{"sweep"}, 
#' 	\code{"group"}, or \code{"nuclear"}. See Details for the description of each penalty type.
#' @param penidx a logical vector indicating the variables to be penalized. The default value is \code{rep(TRUE, p)}, which means all predictors are subject to regularization. If \code{X} contains intercept, use \code{penidx=c(FALSE,rep(TRUE,p-1))}.
#' @param maxiters an optional numeric controlling the maximum number of iterations. The default
#' 	value is maxiters=150.
#' @param ridgedelta an optional numeric controlling the behavior of the Nesterov's accelerated proximal gradient method. The default value is \eqn{\frac{1}{pd}}{1/pd}.
#' @param epsilon an optional numeric controlling the stopping criterion. The algorithm terminates when the relative change in the objective values of two successive iterates is less then \code{epsilon}.
#'   	The default value is \code{epsilon=1e-5}.
#' @param regBeta an optional logical variable used when running negative multinomial regression (\code{dist="NegMN"}).  
#' 	\code{regBeta} controls whether to run regression on the over-dispersion parameter.
#' 	The default is \code{regBeta=FALSE}.
#' @param overdisp an optional numerical variable used only when fitting sparse negative multinomial 
#' 	model \code{dist="NegMN"} and \code{regBeta=FALSE}.  \code{overdisp} gives the over dispersion value
#' 	for all the observations.  The default value is estimated using negative-multinomial regression.  When \code{dist="MN", "DM", "GDM"} or \code{regBeta=TRUE}, the value of \code{overdisp} is ignored. 
#'
#' 	
#' @details 
#' In general, we consider regularization problem
#' \deqn{
#' \min_B h(B) = -l(B)+ J(B),
#' }{
#' min_B h(B) = -l(B)+ J(B),
#' }
#' where \eqn{l(B)} is the loglikelihood function and \eqn{J(B)} is the 
#' regularization function.  
#' 
#' Sparsity in the individual elements of the parameter matrix \eqn{B} is achieved 
#' by the lasso penalty (\code{dist="sweep"})
#' \deqn{
#' J(B) = \lambda \sum_{k\in penidx} \sum_{j=1}^d \|B_{kj}\|
#' }{
#' J(B) = \lambda \sum_{k in penidx} \sum_{j=1}^d ||B_{kj}||
#' }
#' 
#' Sparsity in the rows of the regression parameter matrix \eqn{B} is achieved
#' by the group penalty (\code{dist="group"})
#' \deqn{
#' J(B) = \lambda \sum_{k \in penidx} \|B_{k \cdot}\|_2,
#' }{
#' J(B) = \lambda \sum_{k in penidx} ||B_{k cdot}||_2,
#' }
#' where \eqn{\|v\|_2}{||v||_2} is the \eqn{l_2}{L2} norm of a vector \eqn{v}. In other words, 
#' \eqn{\|B_{k\cdot}\|_2}{||B_{k cdot}||_2} is the \eqn{l_2}{L2} norm of the \eqn{k}-th row of the 
#' parameter matrix \eqn{B}.
#' 
#' Sparsity in the rank of the parameter matrix \eqn{B} is achieved by the nuclear norm penalty (\code{dist="nuclear"})
#' \deqn{
#' J(B) = \lambda \|B\|_*= \lambda \sum_{i=1}^{min(p, d)} \sigma_i(B),
#' }{
#' J(B) = \lambda ||B||_{*}= \lambda \sum_{i=1}^{min(p, d)} \sigma_i(B),
#' }
#' where \eqn{\sigma_i(B)} are the singular values of the parameter matrix \eqn{B}. 
#' The nuclear norm \eqn{\|B\|_*}{||B||_{*}} is a convex relaxation of \eqn{rank(B)=\|\sigma(B)\|_0}{rank(B)=||\sigma(B)||_0}.
#' 
#' See \code{\link{dist}} for details about distributions. 
#' 
#' 
#' @return Returns an object of class \code{"MGLMsparsereg"}. An object of class \code{"MGLMsparsereg"} is a list containing at least the following components:  \itemize{
#' \item{\code{coefficients}}{ the estimated matrix of regression coefficients.}
#' \item{\code{logL}}{ the final loglikelihood value.}
#' \item{\code{AIC}}{ Akaike information criterion.}
#' \item{\code{BIC}}{ Bayesian information criterion.}
#' \item{\code{Dof}}{ degrees of freedom of the estimated parameter.}
#' \item{\code{iter}}{ number of iterations used. }
#' \item{\code{maxlambda}}{ the maximum tuning parameter such that 
#' 	the estimated coefficients are not all zero.  This value is returned only
#' 	when the tuning parameter \code{lambda} given to the function is large enough 
#' 	such that all the parameter estimates are zero; otherwise, \code{maxlambda}
#' 	is not computed.}
#' \item{\code{call}}{ a matched call.}
#' \item{\code{data}}{ the data used to fit the model: a list of the predictor matrix
#' and the response matrix.}
#' \item{\code{penalty}}{ the penalty chosen when running the penalized regression.}
#' }
#'
#' @author Yiwen Zhang and Hua Zhou
#'
#' @examples 
#' ## Generate Dirichlet Multinomial data
#' dist <- "DM"
#' n <- 100
#' p <- 15
#' d <- 5
#' m <- runif(n, min=0, max=25) + 25
#' set.seed(134)
#' X <- matrix(rnorm(n*p),n, p)
#' alpha <- matrix(0, p, d)
#' alpha[c(1,3, 5), ] <- 1
#' Alpha <- exp(X\%*\%alpha)
#' Y <- rdirmn(size=m, alpha=Alpha)
#' 
#' ## Tuning
#' ngridpt <- 10
#' p <- ncol(X)
#' d <- ncol(Y)
#' pen <- 'nuclear'
#' spfit <- MGLMsparsereg(formula=Y~0+X, dist=dist, lambda=Inf, penalty=pen)
#' 
#' @keywords regression 
#'
#' @export
MGLMsparsereg <- function(formula, data, dist, lambda, penalty, weight, init, penidx, 
    maxiters = 150, ridgedelta, epsilon = 1e-05, regBeta = FALSE, overdisp) {
    call <- match.call()
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
    X <- model.matrix(mt, mf, contrasts)
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
    if (dist == "GDM" && d == 2) 
        stop("When d=2, GDM model is equivilant to DM model, please use dist='DM'.")
    if (penalty == "group") 
        penalty <- "group_row"
    if (!penalty %in% c("sweep", "group_row", "nuclear")) 
        stop("Penalty type can only be sweep, group, or nuclear.")
    if (missing(penidx)) 
        penidx <- rep(TRUE, p)
    if (missing(init)) {
        if (dist == "MN") {
            init <- matrix(0, p, (d - 1))
        } else if (dist == "DM") {
            init <- matrix(0, p, d)
        } else if (dist == "GDM") {
            init <- matrix(0, p, 2 * (d - 1))
        } else if (dist == "NegMN") {
            if (regBeta) 
                init <- matrix(0, p, (d + 1)) else init <- matrix(0, p, d)
        }
    }
    if (dist == "NegMN" && regBeta == FALSE && missing(overdisp)) {
        est <- DMD.NegMN.fit(Y)
        overdisp <- est$estimate[d + 1]
    } else {
        overdisp <- NULL
    }
    if (missing(ridgedelta)) 
        ridgedelta <- 1/(p * d)
    est <- eval(call("MGLMsparsereg.fit", Y = Y, X = X, dist = dist, lambda = lambda, 
        penalty = penalty, weight = weight, init = init, penidx = penidx, maxiters = maxiters, 
        ridgedelta = ridgedelta, epsilon = epsilon, regBeta = regBeta, overdisp = overdisp))
    est@call <- match.call()
    
    ## ----------------------------------------## 
    ## Set class
    ## ----------------------------------------## 
    return(est) 

}

## ============================================================## 
## The MGLM sparse reg function 
## ============================================================##
#' @rdname MGLMsparsereg
#' @export 
MGLMsparsereg.fit <- function(Y, X, dist, lambda, penalty, weight, init, penidx, 
    maxiters = 150, ridgedelta, epsilon = 1e-05, regBeta = FALSE, overdisp) {
    
    tol <- 1e-8 
    d <- ncol(Y)
    m <- rowSums(Y)
    p <- ncol(X)
    N <- nrow(Y)
    
    if (dist == "GDM") {
        colOrder <- order(colSums(Y), decreasing = TRUE)
        Y <- Y[, colOrder]
        outOrder <- order(colOrder[-d])
        # Ys <- t( apply(apply( apply(Y,1,rev),2,cumsum),2,rev) )
    }
    
    beta_old <- init
    B <- init
    alpha <- 1
    alpha_iter <- list()
    objval <- Inf
    niter <- 1
    isdescent <- TRUE
    oriRidge <- ridgedelta
    while ((niter < 3) || (!stop)) {
        niter <- niter + 1
        beta_old <- B
        obj1 <- objval
        if (niter <= 2) 
            S <- init
        loss <- MGLM.loss(Y, X, S, dist, weight, regBeta, overdisp)
        loss.S <- loss[[1]]
        loss.D1S <- loss[[2]]
        for (l in 1:50) {
            A <- S - ridgedelta * loss.D1S
            B <- matrix(NA, nrow(B), ncol(B))
            B[!penidx, ] <- A[!penidx, ]
            if (d > 2) {
                Apen <- A[penidx, ]
            } else if (d == 2) 
                Apen <- matrix(A[penidx, ], , 1)
            pen <- matrix_threshold(X = Apen, lambda = ridgedelta * lambda, penalty = penalty)
            B[penidx, ] <- pen[[1]]
            penval <- pen[[2]]
            if (all(abs(B[penidx, ]) < tol)) {
                if (penalty == "sweep") {
                  maxlambda <- max(abs(Apen))/ridgedelta
                } else if (penalty == "group_row") {
                  maxlambda <- max(sqrt(rowSums(Apen^2)))/ridgedelta
                } else if (penalty == "nuclear") {
                  sing.vec <- svd(Apen)$d
                  maxlambda <- sing.vec[1]/ridgedelta
                }
                if (lambda > maxlambda) {
                  B[penidx, ] <- 0
                  penval <- 0
                }
            } else {
                maxlambda <- numeric() # 12/1/17
            }
            objloss <- MGLM.loss(Y, X, B, dist, weight, regBeta, overdisp)
            loss.S2 <- objloss[[1]]
            loss.D1S2 <- objloss[[2]]
            objval <- loss.S2 + penval
            BminusS <- B - S
            surval <- loss.S + sum(loss.D1S * BminusS) + norm(BminusS, type = "F")^2/2/ridgedelta + 
                penval
            if (!is.na(objval) && !is.na(surval)) 
                if (objval <= surval) 
                  break else ridgedelta <- ridgedelta/2
        }
        alpha_old <- alpha
        alpha <- (1 + sqrt(4 + alpha_old^2))/2
        if (!is.na(objval) & objval <= obj1) {
            stop <- abs(obj1 - objval) < epsilon * (abs(obj1) + 1)
            S <- B + (alpha_old - 1)/alpha * (B - beta_old)
        } else {
            objval <- obj1
            if (isdescent) {
                isdescent <- FALSE
                stop <- FALSE
                S <- B + (alpha_old - 1)/alpha * (beta_old - B)
                B <- beta_old
            } else {
                stop <- TRUE
            }
        }
        if (niter >= maxiters) 
            stop <- TRUE
    }
    if (d > 2) 
        Bpen <- B[penidx, ] else if (d == 2) 
        Bpen <- matrix(B[penidx, ], , 1)
    if (penalty == "sweep") {
        Dof <- sum(!penidx) * ncol(B) + sum(Bpen != 0)
    } else if (penalty == "group_row") {
        Dof <- sum(!penidx) * ncol(B) + sum(rowSums(Bpen^2) != 0) + sum((ncol(B) - 
            1) * rowSums(Bpen^2)/rowSums(Apen^2))
    } else if (penalty == "nuclear") {
        Aspectrum <- svd(A[penidx, ])$d
        if (sum(Aspectrum > ridgedelta * lambda) == 0) {
            Dof <- 0
        } else {
            Dof <- 0
            for (i in 1:sum(Aspectrum > ridgedelta * lambda)) {
                Dof <- Dof + 1 + 2 * sum(Aspectrum[i] * (Aspectrum[i] - ridgedelta * 
                  lambda)/(Aspectrum[i]^2 - Aspectrum[-i]^2)) + abs(p - d) * (1 - 
                  ridgedelta * lambda/Aspectrum[i])
            }
            Dof <- Dof + sum(!penidx) * d
        }
    }
    if (dist == "GDM") 
        Dof <- Dof/2 else Dof <- Dof
    logL <- -MGLM.loss(Y, X, B, dist, weight, regBeta, overdisp)[[1]]
    AIC <- -2 * logL + 2 * Dof
    BIC <- -2 * logL + log(N) * Dof
    
    if (dist == "GDM") {
        B <- B[, c(outOrder, outOrder + d - 1)]
    }
    
 
    distribution <- ifelse(dist == "MN", "Multinomial", ifelse(dist == "DM", 
       "Dirichlet Multinomial", ifelse(dist == "GDM", "Generalized Dirichlet Multinomial", 
                                       "Negative Multinomial")))
    penalty_name <- ifelse(penalty == "group_row", "group", ifelse(penalty == "nuclear", 
                                                                  "nuclear", "sweep"))
    new("MGLMsparsereg", call = match.call(), data = list(Y = Y, X = X), 
        coefficients = B, logL = logL, BIC = BIC, AIC = AIC, Dof = Dof, 
        iter = niter, maxlambda = maxlambda, lambda = lambda, 
        distribution = distribution, penalty = penalty_name) # , Beta = est$Beta) 
}

## ============================================================================##
## MGLM loss
## ============================================================================##

MGLM.loss <- function(Y, X, beta, dist, weight, regBeta = FALSE, Beta) {
	
	 if (missing(weight)) 
        weight <- rep(1, nrow(Y))
        
     losslist <- switch(dist,
     					"MN"= DMD.MN.loss(Y, X, beta, dist, weight),
     					"DM" = DMD.DM.loss(Y, X, beta, dist, weight),
     					"GDM" = DMD.GDM.loss(Y, X, beta, dist, weight),
     					"NegMN" = DMD.NegMN.loss(Y, X, beta, dist, weight, regBeta = FALSE, Beta))
	
     return(losslist)
}

## ============================================================================##
## matrix thresholding
## ============================================================================##

matrix_threshold <- function(X, lambda, penalty) {
    
    N <- nrow(X)
    d <- ncol(X)
    B <- matrix(0, N, d)
    if (penalty == "sweep") {
        B <- lsq_thresholding(X, lambda)
        penalty_value <- lambda * sum(abs(B))
    } else if (penalty == "group_row" || penalty == "group") {
        row_12norm <- sqrt(rowSums(X^2))
        vec <- 1 - lambda/row_12norm
        vec[vec < 0] <- 0
        B <- X * vec
        penalty_value <- lambda * sum(sqrt(rowSums(B^2)))
    } else if (penalty == "group_col") {
        row_12norm <- sqrt(colSums(X^2))
        B <- X * max(c(1 - lambda/row_12norm, 0))
        penalty_value <- lambda * sum(sqrt(colSums(B^2)))
    } else if (penalty == "nuclear") {
        decomp <- svt(X, lambda)
        U <- decomp[[1]]
        s <- as.vector(decomp[[2]])
        V <- decomp[[3]]
        if (length(s) == 0) 
            B <- matrix(0, N, d) else B <- U %*% (s * t(V))
        bs <- svd(B)$d
        penalty_value <- lambda * sum(bs)
    }
    
    return(list(B, penalty_value))
}

## ============================================================================##
## lsq_threshold We can only work on lasso for now
## ============================================================================##
lsq_thresholding <- function(b, lambda) {
    if (lambda < 0) 
        stop("Penalty constant lambda should be nonnegative.")
    
    B <- b
    B[abs(b) <= lambda] <- 0
    B[b > lambda] <- B[b > lambda] - lambda
    B[b < -lambda] <- B[b < -lambda] + lambda
    
    return(B)
}

## ============================================================================##
## svt
## ============================================================================##
svt <- function(b, lambda) {
    decomp <- svd(b)
    s <- decomp$d
    if (lambda > 0) {
        s <- lsq_thresholding(as.matrix(s), lambda)
        idx <- s > 0
        s <- s[idx]
        U <- decomp$u[, idx]
        V <- decomp$v[, idx]
    }
    return(list(U = U, s = s, V = V))
} 
