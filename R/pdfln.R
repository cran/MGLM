# Functions to calculate log of pdf of each of the multivariate models.  
# Author: Yiwen Zhang
##============================================================## 

##============================================================## 
## Multinomial model 
##============================================================##
#' @export
dmultn <- function(X, Y, B1, weight) {
    N <- nrow(X)
    d <- ncol(Y)
    batch_sizes <- rowSums(Y)
    B <- cbind(B1, 0)
    alpha <- exp(X %*% B)
    return(sum(Y * (X %*% B) * weight) - sum(Y * log1p(rowSums(alpha)) * weight) + 
        sum(lgamma(batch_sizes + 1) * weight) - sum(lgamma(Y + 1) * weight))
}

#' @name mn 
#' @title The Multinomial Distribution 
#' @description \code{rmn} generates random number vectors given \code{alpha}. 
#' The function \code{rmn(n, size, alpha)} calls \code{rmultinom(n, size, prob)} after converting \code{alpha} to probability. 
#' \code{dmn} computes the log of multinomial probability mass function.
#' 
#' @param Y the multivariate count matrix with dimension \eqn{n \times d}{n x d}, where 
#' \eqn{n = 1,2,\ldots} is number of observations and \eqn{d=2,\ldots} is number of categories.
#' 
#' @param prob the probability parameter of the multinomial distribution.  \code{prob}
#' can be either a vector of length \eqn{d} or a matrix with matching 
#' size of \code{Y}.  If \code{prob} is a vector, it will be replicated \eqn{n} 
#' times to match the dimension of \code{Y}. If the sum(s) of \code{prob} is not 1, it will be automatically scaled to have sum 1.  
#' 
#' @details 
#' A multinomial distribution models the counts of \eqn{d} possible outcomes.
#' The counts of categories are negatively correlated. 
#' \eqn{y=(y_1, \ldots, y_d)} is a \eqn{d} category count vector. 
#' Given the parameter vector \eqn{p = (p_1, \ldots, p_d)}, \eqn{0 < p_j < 1}, 
#' \eqn{\sum_{j=1}^d p_j = 1}{sum_{j=1}^d p_j = 1}, the function calculates the log of the multinomial pmf
#' \deqn{
#'   P(y|p) = C_{y_1, \ldots, y_d}^{m} \prod_{j=1}^{d} p_j^{y_j},
#' }{
#'   P(y|p) = C_{y_1, \ldots, y_d}^{m} prod_{j=1}^{d} p_j^{y_j},
#' }
#' where \eqn{m=\sum_{j=1}^d y_j}{m = sum_{j=1}^d y_j}. Here, \eqn{C_k^n}, often read as "\eqn{n} choose \eqn{k}", 
#' refers the number of \eqn{k} combinations from a set of \eqn{n} elements.
#' 
#' The parameter \eqn{p} can be one vector, like the result from the distribution
#' fitting function; or, \eqn{p} can be a matrix with \eqn{n} rows, like the estimate
#' from the regression function, 
#' \deqn{p_j = \frac{exp(X \beta_j)}{1 + sum_{j'=1}^{d-1} exp(X\beta_{j'})},}{p_j = (exp(X \beta_j)) / (1 + sum_{j'=1}^{d-1} exp(X\beta_{j'})),} where \eqn{j=1,\ldots,d-1}
#' and \eqn{p_d = \frac{1}{1 + \sum_{j'=1}^{d-1} exp(X\beta_{j'})}}{p_d = 1/(1 + sum_{j'=1}^{d-1} exp(X\beta_{j'}))}.
#' The \eqn{d}-th column of the coefficient matrix \eqn{\beta} is set to \eqn{0} to avoid the identifiability issue.
#' 
#' @return The function \code{dmn} returns the value of \eqn{\log(P(y|p))}{logP(y|p)}. 
#' When \code{Y} is a matrix of \eqn{n} rows, the function returns a 
#' vector of length \eqn{n}. 
#' 
#' The function \code{rmn} returns multinomially distributed random number vectors
#' 
#' @author Yiwen Zhang and Hua Zhou
#' @keywords distribution models
#' 
#' @examples m <- 20
#' prob <- c(0.1, 0.2)
#' dm.Y <- rdirmn(n=10, m, prob)	
#' pdfln <- dmn(dm.Y, prob)
NULL
 
#' @rdname mn
#' @export
dmn <- function(Y, prob) {
    
    if (is.vector(prob) && length(prob) > 1) {
        
        if (is.vector(Y) && length(Y) > 1) {
            
            if (length(Y) != length(prob)) {
                stop("size of Y and prob doesn't match.")
            } else {
                Y <- matrix(Y, 1, length(Y))
            }
            
        } else if (is.vector(Y) && length(Y) <= 1) {
            stop("Y can not be a scalar")
        }
        
        prob <- t(matrix(prob, length(prob), dim(Y)[1]))
        
    }
    
    if (any(dim(prob) != dim(Y))) 
        stop("dimensions of prob and Y do not match")
    
    m <- rowSums(Y)
    logl <- lgamma(m + 1) - rowSums(lgamma(Y + 1)) + rowSums(Y * log(prob))
    
    return(logl)
    
}


## ============================================================## 
## Dirichlet multinomial model
## ============================================================## 
## this alpha could be X%*%alpha_hat in a regression estimate
#' @name dirmn
#' @title The Dirichlet Multinomial Distribution
#' @description 
#' \code{ddirmn} computes the log of the Dirichlet multinomial probability mass function.
#' \code{rdirmn} generates Dirichlet multinomially distributed random number vectors. 
#' 
#' @param alpha the parameter of the Dirichlet multinomial distribution. Can be a numerical positive vector or matrix.
#' For \code{ddirmn}, \code{alpha} has to match the size of \code{Y}. If \code{alpha} 
#' is a vector, it will be replicated \eqn{n} times to match the dimension of \code{Y}.
#' 
#' For \code{rdirmn}, if \code{alpha} is a vector, \code{size} must be a scalar, and all the random vectors will
#' be drawn from the same \code{alpha} and \code{size}.  
#' If \code{alpha} is a matrix, the number of rows should match the length of 
#' \code{size}, and each random vector 
#' will be drawn from the corresponding row of \code{alpha} and the corresponding
#' element in the \code{size} vector. See Details below. 
#' @param Y The multivariate count matrix with dimensions \eqn{n \times d}{nxd}, where 
#' \eqn{n = 1,2, \ldots} is the number of observations and \eqn{d=2,3, \ldots} is the number of categories. 
#' 
#' @details 
#' When the multivariate count data exhibits over-dispersion, the traditional 
#' multinomial model is insufficient. Dirichlet multinomial distribution models the
#' probabilities of the categories by a Dirichlet distribution. 
#' Given the parameter vector \eqn{\alpha = (\alpha_1, \ldots, \alpha_d), \alpha_j>0  }, 
#' the probability mass of \eqn{d}-category count vector \eqn{Y=(y_1, \ldots, y_d)}, \eqn{d \ge 2} 
#' under Dirichlet multinomial distribution is
#' \deqn{
#'  P(y|\alpha) = C_{y_1, \ldots, y_d}^{m} \prod_{j=1}^{d} 
#'  \frac{\Gamma(\alpha_j+y_j)}{\Gamma(\alpha_j)}
#'  \frac{\Gamma(\sum_{j'=1}^d \alpha_{j'})}{\Gamma(\sum_{j'=1}^d \alpha_{j'} + \sum_{j'=1}^d y_{j'})},
#'  }{
#'  P(y|\alpha) =
#'  C_{y_1, \ldots, y_d}^{m} prod_{j=1}^d 
#'  {Gamma(\alpha_j+y_j)Gamma(sum_{j'=1}^d \alpha_j')} / {Gamma(\alpha_j)Gamma(sum_{j'=1}^d \alpha_j' + sum_{j'=1}^d y_j')},
#'  }
#' where \eqn{m=\sum_{j=1}^d y_j}{m = sum_{j=1}^d y_j}. Here, \eqn{C_k^n}, often read as "\eqn{n} choose \eqn{k}", 
#' refers the number of \eqn{k} combinations from a set of \eqn{n} elements.
#' 
#'
#' The parameter \eqn{\alpha} can be a vector of length \eqn{d}, 
#' such as the results from the distribution fitting.
#' \eqn{\alpha} can also be a matrix with \eqn{n} rows, such as the inverse link  
#' calculated from the regression parameter estimate \eqn{exp(X\beta)}.
#' 
#' @return For each count vector and each corresponding parameter vector
#' \eqn{\alpha}, the function \code{ddirmn} returns the value \eqn{\log(P(y|\alpha))}{logP(y|\alpha)}. 
#' When \code{Y} is a matrix of \eqn{n} rows, \code{ddirmn} returns a vector of length \eqn{n}.
#'  
#' \code{rdirmn} returns a \eqn{n\times d}{nxd} matrix of the generated random observations.
#' 
#' @examples 
#' m <- 20
#' alpha <- c(0.1, 0.2)
#' dm.Y <- rdirmn(n=10, m, alpha)	
#' pdfln <- ddirmn(dm.Y, alpha)
#' @keywords distribution models
NULL 

#' @rdname dirmn 
#' @export
ddirmn <- function(Y, alpha) {
    
    
    if (is.vector(alpha) && length(alpha) > 1) {
        
        if (is.vector(Y) && length(Y) > 1) {
            if (length(Y) != length(alpha)) {
                stop("size of Y and alpha doesn't match.")
            } else {
                Y <- matrix(Y, 1, length(Y))
            }
            
        } else if (is.vector(Y) && length(Y) <= 1) {
            stop("Y can not be a scalar")
        }
        alpha <- t(matrix(alpha, length(alpha), dim(Y)[1]))
    }
    if (any(dim(alpha) != dim(Y))) {
        stop("dimensions of alpha and Y do not match")
    }
    
    canc <- rowSums(alpha) > 1e+08
    
    ## ----------------------------------------## 
    ## Calculate
    ## ----------------------------------------## 
    alpha_rowsums <- rowSums(alpha)
    m <- rowSums(Y)
    
    logl <- (lgamma(m + 1) + rowSums(lgamma(Y + alpha)) + lgamma(alpha_rowsums)) - 
        (rowSums(lgamma(Y + 1)) + rowSums(lgamma(alpha)) + lgamma(alpha_rowsums + 
        m))
    
    if (sum(canc) > 0) {
        # logMN <- rowSums(Y*log(alpha)) - rowSums(Y*log(alpha_rowsums)) + lgamma(m+1) -
        # rowSums(lgamma(Y+1)) logl[canc] <- logMN[canc]
        logl[canc] <- -Inf
    }
    
    return(logl)
}





## ============================================================## 
## Generalized dirchlet multinomial model
## ============================================================##
#' @name gdirmn
#' @aliases gdirmn dgdirmn rgdirmn 
#' @title The Generalized Dirichlet Multinomial Distribution
#' @description 
#' \code{rgdirmn} generates random observations from the generalized Dirichlet multinomial distribution. 
#' \code{dgdirmn} computes the log of the generalized Dirichlet multinomial probability mass function.
#' 
#' @param Y the multivariate count matrix with dimensions \eqn{n \times d}{nxd}, where 
#' \eqn{n = 1,2, \ldots} is the number of observations and \eqn{d=3,4,\ldots} is the number of categories.
#' 
#' @param alpha the parameter of the generalized Dirichlet multinomial distribution. 
#' \code{alpha} is a numerical positive vector or matrix. 
#' 
#' For \code{gdirmn}, \code{alpha} should match the size of \code{Y}. If \code{alpha} 
#' is a vector, it will be replicated \eqn{n} times to match the dimension of \code{Y}. 
#' 
#' For \code{rdirmn}, if \code{alpha} is a vector, \code{size} must be a scalar.  All the random vectors will
#' be drawn from the same \code{alpha} and \code{size}.  If \code{alpha} is a matrix, the 
#' number of rows should match the length of \code{size}.  Each random vector 
#' will be drawn from the corresponding row of \code{alpha} and the corresponding element of \code{size}. 
#' 
#' @param beta the parameter of the generalized Dirichlet multinomial distribution. \code{beta} should
#' have the same dimension as \code{alpha}.
#' 
#' For \code{rdirm}, if \code{beta} is a vector, \code{size} must be a scalar.  All the random samples will
#' be drawn from the same \code{beta} and \code{size}.  If \code{beta} is a matrix, the 
#' number of rows should match the length of \code{size}.  Each random vector 
#' will be drawn from the corresponding row of \code{beta} and the corresponding element of \code{size}. 
#' 
#' @details 
#' \eqn{Y=(y_1, \ldots, y_d)} are the \eqn{d} category count vectors. Given the parameter vector \eqn{\alpha = (\alpha_1, \ldots, \alpha_{d-1}),
#' \alpha_j>0}, and \eqn{\beta=(\beta_1, \ldots, \beta_{d-1}), \beta_j>0},
#' the generalized Dirichlet multinomial probability mass function is 
#' \deqn{
#'   P(y|\alpha,\beta)
#'   =C_{y_1, \ldots, y_d}^{m} \prod_{j=1}^{d-1} 
#'   \frac{\Gamma(\alpha_j+y_j)}{\Gamma(\alpha_j)}
#'   \frac{\Gamma(\beta_j+z_{j+1})}{\Gamma(\beta_j)}
#'   \frac{\Gamma(\alpha_j+\beta_j)}{\Gamma(\alpha_j+\beta_j+z_j)}  ,
#' }{
#'   P(y|\alpha,\beta)
#'   =C_{y_1, \ldots, y_d}^{m} prod_{j=1}^{d-1} {Gamma(\alpha_j+y_j)Gamma(\beta_j+z_{j+1})Gamma(\alpha_j+\beta_j)} / {Gamma(\alpha_j)Gamma(\beta_j)Gamma(\alpha_j+\beta_j+z_j)},
#' }
#' where \eqn{z_j = \sum_{k=j}^d y_k}{z_j = sum_{k=j}^d y_k} and \eqn{m = \sum_{j=1}^d y_j}{m = sum_{j=1}^d y_j}.
#' Here, \eqn{C_k^n}, often read as "\eqn{n} choose \eqn{k}", 
#' refers the number of \eqn{k} combinations from a set of \eqn{n} elements.
#' 
#' The \eqn{\alpha} and \eqn{\beta} parameters can be vectors, like the results from the 
#' distribution
#' fitting function, or they can be matrices with \eqn{n} rows, 
#' like the estimate
#' from the regression function multiplied by the covariate matrix
#' \eqn{exp(X\alpha)} and \eqn{exp(X\beta)}
#' 
#' @return \code{dgdirmn} returns the value of 
#' \eqn{\log(P(y|\alpha, \beta))}{logP(y|\alpha, \beta)}. 
#' When \code{Y} is a matrix of \eqn{n} rows, the function \code{dgdirmn} returns a vector of length \eqn{n}. 
#' 
#' \code{rgdirmn} returns a \eqn{n\times d}{nxd} matrix of the generated random observations.
#' 
#' @examples 
#' # example 1
#' m <- 20
#' alpha <- c(0.2, 0.5)
#' beta <- c(0.7, 0.4)
#' Y <- rgdirmn(10, m, alpha, beta)
#' dgdirmn(Y, alpha, beta)
#' 
#' # example 2 
#' set.seed(100)
#' alpha <- matrix(abs(rnorm(40)), 10, 4)
#' beta <- matrix(abs(rnorm(40)), 10, 4)
#' size <- rbinom(10, 10, 0.5)
#' GDM.rdm <- rgdirmn(size=size, alpha=alpha, beta=beta)
#' GDM.rdm1 <- rgdirmn(n=20, size=10, alpha=abs(rnorm(4)), beta=abs(rnorm(4)))
#' @keywords distribution models
NULL

#' @rdname gdirmn 
#' @export
dgdirmn <- function(Y, alpha, beta) {
    
    if (is.vector(alpha) && is.vector(beta)) {
        
        if (length(alpha) != length(beta)) 
            stop("the sizes of alpha and beta do not match")
        
        if (is.vector(Y) && length(Y) > 1) {
            
            if (length(Y) != (length(alpha) + 1)) {
                stop("size of Y and alpha doesn't match.")
            } else {
                Y <- matrix(Y, 1, length(Y))
            }
            
        } else if (is.vector(Y) && length(Y) <= 1) {
            stop("Y can not be a scalar")
        }
        
        alpha <- t(matrix(alpha, length(alpha), dim(Y)[1]))
        beta <- t(matrix(beta, length(beta), dim(Y)[1]))
    }
    
    if (nrow(alpha) != nrow(Y) | ncol(alpha) != ncol(Y) - 1) 
        stop("dimensions of the parameters and Y do not match")
    
    ## ----------------------------------------## 
    ## Calculate
    ## ----------------------------------------## 
    
    m <- rowSums(Y)
    d <- ncol(Y)
    
    z <- t(apply(apply(apply(Y, 1, rev), 2, cumsum), 2, rev))
    
    logl <- (lgamma(m + 1) + rowSums(lgamma(Y[, -d] + alpha)) + 
      rowSums(lgamma(z[, -1] + beta)) + rowSums(lgamma(alpha + beta))) - 
      (rowSums(lgamma(Y + 1)) + rowSums(lgamma(alpha)) + rowSums(lgamma(beta)) + 
        rowSums(lgamma(alpha + beta + z[, -d])))
    
    return(logl)
}



## ============================================================## 
## Negative Multinomial model pdf calculation
## ============================================================##
## beta: overdisperson parameter 
#' @name negmn 
#' @title The Negative Multinomial Distribution 
#' @description \code{dnegmn} calculates the log of the negative multinomial probability mass function. 
#' \code{rnegmn} generates random observations from the negative multinomial distribution.
#' @param Y the multivariate response matrix of dimension \eqn{n \times d}{nxd}, 
#' where \eqn{n = 1, 2, \ldots} is number of observations and \eqn{d=2,3,\ldots} is number of categories. 
#' @param prob the probability parameter of the negative multinomial distribution. Should be a numerical non-negative vector or matrix.  
#' 
#' For \code{dnegmn}, \code{prob} can be either
#' a vector of length \eqn{d} \eqn{(d \ge 2)} or a matrix with matching size of \code{Y}.  
#' If \code{prob} is a vector, it will 
#' be replicated \eqn{n} times to match the dimension of \code{Y}.  The sum
#' of each row of \code{prob} should be smaller than 1. 
#' 
#' For \code{rnegmn}, If \code{prob} is a vector, \code{beta}
#' must be a scalar.  All the \code{n} random vectors will be drawn from the
#' same \code{prob} and \code{beta}.  If \code{prob} is a matrix, the number of rows should 
#' match the length of \code{beta}.  Each random vector will be drawn from
#' the corresponding row of \code{prob} and the corresponding element of \code{beta}. Each row of \code{prob} should have sum less than 1.
#' 
#' @param beta the over dispersion parameter of the negative multinomial distribution. \code{beta} can be either a scalar or a vector of length \eqn{n}.
#' 
#' @param alpha an alternative way to specify the probability. Default value is \code{NULL}. See details.
#' 
#' @details 
#' \eqn{y=(y_1, \ldots, y_d)} is a \eqn{d} category vector. Given the parameter vector \eqn{p= (p_1, \ldots, p_d)},
#' \eqn{p_{d+1} = 1/(1 + \sum_{j'=1}^d p_{j'})}{p_{d+1} = 1/(1 + sum_{j'=1}^d p_{j'})}, 
#' and \eqn{\beta}, \eqn{\beta>0}, the negative multinomial probability mass function is 
#' \deqn{
#'   P(y|p,\beta) =  C_{m}^{\beta+m-1}  C_{y_1, \ldots, y_d}^{m} 
#'   \prod_{j=1}^d p_j^{y_j} p_{d+1}^\beta = \frac{\beta_m}{m!}  {m \choose y_1, \ldots, y_d} \prod_{j=1}^d p_j^{y_j} p_{d+1}^\beta,
#' }{
#'   P(y|p,\beta) =  C_{m}^{\beta+m-1}  C_{y_1, \ldots, y_d}^{m} 
#'   prod_{j=1}^d p_j^{y_j} p_{d+1}^\beta = (\beta_m)/(m!) C_{y_1, \ldots, y_d}^{m} prod_{j=1}^d p_j^{y_j} p_{d+1}^\beta,
#' }
#' where \eqn{m = \sum_{j=1}^d y_j}{m = sum_{j=1}^d y_j}. Here, \eqn{C_k^n}, often read as "\eqn{n} choose \eqn{k}", 
#' refers the number of \eqn{k} combinations from a set of \eqn{n} elements.
#' 
#' \code{alpha} is an alternative way to specify the probability: 
#' \deqn{p_j = \frac{\alpha_j}{(1+\sum_{k=1}^{d} \alpha_k)}}{p_j = \alpha_j / (1+sum_{k=1}^{d} \alpha_k)} for \eqn{j=1,\ldots,d} and 
#' \eqn{p_{d+1} = \frac{1}{(1+\sum_{k=1}^{d} \alpha_k)}}{p_{d+1} = 1 / (1+sum_{k=1}^{d} \alpha_k)}. 
#' 
#' The parameter \code{prob} can be a vector and \code{beta} is a scalar; \code{prob} can also
#' be a matrix with \eqn{n} rows, and \code{beta} is a vector of length \eqn{n}
#' like the estimate from the regression function
#' multiplied by the covariate matrix.
#' 
#' @return \code{dnegmn} returns the value of 
#' \eqn{\log(P(y|p, \beta) )}{logP(y|p, \beta)}.  When \code{Y} is a matrix of \eqn{n} rows, the function
#' returns a vector of length \eqn{n}.
#' 
#' \code{rnegmn} returns a \eqn{n\times d}{nxd} matrix of the generated random observations.
#'
#' @examples 
#' ###-----------------------###
#' set.seed(128)
#' n <- 100
#' d <- 4
#' p <- 5
#' a <- -matrix(1,p,d)
#' X <- matrix(runif(n*p), n, p )
#' alpha <- exp(X%*%a)
#' prob <- alpha/(rowSums(alpha)+1)
#' beta <- exp(X%*%matrix(1,p)) 
#' Y <- rnegmn(n, beta, prob)
#' 
#' ###-----------------------###
#' m <- 20
#' n <- 10
#' p <- 5
#' d <- 6
#' a <- -matrix(1,p,d)
#' X <- matrix(runif(n*p), n, p )
#' alpha <- exp(X%*%a)
#' prob <- alpha/(rowSums(alpha)+1)
#' b <- exp(X%*%rep(0.3,p)) 
#' Y <- rnegmn(prob=prob, beta=rep(10, n))
#' dnegmn(Y, b, prob)
#'
#' @author Yiwen Zhang and Hua Zhou
#' @keywords distribution models
NULL
#' @rdname negmn
#' @export
dnegmn <- function(Y, beta, prob = alpha/(rowSums(alpha)+1), alpha = NULL) {
  
  # TO DO: remove the warning after the next version update (current: 0.1.0)
  warning(" note the deprecated argument order;\n dnegmn(Y, prob, beta) and dneg(Y, alpha, beta) from MGLM_0.0.8 have been deprecated;\n use dnegmn(Y, beta, prob = alpha/(rowSums(alpha)+1), alpha=NULL) instead", 
          call. = FALSE)
  
  if (is.null(alpha)) {
    if (is.vector(prob) && length(prob) <= 1) 
      stop("The length of alpha should be larger than 2.")
    
    if (is.vector(prob) && length(prob) > 1) {
      if (length(beta) != 1) 
        stop("The sizes of alpha and beta do not match")
      
      if (is.vector(Y) && length(Y) > 1) {
        
        if (length(Y) != length(prob)) {
          stop("size of Y and alpha doesn't match.")
        } else {
          Y <- matrix(Y, 1, length(Y))
        }
        
      } else if (is.vector(Y) && length(Y) <= 1) {
        stop("Y can not be a scalar")
      }
      
      prob <- matrix(prob, nrow(Y), length(prob), byrow = TRUE)
      beta <- matrix(beta, nrow(Y), 1)
    }
    beta <- matrix(beta, , 1)
    if (nrow(beta) != nrow(Y) | nrow(prob) != nrow(Y)) 
      stop("dimensions of the parameters and Y do not match")
    if (any(rowSums(prob) >= 1)) 
      stop("sum of probability should be smaller than 1")
    ## ----------------------------------------## 
    ## Calculate the log ln
    ## ----------------------------------------## 
    
    m <- rowSums(Y)
    d <- ncol(Y)
    logl <- lgamma(beta + rowSums(Y)) - lgamma(beta) - rowSums(lgamma(Y + 1)) + rowSums(Y * 
                                                                                          log(prob)) + beta * log1p(-rowSums(prob))
    return(logl)
  } else {
    
    if (is.vector(alpha) && length(alpha) <= 1) 
      stop("The length of alpha should be larger than 2.")
    
    if (is.vector(alpha) && length(alpha) > 1) {
      if (length(beta) != 1) 
        stop("The sizes of alpha and beta do not match")
      
      if (is.vector(Y) && length(Y) > 1) {
        
        if (length(Y) != length(alpha)) {
          stop("size of Y and alpha doesn't match.")
        } else {
          Y <- matrix(Y, 1, length(Y))
        }
        
      } else if (is.vector(Y) && length(Y) <= 1) {
        stop("Y can not be a scalar")
      }
      
      alpha <- matrix(alpha, nrow(Y), length(alpha), byrow = TRUE)
      beta <- matrix(beta, nrow(Y), 1)
    }
    beta <- matrix(beta, , 1)
    if (nrow(beta) != nrow(Y) | nrow(alpha) != nrow(Y)) 
      stop("dimensions of the parameters and Y do not match")
    
    ## ----------------------------------------## 
    ## Calculate the log ln
    ## ----------------------------------------## 
    m <- rowSums(Y)
    d <- ncol(Y)
    logl <- lgamma(beta + rowSums(Y)) - lgamma(beta) - rowSums(lgamma(Y + 1)) + rowSums(Y * 
            log(alpha)) - (beta + m) * log1p(rowSums(alpha))
    
    return(logl)
    
  } 
  
}