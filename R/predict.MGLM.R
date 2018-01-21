# Set method to predict 
# Author: Yiwen Zhang
# #' @aliases predict,MGLMreg,ANY-method predict,MGLMreg-method
# #' @usage S4method{predict}{MGLMreg,ANY}(object, newdata)
##============================================================## 

#' @name predict
#' @title Predict method for MGLM Fits 
#' @aliases predict,MGLMreg,ANY-method predict,MGLMreg-method predict,MGLMreg
#' @description Predict using the fitted model from \code{MGLMreg} when given a new set of covariates.
#' @param object model object.
#' @param newdata new covariates data matrix.
#' @return Outputs the probabilities of each category.
#' 
#' This helps answer questions such as whether certain features increase the probability of observing category j.
#' @examples 
#' n <- 200
#' p <- 5
#' d <- 4
#' X <- matrix(runif(p * n), n, p)
#' alpha <- matrix(c(0.6, 0.8, 1), p, d - 1, byrow=TRUE)
#' alpha[c(1, 2),] <- 0
#' Alpha <- exp(X %*% alpha) 
#' beta <- matrix(c(1.2, 1, 0.6), p, d - 1, byrow=TRUE)
#' beta[c(1, 2),] <- 0
#' Beta <- exp(X %*% beta)
#' m <- runif(n, min=0, max=25) + 25
#' Y <- rgdirmn(n, m, Alpha, Beta)
#' gdmReg <- MGLMreg(Y~0+X, dist="GDM")
#' newX <- matrix(runif(1*p), 1, p)
#' pred <- predict(gdmReg, newX)
#' 
#' @importFrom stats predict 
#' @exportMethod predict 
setMethod("predict", signature(object="MGLMreg"), 
          function(object, newdata) {
              beta <- object@coefficients
              d <- ncol(object@data$Y)
              dist <- object@distribution
              
              pred <- switch(dist,
                             "Multinomial"= DMD.MN.predict(beta, d, newdata), 
                             "Dirichlet Multinomial"= DMD.DM.predict(beta, d, newdata),
                             "Generalized Dirichlet Multinomial"= DMD.GDM.predict(beta, d, newdata),
                             "Negative Multinomial" = DMD.NegMN.predict(beta, d, newdata)
              )
              colnames(pred) <- colnames(object@data$Y)
              return(pred) 
})