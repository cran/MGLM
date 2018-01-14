
## ============================================================## 
## Function: fit GDM Reg ; called by MGLMreg 
## ============================================================##

DMD.GDM.reg <- function(Y, init, X, weight, epsilon, maxiters, display, parallel, 
                        cores, cl, sys) {
  
  ## ----------------------------------------## 
  ## Keep some original values
  ## ----------------------------------------## 
  pred <- matrix(NA, nrow(Y), ncol(Y))
  emptyRow <- rowSums(Y) == 0
  Xf <- X
  ow <- getOption("warn")
  
  weight <- weight[rowSums(Y) != 0]
  X <- as.matrix(X[rowSums(Y) != 0, ])
  Y <- as.matrix(Y[rowSums(Y) != 0, colSums(Y) != 0])
  d <- ncol(Y)
  p <- ncol(X)
  m <- rowSums(Y)
  N <- nrow(Y)
  
  colOrder <- order(colSums(Y), decreasing = TRUE)
  Y <- Y[, colOrder]
  outOrder <- order(colOrder[-d])
  Ys <- t(apply(apply(apply(Y, 1, rev), 2, cumsum), 2, rev))
  
  if (is.null(init)) {
    init <- matrix(0, p, 2 * (d - 1))
    options(warn = -1)
    for (j in 1:(d - 1)) {
      den <- Ys[, j]
      den[den == 0] <- m[den == 0]
      fit <- glm.fit(X, Y[, j]/den, family = binomial(link = "logit"))
      init[, c(j, j + d - 1)] <- fit$coefficients
    }
    options(warn = ow)
  }
  
  alpha <- beta <- alpha_se <- beta_se <- matrix(0, p, (d - 1))
  gradients <- matrix(0, 2 * p, (d - 1))
  niter <- rep(0, (d - 1))
  options(warn = -1)
  if (!parallel) {
    for (j in 1:(d - 1)) {
      YGDM <- cbind(Y[, j], Ys[, (j + 1)])
      reg.result <- eval(call("DMD.DM.reg", Y = YGDM, init = init[, c(j, j + 
                                                                        d - 1)], X = X, weight = weight, epsilon = epsilon, maxiters = maxiters, 
                              display = display, parallel = FALSE, cores = NULL, cl = NULL, sys = NULL))
      alpha[, j] <- reg.result$coefficients[, 1]
      beta[, j] <- reg.result$coefficients[, 2]
      alpha_se[, j] <- reg.result$SE[, 1]
      beta_se[, j] <- reg.result$SE[, 2]
      niter[j] <- reg.result$iter
      gradients[, j] <- reg.result$gradient
    }
  } else {
    Y.list <- list()
    init.list <- list()
    for (i in 1:(d - 1)) {
      Y.list[[i]] <- cbind(Y[, i], Ys[, (i + 1)])
      init.list[[i]] <- init[, c(i, i + d - 1)]
    }
    if (sys[1] == "Windows") {
      fit.list <- clusterMap(cl, DMD.DM.reg, Y.list, init.list, .scheduling = "dynamic", 
                             MoreArgs = list(X = X, weight = weight, epsilon = epsilon, maxiters = maxiters, 
                                             display = display, parallel = FALSE))
    } else if (sys[1] != "Windows") {
      fit.list <- mcmapply(cl, DMD.DM.reg, Y.list, init.list, MoreArgs = list(X = X, 
                                                                              weight = weight, epsilon = epsilon, maxiters = maxiters, display = display, 
                                                                              parallel = FALSE), mc.cores = cores, mc.preschedule = FALSE)
    }
    
    for (j in 1:(d - 1)) {
      alpha[, j] <- fit.list[[j]]$coefficients[, 1]
      beta[, j] <- fit.list[[j]]$coefficients[, 2]
      alpha_se[, j] <- fit.list[[j]]$SE[, 1]
      beta_se[, j] <- fit.list[[j]]$SE[, 2]
      niter[j] <- fit.list[[j]]$iter
      gradients[, j] <- fit.list[[j]]$gradient
    }
  }
  options(warn = ow)
  ## ----------------------------------------## 
  ## End of the main loop
  ## ----------------------------------------## 
  
  ll2 <- sum(weight * dgdirmn(Y, exp(X %*% alpha), exp(X %*% beta)), na.rm = TRUE)
  
  ## ----------------------------------------## 
  ## Compute output statistics
  ## ----------------------------------------## 
  BIC <- -2 * ll2 + log(N) * p * (d - 1) * 2
  AIC <- -2 * ll2 + 2 * p * (d - 1) * 2
  gradients <- cbind(gradients[1:p, ], gradients[(p + 1):(2 * p), ])
  wald <- rep(NA, p)
  wald.p <- rep(NA, p)
  SE <- matrix(Inf, p, 2 * (d - 1))
  H <- matrix(NA, 2 * (d - 1), 2 * (d - 1))
  
  ## ----------------------------------------## 
  ## Check the gradients
  ## ----------------------------------------##     
  if (mean(gradients^2) > 1e-04) {
    warning(paste("The algorithm doesn't converge within", sum(niter), "iterations. The norm of the gradient is ", 
                  sum(gradients^2, na.rm = T), 
                  " Please interpret hessian matrix and MLE with caution.", 
                  sep = " "))
  }
  
  ## ----------------------------------------## 
  ## Check whether H is negative definite
  ## ----------------------------------------## 
  B1 <- exp(X %*% alpha)
  B2 <- exp(X %*% beta)
  B12 <- B1 + B2 
  a1 <- B1 * (digamma(B1 + Y[, -d]) - digamma(B1)) - B1^2 * (-trigamma(B1 + Y[, 
                                                                              -d]) + trigamma(B1))
  a2 <- B1 * (digamma(B12 + Ys[, -d]) - digamma(B12)) - B1^2 * (-trigamma(B12 + 
                                                                            Ys[, -d]) + trigamma(B12))
  b <- B1 * B2 * (-trigamma(B12 + Ys[, -d]) + trigamma(B12))
  d1 <- B2 * (digamma(B2 + Ys[, -1]) - digamma(B2)) - B2^2 * (-trigamma(B2 + Ys[, 
                                                                                -1]) + trigamma(B2))
  d2 <- B2 * (digamma(B12 + Ys[, -d]) - digamma(B12)) - B2^2 * (-trigamma(B12 + 
                                                                            Ys[, -d]) + trigamma(B12))
  
  if (any(is.nan(a1), is.nan(a2), is.nan(b), is.nan(d1), is.nan(d2))) {
    warning(paste("Out of range of trigamma. \n\t\t\tNo standard error or tests results reported.\n\t\t\tRegression parameters diverge. Recommend multinomial logit model.", 
                  "\n", sep = ""))
    SE <- matrix(NA, p, 2 * (d - 1))
    wald <- rep(NA, p)
    wald.p <- rep(NA, p)
    H <- matrix(NA, p * 2 * (d - 1), p * 2 * (d - 1))
  } else {
    Ha <- matrix(0, p * (d - 1), p * (d - 1))
    for (i in 1:(d - 1)) {
      id <- (i - 1) * p + c(1:p)
      Ha[id, id] <- t(X) %*% (weight * (a1[, i] - a2[, i]) * X)
    }
    Hb <- matrix(0, p * (d - 1), p * (d - 1))
    for (i in 1:(d - 1)) {
      id <- (i - 1) * p + c(1:p)
      Hb[id, id] <- t(X) %*% (weight * b[, i] * X)
    }
    Hd <- matrix(0, p * (d - 1), p * (d - 1))
    for (i in 1:(d - 1)) {
      id <- (i - 1) * p + c(1:p)
      Hd[id, id] <- t(X) %*% (weight * (d1[, i] - d2[, i]) * X)
    }
    H <- rbind(cbind(Ha, Hb), cbind(Hb, Hd))
    eig <- eigen(H)$values
    
    if (any(is.complex(eig)) || any(eig > 0)) {
      warning(paste("The estimate is a saddle point.\n\t\t\tNo standard error estimate or test results reported.", 
                    sep = ""))
    } else if (any(eig == 0)) {
      warning(paste("The hessian matrix is almost singular.\n\t\t\tNo standard error estimate or test results reported.", 
                    sep = ""))
    } else if (all(eig < 0)) {
      ## ----------------------------------------------## 
      ## Calculate the Wald statistic
      ## ----------------------------------------------## 
      ## Hinv <- chol2inv(chol(-H) )
      Hinv <- solve(-H)
      SE <- matrix(sqrt(diag(Hinv)), p, 2 * (d - 1))
      for (j in 1:p) {
        id <- c(0:(2 * (d - 1) - 1)) * p + j
        invH <- NULL
        try(invH <- solve(Hinv[id, id], c(alpha[j, ], beta[j, ])), silent = TRUE)
        if (is.numeric(invH)) {
          wald[j] <- sum(c(alpha[j, ], beta[j, ]) * invH)
          wald.p[j] <- pchisq(wald[j], 2 * (d - 1), lower.tail = FALSE)
        }
      }
    }
  }
  
  tmpalpha1 <- exp(Xf %*% alpha[, 1])
  tmpalpha2 <- exp(Xf %*% alpha[, 2])
  
  pred[, 1] <- tmpalpha1 / (tmpalpha1 + exp(Xf %*% beta[, 1]))
  pred[, 2] <- (1 - pred[, 1]) * tmpalpha2 / (tmpalpha2 + exp(Xf %*% beta[, 2]))
  if (d > 3) {
    for (f in 3:(d - 1)) {
      tmpalphaf <- exp(Xf %*% alpha[, f])
      pred[, f] <- (1 - rowSums(pred[, 1:(f - 1)])) * tmpalphaf / (tmpalphaf +
                                                                     exp(Xf %*% beta[, f]))
    }
  } else if (d == 3) {
    pred[, d] <- (1 - rowSums(pred[, 1:(d - 1)]))
  }
  pred[emptyRow, ] <- 0
  
  ## ----------------------------------------## 
  ## Clean up the results
  ## ----------------------------------------## 
  coefficients = cbind(alpha[, outOrder], beta[, outOrder])
  SE = SE[, c(outOrder, outOrder + d - 1)]
  
  colnames(coefficients) <- c(paste("alpha_", colnames(Y)[1:(d - 1)], sep = ""), 
                              paste("beta_", colnames(Y)[1:(d - 1)], sep = ""))
  colnames(SE) <- c(paste("alpha", colnames(Y)[1:(d - 1)], sep = ""), 
                    paste("beta_", colnames(Y)[1:(d - 1)], sep = ""))
  colnames(gradients) <- c(paste("alpha_", colnames(Y)[1:(d - 1)], sep = ""), 
                           paste("beta_", colnames(Y)[1:(d - 1)], sep = ""))
  
  ## ----------------------------------------## 
  list(coefficients = coefficients, SE = SE, 
       Hessian = H[c(outOrder, outOrder + d - 1), c(outOrder, outOrder + d - 1)], 
       wald.value = wald, wald.p = wald.p, logL = ll2, BIC = BIC, AIC = AIC, 
       iter = sum(niter), DoF = 2 * p * (d - 1), gradient = gradients, 
       fitted = pred[, order(colOrder)], distribution = "Generalized Dirichlet Multinomial")
  
}


##============================================================## 
## Function: fit Generalized DM ; called by MGLMfit 
##============================================================##
#' @rdname MGLMfit 
DMD.GDM.fit <- function(data, init, weight, epsilon = 1e-08, maxiters = 150, 
                        display = FALSE) {
  
  ow <- getOption("warn")
  ##----------------------------------------## 
  ## Remove the zero rows and columns
  ##----------------------------------------## 
  if (!missing(weight)) 
    weight <- weight[rowSums(data) != 0]
  
  data <- data[rowSums(data) != 0, colSums(data) != 0]
  N <- nrow(data)
  d <- ncol(data)
  m <- rowSums(data)
  max <- max(m)
  k <- c(0:(max - 1))
  Y <- t(apply(apply(apply(data, 1, rev), 2, cumsum), 2, rev))
  Y2 <- as.matrix(Y[, (2:d)])
  
  ##----------------------------------------## 
  ## Give default value to the missing variables
  ##----------------------------------------## 
  if (missing(weight)) 
    weight <- rep(1, N)
  if (missing(init)) {
    y <- as.matrix(Y[apply(Y, 1, min) != 0, 1:(d - 1)])
    y2 <- as.matrix(Y[apply(Y, 1, min) != 0, 2:d])
    x <- as.matrix(Y[apply(Y, 1, min) != 0, 1:(d - 1)])
    rho <- (colSums((x/y)^2)/colSums(x/y)) + colSums((y2/y)^2)/colSums((y2/y))
    init_alpha <- rep(min((1/N) * (colSums(x/y)) * (2 - rho)/(rho - 1)), (d - 
                                                                            1))
    init_alpha[rho == 2] <- 1e-06
    init_beta <- init_alpha
  } else {
    init_alpha <- init[, 1:(d - 1)]
    init_beta <- init[, d:(2 * d - 2)]
  }
  ##----------------------------------------## 
  ## Get prepared for the loop
  ##----------------------------------------## 
  alpha_hat <- rep(0, (d - 1))
  beta_hat <- rep(0, (d - 1))
  gradients <- matrix(0, 2, d - 1)
  LL <- rep(0, (d - 1))
  SE <- matrix(0, 2, (d - 1))
  niter <- rep(0, (d - 1))
  
  ##----------------------------------------## 
  ## Distribution fitting
  ##----------------------------------------## 
  options(warn = -1)
  for (i in 1:(d - 1)) {
    fit <- NULL
    data.fit <- cbind(data[, i], Y2[, i])
    init.fit <- c(init_alpha[i], init_beta[i])
    try(fit <- DMD.DM.fit(data = data.fit, weight = weight, init = init.fit, 
                          epsilon = epsilon, maxiters = maxiters, display), silent = TRUE)
    if (is.null(fit)) {
      stop("GDM model is not suitable for this dataset. \n\t\t\t\tPlease use anoter model or privide initial value.")
    }
    alpha_hat[i] <- fit$estimate[1]
    beta_hat[i] <- fit$estimate[2]
    gradients[, i] <- fit$gradient
    LL[i] <- fit$logL
    SE[, i] <- fit$SE
    niter[i] <- fit$iter
  }
  options(warn = ow)
  
  ##----------------------------------------## 
  ##Check the gradients
  ##----------------------------------------## 
  if (mean(gradients^2) > 1e-04) 
    warning(paste("The algorithm doesn't converge within", niter, "iterations.", sep = " "))
  
  ##----------------------------------------## 
  ## Compute output statistics
  ##----------------------------------------## 
  p_MN <- matrix(apply(data/m, 2, mean), N, d, byrow = TRUE)
  logL_MN <- sum(weight * data * log(p_MN)) + sum(lgamma(m + 1)) - sum(lgamma(data + 
                                                                                1))
  LRT <- 2 * (sum(LL) - logL_MN)
  p_value <- pchisq(q = LRT, df = (d - 1), ncp = 0, lower.tail = FALSE, log.p = FALSE)
  BIC <- -2 * sum(LL) + log(N) * (d - 1) * 2
  AIC <- -2 * sum(LL) + 2 * (d - 1) * 2
  fitted <- alpha_hat / (alpha_hat + beta_hat)
  for (f in 2:(d - 1)) {
    fitted[f] <- (1 - sum(fitted[1:(f - 1)])) * fitted[f]
  }
  fitted[d] <- (1 - sum(fitted[1:(d - 1)]))
  DoF <- 2 * (d - 1)
  
  ##----------------------------------------## 
  ## Clean up the results
  ##----------------------------------------## 
  estimate <- c(alpha_hat, beta_hat)
  if (!is.null(colnames(data))) {
    names(estimate) <- c(paste("alpha", colnames(data)[1:(d - 1)], sep = "_"), 
                        paste("beta", colnames(data)[1:(d - 1)], sep = "_"))
  } else {
    names(estimate) <- c(paste("alpha", 1:(d - 1), sep = "_"), 
                         paste("beta", 1:(d - 1), sep = "_"))
  }
  
 
  ##----------------------------------------## 
  list(estimate = estimate, SE = c(SE), gradient = gradients, logL = sum(LL), 
       BIC = BIC, AIC = AIC, iter = sum(niter), LRT = LRT, LRTpvalue = p_value, 
       fitted = fitted, DoF = DoF, distribution = "Generalized Dirichlet Multinomial", 
       vcov = matrix())
}





##============================================================## 
## Function: predict ; called by predict 
##============================================================##

DMD.GDM.predict <- function(beta, d, newdata) {
		alpha <- beta[, 1:(d - 1)]
    beta <- beta[, d:ncol(beta)]
    pred <- matrix(0, nrow(newdata), d)
    pred[, 1:(d - 1)] <- exp(newdata %*% alpha)/(exp(newdata %*% alpha) + exp(newdata %*% 
            beta))
    pred[, 2] <- (1 - pred[, 1]) * pred[, 2]
    if (nrow(newdata) > 1) {
            for (f in 3:(d - 1)) {
                pred[, f] <- (1 - rowSums(pred[, 1:(f - 1)])) * pred[, f]
            }
            pred[, d] <- (1 - rowSums(pred[, 1:(d - 1)]))
    } else {
            for (f in 3:(d - 1)) {
                pred[, f] <- (1 - sum(pred[, 1:(f - 1)])) * pred[, f]
            }
            pred[, d] <- (1 - sum(pred[, 1:(d - 1)]))
    }
    return(pred)
}


##============================================================## 
## Function: loss ; called by MGLMsparsereg 
##============================================================##
DMD.GDM.loss <- function(Y, X, beta, dist, weight){
  N <- nrow(Y)
  d <- ncol(Y)
  p <- ncol(X)
  m <- rowSums(Y)
  
  Ys <- t(apply(apply(apply(Y, 1, rev), 2, cumsum), 2, rev))
  alpha <- exp(X %*% beta)
  A <- alpha[, 1:(d - 1)]
  B <- alpha[, d:(2 * (d - 1))]
  loss <- -sum(weight * dgdirmn(Y, A, B))
  dalpha1 <- digamma(A + Y[, -d]) + digamma(A + B) - digamma(A) - digamma(A + B + Ys[, -d])
  dalpha2 <- digamma(B + Ys[, -1]) + digamma(A + B) - digamma(B) - digamma(A + B + Ys[, -d])
  dalpha <- cbind(dalpha1, dalpha2)
  lossD1 <- -rowSums(sapply(1:nrow(Y), 
         function(i, A, B, w) return(w[i] * A[i, ] %x% B[i, ]), alpha * dalpha, X, weight))
  lossD1 <- matrix(lossD1, p, 2 * (d - 1))
  
  return(list(loss, lossD1))
  
}
