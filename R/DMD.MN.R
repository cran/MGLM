
##============================================================## 
## Function: fit Multinomial 
##============================================================##
DMD.MN.fit <- function(data, weight) {
  
  N <- nrow(data)
  d <- ncol(data)
  
  t <- sum(data)
  m <- rowSums(data)
  if (missing(weight)) weight <- rep(1, N)
  estimate <-apply(data/m, 2, mean) 
  SE <- sqrt(estimate * (1-estimate)/N)
  p_MN <- t(matrix(estimate, d, N))
  logL <- sum(weight * data * log(p_MN)) + sum(lgamma(m + 1)) - sum(lgamma(data + 1))
  BIC <- -2 * logL + log(N) * (d - 1)
  AIC <- -2 * logL + 2 * (d-1)
  
  
  ##----------------------------------------## 
  ## Clean up the results
  ##----------------------------------------## 
  estimate = apply(data/m, 2, mean)
  if (!is.null(colnames(data))) {
    names(estimate) <- paste("alpha", colnames(data), sep = "_")
  } else {
    names(estimate) <- paste("alpha", 1:d, sep = "_")
  }
  
 
  ##----------------------------------------## 
  list(estimate = estimate, SE = SE, logL = logL, BIC = BIC,
       AIC = AIC, distribution = "Multinomial", LRT = numeric(), 
       LRTpvalue = numeric(), iter = numeric(),
       vcov = matrix(), gradient = matrix(), fitted = vector())
}



## ============================================================## 
## Function: fit Multinomial Reg The dimension of Y is n by d (d>=3)
## ============================================================##
DMD.MN.reg <- function(Y, init, X, weight, epsilon, maxiters, display, parallel, 
                       cores, cl, sys) {
  
  ## ----------------------------------------## 
  ## Function to calculate log-likelihood
  ## ----------------------------------------## 
  dmultn_L <- function(B1) {
    B <- cbind(B1, 0)
    alpha <- exp(X %*% B)
    return(sum(Y * (X %*% B) * weight) - sum(Y * log(rowSums(alpha)) * weight) + 
             sum(lgamma(m + 1) * weight) - sum(lgamma(Y + 1) * weight))
  }
  
  ##----------------------------------------## 
  ## Function to calculate gradients
  ##----------------------------------------## 
  dl.MN_fctn <- function(i, beta) {
    x <- (weight * X)[i, ]
    y <- Y[i, ]
    prob <- exp(colSums(x * beta))/(sum(exp(colSums(x * beta))) + 1)
    dl <- y[-ncol(Y)] - sum(y) * prob
    return(matrix(c(outer(x, dl)), (length(y) - 1) * length(x), 1))
  }
  
  ## ----------------------------------------## 
  ## Function to calculate hessian matrix
  ## ----------------------------------------## 
  H.MN_fctn <- function(i, beta) {
    x <- (weight * X)[i, ]
    y <- Y[i, ]
    prob <- exp(colSums(x * beta))/(sum(exp(colSums(x * beta))) + 1)
    if (d > 2) 
      dprob <- diag(prob) else if (d == 2) 
        dprob <- prob
    d2l <- sum(y) * (outer(prob, prob) - dprob)
    o.x <- outer(x, x)
    H <- d2l %x% o.x
    return(H)
  }
  
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
  m <- rowSums(Y)
  p <- ncol(X)
  N <- nrow(Y)
  beta <- init
  lliter <- rep(0, maxiters)
  lliter[1] <- dmultn_L(init)
  niter <- 1
  options(warn = -1)
  
  ## ----------------------------------------## 
  ## Begin the main loop
  ## ----------------------------------------## 
  while (((niter <= 2) || ((ll2 - ll1)/(abs(ll1) + 1) > epsilon)) & (niter < maxiters)) {
    
    niter <- niter + 1
    ll1 <- lliter[niter - 1]
    
    ## ----------------------------------------## 
    ## Newton update
    ## ----------------------------------------## 
    P <- exp(X %*% beta)
    P <- P/(rowSums(P) + 1)
    dl <- Y[, -d] - P * m
    if (d > 2) {
      dl <- colSums(kr(dl, X, weight))
    } else if (d == 2) {
      dl <- colSums(dl * weight * X)
    }
    # dl <- Reduce('+', lapply(1:N, dl.MN_fctn, beta) )
    H <- Reduce("+", lapply(1:N, H.MN_fctn, beta))
    update <- NULL
    try(update <- solve(H, dl), silent = TRUE)
    if (is.null(update)) {
      ll.Newton <- NA
    } else if (is.numeric(update)) {
      beta_Newton <- beta - matrix(update, p, d - 1)
      ll.Newton <- dmultn_L(beta_Newton)
      ## ----------------------------------------## 
      ## Half stepping
      ## ----------------------------------------## 
      if (is.nan(ll.Newton) || ll.Newton >= 0) {
        ll.Newton <- NA
      } else if (!is.na(ll.Newton) & ll1 >= ll.Newton) {
        for (step in 1:20) {
          beta_N <- beta - matrix(update/(2^step), p, d - 1)
          llnew <- dmultn_L(beta_N)
          if (is.nan(llnew) || llnew > 0) {
            next
          } else if (llnew > ll.Newton) {
            ll.Newton <- llnew
            beta_Newton <- beta_N
          }
          if (is.na(llnew) | llnew > ll1) {
            break
          }
        }
      }
      if (is.nan(ll.Newton) | ll.Newton >= 0) 
        ll.Newton <- NA
    } else {
      ll.Newton <- NA
    }
    
    if (is.na(ll.Newton) || ll.Newton < ll1) {
      ## ----------------------------------------## MM update
      denominator <- 1 + rowSums(exp(X %*% beta))
      a <- m/denominator
      beta_MM <- matrix(0, p, (d - 1))
      if (!parallel) {
        for (i in 1:(d - 1)) {
          beta_MM[, i] <- glm.fit(X, Y[, i]/a, weights = weight * a, family = poisson(link = "log"))$coefficients
        }
      } else {
        y.list <- split(Y/a, rep(1:d, each = nrow(Y)))
        y.list[[d]] <- NULL
        if (sys == "Windows") {
          fit.list <- clusterMap(cl, glm.private, y.list, .scheduling = "dynamic", 
                                 MoreArgs = list(weights = weight * a, X = X, family = poisson(link = "log")))
        } else if (sys != "Windows") {
          fit.list <- mcmapply(cl, glm.private, y.list, MoreArgs = list(weights = weight * 
                                                                          a, X = X, family = poisson(link = "log")), mc.cores = cores, 
                               mc.preschedule = FALSE)
        }
        beta_MM <- do.call("cbind", fit.list)
      }
      ll.MM <- dmultn_L(beta_MM)
      ## ----------------------------------------## Choose the update
      if (is.na(ll.Newton) | (ll.MM < 0 & ll.MM > ll1)) {
        if (display) 
          print(paste("Iteration ", niter, " MM update", ll.MM, sep = ""))
        beta <- beta_MM
        lliter[niter] <- ll.MM
        ll2 <- ll.MM
      }
    } else {
      if (display) 
        print(paste("Iteration ", niter, " Newton's update", ll.Newton, sep = ""))
      beta <- beta_Newton
      lliter[niter] <- ll.Newton
      ll2 <- ll.Newton
    }
  }
  ## ----------------------------------------## 
  ## End of the main loop
  ## ----------------------------------------## 
  options(warn = ow)
  
  ## ----------------------------------------## 
  ## Compute output statistics
  ## ----------------------------------------## 
  BIC <- -2 * lliter[niter] + log(N) * p * (d - 1)
  AIC <- 2 * p * (d - 1) - 2 * lliter[niter]
  P <- exp(X %*% beta)
  P <- cbind(P, 1)
  P <- P/rowSums(P)
  if (d > 2) {
    H <- kr(P[, 1:(d - 1)], X)
  } else if (d == 2) {
    H <- P[, 1] * X
  }
  H <- t(H) %*% (weight * m * H)
  for (i in 1:(d - 1)) {
    id <- (i - 1) * p + (1:p)
    H[id, id] = H[id, id] - t(X) %*% (weight * m * P[, i] * X)
  }
  
  SE <- matrix(NA, p, (d - 1))
  wald <- rep(NA, p)
  wald.p <- rep(NA, p)
  
  
  ## ----------------------------------------## 
  ## Check the estimate
  ## ----------------------------------------## 
  dl <- Reduce("+", lapply(1:N, dl.MN_fctn, beta))
  if (mean(dl^2) > 0.001) {
    warning(paste("The algorithm doesn't converge within", niter, "iterations. Please check gradient.", 
                  "\n", sep = " "))
  }
  eig <- eigen(H)$values
  if (any(is.complex(eig)) || any(eig > 0)) {
    warning("The estimate is a saddle point.")
  } else if (any(eig == 0)) {
    warning(paste("The Hessian matrix is almost singular.", "\n", sep = ""))
  } else {
    ## ----------------------------------------## 
    ## If Hessian is negative definite, estimate SE and Wald
    ## ----------------------------------------## 
    Hinv <- chol2inv(chol(-H))
    SE <- matrix(sqrt(diag(Hinv)), p, (d - 1))
    wald <- rep(0, p)
    wald.p <- rep(0, p)
    for (j in 1:p) {
      id <- c(0:(d - 2)) * p + j
      wald[j] <- beta[j, ] %*% solve(Hinv[id, id]) %*% beta[j, ]
      wald.p[j] <- pchisq(wald[j], (d - 1), lower.tail = FALSE)
    }
  }
  
  
  ## ----------------------------------------## 
  ## Fitted value
  ## ----------------------------------------## 
  tmppred <- exp(Xf %*% beta)
  pred <- tmppred / (rowSums(tmppred) + 1)
  pred <- cbind(pred, (1 - rowSums(pred)))
  pred[emptyRow, ] <- 0
  
  ## ----------------------------------------## 
  ## Clean up the results
  ## ----------------------------------------##
  colnames(beta) <- colnames(Y)[1:(ncol(Y) - 1)]
  colnames(SE) <- colnames(Y)[1:(ncol(Y) - 1)]
  
  ## ----------------------------------------##
  list(coefficients = beta, SE = SE, Hessian = H, wald.value = wald, wald.p = wald.p, 
       DoF = p * (d - 1), logL = lliter[niter], BIC = BIC, AIC = AIC, iter = niter, 
       gradient = dl, fitted = pred, cl = cl, distribution = "Multinomial")
  
}



##============================================================## 
## Function: predict 
##============================================================##

DMD.MN.predict <- function(beta, d, newdata) {
	pred <- exp(newdata %*% beta)/(rowSums(exp(newdata %*% beta)) + 1)
  pred <- cbind(pred, (1 - rowSums(pred)))
	return(pred)
}


##============================================================## 
## Function: loss ; called by MGLMsparsereg 
##============================================================##
DMD.MN.loss <- function(Y, X, beta, dist, weight){
  
  N <- nrow(Y)
  d <- ncol(Y)
  p <- ncol(X)
  m <- rowSums(Y)
  
  P <- matrix(NA, N, d)
  P[, d] <- rep(1, N)
  P[, 1:(d - 1)] <- exp(X %*% beta)
  P <- P/rowSums(P)
  loss <- -sum(weight * dmn(Y, P))
  kr1 <- Y[, -d] - P[, 1:(d - 1)] * m
  if (d > 2) {
    lossD1 <- -colSums(kr(kr1, X, weight))
  } else if (d == 2) {
    lossD1 <- -colSums(kr1 * weight * X)
  }
  lossD1 <- matrix(lossD1, p, (d - 1))
  
  return(list(loss, lossD1))
}