## ============================================================## 
## Function: Dirichlet Multinomial regression
## ============================================================##
objfun <- function(alpha, x, y, d, p) {
  alpha <- matrix(alpha, p, d)
  alpha <- exp(x %*% alpha)
  m <- rowSums(y)
  
  logl <- (lgamma(m + 1) + rowSums(lgamma(y + alpha) - lgamma(alpha)) + lgamma(rowSums(alpha))) - 
    (lgamma(rowSums(alpha) + m) + rowSums(lgamma(y + 1)))
  return(-sum(logl))
}

objfun.grad <- function(alpha, x, y, d, p) {
  alpha <- matrix(alpha, p, d)
  Beta <- exp(x %*% alpha)
  m <- rowSums(y)
  tmpvector <- digamma(rowSums(Beta) + m) - digamma(rowSums(Beta))
  tmpvector[is.nan(tmpvector)] <- 0
  tmpmatrix <- digamma(Beta + y) - digamma(Beta)
  tmpvector2 <- trigamma(rowSums(Beta)) - trigamma(m + rowSums(Beta))
  tmpmatrix2 <- trigamma(Beta) - trigamma(Beta + y)
  tmpmatrix2 <- Beta * tmpmatrix - Beta^2 * tmpmatrix2
  dalpha <- Beta * (tmpmatrix - tmpvector)
  
  expr <- paste("rbind(", paste(rep("dalpha", p), collapse = ","), ")", sep = "")
  dalpha <- eval(parse(text = expr))
  dalpha <- matrix(c(dalpha), nrow(x), ncol = p * d)
  expr2 <- paste("cbind(", paste(rep("x", d), collapse = ","), ")", sep = "")
  x <- eval(parse(text = expr2))
  dl <- colSums(dalpha * x)
  return(-dl)
}

objfun.hessian <- function(alpha, x, y, d, p) {
  alpha <- matrix(alpha, p, d)
  Beta <- exp(x %*% alpha)
  m <- rowSums(y)
  tmpvector <- digamma(rowSums(Beta) + m) - digamma(rowSums(Beta))
  tmpvector[is.nan(tmpvector)] <- 0
  tmpmatrix <- digamma(Beta + y) - digamma(Beta)
  tmpvector2 <- trigamma(rowSums(Beta)) - trigamma(m + rowSums(Beta))
  tmpmatrix2 <- trigamma(Beta) - trigamma(Beta + y)
  tmpmatrix2 <- Beta * tmpmatrix - Beta^2 * tmpmatrix2
  
  ## kr
  expr <- paste("rbind(", paste(rep("Beta", p), collapse = ","), ")", sep = "")
  Beta1 <- eval(parse(text = expr))
  Beta1 <- matrix(c(Beta1), nrow(x), ncol = p * d)
  expr2 <- paste("cbind(", paste(rep("x", d), collapse = ","), ")", sep = "")
  x1 <- eval(parse(text = expr2))
  Hessian <- Beta1 * x1
  Hessian <- t(Hessian) %*% (tmpvector2 * Hessian)
  
  for (i in 1:d) {
    idx <- (i - 1) * p + (1:p)
    Hessian[idx, idx] <- Hessian[idx, idx] - t(x) %*% ((tmpvector * Beta[, i] - 
                                                          tmpmatrix2[, i]) * x)
  }
  return(-Hessian)
}


DMD.DM.reg <- function(Y, init, X, weight, epsilon, maxiters, display, parallel, 
                       cores, cl, sys) {
  
  ## ----------------------------------------## 
  ## Keep some original values
  ## ----------------------------------------## 
  ow <- getOption("warn")
  pred <- matrix(NA, nrow(Y), ncol(Y))
  emptyRow <- rowSums(Y) == 0
  Xf <- X
  
  weight <- weight[rowSums(Y) != 0]
  X <- as.matrix(X[rowSums(Y) != 0, ])
  Y <- as.matrix(Y[rowSums(Y) != 0, colSums(Y) != 0])
  d <- ncol(Y)
  p <- ncol(X)
  m <- rowSums(Y)
  N <- nrow(Y)
  beta <- init
  Beta <- exp(X %*% init)
  lliter <- rep(0, maxiters)
  lliter[1] <- sum(ddirmn(Y, Beta) * weight, na.rm = TRUE)
  ll2 <- lliter[1]
  options(warn = -1)
  niter <- 1
  
  ## ----------------------------------------##
  ## Begin the main loop
  ## ----------------------------------------##
  while ((niter <= 2 || abs(ll2 - ll1)/(abs(ll1) + 1) > epsilon) & (niter < maxiters)) {
    niter <- niter + 1
    ll1 <- lliter[niter - 1]
    tmpvector <- digamma(rowSums(Beta) + m) - digamma(rowSums(Beta))
    tmpvector[is.nan(tmpvector)] <- 0
    tmpmatrix <- digamma(Beta + Y) - digamma(Beta)
    ## ----------------------------------------## 
    ## Newton update tmpvector2 <-
    ## trigamma(rowSums(Beta)) - trigamma(m+rowSums(Beta)) tmpmatrix2 <-
    ## trigamma(Beta) - trigamma(Beta+Y) tmpmatrix2 <- Beta*tmpmatrix -
    ## Beta^2*tmpmatrix2 Hessian <- kr(Beta, X) Hessian <- t(Hessian)%*%(
    ## (tmpvector2*weight)*Hessian ) for(i in 1:d){ idx <- (i-1)*p+(1:p) Hessian[idx,
    ## idx] <- Hessian[idx, idx] - t(X) %*% ( weight*((tmpvector*Beta[,i])-
    ## tmpmatrix2[,i])*X) } dalpha <- Beta*(tmpmatrix - tmpvector) dl <- colSums(
    ## kr(dalpha, X, weight) )
    Hessian <- -objfun.hessian(c(beta), X, Y, d, p)
    dl <- -objfun.grad(c(beta), X, Y, d, p)
    if (all(!is.na(dl)) & mean(dl^2) < 1e-04) 
      break
    temp.try <- NULL
    try(temp.try <- solve(Hessian, dl), silent = TRUE)
    if (is.null(temp.try) | any(is.nan(temp.try))) {
      ll.Newton <- NA
    } else if (is.numeric(temp.try)) {
      beta_Newton <- beta - matrix(temp.try, p, d)
      ll.Newton <- sum(ddirmn(Y, exp(X %*% beta_Newton)) * weight, na.rm = TRUE)
      ## ----------------------------------------## Half stepping
      if (is.nan(ll.Newton) || ll.Newton >= 0) {
        ll.Newton <- NA
      } else if (!is.na(ll.Newton) & ll1 >= ll.Newton) {
        for (st in 1:20) {
          beta_N <- beta - matrix(temp.try * (0.5^st), p, d)
          llnew <- sum(ddirmn(Y, exp(X %*% beta_N)) * weight)
          if (is.na(llnew) | is.nan(llnew) | llnew > 0) {
            next
          } else if (llnew > ll.Newton) {
            ll.Newton <- llnew
            beta_Newton <- beta_N
          }
          if (!is.na(llnew) & llnew > ll1) {
            break
          }
        }
      }
    } else {
      ll.Newton <- NA
    }
    if (is.na(ll.Newton) || ll.Newton < ll1) {
      ## ----------------------------------------## 
      ## MM update
      ## ----------------------------------------## 
      beta_MM <- beta
      weight.fit <- weight * tmpvector
      wnz <- weight.fit != 0 & weight.fit != Inf
      weight.fit <- weight.fit[wnz]
      X1 <- X[wnz, ]
      Y_new <- Beta * tmpmatrix
      Y_new <- Y_new[wnz, ]/weight.fit
      if (!parallel) {
        for (j in 1:d) {
          ## Surrogate 1 Poisson Regression
          wy <- Y_new[, j]
          wy[is.na(wy)] <- Y_new[is.na(wy), j]
          ff <- glm.fit(X1, Y_new[, j], weights = weight.fit, family = poisson(link = "log"), 
                        control = list(epsilon = 1e-08))
          # a <- b <- NA ff <- nlminb(rep(0.1, p), ll.obj, gradient=ll.grad,
          # #hessian=ll.hessian, y=wy, X=X1, w=weight.fit) b <- -ff$objective a <-
          # -ll.obj(beta[,j], wy, X1, weight.fit) b <- -ll.obj(ff$coefficients, wy, X1,
          # weight.fit)
          if (ff$converged) 
            beta_MM[, j] <- ff$coefficients  #par
          # print(paste( ff$converged, sep=' *** '))
        }
      } else {
        weight.fit[weight.fit == 0] <- 1
        y.list <- split(Y_new/weight.fit, rep(1:d, each = nrow(Y)))
        if (sys[1] == "Windows") {
          fit.list <- clusterMap(cl, glm.private, y.list, .scheduling = "dynamic", 
                                 MoreArgs = list(weights = weight.fit, X = X, family = poisson(link = "log")))
        } else if (sys[1] != "Windows") {
          fit.list <- mcmapply(cl, glm.private, y.list, MoreArgs = list(weights = weight.fit, 
                                                                        X = X, family = poisson(link = "log")), mc.cores = cores, mc.preschedule = FALSE)
        }
        beta_MM <- do.call("cbind", fit.list)
      }
      ll.MM <- sum(ddirmn(Y, exp(X %*% beta_MM)) * weight, na.rm = TRUE)
      # print(paste(ll1, ll.MM, ll1-ll.MM, sep=' '))
      ## ----------------------------------------## 
      ## Choose the update
      ## ----------------------------------------## 
      if (is.na(ll.Newton) | (ll.MM < 0 & ll.MM > ll1)) {
        if (display) 
          print(paste("Iteration ", niter, " MM update, log-likelihood ", 
                      ll.MM, sep = ""))
        beta <- beta_MM
        Beta <- exp(X %*% beta_MM)
        lliter[niter] <- ll.MM
        ll2 <- ll.MM
      }
    } else {
      if (display) 
        print(paste("Iteration ", niter, " Newton's update, log-likelihood", 
                    ll.Newton, sep = ""))
      beta <- beta_Newton
      Beta <- exp(X %*% beta_Newton)
      lliter[niter] <- ll.Newton
      ll2 <- ll.Newton
    }
  }
  ## ----------------------------------------## 
  ## End of the main loop
  ## ----------------------------------------## 
  options(warn = ow)
  
  ## ----------------------------------------## 
  ## Compute some output statistics
  ## ----------------------------------------## 
  BIC <- -2 * ll2 + log(N) * p * d
  AIC <- -2 * ll2 + 2 * p * d
  tmpvector <- digamma(rowSums(Beta) + m) - digamma(rowSums(Beta))
  tmpmatrix <- digamma(Beta + Y) - digamma(Beta)
  tmpvector2 <- trigamma(rowSums(Beta)) - trigamma(m + rowSums(Beta))
  tmpmatrix2 <- trigamma(Beta) - trigamma(Beta + Y)
  tmpmatrix2 <- Beta * tmpmatrix - Beta^2 * tmpmatrix2
  SE <- matrix(NA, p, d)
  wald <- rep(NA, p)
  wald.p <- rep(NA, p)
  H <- matrix(NA, p * d, p * d)
  
  ## ---------------------------------------------------------------## 
  ## Check diverge
  ## ---------------------------------------------------------------## 
  if (any(Beta == Inf, is.nan(Beta), is.nan(tmpvector2), is.nan(tmpvector), is.nan(tmpmatrix2))) {
    warning(paste("Out of range of trigamma(). \n No standard error or tests results reported.\n
                  Regression parameters diverge.  Recommend multinomial logit model.", 
                  "\n", sep = ""))
  } else {
    ## ----------------------------------------## 
    ## Check gradients
    ## ----------------------------------------## 
    dalpha <- Beta * (tmpmatrix - tmpvector)
    dl <- colSums(kr(dalpha, X, weight))
    if (mean(dl^2) > 1e-04) {
      warning(paste("The algorithm doesn't converge within", niter, "iterations. The norm of the gradient is", 
                    sqrt(sum(dl^2)), ". Please interpret hessian matrix and MLE with caution.", 
                    "\n", sep = " "))
    }
    ## -----------------------------------## 
    ## Check whether H is negative definite, H
    ## -----------------------------------## 
    ## <- sapply(1:N, function(i, a, b) return(a[i, ]%x%b[i, ]),Beta, X)
    H <- kr(Beta, X)
    H <- t(H) %*% (H * (tmpvector2))
    for (i in 1:d) {
      id <- (i - 1) * p + c(1:p)
      H[id, id] <- H[id, id] + t(X) %*% ((-tmpvector * Beta[, i] + tmpmatrix2[, 
                                                                              i]) * X)
    }
    eig <- eigen(H)$values
    if (any(is.complex(eig)) || any(eig > 0) || any(abs(eig) < 1e-04)) {
      warning("The estimate is a saddle point. The hessian matrix is not negative definite.")
    } else if (all(eig < 0)) {
      ## --------------------------------------## 
      ## Calculate SE and wald
      ## --------------------------------------## 
      Hinv <- chol2inv(chol(-H))
      SE <- matrix(sqrt(diag(Hinv)), p, d)
      wald <- rep(NA, p)
      wald.p <- rep(NA, p)
      for (j in 1:p) {
        id <- c(0:(d - 1)) * p + j
        invH <- NULL
        try(invH <- solve(Hinv[id, id], beta[j, ]), silent = TRUE)
        if (is.numeric(invH)) {
          wald[j] <- sum(beta[j, ] * invH)
          wald.p[j] <- pchisq(wald[j], d, lower.tail = FALSE)
        }
      }
    }
  }
  
  tmppred <- exp(Xf %*% beta)
  pred <- tmppred /rowSums(tmppred)
  pred[emptyRow, ] <- 0
  
  ## ----------------------------------------## 
  ## Clean up the results
  ## ----------------------------------------## 
  colnames(beta) <- colnames(Y)
  colnames(SE) <- colnames(Y)
  
  ## ----------------------------------------## 
  list(coefficients = beta, SE = SE, Hessian = H, wald.value = wald, wald.p = wald.p, 
       DoF = p * d, logL = ll2, BIC = BIC, AIC = AIC, iter = niter, gradient = dl, 
       fitted = pred, logLiter = lliter, distribution = "Dirichlet Multinomial")
}


##============================================================## 
## Function: fit Direchelet Multinomial
##============================================================##
#' @rdname MGLMfit 
DMD.DM.fit <- function(data, init, weight, epsilon = 1e-08, maxiters = 150, display = FALSE) {
  
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
  
  ##----------------------------------------## 
  ## Give default value to the missing variables
  ##----------------------------------------## 
  if (missing(weight)) {
    weight <- rep(1, N)
  }
  if (missing(init)) {
    
    rho <- sum((colSums(weight * (data/m)^2))/(colSums(weight * data/m)))
    if (rho == d) {
      init <- rep(1e-06, d)
    } else {
      init <- as.vector((1/N) * (colSums(weight * data/m)) * (d - rho)/(rho - 
                                                                          1))
    }
  }
  
  ##----------------------------------------## 
  ## Get prepared for the loop
  ##----------------------------------------## 
  Sj <- function(xj, k, weight) Sj <- colSums(weight * outer(xj, k, ">"))
  s <- apply(data, 2, Sj, k = k, weight = weight)
  r <- Sj(m, k, weight = weight)
  alpha_hat <- init
  niter <- 1
  log_term <- sum(lgamma(m + 1)) - sum(lgamma(data + 1))
  LL <- -sum(r * log(sum(alpha_hat) + k)) + sum(s * log(outer(k, alpha_hat, "+"))) + 
    log_term
  alpha_hat <- init
  DM.LL_iter <- rep(0, maxiters)
  DM.LL_iter[1] <- LL
  
  ##----------------------------------------## 
  ## The main loop
  ##----------------------------------------## 
  while (((niter <= 2) || ((DM.LL2 - DM.LL1)/(abs(DM.LL1) + 1) > epsilon)) & (niter < 
                                                                              maxiters)) {
    
    niter <- niter + 1
    DM.LL1 <- -sum(r * log(sum(alpha_hat) + k)) + sum(s * log(outer(k, alpha_hat, 
                                                                    "+"))) + log_term
    
    ##----------------------------------------## 
    ## MM update
    numerator <- colSums(s/(outer(rep(1, length(k)), alpha_hat) + outer(k, rep(1, 
                                                                               d))))
    denominator <- sum(r/(sum(alpha_hat) + k))
    alpha_MM <- alpha_hat * numerator/denominator
    
    ##----------------------------------------## 
    ## Newton update
    a <- sum(r/(sum(alpha_hat) + k)^2)
    b <- colSums(s/(outer(rep(1, length(k)), alpha_hat) + outer(k, rep(1, d)))^2)
    dl <- (colSums(s/(outer(rep(1, max), alpha_hat) + outer(k, rep(1, d)))) - 
             sum(r/(sum(alpha_hat) + k)))
    alpha_Newton <- alpha_hat + (1/b) * dl + (a * sum((1/b) * dl)/(1 - a * sum(1/b))) * 
      (1/b)
    
    ##----------------------------------------## 
    ## Choose the update
    if (any(is.na(alpha_Newton)) | any(alpha_Newton < 0)) {
      if (is.nan(alpha_MM) || is.na(alpha_MM) || any(alpha_MM == Inf)) {
        stop("DM model is not suitable for this dataset. \n\t\t\t\tPlease use anoter model or privide initial value.")
      } else {
        alpha_hat <- alpha_MM
        if (display) 
          print(paste("Iteration", niter, "MM update", sep = " "))
      }
    } else {
      LL.MM <- -sum(r * log(sum(alpha_MM) + k)) + sum(s * log(outer(k, alpha_MM, 
                                                                    "+"))) + log_term
      LL.Newton <- -sum(r * log(sum(alpha_Newton) + k)) + sum(s * log(outer(k, 
                                                                            alpha_Newton, "+"))) + log_term
      if (LL.MM > LL.Newton) {
        alpha_hat <- alpha_MM
        if (display) 
          print(paste("Iteration", niter, "MM update", sep = " "))
      } else {
        alpha_hat <- alpha_Newton
        if (display) 
          print(paste("Iteration", niter, "Newton update", sep = " "))
      }
    }
    DM.LL2 <- -sum(r * log(sum(alpha_hat) + k)) + sum(s * log(outer(k, alpha_hat, 
                                                                    "+"))) + log_term
    DM.LL_iter[niter] <- DM.LL2
  }
  ##----------------------------------------## 
  ## End of the main loop
  ##----------------------------------------## 
  
  ##----------------------------------------## 
  ## Check the gradients
  ##----------------------------------------## 
  a <- sum(r/(sum(alpha_hat) + k)^2)
  b <- colSums(s/(outer(rep(1, max), alpha_hat) + outer(k, rep(1, d)))^2)
  dl <- (colSums(s/(outer(rep(1, max), alpha_hat) + outer(k, rep(1, d)))) - sum(r/(sum(alpha_hat) + 
                                                                                     k)))
  if (mean(dl^2) > 1e-04) 
    warning(paste("The algorithm doesn't converge within", niter, "iterations.", sep = " "))
  
  ##----------------------------------------## 
  ## Compute output statistics 1)SE 
  ## 2) LRT test against the MN model 3) BIC
  ##----------------------------------------## 
  invI <- diag(1/b) + a/(1 - a * sum(1/b)) * outer(1/b, 1/b, "*")
  SE <- sqrt(diag(invI))
  p_MN <- t(matrix(apply(data/m, 2, mean), d, N))
  logL_MN <- sum(weight * data * log(p_MN)) + sum(lgamma(m + 1)) - sum(lgamma(data + 
                                                                                1))
  LRT <- 2 * (DM.LL2 - logL_MN)
  p_value <- pchisq(q = LRT, df = 1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
  BIC <- -2 * DM.LL2 + log(N) * d
  AIC <- -2 * DM.LL2 + 2 * d
  DoF <- d
  fitted <- alpha_hat/sum(alpha_hat)
  
  ##----------------------------------------## 
  ## Clean up the results
  ##----------------------------------------## 
  estimate <- alpha_hat
  if (!is.null(colnames(data))) {
    names(estimate) <- paste("alpha", colnames(data), sep = "_")
  } else {
    names(estimate) <- paste("alpha", 1:d, sep = "_")
  }
  

  ##----------------------------------------## 
  list(estimate = estimate, SE = SE, vcov = invI, gradient = dl, logL = DM.LL2, 
       iter = niter, BIC = BIC, AIC = AIC, LRT = LRT, LRTpvalue = p_value, fitted = fitted, 
       DoF = DoF, distribution = "Dirichlet Multinomial")
}



##============================================================## 
## Function: predict 
##============================================================##

DMD.DM.predict <- function(beta, d, newdata) {
	 pred <- exp(newdata %*% beta)/rowSums(exp(newdata %*% beta))
	 return(pred)
}


##============================================================## 
## Function: loss ; called by MGLMsparsereg 
##============================================================##
DMD.DM.loss <- function(Y, X, beta, dist, weight){
  
  N <- nrow(Y)
  d <- ncol(Y)
  p <- ncol(X)
  m <- rowSums(Y)
  
  alpha <- exp(X %*% beta)
  loss <- -sum(weight * ddirmn(Y, alpha))
  tmpvector <- digamma(rowSums(alpha) + m) - digamma(rowSums(alpha))
  tmpmatrix <- digamma(alpha + Y) - digamma(alpha)
  dalpha <- tmpmatrix - tmpvector
  lossD1 <- -rowSums(sapply(1:nrow(Y), function(i, A, B, w) return(w[i] * A[i, 
                                      ] %x% B[i, ]), alpha * dalpha, X, weight))
  lossD1 <- matrix(lossD1, p, d)
  
  return(list(loss, lossD1))
}