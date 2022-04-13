##============================================================## 
## Function: fit Negative Multinomial
##============================================================##
#' @rdname MGLMfit 
DMD.NegMN.fit <- function(data, init, weight, epsilon = 1e-08, maxiters = 150, display = FALSE) {
  
  N <- nrow(data)
  d <- ncol(data)
  m <- rowSums(data)
  mm <- mean(m)
  max <- max(m)
  k <- c(0:(max - 1))
  cat_count <- colSums(data)
  data <- as.matrix(data)
  ## ----------------------------------------## 
  ## Give default value to the missing variables
  ## ----------------------------------------## 
  if (missing(weight)) 
    weight <- rep(1, N)
  Sj <- function(xj, k, weight) Sj <- colSums(weight * outer(xj, k, ">"))
  Xj <- colSums(data * weight)
  r <- Sj(m, k, weight = weight)
  mbar <- mean(m * weight)
  s2 <- sd(m)^2 * (N - 1)/N
  if (missing(init)) {
    if ((s2 - mbar) > 0) {
      init_beta <- mbar^2/(s2 - mbar)
      init_pi <- rep(0, (d + 1))
      init_pi[d + 1] <- mbar/s2
      init_pi <- c((init_pi[d + 1] * Xj)/(init_beta * N), mbar/s2)
    } else {
      warning("The data shows underdispersion, better use other models.\n\t\t\t\tbeta_init=1 is assigned here.")
      init_beta <- 1
      init_pi <- rep(1/(d + 1), d + 1)
    }
  } else {
    init_beta <- init[(d + 2)]
    init_pi <- init[1:(d + 1)]
  }
  
  ##----------------------------------------## 
  ## Get prepared for the main loop
  ##----------------------------------------## 
  beta_hat <- init_beta
  pi_hat <- init_pi[1:d]
  pi2_hat <- init_pi[(d + 1)]
  LL1 <- sum(lgamma(beta_hat + m)) + sum(data %*% log(pi_hat)) + N * beta_hat * 
    log(pi2_hat) - sum(lgamma(data + 1))
  LL_iter <- rep(NA, maxiters)
  LL_iter[1] <- LL1
  niter <- 1
  
  ##----------------------------------------## 
  ## The main loop
  ##----------------------------------------## 
  while ((niter <= 2) || (((LL2 - LL1)/(abs(LL1) + 1) > epsilon) & (niter < maxiters))) {
    
    LL1 <- LL_iter[niter]
    niter <- niter + 1
    
    ##----------------------------------------## 
    ## MM Update
    beta_MM <- -(sum(r * beta_hat/(beta_hat + k))/(N * log(pi2_hat)))
    pi2_MM <- (N * beta_MM)/(sum(Xj) + N * beta_MM)
    pi_MM <- Xj/(sum(Xj) + N * beta_MM)
    LL_MM <- sum(r * log(beta_MM + k)) + sum(Xj * log(pi_MM)) + N * beta_MM * 
      log(pi2_MM) - sum(lgamma(data + 1))
    
    ##----------------------------------------## 
    ## Newton's Update
    score <- c(Xj/pi_hat - sum(weight) * beta_hat/pi2_hat, sum(r/(beta_hat + 
                                                                    c(0:(max - 1)))) + sum(weight) * log(pi2_hat))
    
    diaginv <- c(pi_hat^2/Xj, 1/(sum(r/(beta_hat + c(0:(max - 1)))^2) - sum(weight)/beta_hat))
    rank1v <- rep(sqrt(beta_hat * sum(weight))/pi2_hat, (d + 1))
    rank1v[d + 1] <- sqrt(sum(weight)/beta_hat)
    newton_iterate <- c(pi_hat, beta_hat) + diaginv * score - sum(diaginv * rank1v * 
                                                                    score)/(1 + sum(diaginv * rank1v^2)) * (diaginv * rank1v)
    pi_Newton <- newton_iterate[1:d]
    pi2_Newton <- 1 - sum(pi_Newton)
    beta_Newton <- newton_iterate[length(newton_iterate)]
    
    ##----------------------------------------## 
    ## Choose update
    check <- all(beta_Newton > 0) && all(pi_Newton > 0) && all(pi_Newton < 1) && 
      all(pi2_Newton < 1) && all(pi2_Newton > 0) && all(!is.nan(beta_Newton)) && 
      all(!is.nan(pi_Newton))
    if (check) {
      LL_Newton <- sum(r * log(beta_Newton + k)) + sum(Xj * log(pi_Newton)) + 
        N * beta_Newton * log(pi2_Newton) - sum(lgamma(data + 1))
      beta_hat <- beta_MM * (LL_MM >= LL_Newton) + beta_Newton * (1 - (LL_MM >= 
                                                                         LL_Newton))
      pi_hat <- pi_MM * (LL_MM >= LL_Newton) + pi_Newton * (1 - (LL_MM >= LL_Newton))
      pi2_hat <- pi2_MM * (LL_MM >= LL_Newton) + pi2_Newton * (1 - (LL_MM >= 
                                                                      LL_Newton))
      LL2 <- LL_MM * (LL_MM >= LL_Newton) + LL_Newton * (1 - (LL_MM >= LL_Newton))
    } else {
      beta_hat <- beta_MM
      pi_hat <- pi_MM
      pi2_hat <- pi2_MM
      LL2 <- LL_MM
    }
    
    ##----------------------------------------## 
    ## Display the method used, if asked
    if (display) {
      if (!check) 
        print(paste("iteration", niter, ", MM Algorithm, logL=", LL_MM, sep = "")) else if (LL_MM >= LL_Newton) 
          print(paste("iteration", niter, ", MM Algorithm, logL=", LL_MM, sep = "")) else print(paste("iteration", niter, ", Newton's Method, logL=", LL_Newton, 
                                                                                                      sep = ""))
    }
    LL_iter[niter] <- LL2
    
  }
  ##----------------------------------------## 
  ## End of the main iteration
  ##----------------------------------------## 
  
  ##----------------------------------------## 
  ## Check gradients
  ##----------------------------------------## 
  score <- c(Xj/pi_hat - sum(weight) * beta_hat/pi2_hat, sum(r/(beta_hat + c(0:(max - 
                                                                                  1)))) + sum(weight) * log(pi2_hat))
  
  if (mean(score^2) > 1e-04) 
    warning(paste("The algorithm doesn't converge within ", niter, " iterations", 
                  sep = ""))
  
  ##----------------------------------------## 
  ## Compute output statistics
  ##----------------------------------------## 
  BIC <- -2 * LL2 + log(N) * (d + 1)
  AIC <- -2 * LL2 + 2 * (d + 1)
  dinv <- c(pi_hat^2/Xj, 1/(sum(r/(beta_hat + k)^2) - sum(weight)/beta_hat))
  r1v <- rep(sqrt(beta_hat * sum(weight))/pi2_hat, d + 1)
  r1v[length(r1v)] <- sqrt(sum(weight)/beta_hat)
  SE <- sqrt(dinv - (dinv * r1v)^2/(1 + sum(dinv * r1v^2)))
  fitted <- beta_hat * pi_hat/sum(c(pi_hat, pi2_hat))
  
  ##----------------------------------------## 
  ## Clean up the results
  ##----------------------------------------## 
  estimate = c(pi_hat, beta_hat)
  if (!is.null(colnames(data))) {
    names(estimate) <- c(paste("p", colnames(data), sep = "_"), "phi")
  } else {
    names(estimate) <- c(paste("p", 1:d, sep = "_"), "phi")
  }
  
  
  ##----------------------------------------##
  list(estimate = estimate, DoF = (d + 1), gradient = score, SE = SE, 
       logL = LL2, lliter = LL_iter, BIC = BIC, AIC = AIC, itern = niter, LRT = numeric(), 
       LRTpvalue = numeric(), fitted = fitted, vcov = matrix(),
       distribution = "Negative Multinomial")
} 






## ============================================================## 
## Function: fit NegMN Reg 
## ============================================================##

DMD.NegMN.reg <- function(Y, init, X, weight, epsilon, maxiters, display, parallel, 
                          cores, cl, sys) {
  ow <- getOption("warn")
  pred <- matrix(NA, nrow(Y), ncol(Y))
  emptyRow <- rowSums(Y) == 0
  Ys <- t(apply(apply(apply(Y, 1, rev), 2, cumsum), 2, rev))
  d <- ncol(Y)
  p <- ncol(X)
  m <- rowSums(Y)
  N <- nrow(Y)
  alpha <- init[, 1:d]
  beta <- init[, (d + 1)]
  Alpha <- exp(X %*% alpha)
  rowsum_Alpha <- rowSums(Alpha) + 1
  Beta <- c(exp(X %*% beta))
  lliter <- rep(NA, maxiters)
  lliter[1] <- sum(weight * dnegmn(Y, Beta, alpha = Alpha), na.rm = TRUE)
  niter <- 1
  options(warn = -1)
  while (((niter <= 2) || ((ll2 - ll1)/(abs(ll1) + 1) > epsilon)) & (niter < maxiters)) {
    niter <- niter + 1
    ll1 <- lliter[niter - 1]
    Alpha <- exp(X %*% alpha)
    rowsum_Alpha <- rowSums(Alpha) + 1
    Beta <- c(exp(X %*% beta))
    Prob <- cbind(Alpha/rowsum_Alpha, 1/rowsum_Alpha)
    tmpBeta <- digamma(Beta + m) - digamma(Beta)
    tmpBeta[is.nan(tmpBeta)] <- 0
    w_beta <- log(rowsum_Alpha)
    dlbeta <- colSums((tmpBeta - w_beta) * Beta * X)
    hbeta_w <- (trigamma(Beta + m) - trigamma(Beta) + tmpBeta - w_beta) * Beta
    hbeta <- sapply(1:N, function(i, A, B, w) return(w[i] * A[i, ] %x% B[i, ]), 
                    X, X, hbeta_w)
    hbeta <- matrix(rowSums(hbeta), p, p)
    if (all(eigen(hbeta)$value < 0)) {
      beta_MM <- beta - solve(hbeta, dlbeta)
      Beta_MM <- c(exp(X %*% beta_MM))
      lltemp <- sum(weight * dnegmn(Y, Beta_MM, alpha = Alpha), na.rm = TRUE)
      if (lltemp < ll1) {
        Y_reg <- weight * Beta * tmpBeta/w_beta
        beta_MM <- glm.fit(X, Y_reg, weights = w_beta, family = poisson(link = "log"), 
                           start = beta)$coefficients
      }
    } else {
      Y_reg <- weight * Beta * tmpBeta/w_beta
      beta_MM <- glm.fit(X, Y_reg, weights = w_beta, family = poisson(link = "log"), 
                         start = beta)$coefficients
    }
    Beta_MM <- c(exp(X %*% beta_MM))
    alpha_MM <- matrix(0, p, d)
    w_alpha <- (Beta_MM + m)/(rowSums(Alpha) + 1)
    if (!parallel) {
      for (j in 1:d) {
        alpha_MM[, j] <- glm.fit(X, weight * Y[, j]/w_alpha, weights = w_alpha, 
                                 family = poisson(link = "log"), start = alpha[, j])$coefficients
      }
    } else {
      y.list <- split(weight * Y/w_alpha, rep(1:d, each = nrow(Y)))
      start.list <- split(alpha, rep(1:d, each = nrow(alpha)))
      if (sys[1] == "Windows") {
        fit.list <- clusterMap(cl, glm.private, y.list, start.list, .scheduling = "dynamic", 
                               MoreArgs = list(X = X, family = poisson(link = "log"), weights = w_alpha))
      } else if (sys[1] != "Windows") {
        fit.list <- mcmapply(cl, glm.private, y.list, start.list, MoreArgs = list(X = X, 
                                                                                  family = poisson(link = "log"), weights = w_alpha), mc.cores = cores, 
                             mc.preschedule = FALSE)
      }
      alpha_MM <- do.call("cbind", fit.list)
    }
    Alpha_MM <- exp(X %*% alpha_MM)
    Prob_MM <- cbind(Alpha_MM/(1 + rowSums(Alpha_MM)), 1/(1 + rowSums(Alpha_MM)))
    Beta_MM <- c(exp(X %*% beta_MM))
    ll.MM <- sum(weight * dnegmn(Y, Beta_MM, alpha = Alpha_MM), na.rm = TRUE)
    deta <- matrix(0, N, (d + 1))
    deta[, 1:d] <- Y - Alpha * Beta - Prob[, 1:d] * (m - Beta * (rowsum_Alpha - 
                                                                   1))
    deta[, (d + 1)] <- Beta * (digamma(Beta + m) - digamma(Beta) - log(rowsum_Alpha))
    score <- rowSums(sapply(1:N, function(i, A, B, weight) return(weight[i] * 
                                                                    A[i, ] %x% B[i, ]), deta, X, weight))
    hessian <- matrix(0, p * (1 + d), p * (1 + d))
    upleft <- sapply(1:N, function(i, A, B) return(A[i, ] %x% B[i, ]), cbind(Prob[, 
                                                                                  1:d], -Beta/(Beta + m)), X)
    hessian <- t(t(upleft) * (weight * (Beta + m))) %*% t(upleft)
    for (j in 1:d) {
      idx <- (j - 1) * p + (1:p)
      hessian[idx, idx] <- hessian[idx, idx] - t(X) %*% (X * (weight * (Beta + 
                                                                          m) * Prob[, j]))
    }
    tmpvector2 <- as.vector(Beta * (tmpBeta + Beta * (trigamma(Beta + m) - trigamma(Beta)) - 
                                      log(rowsum_Alpha) - Beta/(Beta + m)))
    idx <- d * p + c(1:p)
    hessian[idx, idx] = hessian[idx, idx] + t(X) %*% (X * (weight * tmpvector2))
    temp.try <- NULL
    try(temp.try <- solve(hessian, score), silent = TRUE)
    if (is.null(temp.try)) {
      ll.Newton <- NA
    } else if (is.numeric(temp.try)) {
      B_Newton <- cbind(alpha, beta) - matrix(temp.try, p, (d + 1))
      B <- exp(X %*% B_Newton)
      Alpha_Newton <- B[, 1:d]
      Beta_Newton <- B[, (d + 1)]
      ll.Newton <- sum(weight * dnegmn(Y, Beta_Newton, alpha = Alpha_Newton), na.rm = TRUE)
      if (!is.na(ll.Newton) & ll.MM > ll.Newton) {
        for (step in 1:20) {
          B_N <- cbind(alpha, beta) - matrix(temp.try/(2^step), p, d + 1)
          B <- exp(X %*% B_N)
          Alpha_N <- B[, 1:d]
          Beta_N <- B[, (d + 1)]
          llnew <- sum(weight * dnegmn(Y, Beta_N, alpha = Alpha_N), na.rm = T)
          if (is.na(llnew) | is.nan(llnew) | llnew > 0) {
            next
          } else if (llnew > ll.Newton) {
            ll.Newton <- llnew
            B_Newton <- B_N
          }
          if (!is.na(llnew) & llnew > ll1) {
            break
          }
        }
      }
      if (is.nan(ll.Newton)) 
        ll.Newton <- NA else if (ll.Newton >= 0) 
          ll.Newton <- NA
    } else {
      ll.Newton <- NA
    }
    if (is.na(ll.Newton) | ll.MM > ll.Newton) {
      if (display) 
        print(paste("Iteration ", niter, " MM update", sep = ""))
      alpha <- alpha_MM
      beta <- beta_MM
      ll2 <- ll.MM
    } else {
      if (display) 
        print(paste("Iteration ", niter, " Newton's update", sep = ""))
      alpha <- B_Newton[, 1:d]
      beta <- B_Newton[, (d + 1)]
      ll2 <- ll.Newton
    }
    lliter[niter] <- ll2
  }
  options(warn = ow)
  BIC <- -2 * ll2 + log(N) * p * (d + 1)
  AIC <- -2 * ll2 + 2 * p * (d + 1)
  A <- exp(X %*% cbind(alpha, beta))
  Alpha <- A[, 1:d]
  Beta <- A[, (d + 1)]
  Prob <- cbind(Alpha, 1)
  Prob <- Prob/(rowSums(Prob))
  tmpv2 <- Beta + m
  tmpv1 <- digamma(tmpv2) - digamma(Beta)
  SE <- matrix(NA, p, d + 1)
  wald <- rep(NA, p)
  wald.p <- rep(NA, p)
  H <- matrix(NA, p * (d + 1), p * (d + 1))
  if (any(A == Inf, is.nan(A))) {
    warning("Out of range of trigamma.  No SE or tests results reported.\n\t\t\t\tRegression parameters diverge.Recommend multinomial logit model")
  } else {
    deta <- matrix(0, N, (d + 1))
    deta[, 1:d] <- Y - Alpha * Beta - Prob[, 1:d] * (m - Beta * (rowsum_Alpha - 
                                                                   1))
    deta[, (d + 1)] <- Beta * (digamma(Beta + m) - digamma(Beta) + log(Prob[, 
                                                                            d + 1]))
    score <- rowSums(sapply(1:N, function(i, A, B, weight) return(weight[i] * 
                                                                    A[i, ] %x% B[i, ]), deta, X, weight))
    H <- sapply(1:N, function(i, a, x) return(a[i, ] %x% x[i, ]), cbind(Prob[, 
                                                                             1:d], -Beta/tmpv2), X)
    H <- H %*% (weight * tmpv2 * t(H))
    for (i in 1:d) {
      id <- (i - 1) * p + c(1:p)
      H[id, id] <- H[id, id] - t(X) %*% (weight * tmpv2 * Prob[, i] * X)
    }
    id <- d * p + c(1:p)
    tmpv3 <- Beta * (tmpv1 + Beta * (trigamma(tmpv2) - trigamma(Beta)) + log(Prob[, 
                                                                                  (d + 1)]) - Beta/tmpv2)
    H[id, id] <- H[id, id] + t(X) %*% (weight * tmpv3 * X)
    if (mean(score^2) > 1e-04) {
      warning(paste("The algorithm doesn't converge within", niter, "iterations. The norm of the gradient is ", 
                    sum(score^2), " Please interpret hessian matrix and MLE with caution.", 
                    sep = " "))
    }
    eig <- eigen(H)$values
    if (any(eig > 0)) {
      warning("The estimate is a saddle point.")
    } else if (any(eig == 0)) {
      warning("The hessian matrix is almost singular.")
    } else if (all(eig < 0)) {
      Hinv <- chol2inv(chol(-H))
      SE <- matrix(sqrt(diag(Hinv)), p, (d + 1))
      wald <- rep(0, p)
      wald.p <- rep(0, p)
      for (j in 1:p) {
        id <- c(0:d) * p + j
        wald[j] <- c(alpha[j, ], beta[j]) %*% chol2inv(chol(Hinv[id, id])) %*% 
          c(alpha[j, ], beta[j])
        wald.p[j] <- pchisq(wald[j], (d + 1), lower.tail = FALSE)
      }
    }
  }
  tmpalpha <- exp(X %*% alpha)
  pred <- c(exp(X %*% beta)) * tmpalpha / rowSums(tmpalpha)
  
  ## ----------------------------------------## 
  ## Clean up the results
  ## ----------------------------------------## 
  coefficients <- cbind(alpha, beta)
  colnames(coefficients) <- c(colnames(Y), "phi")
  colnames(SE) <- c(colnames(Y), "phi")
  
  rownames(coefficients) <- colnames(X)
  rownames(SE) <- colnames(X)
  
  ## ----------------------------------------## 
  list(coefficients = coefficients, SE = SE, Hessian = H, BIC = BIC, AIC = AIC, 
       wald.value = wald, wald.p = wald.p, logL = lliter[niter], iter = (niter), 
       gradient = score, fitted = pred, distribution = "Negative Multinomial")
}

# DMD.NegMN.reg <- function(Y, init, X, weight, epsilon, maxiters, display,
# parallel, cores, cl, sys){ ##----------------------------------------## ## Keep
# some original values ow <- getOption('warn') fitted <- matrix(NA, nrow(Y),
# ncol(Y)) emptyRow <- rowSums(Y)==0 Ys <- t( apply(apply(
# apply(Y,1,rev),2,cumsum),2,rev) ) d <- ncol(Y) p <- ncol(X) m <- rowSums(Y) N
# <- nrow(Y) alpha <- init[, 1:d] beta <- init[, (d+1)] Alpha <- exp(X%*%alpha)
# rowsum_Alpha <- rowSums(Alpha)+1 Beta <- c(exp(X%*%beta)) lliter <- rep(NA,
# maxiters) lliter[1] <- sum(weight*dnegmn(Y, Beta, alpha=Alpha), na.rm=TRUE) ll2 <-
# lliter[1] niter <- 1 options(warn=-1)
# ##---------------------------------------## ## Begin the main loop while(
# ((niter <=2)|| ((ll2-ll1)/(abs(ll1)+1) > epsilon))&(niter<maxiters) ){ niter <-
# niter+1 ll1 <- lliter[niter-1] Alpha <- exp(X%*%alpha) rowsum_Alpha <-
# rowSums(Alpha)+1 Beta <- c(exp(X%*%beta)) Prob <-
# cbind(Alpha/rowsum_Alpha,1/rowsum_Alpha) tmpBeta <-
# digamma(Beta+m)-digamma(Beta) tmpBeta[is.nan(tmpBeta)] <- 0
# ##----------------------------------------## ## Newton Update deta <- matrix(0,
# N, (d+1)) deta[, 1:d] <- Y - Alpha*Beta - Prob[, 1:d]*(m-Beta*(rowsum_Alpha-1)
# ) deta[, (d+1)] <- Beta*(digamma(Beta+m)-digamma(Beta)-log(rowsum_Alpha)) score
# <- colSums( kr(deta, X, weight) ) hessian <- matrix(0, p*(1+d), p*(1+d)) upleft
# <- kr(cbind(Prob[,1:d], -Beta/(Beta+m)),X) hessian <-
# t(upleft*(weight*(Beta+m)))%*%upleft for(j in 1:d){ idx <- (j-1)*p + (1:p)
# hessian[idx, idx] <- hessian[idx, idx]- t(X)%*%(X*(weight*(Beta+m)*Prob[,j])) }
# tmpvector2 <- as.vector(Beta*(tmpBeta+Beta*(trigamma(Beta+m)-trigamma(Beta))-
# log(rowsum_Alpha)- Beta/(Beta+m) )) idx <- d*p + c(1:p) hessian[idx, idx]=
# hessian[idx, idx]+t(X)%*%(X*(weight*tmpvector2)) temp.try <- NULL try(temp.try
# <- solve(hessian,score), silent=TRUE) if(is.null(temp.try)){ ll.Newton <- NA
# }else if(is.numeric(temp.try)){ B_Newton <- cbind(alpha, beta)-matrix(temp.try,
# p, (d+1)) B <- exp(X%*%B_Newton) Alpha_Newton <- B[, 1:d] Beta_Newton <- B[,
# (d+1)] ll.Newton <- sum(weight*dnegmn(Y,Beta_Newton, alpha=Alpha_Newton), na.rm=TRUE)
# ## ----------------------------------------## ## Half stepping
# if(is.nan(ll.Newton) || ll.Newton>0){ ll.Newton <- NA }else
# if(!is.na(ll.Newton)&ll1 >= ll.Newton){ for(step in 1:40){ B_N <- cbind(alpha,
# beta) - matrix(temp.try*(0.5^step), p, d+1) B <- exp(X%*%B_N) Alpha_N <- B[,
# 1:d] Beta_N <- B[, (d+1)] llnew <- sum(weight*dnegmn(Y, Beta_N, alpha=Alpha_N),na.rm=T)
# if(is.nan(llnew) | is.na(llnew) | llnew>0){ next }else if(llnew > ll.Newton){
# ll.Newton <- llnew B_Newton <- B_N } if(is.na(llnew) |llnew>ll1){ break } } }
# if(is.nan(ll.Newton)) ll.Newton <- NA else if(ll.Newton>=0) ll.Newton <- NA
# }else{ ll.Newton <- NA } if(is.na(ll.Newton) || ll.Newton<ll1){ ##
# ----------------------------------------## ## MM Update w_beta <-
# log(rowsum_Alpha) dlbeta <- colSums((tmpBeta - w_beta)*Beta*X) hbeta_w <-
# (trigamma(Beta+m)-trigamma(Beta)+tmpBeta-w_beta)*Beta hbeta <- kr(X, X,
# hbeta_w) hbeta <- matrix(colSums(hbeta), p, p) if( all(eigen(hbeta)$value<0) ){
# beta_MM <- beta - solve(hbeta, dlbeta) Beta_MM <- c(exp(X%*%beta_MM)) lltemp <-
# sum(weight*dnegmn(Y, Beta_MM, alpha=Alpha), na.rm=TRUE) if(lltemp < ll1){ Y_reg <-
# weight*Beta*tmpBeta/w_beta beta_MM <- glm.fit(X, Y_reg, weights=w_beta,
# family=poisson(link='log'), start=beta)$coefficients } }else{ Y_reg <-
# weight*Beta*tmpBeta/w_beta beta_MM <- glm.fit(X, Y_reg, weights=w_beta,
# family=poisson(link='log'), start=beta)$coefficients } Beta_MM <-
# c(exp(X%*%beta_MM)) alpha_MM <- matrix(0,p,d) w_alpha
# <-(Beta_MM+m)/(rowSums(Alpha)+1) if(!parallel){ for(j in 1:d){ alpha_MM[,j] <-
# glm.fit(X, weight*Y[,j]/w_alpha, weights=w_alpha, family=poisson(link='log'),
# start=alpha[,j])$coefficients } }else{ y.list <- split(weight*Y/w_alpha,
# rep(1:d, each=nrow(Y))) start.list <- split(alpha, rep(1:d, each=nrow(alpha)))
# if(sys[1]=='Windows'){ fit.list <- clusterMap(cl, glm.private, y.list,
# start.list, .scheduling='dynamic', MoreArgs=list(X=X, family=poisson(link
# ='log'), weights=w_alpha)) }else if(sys[1]!='Windows'){ fit.list <-
# mcmapply(cl, glm.private, y.list, start.list, MoreArgs=list(X=X,
# family=poisson(link ='log'), weights=w_alpha), mc.cores = cores, mc.preschedule
# = FALSE) } alpha_MM <- do.call('cbind', fit.list) } Alpha_MM <-
# exp(X%*%alpha_MM) Prob_MM
# <-cbind(Alpha_MM/(1+rowSums(Alpha_MM)),1/(1+rowSums(Alpha_MM))) Beta_MM <-
# c(exp(X%*%beta_MM)) ll.MM <- sum(weight*dnegmn(Y, Beta_MM, alpha=Alpha_MM), na.rm=TRUE)
# ## ----------------------------------------## ## Choose update if(
# is.na(ll.Newton)|ll.MM<0&ll.MM>ll1 ){ if(display) print(paste('Iteration ',
# niter, ' MM update', ll.MM, sep='')) alpha <- alpha_MM beta <- beta_MM
# lliter[niter] <- ll.MM ll2 <- ll.MM } }else{ if(display) print(paste('Iteration
# ', niter, ' Newton's update', ll.Newton, sep='')) alpha <- B_Newton[,1:d] beta
# <- B_Newton[, (d+1)] lliter[niter] <- ll.Newton ll2 <- ll.Newton } } ##
# ----------------------------------------## ## End of the main loop
# options(warn=ow) ##----------------------------------------## ## Compute output
# statistics BIC <- -2*ll2 + log(N)*p*(d+1) AIC <- -2*ll2 + 2*p*(d+1) A <-
# exp(X%*%cbind(alpha, beta)) Alpha <- A[, 1:d] Beta <- A[,(d+1)] Prob <-
# cbind(Alpha, 1) Prob <- Prob/(rowSums(Prob)) tmpv2 <- Beta+m tmpv1 <-
# digamma(tmpv2) - digamma(Beta) SE <- matrix(NA, p,d+1) wald <- rep(NA, p)
# wald.p <- rep(NA, p) H <- matrix(NA, p*(d+1), p*(d+1) )
# ##----------------------------------------## ## Check diverge if(any(A==Inf,
# is.nan(A))){ warning(paste('Out of range of trigamma.  No SE or tests results
# reported.  Regression parameters diverge.Recommend multinomial logit
# model','\n', sep='')) }else{ ##----------------------------------------## ##
# Calculate dl deta <- matrix(0, N, (d+1)) deta[, 1:d] <- Y - Alpha*Beta - Prob[,
# 1:d]*(m-Beta*(rowsum_Alpha-1) ) deta[, (d+1)] <-
# Beta*(digamma(Beta+m)-digamma(Beta)+log(Prob[, d+1])) score <- colSums(kr(deta,
# X, weight)) ##----------------------------------------## ## Calculate H H <-
# kr(cbind(Prob[, 1:d],-Beta/tmpv2), X) H <- t(H)%*%(weight*tmpv2*H) for(i in
# 1:d){ id <- (i-1)*p + c(1:p) H[id, id] <- H[id, id]-t(X)%*%(weight*tmpv2*Prob[,
# i]*X) } id <- d*p+c(1:p) tmpv3 <- Beta*( tmpv1+
# Beta*(trigamma(tmpv2)-trigamma(Beta))+ log(Prob[, (d+1)])-Beta/tmpv2) H[id, id]
# <- H[id, id]+t(X)%*%(weight*tmpv3*X)
# ##----------------------------------------## ## Check gradients
# if(mean(score^2)>1e-4){ warning(paste('The algorithm doesn't converge within',
# niter, 'iterations. The norm of the gradient is ', sum(score^2), ' Please
# interpret hessian matrix and MLE with caution.', '\n', sep=' ') ) }
# ##-----------------------------------## ## Check whether H is negative
# definite, eig <- eigen(H)$values if(any(is.complex(eig)) || any(eig>0)){
# warning(paste('The estimate is a saddle point.', '\n', sep='')) }else
# if(any(eig==0)){ warning(paste('The hessian matrix is almost singular.','\n',
# sep='')) }else if(all(eig<0)){ ##--------------------------------------## ##
# Calculate SE and wald Hinv <-chol2inv(chol(-H) ) SE <-matrix( sqrt(diag(Hinv )
# ), p, (d+1)) wald <- rep(0, p) wald.p<- rep(0, p) for(j in 1:p){ id<- c(0:
# d)*p+j wald[j] <- c(alpha[j,] ,beta[j])%*%chol2inv(chol(Hinv[id,id]))%*%
# c(alpha[j,] ,beta[j]) wald.p[j] <- pchisq(wald[j], (d+1), lower.tail=FALSE) } }
# } fitted <- c(exp(X%*%beta)) * exp(X%*%alpha)/rowSums(exp(X%*%alpha))
# list(coefficients=cbind(alpha, beta), SE=SE, Hessian=H, BIC=BIC, AIC=AIC,
# wald.value=wald, wald.p=wald.p, logL=lliter[niter], iter=(niter),
# gradients=score, fitted=fitted) }


## ============================================================## 
## Function: fit NegMN Reg, but do not link the over-dispersion parameter
## ============================================================##

DMD.NegMN.Alpha.reg <- function(Y, init, X, weight, epsilon, maxiters, display, parallel, 
                                cores, cl, sys) {
  
  ## ----------------------------------------## 
  ## Keep some original values
  ## ----------------------------------------## 
  ow <- getOption("warn")
  pred <- matrix(NA, nrow(Y), ncol(Y))
  emptyRow <- rowSums(Y) == 0
  
  Ys <- t(apply(apply(apply(Y, 1, rev), 2, cumsum), 2, rev))
  d <- ncol(Y)
  p <- ncol(X)
  m <- rowSums(Y)
  N <- nrow(Y)
  
  alpha <- init
  Alpha <- exp(X %*% alpha)
  rowsum_Alpha <- rowSums(Alpha) + 1
  beta <- DMD.NegMN.fit(Y, weight = weight)
  beta <- beta$estimate[(d + 1)]
  if (beta < 0) 
    beta <- 1
  Beta <- rep(beta, N)
  lliter <- rep(NA, maxiters)
  lliter[1] <- sum(weight * dnegmn(Y, Beta, alpha=Alpha), na.rm = TRUE)
  ll2 <- lliter[1]
  betaiter <- rep(NA, maxiters)
  betaiter[1] <- beta
  niter <- 1
  div <- FALSE
  options(warn = -1)
  ## ----------------------------------------## 
  ## Begin the main loop
  ## ----------------------------------------## 
  while (((niter <= 2) || ((ll2 - ll1)/(abs(ll1) + 1) > epsilon)) & (niter < maxiters)) {
    
    niter <- niter + 1
    ll1 <- lliter[niter - 1]
    Alpha <- exp(X %*% alpha)
    rowsum_Alpha <- rowSums(Alpha) + 1
    Prob <- cbind(Alpha/rowsum_Alpha, 1/rowsum_Alpha)
    tmpBeta <- digamma(beta + m) - digamma(beta)
    tmpBeta[is.nan(tmpBeta)] <- 0
    
    dlbeta <- sum(tmpBeta - log(rowsum_Alpha))
    hbeta <- sum(trigamma(beta + m) - trigamma(beta))
    ## ----------------------------------------## 
    ## Newton Update
    deta <- Y - (beta + m) * Prob[, 1:d]
    score <- colSums(kr(deta, X, weight))
    score <- c(score, dlbeta)
    upleft <- kr(Prob[, 1:d], X)
    upright <- -colSums(upleft)
    upleft <- t(upleft) %*% (upleft * (weight * (beta + m)))
    for (j in 1:d) {
      idx <- (j - 1) * p + 1:p
      upleft[idx, idx] <- upleft[idx, idx] - t(X) %*% (X * (weight * (beta + 
                                                                        m) * Prob[, j]))
    }
    hessian <- rbind(cbind(upleft, upright), c(upright, hbeta))
    temp.try <- NULL
    try(temp.try <- solve(hessian, score), silent = TRUE)
    if (is.null(temp.try)) {
      ll.Newton <- NA
    } else if (is.numeric(temp.try)) {
      beta_Newton <- beta - temp.try[length(temp.try)]
      B_Newton <- alpha - matrix(temp.try[1:p * d], p, d)
      Alpha_Newton <- exp(X %*% B_Newton)
      ll.Newton <- sum(weight * dnegmn(Y, rep(beta_Newton, N), alpha=Alpha_Newton), 
                       na.rm = TRUE)
    }
    ## ----------------------------------------## 
    ## Half stepping
    if (is.nan(ll.Newton) || ll.Newton >= 0) {
      ll.Newton <- NA
    } else if (!is.na(ll.Newton) & ll1 >= ll.Newton) {
      for (step in 1:20) {
        temp <- temp.try/(2^step)
        B_N <- alpha - matrix(temp[1:(p * d)], p, d)
        Alpha_N <- exp(X %*% B_N)
        beta_N <- beta - temp[p * d + 1]
        llnew <- sum(weight * dnegmn(Y, rep(beta_N, N), alpha=Alpha_N), na.rm = T)
        if (is.na(llnew) | is.nan(llnew) | llnew > 0) {
          next
        } else if (llnew > ll.Newton) {
          ll.Newton <- llnew
          B_Newton <- B_N
          beta_Newton <- beta_N
        }
        if (!is.na(llnew) & llnew > ll1) {
          break
        }
        if (is.nan(ll.Newton) | ll.Newton >= 0) 
          ll.Newton <- NA
      }
      if (is.nan(ll.Newton)) 
        ll.Newton <- NA else if (ll.Newton >= 0) 
          ll.Newton <- NA
    }
    
    if (is.na(ll.Newton) || ll.Newton < ll1) {
      ## ----------------------------------------## 
      ## MM Update
      dlbeta <- sum(tmpBeta - log(rowsum_Alpha))
      hbeta <- sum(trigamma(beta + m) - trigamma(beta))
      beta_MM <- beta
      if (hbeta != 0) {
        temp_MM <- beta - dlbeta/hbeta
        ## check ll increase
        if (!is.na(temp_MM) && !is.nan(temp_MM) && temp_MM > 0) 
          beta_MM <- temp_MM
      }
      alpha_MM <- matrix(0, p, d)
      w_alpha <- (beta_MM + m)/rowsum_Alpha
      if (sum(w_alpha^2) == 0) 
        break else w_alpha[w_alpha == 0] <- 1
      if (any(is.na(w_alpha))) 
        stop("The algorithm diverged. Please try other model.")
      if (!parallel) {
        for (j in 1:d) {
          alpha_MM[, j] <- glm.fit(X, weight * Y[, j]/w_alpha, weights = w_alpha, 
                                   family = poisson(link = "log"), start = alpha[, j])$coefficients
        }
      } else {
        y.list <- split(weight * Y/w_alpha, rep(1:d, each = nrow(Y)))
        start.list <- split(alpha, rep(1:d, each = nrow(alpha)))
        if (sys[1] == "Windows") {
          fit.list <- clusterMap(cl, glm.private, y.list, start.list, .scheduling = "dynamic", 
                                 MoreArgs = list(X = X, family = poisson(link = "log"), weights = w_alpha))
        } else if (sys[1] != "Windows") {
          fit.list <- mcmapply(cl, glm.private, y.list, start.list, MoreArgs = list(X = X, 
                                                                                    family = poisson(link = "log"), weights = w_alpha), mc.cores = cores, 
                               mc.preschedule = FALSE)
        }
        alpha_MM <- do.call("cbind", fit.list)
      }
      Alpha_MM <- exp(X %*% alpha_MM)
      Prob_MM <- cbind(Alpha_MM/(1 + rowSums(Alpha_MM)), 1/(1 + rowSums(Alpha_MM)))
      ll.MM <- sum(weight * dnegmn(Y, rep(beta_MM, N), alpha=Alpha_MM), na.rm = TRUE)
      
      ## ----------------------------------------## 
      ## Choose update
      ## ----------------------------------------## 
      if (is.na(ll.Newton) | ll.MM > ll1) {
        if (beta_MM <= 0) {
          warning(paste("The estimate of overdispersion parameter is ,", 
                        beta, ". It is smaller or eqal to zero.  Please consider other models.", 
                        sep = " "))
          div <- TRUE
          break
        }
        beta <- beta_MM
        alpha <- alpha_MM
        lliter[niter] <- ll.MM
        ll2 <- ll.MM
        if (display) 
          print(paste("Iteration ", niter, " MM update, log-likelihood ", 
                      ll.MM, sep = ""))
      }
    } else {
      if (beta_Newton <= 0) {
        warning(paste("The estimate of overdispersion parameter is ,", beta, 
                      ". It is smaller or eqal to zero.  Please consider other models.", 
                      sep = " "))
        div <- TRUE
        break
      }
      alpha <- B_Newton
      beta <- beta_Newton
      lliter[niter] <- ll.Newton
      ll2 <- ll.Newton
      if (display) 
        print(paste("Iteration ", niter, " Newton's update, log-likelihood ", 
                    ll.Newton, sep = ""))
    }
    if (div) 
      break
  }
  ## ----------------------------------------## 
  ## End of the main loop
  ## ----------------------------------------## 
  options(warn = ow)
  
  ## ----------------------------------------## 
  ## Compute output statistics
  ## ----------------------------------------## 
  BIC <- -2 * ll2 + log(N) * (p * d + 1)
  AIC <- -2 * ll2 + 2 * (p * d + 1)
  A <- exp(X %*% alpha)
  Alpha <- A[, 1:d]
  Prob <- cbind(Alpha, 1)
  Prob <- Prob/(rowSums(Prob))
  tmpv1 <- digamma(beta + m) - digamma(beta)
  SE <- matrix(NA, p, d)
  SE_beta <- NA
  wald <- rep(NA, p)
  wald.p <- rep(NA, p)
  H <- matrix(NA, p * (d + 1), p * (d + 1))
  
  ## ----------------------------------------## 
  ## Check diverge
  ## ----------------------------------------## 
  if (any(A == Inf, is.nan(A), is.nan(beta), is.nan(beta))) {
    warning(paste("Out of range of trigamma.  No SE or tests results reported.\n\t\t\t\tRegression parameters diverge.Recommend multinomial logit model", 
                  "\n", sep = ""))
  } else {
    ## ----------------------------------------## 
    ## Calculate dl
    ## ----------------------------------------## 
    dlbeta <- sum(tmpv1 - log1p(rowSums(Alpha))) # use log1p instead of log(1+.) 
    deta <- Y - (beta + m) * Prob[, 1:d]
    # score <- rowSums(sapply(1:N, function(i, A, B, weight)
    # return(weight[i]*A[i,]%x%B[i,]), deta, X, weight))
    score <- colSums(kr(deta, X, weight))
    score <- c(score, dlbeta)
    ## ----------------------------------------## 
    ## Calculate H upleft <- sapply(1:N,
    ## function(i,A,B) return(A[i,]%x%B[i,]), Prob[,1:d], X)
    upleft <- kr(Prob[, 1:d], X)
    upright <- -colSums(kr(Prob[, 1:d], X))
    upleft <- t(upleft * (weight * (beta + m))) %*% upleft
    for (j in 1:d) {
      idx <- (j - 1) * p + 1:p
      upleft[idx, idx] <- upleft[idx, idx] - t(X) %*% (X * (weight * (beta + 
                                                                        m) * Prob[, j]))
    }
    H <- rbind(cbind(upleft, upright), c(upright, hbeta))
    
    ## ----------------------------------------## 
    ## Check gradients
    ## ----------------------------------------## 
    if (mean(score[1:(p * d)]^2) > 0.01) {
      warning(paste("The algorithm doesn't converge within", niter, "iterations. The norm of the gradient is ", 
                    sum(score^2), " Please interpret hessian matrix and MLE with caution.", 
                    "\n", sep = " "))
    }
    ## -----------------------------------## 
    ## Check whether H is negative definite,
    ## -----------------------------------## 
    eig <- eigen(H)$values
    if (any(is.complex(eig)) || any(eig > 0)) {
      warning(paste("The estimate is a saddle point.", "\n", sep = ""))
    } else if (any(eig == 0)) {
      warning(paste("The hessian matrix is almost singular.", "\n", sep = ""))
    } else if (all(eig < 0)) {
      ## --------------------------------------## 
      ## Calculate SE and wald
      ## --------------------------------------## 
      Hinv <- chol2inv(chol(-H))
      SE <- matrix(sqrt(diag(Hinv))[1:p * d], p, d)
      SE_beta <- Hinv[length(Hinv)]
      wald <- rep(0, p)
      wald.p <- rep(0, p)
      for (j in 1:p) {
        id <- c(0:(d - 1)) * p + j
        wald[j] <- alpha[j, ] %*% chol2inv(chol(Hinv[id, id])) %*% alpha[j, 
                                                                         ]
        wald.p[j] <- pchisq(wald[j], d, lower.tail = FALSE)
      }
    }
  }
  tmpalpha <- exp(X %*% alpha)
  pred <- beta * tmpalpha / rowSums(tmpalpha)
  
  ## ----------------------------------------## 
  ## Clean up the results
  ## ----------------------------------------## 
  coefficients <- list(alpha = alpha, phi = beta)
  SE <- list(SE.alpha = SE, SE.beta = SE_beta)
  colnames(coefficients$alpha) <- colnames(Y)
  colnames(SE$SE.alpha) <- colnames(Y)
  
  rownames(coefficients[[1]]) <- colnames(X)
  rownames(SE[[1]]) <- colnames(X)
  
  ## ----------------------------------------##
  list(coefficients = coefficients, SE = SE, Hessian = H, BIC = BIC, 
       AIC = AIC, wald.value = wald, wald.p = wald.p, logL = lliter[niter], 
       iter = (niter), gradient = score, fitted = pred, 
       distribution = "Negative Multinomial")
}



##============================================================## 
## Function: predict 
##============================================================##

DMD.NegMN.predict <- function(beta, d, newdata) {
	  if (is.null(beta$phi)) {
         alpha <- beta[, 1:d]
         beta <- beta[, (d + 1)]
      } else {
         alpha <- beta$alpha
         phi <- beta$phi
      }
      pred <- exp(newdata %*% alpha)/(rowSums(exp(newdata %*% alpha)) + 1)
	  return(pred)
}



##============================================================## 
## Function: loss ; called by MGLMsparsereg 
##============================================================##
DMD.NegMN.loss <- function(Y, X, beta, dist, weight, regBeta = FALSE, Beta) {
  
  N <- nrow(Y)
  d <- ncol(Y)
  p <- ncol(X)
  m <- rowSums(Y)
  
  if (regBeta) {
    P <- matrix(NA, N, (d + 1))
    alpha <- exp(X %*% beta)
    Beta <- alpha[, (d + 1)]
    alpha_rowsums <- rowSums(alpha[, 1:d]) + 1
    P[, (d + 1)] <- 1/alpha_rowsums
    P[, 1:d] <- alpha[, 1:d] * P[, (d + 1)]
    loss <- -sum(weight * dnegmn(Y, Beta, alpha=alpha[, 1:d]))
    deta <- matrix(0, nrow(Y), d + 1)
    deta[, 1:d] <- Y - alpha[, 1:d] * Beta - P[, 1:d] * (m - (alpha_rowsums - 
                                                                1) * Beta)
    deta[, (d + 1)] <- Beta * (digamma(Beta + m) - digamma(Beta) + 
                                 log(P[, (d + 1)]))
    lossD1 <- -rowSums(sapply(1:nrow(Y), function(i, A, B, w) return(w[i] * 
                                          A[i, ] %x% B[i, ]), deta, X, weight))
    lossD1 <- matrix(lossD1, p, (d + 1))
  } else {
    P <- matrix(NA, N, (d + 1))
    alpha <- exp(X %*% beta)
    alpha_rowsums <- rowSums(alpha) + 1
    P[, (d + 1)] <- 1/alpha_rowsums
    P[, 1:d] <- alpha[, 1:d] * P[, (d + 1)]
    loss <- -sum(weight * dnegmn(Y, rep(Beta, N), alpha=alpha))
    deta <- Y - (Beta + m) * P[, 1:d]
    lossD1 <- -rowSums(sapply(1:N, function(i, A, B, w) return(w[i] * A[i, 
                                               ] %x% B[i, ]), deta, X, weight))
    lossD1 <- matrix(lossD1, p, d)
  }
  
  return(list(loss, lossD1))
}