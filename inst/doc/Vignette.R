### R code from vignette source 'Vignette.Rnw'

###################################################
### code chunk number 1: Vignette.Rnw:38-39
###################################################
require(MGLM)


###################################################
### code chunk number 2: Vignette.Rnw:52-58
###################################################
set.seed(123)
n <- 200
d <- 4
alpha <- rep(1,d)
m <- 50
Y <- rmn(m, alpha, n)


###################################################
### code chunk number 3: Vignette.Rnw:61-63
###################################################
mnFit <- MGLMfit(Y, dist="DM")
print(mnFit)


###################################################
### code chunk number 4: Vignette.Rnw:70-76
###################################################
set.seed(123)
n <- 200
d <- 4
alpha <- rep(1, d)
m <- 50
Y <- rdirm(m, alpha, n)


###################################################
### code chunk number 5: Vignette.Rnw:78-80
###################################################
dmFit <- MGLMfit(Y, dist="DM")
print(dmFit)


###################################################
### code chunk number 6: Vignette.Rnw:87-89
###################################################
gdmFit <- MGLMfit(Y, dist="GDM")
print(gdmFit)


###################################################
### code chunk number 7: Vignette.Rnw:97-106
###################################################
set.seed(124)
n <- 200
d <- 4
alpha <- rep(1, d-1)
beta <- rep(1, d-1)
m <- 50
Y <- rgdirm(m, alpha, beta, n)
gdmFit <- MGLMfit(Y, dist="GDM")
print(gdmFit)


###################################################
### code chunk number 8: Vignette.Rnw:112-121
###################################################
set.seed(1220)
n <- 100
d <- 4
p <- 5
prob <- rep(0.2, d)
beta <- 10
Y <- rnegmn(prob, beta, n)
negmnFit <- MGLMfit(Y, dist="NegMN")
print(negmnFit)


###################################################
### code chunk number 9: Vignette.Rnw:139-152
###################################################
set.seed(1234)
n <- 200
p <- 5
d <- 4
X <- matrix(runif(p*n), n, p)
alpha <- matrix(c(0.6, 0.8, 1), p, d-1, byrow=TRUE)
alpha[c(1,2),] <- 0
Alpha <- exp(X%*%alpha)
beta <- matrix(c(1.2, 1, 0.6), p, d-1, byrow=TRUE)
beta[c(1,2),] <- 0
Beta <- exp(X%*%beta)
m <- runif(n, min=0, max=25) + 25
Y <- rgdirm(m, Alpha, Beta)


###################################################
### code chunk number 10: Vignette.Rnw:158-160
###################################################
mnReg <- MGLMreg(Y~0+X, dist="MN")
print(mnReg)


###################################################
### code chunk number 11: Vignette.Rnw:167-169
###################################################
dmReg <- MGLMreg(Y~0+X, dist="DM")
print(dmReg)


###################################################
### code chunk number 12: Vignette.Rnw:176-178
###################################################
gdmReg <- MGLMreg(Y~0+X, dist="GDM")
print(gdmReg)


###################################################
### code chunk number 13: Vignette.Rnw:185-187
###################################################
negReg <- MGLMreg(Y~0+X, dist="NegMN", regBeta=FALSE)
print(negReg)


###################################################
### code chunk number 14: Vignette.Rnw:196-199
###################################################
newX <- matrix(runif(1*p), 1, p)
pred <- predict(gdmReg, newX)
pred


###################################################
### code chunk number 15: Vignette.Rnw:215-225
###################################################
n <- 100
p <- 10
d <- 5
m <- runif(n, min=0, max=25) + 25
set.seed(134)
X <- matrix(rnorm(n*p),n, p)
alpha <- matrix(0, p, d)
alpha[c(1,3, 5), ] <- 1
Alpha <- exp(X%*%alpha)
Y <- rdirm(size=m, alpha=Alpha)


###################################################
### code chunk number 16: Vignette.Rnw:232-234
###################################################
sweep <- MGLMtune(Y~0+X, dist="DM", penalty="sweep", ngridpt=30)
print(sweep$select)


###################################################
### code chunk number 17: sweeppath
###################################################
p1 <- plot(x=log(sweep$path$Lambda),y=sweep$path$BIC, type="p",
           xlab="log(lambda)", ylab="BIC") +
  lines(x=log(sweep$path$Lambda),y=sweep$path$BIC)


###################################################
### code chunk number 18: Vignette.Rnw:248-249
###################################################
p1 <- plot(x=log(sweep$path$Lambda),y=sweep$path$BIC, type="p",
           xlab="log(lambda)", ylab="BIC") +
  lines(x=log(sweep$path$Lambda),y=sweep$path$BIC)


###################################################
### code chunk number 19: Vignette.Rnw:260-262
###################################################
group <- MGLMtune(Y~0+X, dist="DM", penalty="group", ngridpt=30)
print(group$select)


###################################################
### code chunk number 20: grouppath
###################################################
p2 <- plot(x=log(group$path$Lambda),y=group$path$BIC, type="p",
           xlab="log(lambda)", ylab="BIC") +
  lines(x=log(group$path$Lambda),y=group$path$BIC)


###################################################
### code chunk number 21: Vignette.Rnw:275-276
###################################################
p2 <- plot(x=log(group$path$Lambda),y=group$path$BIC, type="p",
           xlab="log(lambda)", ylab="BIC") +
  lines(x=log(group$path$Lambda),y=group$path$BIC)


###################################################
### code chunk number 22: Vignette.Rnw:288-290
###################################################
nuclear <- MGLMtune(Y~0+X, dist="DM", penalty="nuclear", ngridpt=30)
print(nuclear$select)


###################################################
### code chunk number 23: nuclearpath
###################################################
p3 <- plot(x=log(nuclear$path$Lambda),y=nuclear$path$BIC, type="p",
           xlab="log(lambda)", ylab="BIC") +
  lines(x=log(nuclear$path$Lambda),y=nuclear$path$BIC)


###################################################
### code chunk number 24: Vignette.Rnw:304-305
###################################################
p3 <- plot(x=log(nuclear$path$Lambda),y=nuclear$path$BIC, type="p",
           xlab="log(lambda)", ylab="BIC") +
  lines(x=log(nuclear$path$Lambda),y=nuclear$path$BIC)


