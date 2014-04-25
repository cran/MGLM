### R code from vignette source 'Vignette.Rnw'

###################################################
### code chunk number 1: Vignette.Rnw:38-41
###################################################
require(ggplot2)
require(reshape2)
require(plyr)


###################################################
### code chunk number 2: Vignette.Rnw:43-44
###################################################
require(MGLM)


###################################################
### code chunk number 3: Vignette.Rnw:49-51
###################################################
data(iris)
Y <-iris[, 1:4]


###################################################
### code chunk number 4: scatter1
###################################################
scatter.plot(Y, facet=TRUE, free=TRUE)


###################################################
### code chunk number 5: corr
###################################################
corr.plot(Y)


###################################################
### code chunk number 6: scatter1
###################################################
scatter.plot(Y, facet=TRUE, free=TRUE)


###################################################
### code chunk number 7: corr
###################################################
corr.plot(Y)


###################################################
### code chunk number 8: Vignette.Rnw:94-100
###################################################
set.seed(123)
n <- 200
d <- 4
alpha <- rep(1,d)
m <- 50
Y <- rmn(m, alpha, n)


###################################################
### code chunk number 9: Vignette.Rnw:103-105
###################################################
mnFit <- MGLMfit(Y, dist="DM")
print(mnFit)


###################################################
### code chunk number 10: Vignette.Rnw:112-118
###################################################
set.seed(123)
n <- 200
d <- 4
alpha <- rep(1, d)
m <- 50
Y <- rdirm(m, alpha, n)


###################################################
### code chunk number 11: Vignette.Rnw:120-122
###################################################
dmFit <- MGLMfit(Y, dist="DM")
print(dmFit)


###################################################
### code chunk number 12: Vignette.Rnw:129-131
###################################################
gdmFit <- MGLMfit(Y, dist="GDM")
print(gdmFit)


###################################################
### code chunk number 13: Vignette.Rnw:139-148
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
### code chunk number 14: Vignette.Rnw:154-163
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
### code chunk number 15: Vignette.Rnw:181-194
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
### code chunk number 16: Vignette.Rnw:200-202
###################################################
mnReg <- MGLMreg(Y~0+X, dist="MN")
print(mnReg)


###################################################
### code chunk number 17: Vignette.Rnw:209-211
###################################################
dmReg <- MGLMreg(Y~0+X, dist="DM")
print(dmReg)


###################################################
### code chunk number 18: Vignette.Rnw:218-220
###################################################
gdmReg <- MGLMreg(Y~0+X, dist="GDM")
print(gdmReg)


###################################################
### code chunk number 19: Vignette.Rnw:227-229
###################################################
negReg <- MGLMreg(Y~0+X, dist="NegMN", regBeta=FALSE)
print(negReg)


###################################################
### code chunk number 20: fit1
###################################################
plot(gdmReg, facet=TRUE, free=TRUE)


###################################################
### code chunk number 21: Vignette.Rnw:243-246
###################################################
newX <- matrix(runif(1*p), 1, p)
pred <- predict(gdmReg, newX)
pred


###################################################
### code chunk number 22: Vignette.Rnw:262-272
###################################################
set.seed(118)
n <- 100
p <- 10
d <- 5
m <- rbinom(n, 200, 0.8)
X <- matrix(rnorm(n*p),n, p)
alpha <- matrix(0, p, d)
alpha[c(1,3, 5), ] <- 1
Alpha <- exp(X%*%alpha)
Y <- rdirm(size=m, alpha=Alpha)


###################################################
### code chunk number 23: Vignette.Rnw:279-281
###################################################
sweep <- MGLMtune(Y~0+X, dist="DM", penalty="sweep", ngridpt=30)
print(sweep$select)


###################################################
### code chunk number 24: sweeppath
###################################################
plot(sweep)


###################################################
### code chunk number 25: Vignette.Rnw:293-294
###################################################
plot(sweep)


###################################################
### code chunk number 26: Vignette.Rnw:305-307
###################################################
group <- MGLMtune(Y~0+X, dist="DM", penalty="group", ngridpt=30)
print(group$select)


###################################################
### code chunk number 27: grouppath
###################################################
plot(group)


###################################################
### code chunk number 28: Vignette.Rnw:318-319
###################################################
plot(group)


###################################################
### code chunk number 29: Vignette.Rnw:331-333
###################################################
nuclear <- MGLMtune(Y~0+X, dist="DM", penalty="nuclear", ngridpt=30, warm.start=FALSE)
print(nuclear$select)


###################################################
### code chunk number 30: nuclearpath
###################################################
plot(nuclear)


###################################################
### code chunk number 31: Vignette.Rnw:345-346
###################################################
plot(nuclear)


