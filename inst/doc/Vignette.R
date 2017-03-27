### R code from vignette source 'Vignette.Rnw'

###################################################
### code chunk number 1: Vignette.Rnw:46-53
###################################################
require(MGLM)
set.seed(123)
n <- 200
d <- 4
alpha <- rep(1, d)
m <- 50
Y <- rmn(n, m, alpha)


###################################################
### code chunk number 2: Vignette.Rnw:56-58
###################################################
mnFit <- MGLMfit(Y, dist="MN")
print(mnFit)


###################################################
### code chunk number 3: Vignette.Rnw:62-64
###################################################
compareFit <- MGLMfit(Y, dist="DM")
print(compareFit)


###################################################
### code chunk number 4: Vignette.Rnw:71-77
###################################################
set.seed(123)
n <- 200
d <- 4
alpha <- rep(1, d)
m <- 50
Y <- rdirm(n, m, alpha)


###################################################
### code chunk number 5: Vignette.Rnw:79-81
###################################################
dmFit <- MGLMfit(Y, dist="DM")
print(dmFit)


###################################################
### code chunk number 6: Vignette.Rnw:88-90
###################################################
gdmFit <- MGLMfit(Y, dist="GDM")
print(gdmFit)


###################################################
### code chunk number 7: Vignette.Rnw:98-107
###################################################
set.seed(124)
n <- 200
d <- 4
alpha <- rep(1, d-1)
beta <- rep(1, d-1)
m <- 50
Y <- rgdirm(n, m, alpha, beta)
gdmFit <- MGLMfit(Y, dist="GDM")
print(gdmFit)


###################################################
### code chunk number 8: Vignette.Rnw:113-122
###################################################
set.seed(1220)
n <- 100
d <- 4
p <- 5
prob <- rep(0.2, d)
beta <- 10
Y <- rnegmn(n, prob, beta)
negmnFit <- MGLMfit(Y, dist="NegMN")
print(negmnFit)


###################################################
### code chunk number 9: Vignette.Rnw:140-153
###################################################
set.seed(1234)
n <- 200
p <- 5
d <- 4
X <- matrix(runif(p * n), n, p)
alpha <- matrix(c(0.6, 0.8, 1), p, d - 1, byrow=TRUE)
alpha[c(1, 2),] <- 0
Alpha <- exp(X %*% alpha)
beta <- matrix(c(1.2, 1, 0.6), p, d - 1, byrow=TRUE)
beta[c(1, 2),] <- 0
Beta <- exp(X %*% beta)
m <- runif(n, min=0, max=25) + 25
Y <- rgdirm(n, m, Alpha, Beta)


###################################################
### code chunk number 10: Vignette.Rnw:159-161
###################################################
mnReg <- MGLMreg(Y~0+X, dist="MN")
print(mnReg)


###################################################
### code chunk number 11: Vignette.Rnw:167-169
###################################################
dmReg <- MGLMreg(Y~0+X, dist="DM")
print(dmReg)


###################################################
### code chunk number 12: Vignette.Rnw:175-177
###################################################
gdmReg <- MGLMreg(Y~0+X, dist="GDM")
print(gdmReg)


###################################################
### code chunk number 13: Vignette.Rnw:183-185
###################################################
negReg <- MGLMreg(Y~0+X, dist="NegMN", regBeta=FALSE)
print(negReg)


###################################################
### code chunk number 14: Vignette.Rnw:198-201
###################################################
newX <- matrix(runif(1*p), 1, p)
pred <- predict(gdmReg, newX)
pred


###################################################
### code chunk number 15: Vignette.Rnw:217-227
###################################################
set.seed(118)
n <- 100
p <- 10
d <- 5
m <- rbinom(n, 200, 0.8)
X <- matrix(rnorm(n * p), n, p)
alpha <- matrix(0, p, d)
alpha[c(1, 3, 5), ] <- 1
Alpha <- exp(X %*% alpha)
Y <- rdirm(size=m, alpha=Alpha)


###################################################
### code chunk number 16: Vignette.Rnw:232-234
###################################################
sweep <- MGLMtune(Y ~ 0 + X, dist="DM", penalty="sweep", ngridpt=30)
print(sweep$select)


###################################################
### code chunk number 17: Vignette.Rnw:256-258
###################################################
group <- MGLMtune(Y ~ 0 + X, dist="DM", penalty="group", ngridpt=30)
print(group$select)


###################################################
### code chunk number 18: Vignette.Rnw:282-284
###################################################
nuclear <- MGLMtune(Y ~ 0 + X, dist="DM", penalty="nuclear", ngridpt=30, warm.start=FALSE)
print(nuclear$select)


