\name{internals}
\alias{dmultn}
\alias{glm.private}
\alias{matrix_threshold}
\alias{lsq_threshold}
\alias{svt}
\alias{lsq_thresholding}
\alias{objfun}
\alias{objfun.grad}
\alias{objfun.hessian}
\alias{MGLM.loss}

\title{Internal Functions}
\description{
These functions are for internal use only or not yet documented.
}
\usage{
dmultn(X,Y,B1, weight)
glm.private(Y, start=NULL, weights, X, family)
matrix_threshold(X, lambda, penalty)
lsq_thresholding(b, lambda)
svt(b, lambda)
objfun(alpha, x, y, d, p)
objfun.grad(alpha, x, y, d, p)
objfun.hessian(alpha, x, y, d, p)
MGLM.loss(Y, X, beta, dist, weight, regBeta = FALSE, Beta) 
}

\author{Yiwen Zhang and Hua Zhou}

\keyword{internal}

