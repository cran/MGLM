\name{internals}
\alias{MGLM.loss}
\alias{dmultn}
\alias{dneg}
\alias{glm.private}
\alias{matrix_threshold}
\alias{lsq_threshold}
\alias{rmn}
\alias{svt}
\alias{lsq_thresholding}

\title{Internal Functions}
\description{
These functions are for internal use only or not yet documented.
}
\usage{
MGLM.loss(Y, X, beta, dist, weight, regBeta=FALSE, Beta)
dmultn(X,Y,B1, weight)
dneg(Y, alpha, beta)
glm.private(Y, start=NULL, weights, X, family)
matrix_threshold(X, lambda, penalty)
lsq_thresholding(b, lambda)
rmn(size, alpha, n)
svt(b, lambda)
}

\author{Yiwen Zhang and Hua Zhou}

\keyword{internal}

