\name{dgdirm}
\alias{dgdirm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the log of the generalized Dirichlet multinomial probability mass function
}
\description{
Calculate the log of the generalized Dirichlet multinomial probability mass function.
}
\usage{
dgdirm(Y, alpha, beta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
The multivariate count matrix with dimensions \eqn{n \times d}, where 
\eqn{n = 1,2, \ldots} is the number of observations and \eqn{d=3,4,\ldots} is the number of categories.
}
  \item{alpha}{
The parameter of the generalized Dirichlet multinomial distribution.  \code{alpha} can
either be a vector or a matrix which matches the size of \code{Y}.  If \code{alpha} 
is a vector, it will be replicated \eqn{n} times to match the dimension of \code{Y}. 
}
  \item{beta}{
The parameter of the generalized Dirichlet multinomial distribution. \code{beta} should
have the same dimension as \code{alpha}.
}
}

\details{
\eqn{Y=(y_1, \ldots, y_d)} are the \eqn{d} category count vectors. Given the parameter vector \eqn{\alpha = (\alpha_1, \ldots, \alpha_{d-1}),
 \alpha_j>0}, and \eqn{\beta=(\beta_1, \ldots, \beta_{d-1}), \beta_j>0},
 it calculates the generalized Dirichlet multinomial pmf
 \deqn{
P(y|\alpha,\beta)
={m \choose y_1, \ldots, y_d} \prod_{j=1}^{d-1} 
\frac{\Gamma(\alpha_j+y_j)}{\Gamma(\alpha_j)}
\frac{\Gamma(\beta_j+z_{j+1})}{\Gamma(\beta_j)}
\frac{\Gamma(\alpha_j+\beta_j)}{\Gamma(\alpha_j+\beta_j+z_j)}  ,
}
where \eqn{z_j = \sum_{k=j}^d y_k} and \eqn{m = \sum_{j=1}^d y_j}.

For each count vector and each corresponding parameter vector
\eqn{\alpha} and \eqn{\beta}, the function \code{dgdirm} returns the value of 
\eqn{\log(P(y|\alpha, \beta))}. 
When \code{Y} is a matrix of \eqn{n} rows, the function returns a vector of length \eqn{n}. 

The \eqn{\alpha} and \eqn{\beta} parameters can be vectors, like the results from the 
distribution
fitting function, or they can be matrices with \eqn{n} rows, 
like the estimate
from the regression function multiplied by the covariate matrix
\eqn{exp(X\alpha)} and \eqn{exp(X\beta)}
}
\value{
A vector of length \eqn{n} with value \eqn{\log(P(y| alpha, beta)}
}
\author{
Yiwen Zhang and Hua Zhou
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
m <- 20
alpha <- c(0.2, 0.5)
beta <- c(0.7, 0.4)
Y <- rgdirm(10, m, alpha, beta)
dgdirm(Y, alpha, beta)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ log pdf }% __ONLY ONE__ keyword per line