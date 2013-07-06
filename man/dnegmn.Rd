\name{dnegmn}
\alias{dnegmn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate log of negative multinomial probability mass function
}
\description{
Calculate log of negative multinomial probability mass function.
}
\usage{
dnegmn(Y, prob, beta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
the multivariate response matrix of dimension \eqn{n \times d}, 
where \eqn{n = 1, 2, \ldots} is number of observations and \eqn{d=2,3,\ldots} is number of categories. 
}
  \item{prob}{
the probability parameter of the 
negative multinomial distribution.  \code{prob} can be either
a vector of length \eqn{d} \eqn{(d \ge 2)} or a matrix with matching size of \code{Y}.  
If \code{prob} is a vector, it will 
be replicated \eqn{n} times to match the dimension of \code{Y}.  The sum
of each row of \code{prob} should be smaller than 1. 
}
  \item{beta}{
the over dispersion parameter of the negative multinomial distribution. 
\code{beta} can be either a scalar or a vector of length \eqn{n}.
}
}
\details{
\eqn{y=(y_1, \ldots, y_d)} is a \eqn{d} category vector. Given the parameter vector \eqn{p= (p_1, \ldots, p_d)},
\eqn{p_{d+1} = 1/(1 + \sum_{j'=1}^d p_{j'})}, 
 and \eqn{\beta}, \eqn{\beta>0}, it calculates the log of the negative multinomial pmf
\deqn{
P(y|p,\beta) =  {\beta+m-1 \choose m}  {m \choose y_1, \ldots, y_d} 
\prod_{j=1}^d p_j^{y_j} p_{d+1}^\beta = \frac{\beta_m}{m!}  {m \choose y_1, \ldots, y_d} \prod_{j=1}^d p_j^{y_j} p_{d+1}^\beta,
}
where \eqn{m = \sum_{j=1}^d y_j}.
  
For each negative multinomial sample, and each corresponding parameter \eqn{\alpha}
and \eqn{\beta}, the function \code{dnegmn} returns the value of 
\eqn{\log(P(y|p, \beta) )}.  When \code{Y} is a matrix of \eqn{n} rows, the function
returns a vector of length \eqn{n}. 
 
The parameter \code{prob} can be a vector and \code{beta} is a scalar; or, \code{prob} can
be a matrix with \eqn{n} rows, and \code{beta} is a vector of length \eqn{n}
like the estimate from the regression function
multiplied by the covariate matrix
}
\value{
A vector of length \eqn{n} with value \eqn{\log(P(y|p, \beta))}
}
\author{
Yiwen Zhang and Hua Zhou
}
\examples{
\donttest{
m <- 20
n <- 10
p <- 5
d <- 6
a <- -matrix(1,p,d)
X <- matrix(runif(n*p), n, p )
alpha <- exp(X\%*\%a)
prob <- alpha/(rowSums(alpha)+1)
b <- exp(X\%*\%rep(0.3,p)) 
Y <- rnegmn(prob=prob, beta=rep(10, n))
dnegmn(Y, prob, b)
}
}
\keyword{ log pdf }
