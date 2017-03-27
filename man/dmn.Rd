\name{dmn}
\alias{dmn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the log of the multinomial probability mass function
}
\description{
Calculate the log of the multinomial probability mass function.
}
\usage{
dmn(Y, prob)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
The multivariate count matrix with dimension \eqn{n \times d}, where 
\eqn{n = 1,2, \ldots} is number of observations and \eqn{d=2,\ldots} is number of categories.  
}
  \item{prob}{
The probability parameter of the multinomial distribution.  \code{prob}
can be either a vector of length \eqn{d} or a matrix with matching 
size of \code{Y}.  If \code{prob} is a vector, it will be replicated \eqn{n} 
times to match the dimension of \code{Y}. If the sum(s) of \code{prob} is not 1, it will be automatically scaled to have sum 1.
}
}
\details{

\eqn{y=(y_1, \ldots, y_d)} is a \eqn{d} category count vector. Given the parameter vector \eqn{p = (p_1, \ldots, p_d)}, \eqn{0 < p_j < 1}, \eqn{\sum_{j=1}^d p_j = 1}, the function calculates the log of the multinomial pmf
\deqn{
P(y|p) = {m \choose y_1, \ldots, y_d} \prod_{j=1}^{d} p_j^{y_j},
}
where \eqn{m=\sum_{j=1}^d y_j}.
 
For each count vector and each corresponding parameter vector
\eqn{p}, the function \code{dmn} returns the value of \eqn{\log(P(y|p))}. 
When \code{Y} is a matrix of \eqn{n} rows, the function returns a 
vector of length \eqn{n}. 
 
The parameter \eqn{p} can be one vector, like the result from the distribution
fitting function; or, \eqn{p} can be a matrix with \eqn{n} rows, like the estimate
from the regression function, 
\eqn{p_j = \frac{exp(X \beta_j)}{1 + \sum_{j'=1}^{d-1} exp(X\beta_{j'})}}, \eqn{j=1,\ldots,d-1},
\eqn{p_d = \frac{1}{1 + \sum_{j'=1}^{d-1} exp(X\beta_{j'})}}.
The \eqn{d}-th column of the coefficient matrix \eqn{\beta} is set to \eqn{0} to avoid the identifiability issue.
}

\value{
A vector of length \eqn{n} with value \eqn{\log(P(y|p))}
}
\author{
Yiwen Zhang and Hua Zhou
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
m <- 20
prob <- c(0.1, 0.2)
dm.Y <- rdirm(n=10, m, prob)	
pdfln <- dmn(dm.Y, prob)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{log pdf}
