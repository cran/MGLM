\name{ddirm}
\alias{ddirm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate log of Dirichlet multinomial probability mass function
}
\description{
Calculate log of Dirichlet multinomial probability mass function.
}
\usage{
ddirm(Y, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
the multivariate count matrix with dimension \eqn{n \times d}, where 
\eqn{n = 1,2, \ldots} is number of observations and \eqn{d=2,3, \ldots} is number of categories. 
}
  \item{alpha}{
the parameter of the Dirichlet multinomial distribution.  \code{alpha}
can either be a vector or a matrix with matching size of \code{Y}.  If \code{alpha} 
is a vector, it will be replicated \eqn{n} times to match the dimension of \code{Y}. 
}
}
\details{
\eqn{Y=(y_1, \ldots, y_d)}, \eqn{d \ge 2}, is \eqn{d} category count vector(s). Given the parameter vector \eqn{\alpha = (\alpha_1, \ldots, \alpha_d), \alpha_j>0  }, it calculates the log of the Dirichlet multinomial pmf
\deqn{
P(y|\alpha) =
 {m \choose y_1, \ldots, y_d} \prod_{j=1}^{d} 
 \frac{\Gamma(\alpha_j+y_j)}{\Gamma(\alpha_j)}
  \frac{\Gamma(\sum_{j'=1}^d \alpha_j)}{\Gamma(\sum_{j'=1}^d \alpha_j + \sum_{j'=1}^d y_j)},
  }
where \eqn{m=\sum_{j=1}^d y_j}.
 
For each count vector and each corresponding parameter vector
\eqn{\alpha}, the function \code{ddirm} returns the value of \eqn{\log(P(y|\alpha))}. 
When \code{Y} is a matrix of \eqn{n} rows, the function returns a vector of length \eqn{n}. 

The parameter \eqn{\alpha} can be a vector of length \eqn{d}, 
such as the results from the distribution fitting; 
or, \eqn{\alpha} can be a matrix with \eqn{n} rows, such as the inverse link  
calculated from the regression parameter estimate \eqn{exp(X \beta) }
}
\value{
A vector of length \eqn{n} with value \eqn{\log(P(y| \alpha))}
}
\author{
Yiwen Zhang and Hua Zhou
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
m <- 20
alpha <- c(0.1, 0.2)
dm.Y <- rdirm(m, alpha, n=10)	
pdfln <- ddirm(dm.Y, alpha)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{log pdf}
