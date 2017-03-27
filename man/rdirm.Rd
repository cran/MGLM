\name{rdirm}
\alias{rdirm}
\title{
Generate random observations from the Dirichlet multinomial distribution
}
\description{
Generate random observations from the Dirichlet multinomial distribution.
}
\usage{
rdirm(n, size, alpha)
}
\arguments{
  \item{n}{
  Number of random vectors to generate.  When \code{size} is a scalar and \code{alpha} is a vector, 
  must specify \code{n}.  When \code{size} is a vector and \code{alpha} is a matrix, \code{n} is optional.
  The default value of \code{n} is the length of \code{size}. If given, \code{n} should be equal to 
  the length of \code{size}. 
}
  \item{size}{
  A number or vector specifying the total number of objects that are put
  into d categories in the Dirichlet multinomial distribution. 
}
  \item{alpha}{
  Numerical positive vector or matrix.  
  If \code{alpha} is a vector, \code{size} must be a scalar, and all the random vectors will
  be drawn from the same \code{alpha} and \code{size}.  
  If \code{alpha} is a matrix, the number of rows should match the length of 
  \code{size}, and each random vector 
  will be drawn from the corresponding row of \code{alpha} and the corresponding
  element in the \code{size} vector. See Details below. 
}
}

\details{
The random vector is generated from the Dirichlet multinomial
distribution with probability mass function 
\deqn{
P(Y_1=y_1, \ldots, Y_d=y_d) = 
{m \choose y_1, \ldots, y_d}  \prod_{j=1}^{d} 
\frac{\Gamma(\alpha_j+y_j)}{\Gamma(\alpha_j)} \frac{\Gamma(\sum_{j'=1}^d \alpha_{j'})}{\Gamma(\sum_{j'=1}^d \alpha_{j'} + \sum_{j'=1}^d y_{j'})}
}
where \eqn{m = \sum_{j=1}^d y_j}.
}
\value{
a \code{n*d} matrix of the generated random observations.
}
\author{
Yiwen Zhang and Hua Zhou
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
set.seed(100)
alpha <- matrix(abs(rnorm(50)), 10, 5)
size <- rbinom(10, 10, 0.5)
DM.rdm <- rdirm(size=size, alpha=alpha)
DM.rdm1 <- rdirm(n=20, size=20, alpha=c(0.2, 0.1, 0.5))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Random data generator}
