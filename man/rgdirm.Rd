\name{rgdirm}
\alias{rgdirm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate random observations from the generalized Dirichlet multinomial distribution
}
\description{
Generate random observations from the generalized Dirichlet multinomial distribution. 
}
\usage{
rgdirm(n, size, alpha, beta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
  Number of random vectors to generate.  When \code{size} is a scalar and \code{alpha} is a vector, 
  must specify \code{n}.  When \code{size} is a vector and \code{alpha} is a matrix, \code{n} is optional.
  The default value of \code{n} is the length of \code{size}. If given, \code{n} should be equal to 
  the length of \code{size}. 
}
  \item{size}{
  A number or vector specifying the total number of objects that are put
  into d categories in the generalized Dirichlet multinomial distribution. 
}
  \item{alpha}{
  Numerical positive vector or matrix.  
  If \code{alpha} is a vector, \code{size} must be a scalar.  All the random vectors will
  be drawn from the same \code{alpha} and \code{size}.  If \code{alpha} is a matrix, the 
  number of rows should match the length of \code{size}.  Each random vector 
  will be drawn from the corresponding row of \code{alpha} and the corresponding element of \code{size}. 
}
  \item{beta}{
  Numerical positive vector or matrix.  
  If \code{beta} is a vector, \code{size} must be a scalar.  All the random samples will
  be drawn from the same \code{beta} and \code{size}.  If \code{beta} is a matrix, the 
  number of rows should match the length of \code{size}.  Each random vector 
  will be drawn from the corresponding row of \code{beta} and the corresponding element of \code{size}. 
  See Details below. 
}
}

\value{
A \code{n*d} matrix of the generated random observations.
}

\details{
The random vectors are generated from the generalized Direchlet multinomial
distribution with probability mass function
\deqn{
P(Y_1=y_1, \ldots, Y_d=y_d) =
{m \choose y_1, \ldots, y_d} \prod_{j=1}^{d-1} \frac{\Gamma(\alpha_j+y_j)}{\Gamma(\alpha_j)}
\frac{\Gamma(\beta_j+z_{j+1})}{\Gamma(\beta_j)}\frac{\Gamma(\alpha_j+\beta_j)}{\Gamma(\alpha_j+\beta_j+z_j)},
}
where \eqn{z_j =\sum_{k=j}^d y_k} and \eqn{m = \sum_{j=1}^d y_j}.
}

\author{
Yiwen Zhang and Hua Zhou
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
set.seed(100)
alpha <- matrix(abs(rnorm(40)), 10, 4)
beta <- matrix(abs(rnorm(40)), 10, 4)
size <- rbinom(10, 10, 0.5)
GDM.rdm <- rgdirm(size=size, alpha=alpha, beta=beta)
GDM.rdm1 <- rgdirm(n=20, size=10, alpha=abs(rnorm(4)), beta=abs(rnorm(4)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Random data generator}
