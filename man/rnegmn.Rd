\name{rnegmn}
\alias{rnegmn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate random observations from the negative multinomial distribution
}
\description{
Generate random observations from the negative multinomial distribution.
}
\usage{
rnegmn(n, prob, beta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
  Number of random vectors to generate.  When \code{beta} is a scalar and \code{prob} is a vector, 
  must specify \code{n}.  When \code{beta} is a vector and \code{prob} is a matrix, \code{n} is optional.
  The default value of \code{n} is the length of \code{beta}. If given, \code{n} should be equal to 
  the length of \code{beta}. 
}
  \item{prob}{
Numerical non-negative vector or matrix.  If \code{prob} is a vector, \code{beta}
must be a scalar.  All the \code{n} random vectors will be drawn from the
same \code{prob} and \code{beta}.  If \code{prob} is a matrix, the number of rows should 
match the length of \code{beta}.  Each random vector will be drawn from
the corresponding row of \code{prob} and the corresponding element of \code{beta}. Each row of \code{prob} should have sum less than 1.
}
  \item{beta}{
Number or vector specifying the over despersion in negative multinomial
distribution. 
}
}


\details{
The random vectors are generated from the negative multinomial
distribution NegMN\eqn{(p_1, p_2, \ldots, p_d, \beta)} with probability mass function
\deqn{
P(Y_1=y_1, \ldots, Y_d=y_d) =
{\beta +m-1 \choose m} {m \choose y_1, \ldots, y_d} \prod_{j=1}^d p_j^{y_j} p_{d+1}^{\beta} 
}
where \eqn{p_{d+1}=1 - \sum_{j=1}^d p_j} and \eqn{m =\sum_{j=1}^d y_j}.
}

\value{
A \code{n*d} matrix of the generated random observations.
}

\author{
Yiwen Zhang and Hua Zhou
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
set.seed(128)
n <- 100
d <- 4
p <- 5
a <- -matrix(1,p,d)
X <- matrix(runif(n*p), n, p )
alpha <- exp(X\%*\%a)
prob <- alpha/(rowSums(alpha)+1)
beta <- exp(X\%*\%matrix(1,p)) 
Y <- rnegmn(n, prob, beta)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Random data generator}

