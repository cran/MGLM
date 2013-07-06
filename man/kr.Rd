\name{kr}
\alias{kr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Khatri-Rao product of two matrices
}
\description{
Return the Khatri-Rao product of two matrices.
}
\usage{
kr(A, B)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A, B}{
  both A and B are matrices, with the same number of columns.  
}
}
\details{
The column-wise Kronecker product. 
}
\value{
A marix of the Khatri-Rao product.
}
\author{
Yiwen Zhang and Hua Zhou
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
X <- matrix(rnorm(30), 3, 10)
Y <- matrix(runif(50), 5, 10)
C <- kr(X, Y)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
