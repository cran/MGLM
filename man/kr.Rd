\name{kr}
\alias{kr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Khatri-Rao product of two matrices
}
\description{
Return the Khatri-Rao product of two matrices, which is a column-wise Kronecker product.
}
\usage{
kr(A, B, w, byrow)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A, B}{
  both A and B are matrices. The two matrices should have the same number of columns.  We also give the user an option to do row-wise Kronecker product, to avoid transpose.  When doing row-wise Kronecker product, the number of rows of A and B should be the same.   
}
  \item{w}{
  the weights vector.  The length of the vector should match with the dimension of the matrices.  If performing column-wise Kronecker product, the length of w should be the same as the column number of A and B.  If performing row-wise Kronecker prodoct, the length of w should be the same as the row number of A and B. The default is a vector of 1 if no value provided.
}
  \item{byrow}{
  a logical variable controling whether to perform row/column-wise Kronecker product.  The default is byrow=TRUE.
}
}
\details{
The column/row-wise Kronecker product. 
}
\value{
A marix of the Khatri-Rao product.
}
\author{
Yiwen Zhang and Hua Zhou
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
X <- matrix(rnorm(30), 10, 3)
Y <- matrix(runif(50), 10, 5)
C <- kr(X, Y)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
