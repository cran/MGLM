\name{corr.plot}
\alias{corr.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pairwise correlation heatmap
}
\description{
A funciton creates pairwise correlation plot.
}
\usage{
corr.plot(Y, main = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
	a data frame or matrix of which the pairwise scatter plots of every two columns 
	are produced. 
}
  \item{main}{
	an optional character string specifying the title for the plot.
}
}
\value{
A plot object of class \code{ggplot2}.
}
\references{
http://docs.ggplot2.org/current/
}
\author{
Yiwen Zhang and Hua Zhou
}

\examples{
data(iris)
plot <- corr.plot(iris[, 1:4])
print(plot)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ plot }
