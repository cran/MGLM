\name{scatter.plot}
\alias{scatter.plot}
\title{
Pairwise scatter plot matrices
}
\description{
A function produces pairwise scatter plot matrices.
}
\usage{
scatter.plot(Y, main = NULL, free = TRUE, ...)
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
  \item{free}{
    an optional logical variable controling whether the scales shared across all
	small plot blocks.  The default is \code{free=TRUE}. 
}
  \item{\dots}{
	other arguments passed on to \code{geom_point}.
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
plot <- scatter.plot(iris[, 1:4])
print(plot)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ plot }
