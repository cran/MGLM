\name{dist}
\alias{dist}
\title{
Details of the distributions
}
\description{
An object that specifies the distribution to be fitted by the \code{MGLMfit} function, or the regression model to be fitted by the \code{MGLMreg} or \code{MGLMsparsereg} functions.
Can be chosen from \code{"MN"}, \code{"DM"}, \code{"NegMN"}, or \code{"GDM"}.
}

%- maybe also 'usage' for other objects documented here.
\details{

\subsection{"MN": Multinomial distribution}{
A multinomial distribution models the counts of \eqn{d} possible outcomes.
The counts of categories are negatively correlated. 
The density of a \eqn{d} category count vector \eqn{y} with parameter 
\eqn{p=(p_1, \ldots, p_d)} is
\deqn{
P(y|p) = C_{y_1, \ldots, y_d}^{m} \prod_{j=1}^{d} p_j^{y_j},
}{
P(y|p) = C_{y_1, \ldots, y_d}^{m} prod_{j=1}^{d} p_j^{y_j},
}
where \eqn{m = \sum_{j=1}^d y_j}{m = sum_{j=1}^d y_j}, \eqn{0 < p_j < 1}, and \eqn{\sum_{j=1}^d p_j = 1}{sum_{j=1}^d p_j = 1}. 
Here, \eqn{C_k^n}, often read as "\eqn{n} choose \eqn{k}", refers the number of \eqn{k} combinations from a set of \eqn{n} elements.

The \code{\link{MGLMreg}} function with \code{dist="MN"} calculates the MLE of regression coefficients \eqn{\beta_j} of the multinomial logit model, which has link function \eqn{p_j = exp(X\beta_j)/(1 + \sum_{j=1}^{d-1} exp(X\beta_j))}{p_j = exp(X\beta_j) / (1 + sum_{j=1}^{d-1} exp(X\beta_j))}, \eqn{j=1,\ldots,d-1}. The \code{\link{MGLMsparsereg}} function with \code{dist="MN"} fits regularized multinomial logit model.
}

\subsection{"DM": Dirichlet multinomial distribution}{

When the multivariate count data exhibits over-dispersion, the traditional 
multinomial model is insufficient.  Dirichlet multinomial distribution models the
probabilities of the categories by a Dirichlet distribution.  
The density of a \eqn{d} category count vector \eqn{y}, with 
parameter \eqn{\alpha = (\alpha_1, \ldots, \alpha_d)},
\eqn{\alpha_j > 0}, is
\deqn{ 
  P(y|\alpha) = C_{y_1, \ldots, y_d}^{m} \prod_{j=1}^{d} 
  \frac{\Gamma(\alpha_j+y_j)}{\Gamma(\alpha_j)}
  \frac{\Gamma(\sum_{j'=1}^d \alpha_{j'})}{\Gamma(\sum_{j'=1}^d \alpha_{j'} + \sum_{j'=1}^d y_{j'})},
  }{
  P(y|\alpha) =
  C_{y_1, \ldots, y_d}^{m} prod_{j=1}^d 
  {Gamma(\alpha_j+y_j)Gamma(sum_{j'=1}^d \alpha_j')} / {Gamma(\alpha_j)Gamma(sum_{j'=1}^d \alpha_j' + sum_{j'=1}^d y_j')},
  }
 where \eqn{m=\sum_{j=1}^d y_j}{m = sum_{j=1}^d y_j}. Here, \eqn{C_k^n}, often read as "\eqn{n} choose \eqn{k}", 
 refers the number of \eqn{k} combinations from a set of \eqn{n} elements.

The \code{\link{MGLMfit}} function with \code{dist="DM"} calculates the maximum likelihood estimate (MLE) of \eqn{(\alpha_1, \ldots, \alpha_d)}. The \code{\link{MGLMreg}} function with \code{dist="DM"} calculates the MLE of regression coefficients \eqn{\beta_j} of the Dirichlet multinomial regression model, which has link function \eqn{\alpha_j = exp(X\beta_j)}, \eqn{j=1,\ldots,d}. The \code{\link{MGLMsparsereg}} function with \code{dist="DM"} fits regularized Dirichlet multinomial regression model.
}

\subsection{"GDM": Generalized Dirichlet multinomial distribution}{

The more flexible Generalized Dirichlet multinomial model can be used when the counts of categories have both positive and negative correlations. 
The probability mass of a count vector \eqn{y} over \eqn{m} trials with parameter
\eqn{(\alpha, \beta)=(\alpha_1, \ldots, \alpha_{d-1}, \beta_1, \ldots, \beta_{d-1})},
\eqn{\alpha_j, \beta_j > 0}, is
\deqn{
P(y|\alpha,\beta)
=C_{y_1, \ldots, y_d}^{m} \prod_{j=1}^{d-1} 
\frac{\Gamma(\alpha_j+y_j)}{\Gamma(\alpha_j)}
\frac{\Gamma(\beta_j+z_{j+1})}{\Gamma(\beta_j)}
\frac{\Gamma(\alpha_j+\beta_j)}{\Gamma(\alpha_j+\beta_j+z_j)}  ,
}{
   P(y|\alpha,\beta)
   =C_{y_1, \ldots, y_d}^{m} prod_{j=1}^{d-1} {Gamma(\alpha_j+y_j)Gamma(\beta_j+z_{j+1})Gamma(\alpha_j+\beta_j)} / {Gamma(\alpha_j)Gamma(\beta_j)Gamma(\alpha_j+\beta_j+z_j)},
 }
where \eqn{z_j = \sum_{k=j}^d y_k}{z_j = sum_{k=j}^d y_k} and \eqn{m=\sum_{j=1}^d y_j}{m = sum_{j=1}^d y_j}. Here, \eqn{C_k^n}, often read as "\eqn{n} choose \eqn{k}", 
#' refers the number of \eqn{k} combinations from a set of \eqn{n} elements.

The \code{\link{MGLMfit}} with \code{dist="GDM"} calculates the MLE of \eqn{(\alpha, \beta)=(\alpha_1, \ldots, \alpha_{d-1}, \beta_1, \ldots, \beta_{d-1})}. The \code{\link{MGLMreg}} function with \code{dist="GDM"} calculates the MLE of regression coefficients \eqn{\alpha_j, \beta_j} of the generalized Dirichlet multinomial regression model, which has link functions \eqn{\alpha_j=exp(X\alpha_j)} and \eqn{\beta_j=exp(X\beta_j)}, \eqn{j=1, \ldots, d-1}. The \code{\link{MGLMsparsereg}} function with \code{dist="GDM"} fits regularized generalized Dirichlet multinomial regression model.
}

\subsection{"NegMN": Negative multinomial distribution}{

Both the multinomial distribution and Dirichlet multinomial distribution are good for 
negatively correlated counts.  When the counts of categories are positively 
correlated, the negative multinomial distribution is preferred.  
The probability mass function of a \eqn{d} category count vector \eqn{y} with parameter
 \eqn{(p_1, \ldots, p_{d+1}, \beta)}, \eqn{\sum_{j=1}^{d+1} p_j=1}{sum_{j=1}^{d+1} p_j = 1}, \eqn{p_j > 0}, \eqn{\beta > 0}, is
\deqn{
P(y|p,\beta) =  C_{m}^{\beta+m-1}  C_{y_1, \ldots, y_d}^{m} 
\prod_{j=1}^d p_j^{y_j} p_{d+1}^\beta \\
= \frac{\beta_m}{m!}  C_{y_1, \ldots, y_d}^{m}  
\prod_{j=1}^d p_j^{y_j} p_{d+1}^\beta,
}{
   P(y|p,\beta) =  C_{m}^{\beta+m-1}  C_{y_1, \ldots, y_d}^{m} 
   prod_{j=1}^d p_j^{y_j} p_{d+1}^\beta = (\beta_m)/(m!) C_{y_1, \ldots, y_d}^{m} prod_{j=1}^d p_j^{y_j} p_{d+1}^\beta,
 }
where \eqn{m = \sum_{j=1}^d y_j}{m = sum_{j=1}^d y_j}. Here, \eqn{C_k^n}, often read as "\eqn{n} choose \eqn{k}", refers the number of \eqn{k} combinations from a set of \eqn{n} elements.

The \code{\link{MGLMfit}} function with \code{dist="NegMN"} calculates the MLE of \eqn{(p_1, \ldots, p_{d+1}, \beta)}. The \code{\link{MGLMreg}} function with \code{dist="NegMN"} and \code{regBeta=FALSE} calculates the MLE of regression coefficients \eqn{(\alpha_1,\ldots,\alpha_d, \beta)} of the negative multinomial regression model, which has link function \eqn{p_{d+1} = 1/(1 + \sum_{j=1}^d exp(X\alpha_j))}{p_{d+1} = 1/(1 + sum_{j=1}^d exp(X\alpha_j))}, \eqn{p_j = exp(X\alpha_j) p_{d+1}}, \eqn{j=1, \ldots, d}. When \code{dist="NegMN"} and \code{regBeta=TRUE}, the overdispersion parameter is linked to covariates via \eqn{\beta=exp(X\alpha_{d+1})}, and the 
function \code{MGLMreg} outputs an estimated matrix of 
\eqn{(\alpha_1, \ldots, \alpha_{d+1})}. The \code{\link{MGLMsparsereg}} function with \code{dist="NegMN"} fits regularized negative multinomial regression model.
}
}
\author{
Yiwen Zhang and Hua Zhou
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{MGLMfit}}, \code{\link{MGLMreg}}, \code{\link{MGLMsparsereg}},
\code{\link{dmn}}, \code{\link{ddirmn}}, \code{\link{dgdirmn}}, \code{\link{dnegmn}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Models}
\keyword{Regression}
\keyword{Distribution fitting}
