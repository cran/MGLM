## All classes & generics 

setClassUnion("numericOrMatrix", c("numeric", "matrix"))
setClassUnion("listOrMatrix", c("list", "matrix"))

#' @name MGLMfit-class
#' @aliases MGLMfit-class  
#' @title Class \code{"MGLMfit"}
#' @description A class containing the model fitting results from the \code{MGLMfit}.
#' @docType class 
#' 
#' @slot 	estimate object of class \code{"vector"}, containing the parameter estimates.
#' @slot 	SE object of class \code{"vector"},
#'   containing the standard errors of the estimates.
#' @slot 	vcov object of class \code{"matrix"},
#'   the variance covariance matrix of the parameter estimates.
#' @slot 	logL object of class \code{"numeric"}, 
#'   the fitted log likelihood. 
#' @slot 	BIC object of class \code{"numeric"}, 
#'   Bayesian information criterion.
#' @slot 	AIC object of class \code{"numeric"},
#'   Akaike information criterion.
#' @slot 	LRTpvalue object of class \code{"numeric"},
#'   likelihood ratio test p value.
#' @slot 	gradient object of class \code{"numeric"} or \code{"matrix"},
#'   containing the gradient.
#' @slot 	iter object of class \code{"numeric"}, 
#'   number of iteration used.
#' @slot 	distribution object of class \code{"character"},
#'   the distribution fitted.
#' @slot 	fitted object of class \code{"vector"},
#'   the fitted mean of each category.
#' @slot LRT object of class \code{"numeric"}, 
#' the likelihood ratio test statistic. 
#' 
#' @examples 
#' showClass("MGLMfit")
#' 
#' @exportClass MGLMfit 
#' @author Yiwen Zhang and Hua Zhou 
#' @keywords classes 
setClass("MGLMfit", representation(estimate = "vector", SE = "vector", vcov = "matrix", 
                                   logL = "numeric", BIC = "numeric", AIC = "numeric", LRT = "numeric",
                                   LRTpvalue = "numeric", iter = "numeric", distribution = "character", 
                                   gradient = "numericOrMatrix", fitted = "vector"))



#' @name MGLMreg-class
#' @aliases MGLMreg-class
#' @docType class
#' @title Class \code{"MGLMreg"}
#' @description Objects can be created by calls of the form \code{new("MGLMreg", ...)}.
#' 
#' @slot call object of class \code{"call"}.
#' @slot data object of class \code{"list"} ,
#' consists of both the predictor matrix and the response matrix. 
#' @slot coefficients object of class \code{"list"} or \code{"matrix"},
#' the estimated parameters.
#' @slot SE object of class \code{"list"} or \code{"matrix"},
#' the standard errors of the parameters.
#' @slot test object of class \code{"matrix"},
#' the test statistics and p-values.
#' @slot Hessian object of class \code{"matrix"},
#' the Hessian matrix.
#' @slot logL object of class \code{"numeric"},
#' the loglikelihood.
#' @slot BIC object of class \code{"numeric"},
#" Bayesian information criterion.
#' @slot AIC object of class \code{"numeric"},
#' Akaike information criterion.
#' @slot iter object of class \code{"numeric"}, 
#'   the number of iteration used.
#' @slot distribution object of class \code{"character"},
#'   the distribution fitted.
#' @slot fitted object of class \code{"vector"},
#'   the fitted value.
#' @slot gradient object of class \code{"numeric"} or \code{"matrix"},
#'  the gradient at the estimated parameter values.
#' @slot wald.value object of class \code{"numeric"},
#'  the Wald statistics.
#' @slot wald.p object of class \code{"numeric"},
#'   the p values of Wald test.
#' @slot Dof object of class \code{"numeric"}, 
#'   the degrees of freedom. 
#' 
#' @examples 
#' showClass("MGLMreg")
#' 
#' @exportClass MGLMreg 
#' @author Yiwen Zhang and Hua Zhou 
#' @keywords classes 
setClass("MGLMreg", representation(coefficients = "listOrMatrix", SE = "listOrMatrix", 
                                   Hessian = "matrix", gradient = "numericOrMatrix", wald.value = "numeric", 
                                   wald.p = "numeric", test = "matrix", logL = "numeric", BIC = "numeric", 
                                   AIC = "numeric", fitted = "matrix", call = "call", iter = "numeric",
                                   distribution = "character", data = "list", Dof = "numeric"))



#' @name MGLMsparsereg-class 
#' @title Class \code{"MGLMsparsereg"}
#' @aliases MGLMsparsereg-class 
#' @description A class containing the results from the \code{MGLMsparsereg}.
#' @docType class 
#' 
#' @slot call object of class \code{"call"}.
#' @slot data object of class \code{"list"} ,
#' consists of both the predictor matrix and the response matrix. 
#' @slot coefficients object of class \code{"matrix"},
#' the estimated parameters.
#' @slot logL object of class \code{"numeric"},
#' the loglikelihood.
#' @slot BIC object of class \code{"numeric"},
#" Bayesian information criterion.
#' @slot AIC object of class \code{"numeric"},
#' Akaike information criterion.
#' @slot Dof object of class \code{"numeric"},
#'  the degrees of freedom.
#' @slot 	iter object of class \code{"numeric"}, 
#'   the number of iteration used.
#' @slot maxlambda object of class \code{"numeric"},
#' the maximum tuning parameter that ensures the estimated regression coefficients are not all zero.
#' @slot lambda object of class \code{"numeric"},
#' the tuning parameter used. 
#' @slot 	distribution object of class \code{"character"},
#'   the distribution fitted.
#' @slot penalty Object of class \code{"character"},
#' the chosen penalty when running penalized regression.
#' 
#' @examples 
#' showClass("MGLMsparsereg")
#' @author Yiwen Zhang and Hua Zhou
#' @exportClass MGLMsparsereg 
#' @keywords classes 
setClass("MGLMsparsereg", representation(call = "call", data = "list", coefficients = "listOrMatrix", 
                                         logL = "numeric", BIC = "numeric", AIC = "numeric", Dof = "numeric", iter = "numeric", 
                                         maxlambda = "numeric", lambda = "numeric", distribution = "character", penalty = "character"))
#     Beta = "numeric"))




#' @name MGLMtune-class
#' @aliases MGLMtune-class  
#' @title Class \code{"MGLMtune"}
#' @description A class containing the results from the \code{MGLMtune}.
#' 
#' @slot call object of class \code{"call"}.
#' @slot select object of class \code{"MGLMsparsereg"},
#' regularized regression results given by the optimal tuning parameter.
#' @slot path object of class \code{"data.frame"},
#' the BIC, AIC, log-likelihood and degrees of freedom given each tuning parameter.
#' @slot select.list object of class \code{"list"},
#' the regularized regression results at each tuning grid point.
#'
#' @examples 
#' showClass("MGLMtune")
#' 
#' @author Yiwen Zhang and Hua Zhou
#' @keywords classes
#' @exportClass MGLMtune 
setClass("MGLMtune", representation(call = "call", select = "MGLMsparsereg", 
                                    path = "data.frame", select.list = "list"))
