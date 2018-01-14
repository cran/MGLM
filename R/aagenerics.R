if (!isGeneric("logLik"))
  setGeneric("logLik", function(object, ...)
    standardGeneric("logLik"), package = "MGLM")


if (!isGeneric("predict")) 
  setGeneric("predict", function(object, newdata) 
    standardGeneric("predict"), package = "MGLM")


if (!isGeneric("BIC"))
  setGeneric("BIC", function(object, ...)
    standardGeneric("BIC"), package = "MGLM")


if (!isGeneric("AIC"))
  setGeneric("AIC", function(object, ...)
    standardGeneric("AIC"), package = "MGLM")


if (!isGeneric("coef"))
  setGeneric("coef", function(object, ...)
    standardGeneric("coef"), package = "MGLM")


if (!isGeneric("path"))
  setGeneric("path", function(object, ...)
    standardGeneric("path"), package = "MGLM")


if (!isGeneric("maxlambda"))
  setGeneric("maxlambda", function(object, ...)
    standardGeneric("maxlambda"), package = "MGLM")


if (!isGeneric("dof"))
  setGeneric("dof", function(object, ...)
    standardGeneric("dof"), package = "MGLM")
