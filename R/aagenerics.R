if (!isGeneric("path"))
  setGeneric("path", function(object, ...)
    standardGeneric("path"), package = "MGLM")

if (!isGeneric("maxlambda"))
  setGeneric("maxlambda", function(object, ...)
    standardGeneric("maxlambda"), package = "MGLM")

if (!isGeneric("dof"))
  setGeneric("dof", function(object, ...)
    standardGeneric("dof"), package = "MGLM")
