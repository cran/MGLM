# Kronecker product 
# Author: Yiwen Zhang Date: Modified on 02/02/2014
##============================================================## 

#' @name kr 
#' @aliases kr
#' @title Khatri-Rao product of two matrices
#' @description Return the Khatri-Rao product of two matrices, which is a column-wise Kronecker product.
#' @param A,B matrices. The two matrices \code{A} and \code{B} should have the same number of columns.  
#' We also give the user an option to do row-wise Kronecker product, to avoid transpose.  
#' When doing row-wise Kronecker product, the number of rows of A and B should be the same.
#' @param w the weights vector. The length of the vector should match with the dimension of the matrices.  
#' If performing column-wise Kronecker product, the length of w should be the same as the column number of A and B.  
#' If performing row-wise Kronecker prodoct, the length of w should be the same as the row number of A and B. 
#' The default is a vector of 1 if no value provided.
#' @param byrow a logical variable controlling whether to perform row/column-wise Kronecker product.  
#' The default is \code{byrow}=TRUE. 
#' @details The column/row-wise Kronecker product. 
#' @return A matrix of the Khatri-Rao product.
#' @author Yiwen Zhang and Hua Zhou
#' 
#' @examples 
#' X <- matrix(rnorm(30), 10, 3)
#' Y <- matrix(runif(50), 10, 5)
#' C <- kr(X, Y)
#' 
#' @export
kr <- function(A, B, w, byrow = TRUE) {
    
    if (byrow) {
        if (nrow(A) != nrow(B)) 
            stop("Dimensions of the matrices do not match.")
        if (missing(w)) 
            w <- rep(1, nrow(A))
        if (nrow(A) != length(w)) 
            stop("Length of the weight does not match with the dimension of the\n             matrices.")
        cola <- ncol(A)
        colb <- ncol(B)
        colab <- cola * colb
        expr <- paste("rbind(", paste(rep("A", colb), collapse = ","), ")", sep = "")
        A <- eval(parse(text = expr))
        A <- matrix(c(A), nrow(B), ncol = colab)
        A <- w * A
        expr2 <- paste("cbind(", paste(rep("B", cola), collapse = ","), ")", sep = "")
        B <- eval(parse(text = expr2))
    } else {
        if (ncol(A) != ncol(B)) 
            stop("Dimensions of the matrices do not match.")
        if (missing(w)) 
            w <- rep(1, ncol(A))
        if (ncol(A) != length(w)) 
            stop("Length of the weight does not match with the dimension of the\n             matrices.")
        rowa <- nrow(A)
        rowb <- nrow(B)
        rowab <- rowa * rowb
        A <- matrix(rep(A, each = rowb), rowab, )
        A <- A * matrix(w, nrow(A), ncol(A), byrow = TRUE)
        expr <- paste("rbind(", paste(rep("B", rowa), collapse = ","), ")", sep = "")
        B <- eval(parse(text = expr))
    }
    
    return(A * B)
} 
