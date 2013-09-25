# Kronecker product
# 
# Author: Yiwen Zhang
###############################################################################

kr <- function(A, B){
	if(ncol(A)!=ncol(B)) stop("Error: dimensions of the matrixes do not match")
	n <- ncol(A)
	kr <- sapply(1:n, function(i, A, B) return(A[,i]%x%B[,i]), A, B)
	return(kr)
}


