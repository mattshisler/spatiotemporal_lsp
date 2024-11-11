make_graph_laplace <- function(length, height, buff, rho, eps = getOption("spam.eps")){
  require(spam)
  
  dims <- c(length + buff, height + buff)
  
  if(rho==0){
    return(spam_diag(dims[1]*dims[2]))
  }
  
  x <- numeric(dims[1])
  x[1:2] <- c(1, -rho[1])
  y <- numeric(prod(dims))
  y[dims[1] + 1] <- -rho[1]
  
  L <- kronecker(diag.spam(dims[2]), toeplitz.spam(x, eps = eps)) + toeplitz.spam(y, eps = eps)
  diag(L) <- rowSums(L!=0)-1
  return(L)
}