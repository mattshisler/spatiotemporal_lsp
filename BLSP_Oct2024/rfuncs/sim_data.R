sim_data <- function(eta,
                     q = 2,
                     Omega.sd,
                     Omega.K = 0.1,
                     Lambda.sd,
                     Lambda.K = 0.1,
                     lambda = c(1,100),
                     B = NULL,
                     dim.domain,
                     start.date,
                     end.date,
                     return.period = 16,
                     missing.image.prob = 0,
                     sigma = 0.05,
                     rho = c(0.99, 0.8),
                     mean.func){
  
  require(lubridate)
  require(MCMCpack)
  # source("/Users/matthewshisler/Spatial BLSP - 20 May 2024/double_logis7.R")
  # source("/Users/matthewshisler/Spatial BLSP - 20 May 2024/make_graph_laplace.R")
  source("make_graph_laplace.R")
  
  length <- dim.domain[1]
  height <- dim.domain[2]
  n      <- length*height
  m      <- length(lubridate::year(start.date):lubridate::year(end.date))
  
  date <- seq(start.date, end.date, by = return.period)
  date <- date[rbinom(length(date), 1, missing.image.prob)==0]
  doy  <- split(date, year(date))
  doy  <- lapply(doy, yday)
  
  doy.list <- doy
  doy.vec <- unlist(doy)
  date.vec <- unlist(date)
  years <- unique(year(date))
  p <- length(eta[[1]])
  
  # Omega.cor <- cov2cor(MCMCpack::riwish((p+1) + Omega.K, (Omega.K)*diag(p)))
  
  ar1_cor <- function(n, rho) {
    exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                      (1:n - 1))
    rho^exponent
  }
  
  Omega.cor <- cov2cor(ar1_cor(p, Omega.K))
  Omega <- diag(Omega.sd)%*%Omega.cor%*%diag(Omega.sd)
  Omega.inv <- solve(Omega)
  
  L1          <- make_graph_laplace(length, height, buff = 0, rho = rho[1])
  chol.L1     <- chol(as.matrix(L1))
  L1.inv      <- solve(L1)
  chol.L1.inv <- chol(as.matrix(L1.inv))
  
  if (is.null(B)){
    lambda <- 1/c(1,100) 
    Lambda.cor <- list()
    # Lambda.cor[[1]] <- cov2cor(MCMCpack::riwish((p+1) + Lambda.K, (Lambda.K)*diag(p)))
    # Lambda.cor[[2]] <- cov2cor(MCMCpack::riwish((p+1) + Lambda.K, (Lambda.K)*diag(p)))
    Lambda.cor[[1]] <- cov2cor(ar1_cor(p, Lambda.K))
    Lambda.cor[[2]] <- cov2cor(ar1_cor(p, Lambda.K))
    Lambda <- list()
    Lambda[[1]] <- (lambda[1])*(diag(Lambda.sd)%*%Lambda.cor[[1]]%*%diag(Lambda.sd))
    Lambda[[2]] <- (lambda[2])*(diag(Lambda.sd)%*%Lambda.cor[[2]]%*%diag(Lambda.sd))
    
    L2          <- make_graph_laplace(length, height, buff = 0, rho = rho[2])
    chol.L2     <- chol(as.matrix(L2))
    L2.inv      <- solve(L2)
    chol.L2.inv <- chol(as.matrix(L2.inv))
  
    B <- list()
    B[[1]] <- matrix(eta[[1]], nrow = p, ncol = n) + t(chol(Lambda[[1]]))%*%matrix(rnorm(p*n, 0, 1), nrow = p, ncol = n)%*%chol.L2.inv
    B[[2]] <- matrix(eta[[2]], nrow = p, ncol = n) + t(chol(Lambda[[2]]))%*%matrix(rnorm(p*n, 0, 1), nrow = p, ncol = n)%*%chol.L2.inv
  }
  
  I.n <- spam_diag(n)
  
  # spat.domain <- expand.grid(x = 1:length, y = 1:height)
  # spat.domain$label <- 1:n

  Z <- matrix(c(rep(1, m), scale(seq(1:m), scale = F)), ncol = q)
  
  Theta <- array(dim = c(p, m, n))
  for(i in 1:m){
    Theta[,i,] <- B[[1]] + B[[2]]*Z[i,2] + t(chol(Omega))%*%matrix(rnorm(p*n, 0, 1), nrow = p, ncol = n)%*%chol.L1.inv
  }
  
  tVI <- lapply(seq_len(m),  function(i) lapply(doy[[i]], function(j)  apply(Theta[,i,], MARGIN = 2, function(x) {mean.func(j, x)})))
  mtVI <- unlist(lapply(tVI, function(x) lapply(x, function(y) matrix(y, nrow = length))), recursive = F)
  atVI <- simplify2array(mtVI)
  otVI <- atVI + array(rnorm(prod(dim(atVI)), 0, sd = sigma), dim = dim(atVI))
  
  spat.domain <- expand.grid(x = 1:length, y = 1:height)
  spat.domain$label <- 1:n
  
  return(list(domain.dim = c(length, height),
         spat_domain = expand.grid(x = 1:length, y = 1:height),
         otVI = otVI,
         atVI = atVI,
         n.years = m,
         date.vec = date,
         doy.vec = doy.vec,
         doy.list = doy,
         eta = eta,
         B = B,
         Theta = Theta,
         Omega = Omega,
         Lambda = Lambda,
         Z = Z,
         rho = rho))
  
}
