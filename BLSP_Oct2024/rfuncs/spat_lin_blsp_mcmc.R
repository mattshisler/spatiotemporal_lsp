spat_lin_blsp_mcmc <- function(Theta.init = NULL,
                               B.init = NULL,
                               Omega.init = NULL,
                               Lambda.init = NULL,
                               sigma2.fixed = 0.05^2,
                               X.array,
                               R,
                               n.day.per.year,
                               sumRtX,
                               sumXtX,
                               Z,
                               burn = 1000,
                               n.iters = 5000,
                               thin = max(1, floor(n.iters*1e-3)),
                               rho = c(0.99, 0.99),
                               domain.dim,
                               missing.vals = NULL){
  
  p <- dim(Theta.init)[1] # dimension of the response
  m <- dim(Theta.init)[2] # number of years
  n <- dim(Theta.init)[3] # number of sites
  q <- ncol(Z)            # number of covariates
  total.days     <- sum(n.day.per.year)
  expand.theta   <- rep(1:m, times = n.day.per.year)
  
  length <- domain.dim[1]
  height <- domain.dim[2]
  rho1 <- rho[1]
  rho2 <- rho[2]
  L1  <- make_graph_laplace(length, height, buff = 0, rho = rho[1])
  L2  <- make_graph_laplace(length, height, buff = 0, rho = rho[2])
  
  
  # initial values
  Omega.mcmc <- Omega.init
  Omega.inv.mcmc <- solve(Omega.mcmc)

  sigma2.mcmc <- sigma2.fixed
  Lambda.mcmc <- Lambda.init
  Lambda.inv.mcmc <- lapply(Lambda.mcmc, solve)
  B.mcmc <- B.init
  Theta.mcmc <- Theta.init
  eta.mcmc <- lapply(seq_len(q), function(x) matrix(rowMeans(B.init[[x]]), nrow = p, ncol = n))
  
  rho2.mcmc <- rho[2]
  
  
  # Some more precomputed terms
  ZtZ <- t(Z)%*%Z
  w.s.plus <- diag(L2)
  sum.w.s.plus <- sum(w.s.plus)
  sum.Bz <- list()
  for (i in 1:m){
    sum.Bz[[i]] <- Reduce("+", Map("*", Z[i,], B.mcmc))  
  }
  
  # storage
  keep.sigma2 <- rep(NA, n.iters/thin)
  keep.Theta <- array(NA, dim = c(p, m, n, n.iters/thin))
  keep.B     <- list(array(NA, dim = c(p, n, n.iters/thin)), array(NA, dim = c(p, n, n.iters/thin)))
  keep.Omega <- array(NA, dim = c(p, p, n.iters/thin))
  keep.Lambda <- list()
  keep.eta <- list()
  for(j in 1:q){
    keep.Lambda[[j]] <- array(NA, dim = c(p, p, n.iters/thin))
    keep.eta[[j]]    <- matrix(0, nrow = p, ncol = n.iters/thin) 
  }
  
  # keep initial values
  keep.Theta[,,,1] <- Theta.mcmc
  for(j in 1:q){
    keep.B[[j]][,,1] <- B.mcmc[[j]] 
    keep.Lambda[[j]][,,1] <- Lambda.mcmc[[j]]
  }
  keep.Omega[,,1]  <- Omega.mcmc
  store.idx <- 1
  
  Q.theta <- kronecker.spam(L1, Omega.inv.mcmc)
  Q.beta <- list()
  Q.beta.prime <- list()
  for(j in 1:q){
    Q.beta[[j]]  <- kronecker.spam(L2, Lambda.inv.mcmc[[j]])
    Q.beta.prime[[j]] <- ZtZ[j,j]*Q.theta + kronecker.spam(L2, Lambda.inv.mcmc[[j]])
  }
  
  Q.theta.prime.chol.struct <- spam::chol(Q.theta + (1/sigma2.mcmc)*sumXtX[[1]])
  Q.beta.prime.chol.struct  <- spam::chol(Q.beta.prime[[1]])
  
  psi <- p + 3 + 0.000001
  G  <- spam_diag(p)*1
  Sigma.inv <- diag(p)*1e-7
  sum.Bz.mat <- t(do.call(rbind, lapply(sum.Bz, function(x) c(x))))
  a <- 0.01
  b <- 0.01
  
  for(iter in 2:n.iters){
    
    # Update Theta
    for(i in 1:m){
      Q.theta.prime <- Q.theta + (1/sigma2.mcmc)*sumXtX[[i]]
      b.theta.prime <- Q.theta%*%c(sum.Bz[[i]]) + (1/sigma2.mcmc)*sumRtX[[i]]
      chol.Q.theta.prime <- spam::update.spam.chol.NgPeyton(Q.theta.prime.chol.struct, Q.theta.prime)
      
      mu <- drop(solve.spam(chol.Q.theta.prime, b.theta.prime))
      nu <- backsolve(chol.Q.theta.prime, array(rnorm(n*p), c(n*p, 1)), k = n*p)
      Theta.mcmc[,i,] <- matrix((nu + mu), nrow = p, ncol = n)
      if(iter%%thin == 0 & iter > burn) keep.Theta[ , i, , store.idx] <- Theta.mcmc[,i,]
    }
    
    # Update B, Lambda, and eta
    for(j in 1:q){
      # update B
      phi <- list()
      for (i in 1:m){
        phi[[i]] <- Theta.mcmc[,i,] - sum.Bz[[i]] + B.mcmc[[j]]*Z[i,j]
      }
      sum.phitZ <- Reduce("+", Map("*", phi, Z[,j]))
      
      Q.beta.prime[[j]] <- ZtZ[j,j]*Q.theta + Q.beta[[j]]
      
      b.beta.prime <- Q.theta%*%c(sum.phitZ) + Q.beta[[j]]%*%(c(eta.mcmc[[j]]))
      chol.Q.beta.prime <- spam::update.spam.chol.NgPeyton(Q.beta.prime.chol.struct, Q.beta.prime[[j]])
      
      mu <- drop(solve.spam(chol.Q.beta.prime, b.beta.prime))
      nu <- backsolve(chol.Q.beta.prime, array(rnorm(n*p), c(n*p, 1)), k = n*p)
      B.mcmc[[j]] <- matrix((nu + mu), nrow = p, ncol = n)
      if(iter%%thin == 0 & iter > burn) keep.B[[j]][ , , store.idx] <- B.mcmc[[j]]
      
      for (i in 1:m){
        sum.Bz[[i]] <- Reduce("+", Map("*", Z[i,], B.mcmc))  
      }
      
      # update Lambda
      Resid.2 <- B.mcmc[[j]] - eta.mcmc[[j]]
      H.beta <- (Resid.2)%*%(L2%*%t(Resid.2))
      Lambda.inv.mcmc[[j]] <- MCMCpack::rwish(psi + n, solve(G + H.beta))
      if(iter%%thin == 0 & iter > burn) keep.Lambda[[j]][,,store.idx] <- solve(Lambda.inv.mcmc[[j]])
      Q.beta[[j]] <- kronecker.spam(L2, Lambda.inv.mcmc[[j]])
      
      # Update eta
      V.eta.inv <- solve(((1-rho2.mcmc)^2)*sum.w.s.plus*Lambda.inv.mcmc[[j]] + Sigma.inv)
      M.eta <- ((1-rho2.mcmc)^2)*Lambda.inv.mcmc[[j]]%*%rowSums(w.s.plus*B.mcmc[[j]])
      eta.mcmc[[j]] <-  matrix(spam::rmvnorm(1, mu = V.eta.inv%*%M.eta, Sigma = V.eta.inv), nrow = p, ncol = n)
      if(iter%%thin == 0 & iter > burn) keep.eta[[j]][ , store.idx] <- eta.mcmc[[j]][,1] 
    }
    

    # Update Omega
    Resid.1 <- lapply(seq_len(m), function(x) Theta.mcmc[,x,]-sum.Bz[[x]])
    H.Omega <- Reduce("+", lapply(seq_len(m), function(x) (Resid.1[[x]])%*%(L1%*%t(Resid.1[[x]]))))
    
    Omega.mcmc <- MCMCpack::riwish(psi + m*n, G + H.Omega)
    # keep.Omega[,,iter] <- Omega.mcmc
    Omega.inv.mcmc <- solve(Omega.mcmc)
    
    if(iter%%thin == 0 & iter > burn) keep.Omega[ , , store.idx] <- Omega.mcmc
    Q.theta <- kronecker.spam(L1, Omega.inv.mcmc)
    
    # Update sigma2
    XTheta <- apply(X.array*Theta.init[, expand.theta, ], MARGIN = c(2,3), sum)
    SSE2   <- sum((R - t(XTheta))^2)
    sigma2.mcmc <- 1/rgamma(1, a + total.days*n/2, b + SSE2/2)
    if(iter%%thin == 0 & iter > burn) keep.sigma2[store.idx] <- sigma2.mcmc
    
    if(iter%%thin == 0 & iter > burn) store.idx <- store.idx + 1
    
    
  }
  
  return(list(n.iters = n.iters,
              thin    = thin,
              burn    = burn,
              keep.Theta  = keep.Theta,
              keep.B      = keep.B,
              keep.Omega  = keep.Omega,
              keep.Lambda = keep.Lambda,
              keep.eta    = keep.eta,
              keep.sigma2 = keep.sigma2,
              RtX = RtX,
              XtX = XtX,
              Z = Z,
              domain.dim = domain.dim,
              Theta.init = Theta.init,
              Omega.init = Omega.init,
              Lambda.init = Lambda.init,
              B.init = B.init,      
              sigma2.fixed = sigma2.fixed,
              rho = rho))
  
}