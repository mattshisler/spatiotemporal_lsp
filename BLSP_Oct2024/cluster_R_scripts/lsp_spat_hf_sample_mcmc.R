args <- as.numeric(commandArgs(trailingOnly = T))
# rho1, rho2, pick-up-from-iteration-#, n.iters-to-sample, thinning 

# print(args)
print("Running mcmc sampler with:")
#print(paste0("rho equal to ", args[1], " ", args[2]))
print(paste0("picking up from iteration ", args[3]))
print(paste0("will run until iteration ", args[4], " iterations"))
print(paste0("thinning factor of ", args[5]))
print(paste0("number of cores ", args[6]))
print(paste0("iter print frequency ", args[7]))
print(paste0("num samples per file ", args[8]))

# load packages
library(spam)
library(lubridate)
library(minpack.lm)
library(MCMCpack)
library(collapse)
library(tidyr)
library(future.apply)

# load functions
source("./rfuncs/make_graph_laplace.R")
source("./rfuncs/nls_fit.R")
source("./rfuncs/double_logis7.R")

# load prepped workspace for hf application
load("./workspace/streamlined_workspace_753_img.RData")


rho1.list <- c(80,80,80,80,
               90,90,90,90,
               95,95,95,95,
               99,99,99,99)

rho2.list <- c(80,90,95,99,
               80,90,95,99,
               80,90,95,99,
               80,90,95,99)

rho.text <- c(rho1.list[args[1]], rho2.list[args[1]])

print(paste0("rho equal to ", rho.text[1], " ", rho.text[2]))

# parameters to vary for different runs.
iter.pickup <- args[3]
mcmc.pickup.file <- paste0("./burnin/mcmc_samples_rho1_", rho.text[1], "_rho2_", rho.text[2], "_iter_", iter.pickup, ".RDS")
mcmc.pickup <- readRDS(mcmc.pickup.file)
rho <- rho.text/100

#rho <- mcmc.pickup$rho
iter.pickup <- mcmc.pickup$iter

end.iter <- args[4]
n.iters <- end.iter - iter.pickup + 1
burn <- 0
thin <- args[5]      # only save every nth observation
n.workers <- args[6] # number of workers
print.freq <- args[7]
num.samples.per.object <- args[8]


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

# pickup the mcmc from where we left off
Omega.mcmc      <- mcmc.pickup$Omega.mcmc
Omega.inv.mcmc  <- solve(Omega.mcmc)
sigma2.mcmc     <- mcmc.pickup$sigma2.mcmc
Lambda.inv.mcmc <- mcmc.pickup$Lambda.inv.mcmc
Lambda.mcmc     <- lapply(Lambda.inv.mcmc, solve)
B.mcmc          <- mcmc.pickup$B.mcmc
Theta.mcmc      <- mcmc.pickup$Theta.mcmc
eta.mcmc        <- mcmc.pickup$eta.mcmc


# Some more precomputed terms
ZtZ <- t(Z)%*%Z
w.s.plus <- diag(L2)
sum.w.s.plus <- sum(w.s.plus)
sum.Bz <- list()
for (i in 1:m){
  sum.Bz[[i]] <- Reduce("+", Map("*", Z[i,], B.mcmc))  
}

# storage
keep.sigma2 <- rep(NA, num.samples.per.object)
keep.Theta <- array(NA, dim = c(p, m, n, num.samples.per.object))
# keep.B     <- list(array(NA, dim = c(p, n, num.samples.per.object)), array(NA, dim = c(p, n, num.samples.per.object)))
keep.Omega <- array(NA, dim = c(p, p, num.samples.per.object))
keep.B <- list()
keep.Lambda <- list()
keep.eta <- list()
for(j in 1:q){
 keep.B[[j]]      <- array(NA, dim = c(p, n, num.samples.per.object))
 keep.Lambda[[j]] <- array(NA, dim = c(p, p, num.samples.per.object))
 keep.eta[[j]]    <- matrix(0, nrow = p, ncol = num.samples.per.object)
}
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

psi <- p + 3 + 1e-6
G  <- spam_diag(p)*1
Sigma.inv <- diag(p)*1e-7
sum.Bz.mat <- t(do.call(rbind, lapply(sum.Bz, function(x) c(x))))
a <- 0.01
b <- 0.01

# function for parallelizing over years.
sample_Theta <- function(Q.theta,
                         sigma2.mcmc,
                         sumXtX,
                         sum.Bz,
                         sumRtX,
                         Q.theta.prime.chol.struct,
                         p, n){
  
  Q.theta.prime <- Q.theta + (1/sigma2.mcmc)*sumXtX
  b.theta.prime <- Q.theta%*%c(sum.Bz) + (1/sigma2.mcmc)*sumRtX
  chol.Q.theta.prime <- spam::update.spam.chol.NgPeyton(Q.theta.prime.chol.struct, Q.theta.prime)
  mu <- drop(solve.spam(chol.Q.theta.prime, b.theta.prime))
  nu <- backsolve(chol.Q.theta.prime, array(rnorm(n*p), c(n*p, 1)), k = n*p)
  return(matrix((nu + mu), nrow = p, ncol = n))
}


for(iter in (iter.pickup+1):end.iter){
  
  plan(multicore, workers = n.workers)
  Theta.mcmc <- future_lapply(1:m, function(x) sample_Theta(Q.theta,
                                                      sigma2.mcmc,
                                                      sumXtX[[x]],
                                                      sum.Bz[[x]],
                                                      sumRtX[[x]],
                                                      Q.theta.prime.chol.struct,
                                                      p, n), future.seed = TRUE)
  
  Theta.mcmc <- aperm(simplify2array(Theta.mcmc), c(1,3,2))
  if(iter%%thin == 0 & iter > burn) keep.Theta[ , , , store.idx] <- Theta.mcmc
  
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
    V.eta.inv <- solve(((1-rho2)^2)*sum.w.s.plus*Lambda.inv.mcmc[[j]] + Sigma.inv)
    M.eta <- ((1-rho2)^2)*Lambda.inv.mcmc[[j]]%*%rowSums(w.s.plus*B.mcmc[[j]])
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
  
  # Draw missing values - update RtX
  for(im in 1:length(doy.miss)){
    Theta.miss <- as.matrix(Theta.mcmc[ , year.miss[im], missing.idx.split[[im]]])
    R[missing.idx.split[[im]], image.miss.idx[im]] <- rnorm(n.missing.per.im[im], mean = colSums(Theta.miss*Xs.miss[[im]]), sd = sqrt(sigma2.mcmc))
  }
  
  RtX <- lapply(1:n.image, function(d) c(TRA(Xs[[d]], R[,d], "*")))
  RtX <- do.call("cbind", RtX)
  sumRtX <- unname(lapply(d.idx, function(x) matrix(rowSums(RtX[,x]), nrow = n*p, ncol = 1)))
  
  # Update sigma2
  XTheta <- apply(X.array*Theta.mcmc[, expand.theta, ], MARGIN = c(2,3), sum)
  SSE2   <- sum((R - t(XTheta))^2)
  sigma2.mcmc <- 1/rgamma(1, a + total.days*n/2, b + SSE2/2)
  if(iter%%thin == 0 & iter > burn) keep.sigma2[store.idx] <- sigma2.mcmc
  
  
  if(iter%%thin == 0 & iter > burn) store.idx <- store.idx + 1
  
  if(iter%%print.freq == 0 & iter > burn) print(paste0(c(rho.text[1], rho.text[2], "mcmc-iter", iter), collapse = "-"))
  
  if(store.idx > num.samples.per.object){
    iter.start <- iter - (num.samples.per.object*thin) + 1
    mcmc.samples <- list(iter.start = iter.start, 
                         iter.end = iter,
                         thin    = thin,
                         burn    = burn,
                         keep.Theta  = keep.Theta,
                         keep.B      = keep.B,
                         keep.Omega  = keep.Omega,
                         keep.Lambda = keep.Lambda,
                         keep.eta    = keep.eta,
                         keep.sigma2 = keep.sigma2)
    saveRDS(mcmc.samples, paste0("./samples/rho1_", rho.text[1], "_rho2_", rho.text[2], "/mcmc_samples_rho1_", rho.text[1] ,"_rho2_", rho.text[2], "_iters_", iter.start, "_", iter, ".RDS"))
    store.idx <- 1
  }
  
}
     #save the state of the chain.
 mcmc.state <- list(iter    = iter,
                    thin    = thin,
                    burn    = burn,
	                  rho     = rho,
                    Theta.mcmc  = Theta.mcmc,
                    B.mcmc      = B.mcmc,
                    Omega.mcmc  = Omega.mcmc,
                    Lambda.inv.mcmc = Lambda.inv.mcmc,
                    eta.mcmc    = eta.mcmc,
                    sigma2.mcmc = sigma2.mcmc)
 saveRDS(mcmc.state, paste0("./burnin/mcmc_samples_rho1_", args[1] ,"_rho2_", args[2], "_iter_", iter,".RDS"))

#mcmc.samples <- list(n.iters = n.iters,
#                     thin    = thin,
#                     burn    = burn,
#                     keep.Theta  = keep.Theta,
#                     keep.B      = keep.B,
#                     keep.Omega  = keep.Omega,
#                     keep.Lambda = keep.Lambda,
#                     keep.eta    = keep.eta,
#                     keep.sigma2 = keep.sigma2)

#saveRDS(mcmc.samples, paste0("./mcmc_samples_rho1_", 100*rho[1] ,"_rho2_", 100*rho[2], "_future.RDS"))

