args <- commandArgs(trailingOnly = T)

library(collapse)
library(abind)

rho1.list <- c(80,80,80,80,
               90,90,90,90,
               95,95,95,95,
               99,99,99,99)

rho2.list <- c(80,90,95,99,
               80,90,95,99,
               80,90,95,99,
               80,90,95,99)

rho.text <- c(rho1.list[args[1]], rho2.list[args[1]])

sample.files <- list.files(paste0("./samples/rho1_", rho.text[1], "_rho2_", rho.text[2]), full.names = T)
mcmc.samples <- readRDS(sample.files[1])

keep.Theta  <- mcmc.samples$keep.Theta
keep.B      <- mcmc.samples$keep.B
keep.sigma2 <- mcmc.samples$keep.sigma2
keep.eta    <- mcmc.samples$keep.eta
keep.Lambda <- mcmc.samples$keep.Lambda
keep.Omega  <- mcmc.samples$keep.Omega

p <- dim(keep.Theta)[1]
m <- dim(keep.Theta)[2]
n <- dim(keep.Theta)[3]
q <- length(keep.B)

for(file in sample.files[-1]){
  temp <- readRDS(file)
  keep.Theta  <- abind(keep.Theta, temp$keep.Theta, along = 4)
  keep.sigma2 <- c(keep.sigma2, temp$keep.sigma2)
  keep.Omega  <- abind(keep.Omega, temp$keep.Omega, along = 3)
  for(j in 1:q){
    keep.B[[j]]      <- abind(keep.B[[j]],      temp$keep.B[[j]], along = 3)
    keep.eta[[j]]    <- cbind(keep.eta[[j]],    temp$keep.eta[[j]]) 
    keep.Lambda[[j]] <- abind(keep.Lambda[[j]], temp$keep.Lambda[[j]], along = 3)
  }
  rm(temp)
}

n.iters <- length(keep.sigma2)

# Compute summaries
Theta.mean   <- apply(keep.Theta, MARGIN = c(1,2,3), mean, na.rm = T)
Theta.sd     <- apply(keep.Theta, MARGIN = c(1,2,3), sd,   na.rm = T)
Theta.lower  <- apply(keep.Theta, MARGIN = c(1,2,3), fquantile, probs = 0.05, na.rm = T, names = F)
Theta.upper  <- apply(keep.Theta, MARGIN = c(1,2,3), fquantile, probs = 0.95, na.rm = T, names = F)

Theta <- list(Theta.mean  = Theta.mean, 
              Theta.sd    = Theta.sd, 
              Theta.lower = Theta.lower, 
              Theta.upper = Theta.upper)

sigma2.mean  <- mean(keep.sigma2, na.rm = T)
sigma2.sd    <- sd(keep.sigma2,   na.rm = T)
sigma2.lower <- fquantile(keep.sigma2, probs = 0.05, na.rm = T, names = F)
sigma2.upper <- fquantile(keep.sigma2, probs = 0.95, na.rm = T, names = F)

sigma2 <- list(sigma2.mean  = sigma2.mean,
               sigma2.sd    = sigma2.sd,
               sigma2.lower = sigma2.lower,
               sigma2.upper = sigma2.upper)

Omega.mean  <- apply(keep.Omega, MARGIN = c(1,2), mean, na.rm = T)
Omega.sd    <- apply(keep.Omega, MARGIN = c(1,2), sd, na.rm = T)
Omega.lower <- apply(keep.Omega, MARGIN = c(1,2), fquantile, probs = 0.05, na.rm = T, names = F)
Omega.upper <- apply(keep.Omega, MARGIN = c(1,2), fquantile, probs = 0.95, na.rm = T, names = F)

Omega <- list(Omega.mean  = Omega.mean,
              Omega.sd    = Omega.sd,
              Omega.lower = Omega.lower,
              Omega.upper = Omega.upper)

B.mean      <- B.sd      <- B.lower      <- B.upper      <- list()
eta.mean    <- eta.sd    <- eta.lower    <- eta.upper    <- list()
Lambda.mean <- Lambda.sd <- Lambda.lower <- Lambda.upper <- list()   

for(j in 1:q){
  B.mean[[j]]   <- apply(keep.B[[j]], MARGIN = c(1,2), mean, na.rm = T)
  B.sd[[j]]     <- apply(keep.B[[j]], MARGIN = c(1,2), sd,   na.rm = T)
  B.lower[[j]]  <- apply(keep.B[[j]], MARGIN = c(1,2), fquantile, probs = 0.05, na.rm = T, names = F)
  B.upper[[j]]  <- apply(keep.B[[j]], MARGIN = c(1,2), fquantile, probs = 0.95, na.rm = T, names = F)

  eta.mean[[j]]   <- rowMeans(keep.eta[[j]], na.rm = T)
  eta.sd[[j]]     <- apply(keep.eta[[j]], MARGIN = 1, sd,   na.rm = T)
  eta.lower[[j]]  <- apply(keep.eta[[j]], MARGIN = 1, fquantile, probs = 0.05, na.rm = T, names = F)
  eta.upper[[j]]  <- apply(keep.eta[[j]], MARGIN = 1, fquantile, probs = 0.95, na.rm = T, names = F)
  
  Lambda.mean[[j]]   <- apply(keep.Lambda[[j]], MARGIN = c(1,2), mean, na.rm = T)
  Lambda.sd[[j]]     <- apply(keep.Lambda[[j]], MARGIN = c(1,2), sd,   na.rm = T)
  Lambda.lower[[j]]  <- apply(keep.Lambda[[j]], MARGIN = c(1,2), fquantile, probs = 0.05, na.rm = T, names = F)
  Lambda.upper[[j]]  <- apply(keep.Lambda[[j]], MARGIN = c(1,2), fquantile, probs = 0.95, na.rm = T, names = F)
  
}

B <- list(B.mean  = B.mean,
          B.sd    = B.sd,
          B.lower = B.lower,
          B.upper = B.upper)

eta <- list(eta.mean  = eta.mean,
            eta.sd    = eta.sd,
            eta.lower = eta.lower,
            eta.upper = eta.upper)

Lambda <- list(Lambda.mean  = Lambda.mean,
               Lambda.sd    = Lambda.sd,
               Lambda.lower = Lambda.lower,
               Lambda.upper = Lambda.upper)

mcmc.summaries <- list(n.iters = n.iters,
                       Theta  = Theta,
                       Omega  = Omega,
                       sigma2 = sigma2,
                       B      = B,
                       eta    = eta,
                       Lambda = Lambda)

saveRDS(mcmc.summaries, paste0("./summaries/mcmc_summaries_rho1_", rho.text[1], "_rho2_", rho.text[1], ".RDS"))