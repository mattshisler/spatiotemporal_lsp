args <- commandArgs(trailingOnly = T)


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

sample.files <- list.files(paste0("./samples/rho1_", rho.text[1], "_rho2_", rho.text[2], "_samples"), full.names = T)
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
n.iters <- length(keep.sigma2)

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

# Compute Diagnostics

# Theta Diagnostics
Theta.diag <- array(NA, c(p,m,n))
for(k in 1:p){
  for(i in 1:m){
    Theta.diag[k,i,] <- apply(keep.Theta[k,i,,], MARGIN = 1, function(x) coda::geweke.diag(x)$z)
  }
}

# Beta Diagnostics
Beta.diag <- replicate(q, matrix(NA, nrow = p, ncol = n), simplify = F)
for(j in 1:q){
  for(k in 1:p){
    Beta.diag[[j]][k,] <- apply(keep.beta[[j]][k,,], MARGIN = 1, function(x) coda::geweke.diag(x)$z)
  }
}

mcmc.diag <- list(n.iters    = n.iters,
                  Theta.diag = Theta.diag,
                  Beta.diag  = Beta.diag)

saveRDS(mcmc.diag, paste0("./diagnostics/mcmc_samples_rho1_", rho.text[1] ,"_rho2_", rho.text[2], "_diagnostics.RDS"))


