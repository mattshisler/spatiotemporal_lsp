args <- commandArgs(trailingOnly = T)

library(lubridate)
library(spam)

source("./rfuncs/double_logis7.R")

load("./prepped_workspace_for_hf_lsp.RData")

sample.files <- list.files(paste0("./rho1_", args[1], "_rho2_", args[2], "_samples"), full.names = T)
mcmc.samples <- readRDS(sample.files[1])

keep.Theta  <- mcmc.samples$keep.Theta[,,args[3],]

p <- dim(keep.Theta)[1]
m <- dim(keep.Theta)[2]
n <- dim(keep.Theta)[3]
q <- length( mcmc.samples$keep.B)

for(file in sample.files[-1]){
  temp <- readRDS(file)
  keep.Theta  <- abind(keep.Theta, temp$keep.Theta, along = 3)
  rm(temp)
}

date.vec.full <- seq(ymd("1984-01-01"),ymd("2020-12-31"),by='days')
doy.vec.full <- yday(date.vec.full)
mu.full <- double_logis7(doy.vec.full, point.agg.theta[, args[3]]) - X[[args[3]]][doy.vec.full,]%*%point.agg.theta[, args[3]]
year.vec.full <- year(date.vec.full) - 1983

mean.vi.full <- matrix(nrow = length(date.vec.full), ncol = dim(keep.Theta)[3])
for(iter in 1:dim(keep.Theta)[3]){
  mean.vi.full[,iter] <- rowSums(X[[args[3]]][doy.vec.full,]*t(pixel.samples[,year.vec.full,iter])) + mu.full
}

ci.vi.full <- apply(mean.vi.full, MARGIN = 1, quantile, probs = c(0.05, 0.95), names = F)
saveRDS(ci.vi.full, paste0("./ci_vi_full_", args[1], "_", args[2], ".RDS"))




