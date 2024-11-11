double_logis7 <- function(d, theta){
  # double logistic function.
  # theta1 is transformed using the logistic function.
  # theta2 is transformed using the 
  # theta4 is transformed using the 
  # theta6 is transformed using the 
  # This allows for all parameters to follow a gaussian distribution
  
  theta[1] <- plogis(theta[1])
  theta[2] <- exp(theta[2])
  theta[4] <- exp(theta[4])
  theta[6] <- exp(theta[6])
  theta[7] <- theta[7]/(1e4)
  
  d1 <- 1 + exp((theta[3] - d)/theta[4])
  
  d2 <- 1 + exp((theta[5] - d)/theta[6])
  
  out <- theta[1] + (theta[2] + (theta[3] + (theta[5] - theta[3])/2  - d)*theta[7] - theta[1])*(1/d1 - 1/d2)
  
  
  return(out)
}