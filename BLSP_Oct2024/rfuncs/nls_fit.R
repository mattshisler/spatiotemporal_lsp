nls_fit <- function(doy, vi, 
                           model_form, 
                           init_val = c(-1.25, -0.5, 100, 2.5, 265, 2.5, 5), 
                           lower_bound = c(-10, -3, 1, -10, 100, -10, 4), 
                           upper_bound = c(10, 2, 200, 10, 370, 10, 6)){
  require(minpack.lm)
  
  start_list <- as.list(init_val)
  
  names(start_list) <- paste0("theta", 1:length(start_list))
  
  avg_fit <- nlsLM(model_form,
                   data = list(y = vi, 
                               t = doy),
                   weights = rep(1, length(t)),
                   start = start_list,
                   lower = lower_bound,
                   upper = upper_bound,
                   control = list("maxiter" = 1000))
  
  return(avg_fit)
  
}