#' Plot the posterior of abundance for inlabru point process model
#' Does some rough things to get 
#' 
#' @param model  fitted inlabru model
#' @param predpts the thing that ipoints() returns, or something similar, must have weight column
#' @param ... other args passed on to predict()
#' @param n.sd how many standard deviations above and below the mean to consider (for domain of N)
#' @param n.grid roughly the number of discrete N values to consider

AS.Nposterior <- function(model, predpts, n.sd = 3, n.grid = 100, ...){
  
  # First get posterior for 
  Lambda <- predict(model, predpts, ~ sum(weight * exp(grf + Intercept))) 
  Lambda
  
  Nlow <- round(Lambda$mean - n.sd*Lambda$sd)
  Nhigh <- round(Lambda$mean + n.sd*Lambda$sd)
  by <- round((Nhigh - Nlow)/n.grid)   # should give roughly n.grid 
  
  Nrange <- seq(Nlow, Nhigh, by = by)
  Npost = predict(model, predpts, 
                 ~ data.frame(N = Nrange, 
                              dpois(Nrange, lambda = sum(weight * exp(grf + Intercept)))),
                 ...)
  
  return(list(Npost = Npost, Nexp = Lambda))
}



