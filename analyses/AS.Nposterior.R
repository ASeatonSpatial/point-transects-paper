#' Plot the posterior of abundance for inlabru point process model
#' Does some rough things to get 
#' 
#' @param model  fitted inlabru model
#' @param predpts the thing that ipoints() returns, or something similar, must have weight column
#' @param ... other args passed on to predict()
#' 

AS.Nposterior <- function(model, predpts, ...){
  
  # First get posterior for 
  Lambda <- predict(model, predpts, ~ sum(weight * exp(grf + Intercept))) 
  Lambda
  
  Nlow <- round(Lambda$mean - 4*Lambda$sd)
  Nhigh <- round(Lambda$mean + 4*Lambda$sd)
  by <- round((Nhigh - Nlow)/500)   # should give roughly 500 
  
  Nrange <- seq(Nlow, Nhigh, by = by)
  Npost = predict(model, predpts, 
                 ~ data.frame(N = Nrange, 
                              dpois(Nrange, lambda = sum(weight * exp(grf + Intercept)))))
  
  return(list(Npost = Npost, Nexp = Lambda))
}



