# Log of the Half-normal detection function
log_hn <- function(distance, sig) {
  -0.5 * (distance / sig)^2
}

# The Half-normal detection function
hn <- function(distance, sig) {
  exp(log_hn(distance, sig))
}

# Log of the Hazard Rate detection function
log_hr <- function(distance, sig, gam) {
  log1p(-exp(-(distance / sig)^(-gam)))
}

# The Hazard Rate detection function
hr <- function(distance, sig, gam) {
  1 - exp(-(distance / sig)^(-gam))
}
