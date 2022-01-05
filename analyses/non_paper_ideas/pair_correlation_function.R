# Try to get the pair correlation function

library(INLA)
library(inlabru)

# set ggplot theme
theme_set(theme_minimal())

fit = readRDS("akepa_2002_fitted_model.RDS")
study_area = readRDS("data/study_area_extended_no_crs.RDS")
samplers = readRDS("data/samplers_extended_no_crs.RDS")
mesh = readRDS("data/mesh_extended_no_crs.RDS")

summary(fit)

# inla.matern.cov needs kappa
rho.mean = fit$summary.hyperpar[1, "mean"]
kappa.mean = sqrt(8) / rho.mean
sigma.mean = fit$summary.hyperpar[2, "mean"]
var.mean = sigma.mean^2

h = seq(0.01, 20000, length.out = 250)

llam_cov = var.mean * inla.matern.cov(nu = 1,
                                      kappa = kappa.mean,
                                      x = h,
                                      d = 2)

df = data.frame(h = h,
                mat_cov = llam_cov)

ggplot() +
  geom_line(aes(x = h, y = mat_cov))

lam_cov = exp(

ggplot() +
  geom_line(aes(x = h, y = pcf))
