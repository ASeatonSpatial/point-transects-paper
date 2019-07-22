# Evaluate fitted model
library(INLA)
library(inlabru)

#### Fitted model  ####
file_name = "akepa_2002_fitted_model.RDS"
fit = readRDS(file_name)

#### Specify detection functions  ####

# Log half-normal
log_hn = function(distance, lsig){ 
  -0.5*(distance/exp(lsig))^2
}

# Half-normal
hn <- function(distance, lsig) exp(log_hn(distance, lsig))

# Include r for switching to polar coords
dsamp = function(distance, lsig){ 
  log(distance) + log_hn(distance, lsig)
}


#### Model evaluation ####

summary(fit)

# half-normal

W = 60   # transect radius

hnpred <- predict(fit, distdf, ~hn(distance,lsig))
ggplot() + 
  gg(hnpred)

ggsave(filename = "../figures/halfnormal.png")

# half-normal polar coords
distdf <- data.frame(distance = seq(1,W,length=100))
dfun <- predict(fit, distdf, ~ dsamp(distance,lsig))
ggplot() + 
  gg(dfun)

# intensity
pxl = pixels(mesh, nx = 100, ny = 100, mask = study_area)
pr.int <- predict(fit, pxl, ~ exp(grf + Intercept))

ggplot() + 
  gg(study_area) + 
  gg(pr.int) +
  scale_fill_viridis_c()

ggsave(filename = "../figures/intensity_mean.png")

ggplot() +
  gg(study_area)  +
  gg(pr.int["cv"]) +
  gg(samplers, colour = "red") +
  scale_fill_viridis_c()

ggsave(filename = "../figures/intensity_cv.png")

# GMRF params
spde.range = spde.posterior(fit, "grf", what = "range")
spde.logvar = spde.posterior(fit, "grf", what = "log.variance")

plot(spde.range)
ggsave(filename = "../figures/range_posterior.png")

plot(spde.logvar)
ggsave(filename = "../figures/spde_logvar_posterior.png")

# Abundance posterior
source("AS.Nposterior.R") 
predpts = ipoints(study_area, mesh)
abnc = AS.Nposterior(fit, predpts)   # n.samples = 1000 default

# Expected Abundance
abnc$Nexp

# Abundance posterior
plot(abnc$Npost)

ggplot() + 
  gg(abnc$Npost) +
  ggtitle(paste("Expected Abundance:  ", abnc$Nexp$mean, "\n SD: ", abnc$Nexp$sd))

ggsave("../figures/abundance_posterior.png")