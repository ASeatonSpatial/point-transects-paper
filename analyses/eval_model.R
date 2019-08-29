# Evaluate fitted model
library(INLA)
library(inlabru)

# set ggplot theme
theme_set(theme_minimal())

set.seed(9701071)

#### WARNING - as it stands this script will overwrite figures for the paper

#### Fitted model  ####
file_name = "akepa_2002_fitted_model.RDS"
fit = readRDS(file_name)

study_area = readRDS("data/study_area_extended_no_crs.RDS")
samplers = readRDS("data/samplers_extended_no_crs.RDS")
mesh = readRDS("data/mesh_extended_no_crs.RDS")

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
distdf <- data.frame(distance = seq(1,W,length=100))

hnpred <- predict(fit, distdf, ~hn(distance,lsig))
ggplot() + 
  gg(hnpred) +
  ggtitle("Posterior half-normal detection function")

ggsave(filename = "../figures/halfnormal.png",
       width = 5, height = 5, units = "in")

# half-normal polar coords?  Is that what this is?  
# distdf <- data.frame(distance = seq(1,W,length=100))
# dfun <- predict(fit, distdf, ~ dsamp(distance,lsig))
# ggplot() + 
#   gg(dfun)

# intensity
pxl = pixels(mesh, nx = 100, ny = 100, mask = study_area)
pr.int <- predict(fit, pxl, ~ exp(grf + Intercept))

p1 = ggplot() + 
  gg(study_area) +
  gg(pr.int) +
  scale_fill_viridis_c() +
  ggtitle("Posterior mean intensity") +  
  xlab("Easting") +
  ylab("Northing")

p2 = ggplot() +
  gg(study_area)  +
  gg(pr.int["cv"]) +
  gg(samplers, colour = "red") +
  scale_fill_viridis_c() +
  ggtitle("Posterior CV of intensity") +  
  xlab("Easting") +
  ylab("Northing")

png(filename = "../figures/intensity_mean_cv.png",
    width = 10, height = 6, units = "in", res = 100)
multiplot(p1, p2, cols = 2)
dev.off()

# lower and upper quantiles

v = c(pr.int["mean"]@data,
      pr.int["q0.025"]@data,
      pr.int["q0.975"]@data)

viridisscale =   scale_fill_viridis_c(limits = range(v))

p1 = ggplot() + 
  gg(study_area) + 
  gg(pr.int) +
  viridisscale +
  xlab("Easting") +
  ylab("Northing")

p2 = ggplot() +
  gg(study_area) +
  gg(pr.int["q0.025"]) +
  viridisscale +
  xlab("Easting") +
  ylab("Northing")

p3 = ggplot() +
  gg(study_area) +
  gg(pr.int["q0.975"]) +
  viridisscale +
  xlab("Easting") +
  ylab("Northing")

layout_matrix <- matrix(c(4, 1, 1, 4, 2, 2, 3, 3), nrow = 2, byrow = TRUE)  # triangle arrangement

png(filename = "../figures/intensity_quantiles.png",
    width = 10, height = 10, units = "in", res = 100)
multiplot(p1, p2, p3, layout = layout_matrix, byrow = TRUE, ncol = 2)
dev.off()

# show that these maps are misleading:
str(pr.int)
cell_area = as.numeric(pr.int@grid@cellsize["x"]*pr.int@grid@cellsize["y"])
sum(pr.int@data["mean"]*cell_area)
sum(pr.int@data["q0.025"]*cell_area)
sum(pr.int@data["q0.975"]*cell_area)

# plot "realized intensity" a few times
draw1 = predict(fit, pxl, ~ exp(grf + Intercept), n.samples = 1)
draw2 = predict(fit, pxl, ~ exp(grf + Intercept), n.samples = 1)
draw3 = predict(fit, pxl, ~ exp(grf + Intercept), n.samples = 1)

draw_all = c(draw1["mean"]@data,
      draw2["mean"]@data,
      draw3["mean"]@data)

viridisscale =   scale_fill_viridis_c(limits = range(draw_all))

p1 = ggplot() + 
  gg(study_area) + 
  gg(draw1["mean"]) +
  viridisscale +
  xlab("Easting") +
  ylab("Northing")

p2 = ggplot() + 
  gg(study_area) + 
  gg(draw2["mean"]) +
  viridisscale +
  xlab("Easting") +
  ylab("Northing")

p3 = ggplot() + 
  gg(study_area) + 
  gg(draw3["mean"]) +
  viridisscale +
  xlab("Easting") +
  ylab("Northing")

png(filename = "../figures/intensity_realized.png",
    width = 15, height = 6, units = "in", res = 100)
multiplot(p1,p2,p3, cols = 3) 
dev.off()

# GMRF params
spde.range = spde.posterior(fit, "grf", what = "range")
spde.logvar = spde.posterior(fit, "grf", what = "log.variance")

plot(spde.range)
ggsave(filename = "../figures/range_posterior.png",
       width = 5, height = 5, units = "in")

plot(spde.logvar)
ggsave(filename = "../figures/spde_logvar_posterior.png",
       width = 5, height = 5, units = "in")

# Abundance posterior
source("AS.Nposterior.R") 
predpts = ipoints(study_area, mesh)
abnc = AS.Nposterior(fit, predpts, n.sd = 5, n.grid = 250, n.samples = 20000)   # n.samples = 1000 default

# Expected Abundance posterior
abnc$Nexp

# Realized Abundance posterior
# plot(abnc$Npost)   # assymmetrical about mean?

# ggplot() + 
#   gg(abnc$Npost) +
#   ggtitle(paste("Expected Abundance:  ", abnc$Nexp$mean, "\n SD: ", abnc$Nexp$sd))

Npost = abnc$Npost

ggplot() +
  geom_line(data = Npost, aes(x = N, y = mean)) +
  geom_vline(aes(xintercept = abnc$Nexp$mean), colour = "blue") +
  ggtitle("Posterior of realized abundance")

ggsave("../figures/realized_abundance_posterior.png",
       width = 5, height = 5, units = "in")

Npost$Nexp = dpois(Npost$N, lambda = abnc$Nexp$mean)
ggplot() +
  geom_line(data = Npost, aes(x = N, y = mean)) +
  geom_line(data = Npost, aes(x = N, y = Nexp), linetype = "dashed")

ggsave("../figures/realized_abundance_vs_exp.png",
       width = 5, height = 5, units = "in")
