#### Try excursions package ####
library(INLA)
library(inlabru)
library(excursions)


#### Fitted model  ####
file_name = "akepa_2002_fitted_model.RDS"
fit = readRDS(file_name)

study_area = readRDS("data/study_area_no_crs.RDS")
samplers = readRDS("data/samplers_no_crs.RDS")
mesh = readRDS("data/mesh_no_crs.RDS")

#### monte-carlo samples of intensity ####

# Does excursions work with an inlabru object?

# There must be a way to get the actual values from
# predict with n.samples > 1 ?

n.mc = 500
pxl = pixels(mesh, nx = 100, ny = 100, mask = study_area)
X = matrix(NA, nrow = nrow(pxl@coords), ncol = n.mc)
for (i in 1:n.mc){
  int_draw <- predict(fit, pxl, ~ exp(grf + Intercept), n.samples = 1)
  X[,i] = int_draw@data$mean
}

#### Try excursions ####

# re-scale X.  Currently per m^2
# per km^2:
Xkm = X*1000000

# per hectare:
Xhec = X*10000

# per hectare threshold:
u = 0.5
alpha = 0.05    # Prob(intensity > u) >= 1 - alpha for all s in E

ex = excursions.mc(Xhec, u = u, type = ">", alpha = alpha)

# plot
Fpxl = pxl
Fpxl$F = ex$F
ggplot() + 
  gg(Fpxl) +
  scale_fill_viridis_c()

Epxl = pxl
Epxl$E = ex$E
ggplot() +
  gg(Epxl)

# contour map - how to visualise this?  ?
cm = contourmap.mc(Xhec, levels = c(0.1, 0.5, 0.8, 1.5, 2, 4), alpha = alpha)
str(cm)

Gpxl = pxl
Gpxl$G = cm$G
ggplot() +
  gg(Gpxl)
