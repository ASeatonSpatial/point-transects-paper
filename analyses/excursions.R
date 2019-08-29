#### Try excursions package ####
library(INLA)
library(inlabru)
library(excursions)

set.seed(1525)

#### Fitted model  ####
file_name = "akepa_2002_fitted_model.RDS"
fit = readRDS(file_name)

study_area = readRDS("data/study_area_extended_no_crs.RDS")
samplers = readRDS("data/samplers_extended_no_crs.RDS")
mesh = readRDS("data/mesh_extended_no_crs.RDS")

#### monte-carlo samples of intensity ####

# Does excursions work with an inlabru object?

# There must be a way to get the actual values from
# predict with n.samples > 1 ?

n.mc = 500     # lower this when fiddling with things
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
u = 1
alpha = 0.05    # Prob(intensity > u) >= 1 - alpha for all s in E

ex = excursions.mc(Xhec, u = u, type = ">", alpha = alpha)

# plot
Fpxl = pxl
Fpxl$F = ex$F

# png(filename = "../figures/excursion_function.png",
#     width = 4.5, height = 5.5, units = "cm", res = 1000)
ggplot() + 
  gg(Fpxl) +
  scale_fill_viridis_c() +
  ggtitle("Excursion function for > 1 bird per hectare")

ggsave(filename = "../figures/excursion_function.png",
       width = 5, height = 6, units = "in")

Epxl = pxl
Epxl$E = ex$E
ggplot() +
  gg(Epxl) +
  ggtitle("Excursion set for > 1 bird per hectare,\nalpha = 0.05")

ggsave(filename = "../figures/excursion_set.png",
       width = 5, height = 6, units = "in")

# contour map - how to visualise this?  ?
cm = contourmap.mc(Xhec, levels = c(0.1, 0.5, 0.8, 1.5, 2, 4), alpha = alpha)
str(cm)

Gpxl = pxl
Gpxl$G = cm$G
ggplot() +
  gg(Gpxl)

cm = contourmap.mc(Xhec, levels = c(1), alpha = alpha)
Gpxl = pxl
Gpxl$G = cm$G
ggplot() +
  gg(Gpxl)
