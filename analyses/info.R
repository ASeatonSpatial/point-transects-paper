# Evaluate fitted model
library(INLA)
library(inlabru)
library(sp)
library(raster)
library(patchwork)
library(rgeos)

#### Fitted model  ####
model_path = here::here("analyses", "akepa_2002_fitted_model_new_inlabru.RDS")
data_path = here::here("analyses", "data")
fit = readRDS(model_path)

#### Where to save figures ####
fig_path = here::here("figures")

study_area = readRDS(here::here(data_path, "study_area_extended_no_crs.RDS"))
samplers = readRDS(here::here(data_path, "samplers_extended_no_crs.RDS"))
mesh = readRDS(here::here(data_path, "mesh_extended_no_crs.RDS"))

# generate samples from joint posterior of Intercept and grf st dev
n.mc =  1000
ps = generate(fit,
              formula = ~ data.frame(Intercept_latent, Stdev_for_grf),
              n.samples = n.mc)
ps = do.call(rbind, ps)
head(ps)
range(ps$Intercept_latent)
range(ps$Stdev_for_grf)
hist(ps$Stdev_for_grf)

e1 = exp(2*(ps[,1] + ps[,2]))
e2 = exp(ps[,1] + ps[,2]/2)

A0 = 1 / n.mc * sum(e1) - (1 / n.mc * sum(e2))^2
A0

# regular posterior variance
pxl = pixels(mesh, mask = study_area)
pred = predict(fit, pxl, ~ exp(Intercept + grf))
A1 = pred$sd^2
max(A1)

# Info and plot it
Info = 1 - A1 / A0
pxl$Info_var = Info

g1 = ggplot() +
  gg(pxl, mapping = aes(fill = Info_var)) +
  gg(samplers, colour = "red", pch = 1) +
  coord_equal() +
  scale_fill_viridis_c(name = "I_var") +
  theme_minimal()

g1

# This looks like the inverse of the intensity map.  In the north
# sampling effort is lower and almost no birds detected.
# In the south sampling effort is higher and more birds detected

# The above is comparing variances, what about comparing sd's?

Info_sd = 1 - sqrt(A1) / sqrt(A0)
pxl$Info_sd = Info_sd
g2 = ggplot() +
  gg(pxl, mapping = aes(fill = Info_sd)) +
  gg(samplers, colour = "red", pch = 1) +
  coord_equal() +
  scale_fill_viridis_c(name = "I_sd") +
  theme_minimal()

g2

# Coefficient of variation version:
A1_cv = pred$cv
A0_cv = sqrt(A0) / mean(exp(ps[,1]))

Info_cv = 1 - A1_cv / A0_cv
pxl$Info_cv = Info_cv

g3 = ggplot() +
  gg(pxl, mapping = aes(fill = Info_cv)) +
  gg(samplers, colour = "red", pch = 1) +
  coord_equal() +
  scale_fill_viridis_c(name = "I_cv") +
  theme_minimal()

g3

g1 + g2 + g3 + plot_annotation(tag_levels = "A")
ggsave(filename = here::here(fig_path, "info.png"),
       width = 10, height = 5, dpi = 150)

# Regular sd plot
pxl$sd = pred$sd
g4 = ggplot() +
  gg(pxl, mapping = aes(fill = sd)) +
  gg(samplers, colour = "red", pch = 1) +
  coord_equal() +
  scale_fill_viridis_c() +
  theme_minimal()


#  Where am I with this:  possibly include this in the thesis
#  Don't bother putting this in the draft paper for now.

