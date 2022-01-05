# compare two different fitted models
# to investigate effects of different versions of
# INLA

library(INLA)
library(inlabru)
library(ggplot2)
library(rgeos)
library(sf)
library(patchwork)

theme_set(theme_classic())

data_path = here::here("analyses", "data")
study_area = readRDS(here::here(data_path, "study_area_extended_no_crs.RDS"))
samplers = readRDS(here::here(data_path, "samplers_extended_no_crs.RDS"))
mesh = readRDS(here::here(data_path, "mesh_extended_no_crs.RDS"))

#### Fitted model  ####
model_path_new = here::here("analyses", "fitted_model.RDS")
model_path_old = here::here("holding_bay", "akepa_2002_fitted_model_new_inlabru.RDS")
fit_new = readRDS(model_path_new)
fit_old = readRDS(model_path_old)

summary(fit_new)
summary(fit_old)

# Prediction locations
pxl = pixels(mesh, nx = 300, ny = 300, mask = study_area)

# map units are per km^2 = 1,000,000 m^2
# I want intensity per hectare = 10,000 m^2
# so divide by 100
pr.int_new <- predict(fit_new, pxl, ~ exp(grf + Intercept)/100)
pr.int_old <- predict(fit_old, pxl, ~ exp(grf + Intercept)/100)

# mean 
vals = c(pr.int_new$mean, pr.int_old$mean)

# things for colour scale 
lower = min(vals)
upper = max(vals)
mean_breaks = c(lower, upper)
mean_labels = c(signif(lower, digits = 3),
                signif(upper, digits = 3))

p_new = ggplot() +
  gg(study_area) +
  gg(pr.int_new["mean"]) +
  scale_fill_viridis_c(breaks = mean_breaks,
                       labels = mean_labels) +
  xlab("Easting") +
  ylab("Northing") +
  ggtitle("New inlabru") + 
  coord_equal() +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(fill = guide_colourbar(title.vjust = 0.95,
                                title.theme = element_text(size = 18),
                                label.theme = element_text(size = 16)))

p_old = ggplot() +
  gg(study_area) +
  gg(pr.int_old["mean"]) +
  scale_fill_viridis_c(breaks = mean_breaks,
                       labels = mean_labels) +
  xlab("Easting") +
  ylab("Northing") +
  ggtitle("Old inlabru") + 
  coord_equal() +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(fill = guide_colourbar(title.vjust = 0.95,
                                title.theme = element_text(size = 18),
                                label.theme = element_text(size = 16)))

p_new + p_old

# sd
vals = c(pr.int_new$sd, pr.int_old$sd)

# things for colour scale 
lower = min(vals)
upper = max(vals)
mean_breaks = c(lower, upper)
mean_labels = c(signif(lower, digits = 3),
                signif(upper, digits = 3))

p_new = ggplot() +
  gg(study_area) +
  gg(pr.int_new["sd"]) +
  scale_fill_viridis_c(breaks = mean_breaks,
                       labels = mean_labels) +
  xlab("Easting") +
  ylab("Northing") +
  ggtitle("New inlabru") + 
  coord_equal() +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(fill = guide_colourbar(title.vjust = 0.95,
                                title.theme = element_text(size = 18),
                                label.theme = element_text(size = 16)))

p_old = ggplot() +
  gg(study_area) +
  gg(pr.int_old["sd"]) +
  scale_fill_viridis_c(breaks = mean_breaks,
                       labels = mean_labels) +
  xlab("Easting") +
  ylab("Northing") +
  ggtitle("Old inlabru") + 
  coord_equal() +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(fill = guide_colourbar(title.vjust = 0.95,
                                title.theme = element_text(size = 18),
                                label.theme = element_text(size = 16)))

p_new + p_old
