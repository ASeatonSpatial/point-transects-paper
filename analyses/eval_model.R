# Evaluate fitted model
library(INLA)
library(inlabru)
library(ggplot2)
library(scales)
library(cowplot)
library(rgeos)

# set ggplot theme
theme_set(theme_minimal())

set.seed(9701071)

#### WARNING - as it stands this script will overwrite figures for the paper

#### Fitted model  ####
model_path = here::here("analyses", "akepa_2002_fitted_model_new_inlabru.RDS")
fit = readRDS(model_path)

#### Other things ####
data_path = here::here("analyses", "data")
study_area = readRDS(here::here(data_path, "study_area_extended_no_crs.RDS"))
samplers = readRDS(here::here(data_path, "samplers_extended_no_crs.RDS"))
mesh = readRDS(here::here(data_path, "mesh_extended_no_crs.RDS"))

#### Where to save figures ####
fig_path = here::here("figures")

#### ggplot theme ####
theme_set(theme_classic())

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

W = 58/1000
Wm = 58    # transect radius in metres
distdf <- data.frame(distance = seq(.Machine$double.eps, Wm, length=100))

# adjust for change of units (km in model fit, metres here)
hnpred <- predict(fit, distdf, ~ hn(distance, lsig + log(1000)),
                  n.samples = 100)

ghn = ggplot() +
        gg(hnpred) +
        ylab("probability of detection\n") +
        xlab("\ndistance (m)") +
        theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 16))

# Matern plots ##

blah = spde.posterior(fit, "grf", what = "matern.correlation")
# blah$x = blah$x / 1000
gm = ggplot() +
      gg(blah) +
      ylab("correlation\n") +
      xlab("\ndistance (km)") +
      theme(axis.title = element_text(size = 16, vjust = 0.3),
            axis.text = element_text(size = 16))

gm

png(filename = here::here(fig_path, "detfn_and_matern.png"),
    width = 10, height = 5, units = "in", res = 100)
plot_grid(NULL, ghn, NULL, gm, NULL,
          labels = c("", "A", "", "B", ""),
          rel_widths = c(0.05, 1, 0.15, 1, 0.05),
          ncol = 5)
dev.off()

# intensity plots
pxl = pixels(mesh, nx = 300, ny = 300, mask = study_area)
nrow(pxl@coords)

# map units are per km^2 = 1,000,000 m^2
# I want per hectare = 10,000 m^2
# so divide by 100
pr.int <- predict(fit, pxl, ~ exp(grf + Intercept)/100)

# mean
lower = min(pr.int["mean"]$mean)
upper = max(pr.int["mean"]$mean)
mean_breaks = c(lower, upper)
mean_labels = c(signif(lower, digits = 3),
                signif(upper, digits = 3))

p1 = ggplot() +
  gg(study_area) +
  gg(pr.int["mean"]) +
  scale_fill_viridis_c(breaks = mean_breaks,
                       labels = mean_labels) +
  xlab("Easting") +
  ylab("Northing") +
  coord_equal() +
  theme_void() +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(fill = guide_colourbar(title.vjust = 0.95,
                                title.theme = element_text(size = 18),
                                label.theme = element_text(size = 16)))

p1

# cv
lower = min(pr.int$cv)
upper = max(pr.int$cv)
cv_breaks = c(lower, upper)
cv_labels = c(signif(lower, digits = 3),
              signif(upper, digits = 3))

p2 = ggplot() +
  gg(study_area)  +
  gg(pr.int["cv"]) +
  scale_fill_viridis_c(breaks = cv_breaks,
                       labels = cv_labels) +
  xlab("Easting") +
  ylab("Northing") +
  coord_equal() +
  theme_void() +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(fill = guide_colourbar(title.vjust = 0.95,
                                title.theme = element_text(size = 18),
                                label.theme = element_text(size = 16)))

p2

# sd
lower = min(pr.int$sd)
upper = max(pr.int$sd)
sd_breaks = c(lower, upper)
sd_labels = format(c(signif(lower, digits = 3),
                     signif(upper, digits = 3)),
                   scientific = TRUE)


p3 = ggplot() +
  gg(study_area)  +
  gg(pr.int["sd"]) +
  scale_fill_viridis_c(breaks = sd_breaks,
                       labels = sd_labels) +
  xlab("Easting") +
  ylab("Northing") +
  coord_equal() +
  theme_void() +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(fill = guide_colourbar(title.vjust = 0.95,
                                title.theme = element_text(size = 18),
                                label.theme = element_text(size = 16)))

p3

png(filename = here::here(fig_path, "intensity_mean_cv_sd.png"),
    width = 10, height = 5, units = "in", res = 100)
plot_grid(NULL, p1, NULL, p2, NULL, p3, NULL,
          labels = c("", "A", "", "B", "", "C", ""),
          rel_widths = c(0.1, 1, 0.4, 1, 0.4, 1, 0.1),
          ncol = 7)
dev.off()

# png(filename = "../figures/intensity_mean_cv.png",
#     width = 10, height = 6, units = "in", res = 100)
# multiplot(p1, p2, p3, cols = 3)
# dev.off()

# png(filename = here::here(fig_path, "intensity_mean.png"),
#     width = 10, height = 6, units = "in", res = 100)
# p1
# dev.off()

# png(filename = here::here(fig_path, "figures/akepa_transects.png"),
#     width = 10, height = 6, units = "in", res = 100)
# p3
# dev.off()

# lower and upper quantiles

v = c(pr.int$q0.025, pr.int$q0.975)

# scale legend breaks
lower = min(v)
upper = max(v)
q_breaks = c(lower, upper)
q_labels = sd_labels = format(c(signif(lower, digits = 3),
                                signif(upper, digits = 3)),
                              scientific = TRUE)

p2 = ggplot() +
  gg(study_area) +
  gg(pr.int["q0.025"]) +
  scale_fill_viridis_c(limits = q_breaks,
                       breaks = q_breaks,
                       labels = q_labels)+
  coord_equal() +
  xlab("Easting") +
  ylab("Northing") +
  theme_void() +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(fill = guide_colourbar(title = "",
                                title.vjust = 0.95,
                                title.theme = element_text(size = 18),
                                label.theme = element_text(size = 16)))

p2

p3 = ggplot() +
  gg(study_area) +
  gg(pr.int["q0.975"]) +
  scale_fill_viridis_c(limits = q_breaks,
                       breaks = q_breaks,
                       labels = q_labels)+
  coord_equal() +
  xlab("Easting") +
  ylab("Northing") +
  theme_void() +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(fill = guide_colourbar(title = "",
                                title.vjust = 0.95,
                                title.theme = element_text(size = 18),
                                label.theme = element_text(size = 16)))

p3

png(filename = here::here(fig_path, "intensity_quantiles.png"),
    width = 10, height = 5, units = "in", res = 100)
plot_grid(NULL, p2, NULL, p3, NULL,
          labels = c("", "A", "", "B", ""),
          rel_widths = c(0.05, 1, 0.15, 1, 0.05),
          ncol = 5)
dev.off()

# show that these maps are misleading:
str(pr.int)
cell_area = as.numeric(pr.int@grid@cellsize["x"]*pr.int@grid@cellsize["y"]*100)
sum(pr.int@data["mean"]*cell_area)
sum(pr.int@data["q0.025"]*cell_area)
sum(pr.int@data["q0.975"]*cell_area)

# plot "realized intensity" a few times
set.seed(1989)
draw1 = predict(fit, pxl, ~ exp(grf + Intercept), n.samples = 1)
draw2 = predict(fit, pxl, ~ exp(grf + Intercept), n.samples = 1)
draw3 = predict(fit, pxl, ~ exp(grf + Intercept), n.samples = 1)

draw_all = c(draw1["mean"]@data,
      draw2["mean"]@data,
      draw3["mean"]@data)

viridisscale = scale_fill_viridis_c(limits = range(draw_all))

p1 = ggplot() +
  gg(study_area) +
  gg(draw1["mean"]) +
  viridisscale +
  xlab("Easting") +
  ylab("Northing") +
  theme_void() +
  theme(legend.position = "none")

p2 = ggplot() +
  gg(study_area) +
  gg(draw2["mean"]) +
  viridisscale +
  xlab("Easting") +
  ylab("Northing") +
  theme_void() +
  theme(legend.position = "none")

p3 = ggplot() +
  gg(study_area) +
  gg(draw3["mean"]) +
  viridisscale +
  xlab("Easting") +
  ylab("Northing") +
  theme_void() +
  theme(legend.position = "none")

png(filename = here::here(fig_path, "intensity_realized.png"),
    width = 10, height = 5, units = "in", res = 100)
plot_grid(NULL, p1, NULL, p2, NULL, p3, NULL,
          labels = c("", "A", "", "B", "", "C", ""),
          rel_widths = c(0.1, 1, 0.4, 1, 0.4, 1, 0.1),
          ncol = 7)

dev.off()

# GMRF params
spde.range = spde.posterior(fit, "grf", what = "range")
spde.logvar = spde.posterior(fit, "grf", what = "log.variance")

plot(spde.range)
ggsave(filename = here::here(fig_path, "range_posterior.png"),
       width = 5, height = 5, units = "in")

plot(spde.logvar)
ggsave(filename = here::here(fig_path, "spde_logvar_posterior.png"),
       width = 5, height = 5, units = "in")


