# Evaluate fitted model
library(INLA)
library(inlabru)
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

W = 58   # transect radius
distdf <- data.frame(distance = seq(1,W,length=100))

hnpred <- predict(fit, distdf, ~hn(distance,lsig))
ghn = ggplot() +
        gg(hnpred) +
        ylab("probability of detection\n") +
        xlab("\ndistance (m)") +
        theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 16))

# detectability per transect
df = ipoints(c(0, W), 50, name = "distance")
dpred = predict(fit, df, ~ 2/W^2*sum(weight*exp(dsamp(distance, lsig))))

# Matern plots ##

blah = spde.posterior(fit, "grf", what = "matern.correlation")
blah$x = blah$x / 1000
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

# ggsave(filename = "../figures/matern_posterior.png",
#        width = 5, height = 5, units = "in")


# plot
# ggplot() +
#   gg(study_area) +
#   gg(mesh) +
#   gg(samplers, colour = "red") +
#   scale_fill_viridis_c() +
#   xlab("Easting") +
#   ylab("Northing") +
#   coord_equal()
#
# plot(mesh, lwd = 0.2, asp =1, main = "")

# intensity plots
pxl = pixels(mesh, nx = 100, ny = 100, mask = study_area)
pr.int <- predict(fit, pxl, ~ exp(grf + Intercept))

# mean
lower = min(pr.int["mean"]$mean)
lower = signif(lower, digits = 3)
lower

upper = max(pr.int["mean"]$mean)
upper = signif(upper, digits = 3)
upper

mean_breaks = c(lower, upper)
p1 = ggplot() +
  gg(study_area) +
  gg(pr.int["mean"]) +
  scale_fill_viridis_c(breaks = mean_breaks, labels = mean_breaks) +
  # ggtitle("A") +
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
lower = signif(lower, digits = 3)
lower

upper = max(pr.int$cv)
upper = signif(upper, digits = 3)
upper

cv_breaks = c(lower, upper)

p2 = ggplot() +
  gg(study_area)  +
  gg(pr.int["cv"]) +
  scale_fill_viridis_c(breaks = cv_breaks, labels = cv_breaks) +
  # ggtitle("B") +
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
lower = signif(lower, digits = 3)
lower = 2.23e-06

upper = max(pr.int$sd)
upper = signif(upper, digits = 3)
upper

sd_breaks = c(lower, upper)
sd_labels = format(sd_breaks, scientific = TRUE)

p3 = ggplot() +
  gg(study_area)  +
  gg(pr.int["sd"]) +
  scale_fill_viridis_c(breaks = sd_breaks, labels = sd_labels) +
  # ggtitle("C") +
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

png(filename = here::here(fig_path, "intensity_mean.png"),
    width = 10, height = 6, units = "in", res = 100)
p1
dev.off()

png(filename = here::here(fig_path, "figures/akepa_transects.png"),
    width = 10, height = 6, units = "in", res = 100)
p3
dev.off()

# lower and upper quantiles

v = c(pr.int$q0.025, pr.int$q0.975)

# scale legend breaks
lower = min(v)
lower = signif(lower, digits = 3)
lower = 8.32e-08

upper = max(v)
upper = signif(upper, digits = 3)
upper = 0.000943

q_breaks = c(lower, upper)
q_labels = format(q_breaks, scientific = TRUE)

p2 = ggplot() +
  gg(study_area) +
  gg(pr.int["q0.025"]) +
  scale_fill_viridis_c(limits = range(v), breaks = q_breaks, labels = q_labels)+
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
  scale_fill_viridis_c(limits = range(v), breaks = q_breaks, labels = q_labels)+
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


