# Produce summary figures of the posterior intensity field

library(INLA)
library(inlabru)
library(ggplot2)
library(scales)
library(cowplot)
library(sf)

set.seed(9701071)

#### Fitted model  ####
model_path = here::here("R",
                        "fitted_models",
                        "fitted_model.RDS")
fit = readRDS(model_path)

#### Other things ####
data_path = here::here("R", "data")
study_area = readRDS(here::here(data_path, "study_area.RDS"))
samplers = readRDS(here::here(data_path, "samplers.RDS"))
mesh = readRDS(here::here(data_path, "mesh.RDS"))

#### Where to save figures ####
fig_path = here::here("figures")

#### ggplot theme ####
# theme_set(theme_minimal())
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

#### Posterior half-normal plot ####

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
ghn

##### Posterior Matern plot ####
blah = spde.posterior(fit, "grf", what = "matern.correlation")
gmat = ggplot() +
  gg(blah) +
  ylab("correlation\n") +
  xlab("\ndistance (km)") +
  theme(axis.title = element_text(size = 16, vjust = 0.3),
        axis.text = element_text(size = 16))

gmat

png(filename = here::here(fig_path, "detfn_and_matern.png"),
    width = 10, height = 5, units = "in", res = 100)
plot_grid(NULL, ghn, NULL, gmat, NULL,
          labels = c("", "A", "", "B", ""),
          rel_widths = c(0.05, 1, 0.15, 1, 0.05),
          ncol = 5)
dev.off()

##### Posterior intensity plots #####

# Prediction locations
pxl = fm_pixels(mesh,
                dims = c(300, 300),
                mask = study_area)

# map units are per km^2 = 1,000,000 m^2
# I want intensity per hectare = 10,000 m^2
# so divide by 100
pr.int <- predict(fit, pxl, ~ exp(grf + Intercept)/100,
                  n.samples = 1000)

# add x and y columns for geom_tile
pr.int.coords = st_coordinates(pr.int)
pr.int$x = pr.int.coords[,1]
pr.int$y = pr.int.coords[,2]

# mean
lower = min(pr.int$mean)
upper = max(pr.int$mean)
mean_breaks = c(lower, upper)
mean_labels = c(signif(lower, digits = 3),
                signif(upper, digits = 3))

p1 = ggplot() +
  geom_sf(data = study_area) +
  geom_tile(data = pr.int,
            mappin = aes(x = x,
                         y = y,
                         fill = mean)) +
  scale_fill_viridis_c(breaks = mean_breaks,
                       labels = mean_labels) +
  xlab("Easting") +
  ylab("Northing") +
  theme_void() +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(fill = guide_colourbar(title.vjust = 0.95,
                                title.theme = element_text(size = 18),
                                label.theme = element_text(size = 16)))

p1

# cv
pr.int$cv = pr.int$sd / pr.int$mean
lower = min(pr.int$cv)
upper = max(pr.int$cv)
cv_breaks = c(lower, upper)
cv_labels = c(signif(lower, digits = 3),
              signif(upper, digits = 3))

p2 = ggplot() +
  geom_sf(data = study_area) +
  geom_tile(data = pr.int,
            mapping = aes(x = x,
                         y = y,
                         fill = cv)) +
  scale_fill_viridis_c(breaks = cv_breaks,
                       labels = cv_labels) +
  xlab("Easting") +
  ylab("Northing") +
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
  geom_sf(data = study_area) +
  geom_tile(data = pr.int,
            mapping = aes(x = x,
                         y = y,
                         fill = sd)) +
  scale_fill_viridis_c(breaks = sd_breaks,
                       labels = sd_labels) +
  xlab("Easting") +
  ylab("Northing") +
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


##### lower and upper quantiles ####

# scale legend breaks
v = c(pr.int$q0.025, pr.int$q0.975)
lower = min(v)
upper = max(v)
q_breaks = c(lower, upper)
q_labels = sd_labels = format(c(signif(lower, digits = 3),
                                signif(upper, digits = 3)),
                              scientific = TRUE)

p2 = ggplot() +
  geom_sf(data = study_area) +
  geom_tile(data = pr.int,
            mapping = aes(x = x,
                         y = y,
                         fill = q0.025)) +
  scale_fill_viridis_c(limits = q_breaks,
                       breaks = q_breaks,
                       labels = q_labels)+
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
  geom_sf(data = study_area) +
  geom_tile(data = pr.int,
            mapping = aes(x = x,
                         y = y,
                         fill = q0.975)) +
  scale_fill_viridis_c(limits = q_breaks,
                       breaks = q_breaks,
                       labels = q_labels)+
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

# show that these maps are misleading by interpreting
# them as an intensity
cell_area = st_area(study_area)/nrow(pr.int)*100
sum(pr.int$q0.025*cell_area)
sum(pr.int$q0.975*cell_area)
# these are not supported by posterior for abundance (see script N_posterior.R)

### plot three realisations of posterior intensity field ####

set.seed(1989)
draw1 = predict(fit, pxl, ~ exp(grf + Intercept), n.samples = 1)
draw2 = predict(fit, pxl, ~ exp(grf + Intercept), n.samples = 1)
draw3 = predict(fit, pxl, ~ exp(grf + Intercept), n.samples = 1)
draw1$x = pr.int$x
draw1$y = pr.int$y
draw2$x = pr.int$x
draw2$y = pr.int$y
draw3$x = pr.int$x
draw3$y = pr.int$y

draw_all = c(draw1$mean,
      draw2$mean,
      draw3$mean)

viridisscale = scale_fill_viridis_c(limits = range(draw_all))

p1 = ggplot() +
  geom_sf(data = study_area) +
  geom_tile(data = draw1,
            mapping = aes(x = x,
                         y = y,
                         fill = mean)) +
  viridisscale +
  xlab("Easting") +
  ylab("Northing") +
  theme_void() +
  theme(legend.position = "none")

p2 = ggplot() +
  geom_sf(data = study_area) +
  geom_tile(data = draw2,
            mapping = aes(x = x,
                         y = y,
                         fill = mean)) +
  viridisscale +
  xlab("Easting") +
  ylab("Northing") +
  theme_void() +
  theme(legend.position = "none")

p3 = ggplot() +
  geom_sf(data = study_area) +
  geom_tile(data = draw3,
            mapping = aes(x = x,
                         y = y,
                         fill = mean)) +
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


