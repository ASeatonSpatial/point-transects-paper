# Produce summary figures of the posterior intensity field

library(INLA)
library(inlabru)
library(ggplot2)
library(scales)
library(cowplot)
library(sf)

set.seed(9701071)

# statsoc doc class textwidth in cm
tw <- 14.69785
twi <- 5.7865551181 # pdf() needs inches

# axis title text size
ax.size <- 12

#### Fitted model  ####
model_path <- here::here(
  "R",
  "fitted_models",
  "fitted_model.RDS"
)
fit <- readRDS(model_path)

#### Other things ####
data_path <- here::here("R", "data")
study_area <- readRDS(here::here(data_path, "study_area.RDS"))
samplers <- readRDS(here::here(data_path, "samplers.RDS"))
mesh <- readRDS(here::here(data_path, "mesh.RDS"))

#### Where to save figures ####
fig_path <- here::here("figures")

#### ggplot theme ####
# theme_set(theme_minimal())
theme_set(theme_classic())

#### Specify detection functions  ####

# Log half-normal
log_hn <- function(distance, lsig) {
  -0.5 * (distance / exp(lsig))^2
}

# Half-normal
hn <- function(distance, lsig) exp(log_hn(distance, lsig))

# Include r for switching to polar coords
dsamp <- function(distance, lsig) {
  log(distance) + log_hn(distance, lsig)
}

#### Model evaluation ####
summary(fit)

##### Posterior intensity plots #####

# Prediction locations
pxl <- fm_pixels(mesh,
  dims = c(300, 300),
  mask = study_area
)

# map units are per km^2 = 1,000,000 m^2
# I want intensity per hectare = 10,000 m^2
# so divide by 100
pr.int <- predict(fit, pxl, ~ exp(grf + Intercept) / 100,
  n.samples = 1000
)

# mean
lower <- min(pr.int$mean)
upper <- max(pr.int$mean)
mean_breaks <- c(lower, upper)
mean_labels <- c(
  signif(lower, digits = 3),
  signif(upper, digits = 3)
)

p1 <- ggplot() +
  geom_sf(data = study_area) +
  gg(
    data = pr.int,
    mapping = aes(
      fill = mean,
      colour = mean
    ),
    geom = "tile"
  ) +
  scale_fill_viridis_c(
    breaks = mean_breaks,
    labels = mean_labels
  ) +
  scale_colour_viridis_c(
    breaks = mean_breaks,
    labels = mean_labels
  ) +
  coord_sf() +
  xlab("Easting") +
  ylab("Northing") +
  theme_void() +
  #  ggtitle("Intensity [count per hectare]") +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(fill = guide_colourbar(
    title.vjust = 0.85,
    title.theme = element_text(size = ax.size - 2),
    label.theme = element_text(size = ax.size - 2)
  ))

p1

# cv
pr.int$cv <- pr.int$sd / pr.int$mean
lower <- min(pr.int$cv)
upper <- max(pr.int$cv)
cv_breaks <- c(lower, upper)
cv_labels <- c(
  signif(lower, digits = 3),
  signif(upper, digits = 3)
)

p2 <- ggplot() +
  geom_sf(data = study_area) +
  gg(
    data = pr.int,
    mapping = aes(
      fill = cv,
      colour = cv
    ),
    geom = "tile"
  ) +
  scale_fill_viridis_c(
    breaks = cv_breaks,
    labels = cv_labels
  ) +
  scale_colour_viridis_c(
    breaks = cv_breaks,
    labels = cv_labels
  ) +
  coord_sf() +
  xlab("Easting") +
  ylab("Northing") +
  theme_void() +
  #  ggtitle("Coefficient of variation [sd/mean]") +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(fill = guide_colourbar(
    title.vjust = 0.85,
    title.theme = element_text(size = ax.size - 2),
    label.theme = element_text(size = ax.size - 2)
  ))

p2

# sd
lower <- min(pr.int$sd)
upper <- max(pr.int$sd)
sd_breaks <- c(lower, upper)
sd_labels <- format(c(
  signif(lower, digits = 3),
  signif(upper, digits = 3)
))

p3 <- ggplot() +
  geom_sf(data = study_area) +
  gg(
    data = pr.int,
    mapping = aes(
      fill = sd,
      colour = sd
    ),
    geom = "tile"
  ) +
  scale_fill_viridis_c(
    breaks = sd_breaks,
    labels = sd_labels
  ) +
  scale_colour_viridis_c(
    breaks = sd_breaks,
    labels = sd_labels
  ) +
  coord_sf() +
  xlab("Easting") +
  ylab("Northing") +
  theme_void() +
  #  ggtitle("Standard deviation [counts per hectare]") +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(fill = guide_colourbar(
    title.vjust = 0.85,
    title.theme = element_text(size = ax.size - 2),
    label.theme = element_text(size = ax.size - 2)
  ))

p3

# Adjust each plot's legend guide to include some margin
# Numbers getting clipped otherwise.

# Adjust each plot's legend guide to avoid stretching the color ramps
p1 <- p1 +
  theme(legend.margin = margin(t = 0, r = 10, b = 0, l = 10)) + # Add margin to legend
  guides(fill = guide_colourbar(
    title.vjust = 0.85,
    title.theme = element_text(size = ax.size - 2),
    label.theme = element_text(size = ax.size - 2),
    barwidth = unit(2, "cm"), # Set a fixed width for the color bar
    barheight = unit(0.5, "cm") # Set a fixed height for the color bar
  ))

p2 <- p2 +
  theme(legend.margin = margin(t = 0, r = 10, b = 0, l = 10)) +
  guides(fill = guide_colourbar(
    title.vjust = 0.85,
    title.theme = element_text(size = ax.size - 2),
    label.theme = element_text(size = ax.size - 2),
    barwidth = unit(2, "cm"), # Consistent width for the color bar
    barheight = unit(0.5, "cm") # Consistent height for the color bar
  ))

p3 <- p3 +
  theme(legend.margin = margin(t = 0, r = 10, b = 0, l = 10)) +
  guides(fill = guide_colourbar(
    title.vjust = 0.85,
    title.theme = element_text(size = ax.size - 2),
    label.theme = element_text(size = ax.size - 2),
    barwidth = unit(2, "cm"), # Consistent width for the color bar
    barheight = unit(0.5, "cm") # Consistent height for the color bar
  ))

pdf(
  file = here::here(fig_path, "intensity_mean_cv_sd.pdf"),
  width = twi, height = twi / 2
)
plot_grid(NULL, p1, NULL, p2, NULL, p3, NULL,
  labels = c("", "A", "", "B", "", "C", ""),
  rel_widths = c(0.1, 1, 0.4, 1, 0.4, 1, 0.1),
  ncol = 7
)
dev.off()


##### lower and upper quantiles ####

# scale legend breaks
v <- c(pr.int$q0.025, pr.int$q0.975)
lower <- min(v)
upper <- max(v)
q_breaks <- c(lower, upper)
q_labels <- sd_labels <- format(
  c(
    signif(lower, digits = 3),
    signif(upper, digits = 3)
  ),
  scientific = TRUE
)

p2 <- ggplot() +
  geom_sf(data = study_area) +
  gg(
    data = pr.int,
    mapping = aes(
      fill = q0.025,
      colour = q0.025
    ),
    geom = "tile"
  ) +
  coord_sf() +
  scale_fill_viridis_c("",
    limits = q_breaks,
    breaks = q_breaks,
    labels = q_labels
  ) +
  scale_colour_viridis_c(
    limits = q_breaks,
    breaks = q_breaks,
    labels = q_labels,
    guide = "none"
  ) +
  xlab("Easting") +
  ylab("Northing") +
  theme_void() +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(fill = guide_colourbar(
    title = "",
    title.vjust = 0.95,
    title.theme = element_text(size = 18),
    label.theme = element_text(size = 16)
  ))

p2

p3 <- ggplot() +
  geom_sf(data = study_area) +
  gg(
    data = pr.int,
    mapping = aes(
      fill = q0.975,
      colour = q0.975
    ),
    geom = "tile"
  ) +
  coord_sf() +
  scale_fill_viridis_c("",
    limits = q_breaks,
    breaks = q_breaks,
    labels = q_labels
  ) +
  scale_colour_viridis_c(
    limits = q_breaks,
    breaks = q_breaks,
    labels = q_labels,
    guide = "none"
  ) +
  xlab("Easting") +
  ylab("Northing") +
  theme_void() +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(fill = guide_colourbar(
    title = "",
    title.vjust = 0.95,
    title.theme = element_text(size = 18),
    label.theme = element_text(size = 16)
  ))

p3

p2 <- p2 +
  theme(legend.margin = margin(t = 0, r = 10, b = 0, l = 10)) +
  guides(fill = guide_colourbar(
    title.vjust = 0.85,
    title.theme = element_text(size = ax.size - 2),
    label.theme = element_text(size = ax.size - 2),
    barwidth = unit(2, "cm"), # Consistent width for the color bar
    barheight = unit(0.5, "cm") # Consistent height for the color bar
  ))

p3 <- p3 +
  theme(legend.margin = margin(t = 0, r = 10, b = 0, l = 10)) +
  guides(fill = guide_colourbar(
    title.vjust = 0.85,
    title.theme = element_text(size = ax.size - 2),
    label.theme = element_text(size = ax.size - 2),
    barwidth = unit(2, "cm"), # Consistent width for the color bar
    barheight = unit(0.5, "cm") # Consistent height for the color bar
  ))

library(patchwork)
pdf(
  file = here::here(fig_path, "intensity_quantiles.pdf"),
  width = twi / 1.5, height = twi / 2
)
# plot_grid(NULL, p2, NULL, p3, NULL,
#           labels = c("", "A", "", "B", ""),
#           rel_widths = c(0.05, 1, 0.15, 1, 0.05),
#           ncol = 5)
wrap_plots(p2, p3) +
  plot_annotation(tag_levels = "A") +
  plot_layout(
    heights = c(1, 1), # Adjust heights for equal space
    nrow = 1
  ) + # 2 rows for vertical stacking
  theme(plot.margin = margin(10, 10, 10, 10)) &
  theme(legend.justification = c(0.5, 0)) # Adjust margins for more whitespace
dev.off()

# show that these maps are misleading by interpreting
# them as an intensity
cell_area <- st_area(study_area) / nrow(pr.int) * 100
sum(pr.int$q0.025 * cell_area)
sum(pr.int$q0.975 * cell_area)
# these are not supported by posterior for abundance (see script N_posterior.R)

### plot three realisations of posterior intensity field ####

set.seed(19822139)
draws <- generate(fit, pxl, ~ exp(grf + Intercept), n.samples = 3)
draws <- cbind(rbind(pxl, pxl, pxl),
  value = as.vector(draws),
  draw = rep(1:3, each = nrow(pxl))
)
draw_all <- draws$value

viridisscale <- scale_fill_viridis_c("",
  limits = range(draw_all)
)

plts <- list()
for (the_draw in 1:3) {
  plts[[the_draw]] <-
    ggplot() +
    geom_sf(data = study_area) +
    gg(
      data = draws |> dplyr::filter(draw == the_draw),
      mapping = aes(
        fill = value,
        colour = value
      ),
      geom = "tile"
    ) +
    viridisscale +
    scale_colour_viridis_c("",
      limits = range(draw_all),
      guide = "none"
    ) +
    xlab("Easting") +
    ylab("Northing") +
    theme_void() +
    theme(legend.position = "none")
}


pdf(
  file = here::here(fig_path, "intensity_realized.pdf"),
  width = twi, height = twi / 2
)
plot_grid(NULL, plts[[1]], NULL, plts[[2]], NULL, plts[[3]], NULL,
  labels = c("", "A", "", "B", "", "C", ""),
  rel_widths = c(0.1, 1, 0.4, 1, 0.4, 1, 0.1),
  ncol = 7
)
dev.off()

#### Further plots are in eval_spde.R ####
