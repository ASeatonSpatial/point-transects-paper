## Fit model and save model object

# Load required packages
library(ggplot2)
library(sf)
library(inlabru)
library(INLA)

#### Save Fitted Model Object? ####
#### WARNING: THIS WILL OVERWRITE EXISTING FITTED MODEL OBJECT
to_save <- TRUE
model_path <- here::here(
  "R",
  "fitted_models",
  "fitted_model.RDS"
)

#### Produce intermediate plots? ####
to_plot <- TRUE

#### Load Data ####
data_path <- here::here("R", "data")
obs <- readRDS(here::here(data_path, "obs.RDS"))
study_area <- readRDS(here::here(data_path, "study_area.RDS"))
samplers <- readRDS(here::here(data_path, "samplers.RDS"))
mesh_path <- here::here(data_path, "mesh.RDS")
mesh <- readRDS(mesh_path)

#### Specify detection function  ####

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

# Plot the data
if (to_plot) {
  g1 <- ggplot() +
    gg(mesh) +
    geom_sf(data = study_area) +
    geom_sf(data = samplers) +
    geom_sf(data = obs, colour = "green") +
    coord_sf()
  g1
}

#### Specify and fit model ####

# SPDE:
matern <- inla.spde2.pcmatern(mesh,
  prior.sigma = c(2, 0.01),
  prior.range = c(300 / 1000, 0.01)
)

cmp <- ~ grf(main = geometry, model = matern) +
  lsig(1) + Intercept(1)

# Predictor formula:
fml <- geometry + distance ~ grf +
  dsamp(distance, lsig) +
  log(2 * pi) + # 2*pi offset for not knowing angle theta
  Intercept

W <- 58 / 1000 # transect radius, units km
distance_domain <- fm_mesh_1d(W * seq(.Machine$double.eps, 1, length.out = 30)^2)
starting_values <- list(lsig = 3.36 - log(1000))

lik <- bru_obs(
  data = obs,
  family = "cp",
  samplers = samplers,
  domain = list(
    geometry = mesh,
    distance = distance_domain
  ),
  formula = fml
)

fit <- bru(
  components = cmp,
  lik,
  options = list(
    bru_max_iter = 40,
    bru_initial = starting_values,
    bru_verbose = 3,
    control.inla = list(int.strategy = "eb")
  )
)

#### Save model object ####

if (to_save) {
  if (!dir.exists(dirname(model_path))) {
    dir.create(dirname(model_path))
  }
  saveRDS(object = fit, file = model_path)
}
