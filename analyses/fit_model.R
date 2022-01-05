## Fit model and save model object

# Load required packages
library(ggplot2)
library(sp)
library(inlabru)
library(INLA)  

#### Save Fitted Model Object? ####
#### WARNING: THIS WILL OVERWRITE EXISTING FITTED MODEL OBJECT 
to_save = TRUE
model_path = here::here("analyses", "fitted_model.RDS")

#### Produce intermediate plots? ####
to_plot = FALSE

#### Load Data ####
data_path = here::here("analyses", "data")
realobs <- readRDS(here::here(data_path, "obs_extended_no_crs.RDS"))
study_area = readRDS(here::here(data_path, "study_area_extended_no_crs.RDS"))
samplers = readRDS(here::here(data_path, "samplers_extended_no_crs.RDS"))
mesh = readRDS(here::here(data_path, "mesh_extended_no_crs.RDS"))

#### Specify detection function  ####

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

# Plot the data
if (to_plot){
  g1 <- ggplot() +
    gg(mesh) +
    gg(study_area) +
    gg(samplers) +
    gg(realobs, colour = "green") +
    coord_equal()
  g1
}

#### Specify and fit model ####

# SPDE:
matern <- inla.spde2.pcmatern(mesh,
                              prior.sigma = c(2, 0.01),
                              prior.range = c(300/1000, 0.01))

cmp <- ~ grf(main = coordinates, model = matern) +
  lsig(1) + Intercept(1)

# Predictor formula:
fml <- coordinates + distance ~ grf +
  dsamp(distance, lsig) +
  log(2*pi) +   # 2*pi offset for not knowing angle theta
  Intercept

W <- 58/1000   # transect radius, units km 
distance_domain <- inla.mesh.1d(seq(.Machine$double.eps, W, length.out = 30))
starting_values <- list(lsig = 3.36 - log(1000))

fit <-  lgcp(components = cmp,
             data = realobs,
             samplers = samplers,
             domain = list(coordinates = mesh,
                           distance = distance_domain),
             formula = fml,
             options = list(bru_max_iter = 40,
                            bru_result = starting_values,
                            bru_verbose = 3,
#                            inla.mode = "experimental",
                            control.inla = list(int.strategy = "eb")))

#### Save model object ####

if (to_save) saveRDS(object = fit, file = model_path)

