## Fit model and save model object

# Load required packages
library(inlabru)
library(INLA)  # need testing version of inla
# library(rgeos)
# library(raster)  # raster::crs() easy way to remove CRS information
# library(maptools)

init.tutorial()

#### Save Fitted Model Object? ####

#### WARNING!!!! THIS WILL OVERWRITE EXISTING OBJECT!!!!
to_save = TRUE
file_name = "akepa_2002_fitted_model.RDS"


#### Load Data ####

realobs <- readRDS("data/obs_extended_no_crs.RDS")
samplers <- readRDS("data/samplers_extended_no_crs.RDS")
study_area <- readRDS("data/study_area_extended_no_crs.RDS")
mesh <- readRDS("data/mesh_extended_no_crs.RDS")

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
g1 <- ggplot() +
  gg(mesh) + 
  gg(study_area) +
  gg(samplers) +
  gg(realobs, colour = "green") +
  coord_equal()

g1

#### Specify and fit model ####

# SPDE:
# (xrange <- study_area@bbox["x",2] - study_area@bbox["x",1])
# (yrange <- study_area@bbox["y",2] - study_area@bbox["y",1])
# sort(gDistance(pts, byid = TRUE), decreasing = TRUE)
matern <- inla.spde2.pcmatern(mesh, 
                              prior.sigma = c(2, 0.01),    
                              prior.range = c(130, 0.01))

cmp <- ~ grf(map = coordinates, model = matern) + 
  lsig + Intercept

# Predictor formula:
fml <- coordinates + distance ~ grf +
  dsamp(distance, lsig) + 
  log(2*pi) +   # 2*pi offset for not knowing angle theta
  Intercept


starting_values <- data.frame(grf = 0, lsig = 3.36, Intercept = 0)

W <- 60   # transect radius
distance_domain <- seq(.Machine$double.eps, W, length.out = 30)

fit <-  lgcp(components = cmp, 
             data = realobs,
             samplers = samplers,
             # domain = list(coordinates = mesh, distance = distance_domain),
             domain = list(distance = distance_domain),
             formula = fml,
             options = list(result = starting_values, max.iter = 40))

summary(fit)


#### Save model object ####

if (to_save) saveRDS(object = fit, file = file_name)

  