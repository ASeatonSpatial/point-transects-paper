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
to_save = FALSE
file_name = "akepa_2002_fitted_model_new_inlabru.RDS"


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

# What are the shorter distances between transects?
# (xrange <- study_area@bbox["x",2] - study_area@bbox["x",1])
# (yrange <- study_area@bbox["y",2] - study_area@bbox["y",1])
# sort(gDistance(pts, byid = TRUE), decreasing = TRUE)

matern <- inla.spde2.pcmatern(mesh, 
                              prior.sigma = c(2, 0.01),    
                              prior.range = c(300/1000, 0.01))   # was 300 before, now km units so divide

cmp <- ~ grf(main = coordinates, model = matern) + 
  lsig(1) + Intercept(1)

# Predictor formula:
fml <- coordinates + distance ~ grf +
  dsamp(distance, lsig) + 
  log(2*pi) +   # 2*pi offset for not knowing angle theta
  Intercept



W <- 58/1000   # transect radius
distance_domain <- inla.mesh.1d(seq(.Machine$double.eps, W, length.out = 30))


# inlabru:::iinla.setOption(iinla.verbose = TRUE)
# inlabru:::iinla.getOption("iinla.verbose")

starting_values <- list(lsig = 3.36 - log(1000))

fit <-  lgcp(components = cmp, 
             data = realobs,
             samplers = samplers,
             domain = list(coordinates = mesh,
                           distance = distance_domain),
             formula = fml,
             options = list(bru_max_iter = 40,
                            bru_result = starting_values,
                            # bru_method = list(taylor = "legacy"),
                            bru_verbose = 3))

summary(fit)

#### Save model object ####

if (to_save) saveRDS(object = fit, file = file_name)

  