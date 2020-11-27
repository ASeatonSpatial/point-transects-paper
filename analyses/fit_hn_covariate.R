library(inlabru)
library(INLA)
library(dplyr)

init.tutorial()

#### Load Data ####

data_path = here::here("analyses", "data")
realobs <- readRDS(here::here(data_path, "obs_extended_no_crs.RDS"))
study_area = readRDS(here::here(data_path, "study_area_extended_no_crs.RDS"))
samplers = readRDS(here::here(data_path, "samplers_extended_no_crs.RDS"))
mesh = readRDS(here::here(data_path, "mesh_extended_no_crs.RDS"))

### Plot detections for different habitat types

df = realobs@data
df$HabCode = as.factor(df$HabCode)
df$Stratum = as.factor(df$Stratum)

### Plot detections by Stratum
df %>%
  group_by(Stratum) %>%
  summarise(n = n())

df %>%
  ggplot() +
  geom_histogram(aes(x = distance, y = after_stat(density)), bins = 15) +
  facet_wrap(~ Stratum)

# These look very slightly different.  Could try fitting a different detection
# function to each

# create binary variables
realobs$stratum = ifelse(realobs$Stratum == "CF", 0, 1)
samplers$stratum = ifelse(samplers$Stratum == "CF", 0, 1)

realobs$habcode = ifelse(realobs$HabCode == "Mesic", 1, 0)
samplers$habcode = ifelse(samplers$HabCode == "Mesic", 1, 0)

#### Specify detection function  ####

# x is a binary variable (2 level factor)
# lsig_x ajdusts for when x == 1
# lsig_int is for x == 0

# Log half-normal
log_hn = function(distance, x, lsig_int, lsig_x){
  -0.5*(distance/exp(lsig_int + lsig_x*x))^2
}

# Half-normal
hn <- function(distance, x, lsig_int, lsig_x) {
  exp(log_hn(distance, x, lsig_int, lsig_x))
}

# Include r for switching to polar coords
# 2*pi for not knowing angle theta
dsamp = function(distance, x, lsig_int, lsig_x){
  log(2*pi) + log(distance) + log_hn(distance, x, lsig_int, lsig_x)
}


# define distance domain and starting values
# for detection scale param intercept
W <- 58/1000   # transect radius
distance_domain <- inla.mesh.1d(seq(.Machine$double.eps, W, length.out = 30))
starting_values <- list(lsig_int = 3.36 - log(1000))    # lsig_x = 0 implied here

# SPDE:
matern <- inla.spde2.pcmatern(mesh,
                              prior.sigma = c(2, 0.01),
                              prior.range = c(300/1000, 0.01))   # was 300 before, now km units so divide

# model components
cmp <- ~ grf(main = coordinates, model = matern) +
  lsig_int(1) + lsig_x(1) + Intercept(1)

# Predictor formula:
fml <- coordinates + distance ~ grf +
  dsamp(distance, stratum, lsig_int, lsig_x) +
  Intercept


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

INLA:::summary.inla(fit)

### Try it on HabCode variable?

# Not enough Mesic detections?  37 is low for distance sampling but
# it might work

df %>%
  group_by(HabCode) %>%
  summarise(n = n())

df %>%
  ggplot() +
  geom_histogram(aes(x = distance, y  = after_stat(density)), bins = 15) +
  facet_wrap(~ HabCode)

# Predictor formula:
fml <- coordinates + distance ~ grf +
  dsamp(distance, habcode, lsig_int, lsig_x) +
  Intercept

fit2 <-  lgcp(components = cmp,
              data = realobs,
              samplers = samplers,
              domain = list(coordinates = mesh,
                            distance = distance_domain),
              formula = fml,
              options = list(bru_max_iter = 40,
                             bru_result = starting_values,
                             # bru_method = list(taylor = "legacy"),
                             bru_verbose = 3))

INLA:::summary.inla(fit2)

# Plot detection functions
distdf <- data.frame(distance = seq(.Machine$double.eps,W,length.out=100),
                     habcode = 0)   # Wet habitat

hn_pred_wet <- predict(fit, distdf, ~hn(distance, habcode, lsig_int, lsig_x))

distdf <- data.frame(distance = seq(.Machine$double.eps,W,length.out=100),
                     habcode = 1)   # Mesic habitat

hn_pred_mesic <- predict(fit, distdf, ~hn(distance, habcode, lsig_int, lsig_x))

ggplot() +
  geom_line(data = hn_pred_wet,
            mapping = aes(x = distance, y = mean)) +
  geom_ribbon(data = hn_pred_wet,
               mapping = aes(x = distance, ymax = q0.975, ymin = q0.025),
              alpha = 0.25) +
  geom_line(data = hn_pred_mesic,
            mapping = aes(x = distance, y = mean), colour = "blue") +
  geom_ribbon(data = hn_pred_mesic,
              mapping = aes(x = distance, ymax = q0.975, ymin = q0.025),
              alpha = 0.25, fill = "blue") +
  ggtitle("blue = mesic habitat, black = wet habitat")


# Is it surprising that mesic with 37 obs seems to be more
# precise than wet which has many more?


#### Try both at once:

# x1 and x2 are binary variables (2 level factor)
# lsig_x1 adjusts for when x1 == 1
# lsig_x2 adjusts for when x2 == 1
# lsig_int is when x1 == 0 & x2 == 0

# Log half-normal
log_hn = function(distance, x1, x2, lsig_int, lsig_x1, lsig_x2){
  -0.5*(distance/exp(lsig_int + lsig_x1*x1 + lsig_x2*x2))^2
}

# Half-normal
hn <- function(distance, x1, x2, lsig_int, lsig_x1, lsig_x2) {
  exp(log_hn(distance, x1, x2, lsig_int, lsig_x1, lsig_x2))
}

# Include r for switching to polar coords
# 2*pi for not knowing angle theta
dsamp = function(distance, x1, x2, lsig_int, lsig_x1, lsig_x2){
  log(2*pi) + log(distance) + log_hn(distance, x1, x2, lsig_int, lsig_x1, lsig_x2)
}

cmp <- ~ grf(main = coordinates, model = matern) +
  lsig_int(1) + lsig_x1(1) + lsig_x2(1) + Intercept(1)

# lsig_x1 is for habcode
# lsig_x2 is for stratum
fml <- coordinates + distance ~ grf +
  dsamp(distance, habcode, stratum, lsig_int, lsig_x1, lsig_x2) +
  Intercept

fit3 <-  lgcp(components = cmp,
              data = realobs,
              samplers = samplers,
              domain = list(coordinates = mesh,
                            distance = distance_domain),
              formula = fml,
              options = list(bru_max_iter = 40,
                             bru_result = starting_values,
                             # bru_method = list(taylor = "legacy"),
                             bru_verbose = 3))

INLA:::summary.inla(fit3)
