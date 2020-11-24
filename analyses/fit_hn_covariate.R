library(inlabru)
library(INLA)
library(dplyr)
library(patchwork)

init.tutorial()

#### Load Data ####

realobs <- readRDS("data/obs_extended_no_crs.RDS")
samplers <- readRDS("data/samplers_extended_no_crs.RDS")
study_area <- readRDS("data/study_area_extended_no_crs.RDS")
mesh <- readRDS("data/mesh_extended_no_crs.RDS")

### Plot detections for different habitat types

df = realobs@data
df$HabCode = as.factor(df$HabCode)
df$Stratum


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

# create 0-1 version of stratum factor
realobs$stratum = ifelse(realobs$Stratum == "CF", 0, 1)
samplers$stratum = ifelse(samplers$Stratum == "CF", 0, 1)

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

#### Specify and fit model ####

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

W <- 58/1000   # transect radius
distance_domain <- inla.mesh.1d(seq(.Machine$double.eps, W, length.out = 30))
starting_values <- list(lsig_int = 3.36 - log(1000))    # lsig_strat = 0 implied here

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

# create binary variable
realobs$habcode = ifelse(realobs$HabCode == "Mesic", 1, 0)
samplers$habcode = ifelse(samplers$HabCode == "Mesic", 1, 0)

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
  