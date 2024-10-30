# Posterior abundance plot

library(INLA)
library(inlabru)
library(ggplot2)
library(sf)
library(patchwork)
library(dplyr)
source(here::here("R", "gg.R"))   # some small edits to inlabru gg methods

# set ggplot theme
theme_set(theme_minimal())

#### Fitted model  ####
model_path = here::here("R",
                        "fitted_models",
                        "fitted_model.RDS")
fit = readRDS(model_path)

#### Other things ####
data_path = here::here("R", "data")
study_area = readRDS(here::here(data_path, "study_area.RDS"))
samplers = readRDS(here::here(data_path, "samplers.RDS"))
obs = readRDS(here::here(data_path, "obs.RDS"))
mesh = readRDS(here::here(data_path, "mesh.RDS"))

#### Where to save figures ####
fig_path = here::here("figures")

# integration points
ips = fm_int(domain = mesh,
             samplers = study_area)

# test to get range for N posterior
test = generate(fit,
                newdata = ips,
                formula = ~ sum(weight * exp(Intercept + grf)),
                n.samples = 500)

# define range of N
N_min = round(mean(test) - 4.5*sd(test))
N_max = round(mean(test) + 4.5*sd(test))
by = round((N_max - N_min)/100)
N_seq = seq(N_min, N_max, by = by) # should be roughly length 100

n.mc = 20000
N_post = predict(fit,
                 newdata = ips,
                 ~ {
                     c(
                       predict_everywhere = dpois(N_seq,
                                                 lambda = sum(weight * exp(Intercept + grf))
                                                 )
                     )
                 },
                 n.samples = n.mc)

N_post$N = N_seq

# Just the abundance posterior
N_post %>%
  ggplot() +
  geom_line(aes(x = N, y = mean)) +
  geom_ribbon(aes(x = N,
                  ymin = mean - 2 * sd / sqrt(n.mc),
                  ymax = mean + 2 * sd / sqrt(n.mc)),
              alpha = 0.2) +
  ylab("p(N)")

ggsave(filename = here::here(fig_path, "N_posterior.png"),
       width = 5, height = 5, units = "in")

Nweight = (N_max - N_min) / length(N_seq)
sum(Nweight * N_post$mean) # should be close to 1

# get mean density per hectare
# this allows comparison with Figure 2 in Camp et al 2020
whole_area = st_area(study_area)*100   # hectares
N_mode = N_seq[which(N_post$mean == max(N_post$mean))]  # posterior mode
N_mode / whole_area
5500 / whole_area
4200 / whole_area  # some random end points just to get a feeling
7500 / whole_area

# just convert the plot scale to show density
N_post %>%
  mutate(dens = N/whole_area) %>%
  ggplot() +
  geom_line(aes(x = dens, y = mean)) +
  geom_ribbon(aes(x = dens,
                  ymin = mean - 2 * sd / sqrt(n.mc),
                  ymax = mean + 2 * sd / sqrt(n.mc)),
              alpha = 0.2) +
  ylab("posterior probability") +
  xlab("akepa density per hectare")

ggsave(filename = here::here(fig_path, "density_posterior.png"),
       width = 5, height = 5, units = "in")

