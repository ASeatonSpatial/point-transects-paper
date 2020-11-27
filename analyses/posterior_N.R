library(INLA)
library(inlabru)
library(rgeos)

# set ggplot theme
theme_set(theme_minimal())

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

#### A new way: don't predict within transects ####
# Unfinished - need to create my own buffer around the point locations
# gBuffer doesn't create things that are close enough to a circle for my liking

# Total surveyed region as % of study area?
nrow(samplers@data)*pi*58^2 / gArea(study_area) * 100  # 58m radius for point transects

# Buffer samplers to create transects
samplers_buffered = gBuffer(samplers, width = 58, byid = TRUE)
gArea(samplers_buffered) - nrow(samplers@coords)*pi*58^2   # this doesn't look great

# Remove transects from study area
study_area_no_samplers = gDifference(study_area, samplers_buffered)
gArea(study_area) - gArea(study_area_no_samplers) - gArea(samplers_buffered)  # should be zero but isn't



#### The original way of doing it ####
# Posterior abundance

# try the new way:
ips = ipoints(study_area_no_samplers, domain = mesh)
ggplot() +
  gg(study_area) +
  gg(ips) +
  coord_equal()

# how does this compare?
ips2 = ipoints(study_area, domain = mesh)
nrow(ips)
nrow(ips2)

# source("AS.Nposterior.R")
predpts = ipoints(study_area, mesh)
# abnc = AS.Nposterior(fit, predpts,
#                      n.sd = 5,
#                      n.grid = 100,
#                      n.samples = 100)   # n.samples = 1000 default

# First get range using mean + sd in expected abundance
n.sd = 5
n.grid = 100
Lambda <- predict(fit, predpts, ~ sum(weight * exp(grf + Intercept)))
Nlow <- round(Lambda$mean - n.sd*Lambda$sd)
Nhigh <- round(Lambda$mean + n.sd*Lambda$sd)
by <- round((Nhigh - Nlow)/n.grid)   # should give roughly n.grid
Nrange <- seq(Nlow, Nhigh, by = by)

getNpost = function(){
  predict(fit, predpts,
          ~ data.frame(N = Nrange,
                       dpois(Nrange, lambda = sum(weight * exp(grf + Intercept)))),
          n.samples = 100)
}

n.blah = 100
Npost = matrix(NA, nrow = n.grid+1, ncol = n.blah)
for (i in 1:n.blah) {
  print(paste("i = ", i))
  Npost[,i] = getNpost()$mean
  gc()
}

Npost.mean = data.frame(N=Nrange, mean = rowMeans(Npost))
#  saveRDS(Npost.mean, file = here::here("analyses", "Npost.mean.RDS"))

ggplot() +
  geom_line(data = Npost.mean, aes(x = N, y = mean)) +
  ylab("Posterior probability\n") +
  xlab("\nAbundance") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        plot.margin = margin(0.2,0.5,0.2,0.2, "in"))

ggsave(here::here(fig_path, "realized_abundance_posterior.png"),
       width = 5, height = 5, units = "in")

Npost$Nexp = dpois(Npost$N, lambda = abnc$Nexp$mean)
ggplot() +
  geom_line(data = Npost, aes(x = N, y = mean)) +
  geom_line(data = Npost, aes(x = N, y = Nexp), linetype = "dashed")

ggsave(here::here(fig_path,"realized_abundance_vs_exp.png"),
       width = 5, height = 5, units = "in")
