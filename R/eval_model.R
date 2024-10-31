# Produce summary figures of the posterior intensity field

library(INLA)
library(inlabru)
library(ggplot2)
library(scales)
library(cowplot)
library(sf)

set.seed(9701071)

# statsoc doc class textwidth in cm
tw = 14.69785
twi = 5.7865551181  # pdf() needs inches

# axis title text size
ax.size = 12

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
                         fill = mean,
                         colour = mean)) +
  scale_fill_viridis_c(breaks = mean_breaks,
                       labels = mean_labels) +
  scale_colour_viridis_c(breaks = mean_breaks,
                         labels = mean_labels) +
  xlab("Easting") +
  ylab("Northing") +
  theme_void() +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(fill = guide_colourbar(title.vjust = 0.85,
                                title.theme = element_text(size = ax.size-2),
                                label.theme = element_text(size = ax.size-2)))

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
                         fill = cv,
                         colour = cv)) +
  scale_fill_viridis_c(breaks = cv_breaks,
                       labels = cv_labels) +
  scale_colour_viridis_c(breaks = cv_breaks,
                       labels = cv_labels) +
  xlab("Easting") +
  ylab("Northing") +
  theme_void() +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(fill = guide_colourbar(title.vjust = 0.85,
                                title.theme = element_text(size = ax.size-2),
                                label.theme = element_text(size = ax.size-2)))

p2

# sd
lower = min(pr.int$sd)
upper = max(pr.int$sd)
sd_breaks = c(lower, upper)
sd_labels = format(c(signif(lower, digits = 3),
                     signif(upper, digits = 3)))

p3 = ggplot() +
  geom_sf(data = study_area) +
  geom_tile(data = pr.int,
            mapping = aes(x = x,
                         y = y,
                         fill = sd,
                         colour = sd)) +
  scale_fill_viridis_c(breaks = sd_breaks,
                       labels = sd_labels) +
  scale_colour_viridis_c(breaks = sd_breaks,
                       labels = sd_labels) +
  xlab("Easting") +
  ylab("Northing") +
  theme_void() +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(fill = guide_colourbar(title.vjust = 0.85,
                                title.theme = element_text(size = ax.size-2),
                                label.theme = element_text(size = ax.size-2)))

p3

# Adjust each plot's legend guide to include some margin
# Numbers getting clipped otherwise.

# Adjust each plot's legend guide to avoid stretching the color ramps
p1 = p1 +
  theme(legend.margin = margin(t = 0, r = 10, b = 0, l = 10)) +  # Add margin to legend
  guides(fill = guide_colourbar(
    title.vjust = 0.85,
    title.theme = element_text(size = ax.size-2),
    label.theme = element_text(size = ax.size-2),
    barwidth = unit(2, "cm"),  # Set a fixed width for the color bar
    barheight = unit(0.5, "cm") # Set a fixed height for the color bar
  ))

p2 = p2 +
  theme(legend.margin = margin(t = 0, r = 10, b = 0, l = 10)) +
  guides(fill = guide_colourbar(
    title.vjust = 0.85,
    title.theme = element_text(size = ax.size-2),
    label.theme = element_text(size = ax.size-2),
    barwidth = unit(2, "cm"),  # Consistent width for the color bar
    barheight = unit(0.5, "cm") # Consistent height for the color bar
  ))

p3 = p3 +
  theme(legend.margin = margin(t = 0, r = 10, b = 0, l = 10)) +
  guides(fill = guide_colourbar(
    title.vjust = 0.85,
    title.theme = element_text(size = ax.size-2),
    label.theme = element_text(size = ax.size-2),
    barwidth = unit(2, "cm"),  # Consistent width for the color bar
    barheight = unit(0.5, "cm") # Consistent height for the color bar
  ))

pdf(file = here::here(fig_path, "intensity_mean_cv_sd.pdf"),
    width = twi, height = twi/2)
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
                         fill = q0.025,
                         colour = q0.025)) +
  scale_fill_viridis_c("",
                       limits = q_breaks,
                       breaks = q_breaks,
                       labels = q_labels)+
  scale_colour_viridis_c(limits = q_breaks,
                         breaks = q_breaks,
                         labels = q_labels,
                         guide = "none")+
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
                         fill = q0.975,
                         colour = q0.975)) +
  scale_fill_viridis_c("",
                       limits = q_breaks,
                       breaks = q_breaks,
                       labels = q_labels)+
  scale_colour_viridis_c(limits = q_breaks,
                         breaks = q_breaks,
                         labels = q_labels,
                         guide = "none")+
  xlab("Easting") +
  ylab("Northing") +
  theme_void() +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(fill = guide_colourbar(title = "",
                                title.vjust = 0.95,
                                title.theme = element_text(size = 18),
                                label.theme = element_text(size = 16)))

p3

p2 = p2 +
  theme(legend.margin = margin(t = 0, r = 10, b = 0, l = 10)) +
  guides(fill = guide_colourbar(
    title.vjust = 0.85,
    title.theme = element_text(size = ax.size-2),
    label.theme = element_text(size = ax.size-2),
    barwidth = unit(2, "cm"),  # Consistent width for the color bar
    barheight = unit(0.5, "cm") # Consistent height for the color bar
  ))

p3 = p3 +
  theme(legend.margin = margin(t = 0, r = 10, b = 0, l = 10)) +
  guides(fill = guide_colourbar(
    title.vjust = 0.85,
    title.theme = element_text(size = ax.size-2),
    label.theme = element_text(size = ax.size-2),
    barwidth = unit(2, "cm"),  # Consistent width for the color bar
    barheight = unit(0.5, "cm") # Consistent height for the color bar
  ))

library(patchwork)
pdf(file = here::here(fig_path, "intensity_quantiles.pdf"),
    width = twi/1.5, height = twi/2)
# plot_grid(NULL, p2, NULL, p3, NULL,
#           labels = c("", "A", "", "B", ""),
#           rel_widths = c(0.05, 1, 0.15, 1, 0.05),
#           ncol = 5)
wrap_plots(p2, p3) +
  plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(1, 1), # Adjust heights for equal space
            nrow = 1) +  # 2 rows for vertical stacking
  theme(plot.margin = margin(10, 10, 10, 10)) &
  theme(legend.justification = c(0.5, 0))  # Adjust margins for more whitespace
dev.off()

# show that these maps are misleading by interpreting
# them as an intensity
cell_area = st_area(study_area)/nrow(pr.int)*100
sum(pr.int$q0.025*cell_area)
sum(pr.int$q0.975*cell_area)
# these are not supported by posterior for abundance (see script N_posterior.R)

### plot three realisations of posterior intensity field ####

set.seed(19822139)
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

viridisscale = scale_fill_viridis_c("",
                                    limits = range(draw_all))

p1 = ggplot() +
  geom_sf(data = study_area) +
  geom_tile(data = draw1,
            mapping = aes(x = x,
                         y = y,
                         fill = mean,
                         colour = mean)) +
  viridisscale +
  scale_colour_viridis_c("",
                         limits = range(draw_all),
                         guide = "none") +
  xlab("Easting") +
  ylab("Northing") +
  theme_void() +
  theme(legend.position = "none")

p2 = ggplot() +
  geom_sf(data = study_area) +
  geom_tile(data = draw2,
            mapping = aes(x = x,
                         y = y,
                         fill = mean,
                         colour = mean)) +
  viridisscale +
  scale_colour_viridis_c("",
                         limits = range(draw_all),
                         guide = "none") +
  xlab("Easting") +
  ylab("Northing") +
  theme_void() +
  theme(legend.position = "none")

p3 = ggplot() +
  geom_sf(data = study_area) +
  geom_tile(data = draw3,
            mapping = aes(x = x,
                         y = y,
                         fill = mean,
                         colour = mean)) +
  viridisscale +
  scale_colour_viridis_c("",
                         limits = range(draw_all),
                         guide = "none") +
  xlab("Easting") +
  ylab("Northing") +
  theme_void() +
  theme(legend.position = "none")

pdf(file = here::here(fig_path, "intensity_realized.pdf"),
    width = twi, height = twi/2)
plot_grid(NULL, p1, NULL, p2, NULL, p3, NULL,
          labels = c("", "A", "", "B", "", "C", ""),
          rel_widths = c(0.1, 1, 0.4, 1, 0.4, 1, 0.1),
          ncol = 7)
dev.off()

#### Posterior half-normal plot ####

W = 58/1000
Wm = 58    # transect radius in metres
distdf <- data.frame(distance = seq(.Machine$double.eps, Wm, length=100))

# adjust for change of units (km in model fit, metres here)
hnpred <- predict(fit, distdf, ~ hn(distance, lsig + log(1000)),
                  n.samples = 100)

ghn = ggplot() +
        gg(hnpred) +
        ylab("detection probability\n") +
        xlab("\ndistance (m)") +
        theme(axis.title = element_text(size = ax.size),
        axis.text = element_text(size = ax.size))
ghn

##### Posterior Matern plot ####
blah = spde.posterior(fit, "grf", what = "matern.correlation")
gmat = ggplot() +
  gg(blah) +
  ylab("correlation\n") +
  xlab("\ndistance (km)") +
  theme(axis.title = element_text(size = ax.size, vjust = 0.3),
        axis.text = element_text(size = ax.size))

gmat

##### Pairwise distances plot ####
# polygon samplers
samplers$ID = 1:nrow(samplers)
W = 58/1000
samplers_buffered = st_buffer(samplers, dist = W)

### Generate posterior point patterns
n.pp = 500

# simulate intensities and detection functions
# each column of llam is a realisation of the log-intensity at each mesh node
post.sample = generate(fit,
                formula = ~ {c(
                  llam = grf_latent + Intercept_latent,
                  lsig = lsig_latent)
                },
                n.samples = n.pp)

llam = post.sample[1:mesh$n, ]
lsig.post = post.sample["lsig",]

# create distance bins and matrix to store counts
breaks = 0:12
counts = matrix(NA, nrow = n.pp, ncol = length(breaks)-1)

# detection functions
# Log half-normal
log_hn = function(distance, lsig){
  -0.5*(distance/exp(lsig))^2
}

# Half-normal
hn <- function(distance, lsig) exp(log_hn(distance, lsig))

# Need to create a new mesh just in the study area and
# project llam to this. Boundary was giving very large intensities
# when trying with full mesh and samplers argument

## Generated by meshbuilder()

## Build boundary information:
## (fmesher supports SpatialPolygons, but this app is not (yet) intelligent enough for that.)
boundary <- list(
  fm_as_segm(study_area),
  NULL)

## Build the mesh:
inner_mesh <- fm_mesh_2d_inla(boundary=boundary,
                              max.edge=c(0.25, 0.61),
                              min.angle=c(30, 21),
                              max.n=c(48000, 16000), ## Safeguard against large meshes.
                              max.n.strict=c(128000, 128000), ## Don't build a huge mesh!
                              cutoff=0.024, ## Filter away adjacent points.
                              offset=c(0.1, 0.13)) ## Offset for extra boundaries, if needed.

A = inla.spde.make.A(mesh = mesh,
                     loc = inner_mesh$loc)

inner_llam = A %*% llam

# Note:  I think there is a way to use inner_mesh in the generate()
# call above to avoid having to create A and project afterwards.
# Have not tested this though.

for (i in 1:n.pp){

  cat("\nGenerating point pattern ", i, "\n")
  a.pp = sample.lgcp(mesh = inner_mesh,
                     loglambda = inner_llam[,i],
                     samplers = study_area)
  a.pp = st_as_sf(a.pp)
  # keep only detectable points
  a.pp = a.pp[a.pp %>% st_within(samplers_buffered) %>% lengths > 0, ]
  a.lsig = lsig.post[i]

  cat("\nThinning point pattern ", i, "\n")

  # add transect ID for each point
  pp_det = st_join(a.pp,
                   samplers_buffered[,"ID"])

  distances = rep(NA, nrow(pp_det))

  for (j in 1:nrow(pp_det)){
    pt = pp_det[j,]
    transect = subset(samplers, ID == pt$ID)
    distances[j] = as.numeric(st_distance(pt, transect))
  }

  pdet = hn(distances, a.lsig)
  pp_det$distances = distances
  pp_det$pdet = pdet

  detected <- rbinom(nrow(pp_det), 1, pdet)
  a.obs <- subset(pp_det, as.logical(detected))

  # calculate pairwise distances
  a.ds = st_distance(a.obs, byid = TRUE)
  a.ds[upper.tri(a.ds, diag = TRUE)] = NA
  a.ds.vec = as.numeric(a.ds)
  a.ds.vec = a.ds.vec[!is.na(a.ds.vec)]
  counts[i,] = hist(a.ds.vec, breaks = breaks, plot = FALSE)$counts

  cat("\nFinished processing point pattern ", i, "\n")

  # memory leak?
  # rm(a.ds, a.ds.vec, a.obs, a.lsig, a.pp)
  # gc()
}

# pairwise distances for observed data
obs.ds = st_distance(obs)
obs.ds[upper.tri(obs.ds, diag = TRUE)] = NA
obs.ds.vec = as.numeric(obs.ds)
count.obs = hist(obs.ds.vec, breaks = breaks, plot = FALSE)$counts

# Prepare the counts data for ggplot
counts_df <- as.data.frame(counts)
counts_long <- counts_df %>%
  pivot_longer(cols = everything(), names_to = "distance_bin", values_to = "count") %>%
  mutate(distance_bin = as.numeric(gsub("V", "", distance_bin)))

# Prepare observed data as a separate data frame
count_obs_df <- data.frame(
  distance_bin = 1:(length(breaks) - 1),
  count = count.obs
)

# Define x-axis tick labels
xticks <- breaks[-1]

gpp = ggplot() +
  geom_boxplot(data = counts_long, aes(x = factor(distance_bin), y = count),
               fill = "white", outlier.shape = NA) +
  geom_point(data = count_obs_df, aes(x = distance_bin, y = count),
             color = "red", size = 2) +
  scale_x_discrete(labels = xticks) +
  labs(x = "distance (km)", y = "count of pairwise distances") +
  theme_minimal()

# Plot with ggplot2
pdf(file = here::here("figures", "dist_mat_pairs.pdf"),
    width = twi, height = twi)

(ghn + gmat) / gpp + plot_annotation(tag_levels = "A")

dev.off()
