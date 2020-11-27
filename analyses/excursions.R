#### Try excursions package ####
library(INLA)
library(inlabru)
library(excursions)
source("gg.R")   # small edit of gg.SpatialPixelsDataFrame

## todo - replace predict loop with generate() ?

# set ggplot theme
theme_set(theme_minimal())

set.seed(1525)

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

#### monte-carlo samples of intensity ####

# Does excursions work with an inlabru object?

# There must be a way to get the actual values from
# predict with n.samples > 1 ?

# yes update this to use generate

n.mc = 500     # lower this when fiddling with things
pxl = pixels(mesh, nx = 100, ny = 100, mask = study_area)
X = matrix(NA, nrow = nrow(pxl@coords), ncol = n.mc)
for (i in 1:n.mc){
  print(paste("i = ", i))
  int_draw <- predict(fit, pxl, ~ exp(grf + Intercept), n.samples = 1)
  X[,i] = int_draw@data$mean
}

#### Try excursions ####

# re-scale X.  Currently per m^2
# per km^2:
Xkm = X*1000000

# per hectare:
Xhec = X*10000

# per hectare threshold:
u = 1
alpha = 0.05    # Prob(intensity > u) >= 1 - alpha for all s in E

ex = excursions.mc(Xhec, u = u, type = ">", alpha = alpha)

# plot
Fpxl = pxl
Fpxl$F = ex$F

# png(filename = "../figures/excursion_function.png",
#     width = 4.5, height = 5.5, units = "cm", res = 1000)
p2 = ggplot() +
  gg(study_area) +
  gg(Fpxl) +
  scale_fill_viridis_c(limits = c(0,1),
                       breaks = c(0, 0.5, 1),
                       labels = c(0, 0.5, 1))+
  coord_equal() +
  theme_void() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal") +
        #plot.margin = margin(0.1, 0, 0.1, 0, "in")) +
  guides(fill = guide_colourbar(title.vjust = 0.95,
                                title.theme = element_blank(),
                                label.theme = element_text(size = 16, hjust = 0.5)))

p2

# ggsave(filename = "../figures/excursion_function.png",
#        width = 5, height = 6, units = "in")

Epxl = pxl
Epxl$E = as.factor(ex$E)
Epxl = Epxl[Epxl$E == 1,]
p1 = ggplot() +
  gg(study_area) +
  gg(Epxl, aes(group = E)) +
  coord_equal() +
  # ggtitle("Excursion set for > 1 bird per hectare,\nalpha = 0.05") +
  scale_fill_manual(values = c("blue"),
                    labels = c("Excursion set")) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        plot.margin = margin(0.1, 0, 0.3, 0, "in"),
        legend.box.margin = margin(-4,0,0,0)) +
  guides(fill = guide_legend(title.vjust = 0.95,
                             title.theme = element_blank(),
                             label.theme = element_text(size = 16, hjust = 0.5)))

p1

# ggsave(filename = "../figures/excursion_set.png",
#        width = 5, height = 6, units = "in")

png(filename = here::here(fig_path, "excursions.png"),
    width = 10, height = 5, units = "in", res = 100)
plot_grid(NULL, p1, NULL, p2, NULL,
          labels = c("", "A", "", "B", ""),
          rel_widths = c(0.05, 1, 0.15, 1, 0.05),
          ncol = 5)
dev.off()

# contour map - how to visualise this?  ?
cm = contourmap.mc(Xhec, levels = c(0.1, 0.5, 0.8, 1.5, 2, 4), alpha = alpha)
str(cm)

Gpxl = pxl
Gpxl$G = cm$G
ggplot() +
  gg(Gpxl)

cm = contourmap.mc(Xhec, levels = c(1), alpha = alpha)
Gpxl = pxl
Gpxl$G = cm$G
ggplot() +
  gg(Gpxl)
