#### Try excursions package ####
library(INLA)
library(inlabru)
library(ggplot2)
library(excursions)
library(cowplot)
library(sf)
source(here::here("R", "gg.R"))   # small edit of gg.SpatialPixelsDataFrame

# set ggplot theme
theme_set(theme_minimal())

set.seed(1525)

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

#### monte-carlo samples of intensity ####

n.mc = 500     # lower this when fiddling with things
n.mc = 25     # lower this when fiddling with things
pxl = fm_pixels(mesh,
                dims = c(50, 50),
                mask = study_area)

pxl.coords = st_coordinates(pxl)
pxl$x = pxl.coords[,1]
pxl$y = pxl.coords[,2]

X = generate(fit,
             pxl,
             ~ exp(grf + Intercept),
             n.samples = n.mc)

#### excursions ####

# re-scale X (per km) to per hectare:
Xhec = X/100

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
  geom_sf(data = study_area) +
  geom_tile(data = Fpxl,
            mapping = aes(x = x,
                          y = y,
                          fill = F)) +
  scale_fill_viridis_c(limits = c(0,1),
                       breaks = c(0, 0.5, 1),
                       labels = c(0, 0.5, 1))+
  theme_void() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal") +
        #plot.margin = margin(0.1, 0, 0.1, 0, "in")) +
  guides(fill = guide_colourbar(title.vjust = 0.95,
                                title.theme = element_blank(),
                                label.theme = element_text(size = 16, hjust = 0.5)))

p2

Epxl = pxl
Epxl$E = as.factor(ex$E)
Epxl = Epxl[Epxl$E == 1,]
p1 = ggplot() +
  geom_sf(data = study_area) +
  geom_tile(data = Epxl,
            mapping = aes(x = x,
                          y = y,
                          fill = E)) +
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
