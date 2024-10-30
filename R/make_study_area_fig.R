# Evaluate fitted model
library(INLA)
library(inlabru)
library(sf)
library(ggplot2)

theme_set(theme_minimal())

# statsoc doc class textwidth in cm
tw = 14.69785
twi = 5.7865551181  # pdf() needs inches
ax.size = 12 # axis text size

data_path = here::here("R", "data")
study_area = readRDS(here::here(data_path, "study_area.RDS"))
samplers = readRDS(here::here(data_path, "samplers.RDS"))
study_area_boundaries <- sf::st_cast(study_area, "LINESTRING")

pdf(file = here::here("figures", "study_area_design.pdf"),
    width = twi/2, height = twi/2)

ggplot() +
  geom_sf(data = study_area_boundaries,
          # fill = NA,
          size = 2,
          colour = "red") +
  geom_sf(data = samplers) +
  xlab("Easting (km)") +
  ylab("Northing (km)") +
  scale_x_continuous(limits = c(254, 262),
                     breaks = seq(254, 262, by = 2)) +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = ax.size),
    axis.text.x = element_text(angle = -45, vjust = 0, size = ax.size), # Rotate x-axis labels
    axis.title = element_text(size = ax.size),
    legend.title = element_text(size = ax.size, vjust = 0.75),
    legend.text = element_text(size = ax.size),
    legend.key.size = unit(2, "line")
  )

dev.off()
