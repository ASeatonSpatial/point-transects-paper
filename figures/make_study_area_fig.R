# Evaluate fitted model
library(INLA)
library(inlabru)

# set ggplot theme
theme_set(theme_minimal())

study_area = readRDS("../analyses/data/study_area_extended_no_crs.RDS")
samplers = readRDS("../analyses/data/samplers_extended_no_crs.RDS")

png(filename = "study_area_design.png",
    width = 10, height = 6, units = "in", res = 100)

ggplot() +
  gg(study_area) +
  gg(samplers) +
  coord_equal() +
  xlab("Easting") +
  ylab("Northing")

dev.off()
