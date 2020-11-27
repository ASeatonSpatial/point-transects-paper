# Evaluate fitted model
library(INLA)
library(inlabru)

data_path = here::here("analyses", "data")
study_area = readRDS(here::here(data_path, "study_area_extended_no_crs.RDS"))
samplers = readRDS(here::here(data_path, "samplers_extended_no_crs.RDS"))
mesh = readRDS(here::here(data_path, "mesh_extended_no_crs.RDS"))

# number of transects
nrow(samplers@coords)

# number of detections
nrow(obs@coords)

# number of transects with at least one detection
length(unique(obs$SampleLabel))

# range of counts at each sampling unit
library(dplyr)
hmm = data.frame(SampleLabel = obs$SampleLabel)
hmm %<>%
  group_by(SampleLabel) %>%
  summarise(n = n())
max(hmm$n)
