# Some summaries of the data to go in the paper

library(INLA)
library(inlabru)

data_path <- here::here("R", "data")
study_area <- readRDS(here::here(data_path, "study_area.RDS"))
samplers <- readRDS(here::here(data_path, "samplers.RDS"))
mesh <- readRDS(here::here(data_path, "mesh.RDS"))
obs <- readRDS(here::here(data_path, "obs.RDS"))

# number of transects
nrow(samplers)

# number of detections within truncation distance
nrow(obs)

# number of transects with at least one detection
length(unique(obs$SampleLabel))

# range of counts at each sampling unit
library(dplyr)
hmm <- data.frame(SampleLabel = obs$SampleLabel)
hmm %<>%
  group_by(SampleLabel) %>%
  summarise(n = n())
rbind(n = seq_len(max(hmm$n)),
      count = tabulate(hmm$n))

# Comparison of hazard-rate and half-norm detection functions.
plot.ecdf(obs$distance)
# Hand-picked scalings to match most of the empirical CDF
# Using posterior means of sig and gam for hr():
curve(cumsum(hr(x, 0.0398, 4.906) * x) * 0.57,
      0,
      58 / 1000,
      add = TRUE,
      col = 4)
# Hand-picked sig value for hn():
curve(cumsum(hn(x, 0.06) * x) * 0.64, add = TRUE, col = 2)
legend(
  "topleft",
  legend = c("Empirical", "Hazard-rate", "Half-normal"),
  col = c(1, 4, 2),
  lty = 1
)
