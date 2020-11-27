# convert AnalysisDataset2002.csv to appropriate SpatialPointDataFrame objects

library(dplyr)
library(sp)
library(inlabru)   # for quick plotting

data_path = here::here("analyses", "data")
newdata = read.csv(here::here(data_path, "AnalysisDataset2002.csv"),
                   stringsAsFactors = FALSE)

# old_obs = readRDS("obs_2002_no_crs.RDS")
# old_samplers = readRDS("samplers_no_crs.RDS")

# Stratum codes:
# OF = open forest
# CF = closed forest
# RF = reforested    - no akepa here

newdata %>%
  select(SampleLabel, Easting, Northing, Stratum, HabCode) %>%
  filter(Stratum %in% c("OF", "CF")) %>%
  distinct() -> samplers

samplers$weight = 1
samplers$x = samplers$Easting/1000
samplers$y = samplers$Northing/1000
coordinates(samplers) = ~ x + y

obs = newdata %>%
  select(SampleLabel, Distance, Easting, Northing, Stratum, HabCode) %>%
  filter(Stratum %in% c("OF", "CF"), !is.na(Distance)) %>%
  rename(distance = Distance)

obs$x = obs$Easting/1000
obs$y = obs$Northing/1000
obs$distance = obs$distance / 1000

coordinates(obs) = ~ x + y

# check it looks okay

ggplot() +
  gg(samplers) +
  gg(obs, col = "green")

saveRDS(samplers, file = here::here(data_path, "samplers_extended_no_crs.RDS"))
saveRDS(obs, file = here::here(data_path, "obs_extended_no_crs.RDS"))
