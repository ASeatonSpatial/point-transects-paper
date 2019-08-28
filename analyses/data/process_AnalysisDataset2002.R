# convert AnalysisDataset2002.csv to appropriate SpatialPointDataFrame objects

library(dplyr)
library(sp)
library(inlabru)   # for quick plotting

newdata = read.csv("AnalysisDataset2002.csv", stringsAsFactors = FALSE)

old_obs = readRDS("obs_2002_no_crs.RDS")
old_samplers = readRDS("samplers_no_crs.RDS")

# Stratum codes:
# OF = open forest
# CF = closed forest
# RF = reforested    - no akepa here

newdata %>%
  select(SampleLabel, Easting, Northing, Stratum) %>% 
  filter(Stratum %in% c("OF", "CF")) %>%
  distinct() -> samplers

samplers$weight = 1
coordinates(samplers) = ~ Easting + Northing

obs = newdata %>% 
  select(SampleLabel, Distance, Easting, Northing, Stratum) %>% 
  filter(Stratum %in% c("OF", "CF"), !is.na(Distance)) -> obs

coordinates(obs) = ~ Easting + Northing

# check it looks okay

ggplot() +
  gg(samplers) +
  gg(obs, col = "green")

saveRDS(samplers, file = "samplers_extended_no_crs.RDS")
saveRDS(samplers, file = "obs_extended_no_crs.RDS")
