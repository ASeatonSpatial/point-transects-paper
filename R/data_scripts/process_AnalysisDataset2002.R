# process AnalysisDataset2002.csv
# into points sf object

library(dplyr)
library(sf)
library(ggplot2)

data_path = here::here("R", "data_raw")
newdata = read.csv(here::here(data_path, "AnalysisDataset2002.csv"),
                   stringsAsFactors = FALSE)

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
samplers_sf = st_as_sf(samplers,
                       coords = c("x", "y"))

obs = newdata %>%
  select(SampleLabel, Distance, Easting, Northing, Stratum, HabCode) %>%
  filter(Stratum %in% c("OF", "CF"), !is.na(Distance)) %>%
  rename(distance = Distance)

obs$x = obs$Easting/1000
obs$y = obs$Northing/1000
obs$distance = obs$distance / 1000

obs_sf = st_as_sf(obs,
                  coords = c("x", "y"))

# check it looks okay
ggplot() +
  geom_sf(data = samplers_sf) +
  geom_sf(data = obs_sf, col = "green")

# write to disk
out_path = here::here("R", "data")
saveRDS(samplers_sf, file = here::here(out_path, "samplers.RDS"))
saveRDS(obs_sf, file = here::here(out_path, "obs.RDS"))
